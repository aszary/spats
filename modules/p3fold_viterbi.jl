module P3FoldViterbi

using Statistics
using FFTW
using DSP
using SparseArrays
using LinearAlgebra

# Background / design rationale: p3fold-refine-notes.md §3.1-3.2.
#
# psrsalsa's `pfold` refine mode (foldP3 in fold.c) assigns phase block-by-block,
# greedily matching each block only against the sum of *previous* blocks
# (fold.c:158-176). A bad block can derail everything after it, and the
# correlation measure (template²·blockmap²) isn't normalized, so loud blocks
# dominate the phase choice.
#
# Here the phase is assigned per *pulse* (not per block) and globally via
# Viterbi: the whole observation is optimized at once, so one noisy pulse
# can't permanently knock later pulses off phase. The emission score is a
# normalized (zero-mean, unit-power) Pearson correlation, so faint and loud
# pulses are weighted by how well they match the template shape, not by their
# raw amplitude. The transition cost penalizes deviation from the phase
# increment implied by the nominal P3 but never forbids a jump outright, so
# genuine P3 changes can still be tracked (controlled by `continuity_weight`).


"""
    emission_corr(data, template, on) -> Matrix{Float64}

Per-pulse, per-phase-bin normalized correlation matrix (N_pulses × ybins).

`C[i, j]` is the zero-mean Pearson correlation between pulse `i`'s on-pulse
window (`on`) and phase-bin `j`'s on-pulse window in `template`. Pulses or
phase bins with zero variance (e.g. nulled pulses, still-empty template
bins) get correlation 0 against everything, so they neither win nor actively
mislead the Viterbi search.
"""
function emission_corr(data::AbstractMatrix, template::AbstractMatrix, on::UnitRange)
    ybins = size(template, 1)
    N = size(data, 1)

    tmpl_on = template[:, on]
    tmean = vec(mean(tmpl_on, dims=2))
    tc = tmpl_on .- tmean
    tnorm = [sqrt(sum(abs2, @view tc[j, :])) for j in 1:ybins]

    C = zeros(Float64, N, ybins)
    for i in 1:N
        d = @view data[i, on]
        dc = d .- mean(d)
        dn = sqrt(sum(abs2, dc))
        dn == 0 && continue
        for j in 1:ybins
            tnorm[j] == 0 && continue
            C[i, j] = sum(dc .* @view tc[j, :]) / (dn * tnorm[j])
        end
    end
    return C
end


"""
    circ_residual(x, m) -> Float64

Wrap `x` into (-m/2, m/2], i.e. the signed residual of `x` modulo `m`
closest to 0. Used to measure how far an actual phase step is from the
expected one, on a circle of `m` phase bins.
"""
function circ_residual(x::Real, m::Real)
    return mod(x + m / 2, m) - m / 2
end


"""
    viterbi(C, delta, continuity_weight) -> Vector{Int}

Globally optimal phase-bin path (1..ybins per pulse) given the emission
matrix `C` (N_pulses × ybins, see `emission_corr`).

Transition cost between consecutive pulses is
`continuity_weight * circ_residual(j - j', delta, ybins)^2`, i.e. a quadratic
penalty on deviating from the expected phase increment `delta = ybins / p3`.
`continuity_weight = 0` makes the search flat (any phase jump equally
"free", pure per-pulse template matching); larger values increasingly
enforce smooth, P3-like drift. This replaces foldP3's hard one-cycle-per-block
limit with a soft, tunable preference (notes §3.1).
"""
function viterbi(C::AbstractMatrix, delta::Real, continuity_weight::Real)
    N, ybins = size(C)
    cost = fill(Inf, N, ybins)
    backptr = zeros(Int, N, ybins)

    cost[1, :] .= -C[1, :]

    for i in 2:N
        for j in 1:ybins
            bestc = Inf
            bestp = 1
            for jp in 1:ybins
                resid = circ_residual((j - jp) - delta, ybins)
                trans = continuity_weight * resid^2
                c = cost[i-1, jp] + trans
                if c < bestc
                    bestc = c
                    bestp = jp
                end
            end
            cost[i, j] = bestc - C[i, j]
            backptr[i, j] = bestp
        end
    end

    path = zeros(Int, N)
    path[N] = argmin(@view cost[N, :])
    for i in N-1:-1:1
        path[i] = backptr[i+1, path[i+1]]
    end
    return path
end


"""
    build_template(data, path, ybins) -> Matrix{Float64}

Sum pulses into their assigned phase bin, ybins × N_bins (unnormalized sum,
same convention as `Tools.p3fold`).
"""
function build_template(data::AbstractMatrix, path::AbstractVector{Int}, ybins::Int)
    nb = size(data, 2)
    template = zeros(Float64, ybins, nb)
    for i in eachindex(path)
        template[path[i], :] .+= @view data[i, :]
    end
    return template
end


"""
    circ_unwrap_steps(path, ybins, delta) -> Vector{Float64}

Per-pulse phase increment `path[i+1] - path[i]`, unwrapped: each step is
resolved to the representative value closest to the expected nominal
increment `delta`, removing the ±ybins circular ambiguity. Length N-1.
"""
function circ_unwrap_steps(path::AbstractVector{Int}, ybins::Int, delta::Real)
    N = length(path)
    steps = zeros(Float64, N - 1)
    for i in 1:N-1
        raw = path[i+1] - path[i]
        steps[i] = delta + circ_residual(raw - delta, ybins)
    end
    return steps
end


"""
    instantaneous_p3(path, ybins, p3_nominal; window=20) -> Vector{Float64}

Smoothed, per-pulse instantaneous P3 [pulse periods] implied by the Viterbi
phase path — this is "the P3 actually used" for each pulse, as opposed to
the single nominal `p3` fed into `fold`.

The discrete per-pulse phase increments are unwrapped into a continuous
phase track (`circ_unwrap_steps`), then a local linear fit over a sliding
window of `window` pulses gives the local slope dphase/dpulse, from which
`P3 = ybins / slope`. Pulses too close to either end for a full window keep
`p3_nominal`.
"""
function instantaneous_p3(path::AbstractVector{Int}, ybins::Int, p3_nominal::Real; window::Int=20)
    delta = ybins / p3_nominal
    N = length(path)
    steps = circ_unwrap_steps(path, ybins, delta)
    cum = vcat(0.0, cumsum(steps))   # cum[i] = unwrapped phase relative to pulse 1

    p3_inst = fill(Float64(p3_nominal), N)
    halfw = max(1, window ÷ 2)
    for i in 1:N
        lo = max(1, i - halfw)
        hi = min(N, i + halfw)
        hi - lo < 2 && continue
        xs = lo:hi
        ys = @view cum[lo:hi]
        xbar = mean(xs)
        ybar = mean(ys)
        num = sum((x - xbar) * (y - ybar) for (x, y) in zip(xs, ys))
        den = sum((x - xbar)^2 for x in xs)
        slope = den == 0 ? delta : num / den
        p3_inst[i] = slope == 0 ? p3_nominal : ybins / slope
    end
    return p3_inst
end


"""
    fold(data, p3, bin_st, bin_end; ybins, n_iter, continuity_weight) -> NamedTuple

Refined P3-fold via per-pulse Viterbi phase assignment, intended as a
globally-optimized alternative to `pfold -p3fold`'s block-greedy refine
(p3fold-refine-notes.md §3).

Bootstraps from the same fixed-phase fold as `Tools.p3fold` (so `n_iter=0`
reproduces it exactly), then alternates:
  1. score every pulse against every phase bin of the current template
     (`emission_corr`),
  2. find the globally optimal phase path (`viterbi`),
  3. rebuild the template by summing pulses along that path,
for `n_iter` EM-like rounds.

Arguments:
  data     – single-pulse matrix (N_pulses × N_bins), real intensity
  p3       – nominal P3 [pulse periods P0]; only used to seed the bootstrap
             fold and to set the expected phase increment per pulse
             (`delta = ybins / p3`) for the continuity penalty
  bin_st, bin_end – on-pulse window (1-indexed) used for the correlation
             score; the folded output still spans all N_bins
  ybins    – number of P3-phase bins (states), default 10
  n_iter   – number of refit rounds, default 5
  continuity_weight – transition penalty strength; 0 = flat (any phase
             jump equally likely, pure template matching), default 0.2
  p3_window – smoothing window [pulses] for `p3_per_pulse`, default 20

Returns:
  folded       – ybins × N_bins matrix, the refined p3-fold
  phase        – Vector{Int}, assigned phase bin (1..ybins) per pulse
  confidence   – Vector{Float64}, correlation of each pulse against its
                 assigned phase bin (low ⇒ poor match, e.g. nulled pulse)
  margin       – Vector{Float64}, correlation gap between the best and
                 second-best phase bin per pulse (low ⇒ ambiguous phase,
                 notes §3.4)
  p3_per_pulse – Vector{Float64}, smoothed instantaneous P3 [pulse periods]
                 implied by `phase` (see `instantaneous_p3`) — the local P3
                 the Viterbi path is actually tracking at each pulse
"""
function fold(data::AbstractMatrix, p3::Real, bin_st::Int, bin_end::Int;
              ybins::Int=10, n_iter::Int=5, continuity_weight::Real=0.2,
              p3_window::Int=20)
    N = size(data, 1)
    on = bin_st:bin_end
    delta = ybins / p3

    path = [floor(Int, mod(i * delta, ybins)) + 1 for i in 1:N]
    template = build_template(data, path, ybins)

    C = zeros(Float64, N, ybins)
    for _ in 1:n_iter
        C = emission_corr(data, template, on)
        path = viterbi(C, delta, continuity_weight)
        template = build_template(data, path, ybins)
    end

    confidence = zeros(Float64, N)
    margin = zeros(Float64, N)
    if n_iter > 0
        for i in 1:N
            row = @view C[i, :]
            confidence[i] = row[path[i]]
            best, second = partialsort(row, 1:2, rev=true)
            margin[i] = best - second
        end
    end

    p3_per_pulse = instantaneous_p3(path, ybins, p3; window=p3_window)

    return (folded=template, phase=path, confidence=confidence, margin=margin,
            p3_per_pulse=p3_per_pulse)
end


"""
    windowed_slope(y, window) -> Vector{Float64}

Local slope dy/dx of `y` (sampled at integer x=1..length(y)), from a
sliding linear-regression window of `window` samples centered at each
point (the window shrinks near the edges).
"""
function windowed_slope(y::AbstractVector, window::Int)
    N = length(y)
    slope = zeros(Float64, N)
    halfw = max(1, window ÷ 2)
    for i in 1:N
        lo = max(1, i - halfw)
        hi = min(N, i + halfw)
        if hi - lo < 2
            slope[i] = i > 1 ? y[i] - y[i-1] : (N > 1 ? y[2] - y[1] : 0.0)
            continue
        end
        xs = lo:hi
        ys = @view y[lo:hi]
        xbar = mean(xs)
        ybar = mean(ys)
        num = sum((x - xbar) * (yy - ybar) for (x, yy) in zip(xs, ys))
        den = sum((x - xbar)^2 for x in xs)
        slope[i] = den == 0 ? 0.0 : num / den
    end
    return slope
end


"""
    second_difference_matrix(N) -> SparseMatrixCSC

(N-2) × N second-difference operator: `(D*x)[i] = x[i+2] - 2x[i+1] + x[i]`,
i.e. the discrete curvature of `x` at each interior point. Building block
for `whittaker_smooth`.
"""
function second_difference_matrix(N::Int)
    I_idx = Int[]; J_idx = Int[]; V = Float64[]
    for i in 1:N-2
        push!(I_idx, i); push!(J_idx, i);   push!(V, 1.0)
        push!(I_idx, i); push!(J_idx, i+1); push!(V, -2.0)
        push!(I_idx, i); push!(J_idx, i+2); push!(V, 1.0)
    end
    return sparse(I_idx, J_idx, V, N - 2, N)
end


"""
    whittaker_smooth(x, λ) -> Vector{ComplexF64}

Curvature-penalized smoother (Whittaker 1923 / Hodrick-Prescott filter),
used as `coherent_fold`'s cutoff-free alternative to a fixed Butterworth
lowpass. The smoothed series `xs` is the minimiser of

    ‖x - xs‖² + λ·‖D·xs‖²

where `D` is `second_difference_matrix(N)`. Directly penalising curvature
means the penalty is on how much the smoothed track is allowed to *bend*
from sample to sample — matching the physical prior that P3 changes
gradually, no discontinuous jumps — rather than on which Fourier
frequencies survive, which is what a lowpass cutoff controls instead and
doesn't map onto "no jumps" as directly.

`λ = 0` reproduces `x` exactly; `λ → ∞` collapses `xs` to a straight line
(zero curvature). Solved as one sparse linear system
`(I + λD'D) xs = x`, applied to the real and imaginary parts separately.
"""
function whittaker_smooth(x::AbstractVector{<:Complex}, λ::Real)
    N = length(x)
    D = second_difference_matrix(N)
    A = sparse(I, N, N) + λ .* (D' * D)
    xs_re = A \ real.(x)
    xs_im = A \ imag.(x)
    return xs_re .+ im .* xs_im
end


"""
    whittaker_smooth_auto(x; strength=1.0) -> Vector{ComplexF64}

Picks `λ` for `whittaker_smooth` automatically via the discrepancy
principle, instead of scanning a hand-picked grid.

Estimates the noise level directly from `x`'s own spectrum: for a
*complex* time series, low |frequency| content (DC and its near
neighbours) sits at both ends of the DFT bin index range, so the middle
third of the spectrum (farthest from DC, around the Nyquist frequency) is
the part assumed to be noise-dominated (the same "real P3 wobble is slow"
assumption every filter here makes). By Parseval (`sum|X_k|² = N·sum|x_n|²`
for Julia's unnormalized `fft`), the median power in that band is itself
already a direct estimate of the total residual energy `Σ|noise_n|²`
one would expect if `x` were pure noise — this times `strength` is
`target` below.

`λ` is then increased until the smoother's residual `‖x - xs(λ)‖²`
matches `target`: too little smoothing (small λ) leaves genuine noise
unremoved (residual too small), too much (large λ) starts eating real
signal (residual too large) — the crossing point is where the smoothed
track has absorbed the noise and, as much as possible, nothing else.
Found by bisection in log(λ), since the residual grows monotonically
with λ.

`strength` scales that target directly: `strength=1` (default) accepts the
noise-floor estimate as-is; `strength<1` demands a smaller residual, i.e.
less smoothing than the discrepancy principle's own estimate calls for
(lets more real structure — and more noise — through, so the track can
still bend where the noise floor was overestimated); `strength>1` demands
more, i.e. more smoothing than the estimate. Use this to correct for the
estimate being systematically too aggressive or too loose for a given
pulsar's actual noise, without abandoning the automatic λ search entirely.
"""
function whittaker_smooth_auto(x::AbstractVector{<:Complex}; strength::Real=1.0)
    N = length(x)
    N < 9 && return whittaker_smooth(x, 1.0)

    X = fft(x)
    P = abs2.(X)
    lo_band = N ÷ 3
    hi_band = N - lo_band
    target = strength * median(@view P[lo_band+1:hi_band])

    residual(logλ) = sum(abs2, x .- whittaker_smooth(x, exp(logλ))) - target

    lo, hi = -4.0, 14.0
    flo, fhi = residual(lo), residual(hi)
    # widen the bracket if the crossing isn't inside [lo, hi] yet
    while flo > 0 && lo > -20.0
        lo -= 4.0; flo = residual(lo)
    end
    while fhi < 0 && hi < 30.0
        hi += 4.0; fhi = residual(hi)
    end

    for _ in 1:25
        mid = (lo + hi) / 2
        fmid = residual(mid)
        if fmid < 0
            lo = mid
        else
            hi = mid
        end
    end
    return whittaker_smooth(x, exp((lo + hi) / 2))
end


"""
    coherent_fold(data, p3, bin_st, bin_end; ybins, lowpass_cutoff, filter_order, p3_window) -> NamedTuple

Coherent, matched-filter P3-fold — an alternative to `fold` (per-pulse
Viterbi) for pulsars where the per-pulse modulation depth is far below the
noise floor but the *aggregate* signal is highly significant (see
p3fold-refine-notes.md §3.3, "niezależny estymator fazy: transformacja
Hilberta", and the conversation that motivated this: `PhaseDrift.drift_test`
can detect a coherent phase slope at tens of σ even when blind per-pulse
template matching, i.e. `fold`, has essentially nothing to lock onto).

No per-pulse blind matching is attempted here. Instead:
  1. a single, full-length (non-segmented) FFT gives the complex spatial
     template `L_on` at f3 = 1/p3 — the same un-chunked estimate
     `PhaseDrift.drift_test` uses, and for the same reason: chunking (like
     `pspec`'s segmented LRFS, see `Data.twodfs_lrfs`) requires phase
     coherence *between* segments, which P3 wobble destroys; a single
     global FFT only ever needs coherence *within* one frequency bin, which
     survives mild wobble.
  2. every pulse is projected onto `conj(L_on)` — a spatial matched filter
     using *all* on-pulse bins, weighted optimally — giving one high-SNR
     complex number per pulse instead of relying on raw per-pulse shape
     correlation.
  3. that per-pulse series is coherently demodulated at f3 and low-pass
     filtered (`DSP.filtfilt`, zero-phase) — a *sliding*, unsegmented
     analogue of pspec's blocked averaging, so there are no hard block
     boundaries for P3 wobble to decohere across.
  4. the slowly-varying phase left after filtering is added back onto the
     f3 ramp to get each pulse's absolute P3-phase directly, which drives
     the fold — replacing blind per-pulse correlation with a directly
     measured, high-SNR phase.

Arguments:
  data     – single-pulse matrix (N_pulses × N_bins), real intensity
  p3       – nominal P3 [pulse periods]; sets the demodulation frequency f3 = 1/p3
  bin_st, bin_end – on-pulse window (1-indexed)
  ybins    – number of P3-phase bins for the output fold, default 10
  lowpass_cutoff – low-pass cutoff [cycles/pulse] applied after
             demodulation; sets the fastest P3 wobble that can still be
             tracked (higher = more responsive to fast change but noisier;
             lower = smoother but assumes more stable P3), default 1/200
  filter_order – Butterworth filter order for the low-pass, default 4
             (ignored when `auto=true`)
  p3_window – local-linear-regression window [pulses] used to turn the
             (possibly still slightly noisy) smoothed phase into a slope,
             i.e. `p3_per_pulse` — this, not `lowpass_cutoff`/`auto`, is
             what controls how jumpy the *displayed* P3(t) track looks:
             differentiating over a short window amplifies whatever noise
             is left, however clean the underlying folded image is. Larger
             = smoother P3(t) display but slower to follow genuine fast
             changes. Default 60 (was 20 — that visibly amplified noise
             into jumps even when the fold image itself looked clean).
  bin_search – half-width [FFT bins] of the neighbourhood around the
             nominal-p3 bin searched for the true on-pulse power peak
             before building `L_on` (see step 1 below), default 2
  auto     – if true, replace the fixed Butterworth lowpass with the
             curvature-penalized `whittaker_smooth_auto` (directly enforces
             a smooth, non-jumpy P3(t) track, with its own smoothing
             strength λ picked automatically) and ignore `lowpass_cutoff` /
             `filter_order` entirely — no cutoff to pick or scan for,
             default false

Returns:
  folded       – ybins × N_bins matrix, the coherently-refolded p3-fold
  phase        – Vector{Float64}, total unwrapped P3-phase per pulse [rad]
  bin          – Vector{Int}, assigned phase bin (1..ybins) per pulse
  p3_per_pulse – Vector{Float64}, instantaneous P3 [pulse periods] from the
                 local slope of `phase`
  snr          – matched-filter detection significance (≈ the same
                 quantity `drift_test` reports) — a sanity check that there
                 is signal to track at all before trusting the fold
"""

function coherent_fold(data::AbstractMatrix, p3::Real, bin_st::Int, bin_end::Int;
                        ybins::Int=10, lowpass_cutoff::Real=1/200, filter_order::Int=4,
                        p3_window::Int=60, bin_search::Int=2, warn_weak::Bool=true,
                        auto::Bool=false, auto_strength::Real=1.0)
    N = size(data, 1)
    on = bin_st:bin_end
    f3 = 1.0 / p3

    # 1. single, full-length (non-segmented) complex spatial template at f3.
    # The nominal p3 rarely lands the true f3 exactly on the bin
    # `round(N/p3)` predicts — a fraction-of-a-bin mismatch measurably
    # weakens the matched filter (L_on) built from it. Refine by searching
    # the `±bin_search` neighbourhood for the bin with the most on-pulse
    # power, i.e. the one the real signal actually sits in.
    F = fft(data, 1)
    k0 = clamp(round(Int, N / p3), 1, N ÷ 2)
    k_lo = max(1, k0 - bin_search)
    k_hi = min(N ÷ 2, k0 + bin_search)
    k = k_lo + argmax([sum(abs2, @view F[kk+1, on]) for kk in k_lo:k_hi]) - 1
    L = F[k+1, :]
    L_on = L[on]

    off = vcat(1:bin_st-1, bin_end+1:size(data, 2))
    sigma_off = isempty(off) ? 0.0 : std([real.(L[off]); imag.(L[off])])
    snr = sigma_off == 0 ? Inf : sqrt(sum(abs2, L_on)) / (sigma_off * sqrt(length(on)))
    # noise floor: this is an RMS-of-magnitude statistic (unlike drift_test's
    # mean-of-magnitude, floor ~1.25), so pure complex Gaussian noise gives
    # E[RMS|z|]/sigma = sqrt(2) ≈ 1.41, not 1.25 — same idea as drift_test's
    # own snr<3 warning, just referenced to this statistic's own noise floor
    warn_weak && snr < 3 && @warn "Weak f3 feature in coherent_fold (snr = $(round(snr, digits=2)), " *
                     "pure noise gives ~1.41): the fold below is tracking a marginal " *
                     "signal — cross-check against PhaseDrift.drift_test's `significance` " *
                     "(a different, Monte-Carlo-calibrated quantity) before trusting it"

    # 2. spatial matched-filter projection: one high-SNR complex number per pulse.
    # The static (non-modulated) average profile must be removed first — it lives at
    # frequency 0 and, unlike in `drift_test` (which reads a single FFT bin and never
    # mixes frequencies), a time-domain projection like this carries every frequency
    # through, so the huge DC term would otherwise swamp the low-pass filter below.
    on_data = data[:, on]
    on_data_demeaned = on_data .- mean(on_data, dims=1)
    w = conj.(L_on)
    z = on_data_demeaned * w

    # 3. coherent demodulation at f3, then smooth the baseband (sliding, no block
    # edges). `auto=true` replaces the fixed Butterworth lowpass with the
    # curvature-penalized `whittaker_smooth_auto` — no `lowpass_cutoff` or
    # `filter_order` to choose at all.
    n = 1:N
    carrier = exp.((-1im * 2π * f3) .* n)
    baseband = z .* carrier
    if auto
        baseband_smooth = whittaker_smooth_auto(baseband; strength=auto_strength)
    else
        respf = digitalfilter(Lowpass(lowpass_cutoff), Butterworth(filter_order); fs=1.0)
        baseband_smooth = filtfilt(respf, real.(baseband)) .+ im .* filtfilt(respf, imag.(baseband))
    end

    # 4. residual phase -> total phase -> fold-bin assignment
    resid = DSP.unwrap(angle.(baseband_smooth))
    phase_total = (2π * f3) .* n .+ resid

    bin = [Int(floor(mod(phase_total[i] / (2π) * ybins, ybins))) + 1 for i in 1:N]
    template = build_template(data, bin, ybins)

    slope = windowed_slope(phase_total, p3_window)
    p3_per_pulse = (2π) ./ slope

    return (folded=template, phase=phase_total, bin=bin, p3_per_pulse=p3_per_pulse, snr=snr)
end


"""
    coherent_fold_jackknife(data, p3, bin_st, bin_end; ybins, lowpass_cutoff, filter_order,
                             p3_window, n_groups) -> NamedTuple

Empirical error bars on `coherent_fold`'s `p3_per_pulse` and `phase`, via
the bootstrap/jackknife idea from p3fold-refine-notes.md §3.4 ("wielkość
skoku fazy vs typowy rozrzut fazy w spokojnych, niewątpliwych odcinkach").

The on-pulse window is split into `n_groups` disjoint, contiguous longitude
sub-ranges. Different longitude bins carry independent detector noise, so
running `coherent_fold` separately on each sub-range gives `n_groups`
*independent* measurements of the same underlying phase track. Their
spread at each pulse — divided by `√n_groups` — estimates the uncertainty
of the full-bin estimate (the one that uses all the bins together),
analogous to how splitting a sample into subsamples and looking at the
scatter of subsample means estimates the standard error of the full mean.

This is empirical, not a closed-form noise propagation: it automatically
captures whatever correlation the demodulation + low-pass filtering
introduces, without having to model it. The trade-off is that each group
has less on-pulse signal than the full window, so `n_groups` shouldn't be
pushed so high that individual groups have too little signal to track
phase at all (watch the per-group SNR if results look like pure noise).

Arguments: same as `coherent_fold`, plus
  n_groups – number of independent longitude sub-ranges, default 4

Returns: the full-bin `coherent_fold` result, plus
  p3_per_pulse_err – Vector{Float64}, 1σ uncertainty on `p3_per_pulse`
  phase_err        – Vector{Float64}, 1σ uncertainty on `phase` [rad]
                      (sub-range phase tracks are de-meaned first, since
                      each has its own arbitrary unwrap integration
                      constant that carries no physical information)
"""
function coherent_fold_jackknife(data::AbstractMatrix, p3::Real, bin_st::Int, bin_end::Int;
                                  ybins::Int=10, lowpass_cutoff::Real=1/200, filter_order::Int=4,
                                  p3_window::Int=60, n_groups::Int=4, auto::Bool=false,
                                  auto_strength::Real=1.0)
    main = coherent_fold(data, p3, bin_st, bin_end; ybins=ybins, lowpass_cutoff=lowpass_cutoff,
                          filter_order=filter_order, p3_window=p3_window, auto=auto,
                          auto_strength=auto_strength)

    N = size(data, 1)
    edges = round.(Int, range(bin_st, bin_end + 1, length=n_groups + 1))
    group_p3 = fill(NaN, n_groups, N)
    group_phase = fill(NaN, n_groups, N)
    for g in 1:n_groups
        st, en = edges[g], edges[g+1] - 1
        en < st && continue
        r = coherent_fold(data, p3, st, en; ybins=ybins, lowpass_cutoff=lowpass_cutoff,
                           filter_order=filter_order, p3_window=p3_window, warn_weak=false,
                           auto=auto, auto_strength=auto_strength)
        group_p3[g, :] = r.p3_per_pulse
        group_phase[g, :] = r.phase .- mean(r.phase)
    end

    p3_per_pulse_err = [std(@view group_p3[:, i]) / sqrt(n_groups) for i in 1:N]
    phase_err = [std(@view group_phase[:, i]) / sqrt(n_groups) for i in 1:N]

    return (folded=main.folded, phase=main.phase, bin=main.bin, p3_per_pulse=main.p3_per_pulse,
            snr=main.snr, p3_per_pulse_err=p3_per_pulse_err, phase_err=phase_err)
end

"""
    scan_lowpass_cutoff(data, p3, bin_st, bin_end; cutoffs, filter_order, ybins,
                         p3_window, n_groups) -> Vector{NamedTuple}

Diagnostic scan to help pick `lowpass_cutoff` for `coherent_fold` /
`p3fold_coherent` per pulsar, instead of guessing a single fixed value
(the cutoff trades off responsiveness to real P3 wobble against noise
suppression, and the right balance is pulsar-specific).

For each candidate cutoff, splits the on-pulse window into `n_groups`
independent longitude sub-ranges (same construction as
`coherent_fold_jackknife`) and runs `coherent_fold` on each separately.
Reports, per cutoff:

  consistency – mean pairwise Pearson correlation between the `n_groups`
                independent `p3_per_pulse` tracks. Real P3 wobble is a
                single underlying physical signal shared by every
                longitude sub-range, so independent tracks should agree;
                pure noise passed through the same filter does not.
                High consistency ⇒ this cutoff is letting real signal
                through, not just noise-shaped smoothness.
  p3_std      – std of the full-bin `p3_per_pulse` track. Rises with
                looser (higher) cutoffs as more noise leaks through;
                falls with tighter (lower) cutoffs as everything gets
                smoothed out (including real wobble at low cutoff).
                Conflates real wobble and noise, so on its own it can't
                distinguish "stable P3" from "over-smoothed noise".
  signal_std  – p3_std * sqrt(max(consistency, 0)): an estimate of how much
                of p3_std is reproducible *signal* rather than independent
                per-subband noise (derivation in the source). This is what
                `p3fold_cutoff_scan` maximises to pick a cutoff — neither
                consistency alone (can be high while p3_std ≈ 0, i.e. a
                cutoff so tight nothing survives to agree on) nor p3_std
                alone (grows with noise) answers "how much real P3 wobble
                does this cutoff actually resolve".
  snr         – matched-filter detection SNR (from `coherent_fold`;
                independent of the cutoff, included as a per-pulsar
                sanity reference).
  p3_per_pulse – the full-bin instantaneous-P3 track for this cutoff
                (same quantity `coherent_fold` returns), so
                `Plot.p3_tracks_vs_cutoff` can overlay how the tracked
                P3 wobble itself changes shape as the cutoff loosens.

No cutoff is chosen inside this function — `p3fold_cutoff_scan` picks
argmax(signal_std) as a suggestion, but always inspect the plots
(`Plot.lowpass_cutoff_scan`'s three panels) before trusting it.

Arguments: same as `coherent_fold`, plus
  cutoffs  – candidate lowpass_cutoff values [cycles/pulse] to try,
             default a log-spaced sweep from 1/1000 to 1/30
"""
function scan_lowpass_cutoff(data::AbstractMatrix, p3::Real, bin_st::Int, bin_end::Int;
                              cutoffs::AbstractVector{<:Real}=[1/1000, 1/500, 1/300, 1/200,
                                                                1/150, 1/100, 1/60, 1/30],
                              filter_order::Int=4, ybins::Int=10, p3_window::Int=60,
                              n_groups::Int=4)
    N = size(data, 1)
    edges = round.(Int, range(bin_st, bin_end + 1, length=n_groups + 1))

    results = NamedTuple[]
    for co in cutoffs
        main = coherent_fold(data, p3, bin_st, bin_end; ybins=ybins, lowpass_cutoff=co,
                              filter_order=filter_order, p3_window=p3_window, warn_weak=false)

        group_p3 = fill(NaN, n_groups, N)
        for g in 1:n_groups
            st, en = edges[g], edges[g+1] - 1
            en < st && continue
            r = coherent_fold(data, p3, st, en; ybins=ybins, lowpass_cutoff=co,
                               filter_order=filter_order, p3_window=p3_window, warn_weak=false)
            group_p3[g, :] = r.p3_per_pulse
        end

        cors = Float64[]
        for a in 1:n_groups, b in a+1:n_groups
            va, vb = @view(group_p3[a, :]), @view(group_p3[b, :])
            if all(isfinite, va) && all(isfinite, vb) && std(va) > 0 && std(vb) > 0
                push!(cors, cor(va, vb))
            end
        end
        consistency = isempty(cors) ? NaN : mean(cors)
        p3_std_val  = std(main.p3_per_pulse)
        # effective *signal* variability: `consistency` estimates what
        # fraction of the observed p3_std is reproducible (shared between
        # independent subbands) rather than independent per-subband noise —
        # for x = s + n1, y = s + n2 with s, n1, n2 independent, cor(x,y) ≈
        # var(s)/var(x) when n1, n2 have variance comparable to x's noise, so
        # std(s) ≈ p3_std * sqrt(consistency). This is what
        # `Plot.lowpass_cutoff_scan`'s third panel maximises to pick a
        # cutoff, instead of consistency alone (high but on ~0 p3_std) or
        # p3_std alone (grows with noise, not with real signal).
        signal_std = (isnan(consistency) || consistency <= 0) ? 0.0 :
                     p3_std_val * sqrt(consistency)

        push!(results, (cutoff=co, consistency=consistency,
                         p3_std=p3_std_val, signal_std=signal_std, snr=main.snr,
                         p3_per_pulse=main.p3_per_pulse))
    end
    return results
end

end # module P3FoldViterbi
