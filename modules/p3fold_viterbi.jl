module P3FoldViterbi

using Statistics
using FFTW
using DSP
using Random

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
  p3_window – smoothing window [pulses] for `p3_per_pulse`, default 20

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
                        p3_window::Int=20)
    N = size(data, 1)
    on = bin_st:bin_end
    f3 = 1.0 / p3

    # 1. single, full-length (non-segmented) complex spatial template at f3
    F = fft(data, 1)
    k = clamp(round(Int, N / p3), 1, N ÷ 2)
    L = F[k+1, :]
    L_on = L[on]

    off = vcat(1:bin_st-1, bin_end+1:size(data, 2))
    sigma_off = isempty(off) ? 0.0 : std([real.(L[off]); imag.(L[off])])
    snr = sigma_off == 0 ? Inf : sqrt(sum(abs2, L_on)) / (sigma_off * sqrt(length(on)))

    # 2. spatial matched-filter projection: one high-SNR complex number per pulse.
    # The static (non-modulated) average profile must be removed first — it lives at
    # frequency 0 and, unlike in `drift_test` (which reads a single FFT bin and never
    # mixes frequencies), a time-domain projection like this carries every frequency
    # through, so the huge DC term would otherwise swamp the low-pass filter below.
    on_data = data[:, on]
    on_data_demeaned = on_data .- mean(on_data, dims=1)
    w = conj.(L_on)
    z = on_data_demeaned * w

    # 3. coherent demodulation at f3, then low-pass filter (sliding, no block edges)
    n = 1:N
    carrier = exp.((-1im * 2π * f3) .* n)
    baseband = z .* carrier
    respf = digitalfilter(Lowpass(lowpass_cutoff), Butterworth(filter_order); fs=1.0)
    baseband_smooth = filtfilt(respf, real.(baseband)) .+ im .* filtfilt(respf, imag.(baseband))

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
                                  p3_window::Int=20, n_groups::Int=4)
    main = coherent_fold(data, p3, bin_st, bin_end; ybins=ybins, lowpass_cutoff=lowpass_cutoff,
                          filter_order=filter_order, p3_window=p3_window)

    N = size(data, 1)
    edges = round.(Int, range(bin_st, bin_end + 1, length=n_groups + 1))
    group_p3 = fill(NaN, n_groups, N)
    group_phase = fill(NaN, n_groups, N)
    for g in 1:n_groups
        st, en = edges[g], edges[g+1] - 1
        en < st && continue
        r = coherent_fold(data, p3, st, en; ybins=ybins, lowpass_cutoff=lowpass_cutoff,
                           filter_order=filter_order, p3_window=p3_window)
        group_p3[g, :] = r.p3_per_pulse
        group_phase[g, :] = r.phase .- mean(r.phase)
    end

    p3_per_pulse_err = [std(@view group_p3[:, i]) / sqrt(n_groups) for i in 1:N]
    phase_err = [std(@view group_phase[:, i]) / sqrt(n_groups) for i in 1:N]

    return (folded=main.folded, phase=main.phase, bin=main.bin, p3_per_pulse=main.p3_per_pulse,
            snr=main.snr, p3_per_pulse_err=p3_per_pulse_err, phase_err=phase_err)
end


# ---------------------------------------------------------------------------
# Agent variants: automatic low-pass cutoff and a variable-P3 track for
# `coherent_fold` (tests and rationale: docs/coherent_fold_params.md).
# ---------------------------------------------------------------------------

"Spatial matched filter at f3 = 1/p3: conj of the full-length FFT bin k = round(N/p3) (as in `coherent_fold`)."
function _template_weights(X::AbstractMatrix, p3::Real)
    N = size(X, 1)
    F = fft(X, 1)
    k = clamp(round(Int, N / p3), 1, N ÷ 2)
    return conj.(F[k+1, :])
end

function _lowpass(x::AbstractVector, cutoff::Real, order::Int)
    respf = digitalfilter(Lowpass(cutoff), Butterworth(order); fs=1.0)
    return filtfilt(respf, real.(x)) .+ im .* filtfilt(respf, imag.(x))
end

function _interp_lin(x, y, x0)
    j = searchsortedlast(x, x0)
    x[j] == x0 && return y[j]
    return y[j] + (y[j+1] - y[j]) * (x0 - x[j]) / (x[j+1] - x[j])
end

"""
    _carrier_track(X, p3, cutoff; filter_order, niter, threshold) -> (s, carrier)

Demodulate the matched-filter series z = X·w at the carrier phase and low-pass
it: s = LP(z·e^{−i·carrier}); total P3-phase = carrier + arg(s). With
`niter = 0` the carrier is the constant-P3 ramp 2πn/p3 (`coherent_fold`).
Each further pass moves the carrier onto the phase found so far (unwrapped
over pulses with |s| ≥ `threshold`, linear across the rest), so a local P3
far from p3 sits near zero frequency and is not attenuated by the filter.
X must have its column means removed.
"""
function _carrier_track(X::AbstractMatrix, p3::Real, cutoff::Real; filter_order::Int=6,
                        niter::Int=0, threshold=nothing)
    N = size(X, 1)
    n = collect(1:N)
    z = X * _template_weights(X, p3)
    carrier = (2π / p3) .* n
    s = _lowpass(z .* exp.(-1im .* carrier), cutoff, filter_order)
    for _ in 1:niter
        keep = threshold === nothing ? trues(N) : abs.(s) .>= threshold
        idx = findall(keep)
        length(idx) < 3 && break
        θ = DSP.unwrap(angle.(s[idx]))
        carrier = carrier .+ [i ≤ idx[1] ? θ[1] : i ≥ idx[end] ? θ[end] : _interp_lin(idx, θ, i) for i in n]
        s = _lowpass(z .* exp.(-1im .* carrier), cutoff, filter_order)
    end
    return s, carrier
end

"Quantile `q` of |s| for pulse-order shuffles (no P3 periodicity left) — the noise level of the demodulated amplitude."
function _shuffle_level(X, p3, cutoff; filter_order=6, niter=0, threshold=nothing, q=0.5, nshuffle=5, seed=1)
    rng = MersenneTwister(seed)
    N = size(X, 1)
    vals = Float64[]
    for _ in 1:nshuffle
        s, _ = _carrier_track(X[randperm(rng, N), :], p3, cutoff; filter_order=filter_order,
                              niter=niter, threshold=threshold)
        append!(vals, abs.(s))
    end
    return quantile(vals, q)
end

"Fraction of fluctuation variance of Y explained by folding with `phase` into `ybins`, minus the noise bias (ybins−1)/(N−1)."
function _fold_r2(Y::AbstractMatrix, phase::AbstractVector, ybins::Int)
    N = size(Y, 1)
    b = [Int(floor(mod(phase[i] / (2π) * ybins, ybins))) + 1 for i in 1:N]
    between = 0.0
    for y in 1:ybins
        idx = findall(==(y), b)
        isempty(idx) && continue
        between += length(idx) * sum(abs2, mean(Y[idx, :], dims=1))
    end
    return between / sum(abs2, Y) - (ybins - 1) / (N - 1)
end

"Interleaved blocks of 4 on-pulse bins: two halves with independent noise, both covering the whole profile."
_bin_halves(nb::Int) = (A = [j for j in 1:nb if iseven((j - 1) ÷ 4)]; (A, setdiff(1:nb, A)))

function _cv_r2(X, A, B, p3, cutoff, ybins, filter_order)
    sc = 0.0
    for (P, Q) in ((A, B), (B, A))
        s, carrier = _carrier_track(X[:, P], p3, cutoff; filter_order=filter_order)
        sc += _fold_r2(X[:, Q], carrier .+ angle.(s), ybins) / 2
    end
    return sc
end


"""
    auto_cutoff_agent(data, p3, bin_st, bin_end; ybins=10, grid=(1/16, 1/10, 1/8, 1/6, 1/4, 1/3),
                      filter_order=6, nshuffle=5, seed=2) -> NamedTuple

Choose `lowpass_cutoff` for `coherent_fold` from the data. Candidates are
f_c = g·f3 for g in `grid` (f3 = 1/p3): the cutoff must exceed the P3 wander
|1/P3_local − f3| to follow it, and stay ≲ f3/3 — above that the −f3 image and
the intensity modulation (nulls) leak into the phase.

Score: cross-validated fold quality. The phase is estimated from one half of
the on-pulse bins (interleaved blocks of 4), the other half is folded with
it, and vice versa; the score is the fraction of fluctuation variance
explained by the fold. Its mean over `nshuffle` pulse-order shuffles is
subtracted: a fast cutoff lets the phase follow each pulse's own subpulse
position, which "explains" variance in the other half too, but has nothing to
do with the P3 periodicity (shuffling keeps it, destroys the periodicity).

Returns: cutoff [cycles/P], grid (cutoffs tried), score (ΔR² per cutoff).
On 10 pulsars the optimum fell at 0.07–0.33·f3; the fixed 1/300 of
`p3fold_coherent` gave 2–6× lower ΔR² for P3 ≲ 15 (synthetic benchmark:
fold correlation with the truth 0.68 → 0.85).
"""
function auto_cutoff_agent(data::AbstractMatrix, p3::Real, bin_st::Int, bin_end::Int;
                           ybins::Int=10, grid=(1/16, 1/10, 1/8, 1/6, 1/4, 1/3),
                           filter_order::Int=6, nshuffle::Int=5, seed::Int=2)
    X = data[:, bin_st:bin_end]
    X = X .- mean(X, dims=1)
    N = size(X, 1)
    A, B = _bin_halves(size(X, 2))
    cutoffs = [g / p3 for g in grid]
    rng = MersenneTwister(seed)
    perms = [randperm(rng, N) for _ in 1:nshuffle]
    score = map(cutoffs) do fc
        _cv_r2(X, A, B, p3, fc, ybins, filter_order) -
            mean(_cv_r2(X[pp, :], A, B, p3, fc, ybins, filter_order) for pp in perms)
    end
    return (cutoff=cutoffs[argmax(score)], grid=cutoffs, score=score)
end


"""
Local P3 from the slope of the total phase (weighted by |s|²) over ±window/2
pulses. Computed within each continuous run of `keep` pulses only: the phase
is unwrapped and the slope fitted without crossing a gap, because the number
of P3 cycles inside a gap is unknown (a gap of ~P3 pulses makes it
ambiguous, J1750-3503) — across a gap P3(n) may jump, it is not interpolated.
NaN outside `keep`, where a run has too few pulses in the window, and where
the result falls outside `bounds` (a slope ≈ 0 or < 0 in a short stretch gave
P3 ≈ −700, J1919+0134).
"""
function _weighted_p3(s, carrier, keep, window; bounds=(2.0, Inf))
    N = length(s)
    out = fill(NaN, N)
    wt = abs2.(s)
    h = max(2, window ÷ 2)
    i = 1
    while i ≤ N
        if !keep[i]
            i += 1
            continue
        end
        j = i
        while j < N && keep[j+1]
            j += 1
        end
        run = i:j
        ph = carrier[run] .+ DSP.unwrap(angle.(s[run]))
        for (m, t) in enumerate(run)
            lo = max(1, m - h); hi = min(length(run), m + h)
            hi - lo + 1 < max(5, h ÷ 2) && continue
            x = Float64.(run[lo:hi]); y = ph[lo:hi]; w = wt[run[lo:hi]]
            xm = sum(w .* x) / sum(w); ym = sum(w .* y) / sum(w)
            v = 2π / (sum(w .* (x .- xm) .* (y .- ym)) / sum(w .* (x .- xm) .^ 2))
            out[t] = bounds[1] ≤ v ≤ bounds[2] ? v : NaN
        end
        i = j + 1
    end
    return out
end

"""
    _energy_nulls(data, on; smooth=5, minlen=2) -> BitVector

Nulls from pulse energy E(n) = Σ_on I(n,φ), averaged over `smooth` pulses.
Null fraction nf = 2·frac(Ē < 0) (Ritchings 1976: noise is symmetric, so nulls
put as many values below zero as above); if nf < 5% no pulse is flagged,
otherwise the nf lowest Ē, kept only in episodes of ≥ `minlen` pulses.
Averaging first matters for weak pulsars: single pulses of J1750-3503
(energy S/N 1.6, no nulls) fall below zero often enough to fake an 18% null
fraction; nulls last several pulses, noise dips of emitting pulses do not, so
the mean over 5 pulses raises the emission S/N √5× and leaves nulls at zero.
Independent of the low-pass filter, so it catches nulls shorter than the
filter memory, which the |s| threshold misses.
"""
function _energy_nulls(data::AbstractMatrix, on; smooth::Int=5, minlen::Int=2)
    N = size(data, 1)
    E = vec(sum(data[:, on], dims=2))
    h = smooth ÷ 2
    Es = [mean(E[max(1, i - h):min(N, i + h)]) for i in 1:N]
    nf = 2 * mean(Es .< 0)
    out = falses(N)
    nf < 0.05 && return out
    m = Es .< quantile(Es, min(nf, 0.9))
    i = 1
    while i ≤ N
        if m[i]
            j = i
            while j < N && m[j+1]
                j += 1
            end
            j - i + 1 ≥ minlen && (out[i:j] .= true)
            i = j + 1
        else
            i += 1
        end
    end
    return out
end

"Variant C on one on-pulse range: adaptive carrier, shuffle-quantile threshold, edges and `exclude` (nulls) dropped."
function _variant_c(X, p3, cutoff; filter_order=6, niter=2, nshuffle=5, threshold_q=0.5, exclude=nothing)
    N = size(X, 1)
    thr0 = _shuffle_level(X, p3, cutoff; filter_order=filter_order, nshuffle=nshuffle, q=threshold_q)
    s, carrier = _carrier_track(X, p3, cutoff; filter_order=filter_order, niter=niter, threshold=thr0)
    thr = niter == 0 ? thr0 : _shuffle_level(X, p3, cutoff; filter_order=filter_order, niter=niter,
                                             threshold=thr0, nshuffle=nshuffle, q=threshold_q)
    # filtfilt edges: 1/(2 f_c); regression window at least 3·P3 and 30 pulses — at
    # f_c ≈ f3/3 the 1/(2 f_c) ≈ 1.5·P3 window (5–10 P for P3 ≈ 4) made P3(n) jitter
    # pulse to pulse; the 30-P floor (≈ 8·P3 at P3 ≈ 4) halves the P3(n) error for
    # P3 < 8 on synthetic data and changes nothing for P3 ≳ 10
    nedge = max(4, round(Int, 1 / (2cutoff)))
    window = max(nedge, round(Int, 3p3), 30)
    edge = falses(N)
    edge[1:min(N, nedge)] .= true
    edge[max(1, N - nedge + 1):N] .= true
    keep = (abs.(s) .>= thr) .& .!edge
    exclude === nothing || (keep .&= .!exclude)
    return (s=s, carrier=carrier, keep=keep, threshold=thr, window=window,
            p3=_weighted_p3(s, carrier, keep, window; bounds=(2.0, 3p3)))
end


"""
    coherent_fold_agent(data, p3, bin_st, bin_end; ybins=10, lowpass_cutoff=:auto, filter_order=6,
                        niter=2, nshuffle=5, n_groups=4, threshold_q=0.5, split_nulls=true) -> NamedTuple

`coherent_fold` with (1) the low-pass cutoff chosen by `auto_cutoff_agent`
(or a number), and (2) the P3(n) track of "variant C":

  1. matched-filter series z(n) and demodulation at f3, low-pass f_c;
  2. `niter` further passes with the carrier following the phase found so
     far (adaptive carrier): a local P3 far from p3 is no longer attenuated
     by the filter and |s| no longer drops to zero there;
  3. pulses with |s| below the `threshold_q` quantile (default: median) of
     |s| for pulse-order shuffles — the noise level, no periodicity left (no
     usable phase: nulls, weak or ambiguous stretches) — and the first/last
     1/(2 f_c) pulses (`filtfilt` edges) are left out of P3(n) — not out of
     the fold, whose bin assignment is modulo 2π; a lower `threshold_q`
     shortens the gaps at the cost of noisier phase;
     with `split_nulls` also the nulls found from pulse energy
     (`_energy_nulls`) — they can be shorter than the filter memory and
     then never reach the |s| threshold, while |s| dips inside them flip the
     phase (single P3(n) spikes → ∞ on synthetic data);
  4. P3(n) from the slope of the total phase, weighted by |s|², over
     max(1/(2 f_c), 3·P3, 30) pulses (values outside [2, 3·p3] → NaN), within each continuous run of used pulses (no
     unwrapping or fitting across a gap: the number of cycles in a gap is
     unknown, so P3(n) may jump there).

Same on-pulse window and template as `coherent_fold`. Errors: P3(n) recomputed
on `n_groups` contiguous longitude sub-ranges (independent noise), std/√n
over the groups that give a value.

Synthetic benchmark (36 cases, P3 6–45, wander/step, 30% nulls, S/N 1.5–5):
relative P3(n) error 5.6% (median) at 94% coverage; fold correlation with
the truth 0.85 (fixed 1/300: 0.68).

Returns:
  folded         – ybins × N_bins fold with the variant-C phase: MEAN of the pulses
                   in each phase bin (NaN for an empty bin), not the sum
  counts         – number of pulses in each phase bin
  phase, bin     – total P3-phase [rad] and fold bin per pulse
  p3_per_pulse   – P3(n) [P], NaN where not measured
  p3_per_pulse_err – 1σ from the longitude groups (NaN if < 2 groups)
  used           – pulses that enter P3(n)
  nulls          – pulses flagged as nulls from energy (all false if `split_nulls=false`)
  amplitude, threshold – |s(n)| and the shuffle level
  lowpass_cutoff – cutoff used; cutoff_grid, cutoff_score (if :auto)
  snr            – matched-filter detection significance (as `coherent_fold`)
"""
function coherent_fold_agent(data::AbstractMatrix, p3::Real, bin_st::Int, bin_end::Int;
                             ybins::Int=10, lowpass_cutoff=:auto, filter_order::Int=6,
                             niter::Int=2, nshuffle::Int=5, n_groups::Int=4, threshold_q::Real=0.5,
                             split_nulls::Bool=true)
    N = size(data, 1)
    on = bin_st:bin_end
    if lowpass_cutoff === :auto
        ac = auto_cutoff_agent(data, p3, bin_st, bin_end; ybins=ybins, filter_order=filter_order,
                               nshuffle=nshuffle)
        fc, grid, score = ac.cutoff, ac.grid, ac.score
    else
        fc, grid, score = Float64(lowpass_cutoff), nothing, nothing
    end

    X = data[:, on]
    X = X .- mean(X, dims=1)
    nulls = split_nulls ? _energy_nulls(data, on) : falses(N)
    c = _variant_c(X, p3, fc; filter_order=filter_order, niter=niter, nshuffle=nshuffle, threshold_q=threshold_q,
                   exclude=nulls)
    phase = c.carrier .+ angle.(c.s)
    bin = [Int(floor(mod(phase[i] / (2π) * ybins, ybins))) + 1 for i in 1:N]
    # mean per phase bin, not the sum of `build_template`: the data-driven phases fill
    # the bins unevenly (P3 ≈ 2: phases cluster at two values, J1539-4828: 38–155 pulses
    # per bin at ybins = 16), and a sum shows the bin occupancy instead of the emission
    counts = [count(==(y), bin) for y in 1:ybins]
    folded = build_template(data, bin, ybins) ./ counts
    folded[counts .== 0, :] .= NaN

    # detection significance, as in coherent_fold
    F = fft(data, 1)
    k = clamp(round(Int, N / p3), 1, N ÷ 2)
    L = F[k+1, :]
    off = vcat(1:bin_st-1, bin_end+1:size(data, 2))
    sigma_off = isempty(off) ? 0.0 : std([real.(L[off]); imag.(L[off])])
    snr = sigma_off == 0 ? Inf : sqrt(sum(abs2, L[on])) / (sigma_off * sqrt(length(on)))

    # errors from independent longitude sub-ranges
    edges = round.(Int, range(bin_st, bin_end + 1, length=n_groups + 1))
    gp3 = fill(NaN, n_groups, N)
    for g in 1:n_groups
        st, en = edges[g], edges[g+1] - 1
        en - st < 2 && continue
        Xg = data[:, st:en]
        Xg = Xg .- mean(Xg, dims=1)
        gp3[g, :] = _variant_c(Xg, p3, fc; filter_order=filter_order, niter=niter, nshuffle=nshuffle,
                               threshold_q=threshold_q, exclude=nulls).p3
    end
    p3_err = map(1:N) do i
        v = filter(isfinite, @view gp3[:, i])
        length(v) < 2 ? NaN : std(v) / sqrt(length(v))
    end

    return (folded=folded, counts=counts, phase=phase, bin=bin, p3_per_pulse=c.p3, p3_per_pulse_err=p3_err,
            used=c.keep, nulls=nulls, amplitude=abs.(c.s), threshold=c.threshold, lowpass_cutoff=fc,
            cutoff_grid=grid, cutoff_score=score, snr=snr)
end

end # module P3FoldViterbi
