"""
Travel — does the subpulse pattern *move* in longitude, or does it only
brighten and fade in place?

This is the same question `PhaseDrift` (modules/phase_modulation.jl) asks, but
answered without ever going near P3. The whole method rests on one algebraic
identity.

Write the fluctuations (static profile removed) as δI(n, φ) and form the
two-dimensional autocorrelation

    K(Δ, τ) = Σ_{n,φ} δI(n,φ) · δI(n+τ, φ+Δ)

Pure amplitude modulation is *separable*: every longitude follows one and the
same temporal waveform, scaled by a real (possibly negative) factor,
δI(n,φ) = a(φ)·w(n). Then

    K(Δ, τ) = [Σ_φ a(φ)a(φ+Δ)] · [Σ_n w(n)w(n+τ)]

and the second factor is *exactly* even in τ — substituting m = n−τ turns
Σ_n w(n)w(n−τ) into Σ_m w(m+τ)w(m), the same sum over the same pairs. No
stationarity, no periodicity, no P3 is assumed anywhere: w(n) may wobble,
change period, null, or stop being periodic altogether. So

    A(Δ, τ) = K(Δ, τ) − K(Δ, −τ)  ≡  0     under amplitude modulation

identically, for the realised data and not merely in expectation. A travelling
pattern δI = f(φ − v·n) gives K(Δ,τ) = C_f(Δ − vτ), which is not even in τ, and
A picks it up.

Because K(Δ,τ) = K(−Δ,−τ) (relabel the summation indices), A is antisymmetric
in both arguments and A(Δ,0) = A(0,τ) = 0, so all the information sits in the
quadrant Δ ≥ 1, τ ≥ 1 — that is what `antisym_map` returns.

Written without the separability assumption,

    A(Δ, τ) = Σ_φ [ C_{φ,φ+Δ}(τ) − C_{φ+Δ,φ}(τ) ]

i.e. A is the net answer to "does longitude φ lead φ+Δ, or the other way
round?". It vanishes whenever no longitude systematically leads another, which
covers more than the textbook P3-only case: two components in antiphase (a(φ)
changing sign) stay separable and give exactly zero, and two components carrying
*independent* modulations of different P3 give zero in expectation. Zapped
pulses are harmless for the same reason — blanking a pulse multiplies every
longitude by the same factor, so it stays separable.

Relation to the 2DFS. The Fourier transform of K is the 2D power spectrum, and
the part of A that survives is the asymmetry of that spectrum about 1/P2 = 0 —
formally the same channel the classic 2DFS criterion uses. The difference is in
the estimator: the usual procedure takes the *centroid* of the power inside a
hand-drawn box, so the symmetric ridge that stochastic pulse-shape variability
piles up along 1/P2 = 0 enters the estimate and biases it. Here the symmetric
part is projected out algebraically before anything is estimated, there is no
box and no centroid, and the whole (Δ, τ) plane is used instead of one
frequency bin.

What is *not* exactly zero: the data carry noise, so the measured A scatters
around zero and the spread still has to be calibrated (`travel_test` does this
with surrogates). The point is that under H0 the signal contributes exactly
zero — not approximately — so only the noise terms need modelling, and the
systematic that plagues the power centroid has no counterpart here. The one
real vulnerability is a noise field that is *not* time-reversal symmetric
(gain drift, RFI ramping through the observation, a badly subtracted baseline);
`travel_test` therefore runs the identical statistic on an off-pulse strip,
where the answer must be consistent with zero.
"""
module Travel

using FFTW
using Statistics
using Random
using LinearAlgebra


"""
    highpass(X, halfwin) -> Matrix{Float64}

Per-longitude high-pass: subtract from every time series X[:, j] its centred
running mean over 2·halfwin+1 pulses (edge-truncated). Removes the static
profile together with slow intensity drifts (scintillation, gain wander), which
is what protects the test against the one systematic it is vulnerable to.

The filter is the *same* at every longitude, so it maps a separable field onto
a separable field: a(φ)w(n) → a(φ)·[w − runmean(w)](n). The exact-zero property
of `antisym_map` therefore survives the preprocessing untouched.

`halfwin <= 0` falls back to plain mean subtraction (static profile only).
Mirrors `PhaseDrift.subtract_running_mean`, duplicated so this module can be
included on its own.
"""
function highpass(X::AbstractMatrix, halfwin::Int)
    halfwin <= 0 && return X .- mean(X, dims=1)
    N = size(X, 1)
    C = cumsum(X, dims=1)
    Y = Matrix{Float64}(undef, size(X))
    for n in 1:N
        lo, hi = max(1, n - halfwin), min(N, n + halfwin)
        if lo == 1
            Y[n, :] = X[n, :] .- (@view C[hi, :]) ./ hi
        else
            Y[n, :] = X[n, :] .- ((@view C[hi, :]) .- (@view C[lo-1, :])) ./ (hi - lo + 1)
        end
    end
    return Y
end


"""
    corr_map(X, max_lag, max_dphi) -> Matrix{Float64}

Two-dimensional autocorrelation K(Δ, τ) = Σ_{n,φ} X[n,φ]·X[n+τ, φ+Δ] of the
fluctuation matrix X (N_pulses × N_bins), via FFT with zero padding so the
result is the *linear* (not circular) correlation — pulses at the start of the
observation are never paired with pulses at the end, which matters because the
statistic is precisely a test for time asymmetry.

Returned as a (max_lag+1) × (2·max_dphi+1) matrix: row i is τ = i−1, column j is
Δ = j − max_dphi − 1. Negative τ is not returned because K(Δ,−τ) = K(−Δ,τ).

Cost is one real 2D FFT pair, independent of the requested lag ranges.
"""
function corr_map(X::AbstractMatrix, max_lag::Int, max_dphi::Int)
    N, M = size(X)
    max_lag  >= 1 || error("max_lag must be ≥ 1 (got $max_lag)")
    max_dphi >= 1 || error("max_dphi must be ≥ 1 (got $max_dphi)")
    max_lag  < N || error("max_lag ($max_lag) must be < N_pulses ($N)")
    max_dphi < M || error("max_dphi ($max_dphi) must be < N_bins ($M)")

    Np = nextprod((2, 3, 5), N + max_lag)
    Mp = nextprod((2, 3, 5), M + max_dphi)
    Xp = zeros(Float64, Np, Mp)
    Xp[1:N, 1:M] .= X
    Kf = irfft(abs2.(rfft(Xp)), Np)

    K = Matrix{Float64}(undef, max_lag + 1, 2 * max_dphi + 1)
    for (i, tau) in enumerate(0:max_lag), (j, dphi) in enumerate(-max_dphi:max_dphi)
        K[i, j] = Kf[tau + 1, mod(dphi, Mp) + 1]
    end
    return K
end


"""
    corr_map_direct(X, max_lag, max_dphi) -> Matrix{Float64}

Reference implementation of `corr_map` as an explicit quadruple sum. O(N·M·
max_lag·max_dphi) and only for validating the FFT path on small inputs — see
`Travel.selftest`.
"""
function corr_map_direct(X::AbstractMatrix, max_lag::Int, max_dphi::Int)
    N, M = size(X)
    K = zeros(max_lag + 1, 2 * max_dphi + 1)
    for (i, tau) in enumerate(0:max_lag), (j, dphi) in enumerate(-max_dphi:max_dphi)
        s = 0.0
        for n in 1:N-tau, phi in max(1, 1-dphi):min(M, M-dphi)
            s += X[n, phi] * X[n+tau, phi+dphi]
        end
        K[i, j] = s
    end
    return K
end


"""
    antisym_map(K) -> Matrix{Float64}

The travel statistic map A(Δ, τ) = K(Δ, τ) − K(−Δ, τ), for τ = 1…max_lag and
Δ = 1…max_dphi, from the output of `corr_map`. (K(−Δ,τ) = K(Δ,−τ) by the
central symmetry of K, so this is the same thing as K(Δ,τ) − K(Δ,−τ).)

Rows are τ = 1, 2, …; columns are Δ = 1, 2, …. Everything outside this quadrant
is either zero by construction (τ = 0 or Δ = 0) or a sign flip of what is here.
"""
function antisym_map(K::AbstractMatrix)
    max_lag  = size(K, 1) - 1
    max_dphi = (size(K, 2) - 1) ÷ 2
    ctr = max_dphi + 1
    A = Matrix{Float64}(undef, max_lag, max_dphi)
    for i in 1:max_lag, d in 1:max_dphi
        A[i, d] = K[i+1, ctr+d] - K[i+1, ctr-d]
    end
    return A
end


"""
    travel_T(X, max_lag, max_dphi) -> Float64

Omnibus travel statistic T = Σ_{Δ,τ} A(Δ,τ)², assumption-free: no drift rate,
no P3 and no template enter it. Zero in the mean under amplitude modulation,
positive whenever some longitude systematically leads another.
"""
travel_T(X::AbstractMatrix, max_lag::Int, max_dphi::Int) =
    sum(abs2, antisym_map(corr_map(X, max_lag, max_dphi)))


"""
    sym_map(K) -> Matrix{Float64}

The Δ-*even* companion of `antisym_map`: E(Δ,τ) = K(Δ,τ) + K(−Δ,τ), on the same
quadrant (τ = 1…max_lag, Δ = 1…max_dphi).

Everything `antisym_map` throws away lands here — including all of a separable
amplitude modulation. On its own E discriminates nothing, which is the point:
it measures how much coherent modulation the pulsar has at a given (P2, P3)
*whether or not it travels*, and so provides the per-pulsar yardstick that turns
the travel projection into a ratio. See `travel_test`'s `R` field.
"""
function sym_map(K::AbstractMatrix)
    max_lag  = size(K, 1) - 1
    max_dphi = (size(K, 2) - 1) ÷ 2
    ctr = max_dphi + 1
    E = Matrix{Float64}(undef, max_lag, max_dphi)
    for i in 1:max_lag, d in 1:max_dphi
        E[i, d] = K[i+1, ctr+d] + K[i+1, ctr-d]
    end
    return E
end


"""
    _travel_maps(X, max_lag, max_dphi, edges) -> (K00, A, E, blocks)

The travel map of the whole stretch together with the maps of the individual
pulse blocks delimited by `edges`, computed in one place so that the observed
data and the surrogates go through exactly the same reduction. Blocks shorter
than `max_lag` are dropped.
"""
function _travel_maps(X::AbstractMatrix, max_lag::Int, max_dphi::Int,
                      edges::Vector{Int})
    K = corr_map(X, max_lag, max_dphi)
    A = antisym_map(K)
    E = sym_map(K)
    blocks = Matrix{Float64}[]
    for b in 1:length(edges)-1
        rows = edges[b]:edges[b+1]-1
        length(rows) <= max_lag && continue
        push!(blocks, antisym_map(corr_map(@view(X[rows, :]), max_lag, max_dphi)))
    end
    return K[1, max_dphi+1], A, E, blocks
end

_inc_power(blocks) = isempty(blocks) ? 0.0 : sum(Ab -> sum(abs2, Ab), blocks)


"""
    drift_template(max_lag, max_dphi, p2, p3; npulses, non) -> Matrix{Float64}

The map a rigidly drifting pattern makes,

    A(Δ, τ) / K(0,0) = (1 − τ/N)(1 − Δ/M) · 2 · sin(2πΔ/P2) · sin(2πτ/P3)

on the same (τ, Δ) grid `antisym_map` returns. `p2` is in longitude bins and
carries the drift sense in its sign (positive = later longitudes light up later,
Szary+2022); `p3` is in pulse periods.

The leading factor is the triangular taper of a *linear* correlation: `corr_map`
sums over (N−τ)(M−|Δ|) pairs at lag (Δ,τ) against NM pairs at the origin, so the
measured map is tapered even for a perfectly coherent pattern. It is pure
geometry, known exactly, and leaving it out biases the matched projection low —
by 23% at N=600, M=40 over the default lag ranges, which is the difference
between "this pulsar barely drifts" and "it drifts as claimed". Pass `npulses`
(N) and `non` (on-pulse width M) to apply it; omit them for the untapered form.

Projecting the measured map onto this is what turns the test around: instead of
asking "is there any travel", it asks "is the travel the one this pulsar is
claimed to have", which is a hypothesis that can be *rejected* rather than
merely left unconfirmed. See `travel_test`'s `frac` fields.
"""
function drift_template(max_lag::Int, max_dphi::Int, p2::Real, p3::Real;
                        npulses::Union{Int,Nothing}=nothing,
                        non::Union{Int,Nothing}=nothing)
    abs(p2) > 0 || error("p2 must be nonzero")
    p3 > 0      || error("p3 must be positive (got $p3)")
    tap_t = isnothing(npulses) ? ones(max_lag)  : [1 - t / npulses for t in 1:max_lag]
    tap_d = isnothing(non)     ? ones(max_dphi) : [1 - d / non      for d in 1:max_dphi]
    return [2 * tap_t[t] * tap_d[d] * sin(2π * d / p2) * sin(2π * t / p3)
            for t in 1:max_lag, d in 1:max_dphi]
end


"""
    drift_template_even(max_lag, max_dphi, p2, p3; npulses, non) -> Matrix{Float64}

The Δ-even half of the same drift model, to be projected onto `sym_map`:

    E(Δ, τ) / K(0,0) = (1 − τ/N)(1 − Δ/M) · 2 · cos(2πΔ/P2) · cos(2πτ/P3)

Writing the drift's correlation out,

    cos(2π(Δ/P2 − τ/P3)) = cos(2πΔ/P2)·cos(2πτ/P3) + sin(2πΔ/P2)·sin(2πτ/P3)
                           └──── even in Δ ────┘     └──── odd in Δ ────┘

shows the two halves carry *equal* coefficients for a rigidly travelling
pattern. Amplitude modulation, being separable, puts everything in the even half
and nothing in the odd one. The ratio of the two projections is therefore 1 for
a pure drift and 0 for pure amplitude modulation — and, unlike either projection
alone, it needs no external yardstick: modulation strength, harmonic content,
loss of coherence with lag and the triangular taper all multiply the two halves
identically and cancel. See `travel_test`'s `R`.
"""
function drift_template_even(max_lag::Int, max_dphi::Int, p2::Real, p3::Real;
                             npulses::Union{Int,Nothing}=nothing,
                             non::Union{Int,Nothing}=nothing)
    abs(p2) > 0 || error("p2 must be nonzero")
    p3 > 0      || error("p3 must be positive (got $p3)")
    tap_t = isnothing(npulses) ? ones(max_lag)  : [1 - t / npulses for t in 1:max_lag]
    tap_d = isnothing(non)     ? ones(max_dphi) : [1 - d / non      for d in 1:max_dphi]
    return [2 * tap_t[t] * tap_d[d] * cos(2π * d / p2) * cos(2π * t / p3)
            for t in 1:max_lag, d in 1:max_dphi]
end


"""
    separable_signal(X, noise_var) -> Matrix{Float64}

Rank-1 (i.e. exactly separable, exactly H0) stand-in for the pulsar's own
modulation, for the surrogates of `travel_test`: the leading SVD mode of X,
rescaled so its total power equals the noise-debiased fluctuation power
Σ X² − N·M·σ².

Only the *amplitude* of the surrogate signal matters, which is why a crude
rank-1 stand-in is enough. Under H0 the signal contributes exactly zero to A
whatever its shape, so it enters the null only through the signal×noise cross
term, whose size is set by the signal power and its temporal autocorrelation —
both of which the rescaled leading mode reproduces.
"""
function separable_signal(X::AbstractMatrix, noise_var::Real)
    F = svd(X)
    X1 = F.S[1] .* (F.U[:, 1] * F.V[:, 1]')
    p_tot = sum(abs2, X) - length(X) * noise_var
    p1 = sum(abs2, X1)
    (p1 > 0 && p_tot > 0) || return zeros(size(X))
    return X1 .* sqrt(p_tot / p1)
end


"""
    offpulse_starts(bin_st, bin_end, nbins, width) -> Vector{Int}

First columns of every contiguous strip of `width` longitude bins that lies
entirely outside the on-pulse window [bin_st, bin_end]. Empty if the off-pulse
region is narrower than the on-pulse one.
"""
function offpulse_starts(bin_st::Int, bin_end::Int, nbins::Int, width::Int)
    s = Int[]
    for a in 1:(bin_st - width)
        push!(s, a)
    end
    for a in (bin_end + 1):(nbins - width + 1)
        push!(s, a)
    end
    return s
end


"""
    offpulse_strips(bin_st, bin_end, nbins, width) -> (strips, wrapped)

Column index sets for the noise surrogates, `width` bins each, all off-pulse.

Normally these are the contiguous strips of `offpulse_starts`, which keep the
correlation between neighbouring longitude bins exactly as the data has it. A
pulsar whose profile is wider than the longest off-pulse run has no such strip
(seen in the TPA sample for on-pulse windows of ~350 bins); rather than refuse
to analyse it, the strips are then taken cyclically around the off-pulse index
list, which costs one splice point per strip and is flagged by `wrapped`.
"""
function offpulse_strips(bin_st::Int, bin_end::Int, nbins::Int, width::Int)
    s = offpulse_starts(bin_st, bin_end, nbins, width)
    isempty(s) || return ([collect(a:a+width-1) for a in s], false)
    off  = vcat(1:bin_st-1, bin_end+1:nbins)
    noff = length(off)
    noff >= width || error(
        "Off-pulse region ($noff bins) is narrower than the on-pulse window " *
        "($width bins): no noise surrogates possible")
    return ([off[[mod1(a + k, noff) for k in 0:width-1]] for a in 1:noff], true)
end


"""
    noise_block!(dest, Xhp, strips, rng)

Draw one noise realisation into `dest`: a contiguous off-pulse strip of `width`
longitude bins, taken from the already high-passed data `Xhp` and given a random
circular shift along the pulse axis.

Using the pulsar's own off-pulse data rather than a noise model keeps the real
noise level, its temporal correlation *and* the correlation between neighbouring
longitude bins (a contiguous strip, not resampled columns) — all three enter the
signal×noise cross term that sets the width of the null. The circular shift
decorrelates realisations drawn from overlapping strips; it splices the end of
the observation onto the start once, which is negligible over ~10³ pulses.
"""
function noise_block!(dest::AbstractMatrix, Xhp::AbstractMatrix,
                      strips::Vector{Vector{Int}}, rng)
    idx = strips[rand(rng, 1:length(strips))]
    sh  = rand(rng, 0:size(Xhp, 1)-1)
    dest .= circshift(@view(Xhp[:, idx]), (sh, 0))
    return dest
end


"""
    ridge(A) -> NamedTuple

Structure of the travel map, from its leading singular mode.

For a pattern drifting as cos(2π(φ/P2 − n/P3)) the identity
cos(x−y) − cos(x+y) = 2 sin x sin y gives

    A(Δ, τ) = 2 · sin(2πΔ/P2) · sin(2πτ/P3)

— exactly rank 1, a product of one sinusoid in longitude lag and one in pulse
lag. So the leading SVD mode of A is the natural summary: `rank1_frac` says how
well that separable-oscillation form fits (near 1 for a clean coherent drift,
small for noise or for a map built out of unrelated structure), and the two
singular vectors carry the periods.

  * the Δ-profile crosses zero at Δ = P2/2 → `p2` (in bins). No crossing inside
    the searched range means only a lower limit, flagged by `p2_lower_limit` —
    which is the interesting case, since a P2 much wider than the visible
    profile is exactly the drift a 2DFS centroid cannot see.
  * the τ-profile crosses zero at τ = P3/2 → `p3`, an estimate the method
    produces as a by-product rather than requires as an input. NaN when the
    searched lag range is shorter than about P3/2 — raise `max_lag` for a
    long-P3 pulsar if this number is wanted.

Both periods are *indicative*. The factor-of-two relations are exact only for
the single-sinusoid case above; real profiles carry harmonics and the pattern
loses coherence with lag, and both shift the first crossing. Measured against
catalogue P3 the recovered values come out within a few to twenty per cent
(J0820-1350: 4.8 vs 4.78; J1110-5637: 10.4 vs 8.58). Use them to read the map,
not as measurements — and only when the travel is actually detected, since on a
noise map the crossings are meaningless.

`direction` is the sign of A at the smallest lags (positive = later longitudes
light up later = positive drift, early → later longitudes, Szary+2022); it is
read off the rank-1 reconstruction rather than a single noisy matrix element.
"""
function ridge(A::AbstractMatrix)
    F = svd(A)
    s = F.S
    rank1_frac = sum(abs2, s) > 0 ? s[1]^2 / sum(abs2, s) : 0.0
    u = F.U[:, 1]                      # τ-profile
    v = F.V[:, 1]                      # Δ-profile
    # fix the arbitrary joint sign so that the τ-profile starts positive
    if u[1] < 0
        u = -u; v = -v
    end
    first_zero = function (y)
        for i in 1:length(y)-1
            if y[i] != 0 && sign(y[i+1]) != sign(y[i])
                return i + abs(y[i]) / (abs(y[i]) + abs(y[i+1]))
            end
        end
        return NaN
    end
    zt, zd = first_zero(u), first_zero(v)
    return (rank1_frac = rank1_frac,
            tau_mode   = u,
            dphi_mode  = v,
            p3         = isnan(zt) ? NaN : 2 * zt,
            p2         = isnan(zd) ? NaN : 2 * zd,
            p2_lower_limit = isnan(zd) ? 2.0 * length(v) : NaN,
            direction  = sign(u[1] * v[1]))
end


"""
    fit_geometry(E, K00, max_lag, max_dphi, npulses, non; n2, n3) -> NamedTuple

Find the (|P2|, P3) at which the pulsar has the most coherent modulation, by
maximising the *even* projection over a grid.

Fitting on the even half and not the odd one is the whole point: the even half
measures how much coherent modulation sits at a given geometry regardless of
whether it travels, so choosing the geometry this way cannot manufacture travel.
The odd projection is then measured at the geometry the pulsar itself picks, and
its sign — not the grid — decides the drift direction.

This replaces taking P3 from a catalogue, which fails badly: across the TPA
sample the ratio of the map's own P3 to the catalogue value has median 0.37 for
pulsars that then read R < 0.3, against 1.02 for those reading R > 0.7. A
mismatched P3 puts the template's τ oscillation out of phase with the map and
drives the projection *negative*, which looks exactly like "no travel" and is
not.

Evaluated in closed form rather than by building 1600 templates. With
W[t,d] = E[t,d]·(1−t/N)(1−d/M), the numerator factorises into two matrix
products, and ‖template‖² separates exactly into a τ factor times a Δ factor.
"""
function fit_geometry(E::AbstractMatrix, K00::Real, max_lag::Int, max_dphi::Int,
                      npulses::Int, non::Int; n2::Int=48, n3::Int=48,
                      p3_fixed::Union{Real,Nothing}=nothing,
                      p2_cap::Union{Real,Nothing}=nothing)
    p2top = isnothing(p2_cap) ? 8.0 * max_dphi : Float64(p2_cap)
    p2top > 4.0 || error("p2_cap must exceed 4 bins (got $p2top)")
    p2g = exp.(range(log(4.0), log(p2top), length=n2))
    p3g = isnothing(p3_fixed) ?
          exp.(range(log(2.0), log(max(4.0, 2.0 * max_lag)), length=n3)) :
          [Float64(p3_fixed)]
    n3  = length(p3g)
    tapt = [1 - t / npulses for t in 1:max_lag]
    tapd = [1 - d / non      for d in 1:max_dphi]

    W  = E .* tapt .* tapd'                                  # (τ × Δ)
    Cd = [cos(2π * d / p2g[i]) for d in 1:max_dphi, i in 1:n2]
    Ct = [cos(2π * t / p3g[j]) for t in 1:max_lag,  j in 1:n3]
    num = 2 .* (Ct' * (W * Cd))                              # ⟨E, te⟩, (p3 × p2)

    nd = [sum(tapd[d]^2 * cos(2π * d / p2g[i])^2 for d in 1:max_dphi) for i in 1:n2]
    nt = [sum(tapt[t]^2 * cos(2π * t / p3g[j])^2 for t in 1:max_lag)  for j in 1:n3]
    den = 4 .* (nt * nd')                                    # ‖te‖²

    # Selection must not be won by the smooth bulk of E. At large P2 and P3 the
    # even template is nearly constant and simply matches the always-positive
    # body of the correlation, which is how J2053-7200 (true P3 = 3.06) fitted
    # P3 = 58 and read a spurious R = 0.12. So the geometry is chosen by the
    # matched-filter amplitude of the template *orthogonalised against the
    # constant* (taper-only) model, which that corner cannot exploit.
    # Measurement then uses the plain template, keeping the equal-coefficient
    # identity that R relies on; where the template oscillates over many cycles
    # ⟨te,t0⟩ ≈ 0 anyway and the two coincide.
    md0 = [sum(tapd[d]^2 * cos(2π * d / p2g[i]) for d in 1:max_dphi) for i in 1:n2]
    mt0 = [sum(tapt[t]^2 * cos(2π * t / p3g[j]) for t in 1:max_lag)  for j in 1:n3]
    te_t0 = 2 .* (mt0 * md0')                                # ⟨te, t0⟩
    t0_t0 = sum(abs2, tapt) * sum(abs2, tapd)                # ⟨t0, t0⟩
    E_t0  = sum(W)                                           # ⟨E, t0⟩

    num_o  = num .- te_t0 .* (E_t0 / t0_t0)
    den_o  = den .- te_t0 .^ 2 ./ t0_t0
    S = num_o ./ sqrt.(max.(den_o, eps()))                   # matched-filter amplitude
    S[.!isfinite.(S)] .= -Inf
    j, i = Tuple(argmax(S))
    return (p2 = p2g[i], p3 = p3g[j], frac_even = num[j, i] / (K00 * den[j, i]),
            at_bound = i >= n2 - 1)
end


"""
    travel_test(data, bin_st, bin_end; kwargs...) -> NamedTuple

Test whether the subpulse pattern travels in longitude, from the time asymmetry
of the two-dimensional correlation (see the module docstring for why this is
exactly zero under amplitude modulation).

Nothing here needs P3: it is neither an argument nor used internally, so pulsars
whose modulation wobbles, changes period, switches mode or nulls are handled
without special treatment, and there is no FFT bin to miss.

Arguments:
  data      – single-pulse matrix (N_pulses × N_bins), real intensity
  bin_st    – first on-pulse bin (1-indexed)
  bin_end   – last on-pulse bin (1-indexed)

Keywords:
  max_lag    – largest pulse lag τ (default 40). If P3 is known, 2–3·P3 is a
               good choice: lags far beyond the coherence time only add noise
  max_dphi   – largest longitude lag Δ (default: half the on-pulse width)
  hp_halfwin – half-width of the running-mean high-pass (default 50 pulses);
               deliberately mild, it removes gain wander without touching the
               modulation. 0 = static profile only
  nreal      – surrogate realisations (default 500)
  seed       – RNG seed (default 7, nothing for non-reproducible)
  nblocks    – contiguous pulse blocks for `T_inc` and the consistency check
               (default 4). Set it so a block is shorter than the expected
               drift episode, otherwise a reversal hides inside one block
  p2_template, p3_template – claimed P2 (in longitude bins, signed) and P3 (in
               pulse periods). Supplying both switches on the matched
               projection, which is what makes a *demotion* possible; omitting
               them leaves the `frac` fields NaN. Pass `p2_template=:auto` to
               fit both from the pulsar's own map with `fit_geometry` (P3 is
               then ignored) — the right choice for a blind batch, since a
               catalogue P3 that disagrees with the map drives the projection
               negative. Note the geometry is then chosen on the same data the
               projection is measured from, which biases `frac_even` up and so
               `R` slightly *down*; the surrogates use the fixed fitted
               template and do not carry that selection
  pulse_st, pulse_end – analyse only this pulse range (default: all)

Fields of the returned NamedTuple:
  on_bins, taus, dphis – the analysed ranges
  K             – correlation map from `corr_map`
  A             – travel map A(Δ,τ), rows τ = 1…max_lag, columns Δ = 1…max_dphi
  A_norm        – A divided by K(0,0), dimensionless
  T, T_null     – omnibus statistic and its surrogate distribution
  significance  – (T − mean(T_null)) / std(T_null)  [σ]
  p_value       – count(T_null ≥ T)/nreal (0 means < 1/nreal)
  T_inc, T_inc_null, significance_inc, p_value_inc
                – the same statistic summed *incoherently* over pulse blocks,
                  Σ_b Σ A_b². T itself is computed on one global map, so a
                  drifter that spends equal time in each sense cancels and reads
                  as no travel at all — a false demotion waiting to happen.
                  (J1750-3503 survives T only because its episodes are lopsided,
                  28 P one way against 88 P the other.) T_inc adds the blocks'
                  power instead of their maps, so opposite senses accumulate.
                  It pays for that with a higher noise floor — nb blocks of
                  noise instead of one coherent average — so it is the less
                  sensitive of the two for a steady drift. Quote both; taking
                  the better of the two costs a mild trials penalty
  frac, frac_err, frac_sig, frac_limit
                – matched projection onto `drift_template` at the claimed
                  (p2_template, p3_template), as the fraction of the observed
                  modulation power that sits in a coherent drift of exactly
                  that geometry: frac = ⟨A, template⟩ / (K(0,0)·‖template‖²),
                  so frac = 1 would mean the whole modulation is that drift.
                  `frac_err` is its surrogate scatter and `frac_limit` the 3σ
                  upper limit. This is the field that lets a "drift"
                  classification be *rejected* rather than just unconfirmed: if
                  the pulsar is labelled a drifter at some P2 and frac_limit
                  comes out at a few per cent, at most a few per cent of its
                  modulation can be travelling, which the label cannot survive.
                  Model-dependent, unlike T — the template is a single sinusoid
                  in each direction, so harmonics and coherence decay make a
                  real drift project onto it imperfectly and bias frac low.
                  Never demote on frac alone; require T and T_inc to be quiet
                  too, and require the modulation itself to be well detected,
                  otherwise the limit measures sensitivity rather than physics
  frac_even, frac_even_err, R, R_err
                – `frac_even` is the same projection onto the Δ-*even* half of
                  the drift model (`drift_template_even` against `sym_map`),
                  and R = frac / frac_even is the share of the pulsar's coherent
                  modulation at that (P2, P3) which actually travels: 1 for a
                  rigid drift, 0 for amplitude modulation.

                  R exists to break a circularity. Judging `frac` needs to know
                  what a genuine drifter returns, and the obvious reference
                  sample — the pulsars already labelled `drift` — is exactly the
                  set suspected of contamination, so calibrating on it teaches
                  that drifters can have frac ≈ 0 and destroys the power to
                  demote anything. R takes its yardstick from the same pulsar
                  instead: both halves of the drift model carry equal
                  coefficients, so modulation strength, harmonics, coherence
                  loss and the taper cancel in the ratio and no external sample
                  is needed. A pulsar that both drifts and has separable
                  modulation lands between 0 and 1, which is the physically
                  meaningful reading rather than a failure.

                  NaN when the even projection is not itself measured at 3σ —
                  without coherent modulation at the claimed geometry there is
                  no denominator and the question is empty. R_err propagates the
                  two errors ignoring their (positive) correlation, so it is
                  mildly conservative.

                  R is *not* a strict fraction and can overshoot 1. The
                  denominator collects the even projection of the non-travelling
                  modulation too, and that contribution carries no fixed sign:
                  an a(φ) whose longitude autocorrelation happens to project
                  negatively onto cos(2πΔ/P2) shrinks the denominator and pushes
                  R above 1 (seen at R = 1.32 for a synthetic drift plus a
                  monotonic-ramp modulation). The useful reading is therefore
                  ordinal, not literal: R ≈ 0 means the coherent modulation at
                  this geometry does not travel, R of order 1 means it does.
                  Crucially the pathology lives at the high end, while demotion
                  turns on the low end, where the numerator is what approaches
                  zero and the ratio stays well behaved
  rank1_frac, tau_mode, dphi_mode, p2, p2_lower_limit, p3, direction
                – structure of the map, from `ridge`
  block_proj    – leave-one-out projection of each pulse block's own map onto
                  the sum of the others: the time-resolved sign and strength of
                  the travel, the analogue of a slope(t) panel. Zero for noise,
                  one for a travel that repeats identically throughout. A
                  drifter that reverses shows both signs; two components whose
                  modulations are independent (so their relative lag wanders)
                  do not reproduce across blocks and land near zero even when
                  the global T looks significant
  block_consistency – mean of `block_proj`
  T_off, significance_off, p_value_off, T_inc_off, significance_inc_off
                – the identical statistics on an off-pulse strip, which must be
                  consistent with zero. A significant result here means the
                  noise is not time-reversal symmetric (gain drift, RFI, bad
                  baseline) and the on-pulse numbers cannot be trusted
"""
function travel_test(data::AbstractMatrix, bin_st::Int, bin_end::Int;
                     max_lag::Int=40, max_dphi::Union{Int,Nothing}=nothing,
                     hp_halfwin::Int=50, nreal::Int=500,
                     seed::Union{Int,Nothing}=7, nblocks::Int=4,
                     p2_template::Union{Real,Symbol,Nothing}=nothing,
                     p3_template::Union{Real,Nothing}=nothing,
                     p2_cap_frac::Real=8.0, orth_even::Bool=false,
                     pulse_st::Union{Int,Nothing}=nothing,
                     pulse_end::Union{Int,Nothing}=nothing)
    nbins = size(data, 2)
    ps = isnothing(pulse_st)  ? 1 : pulse_st
    pe = isnothing(pulse_end) ? size(data, 1) : pulse_end
    (1 <= ps < pe <= size(data, 1)) ||
        error("bad pulse range $ps:$pe for $(size(data,1)) pulses")
    on  = bin_st:bin_end
    M   = length(on)
    M >= 4 || error("Need at least 4 on-pulse bins (got $M)")
    md  = isnothing(max_dphi) ? max(1, M ÷ 2) : max_dphi

    Xhp = highpass(Float64.(@view data[ps:pe, :]), hp_halfwin)
    X   = Xhp[:, on]
    N   = size(X, 1)

    strips, wrapped = offpulse_strips(bin_st, bin_end, nbins, M)
    noise_var = var(vcat((vec(@view Xhp[:, ix]) for ix in strips[1:min(end, 8)])...))

    # the block split is fixed up front: the incoherent statistic, the
    # consistency check and the surrogates must all use the same one
    nb = clamp(nblocks, 1, max(1, N ÷ (4 * max_lag)))
    edges = round.(Int, range(1, N + 1, length=nb + 1))

    K   = corr_map(X, max_lag, md)
    K00, A, E, blocks = _travel_maps(X, max_lag, md, edges)
    T     = sum(abs2, A)
    E2    = sum(abs2, E)
    T_inc = _inc_power(blocks)

    r = ridge(A)
    # :auto takes the template geometry from the pulsar's own map. Convenient for
    # a blind batch, but the P2 is then chosen on the same data the projection is
    # measured from, which biases frac_odd up; the surrogates use the fixed
    # template and so do not carry that selection. Read R only where T or T_inc
    # actually detects something.
    p2_at_bound = false
    p2t, p3t = if p2_template === :auto
        g = fit_geometry(E, K00, max_lag, md, N, M; p3_fixed=p3_template,
                         p2_cap=p2_cap_frac * md)
        p2_at_bound = g.at_bound
        g.p2, g.p3          # |P2|; the sign is read off the odd projection below
    else
        p2_template, p3_template
    end
    have_t = !(isnothing(p2t) || isnothing(p3t) || isnan(p2t) || isnan(p3t))
    tmpl  = have_t ? drift_template(max_lag, md, p2t, p3t;
                                    npulses=N, non=M) : nothing
    tmple = have_t ? drift_template_even(max_lag, md, p2t, p3t;
                                         npulses=N, non=M) : nothing
    if have_t && orth_even
        # The even template stops oscillating once P2 exceeds the searched Δ
        # range; it then matches the smooth body of E instead of the modulation,
        # inflating the denominator and crushing R. Removing the component along
        # the taper-only ("constant") model kills that channel. It leaves the
        # equal-coefficient identity intact: for E ∝ te exactly,
        # ⟨te, te⊥⟩ = ‖te⊥‖², so a pure drift still measures frac_even = 1.
        t0 = [(1 - t / N) * (1 - d / M) for t in 1:max_lag, d in 1:md]
        tmple = tmple .- (dot(tmple, t0) / sum(abs2, t0)) .* t0
    end
    tnorm  = have_t ? sum(abs2, tmpl)  : 0.0
    tnorme = have_t ? sum(abs2, tmple) : 0.0
    fracof(Am)  = (!have_t || tnorm  <= 0 || K00 <= 0) ? NaN :
                  dot(Am, tmpl)  / (K00 * tnorm)
    fraceof(Em) = (!have_t || tnorme <= 0 || K00 <= 0) ? NaN :
                  dot(Em, tmple) / (K00 * tnorme)
    frac      = fracof(A)
    frac_even = fraceof(E)

    Xsig = separable_signal(X, noise_var)
    rng  = isnothing(seed) ? Random.default_rng() : MersenneTwister(seed)
    buf  = Matrix{Float64}(undef, N, M)
    Xn   = Matrix{Float64}(undef, N, M)

    T_null      = zeros(nreal)
    E2_null     = zeros(nreal)
    T_inc_null  = zeros(nreal)
    frac_null   = zeros(nreal)
    frace_null  = zeros(nreal)
    for i in 1:nreal
        noise_block!(buf, Xhp, strips, rng)
        Xn .= Xsig .+ buf
        _, As, Es, bs = _travel_maps(Xn, max_lag, md, edges)
        T_null[i]     = sum(abs2, As)
        E2_null[i]    = sum(abs2, Es)
        T_inc_null[i] = _inc_power(bs)
        frac_null[i]  = fracof(As)
        frace_null[i] = fraceof(Es)
    end
    significance     = (T - mean(T_null)) / std(T_null)
    p_value          = count(>=(T), T_null) / nreal
    significance_inc = (T_inc - mean(T_inc_null)) / std(T_inc_null)
    p_value_inc      = count(>=(T_inc), T_inc_null) / nreal
    frac_err      = have_t ? std(frac_null)  : NaN
    frac_sig      = have_t ? frac / frac_err : NaN
    frac_limit    = have_t ? max(frac, 0.0) + 3 * frac_err : NaN
    frac_even_err = have_t ? std(frace_null) : NaN
    # travelling share of the coherent modulation; error propagated ignoring the
    # (positive) correlation between numerator and denominator, so mildly
    # conservative. Undefined unless the even projection is itself well measured
    R     = (have_t && frac_even > 3 * frac_even_err) ? frac / frac_even : NaN
    R_err = isnan(R) ? NaN :
            abs(R) * sqrt((frac_err / frac)^2 + (frac_even_err / frac_even)^2)

    # rho — the same discrimination as R, read as an angle instead of a ratio.
    # The drift model gives the two halves equal coefficients, i.e. a direction
    # at 45 degrees in the (odd, even) plane, so
    #     rho = sqrt(2)*odd / hypot(odd, even)
    # is 1 for a rigid drift and 0 for amplitude modulation. Unlike R it is a
    # sine rather than a tangent: bounded by sqrt(2), with no pole when the even
    # projection is small, and always defined — which is what removes the need
    # for a "denominator significant" guard and for any quality thresholds.
    # R is kept alongside because the published relations are stated in it.
    dnorm = have_t ? hypot(frac, frac_even) : 0.0
    rho   = (have_t && dnorm > 0) ? sqrt(2) * frac / dnorm : NaN
    rho_err = isnan(rho) ? NaN :
              sqrt(2 * frac_even^2 * (frac_even^2 * frac_err^2 +
                                      frac^2 * frac_even_err^2)) / dnorm^3

    # off-pulse control: the same two statistics where there is no signal at all
    idx0 = strips[length(strips) ÷ 2 + 1]
    _, Aoff, _, boff = _travel_maps(@view(Xhp[:, idx0]), max_lag, md, edges)
    T_off     = sum(abs2, Aoff)
    T_inc_off = _inc_power(boff)
    T_off_null     = zeros(nreal)
    E2_off_null    = zeros(nreal)
    T_inc_off_null = zeros(nreal)
    for i in 1:nreal
        noise_block!(buf, Xhp, strips, rng)
        _, As, Es, bs = _travel_maps(buf, max_lag, md, edges)
        T_off_null[i]     = sum(abs2, As)
        E2_off_null[i]    = sum(abs2, Es)
        T_inc_off_null[i] = _inc_power(bs)
    end
    significance_off     = (T_off - mean(T_off_null)) / std(T_off_null)
    p_value_off          = count(>=(T_off), T_off_null) / nreal
    significance_inc_off = (T_inc_off - mean(T_inc_off_null)) / std(T_inc_off_null)

    # rho_free — the same discrimination with no template and no geometry at all.
    # A rigid drift puts equal coefficients in both halves of the model, so their
    # *norms* are equal too: ‖A‖² = 4Σtaper²sin²sin² and ‖E‖² = 4Σtaper²cos²cos²,
    # and ⟨sin²sin²⟩ = ⟨cos²cos²⟩ = 1/4 once a few cycles fit in the searched
    # ranges. Hence sqrt(‖A‖²/‖E‖²) = 1 for drift, 0 for amplitude modulation —
    # without fitting P2. That matters because fitting P2 is what made the result
    # depend on the on-pulse window: a generous window drags the fit to large P2,
    # and a too-large P2 cripples the odd channel specifically (sin → 0 at Δ → 0,
    # exactly where the signal sits, while cos → 1). Both norms carry a positive
    # noise bias, removed with the off-pulse maps that the control already builds.
    # Price: no P2, no drift direction, and more variance than a matched filter,
    # since noise-dominated grid cells enter the norms.
    A_noise = mean(T_off_null)
    E_noise = mean(E2_off_null)
    As_sig  = max(T  - A_noise, 0.0)
    Es_sig  = max(E2 - E_noise, 0.0)
    rho_free = Es_sig > 0 ? sqrt(As_sig / Es_sig) : NaN
    rho_free_err = (As_sig > 0 && Es_sig > 0 && !isnan(rho_free)) ?
        0.5 * rho_free * sqrt((std(T_null) / As_sig)^2 + (std(E2_null) / Es_sig)^2) : NaN

    # The equal-norms argument needs ⟨sin²(2πτ/P3)⟩ = ⟨cos²(2πτ/P3)⟩, which holds
    # only once several cycles fit in τ = 1…max_lag. Near the temporal Nyquist
    # limit it fails hard: at P3 = 2.05 the ratio is 0.160, so a genuine drift
    # reads rho_free = 0.41 instead of 1 and would be called amplitude modulation.
    # The factor depends only on P3 and max_lag — both known, no fitting — so it
    # divides out analytically. Verified on synthetic drift at P3 = 2.05…9:
    # raw 0.407/1.202/1.053/1.064/1.044 → corrected 1.017/1.019/1.021/1.020/1.020.
    # The matching Δ factor would need P2 and is deliberately NOT applied; it is
    # what makes rho_free fall off once P2 exceeds the searched longitude range.
    ctau = NaN
    if !isnothing(p3_template) && p3_template > 0
        tp2 = [(1 - t / N)^2 for t in 1:max_lag]
        sn  = sum(tp2[t] * sin(2π * t / p3_template)^2 for t in 1:max_lag)
        cs  = sum(tp2[t] * cos(2π * t / p3_template)^2 for t in 1:max_lag)
        (sn > 0 && cs > 0) && (ctau = sqrt(sn / cs))
    end
    rho_free_corr     = isnan(ctau) ? NaN : rho_free / ctau
    rho_free_corr_err = isnan(ctau) ? NaN : rho_free_err / ctau

    # time-resolved consistency, leave-one-out so that pure noise gives zero:
    # projecting a block onto the *global* map would keep the block's own
    # contribution and yield ~1/sqrt(nblocks) for noise alone
    block_proj = Float64[]
    if length(blocks) >= 2
        Asum = sum(blocks)
        for Ab in blocks
            Arest = Asum .- Ab
            nrm = sqrt(sum(abs2, Ab)) * sqrt(sum(abs2, Arest))
            push!(block_proj, nrm > 0 ? dot(Ab, Arest) / nrm : 0.0)
        end
    end

    return (
        on_bins      = on,
        taus         = 1:max_lag,
        dphis        = 1:md,
        K            = K,
        A            = A,
        A_norm       = A ./ K[1, md+1],
        T            = T,
        T_null       = T_null,
        significance = significance,
        p_value      = p_value,
        T_inc        = T_inc,
        T_inc_null   = T_inc_null,
        significance_inc = significance_inc,
        p_value_inc  = p_value_inc,
        p2_template  = p2t,
        p3_template  = p3t,
        frac         = frac,
        frac_err     = frac_err,
        frac_sig     = frac_sig,
        frac_limit   = frac_limit,
        frac_even     = frac_even,
        frac_even_err = frac_even_err,
        R            = R,
        R_err        = R_err,
        rho          = rho,
        rho_err      = rho_err,
        E2           = E2,
        rho_free     = rho_free,
        rho_free_err = rho_free_err,
        ctau         = ctau,
        rho_free_corr     = rho_free_corr,
        rho_free_corr_err = rho_free_corr_err,
        rank1_frac   = r.rank1_frac,
        tau_mode     = r.tau_mode,
        dphi_mode    = r.dphi_mode,
        p2           = r.p2,
        p2_lower_limit = r.p2_lower_limit,
        p3           = r.p3,
        direction    = r.direction,
        block_proj   = block_proj,
        block_consistency = isempty(block_proj) ? NaN : mean(block_proj),
        T_off        = T_off,
        significance_off = significance_off,
        p_value_off  = p_value_off,
        T_inc_off    = T_inc_off,
        significance_inc_off = significance_inc_off,
        pulse_range  = (ps, pe),
        offpulse_wrapped = wrapped,
        p2_at_bound  = p2_at_bound,
    )
end


"""
    selftest(; verbose=true) -> Bool

Checks the two claims the method stands on, on synthetic data:

  1. `corr_map` (FFT, zero-padded) equals `corr_map_direct` (explicit sum);
  2. A ≡ 0 to machine precision for a separable field a(φ)·w(n) — including a
     w(n) whose period wanders, and an a(φ) that changes sign (antiphase
     amplitude modulation, the case that fools power-based tests);
  3. A is clearly nonzero for a travelling pattern.
"""
function selftest(; verbose::Bool=true)
    rng = MersenneTwister(42)
    ok = true

    X = randn(rng, 40, 12)
    d1 = maximum(abs.(corr_map(X, 5, 4) .- corr_map_direct(X, 5, 4)))
    scale = maximum(abs.(corr_map_direct(X, 5, 4)))
    ok &= d1 / scale < 1e-10
    verbose && println("FFT vs direct sum:        rel. diff $(d1/scale)")

    N, M = 600, 40
    phi = collect(1:M)
    # wandering period, nulls, and a sign flip half way across the profile
    p3t = [12.0 + 5 * sin(2π * n / 350) for n in 1:N]
    w = [(n > 200 && n < 260) ? 0.0 : cos(2π * sum(1 ./ p3t[1:n])) for n in 1:N]
    a = [ph < M / 2 ? 1.0 : -0.7 for ph in phi]
    Xsep = a' .* w
    Asep = antisym_map(corr_map(Xsep, 30, 15))
    rel = maximum(abs.(Asep)) / corr_map(Xsep, 1, 1)[1, 2]
    ok &= rel < 1e-10
    verbose && println("separable field, max|A|:  $(rel) of K(0,0)")

    Xdr = [cos(2π * (ph / 18 - n / 12)) for n in 1:N, ph in phi]
    Adr = antisym_map(corr_map(Xdr, 30, 15))
    rel_dr = maximum(abs.(Adr)) / corr_map(Xdr, 1, 1)[1, 2]
    ok &= rel_dr > 0.1
    verbose && println("drifting pattern, max|A|: $(rel_dr) of K(0,0)")

    r = ridge(Adr)
    verbose && println("recovered P2 = $(round(r.p2, digits=1)) bins (true 18), " *
                       "P3 = $(round(r.p3, digits=1)) (true 12), " *
                       "rank1 = $(round(r.rank1_frac, digits=3))")

    # matched projection recovers the full modulation power at the true geometry
    K00 = corr_map(Xdr, 1, 1)[1, 2]
    tm  = drift_template(30, 15, 18, 12; npulses=N, non=M)
    f   = dot(Adr, tm) / (K00 * sum(abs2, tm))
    ok &= 0.9 < f < 1.1
    verbose && println("matched frac at true P2/P3: $(round(f, digits=3)) (expect ~1)")

    # self-normalising ratio: 1 for a pure drift, 0 for separable modulation,
    # intermediate for a mixture -- and no reference sample anywhere
    tme = drift_template_even(30, 15, 18, 12; npulses=N, non=M)
    Rof = function (Xf)
        Kf = corr_map(Xf, 30, 15)
        k0 = Kf[1, 16]
        (dot(antisym_map(Kf), tm) / (k0 * sum(abs2, tm))) /
        (dot(sym_map(Kf), tme)   / (k0 * sum(abs2, tme)))
    end
    R_drift = Rof(Xdr)
    ok &= 0.9 < R_drift < 1.1
    # travelling wave plus a standing wave of the same periods — the textbook
    # mixture. A standing wave is half a forward plus half a backward traveller,
    # so amplitudes a_f = 1 + s/2 and a_b = s/2 give the exact prediction
    #   R = (a_f^2 - a_b^2) / (a_f^2 + a_b^2)
    s = 2.0
    Xmix = Xdr .+ [s * cos(2π * ph / 18) * cos(2π * n / 12) for n in 1:N, ph in phi]
    af, ab = 1 + s / 2, s / 2
    R_pred = (af^2 - ab^2) / (af^2 + ab^2)
    R_mix  = Rof(Xmix)
    ok &= abs(R_mix - R_pred) < 0.1
    verbose && println("R pure drift: $(round(R_drift, digits=3)) (expect ~1);  " *
                       "drift + standing wave: $(round(R_mix, digits=3)) " *
                       "(predicted $(round(R_pred, digits=3)))")

    # a balanced reverser: the global map cancels, the incoherent block sum does not
    Nr = 800
    Xrev = vcat([cos(2π * (ph / 18 - n / 12)) for n in 1:Nr÷2, ph in phi],
                [cos(2π * (-ph / 18 - n / 12)) for n in 1:Nr÷2, ph in phi])
    edges = [1, Nr÷2 + 1, Nr + 1]
    _, Arev, _, brev = _travel_maps(Xrev, 30, 15, edges)
    ratio = sum(abs2, Arev) / _inc_power(brev)
    ok &= ratio < 0.05
    verbose && println("balanced reverser, T/T_inc: $(round(ratio, sigdigits=2)) " *
                       "(global map cancels, T_inc does not)")

    verbose && println(ok ? "selftest PASSED" : "selftest FAILED")
    return ok
end

end  # module Travel
