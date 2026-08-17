module PhaseDrift

using FFTW
using Statistics
using Random
using LinearAlgebra

"""
    complex_lrfs_at_f3(data, p3) -> Vector{ComplexF64}

FFT along the pulse axis (dim 1) and extract the complex coefficient at
frequency f3 = 1/P3 for every longitude bin.

Frequency index k = round(N_pulses / P3), Julia 1-indexed as k+1.
Result shape: (N_bins,)

Sign convention: the coefficient is conjugated, i.e. the effective kernel is
e^(+2πi f3 n). Positive drift in the publication sense (subpulses moving from
early to later longitudes, Szary+2022) means later longitudes light up later
in time, and with the +i kernel a time delay is a *positive* phase shift —
so dψ/dφ > 0 ⇔ positive drift. With the raw forward-FFT kernel e^(-2πi f3 n)
the sign would come out inverted.
"""
function complex_lrfs_at_f3(data::AbstractMatrix, p3::Real)
    N = size(data, 1)
    F = fft(data, 1)                      # (N_pulses, N_bins)
    k = round(Int, N / p3)
    k = clamp(k, 1, N ÷ 2)
    return conj.(F[k + 1, :])             # +1: Julia 1-indexed; conj: see docstring
end


"""
    off_pulse_sigma(L, bin_st, bin_end) -> Float64

Estimate complex noise sigma from the off-pulse region of a complex LRFS slice.
Off-pulse bins are those outside [bin_st, bin_end] (1-indexed).
σ is the std of the real and imaginary parts pooled together.
"""
function off_pulse_sigma(L::AbstractVector, bin_st::Int, bin_end::Int)
    N_bins = length(L)
    off = vcat(1:bin_st-1, bin_end+1:N_bins)
    isempty(off) && error(
        "No off-pulse bins available (bin_st=$bin_st bin_end=$bin_end N_bins=$N_bins)")
    L_off = L[off]
    return std([real.(L_off); imag.(L_off)])
end


"""
    k_from_p3(N, p3) -> Int

Frequency index k = round(N/P3), clamped to the valid FFT range [1, N÷2].
Same convention as `complex_lrfs_at_f3`.
"""
function k_from_p3(N::Int, p3::Real)
    k = round(Int, N / p3)
    return clamp(k, 1, N ÷ 2)
end


"""
    circ_std(angles) -> Float64

Circular standard deviation of a set of angles [rad], sqrt(-2 log R) with
R the mean resultant length. Safe under phase wrapping, unlike a plain std.
"""
function circ_std(angles::AbstractVector)
    n = length(angles)
    n == 0 && return 0.0
    R = clamp(abs(sum(exp.(1im .* angles))) / n, 1e-12, 1.0)
    return sqrt(-2 * log(R))
end


"""
    p3_error_phases(data, p3, p3_error, bin_st, bin_end; n_samples) -> NamedTuple

Phase profiles ψ(φ) swept densely over the P3 uncertainty interval
(p3 ± p3_error). `drift_test` locks onto a single integer frequency index
k = round(N/P3), so a continuous P3 only ever produces a handful of distinct
k's — but those k's are not equally likely a priori: under a uniform prior
on P3 within its error bar, some k's correspond to a wider slice of that
interval than others. Sampling `n_samples` P3 values evenly across
[p3-p3_error, p3+p3_error] (rather than each distinct k once) reproduces
that weighting naturally, since the more probable k's simply recur more
often among the samples, both in the overlay scatter and in psi_sigma.

Returns (empty arrays if p3_error <= 0):
  k_grid    – k index used by each of the n_samples sweep points
  p3_grid   – the n_samples P3 values themselves
  psi_grid  – phase array, size (n_samples, length(on_bins)) [rad]
  psi_sigma – per-bin circular std across the sweep [rad], the panel-2 error bar
"""
function p3_error_phases(data::AbstractMatrix, p3::Real, p3_error::Real,
                         bin_st::Int, bin_end::Int; n_samples::Int=151)
    on_bins = bin_st:bin_end
    N = size(data, 1)
    if !(p3_error > 0)
        return (k_grid=Int[], p3_grid=Float64[],
                psi_grid=zeros(0, length(on_bins)), psi_sigma=zeros(length(on_bins)))
    end
    p3_grid = collect(range(p3 - p3_error, p3 + p3_error, length=n_samples))
    F = fft(data, 1)
    psi_grid = Array{Float64}(undef, n_samples, length(on_bins))
    k_grid = Vector{Int}(undef, n_samples)
    for (i, p3i) in enumerate(p3_grid)
        k = k_from_p3(N, p3i)
        k_grid[i] = k
        # conj: same e^(+2πi f3 n) sign convention as `complex_lrfs_at_f3`
        psi_grid[i, :] = angle.(conj.(F[k + 1, on_bins]))
    end
    psi_sigma = [circ_std(psi_grid[:, j]) for j in axes(psi_grid, 2)]
    return (k_grid=k_grid, p3_grid=p3_grid, psi_grid=psi_grid, psi_sigma=psi_sigma)
end


"""
    slope_stat(L_on) -> Float64

Coherent estimator of the phase gradient dψ/dφ [rad/bin].

  g = Σ_φ  conj(L[φ]) · L[φ+1]
  slope = arg(g)

Naturally weighted by |L[φ]||L[φ+1]|, robust against phase wrapping.

With `complex_lrfs_at_f3`'s sign convention: slope > 0 ⇔ positive drift,
i.e. subpulses moving from early to later longitudes (Szary+2022).
"""
function slope_stat(L_on::AbstractVector)
    n = length(L_on)
    g = sum(conj(L_on[i]) * L_on[i+1] for i in 1:n-1)
    return angle(g)
end


"""
    wrap_pi(x) -> Float64

Wrap an angle into (-π, π].
"""
wrap_pi(x::Real) = atan(sin(x), cos(x))


"""
    local_slopes(L_on, sigma_off; snr_min) -> NamedTuple

Local phase gradient, longitude bin by longitude bin — the resolved version
of `slope_stat`, which collapses the same quantity into a single number:

  Δψ_j = arg( conj(L[j]) · L[j+1] )     [rad/bin],  j = 1 … n-1

This is what discriminates a *systematic* drift from a phase *jump*. Pure
amplitude modulation gives a ψ(φ) that is piecewise constant with 180° steps
wherever |L| passes through a node, so Δψ is ~0 everywhere except at isolated
spikes; a genuine subpulse drift gives Δψ flat at the global slope.

Errors from noise propagation: for L = A·exp(iψ) + n with n complex and
σ_off per component, σ_ψ ≈ σ_off/A, hence

  σ_Δψ,j = σ_off · sqrt(1/A_j² + 1/A_{j+1}²)

That Gaussian-phase approximation breaks down once A ≲ 3σ_off (the phase
error saturates at the uniform-circle value ~104°), so increments touching
such a bin are flagged as unusable rather than trusted.

Fields:
  dpsi       – local increments Δψ_j [rad], length n-1
  sigma      – 1σ error on each Δψ_j [rad]
  good       – Bool mask, both bins of the increment above snr_min·σ_off
  phase_var  – per-bin phase variance σ_off²/A² [rad²], length n (for the
               covariance matrix in `slope_chi2`)
"""
function local_slopes(L_on::AbstractVector, sigma_off::Real; snr_min::Real=3.0)
    n = length(L_on)
    n < 3 && error("Need at least 3 on-pulse bins for local slopes (got $n)")
    A = abs.(L_on)
    phase_var = (sigma_off ./ max.(A, eps())) .^ 2
    dpsi  = [angle(conj(L_on[j]) * L_on[j+1]) for j in 1:n-1]
    sigma = sqrt.(phase_var[1:n-1] .+ phase_var[2:n])
    good  = [A[j] > snr_min * sigma_off && A[j+1] > snr_min * sigma_off for j in 1:n-1]
    return (dpsi=dpsi, sigma=sigma, good=good, phase_var=phase_var)
end


"""
    slope_chi2(ls, slope) -> NamedTuple

Goodness of fit of a *constant* phase gradient to the local increments of
`local_slopes`, i.e. the "is the drift systematic?" statistic.

Consecutive increments are not independent — Δψ_j and Δψ_{j+1} share the bin
j+1 — so the residuals must be weighted by the full covariance matrix, not by
σ_Δψ alone. With v_j = σ_off²/A_j² the per-bin phase variance:

  Cov(Δψ_j, Δψ_j)   =  v_j + v_{j+1}
  Cov(Δψ_j, Δψ_l)   = -v_{j+1}   for l = j+1  (and symmetric)
  Cov(Δψ_j, Δψ_l)   =  0         otherwise

  χ² = rᵀ C⁻¹ r,   r_j = wrap(Δψ_j - slope)

Only increments flagged `good` enter; the covariance is built from the general
shares-a-bin rule, so a non-contiguous mask is handled correctly. One free
parameter (the slope itself) is subtracted from the dof.

χ²_red ≈ 1  → phase changes systematically (constant gradient fits).
χ²_red ≫ 1  → phase jumps with longitude.

Returns (chi2, dof, chi2_red, n_used); all NaN/0 if fewer than 3 usable
increments survive the S/N mask.
"""
function slope_chi2(ls::NamedTuple, slope::Real)
    idx = findall(ls.good)
    m   = length(idx)
    m < 3 && return (chi2=NaN, dof=0, chi2_red=NaN, n_used=m)
    v = ls.phase_var
    C = zeros(m, m)
    for a in 1:m, b in 1:m
        j, l = idx[a], idx[b]
        if j == l
            C[a, b] = v[j] + v[j+1]
        elseif l == j + 1
            C[a, b] = -v[j+1]
        elseif j == l + 1
            C[a, b] = -v[j]
        end
    end
    r    = [wrap_pi(ls.dpsi[j] - slope) for j in idx]
    chi2 = dot(r, C \ r)
    dof  = m - 1
    return (chi2=chi2, dof=dof, chi2_red=chi2 / dof, n_used=m)
end


"""
    amp_null_slopes(L_on, sigma_off; nreal, seed) -> Vector{Float64}

Amplitude-modulation null distribution.

Model: keep the measured amplitude with the noise contribution removed
(|L|² is chi²-biased upward by 2σ²; using raw |L_on| would make the null
phases artificially stable and inflate the significance when the f3
feature is weak — pure noise then reads as several σ), set phase to zero
(flat), add complex noise N(0, σ_off).
Returns nreal null slope values [rad/bin].
"""
function amp_null_slopes(L_on::AbstractVector, sigma_off::Real;
                         nreal::Int=6000, seed::Union{Int,Nothing}=7)
    rng = isnothing(seed) ? Random.default_rng() : MersenneTwister(seed)
    A = sqrt.(max.(abs2.(L_on) .- 2 * sigma_off^2, 0.0))
    n = length(A)
    slopes = zeros(nreal)
    for i in 1:nreal
        noise = sigma_off .* (randn(rng, n) .+ 1im .* randn(rng, n))
        slopes[i] = slope_stat(A .+ noise)
    end
    return slopes
end


"""
    drift_test(data, p3, bin_st, bin_end; p3_error, nreal, seed) -> NamedTuple

Phase-drift vs amplitude-modulation discriminator.

Computes the complex LRFS at f3=1/P3, measures the coherent phase slope,
and compares it against a null distribution for pure amplitude modulation
(flat phase + off-pulse noise). Returns significance in σ units.

If `p3_error` > 0, also computes the phase profiles ψ(φ) for every FFT
harmonic k consistent with P3 ± p3_error (see `p3_error_phases`); their
spread is the phase error induced by the P3 uncertainty, used for the
panel-2 error bars and faint overlay points.

Arguments:
  data      – single-pulse matrix (N_pulses × N_bins), real intensity
  p3        – P3 period [pulse periods P0]
  bin_st    – first on-pulse bin (1-indexed)
  bin_end   – last on-pulse bin (1-indexed)
  p3_error  – P3 uncertainty [pulse periods P0] (default 0.0, disabled)
  nreal     – number of null realisations (default 6000)
  seed      – RNG seed (default 7, nothing for non-reproducible)

Fields of returned NamedTuple:
  L              – complex LRFS slice at f3, all N_bins
  on_bins        – on-pulse bin UnitRange
  L_on           – L restricted to on-pulse
  sigma_off      – off-pulse noise estimate in complex L
  snr            – f3 feature strength, mean(|L_on|)/σ_off; pure noise gives
                   ~1.25 (Rayleigh mean), a warning is issued below 3
  slope          – measured phase slope [rad/bin]; > 0 = positive drift
                   (early → later longitudes, Szary+2022 convention)
  null           – null slope distribution [nreal]
  significance   – slope / std(null)  [σ]
  p3             – P3 used
  p3_error       – P3 uncertainty used
  p3err_k        – k index used by each densely-sampled P3 sweep point
  p3err_p3       – the swept P3 values themselves
  p3err_psi      – phase array from the P3 sweep [rad], (length(p3err_k) × N_on)
  p3err_sigma    – per-bin circular std across the P3 sweep [rad]
  dpsi           – local phase gradient per bin pair [rad/bin] (`local_slopes`)
  dpsi_sigma     – 1σ noise error on dpsi [rad/bin]
  dpsi_good      – mask: both bins of the increment above snr_min·σ_off
  chi2, chi2_dof, chi2_red, chi2_n
                 – fit of a constant gradient to dpsi (`slope_chi2`);
                   χ²_red ~ 1 systematic drift, ≫ 1 phase jumps
"""
function drift_test(data::AbstractMatrix, p3::Real, bin_st::Int, bin_end::Int;
                    p3_error::Real=0.0, nreal::Int=6000, seed::Union{Int,Nothing}=7,
                    snr_min::Real=3.0)
    on_bins   = bin_st:bin_end
    L         = complex_lrfs_at_f3(data, p3)
    sigma_off = off_pulse_sigma(L, bin_st, bin_end)
    L_on      = L[on_bins]
    snr       = sum(abs.(L_on)) / length(L_on) / sigma_off
    snr < 3 && @warn "Weak f3 feature (SNR = $(round(snr, digits=2)), pure noise gives ~1.25): " *
                     "the drift-vs-modulation significance is not meaningful without a detected P3 feature"
    slope     = slope_stat(L_on)
    null      = amp_null_slopes(L_on, sigma_off; nreal=nreal, seed=seed)
    significance = slope / std(null)
    p3err = p3_error_phases(data, p3, p3_error, bin_st, bin_end)
    ls    = local_slopes(L_on, sigma_off; snr_min=snr_min)
    chi2  = slope_chi2(ls, slope)
    return (
        L            = L,
        on_bins      = on_bins,
        L_on         = L_on,
        sigma_off    = sigma_off,
        snr          = snr,
        slope        = slope,
        null         = null,
        significance = significance,
        p3           = p3,
        p3_error     = p3_error,
        p3err_k      = p3err.k_grid,
        p3err_p3     = p3err.p3_grid,
        p3err_psi    = p3err.psi_grid,
        p3err_sigma  = p3err.psi_sigma,
        dpsi         = ls.dpsi,
        dpsi_sigma   = ls.sigma,
        dpsi_good    = ls.good,
        chi2         = chi2.chi2,
        chi2_dof     = chi2.dof,
        chi2_red     = chi2.chi2_red,
        chi2_n       = chi2.n_used,
    )
end

"""
    windowed_lrfs(data, f3, window, stride) -> NamedTuple

Sliding-window complex LRFS at exactly f3 [cycles/pulse]: for every window
position b the complex amplitude

  L_b(φ) = Σ_{n∈b} data[n,φ] · e^(+2πi f3 n)

is computed for all longitude bins via one carrier multiplication and a
cumulative sum along the pulse axis, so the cost is O(N·N_bins) regardless
of `stride`. There is no quantisation of f3 to an integer FFT index, and
the short window's wide response (main lobe f3 ± 1/window) captures a
wobbling P3 without loss. The +i kernel matches `complex_lrfs_at_f3`'s
sign convention: dψ/dφ > 0 ⇔ positive drift (early → later longitudes,
Szary+2022) — see there.

DC and slow intensity fluctuations must be removed from `data` *before*
calling this (`drift_test_sliding` uses `subtract_running_mean`): for
window < P3 the window's response at f3 overlaps f = 0, so slow power
would leak straight into every L_b. Deliberately there is no *per-window*
demeaning — with under one carrier cycle per window the f3 signal is
largely indistinguishable from an offset, so per-window demeaning would
remove much of the signal along with the DC. Residual leakage has a
longitude-stationary pattern and lands in Re(g) of `drift_test_sliding`'s
statistic, which only uses Im(g) — see there.

Returns (L, starts): L is nwin × N_bins complex, window b spans pulses
starts[b] : starts[b]+window-1.
"""
function windowed_lrfs(data::AbstractMatrix, f3::Real, window::Int, stride::Int)
    N, nbins = size(data)
    window < 2       && error("window must be ≥ 2 (got $window)")
    stride < 1       && error("stride must be ≥ 1 (got $stride)")
    N < window       && error("Need at least window=$window pulses (got $N)")
    carrier = [exp(2im * π * f3 * n) for n in 1:N]
    C = cumsum(data .* carrier, dims=1)
    starts = collect(1:stride:N-window+1)
    L = Array{ComplexF64}(undef, length(starts), nbins)
    for (b, s) in enumerate(starts)
        e = s + window - 1
        if s == 1
            L[b, :] = C[e, :]
        else
            L[b, :] = C[e, :] .- C[s-1, :]
        end
    end
    return (L=L, starts=starts)
end


"""
    subtract_running_mean(X, halfwin) -> Matrix{Float64}

Per-bin high-pass: subtract from every time series X[:, j] its centred
running mean over 2·halfwin+1 pulses (edge-truncated). With halfwin ≈ P3
the running-mean's frequency response has a null at f3 = 1/P3 (its length
is then ≈ 2 cycles of f3), so the subtraction removes the static profile
*and* slow intensity fluctuations (nulling, scintillation) while leaving
the f3 modulation essentially untouched (~1% for P3 ≫ 1). This matters for
short-window LRFS: slow power leaking into a window biases the coherent
slope arg(g) toward zero (it adds a longitude-stationary, flat-phase
component to Re g). halfwin ≤ 0 falls back to plain mean subtraction.
"""
function subtract_running_mean(X::AbstractMatrix, halfwin::Int)
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
    null_amplitudes(L_on, sigma_win; amp_smooth) -> Matrix{Float64}

Per-window amplitude profiles for the flat-phase surrogates of
`drift_test_sliding` / `drift_test_profile`: |L_b(φ)| with the noise
contribution removed, returned transposed (n_on × nwin) for the access
order of the surrogate loops.

The null variance of the drift quadrature is driven by the amplitude
*gradient* Σ_j (A[j-1] - A[j+1])², and the raw per-bin |L| jitters at the
noise level — that jitter would inflate the gradient, and with it the null,
by a large factor. So the power profile is boxcar-smoothed over
±`amp_smooth` bins in longitude before the debias, and then rescaled per
window so the total power still matches the unsmoothed, debiased sum. A
final global rescale removes the bias left by the per-window truncations
at zero (which otherwise keep the upward noise fluctuations of the window
power). At mean(S)/std(S) ~ 100 a few-percent power error is worth several
σ, so none of these corrections is cosmetic.
"""
function null_amplitudes(L_on::AbstractMatrix, sigma_win::AbstractVector;
                         amp_smooth::Int=5)
    nwin, non = size(L_on)
    At    = Matrix{Float64}(undef, non, nwin)
    pow_u = Vector{Float64}(undef, nwin)
    for b in 1:nwin
        p2 = abs2.(@view L_on[b, :])
        pow_u[b] = sum(p2) - 2 * sigma_win[b]^2 * non
        for j in 1:non
            lo, hi = max(1, j - amp_smooth), min(non, j + amp_smooth)
            At[j, b] = sqrt(max(mean(@view p2[lo:hi]) - 2 * sigma_win[b]^2, 0.0))
        end
        pow_t = sum(abs2, @view At[:, b])
        At[:, b] .*= pow_t > 0 ? sqrt(max(pow_u[b], 0.0) / pow_t) : 0.0
    end
    tot_clip = sum(x -> max(x, 0.0), pow_u)
    At .*= tot_clip > 0 ? sqrt(max(sum(pow_u), 0.0) / tot_clip) : 0.0
    return At
end


"""
    null_im_g!(im_g, At, L_offt, pidx, rot, rng, replace) -> im_g

One flat-phase surrogate realisation for `drift_test_sliding`: fills `im_g`
with Im g_b per window, where L_b^null(φ_j) = A_b(φ_j) + N_b(φ_j). A_b are
the measured (debiased, power-matched, real) per-window amplitudes in `At`
(n_on × nwin, transposed for access order). The noise N_b is *bootstrapped
from the pulsar's own off-pulse windowed LRFS* `L_offt` (n_off × nwin):
each on-pulse bin j is assigned a random off-pulse bin (a fresh draw
`pidx`, without replacement when n_off suffices) and a random phase
rotation `rot[j]`, both held fixed across windows within the realisation.
Off-pulse bins carry no signal but went through the identical high-pass +
windowed-DFT pipeline, so N_b has exactly the real noise level, spectrum
and window-to-window correlation — overlapping windows share noise just as
in the data, which is what makes the null spread of S = Σ|Im g_b| honest
for stride < window, with no noise model or calibration constant at all.
"""
function null_im_g!(im_g::Vector{Float64}, At::Matrix{Float64},
                    L_offt::Matrix{ComplexF64}, pidx::Vector{Int},
                    rot::Vector{ComplexF64}, rng, replace::Bool)
    non, nwin = size(At)
    noff = size(L_offt, 1)
    if replace
        rand!(rng, pidx, 1:noff)
    else
        # first non entries of a fresh permutation of the off-pulse bins
        randperm!(rng, pidx)
    end
    for j in 1:non
        rot[j] = cis(2π * rand(rng))
    end
    @inbounds for b in 1:nwin
        g  = zero(ComplexF64)
        Lj = At[1, b] + rot[1] * L_offt[pidx[1], b]
        for j in 2:non
            Lj1 = At[j, b] + rot[j] * L_offt[pidx[j], b]
            g += conj(Lj) * Lj1
            Lj = Lj1
        end
        im_g[b] = imag(g)
    end
    return im_g
end


"""
    drift_test_sliding(data, p3, bin_st, bin_end;
                       window, stride, nreal, seed, sig_min) -> NamedTuple

Reversal-tolerant, time-resolved drift detector — the sliding-window
counterpart of `drift_test`, built for drifters that change drift direction
(e.g. J1750-3503: negative-drift episodes of 28±4 P, positive of 88±15 P,
Szary+2022). For such pulsars the global coherent slope of `drift_test`
mixes episodes of opposite phase gradient and averages to ~zero — "no
drift" at any significance, no matter how strong the drift is — and its
single FFT bin loses most of the feature to P3 wobble on top.

Per window position b (length `window`, step `stride`):

  L_b(φ) = Σ_{n∈b} (data[n,φ] - running mean) · e^(+2πi f3 n)
  g_b    = Σ_φ conj(L_b[φ]) · L_b[φ+1]        (per-window `slope_stat`)

and the detection statistic combines the drift quadratures incoherently:

  S = Σ_b |Im g_b|

Why Im g: any longitude-stationary modulation — pure amplitude modulation,
but also short-window leakage from nulling and slow intensity changes —
contributes c_b·B(φ) with B real, which lands entirely in Re g_b (a common
per-window phase cancels in the conjugate products). A phase *gradient* is
exactly what rotates power into Im g_b, and |·| makes S invariant under
drift-sign reversals between windows: episodes of + and − drift add
instead of cancelling. arg g_b is the usual slope [rad/bin], returned as
slope(t) — for a reversing drifter it flips sign between episodes, which
is the direct measurement `phase_modulation`/`p3fold_coherent` never show.

Null ("same modulation power, no drift"): per-window measured amplitudes
with the phase flattened — so P3 wobble, nulls and envelope evolution are
all kept — plus noise bootstrapped from the pulsar's own off-pulse
windowed LRFS (`null_im_g!`), which shares one noise field across
overlapping windows exactly as the data does. That sharing is essential
for stride < window: overlapping windows have positively correlated
|Im g_b|, which widens the null spread of S — independent per-window noise
would narrow it and inflate the significance. The amplitudes are
noise-debiased per bin (A_b = sqrt(max(|L_b|² − 2σ_b², 0))) and then
rescaled per window so the *total* signal power matches Σ(|L_b|² − 2σ_b²):
without that rescale the truncation at 0 leaves phantom amplitude in weak
bins, and because mean(S)/std(S) is large, a few-percent power mismatch
shifts the significance by several σ.

Preprocessing: the static profile and slow intensity fluctuations are
removed by `subtract_running_mean` with half-width `hp_halfwin` (default
round(p3), which nulls the filter response exactly at f3). Slow power
leaking into short windows is longitude-stationary and flat-phase, so it
cannot fake drift (it lands in Re g), but it dilutes the *measured* slope
arg(g_b) toward zero — the high-pass keeps slope(t) honest.

Caveat: S detects any longitude-dependence of the modulation phase. Two
longitude-separated components modulated with independent temporal phases
(a phase step, not a drift) also put power into Im g; `drift_test`'s
chi2_red and the slope(t) panel distinguish the two cases.

Choosing `window`: about half the shortest expected drift episode. Upper
limit — a window longer than an episode never fits inside one and
partially self-cancels (the global-slope pathology in miniature). Lower
limit — per-window feature SNR scales as √window and the incoherent sum
loses power steeply once it drops below ~1 (the |·| folding buries weak
signal in the Rayleigh floor); a warning is issued when the median window
SNR is < 1.5 (pure noise gives ~1.25). If you scan `window` (8/16/32/64/128
is a useful diagnostic), quote the significance at the pre-chosen default
or apply a trials correction — the maximum of a scan is biased high.

Arguments:
  data     – single-pulse matrix (N_pulses × N_bins), real intensity
  p3       – nominal P3 [P0]; sets f3 = 1/p3 (aliased P3 < 2 not supported)
  bin_st   – first on-pulse bin (1-indexed)
  bin_end  – last on-pulse bin (1-indexed)
  window   – sliding window length [pulses] (default 16)
  stride   – window step [pulses]; 1 = sliding (default), window = disjoint blocks
  nreal    – number of null realisations (default 1000)
  seed     – RNG seed (default 7, nothing for non-reproducible)
  sig_min  – per-window |Im g_b|/σ threshold for the `detected` mask (default 3)
  hp_halfwin – half-width of the running-mean high-pass [pulses];
             nothing (default) → round(p3), 0 → plain mean subtraction
  amp_smooth – half-width [bins] of the longitude boxcar smoothing the null
             amplitude profiles (default 5); the null is driven by the
             amplitude gradient, so raw noisy |L| would inflate it. With the
             defaults the null is deliberately left slightly conservative
             (H0 reads ≈ −1σ in synthetic calibration, never positive)

Fields of returned NamedTuple:
  p3, f3, window, stride – inputs echoed back
  on_bins      – on-pulse bin UnitRange
  starts       – first pulse of each window
  centers      – window centre pulse numbers (x-axis of the time panels)
  L_on         – windowed complex LRFS, nwin × N_on
  Lmap         – |L_on| in units of the per-window noise σ_b (nwin × N_on)
  sigma_win    – per-window off-pulse noise σ_b in L
  snr_win      – per-window feature SNR, mean(|L_b|)/σ_b; noise gives ~1.25
  snr_med      – median of snr_win
  g            – per-window coherent statistic g_b (complex)
  slope        – arg g_b [rad/bin], the time-resolved drift slope; > 0 =
                 positive drift (early → later longitudes, Szary+2022)
  slope_err    – 1σ error on slope from the null Im-spread, min-capped at π
  slope_sig    – |Im g_b| / σ_Im,b, per-window drift significance
  detected     – slope_sig .≥ sig_min
  S            – observed Σ|Im g_b|
  S_null       – null distribution of S [nreal]
  significance – (S − mean(S_null)) / std(S_null)  [σ]
  p_value      – empirical p, count(S_null ≥ S)/nreal (0 means < 1/nreal)
"""
function drift_test_sliding(data::AbstractMatrix, p3::Real, bin_st::Int, bin_end::Int;
                            window::Int=16, stride::Int=1, nreal::Int=1000,
                            seed::Union{Int,Nothing}=7, sig_min::Real=3.0,
                            hp_halfwin::Union{Int,Nothing}=nothing, amp_smooth::Int=5)
    p3 > 2 || error("p3 must be > 2 P0 (got $p3): f3 = 1/p3 would exceed Nyquist")
    nbins = size(data, 2)
    on  = bin_st:bin_end
    non = length(on)
    non < 3 && error("Need at least 3 on-pulse bins (got $non)")
    f3 = 1.0 / p3

    # high-pass: kill DC + slow fluctuations, filter null sits at f3 (see docstrings)
    hp = isnothing(hp_halfwin) ? round(Int, p3) : hp_halfwin
    data_dm = subtract_running_mean(data, hp)
    wl = windowed_lrfs(data_dm, f3, window, stride)
    L, starts = wl.L, wl.starts
    nwin    = length(starts)
    centers = starts .+ (window - 1) / 2

    off = vcat(1:bin_st-1, bin_end+1:nbins)
    isempty(off) && error(
        "No off-pulse bins available (bin_st=$bin_st bin_end=$bin_end N_bins=$nbins)")
    sigma_win = [std([real.(@view L[b, off]); imag.(@view L[b, off])]) for b in 1:nwin]
    L_on = L[:, on]
    snr_win = [mean(abs.(@view L_on[b, :])) / sigma_win[b] for b in 1:nwin]
    snr_med = median(snr_win)
    snr_med < 1.5 && @warn "Weak per-window f3 feature (median SNR = " *
        "$(round(snr_med, digits=2)), pure noise gives ~1.25): detection power " *
        "collapses below SNR ~1 — consider a longer `window`"

    g = Vector{ComplexF64}(undef, nwin)
    @inbounds for b in 1:nwin
        acc = zero(ComplexF64)
        for j in 1:non-1
            acc += conj(L_on[b, j]) * L_on[b, j+1]
        end
        g[b] = acc
    end
    slope = angle.(g)
    S = sum(abs, imag.(g))

    At = null_amplitudes(L_on, sigma_win; amp_smooth=amp_smooth)
    noff = length(off)
    replace = noff < non
    replace && @warn "Fewer off-pulse than on-pulse bins ($noff < $non): " *
        "null noise bootstrap samples with replacement"
    L_offt = collect(permutedims(L[:, off]))
    rng     = isnothing(seed) ? Random.default_rng() : MersenneTwister(seed)
    pidx    = collect(1:(replace ? non : noff))
    rot     = Vector{ComplexF64}(undef, non)
    S_null  = zeros(nreal)
    im_sum  = zeros(nwin)
    im_sum2 = zeros(nwin)
    im_g    = zeros(nwin)
    for i in 1:nreal
        null_im_g!(im_g, At, L_offt, pidx, rot, rng, replace)
        S_null[i] = sum(abs, im_g)
        im_sum  .+= im_g
        im_sum2 .+= im_g .^ 2
    end
    sigma_im = sqrt.(max.(im_sum2 ./ nreal .- (im_sum ./ nreal) .^ 2, eps()))

    significance = (S - mean(S_null)) / std(S_null)
    p_value      = count(>=(S), S_null) / nreal
    slope_sig    = abs.(imag.(g)) ./ sigma_im
    slope_err    = min.(sigma_im ./ abs.(g), Float64(π))
    detected     = slope_sig .>= sig_min

    return (
        p3           = p3,
        f3           = f3,
        window       = window,
        stride       = stride,
        on_bins      = on,
        starts       = starts,
        centers      = centers,
        L_on         = L_on,
        Lmap         = abs.(L_on) ./ sigma_win,
        sigma_win    = sigma_win,
        snr_win      = snr_win,
        snr_med      = snr_med,
        g            = g,
        slope        = slope,
        slope_err    = slope_err,
        slope_sig    = slope_sig,
        detected     = detected,
        S            = S,
        S_null       = S_null,
        significance = significance,
        p_value      = p_value,
    )
end

"""
    null_G!(G, At, L_offt, pidx, rot, rng, replace) -> G

One flat-phase surrogate realisation for `drift_test_profile`: fills `G`
with G_j = Σ_b conj(L_b[j])·L_b[j+1] for the surrogate windowed LRFS
L_b^null(φ_j) = A_b(φ_j) + N_b(φ_j).

Same construction as `null_im_g!` — measured amplitudes `At` from
`null_amplitudes`, noise bootstrapped from the off-pulse windowed LRFS
`L_offt` (one random off-pulse bin plus one random phase rotation per
on-pulse bin, held fixed across windows) — only the reduction differs:
here the products are summed over *windows* at fixed longitude instead of
over longitude at fixed window. Adjacent G_j share the noise of bin j+1
exactly as in the data, so the correlation between neighbouring phase
increments is reproduced without being modelled.
"""
function null_G!(G::Vector{ComplexF64}, At::Matrix{Float64},
                 L_offt::Matrix{ComplexF64}, pidx::Vector{Int},
                 rot::Vector{ComplexF64}, rng, replace::Bool)
    non, nwin = size(At)
    noff = size(L_offt, 1)
    if replace
        rand!(rng, pidx, 1:noff)
    else
        randperm!(rng, pidx)
    end
    for j in 1:non
        rot[j] = cis(2π * rand(rng))
    end
    fill!(G, zero(ComplexF64))
    @inbounds for b in 1:nwin
        Lj = At[1, b] + rot[1] * L_offt[pidx[1], b]
        for j in 2:non
            Lj1 = At[j, b] + rot[j] * L_offt[pidx[j], b]
            G[j-1] += conj(Lj) * Lj1
            Lj = Lj1
        end
    end
    return G
end


"""
    profile_stat(G) -> (slope, resid, T)

Constant-gradient fit to the longitude-resolved phase increments of
`drift_test_profile`. The fitted gradient is the coherent slope
`slope = arg(Σ_j G_j)` (identical estimator to `slope_stat`), and the
residual of bin j is the component of G_j perpendicular to it,

  resid_j = Im( G_j · e^(-i·slope) ),     T = Σ_j resid_j²

`resid_j` vanishes whenever arg(G_j) equals the fitted slope, whatever the
amplitude |G_j|, so T measures *only* departure from a constant gradient.
Its natural weighting is also the correct one: the phase error of arg(G_j)
scales as 1/|G_j|, so Im(G_j e^(-i·slope)) carries noise of roughly the
same size in every bin regardless of how strong that bin is, and an
unweighted sum of squares is already the inverse-variance form.
"""
function profile_stat(G::AbstractVector)
    slope = angle(sum(G))
    resid = [imag(G[j] * cis(-slope)) for j in eachindex(G)]
    return (slope=slope, resid=resid, T=sum(abs2, resid))
end


"""
    drift_test_profile(data, p3, bin_st, bin_end; window, stride, nreal, seed,
                       hp_halfwin, amp_smooth, pulse_st, pulse_end) -> NamedTuple

Longitude-resolved phase gradient from *windowed* LRFS — the test that
separates a genuine subpulse drift from a phase step between two
longitude-separated components, for pulsars where `drift_test`'s global
single-bin LRFS carries no usable signal.

`drift_test_sliding` reports Σ_φ of the pairwise products, which is a
single number per window: it detects that the modulation phase depends on
longitude, but not *how*. A steady drift (ψ rising across the profile) and
a phase *step* inside the profile (ψ flat on either side, jumping once)
both put power into that sum, and both read as a significant "drift" with
a plausible slope. Only the longitude-resolved increments tell them apart,
which is what `drift_test`/`slope_chi2` does — but from the global FFT
bin, and that bin is empty whenever P3 wobbles: the feature smears over
Δk ≈ k·ΔP3/P3 bins with k = N/P3, so a short P3 in a long observation
(J2053-7200: k ≈ 340, ±1.3% wobble ⇒ ~9 bins) leaves nothing to measure.

The fix uses the same pairwise products as `drift_test_sliding`, summed
along the other axis. With P[b,j] = conj(L_b[j])·L_b[j+1]:

  Σ_j P[b,j] = g_b  → slope(t), time-resolved   (`drift_test_sliding`)
  Σ_b P[b,j] = G_j  → dψ/dφ(φ), longitude-resolved   (here)

Summing over windows is legitimate despite P3 wobble because P[b,j] is
invariant to the window's absolute phase: an unknown e^(iθ_b) cancels in
the conjugate product. That is the same invariance the sliding statistic
relies on, so the windowed LRFS keeps the signal that the global FFT bin
loses, while short windows keep their wide response (f3 ± 1/window).

Goodness of fit is Monte-Carlo calibrated rather than analytic: `T` from
`profile_stat` is recomputed on flat-phase surrogates built from the
measured amplitudes plus off-pulse-bootstrap noise (`null_G!`), which
reproduces the correlation between neighbouring increments — they share a
bin — without a covariance model. Reported as

  chi2_red = T / median(T_null)      ≈ 1  constant gradient fits → drift
                                     ≫ 1  phase jumps → not a drift

with an empirical p-value alongside. This is a calibrated ratio, not a
formal reduced χ²; it is built to read like `drift_test`'s `chi2_red` so
the two can be compared directly.

Synthetic calibration (scratchpad tests, W=32) of what it can and cannot
separate:

  * a sharp phase step inside one broad component — the case this exists
    for — is caught decisively: S reads a spurious 25σ "drift" at a
    plausible -2.2°/bin, while chi2_red = 20 with p = 0 and the largest
    residuals land exactly on the step bin;
  * a genuine steady drift reads chi2_red = 0.5-0.9 at S = 33-97σ, i.e.
    the test stays quiet where it should;
  * two *well-separated* components with a phase offset need no test at
    all: the step falls in the low-amplitude gap between them, so it
    barely contributes to the products and S does not detect it either;
  * two *heavily overlapping* components with a phase offset (separation
    ≲ 2 component widths) turn the step into a smooth ramp and read
    chi2_red ≈ 1-1.8. That is a real physical degeneracy rather than a
    failure — across an overlap region "two components offset in
    modulation phase" and "drift" describe the same observable — and it
    cannot be resolved from the f3 phase alone.

chi2_red is only interpretable when the modulation itself is detected:
check `snr_med` here and the significance from `drift_test_sliding` first,
exactly as `drift_test`'s χ² is meaningless once its `snr` approaches the
pure-noise value.

Caveat — the window sum is coherent, so it assumes the drift keeps **one
sense** through the analysed stretch. For a reversing drifter
(J1750-3503) the episodes would cancel exactly as they do in `drift_test`;
restrict the analysis to a single episode with `pulse_st`/`pulse_end`,
which are applied before everything else (high-pass included).

Arguments: as `drift_test_sliding`, plus
  pulse_st, pulse_end – analyse only this pulse range (default: all)

Fields of the returned NamedTuple:
  p3, f3, window, stride, on_bins – inputs echoed back
  nwin        – number of windows summed over
  pulse_range – (first, last) pulse actually analysed
  G           – Σ_b conj(L_b[j])·L_b[j+1], length n_on-1
  dpsi        – arg G_j [rad/bin], the longitude-resolved phase gradient
  dpsi_sigma  – 1σ on dpsi from the surrogates [rad]
  psi_cum     – cumsum(dpsi) [rad]: ψ(φ) rebuilt from the increments —
                a straight line under drift, a staircase under phase steps
  slope       – fitted constant gradient [rad/bin] (> 0 = positive drift)
  resid       – per-bin departure from that gradient (`profile_stat`)
  amp         – mean_b |L_b(φ)| / σ_b, the on-pulse S/N profile
  snr_med     – median per-window feature S/N (as in `drift_test_sliding`)
  T, T_null   – fit statistic and its surrogate distribution
  chi2_red    – T / median(T_null)
  p_value     – count(T_null ≥ T)/nreal (0 means < 1/nreal)
"""
function drift_test_profile(data::AbstractMatrix, p3::Real, bin_st::Int, bin_end::Int;
                            window::Int=32, stride::Int=1, nreal::Int=500,
                            seed::Union{Int,Nothing}=7,
                            hp_halfwin::Union{Int,Nothing}=nothing, amp_smooth::Int=5,
                            pulse_st::Union{Int,Nothing}=nothing,
                            pulse_end::Union{Int,Nothing}=nothing)
    p3 > 2 || error("p3 must be > 2 P0 (got $p3): f3 = 1/p3 would exceed Nyquist")
    nbins = size(data, 2)
    on    = bin_st:bin_end
    non   = length(on)
    non < 3 && error("Need at least 3 on-pulse bins (got $non)")
    f3 = 1.0 / p3

    n_st  = isnothing(pulse_st) ? 1 : pulse_st
    n_end = isnothing(pulse_end) ? size(data, 1) : pulse_end
    (1 <= n_st < n_end <= size(data, 1)) ||
        error("Bad pulse range $n_st:$n_end for $(size(data, 1)) pulses")
    dsub = @view data[n_st:n_end, :]

    hp      = isnothing(hp_halfwin) ? round(Int, p3) : hp_halfwin
    data_dm = subtract_running_mean(dsub, hp)
    wl      = windowed_lrfs(data_dm, f3, window, stride)
    L       = wl.L
    nwin    = size(L, 1)

    off = vcat(1:bin_st-1, bin_end+1:nbins)
    isempty(off) && error(
        "No off-pulse bins available (bin_st=$bin_st bin_end=$bin_end N_bins=$nbins)")
    sigma_win = [std([real.(@view L[b, off]); imag.(@view L[b, off])]) for b in 1:nwin]
    L_on      = L[:, on]
    snr_med   = median([mean(abs.(@view L_on[b, :])) / sigma_win[b] for b in 1:nwin])
    snr_med < 1.5 && @warn "Weak per-window f3 feature (median SNR = " *
        "$(round(snr_med, digits=2)), pure noise gives ~1.25): consider a longer `window`"

    # G_j = Σ_b conj(L_b[j])·L_b[j+1] — the pairwise products of
    # `drift_test_sliding` summed over windows instead of over longitude
    G = zeros(ComplexF64, non - 1)
    @inbounds for b in 1:nwin
        Lj = L_on[b, 1]
        for j in 2:non
            Lj1 = L_on[b, j]
            G[j-1] += conj(Lj) * Lj1
            Lj = Lj1
        end
    end
    st = profile_stat(G)

    At      = null_amplitudes(L_on, sigma_win; amp_smooth=amp_smooth)
    noff    = length(off)
    replace = noff < non
    replace && @warn "Fewer off-pulse than on-pulse bins ($noff < $non): " *
        "null noise bootstrap samples with replacement"
    L_offt = collect(permutedims(L[:, off]))
    rng    = isnothing(seed) ? Random.default_rng() : MersenneTwister(seed)
    pidx   = collect(1:(replace ? non : noff))
    rot    = Vector{ComplexF64}(undef, non)
    Gn     = zeros(ComplexF64, non - 1)
    T_null = zeros(nreal)
    a_sum  = zeros(non - 1)
    a_sum2 = zeros(non - 1)
    for i in 1:nreal
        null_G!(Gn, At, L_offt, pidx, rot, rng, replace)
        sn = profile_stat(Gn)
        T_null[i] = sn.T
        for j in eachindex(Gn)
            d = wrap_pi(angle(Gn[j]) - sn.slope)
            a_sum[j]  += d
            a_sum2[j] += d^2
        end
    end
    dpsi_sigma = sqrt.(max.(a_sum2 ./ nreal .- (a_sum ./ nreal) .^ 2, eps()))

    med_null = median(T_null)
    dpsi     = angle.(G)
    amp      = [mean(abs(L_on[b, j]) / sigma_win[b] for b in 1:nwin) for j in 1:non]

    return (
        p3          = p3,
        f3          = f3,
        window      = window,
        stride      = stride,
        on_bins     = on,
        nwin        = nwin,
        pulse_range = (n_st, n_end),
        G           = G,
        dpsi        = dpsi,
        dpsi_sigma  = dpsi_sigma,
        psi_cum     = cumsum(dpsi),
        slope       = st.slope,
        resid       = st.resid,
        amp         = amp,
        snr_med     = snr_med,
        T           = st.T,
        T_null      = T_null,
        chi2_red    = med_null > 0 ? st.T / med_null : NaN,
        p_value     = count(>=(st.T), T_null) / nreal,
    )
end

end  # module PhaseDrift
