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
"""
function complex_lrfs_at_f3(data::AbstractMatrix, p3::Real)
    N = size(data, 1)
    F = fft(data, 1)                      # (N_pulses, N_bins)
    k = round(Int, N / p3)
    k = clamp(k, 1, N ÷ 2)
    return F[k + 1, :]                    # +1: Julia 1-indexed (DC is at index 1)
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
        psi_grid[i, :] = angle.(F[k + 1, on_bins])
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
  slope          – measured phase slope [rad/bin]
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

end  # module PhaseDrift
