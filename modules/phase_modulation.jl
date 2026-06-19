module PhaseDrift

using FFTW
using Statistics
using Random

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
    amp_null_slopes(L_on, sigma_off; nreal, seed) -> Vector{Float64}

Amplitude-modulation null distribution.

Model: preserve |L_on|, set phase to zero (flat), add complex noise N(0, σ_off).
Returns nreal null slope values [rad/bin].
"""
function amp_null_slopes(L_on::AbstractVector, sigma_off::Real;
                         nreal::Int=6000, seed::Union{Int,Nothing}=7)
    rng = isnothing(seed) ? Random.default_rng() : MersenneTwister(seed)
    A = abs.(L_on)
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
  slope          – measured phase slope [rad/bin]
  null           – null slope distribution [nreal]
  significance   – slope / std(null)  [σ]
  p3             – P3 used
  p3_error       – P3 uncertainty used
  p3err_k        – k index used by each densely-sampled P3 sweep point
  p3err_p3       – the swept P3 values themselves
  p3err_psi      – phase array from the P3 sweep [rad], (length(p3err_k) × N_on)
  p3err_sigma    – per-bin circular std across the P3 sweep [rad]
"""
function drift_test(data::AbstractMatrix, p3::Real, bin_st::Int, bin_end::Int;
                    p3_error::Real=0.0, nreal::Int=6000, seed::Union{Int,Nothing}=7)
    on_bins   = bin_st:bin_end
    L         = complex_lrfs_at_f3(data, p3)
    sigma_off = off_pulse_sigma(L, bin_st, bin_end)
    L_on      = L[on_bins]
    slope     = slope_stat(L_on)
    null      = amp_null_slopes(L_on, sigma_off; nreal=nreal, seed=seed)
    significance = slope / std(null)
    p3err = p3_error_phases(data, p3, p3_error, bin_st, bin_end)
    return (
        L            = L,
        on_bins      = on_bins,
        L_on         = L_on,
        sigma_off    = sigma_off,
        slope        = slope,
        null         = null,
        significance = significance,
        p3           = p3,
        p3_error     = p3_error,
        p3err_k      = p3err.k_grid,
        p3err_p3     = p3err.p3_grid,
        p3err_psi    = p3err.psi_grid,
        p3err_sigma  = p3err.psi_sigma,
    )
end

end  # module PhaseDrift
