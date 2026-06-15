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
    drift_test(data, p3, bin_st, bin_end; nreal, seed) -> NamedTuple

Phase-drift vs amplitude-modulation discriminator.

Computes the complex LRFS at f3=1/P3, measures the coherent phase slope,
and compares it against a null distribution for pure amplitude modulation
(flat phase + off-pulse noise). Returns significance in σ units.

Arguments:
  data    – single-pulse matrix (N_pulses × N_bins), real intensity
  p3      – P3 period [pulse periods P0]
  bin_st  – first on-pulse bin (1-indexed)
  bin_end – last on-pulse bin (1-indexed)
  nreal   – number of null realisations (default 6000)
  seed    – RNG seed (default 7, nothing for non-reproducible)

Fields of returned NamedTuple:
  L            – complex LRFS slice at f3, all N_bins
  on_bins      – on-pulse bin UnitRange
  L_on         – L restricted to on-pulse
  sigma_off    – off-pulse noise estimate in complex L
  slope        – measured phase slope [rad/bin]
  null         – null slope distribution [nreal]
  significance – slope / std(null)  [σ]
  p3           – P3 used
"""
function drift_test(data::AbstractMatrix, p3::Real, bin_st::Int, bin_end::Int;
                    nreal::Int=6000, seed::Union{Int,Nothing}=7)
    on_bins   = bin_st:bin_end
    L         = complex_lrfs_at_f3(data, p3)
    sigma_off = off_pulse_sigma(L, bin_st, bin_end)
    L_on      = L[on_bins]
    slope     = slope_stat(L_on)
    null      = amp_null_slopes(L_on, sigma_off; nreal=nreal, seed=seed)
    significance = slope / std(null)
    return (
        L            = L,
        on_bins      = on_bins,
        L_on         = L_on,
        sigma_off    = sigma_off,
        slope        = slope,
        null         = null,
        significance = significance,
        p3           = p3,
    )
end

end  # module PhaseDrift
