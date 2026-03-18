"""
    FFTCrossCorr

Module for measuring phase offsets between pulsar profiles at two frequencies
using FFT-based cross-correlation methods.

Methods available:
  - `xcorr_fft`        : standard FFT cross-correlation with parabolic sub-bin refinement
  - `xcorr_phase_slope`: phase-slope method (most noise-robust)
  - `xcorr_all`        : run all methods and return a summary

Units: bins, degrees, milliseconds.
"""
module FFT

using FFTW
using Statistics: mean, std
using Printf

export xcorr_fft, xcorr_phase_slope, xcorr_all, bootstrap_uncertainty, print_offset_summary

# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

"""Wrap offset to the interval (-N/2, N/2]."""
_wrap(offset, N) = offset > N / 2 ? offset - N : (offset <= -N / 2 ? offset + N : offset)

"""Convert bins → degrees and milliseconds given profile length and pulsar period."""
function _convert(offset_bins, nbin, period_s)
    deg = offset_bins / nbin * 360.0
    ms  = offset_bins / nbin * period_s * 1000.0
    return deg, ms
end

"""
Parabolic interpolation around the peak of vector `h` at index `idx` (1-based).
Returns sub-bin correction δ ∈ (-0.5, 0.5).
"""
function _parabolic_interp(h, idx)
    N  = length(h)
    y1 = h[mod1(idx - 1, N)]
    y2 = h[idx]
    y3 = h[mod1(idx + 1, N)]
    denom = y1 - 2y2 + y3
    return iszero(denom) ? 0.0 : 0.5 * (y1 - y3) / denom
end

# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────

"""
    xcorr_fft(prof_high, prof_low; phase_corr=false) -> offset_bins

Compute the offset between two pulsar profiles using FFT cross-correlation.

# Arguments
- `prof_high` : intensity profile at higher frequency (length N)
- `prof_low`  : intensity profile at lower frequency (length N)
- `phase_corr`: if `true`, use phase correlation (normalise cross-spectrum by
                amplitude) to produce a sharper, noise-suppressed peak

# Returns
Offset in bins (positive = `prof_high` leads `prof_low`).

# Method
1. Zero-mean both profiles.
2. Compute cross-power spectrum  H = conj(FFT(f)) * FFT(g).
3. Optionally normalise by |H| (phase correlation).
4. IFFT → correlation function h.
5. Locate peak; refine to sub-bin precision via parabolic interpolation.
6. Wrap result to (-N/2, N/2].
"""
function xcorr_fft(prof_high::AbstractVector, prof_low::AbstractVector;
                   phase_corr::Bool=false)
    N = length(prof_high)
    length(prof_low) == N || throw(ArgumentError("Profiles must have equal length."))

    f = prof_high .- mean(prof_high)
    g = prof_low  .- mean(prof_low)

    F = fft(f)
    G = fft(g)
    H = conj.(F) .* G

    if phase_corr
        H = H ./ (abs.(H) .+ 1e-10)   # normalise → sharper peak
    end

    h   = real.(ifft(H))
    idx = argmax(h)
    δ   = _parabolic_interp(h, idx)

    offset = _wrap(Float64(idx - 1) + δ, N)
    return offset
end


"""
    xcorr_phase_slope(prof_high, prof_low; max_freq=nothing, snr_weight=true)
        -> offset_bins

Estimate the offset by fitting a line to the phase of the cross-power spectrum.

A shift τ (in bins) produces a linear phase ramp:
    φ(k) = 2π k τ / N

so τ = slope × N / (2π).

# Arguments
- `max_freq`   : highest Fourier harmonic to include (default: N÷4).
                 Exclude high-frequency noise by keeping this low.
- `snr_weight` : weight each harmonic by its amplitude |H(k)| before
                 fitting — down-weights noisy harmonics automatically.

# Returns
Offset in bins (positive = `prof_high` leads `prof_low`).
"""
function xcorr_phase_slope(prof_high::AbstractVector, prof_low::AbstractVector;
                            max_freq::Union{Int,Nothing}=nothing,
                            snr_weight::Bool=true)
    N = length(prof_high)
    length(prof_low) == N || throw(ArgumentError("Profiles must have equal length."))

    kmax = isnothing(max_freq) ? N ÷ 4 : min(max_freq, N ÷ 2 - 1)

    f = prof_high .- mean(prof_high)
    g = prof_low  .- mean(prof_low)

    H      = conj.(fft(f)) .* fft(g)
    ks     = Float64.(1:kmax)                  # skip DC (k=0)
    phases = angle.(H[2:kmax+1])               # unwrapped not needed for slope fit
    phases = unwrap(phases)                    # unwrap to remove 2π jumps
    amps   = abs.(H[2:kmax+1])

    w = snr_weight ? amps : ones(kmax)

    # Weighted least-squares: φ = slope * k  (no intercept — DC already skipped)
    slope = sum(w .* ks .* phases) / sum(w .* ks .^ 2)

    offset = _wrap(slope * N / (2π), N)
    return offset
end

"""Unwrap a phase vector (remove 2π jumps)."""
function unwrap(phases::AbstractVector)
    out = copy(phases)
    for i in 2:length(out)
        diff = out[i] - out[i-1]
        if diff >  π; out[i] -= 2π * round(diff / (2π)); end
        if diff < -π; out[i] += 2π * round(-diff / (2π)); end
    end
    return out
end


"""
    xcorr_all(prof_high, prof_low; nbin=nothing, period_s=1.0, kw...)
        -> NamedTuple

Run all three cross-correlation variants and return a unified result.

# Keyword arguments
- `nbin`     : number of bins (defaults to `length(prof_high)`)
- `period_s` : pulsar period in seconds (for ms conversion)
- `max_freq` : passed to `xcorr_phase_slope`
- `snr_weight`: passed to `xcorr_phase_slope`

# Returns
NamedTuple with fields:
  `standard`, `phase_corr`, `phase_slope` — each a NamedTuple
  `(offset_bins, offset_deg, offset_ms)`.
"""
function xcorr_all(prof_high::AbstractVector, prof_low::AbstractVector;
                   nbin::Union{Int,Nothing}=nothing,
                   period_s::Float64=1.0,
                   max_freq::Union{Int,Nothing}=nothing,
                   snr_weight::Bool=true)

    N = isnothing(nbin) ? length(prof_high) : nbin

    function result(ob)
        deg, ms = _convert(ob, N, period_s)
        (offset_bins=ob, offset_deg=deg, offset_ms=ms)
    end

    std_offset  = xcorr_fft(prof_high, prof_low; phase_corr=false)
    pc_offset   = xcorr_fft(prof_high, prof_low; phase_corr=true)
    ps_offset   = xcorr_phase_slope(prof_high, prof_low;
                                    max_freq=max_freq, snr_weight=snr_weight)

    return (standard    = result(std_offset),
            phase_corr  = result(pc_offset),
            phase_slope = result(ps_offset))
end


"""
    bootstrap_uncertainty(prof_high, prof_low, method=:standard;
                          n_boot=1000, nbin=nothing, period_s=1.0, kw...)
        -> (offset_bins, σ_bins, σ_deg, σ_ms)

Estimate offset uncertainty via bootstrap resampling.

Adds Gaussian noise to both profiles on each iteration (noise estimated from
the off-pulse RMS), recomputes the offset, and returns the standard deviation.

# Arguments
- `method` : `:standard`, `:phase_corr`, or `:phase_slope`
- `n_boot` : number of bootstrap iterations (default 1000)
"""



function bootstrap_uncertainty(prof_high::AbstractVector, prof_low::AbstractVector,
                                method::Symbol=:standard;
                                n_boot::Int=100,
                                nbin::Union{Int,Nothing}=nothing,
                                period_s::Float64=1.0,
                                kw...)

    N    = length(prof_high)
    Nb   = isnothing(nbin) ? N : nbin

    # Estimate noise from lowest 20% of bins (proxy for off-pulse)
    noise_h = std(sort(prof_high)[1:max(1, N÷5)])
    noise_l = std(sort(prof_low )[1:max(1, N÷5)])

    offsets = Vector{Float64}(undef, n_boot)
    for i in 1:n_boot
        # Random noise bootstrap method
        ph = prof_high .+ noise_h .* randn(N)
        pl = prof_low  .+ noise_l .* randn(N)
        
        # Method
        offsets[i] = if method == :standard
            xcorr_fft(ph, pl; phase_corr=false)
        elseif method == :phase_corr
            xcorr_fft(ph, pl; phase_corr=true)
        elseif method == :phase_slope
            xcorr_phase_slope(ph, pl; kw...)
        else
            error("Unknown method: $method")
        end
    end

    offset_mean = mean(offsets)
    # mean offset to deg
    off_deg, off_ms = _convert(offset_mean, Nb, period_s)
    
    σ_bins = std(offsets)
    #error calculation to deg and ms
    σ_deg, σ_ms = _convert(σ_bins, Nb, period_s)

    
    return (
        offset_bins = offset_mean, 
        offset_deg  = off_deg,     
        offset_ms   = off_ms,      
        σ_bins      = σ_bins, 
        σ_deg       = σ_deg, 
        σ_ms        = σ_ms
    )
end







"""
    print_offset_summary(results, uncertainties=nothing; label="")

Pretty-print the output of `xcorr_all` and optionally `bootstrap_uncertainty`.
"""
function print_offset_summary(results::NamedTuple,
                               uncertainties=nothing;
                               label::String="")
    header = isempty(label) ? "Cross-correlation offsets" : "Cross-correlation offsets — $label"
    println("\n", "─"^60)
    println(header)
    println("─"^60)
    @printf("%-20s  %+8s  %+9s  %+9s\n", "Method", "bins", "deg", "ms")
    println("─"^60)
    for (name, r) in pairs(results)
        @printf("%-20s  %+8.3f  %+9.3f  %+9.3f\n",
                string(name), r.offset_bins, r.offset_deg, r.offset_ms)
    end
    if !isnothing(uncertainties)
        println("─"^60)
        println("Bootstrap uncertainty ($(string(uncertainties.method))):")
        @printf("  σ = %.3f bins  =  %.3f°  =  %.3f ms\n",
                uncertainties.σ_bins, uncertainties.σ_deg, uncertainties.σ_ms)
    end
    println("─"^60, "\n")
end

end # module FFTCrossCorr