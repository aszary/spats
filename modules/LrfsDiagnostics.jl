module LrfsDiagnostics

using FFTW
using Statistics
using LinearAlgebra

export lrfs_phase_track

"""
    _lrfs(X) -> (F::Matrix{ComplexF64}, intensity::Vector{Float64})

Self-contained LRFS computation matching Tools.lrfs logic:
- FFT along the pulse axis for each longitude bin
- Returns the complex spectrum matrix and summed power per frequency row.
"""
function _lrfs(X::AbstractMatrix{<:Real})
    da = collect(transpose(X))        # (bins × pulses)
    bins, pulse_num = size(da)
    half = floor(Int, pulse_num / 2)

    raw = fft(da, 2)[:, 1:half]       # (bins × half_freqs)
    F   = collect(transpose(raw))     # (half_freqs × bins)  — same layout as Tools.lrfs

    intensity = vec(sum(abs.(F), dims=2))  # sum across bins per freq row
    return F, intensity
end

"""
    unwrap!(phase::AbstractVector{<:Real})

Unwraps phase jumps > π to ensure a continuous phase track.
"""
function unwrap!(phase::AbstractVector{<:Real})
    for i in 2:length(phase)
        diff = phase[i] - phase[i-1]
        if diff > pi
            phase[i:end] .-= 2pi
        elseif diff < -pi
            phase[i:end] .+= 2pi
        end
    end
    return phase
end

"""
    lrfs_phase_track(X::AbstractMatrix{<:Real}, bin_st::Int, bin_end::Int)

Determines if a pulsar is P3-only or drifting by finding the dominant P3 frequency
in the LRFS and tracking its complex phase across longitude bins.
"""
function lrfs_phase_track(X::AbstractMatrix{<:Real}, bin_st::Int, bin_end::Int)
    # Clamp bounds and check valid dimensions
    bin_st = max(1, bin_st)
    bin_end = min(size(X, 2), bin_end)
    if bin_st >= bin_end || size(X, 1) < 4
        return (p3_pulses = 0.0, phase_slope = 0.0, classification = :undetermined, phase_track = Float64[])
    end

    # 1. Preprocess: isolate on-pulse and remove static profile
    X_on   = view(X, :, bin_st:bin_end)
    X_prep = X_on .- mean(X_on, dims=1)

    N, M = size(X_prep)

    # 2. Compute LRFS (self-contained, no external dependency)
    F, intensity = _lrfs(X_prep)

    if length(intensity) < 2
        return (p3_pulses = 0.0, phase_slope = 0.0, classification = :undetermined, phase_track = Float64[])
    end

    # 3. Find the dominant P3 frequency — skip DC (index 1)
    peak_idx = argmax(view(intensity, 2:length(intensity))) + 1

    # P3 in pulses
    p3_value = N / (peak_idx - 1)
    
    # 4. Extract phase track for the dominant P3 frequency
    p3_complex_row = vec(F[peak_idx, :])
    phase_track = angle.(p3_complex_row)
    unwrap!(phase_track)
    
    # 5. Calculate phase slope across longitude (weighted linear regression)
    # Weight by the amplitude of the fluctuation to ignore noisy off-pulse bins
    weights = abs.(p3_complex_row)
    sum_w = sum(weights)
    
    if sum_w == 0.0
        return (p3_pulses = p3_value, phase_slope = 0.0, classification = :undetermined, phase_track = phase_track)
    end
    
    # Weighted linear regression for the slope (m)
    x = collect(1:M)
    x_bar = sum(weights .* x) / sum_w
    y_bar = sum(weights .* phase_track) / sum_w
    
    cov_xy = sum(weights .* (x .- x_bar) .* (phase_track .- y_bar))
    var_x = sum(weights .* (x .- x_bar).^2)
    
    slope = var_x > 0 ? cov_xy / var_x : 0.0
    
    # 6. Classification based on Phase Slope
    # Slope is in radians per bin. 
    classification = :undetermined
    if abs(slope) > 0.05
        classification = :pure_drift
    else
        classification = :amplitude_modulation
    end
    
    return (
        p3_pulses = p3_value,
        phase_slope = slope,
        classification = classification,
        phase_track = phase_track
    )
end

end # module
