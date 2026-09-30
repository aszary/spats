module LrfsDiagnostics

using FFTW
using Statistics
using LinearAlgebra

export lrfs_phase_track

# Include the existing tools module
try
    import ..Tools
catch
    try
        include("tools.jl")
        import .Tools
    catch
        @warn "Could not load tools.jl for Tools.lrfs"
    end
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
    # 1. Preprocess: isolate on-pulse and remove static profile
    X_on = view(X, :, bin_st:bin_end)
    X_prep = X_on .- mean(X_on, dims=1)
    
    N, M = size(X_prep)
    
    # 2. Compute LRFS using the preexisting Tools.lrfs
    F = Tools.lrfs(X_prep)
    
    # 3. Find the dominant P3 frequency
    # Sum power across all longitudes for each frequency row
    power_spectrum = vec(sum(abs2.(F), dims=2))
    
    # Ignore DC component (index 1) and find the peak
    peak_idx = argmax(view(power_spectrum, 2:length(power_spectrum))) + 1
    
    # Calculate actual P3 value in units of pulses
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
