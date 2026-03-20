module ProfileMetrics

using FFTW
using Statistics
using Printf

export global_fourier_shift, component_barycenters, component_barycenters_fixed_windows

"""
Calculates the global sub-bin phase offset between two profiles using the 
Fourier Phase Gradient Method (based on Taylor 1992 pulsar timing algorithm).
This methodology is entirely non-parametric, shape-independent, and highly accurate.

It leverages the Shift Theorem to calculate the global delay across the entire 
frequency band of the profile's envelope, weighted by spectral power to minimize 
the impact of high-frequency noise.

Returns the shift in bins. A positive shift indicates that `prof_target` is 
delayed (shifted to the right) relative to `prof_ref`.
"""
function global_fourier_shift(prof_ref::AbstractVector, prof_target::AbstractVector)
    N = length(prof_ref)
    
    # 1. Transform profiles to the frequency domain
    F_ref = rfft(prof_ref)
    F_tar = rfft(prof_target)
    
    # 2. Compute the cross-spectrum: C(f) = F_tar(f) * conj(F_ref(f))
    # According to the Shift Theorem, if f(t-dt) is the shifted signal,
    # its cross spectrum phase will form a linear ramp: phase = -2*pi*f*dt/N
    cross_spec = F_tar .* conj.(F_ref)
    
    freqs = 0:(length(F_ref)-1)
    raw_phases = angle.(cross_spec)
    amplitudes = abs.(cross_spec)
    
    # 3. Phase Unwrapping
    # We must eliminate artificial 2pi phase jumps to get a continuous linear gradient.
    phases = copy(raw_phases)
    for i in 2:length(phases)
        diff = phases[i] - phases[i-1]
        # Map the difference to the principal branch [-pi, pi)
        wrapped_diff = mod(diff + pi, 2pi) - pi
        phases[i] = phases[i-1] + wrapped_diff
    end
    
    # 4. Weighted Linear Regression
    # We fit a line (y = m*x) to the phase gradient. We use the amplitude of the 
    # cross-spectrum as weights because high-power frequencies contain the true signal,
    # while high-frequency noise has low power and should be ignored.
    sum_wxy = 0.0
    sum_wxx = 0.0
    
    for i in 2:length(freqs) # We skip the DC component (f=0)
        w = amplitudes[i]
        x = freqs[i]
        y = phases[i]
        sum_wxy += w * x * y
        sum_wxx += w * x * x
    end
    
    # Calculate the slope of the phase
    slope = sum_wxx > 0 ? (sum_wxy / sum_wxx) : 0.0
    
    # 5. Convert phase slope back to a time-domain bin shift
    # slope = -2 * pi * shift / N  ==>  shift = -slope * N / (2 * pi)
    shift_bins = -slope * N / (2 * pi)
    
    return shift_bins
end

"""
Calculates the Geometric Center of Mass (Barycenter) for individual pulse components.
This completely avoids analytical curve-fitting (like Gaussians) which fails on asymmetric 
profiles. It mathematically isolates the "island" of each component and computes its true 
balance point relative to the local noise floor.

Barycentric calculations provide a robust, non-parametric estimator for the centroid 
of flux, making no assumptions about the intrinsic profile shape (e.g. skewness due 
to scattering or intrinsic emission mechanisms).
"""
function component_barycenters(prof_ref::AbstractVector, prof_target::AbstractVector; threshold=0.10, n_comp=nothing)
    N = length(prof_ref)
    
    # Normalize copies to 0.0 - 1.0 range for robust thresholding (with protection against flat profiles)
    min_ref, max_ref = minimum(prof_ref), maximum(prof_ref)
    min_tar, max_tar = minimum(prof_target), maximum(prof_target)
    
    amp_ref = max_ref - min_ref
    amp_tar = max_tar - min_tar
    
    # If profiles are flat, return empty result
    if amp_ref ≈ 0.0 || amp_tar ≈ 0.0
        return []
    end
    
    p_ref = (prof_ref .- min_ref) ./ amp_ref
    p_tar = (prof_target .- min_tar) ./ amp_tar
    
    # --- NEW: Adaptive Smoothing ---
    # Smooth the reference profile to eliminate noise wiggles that cause false peaks
    # and premature window boundaries. This is used ONLY as a topological map.
    window_size = max(2, div(N, 40))
    smoothed_p_ref = copy(p_ref)
    for i in 1:N
        start_idx = max(1, i - window_size)
        end_idx = min(N, i + window_size)
        smoothed_p_ref[i] = mean(p_ref[start_idx:end_idx])
    end
    
    # 1. Identify peaks using local maxima detection on SMOOTHED profile
    raw_peaks = Int[]
    for i in 2:N-1
        if smoothed_p_ref[i] > smoothed_p_ref[i-1] && smoothed_p_ref[i] > smoothed_p_ref[i+1] && smoothed_p_ref[i] > threshold
            push!(raw_peaks, i)
        end
    end
    
    # Filter peaks: enforce minimum distance to avoid overlapping detections
    min_dist = max(3, div(N, 20)) # Minimum distance is ~5% of the window
    peaks = Int[]
    sort!(raw_peaks, by=p -> smoothed_p_ref[p], rev=true) # Start with highest
    for p in raw_peaks
        if all(abs(p - ep) > min_dist for ep in peaks)
            push!(peaks, p)
        end
    end
    
    max_comps = isnothing(n_comp) ? length(peaks) : min(length(peaks), n_comp)
    selected_peaks = sort(peaks[1:max_comps]) 
    
    results = []
    for p in selected_peaks
        # 2. Define the Component Window Boundaries
        # Using the SMOOTHED profile guarantees we don't stop prematurely on noise wiggles.
        cut_off = 0.02
        
        left = p
        while left > 1 && smoothed_p_ref[left] > cut_off && smoothed_p_ref[left-1] <= smoothed_p_ref[left]
            left -= 1
        end
        
        right = p
        while right < N && smoothed_p_ref[right] > cut_off && smoothed_p_ref[right+1] <= smoothed_p_ref[right]
            right += 1
        end
        
        if right - left < 2; continue; end
        
        # 3. Calculate Barycenters using the extracted window (on ORIGINAL raw data!)
        # We subtract the local background level to isolate only the component's "mass"
        x_win = left:right
        
        bg_ref = min(p_ref[left], p_ref[right])
        vals_ref = max.(p_ref[x_win] .- bg_ref, 0.0)
        mass_ref = sum(vals_ref)
        com_ref = mass_ref > 0 ? sum(x_win .* vals_ref) / mass_ref : Float64(p)
        
        bg_tar = min(p_tar[left], p_tar[right])
        vals_tar = max.(p_tar[x_win] .- bg_tar, 0.0)
        mass_tar = sum(vals_tar)
        com_tar = mass_tar > 0 ? sum(x_win .* vals_tar) / mass_tar : Float64(p)
        
        push!(results, (peak_bin=p, left_bound=left, right_bound=right, com_low=com_ref, com_high=com_tar, offset=com_tar - com_ref))
    end
    
    return results
end

"""
Calculates Center of Mass for a specific phase (row) using strictly pre-defined boundary windows.
This ensures phase-resolved tracking doesn't artificially jump around due to noise fluctuations.

It accurately measures the physical phase separation of flux between low and high 
frequency sub-components without artificially constraining them into symmetrical bell curves.
"""
function component_barycenters_fixed_windows(prof_ref::AbstractVector, prof_target::AbstractVector, windows)
    # Normalize with protection against flat profiles (zero amplitude)
    min_ref, max_ref = minimum(prof_ref), maximum(prof_ref)
    min_tar, max_tar = minimum(prof_target), maximum(prof_target)
    
    amp_ref = max_ref - min_ref
    amp_tar = max_tar - min_tar
    
    # If profile is flat, return NaN for all windows (no signal to measure)
    if amp_ref ≈ 0.0 || amp_tar ≈ 0.0
        return [(com_ref=NaN, offset=NaN) for _ in windows]
    end
    
    p_ref = (prof_ref .- min_ref) ./ amp_ref
    p_tar = (prof_target .- min_tar) ./ amp_tar
    
    return map(windows) do (left, right)
        x_win = left:right
        
        # Use the local minimum of the window as the baseline instead of edge boundaries.
        # This prevents discarding valid noisy subpulses that happen to sit on a slope.
        bg_ref = minimum(p_ref[x_win])
        bg_tar = minimum(p_tar[x_win])
        
        vals_ref = p_ref[x_win] .- bg_ref
        vals_tar = p_tar[x_win] .- bg_tar
        m_ref = sum(vals_ref)
        m_tar = sum(vals_tar)
        com_ref = m_ref > 0 ? sum(x_win .* vals_ref) / m_ref : NaN
        com_tar = m_tar > 0 ? sum(x_win .* vals_tar) / m_tar : NaN
        return (com_ref=com_ref, offset=com_tar - com_ref)
    end
end

end # module