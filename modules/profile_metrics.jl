module ProfileMetrics

using FFTW
using Statistics
using Printf

export global_fourier_shift, component_barycenters, component_barycenters_fixed_windows

"""
Calculates the global sub-bin phase offset between two profiles using the 
Fourier Phase Gradient Method (based on Taylor 1992 pulsar timing algorithm).
This methodology is entirely non-parametric, shape-independent, and highly accurate.

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
"""
function component_barycenters(prof_ref::AbstractVector, prof_target::AbstractVector; threshold=0.15, n_comp=3)
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
    
    # 1. Identify peaks using local maxima detection
    peaks = Int[]
    for i in 2:N-1
        if p_ref[i] > p_ref[i-1] && p_ref[i] > p_ref[i+1] && p_ref[i] > threshold
            push!(peaks, i)
        end
    end
    
    # Process up to n_comp strongest components
    sort!(peaks, by=p -> p_ref[p], rev=true)
    selected_peaks = sort(peaks[1:min(length(peaks), n_comp)]) 
    
    results = []
    for p in selected_peaks
        # 2. Define the Component Window Boundaries
        # Stop scanning outwards when we hit a local minimum or signal drops below 30% of peak
        cut_off = p_ref[p] * 0.3
        
        left = p
        while left > 1 && p_ref[left] > cut_off && p_ref[left-1] <= p_ref[left]
            left -= 1
        end
        
        right = p
        while right < N && p_ref[right] > cut_off && p_ref[right+1] <= p_ref[right]
            right += 1
        end
        
        if right - left < 2; continue; end
        
        # 3. Calculate Barycenters using the extracted window
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
        bg_ref = min(p_ref[left], p_ref[right])
        bg_tar = min(p_tar[left], p_tar[right])
        vals_ref = max.(p_ref[x_win] .- bg_ref, 0.0)
        vals_tar = max.(p_tar[x_win] .- bg_tar, 0.0)
        m_ref = sum(vals_ref)
        m_tar = sum(vals_tar)
        com_ref = m_ref > 0 ? sum(x_win .* vals_ref) / m_ref : NaN
        com_tar = m_tar > 0 ? sum(x_win .* vals_tar) / m_tar : NaN
        return (com_ref=com_ref, offset=com_tar - com_ref)
    end
end

end # module