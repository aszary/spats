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
    
    min_ref, max_ref = minimum(prof_ref), maximum(prof_ref)
    min_tar, max_tar = minimum(prof_target), maximum(prof_target)
    
    if (max_ref - min_ref) ≈ 0.0 || (max_tar - min_tar) ≈ 0.0
        return []
    end
    
    p_ref = (prof_ref .- min_ref) ./ (max_ref - min_ref + 1e-9)
    p_tar = (prof_target .- min_tar) ./ (max_tar - min_tar + 1e-9)
    
    # 1. Create a Topological Map (Heavy smoothing to ignore noise wiggles, ONLY for boundary logic)
    w_size = max(2, div(N, 100)) # Reduced smoothing to separate very close subcomponents
    smooth_ref = copy(p_ref)
    for i in 1:N
        start_idx = max(1, i - w_size)
        end_idx = min(N, i + w_size)
        smooth_ref[i] = mean(p_ref[start_idx:end_idx])
    end
    
    # 2. Identify structural peaks on the topological map
    raw_peaks = Int[]
    for i in 2:N-1
        if smooth_ref[i] > smooth_ref[i-1] && smooth_ref[i] > smooth_ref[i+1] && smooth_ref[i] > threshold
            push!(raw_peaks, i)
        end
    end
    
    # Filter peaks: extremely close overlaps
    min_dist = max(2, div(N, 100)) # Significantly reduced to allow very close peaks
    peaks = Int[]
    sort!(raw_peaks, by=p -> smooth_ref[p], rev=true) 
    for p in raw_peaks
        if all(abs(p - ep) > min_dist for ep in peaks)
            push!(peaks, p)
        end
    end
    
    # Enforce a hard limit (literature standard usually evaluates up to 6 distinct subcomponents)
    limit = isnothing(n_comp) ? 6 : n_comp
    max_comps = min(length(peaks), limit)
    selected_peaks = sort(peaks[1:max_comps]) 
    
    results = []
    for i in 1:length(selected_peaks)
        p = selected_peaks[i]
        
        # 3. Define Watershed Boundaries (Local Minima Valleys)
        left_limit = i > 1 ? selected_peaks[i-1] : 1
        right_limit = i < length(selected_peaks) ? selected_peaks[i+1] : N
        
        left = p
        while left > left_limit + 1 && smooth_ref[left-1] <= smooth_ref[left]
            left -= 1
        end
        
        right = p
        while right < right_limit - 1 && smooth_ref[right+1] <= smooth_ref[right]
            right += 1
        end
        
        # Safely expand to ensure we capture the tail, avoiding overlapping
        left = max(1, left - div(N, 100))
        right = min(N, right + div(N, 100))
        
        x_win = left:right
        
        bg_ref = min(p_ref[left], p_ref[right])
        m_ref = sum(max.(p_ref[x_win] .- bg_ref, 0.0))
        
        bg_tar = min(p_tar[left], p_tar[right])
        m_tar = sum(max.(p_tar[x_win] .- bg_tar, 0.0))
        
        com_ref = m_ref > 0 ? sum(x_win .* max.(p_ref[x_win] .- bg_ref, 0.0)) / m_ref : Float64(p)
        com_tar = m_tar > 0 ? sum(x_win .* max.(p_tar[x_win] .- bg_tar, 0.0)) / m_tar : Float64(p)
        
        push!(results, (peak_bin=p, left_bound=left, right_bound=right, com_low=com_ref, com_high=com_tar, offset=com_tar - com_ref))
    end
    
    return results
end

"""
Calculates high-precision sub-bin Phase-Resolved Offsets using the Windowed Cross-Correlation Function (CCF).
Based on traditional pulsar timing techniques, it cross-correlates the low and high frequency 
subpulses within the mathematically constrained window and uses a parabolic interpolation 
at the peak of the correlation function to extract shape-independent sub-bin precision shifts.
"""
function component_barycenters_fixed_windows(prof_ref::AbstractVector, prof_target::AbstractVector, windows)
    # Normalize with protection against flat profiles (zero amplitude)
    min_ref, max_ref = minimum(prof_ref), maximum(prof_ref)
    min_tar, max_tar = minimum(prof_target), maximum(prof_target)
    
    amp_ref = max_ref - min_ref
    amp_tar = max_tar - min_tar
    
    if amp_ref ≈ 0.0 || amp_tar ≈ 0.0
        return [(com_ref=NaN, offset=NaN) for _ in windows]
    end
    
    p_ref = (prof_ref .- min_ref) ./ (amp_ref + 1e-9)
    p_tar = (prof_target .- min_tar) ./ (amp_tar + 1e-9)
    
    return map(windows) do (left, right)
        x_win = left:right
        r_vals = p_ref[x_win]
        t_vals = p_tar[x_win]
        
        # Define maximum expected physical shift
        max_shift = max(1, min(length(x_win) - 1, div(length(prof_ref), 10)))
        shifts = -max_shift:max_shift
        cc_array = zeros(length(shifts))
        
        # Time-domain Cross-Correlation
        for (idx, s) in enumerate(shifts)
            cc = 0.0
            for i in 1:length(x_win)
                i_shifted = i - s
                if i_shifted >= 1 && i_shifted <= length(x_win)
                    cc += r_vals[i] * t_vals[i_shifted]
                end
            end
            cc_array[idx] = cc
        end
        
        # Find Peak and apply Parabolic Interpolation for Sub-bin Precision
        max_idx = argmax(cc_array)
        if max_idx > 1 && max_idx < length(cc_array)
            y1, y2, y3 = cc_array[max_idx-1], cc_array[max_idx], cc_array[max_idx+1]
            denom = (y1 - 2*y2 + y3)
            sub_shift = denom != 0 ? 0.5 * (y1 - y3) / denom : 0.0
            offset = shifts[max_idx] + sub_shift
        else
            offset = Float64(shifts[max_idx])
        end
        
        # Calculate structural Barycenter of reference exclusively for X-axis plotting
        bg_ref = minimum(r_vals)
        m_ref = sum(max.(r_vals .- bg_ref, 0.0))
        com_ref = m_ref > 0 ? sum(x_win .* max.(r_vals .- bg_ref, 0.0)) / m_ref : (left+right)/2.0
        
        return (com_ref=com_ref, offset=offset)
    end
end

end # module