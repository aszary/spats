module ProfileMetrics

using FFTW
using Statistics
using Printf

export global_fourier_shift, component_barycenters, component_barycenters_fixed_windows, component_fourier_offsets_fixed_windows, dynamic_component_fourier_offsets

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
    
    # Obliczenie BŁĘDU FORMALNEGO z WLS (Weighted Least Squares)
    N_points = length(freqs) - 1
    if N_points > 1 && sum_wxx > 0
        residuals = phases .- slope .* freqs
        # Poprawna fizycznie estymacja wariancji dla modelu bez wyrazu wolnego (y = mx)
        # Wagi określają proporcje, absolutną skalę szumu odtwarzamy z reszt (Reduced Chi-Square scaling)
        s2 = sum(amplitudes[2:end] .* residuals[2:end].^2) / (N_points - 1)
        slope_err = sqrt(s2 / sum_wxx)
        shift_error = slope_err * N / (2 * pi)
        
        # Limit fizyczny błędu (np. limit rozdzielczości binów) chroniący przed overfittingiem 
        # dla bardzo małej liczby złapanych harmonicznych na idealnym szumie.
        shift_error = max(shift_error, 0.01)
    else
        shift_error = 0.0
    end
    
    # 5. Convert phase slope back to a time-domain bin shift
    # slope = -2 * pi * shift / N  ==>  shift = -slope * N / (2 * pi)
    shift_bins = -slope * N / (2 * pi)
    
    return (shift=shift_bins, error=shift_error)
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
    
    # Wygładzamy profil w celu stabilniejszej detekcji szczytów i granic (aby okna nie były zbyt małe przez szum)
    smoothed_ref = copy(p_ref)
    for i in 3:N-2
        smoothed_ref[i] = mean(p_ref[i-2:i+2])
    end
    
    # 1. Identify peaks using local maxima detection
    peaks = Int[]
    for i in 2:N-1
        if smoothed_ref[i] > smoothed_ref[i-1] && smoothed_ref[i] > smoothed_ref[i+1] && smoothed_ref[i] > threshold
            push!(peaks, i)
        end
    end
    
    # Process up to n_comp strongest components
    sort!(peaks, by=p -> smoothed_ref[p], rev=true)
    selected_peaks = sort(peaks[1:min(length(peaks), n_comp)]) 
    
    results = []
    for p in selected_peaks
        # 2. Define the Component Window Boundaries
        # Stop scanning outwards when we hit a local minimum or signal drops below the 5% noise floor.
        # This ensures we integrate the *entire* physical component, not just the tip.
        cut_off = 0.02  # Zmniejszono z 0.05, aby lepiej "chwytać" ogony komponentu
        
        left = p
        while left > 1 && smoothed_ref[left] > cut_off && smoothed_ref[left-1] <= smoothed_ref[left]
            left -= 1
        end
        
        right = p
        while right < N && smoothed_ref[right] > cut_off && smoothed_ref[right+1] <= smoothed_ref[right]
            right += 1
        end
        
        # Rozszerzamy sztucznie okno na boki, co chroni przed wycinaniem sygnału pod oknem Fouriera
        left = max(1, left - 2)
        right = min(N, right + 2)
        
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

"""
Calculates a localized Fourier Phase Gradient shift within a specific window.
Uses a Hann taper to mitigate spectral leakage at window boundaries.
"""
function windowed_fourier_shift(prof_ref::AbstractVector, prof_target::AbstractVector, left::Int, right::Int)
    N_full = length(prof_ref)
    left = max(1, left)
    right = min(N_full, right)
    
    w_ref = prof_ref[left:right]
    w_tar = prof_target[left:right]
    
    N = length(w_ref)
    if N < 5
        return NaN
    end
    
    # Mnożymy profil w oknie przez funkcję Hanna, aby zachowywał się jak lokalny rozkład ucinający szum (podobnie do fitu Gaussa)
    taper = 0.5 .* (1 .- cos.(2 * pi .* (0:N-1) ./ (N - 1)))
    
    bg_ref = min(w_ref[1], w_ref[end])
    bg_tar = min(w_tar[1], w_tar[end])
    
    w_ref = (w_ref .- bg_ref) .* taper
    w_tar = (w_tar .- bg_tar) .* taper
    
    shift_res = global_fourier_shift(w_ref, w_tar)
    
    if abs(shift_res.shift) > N / 2
        return (shift=NaN, error=NaN)
    end
    
    return shift_res
end

"""
Calculates offsets per component using localized Fourier Phase Gradients 
instead of purely relying on Barycenters. Provides robustness similar to Gaussian fitting.
"""
function component_fourier_offsets_fixed_windows(prof_ref::AbstractVector, prof_target::AbstractVector, windows)
    min_ref, max_ref = minimum(prof_ref), maximum(prof_ref)
    min_tar, max_tar = minimum(prof_target), maximum(prof_target)
    
    amp_ref = max_ref - min_ref
    amp_tar = max_tar - min_tar
    
    if amp_ref ≈ 0.0 || amp_tar ≈ 0.0
        return [(com_ref=NaN, offset=NaN) for _ in windows]
    end
    
    p_ref = (prof_ref .- min_ref) ./ amp_ref
    p_tar = (prof_target .- min_tar) ./ amp_tar
    
    return map(windows) do (left, right)
        x_win = left:right
        
        # Do ustalenia horyzontalnego piku (Długość geograficzna fazy/Longitude) używamy środka masy
        bg_ref = min(p_ref[left], p_ref[right])
        vals_ref = max.(p_ref[x_win] .- bg_ref, 0.0)
        m_ref = sum(vals_ref)
        com_ref = m_ref > 0 ? sum(x_win .* vals_ref) / m_ref : NaN
        
        # Do precyzyjnego liczenia offsetu między częstotliwościami używamy lokalnego Fouriera
        offset_res = windowed_fourier_shift(p_ref, p_tar, left, right)
        
        return (com_ref=com_ref, offset=offset_res.shift, error=offset_res.error)
    end
end

"""
Dynamically identifies signal islands (components) in the current profile 
and calculates their localized Fourier Phase Gradients. 
Automatically adapts to however many sub-pulses appear in the specific phase.
"""
function dynamic_component_fourier_offsets(prof_ref::AbstractVector, prof_target::AbstractVector; threshold=0.03)
    N = length(prof_ref)
    
    min_ref, max_ref = minimum(prof_ref), maximum(prof_ref)
    min_tar, max_tar = minimum(prof_target), maximum(prof_target)
    
    amp_ref = max_ref - min_ref
    amp_tar = max_tar - min_tar
    
    if amp_ref ≈ 0.0 || amp_tar ≈ 0.0
        return NamedTuple{(:com_ref, :offset, :left, :right), Tuple{Float64, Float64, Int, Int}}[]
    end
    
    p_ref = (prof_ref .- min_ref) ./ amp_ref
    p_tar = (prof_target .- min_tar) ./ amp_tar
    
    # Wygładzamy profil by nie łapać szumowych igieł jako osobnych komponentów
    smoothed_ref = copy(p_ref)
    for i in 3:N-2
        smoothed_ref[i] = mean(p_ref[i-2:i+2])
    end
    
    # Szukamy lokalnych maksimów
    peaks = Int[]
    for i in 2:N-1
        if smoothed_ref[i] > smoothed_ref[i-1] && smoothed_ref[i] > smoothed_ref[i+1] && smoothed_ref[i] > threshold
            push!(peaks, i)
        end
    end
    
    # Sortujemy najsilniejsze najpierw
    sort!(peaks, by=p -> smoothed_ref[p], rev=true)
    
    results = NamedTuple{(:com_ref, :offset, :error, :left, :right), Tuple{Float64, Float64, Float64, Int, Int}}[]
    visited = falses(N)
    
    for p in peaks
        if visited[p]
            continue
        end
        
        cut_off = threshold * 0.5 # Schodzimy niżej niż próg detekcji, by złapać ogony
        
        left = p
        while left > 1 && smoothed_ref[left] > cut_off && smoothed_ref[left-1] <= smoothed_ref[left]
            left -= 1
        end
        
        right = p
        while right < N && smoothed_ref[right] > cut_off && smoothed_ref[right+1] <= smoothed_ref[right]
            right += 1
        end
        
        visited[left:right] .= true
        
        left = max(1, left - 3)
        right = min(N, right + 3)
        
        if right - left < 5; continue; end
        
        x_win = left:right
        bg_ref = min(p_ref[left], p_ref[right])
        vals_ref = max.(p_ref[x_win] .- bg_ref, 0.0)
        m_ref = sum(vals_ref)
        com_ref = m_ref > 0 ? sum(x_win .* vals_ref) / m_ref : NaN
        
        offset_res = windowed_fourier_shift(p_ref, p_tar, left, right)
        
        if !isnan(com_ref) && !isnan(offset_res.shift)
            push!(results, (com_ref=com_ref, offset=offset_res.shift, error=offset_res.error, left=left, right=right))
        end
    end
    
    return results
end

end # module