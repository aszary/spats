module LrfsDiagnostics

using FFTW
using Statistics
using LinearAlgebra

export lrfs_phase_track

"""
    _lrfs(X) -> (F::Matrix{ComplexF64}, intensity::Vector{Float64})

Self-contained LRFS computation matching Tools.lrfs logic:
- Real FFT (rfft) along the pulse axis for each longitude bin, which includes the
  Nyquist bin and avoids the half-range truncation of the old `fft(…)[:,1:half]` approach.
- Returns the complex spectrum matrix F (half_freqs × bins) and summed power per row.
"""
function _lrfs(X::AbstractMatrix{<:Real})
    da = collect(transpose(X))   # (bins × pulses)

    # rfft along the pulse dimension: output shape (bins × (pulse_num÷2 + 1))
    raw = rfft(da, 2)            # includes Nyquist; no manual slicing needed
    F   = collect(transpose(raw))  # (half_freqs × bins) — same layout as Tools.lrfs

    intensity = vec(sum(abs2.(F), dims=2))  # sum power (abs²) across bins per freq row
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
    lrfs_phase_track(X, bin_st, bin_end; known_p3=nothing)

Determines if a pulsar is P3-only or drifting by tracking the complex phase of
a specific P3 frequency across longitude bins.

If `known_p3` (in pulses) is provided, the corresponding frequency row is used
directly — this is strongly preferred when the P3 value is known from the
catalogue (e.g. from the input list file) because the auto-detected dominant
peak can land on a harmonic or a noise spike for weak/noisy pulsars.

If `known_p3` is `nothing`, the dominant power peak in the LRFS intensity is
used as a fallback (auto-detection).
"""
function lrfs_phase_track(X::AbstractMatrix{<:Real}, bin_st::Int, bin_end::Int;
                          known_p3::Union{Float64,Nothing}=nothing)
    # Clamp bounds and check valid dimensions
    bin_st  = max(1, bin_st)
    bin_end = min(size(X, 2), bin_end)
    if bin_st >= bin_end || size(X, 1) < 4
        return (p3_pulses = 0.0, phase_slope = 0.0,
                classification = :undetermined, phase_track = Float64[],
                p3_source = :invalid)
    end

    # 1. Preprocess: isolate on-pulse and remove static profile
    X_on   = view(X, :, bin_st:bin_end)
    X_prep = X_on .- mean(X_on, dims=1)

    N, M = size(X_prep)

    # 2. Compute LRFS (self-contained, no external dependency)
    F, intensity = _lrfs(X_prep)

    if length(intensity) < 2
        return (p3_pulses = 0.0, phase_slope = 0.0,
                classification = :undetermined, phase_track = Float64[],
                p3_source = :invalid)
    end

    # 3. Determine which frequency row to analyse
    peak_idx, p3_source = if !isnothing(known_p3) && known_p3 > 1.0
        # Convert known P3 (in pulses) to the nearest LRFS frequency index.
        # The k-th bin of rfft output (1-based) corresponds to frequency (k-1)/N
        # cycles/pulse, so P3 = N/(k-1).  Solving: k = round(N/known_p3) + 1.
        # Clamp to valid range [2, length(intensity)] (index 1 is DC).
        k = clamp(round(Int, N / known_p3) + 1, 2, length(intensity))
        k, :known_p3
    else
        # Fallback: find dominant power peak, skipping DC (index 1)
        argmax(view(intensity, 2:length(intensity))) + 1, :auto_detected
    end

    # P3 actually used (may differ slightly from known_p3 due to frequency quantisation).
    # rfft bin k (1-based) ↔ frequency (k-1)/N → P3 = N/(k-1).
    p3_value = N / (peak_idx - 1)

    # 4. Extract the complex row at the chosen frequency.
    # Conjugate to match the PhaseDrift sign convention (e^{+2πi f₃ n}):
    # positive phase slope → positive drift (later subpulses arrive at later longitudes).
    Z = conj.(vec(F[peak_idx, :]))

    # Phase track for diagnostic output (unwrap from bin 1 outward)
    phase_track = angle.(Z)
    unwrap!(phase_track)

    # 5. Coherent slope estimator — robust against phase-unwrap failures in
    # low-SNR off-pulse bins.  g = Σ conj(Z[j])·Z[j+1]; slope = arg(g) in rad/bin.
    # This is equivalent to the amplitude-weighted mean phase increment and does
    # not rely on unwrapping.
    if length(Z) < 2
        return (p3_pulses = p3_value, phase_slope = 0.0,
                classification = :undetermined, phase_track = phase_track,
                p3_source = p3_source)
    end

    g = sum(conj(Z[j]) * Z[j+1] for j in 1:(length(Z)-1))

    slope = abs(g) > 0 ? angle(g) : 0.0   # radians per longitude bin

    # 6. Classification: flat phase → amplitude modulation; sloped → drift
    classification = abs(slope) > 0.05 ? :pure_drift : :amplitude_modulation

    return (
        p3_pulses      = p3_value,
        phase_slope    = slope,
        classification = classification,
        phase_track    = phase_track,
        p3_source      = p3_source   # :known_p3 or :auto_detected
    )
end

end # module
