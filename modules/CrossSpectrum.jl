module CrossSpectrum

using FFTW
using Statistics
using LinearAlgebra

export cross_spectrum, phase_gradient

"""
    cross_spectrum(X::AbstractMatrix{<:Real}, bin1::Int, bin2::Int)

Computes the complex cross-spectral density S_12(f) = F{X_1} * F{X_2}^* 
between two longitude bins of the single-pulse matrix `X`.
"""
function cross_spectrum(X::AbstractMatrix{<:Real}, bin1::Int, bin2::Int)
    F1 = rfft(view(X, :, bin1))
    F2 = rfft(view(X, :, bin2))
    return F1 .* conj.(F2)
end

"""
    unwrap!(phase::AbstractVector{<:Real})

Unwraps phase jumps to ensure continuity for linear regression.
"""
function unwrap!(phase::AbstractVector{<:Real})
    @inbounds for i in 2:length(phase)
        diff = phase[i] - phase[i-1]
        if diff > π
            phase[i:end] .-= 2π
        elseif diff < -π
            phase[i:end] .+= 2π
        end
    end
    return phase
end

"""
    phase_gradient(X::AbstractMatrix{<:Real})

Computes pairwise cross-spectral phase differences across adjacent longitude bins.
Optimized to compute 1D rfft across all columns in a single batched FFTW operation.
Returns the phase spectrum, linear regression fit dθ/df, and cross-coherence.
"""
function phase_gradient(X::AbstractMatrix{<:Real})
    N, M = size(X)
    if M < 2 || N < 4
        return (phase = Float64[], gradient = 0.0, coherence = 0.0)
    end
    
    freqs = rfftfreq(N)
    
    # Fast batched 1D rfft along time axis (dim 1) for all longitude bins at once
    FX = rfft(X, 1)  # size: (N ÷ 2 + 1, M)
    K = size(FX, 1)
    
    avg_cross_spec = zeros(ComplexF64, K)
    
    # Vectorized accumulation across adjacent bins without allocating intermediate arrays
    @inbounds for i in 2:M
        for k in 1:K
            avg_cross_spec[k] += FX[k, i] * conj(FX[k, i-1])
        end
    end
    
    phase = angle.(avg_cross_spec)
    unwrap!(phase)
    
    mag = abs.(avg_cross_spec)
    
    # Ignore DC component for regression
    sum_w = sum(@view mag[2:end])
    if sum_w == 0.0
        return (phase = phase, gradient = 0.0, coherence = 0.0)
    end
    
    # Weighted linear regression: y = m*x (forcing through origin for phase delay)
    x = @view freqs[2:end]
    y = @view phase[2:end]
    w = @view mag[2:end]
    
    # Calculate gradient (m)
    gradient = sum(w .* x .* y) / sum(w .* x.^2)
    
    # Total auto power using already-computed FX
    total_auto_power = sum(abs2, FX) / M
    coherence = sum_w / (total_auto_power + 1e-12)
    
    return (phase = phase, gradient = gradient, coherence = coherence)
end

end # module
