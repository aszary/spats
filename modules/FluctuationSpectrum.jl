module FluctuationSpectrum

using FFTW
using Statistics

export lrfs, tdfs

"""
    lrfs(X::AbstractMatrix{<:Real})

Computes the Longitude-Resolved Fluctuation Spectrum (LRFS).
Takes an N_pulses x M_bins matrix `X` and performs a 1D real FFT along the pulses dimension (dim 1).
Returns the power spectrum matrix.
"""
function lrfs(X::AbstractMatrix{<:Real})
    F = rfft(X, 1)
    return abs2.(F)
end

"""
    tdfs(X::AbstractMatrix{<:Real})

Computes the 2D Fluctuation Spectrum (2DFS) across both time and spatial axes.
Takes an N_pulses x M_bins matrix `X`. 
Returns a NamedTuple containing the shifted power spectrum, P3 and P2 centroids, 
and the asymmetry ratio of power across spatial frequency 1/P2 = 0.
"""
function tdfs(X::AbstractMatrix{<:Real})
    N, M = size(X)
    if N < 2 || M < 2
        return (
            power = zeros(Float64, N, M),
            p3_centroid = 0.0,
            p2_centroid = 0.0,
            asymmetry_ratio = 0.0,
            symmetric_fraction = 1.0
        )
    end
    
    # 2D FFT
    F = fft(X)
    power = abs2.(F)
    
    # Shift zero frequency to center
    power_shifted = fftshift(power)
    
    # Compute asymmetry across spatial frequency (dim 2)
    center_M = M ÷ 2 + 1
    
    left_power = sum(@view power_shifted[:, 1:center_M-1])
    right_power = sum(@view power_shifted[:, center_M+1:M])
    
    total_power = left_power + right_power
    asymmetry_ratio = total_power > 0 ? (left_power - right_power) / total_power : 0.0
    
    # Centroid estimation (center of mass)
    row_sum = vec(sum(power_shifted, dims=2))
    col_sum = vec(sum(power_shifted, dims=1))
    
    s_row = sum(row_sum)
    s_col = sum(col_sum)
    
    p3_idx_centroid = s_row > 0 ? sum((1:N) .* row_sum) / s_row : (N ÷ 2 + 1)
    p2_idx_centroid = s_col > 0 ? sum((1:M) .* col_sum) / s_col : center_M
    
    # Map index back to relative frequency (cycles per bin/pulse)
    p3_centroid = (p3_idx_centroid - (N ÷ 2 + 1)) / N
    p2_centroid = (p2_idx_centroid - center_M) / M
    
    symmetric_fraction = 1.0 - abs(asymmetry_ratio)

    return (
        power = power_shifted, 
        p3_centroid = p3_centroid, 
        p2_centroid = p2_centroid, 
        asymmetry_ratio = asymmetry_ratio,
        symmetric_fraction = symmetric_fraction
    )
end

end # module
