module CrossCorrelation

using Statistics
using DSP

export ccf_map, peak_lags

"""
    ccf_map(X::AbstractMatrix{<:Real}, max_lag::Int)

Computes pair-wise temporal cross-correlation between all longitudes up to `max_lag`.
Returns a 3D matrix (Lags × Bin1 × Bin2) and the lag vector.
"""
function ccf_map(X::AbstractMatrix{<:Real}, max_lag::Int)
    N, M = size(X)
    lags = -max_lag:max_lag
    L = length(lags)
    
    ccf_matrix = zeros(Float64, L, M, M)
    
    for i in 1:M
        ts_i = view(X, :, i) .- mean(view(X, :, i))
        for j in 1:M
            ts_j = view(X, :, j) .- mean(view(X, :, j))
            
            # Full cross correlation
            c = xcorr(ts_i, ts_j)
            
            # xcorr length is 2N-1, center is at index N (which corresponds to lag 0)
            center = N
            for (k, lag) in enumerate(lags)
                idx = center + lag
                if 1 <= idx <= length(c)
                    ccf_matrix[k, i, j] = c[idx]
                end
            end
        end
    end
    
    return ccf_matrix, lags
end

"""
    peak_lags(X::AbstractMatrix{<:Real}, max_lag::Int)

Extracts the temporal lag tau_peak at which correlation peaks for each longitude pair.
Returns an M x M matrix of peak lags, and estimated mean drift velocity.
"""
function peak_lags(X::AbstractMatrix{<:Real}, max_lag::Int)
    N, M = size(X)
    ccf_matrix, lags = ccf_map(X, max_lag)
    
    peak_lag_mat = zeros(Int, M, M)
    
    for i in 1:M
        for j in 1:M
            # Find the index of the maximum correlation value
            max_val, max_idx = findmax(view(ccf_matrix, :, i, j))
            peak_lag_mat[i, j] = lags[max_idx]
        end
    end
    
    # Estimate drift velocity = Delta phi / tau_peak
    # Average across adjacent bins (Delta phi = 1)
    velocity_estimates = Float64[]
    for i in 2:M
        tau = peak_lag_mat[i, i-1]
        if tau != 0
            push!(velocity_estimates, 1.0 / tau)
        end
    end
    
    mean_velocity = isempty(velocity_estimates) ? 0.0 : mean(velocity_estimates)
    
    return (peak_lag_mat = peak_lag_mat, mean_velocity = mean_velocity)
end

end # module
