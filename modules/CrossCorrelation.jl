module CrossCorrelation

using Statistics
using FFTW

export ccf_map, peak_lags

"""
    ccf_adjacent(X::AbstractMatrix{<:Real}, max_lag::Int)

Computes temporal cross-correlation only between **adjacent** longitude bins
(i.e., bin i and bin i+1). O(M) complexity with FFT reuse across iterations.
"""
function ccf_adjacent(X::AbstractMatrix{<:Real}, max_lag::Int)
    N, M = size(X)
    lags = collect(-max_lag:max_lag)
    L    = length(lags)

    if M < 2 || N < 4
        return zeros(Float64, max(0, M - 1), L), lags
    end

    # Pre-plan FFT
    nfft = nextpow(2, 2N)
    plan_fwd = plan_rfft(zeros(Float64, nfft))
    plan_inv = plan_irfft(zeros(ComplexF64, nfft ÷ 2 + 1), nfft)

    ccf_mat = zeros(Float64, M - 1, L)

    # Reusable buffers
    pad_buf = zeros(Float64, nfft)
    
    # Compute first column's FFT
    v1 = @view X[:, 1]
    m1 = mean(v1)
    @inbounds for k in 1:N
        pad_buf[k] = v1[k] - m1
    end
    @inbounds for k in N+1:nfft
        pad_buf[k] = 0.0
    end
    FA = plan_fwd * pad_buf

    for i in 1:(M - 1)
        # Compute column i+1's FFT into FB
        v2 = @view X[:, i+1]
        m2 = mean(v2)
        @inbounds for k in 1:N
            pad_buf[k] = v2[k] - m2
        end
        @inbounds for k in N+1:nfft
            pad_buf[k] = 0.0
        end
        FB = plan_fwd * pad_buf

        # Circular cross-correlation via product in frequency domain
        C = FA .* conj.(FB)
        r = plan_inv * C
        r ./= N

        # Extract lags
        for (li, lag) in enumerate(lags)
            idx = lag >= 0 ? lag + 1 : nfft + lag + 1
            ccf_mat[i, li] = r[idx]
        end

        # For the next iteration, column i+1 becomes column i: reuse FB as FA!
        FA = FB
    end

    return ccf_mat, lags
end

"""
    ccf_map(X::AbstractMatrix{<:Real}, max_lag::Int)

For backward compatibility: calls ccf_adjacent and returns a 3D matrix
(Lags × M × M) with super-diagonal and sub-diagonal filled.
"""
function ccf_map(X::AbstractMatrix{<:Real}, max_lag::Int)
    N, M = size(X)
    lags = collect(-max_lag:max_lag)
    L    = length(lags)
    ccf_matrix = zeros(Float64, L, M, M)
    adj, _ = ccf_adjacent(X, max_lag)
    for i in 1:(M-1)
        ccf_matrix[:, i, i+1] .= view(adj, i, :)
        ccf_matrix[:, i+1, i] .= reverse(view(adj, i, :))
    end
    return ccf_matrix, lags
end

"""
    peak_lags(X::AbstractMatrix{<:Real}, max_lag::Int)

Extracts the peak lag between adjacent longitude bins and estimates drift velocity.
Runs in O(M) time.
"""
function peak_lags(X::AbstractMatrix{<:Real}, max_lag::Int)
    N, M = size(X)
    if M < 2 || N < 4
        return (peak_lag_mat = zeros(Int, M, M), mean_velocity = 0.0)
    end

    adj_ccf, lags = ccf_adjacent(X, max_lag)

    adj_peak_lags = Vector{Int}(undef, M - 1)
    for i in 1:(M - 1)
        adj_peak_lags[i] = lags[argmax(view(adj_ccf, i, :))]
    end

    peak_lag_mat = zeros(Int, M, M)
    for i in 1:(M - 1)
        peak_lag_mat[i, i+1] =  adj_peak_lags[i]
        peak_lag_mat[i+1, i] = -adj_peak_lags[i]
    end

    velocity_estimates = [1.0 / τ for τ in adj_peak_lags if τ != 0]
    mean_velocity = isempty(velocity_estimates) ? 0.0 : mean(velocity_estimates)

    return (peak_lag_mat = peak_lag_mat, mean_velocity = mean_velocity)
end

end # module
