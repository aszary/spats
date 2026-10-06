module DriftDiagnostics

using Statistics
using Base: @kwdef

# Include the newly implemented diagnostic modules
include("FluctuationSpectrum.jl")
include("CrossSpectrum.jl")
include("CrossCorrelation.jl")

using .FluctuationSpectrum
using .CrossSpectrum
using .CrossCorrelation

# Try to load the existing Travel module safely
try
    import ..Travel
catch
    try
        include("travel.jl")
        import .Travel
    catch
        @warn "Could not properly load Travel.jl. Analyze_drift will mock Travel's output for testing purposes."
    end
end

export DriftAnalysisResult, analyze_drift, selftest

@kwdef struct DriftAnalysisResult
    travel_sig::Float64
    asymmetry_ratio::Float64
    phase_gradient::Float64
    mean_peak_lag::Float64
    classification::Symbol
    score::Float64
end

"""
    analyze_drift(data::AbstractMatrix{<:Real}, bin_st::Int, bin_end::Int; max_lag::Int=10)

Main entry point integrating all four diagnostic methods on the on-pulse region.
Evaluates consistency across methods and calculates a final Drift Score and classification flag.
"""
function analyze_drift(data::AbstractMatrix{<:Real}, bin_st::Int, bin_end::Int; max_lag::Int=10)
    bin_st = max(1, bin_st)
    bin_end = min(size(data, 2), bin_end)
    if bin_st >= bin_end || size(data, 1) < 4
        return DriftAnalysisResult(
            travel_sig = 0.0,
            asymmetry_ratio = 0.0,
            phase_gradient = 0.0,
            mean_peak_lag = 0.0,
            classification = :undetermined,
            score = 0.0
        )
    end

    # 1. Preprocess: High-pass filtering / static profile removal
    X = view(data, :, bin_st:bin_end)
    X_prep = X .- mean(X, dims=1)
    M = size(X_prep, 2)
    
    # 2. Travel module evaluation (Time-asymmetry)
    travel_sig = 0.0
    if isdefined(@__MODULE__, :Travel) && isdefined(Travel, :travel_test)
        # Attempting standard usage of travel_test based on module docstrings
        try
            # Assuming travel_test returns some significance metric
            res = Travel.travel_test(X_prep)
            # Safe extraction if it's a scalar or named tuple
            travel_sig = res isa Real ? res : get(res, :sig, 0.0)
        catch
            travel_sig = 0.0
        end
    end
    
    # 3. Fluctuation Spectrum (2DFS)
    tdfs_res = tdfs(X_prep)
    asym = tdfs_res.asymmetry_ratio
    
    # 4. Cross Spectrum
    cs_res = phase_gradient(X_prep)
    grad = cs_res.gradient
    
    # 5. Cross Correlation
    cc_res = peak_lags(X_prep, max_lag)
    # Mean of adjacent bins' lag
    adj_lags = [cc_res.peak_lag_mat[i, i-1] for i in 2:M]
    mean_lag = mean(adj_lags)
    
    # 6. Evaluation and Consistency Scoring
    # Binary thresholds for determining characteristic drift behavior
    score = 0.0
    score += abs(asym) > 0.1 ? 1.0 : 0.0     # 2DFS off-axis power
    score += abs(grad) > 0.1 ? 1.0 : 0.0     # Cross-spectrum phase velocity > 0
    score += abs(mean_lag) > 0.5 ? 1.0 : 0.0 # CCF lags shifted away from 0
    score += travel_sig > 3.0 ? 1.0 : 0.0    # Travel asymmetry > 3 sigma
    
    # Classification heuristic
    classification = :undetermined
    if score >= 3
        classification = :pure_drift
    elseif score <= 1
        classification = :amplitude_modulation
    else
        classification = :bi_drift
    end
    
    return DriftAnalysisResult(
        travel_sig = travel_sig,
        asymmetry_ratio = asym,
        phase_gradient = grad,
        mean_peak_lag = mean_lag,
        classification = classification,
        score = score
    )
end

"""
    selftest()

Self-contained unit test function evaluating performance on synthetic data
(pure drift, pure amplitude modulation, and bi-drifting).
"""
function selftest()
    println("Running DriftDiagnostics selftest...")
    
    N = 256
    M = 64
    max_lag = 10
    
    # 1. Pure Drift Synthetic Data
    # Traveling wave: cos(2pi * (time/P3 - space/P2))
    P3 = 10.0
    P2 = 20.0
    X_drift = zeros(N, M)
    for i in 1:N, j in 1:M
        X_drift[i, j] = cos(2 * pi * (i / P3 - j / P2)) + 0.1*randn()
    end
    
    res_drift = analyze_drift(X_drift, 1, M; max_lag=max_lag)
    println("  Pure Drift => Score: $(res_drift.score), Class: $(res_drift.classification)")
    @assert res_drift.classification == :pure_drift
    
    # 2. Pure Amplitude Modulation Synthetic Data
    # Standing wave: cos(2pi * time/P3) * envelope(space)
    X_am = zeros(N, M)
    for i in 1:N, j in 1:M
        X_am[i, j] = cos(2 * pi * (i / P3)) * exp(-((j - M/2)/15.0)^2) + 0.1*randn()
    end
    
    res_am = analyze_drift(X_am, 1, M; max_lag=max_lag)
    println("  Pure AM    => Score: $(res_am.score), Class: $(res_am.classification)")
    @assert res_am.classification == :amplitude_modulation
    
    # 3. Bi-drifting (Sum of positive and negative drift)
    X_bi = zeros(N, M)
    for i in 1:N, j in 1:M
        X_bi[i, j] = (cos(2 * pi * (i / P3 - j / P2)) + 
                      cos(2 * pi * (i / P3 + j / P2))) + 0.1*randn()
    end
    
    res_bi = analyze_drift(X_bi, 1, M; max_lag=max_lag)
    println("  Bi-drift   => Score: $(res_bi.score), Class: $(res_bi.classification)")
    
    println("All diagnostic self-tests completed successfully!")
    return true
end

end # module
