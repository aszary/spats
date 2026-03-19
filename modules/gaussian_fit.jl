module GaussianFit

using LsqFit
using LinearAlgebra: diag, svd, Diagonal
using Printf
using Statistics: mean

export fit_gaussians, best_ngaussians, component_offsets, print_fit_summary, evaluate_model, evaluate_components

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

"""Evaluate a single Gaussian"""
_gauss(x, A, mu, sigma) = A * exp(-0.5 * ((x - mu) / sigma)^2)

"""Evaluate a sum of `n` Gaussians plus a constant baseline"""
function _multi_gauss(x::AbstractVector, p::AbstractVector, n::Int)
    result = fill(p[1], length(x))
    for i in 1:n
        A, mu, sigma = p[2 + 3*(i-1)], p[3 + 3*(i-1)], p[4 + 3*(i-1)]
        result .+= _gauss.(x, A, mu, sigma)
    end
    return result
end

"""Moving average to smooth residuals before finding peaks for initial guesses"""
function _smooth(y::AbstractVector, window::Int=5)
    half = window ÷ 2
    len = length(y)
    return [mean(y[max(1, i-half):min(len, i+half)]) for i in 1:len]
end

"""Estimate initial parameters for `n` Gaussians by finding the largest residual peak"""
function _auto_p0(x::AbstractVector, y::AbstractVector, n::Int)
    baseline_est = minimum(y)
    p0 = [baseline_est]
    residual = max.(copy(y) .- baseline_est, 0.0)
    x_range  = maximum(x) - minimum(x)
    
    for _ in 1:n
        smoothed_res = _smooth(residual, 5) # Smooth to avoid noise spikes
        idx    = argmax(smoothed_res)
        A_est  = residual[idx]
        mu_est = x[idx]

        # Estimate sigma from FWHM 
        half  = A_est / 2.0
        left  = findlast(i -> smoothed_res[i] < half, 1:idx)
        right_rel = findfirst(i -> smoothed_res[i] < half, idx:length(smoothed_res))
        right = isnothing(right_rel) ? nothing : idx - 1 + right_rel
        
        if !isnothing(left) && !isnothing(right) && right > left
            fwhm    = x[right] - x[left]
            sig_est = max(fwhm / 2.355, 1.0)
        else
            sig_est = x_range / (4 * n)
        end

        append!(p0, [A_est, mu_est, sig_est])
        residual .-= _gauss.(x, A_est, mu_est, sig_est)
        residual   = max.(residual, 0.0)
    end
    return p0
end

function _default_bounds(x::AbstractVector, y::AbstractVector, n::Int)
    xmin, xmax = minimum(x), maximum(x)
    ymax       = maximum(y)
    sigma_max  = (xmax - xmin) / (1.5 * n)
    lower = vcat([0.0],    repeat([0.0,    xmin, 0.5      ], n))
    upper = vcat([ymax],   repeat([2*ymax, xmax, sigma_max], n))
    return lower, upper
end

function _safe_errors(fit)
    try
        cov = estimate_covar(fit)
        vars = diag(cov)
        if all(isfinite, vars) && all(>=(0), vars)
            return sqrt.(vars)
        end
    catch
    end

    try
        J   = fit.jacobian
        JtJ = J' * J
        F   = svd(JtJ)
        threshold = maximum(F.S) * length(F.S) * eps(eltype(F.S))
        S_inv = [s > threshold ? 1/s : 0.0 for s in F.S]
        cov_pinv = F.V * Diagonal(S_inv) * F.Vt
        vars = diag(cov_pinv)
        if all(isfinite, vars) && all(>=(0), vars)
            return sqrt.(vars)
        end
    catch
    end
    return fill(NaN, length(coef(fit)))
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

function fit_gaussians(x::AbstractVector, y::AbstractVector, n::Int;
                       p0::Union{AbstractVector, Nothing}    = nothing,
                       lower::Union{AbstractVector, Nothing} = nothing,
                       upper::Union{AbstractVector, Nothing} = nothing)

    p0_val  = isnothing(p0)    ? _auto_p0(x, y, n)        : p0
    lo, hi  = isnothing(lower) ? _default_bounds(x, y, n) : (lower, upper)

    model(xdata, p) = _multi_gauss(xdata, p, n)

    try
        fit    = curve_fit(model, x, y, p0_val; lower=lo, upper=hi)
        params = coef(fit)
        errors = _safe_errors(fit)

        yfit      = model(x, params)
        residuals = y .- yfit
        rms       = sqrt(sum(residuals.^2) / length(residuals))
        k         = 1 + 3 * n
        N         = length(x)
        aic       = N * log(rms^2) + 2 * k
        bic       = N * log(rms^2) + k * log(N)

        # Sort components by mu (position) left to right
        comps = [(A=params[2+3*(i-1)], mu=params[3+3*(i-1)], sigma=params[4+3*(i-1)],
                  A_err=errors[2+3*(i-1)], mu_err=errors[3+3*(i-1)], sigma_err=errors[4+3*(i-1)]) for i in 1:n]
        sort!(comps, by = c -> c.mu)
        
        # Re-sync params array if order changed (helpful when p0 is passed next time)
        sorted_params = [params[1]]
        for c in comps; append!(sorted_params, [c.A, c.mu, c.sigma]); end

        return (params = sorted_params, errors = errors, baseline = params[1], components = comps,
                yfit = yfit, residuals = residuals, rms = rms, aic = aic, bic = bic, converged = true)
    catch e
        return (params=nothing, errors=nothing, baseline=nothing, components=nothing,
                yfit=nothing, residuals=nothing, rms=Inf, aic=Inf, bic=Inf, converged=false)
    end
end

function best_ngaussians(x::AbstractVector, y::AbstractVector; max_n::Int = 3)
    results = [fit_gaussians(x, y, n) for n in 1:max_n]
    valid_results = filter(r -> r.converged, results)
    isempty(valid_results) && return results[1], 1, [Inf]
    
    aics    = [r.aic for r in valid_results]
    best_idx = argmin(aics)
    return valid_results[best_idx], best_idx, aics
end

function component_offsets(fit_high, fit_low; nbin::Int = 1024)
    nh = length(fit_high.components)
    nl = length(fit_low.components)
    @assert nh == nl "Number of components must match ($nh vs $nl)"

    return [(
        component      = i,
        mu_high        = fit_high.components[i].mu,
        mu_low         = fit_low.components[i].mu,
        offset_bins    = fit_high.components[i].mu - fit_low.components[i].mu,
        offset_err     = sqrt(fit_high.components[i].mu_err^2 + fit_low.components[i].mu_err^2)
    ) for i in 1:nh]
end

function evaluate_model(x, fit)
    return fit.converged ? fit.yfit : fill(NaN, length(x))
end

function evaluate_components(x, fit)
    return fit.converged ? [_gauss.(x, c.A, c.mu, c.sigma) for c in fit.components] : []
end

function print_fit_summary(fit, n::Int; label::String = "")
    isempty(label) || println("=== $label ===")
    if !fit.converged
        println("Fit failed to converge.")
        return
    end
    @printf("  RMS: %.5f | AIC: %.1f | Baseline: %.4f\n", fit.rms, fit.aic, fit.baseline)
    for (i, c) in enumerate(fit.components)
        @printf("  G%-2d A = %.4f | mu = %.2f ± %.2f bins | sigma = %.2f\n", i, c.A, c.mu, c.mu_err, c.sigma)
    end
end

end # module