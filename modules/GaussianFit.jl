module GaussianFit

using LsqFit
using LinearAlgebra: diag, svd
using Printf

export fit_gaussians, best_ngaussians, component_offsets, print_fit_summary

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

"""
    _gauss(x, A, mu, sigma)

Evaluate a single Gaussian with amplitude `A`, centre `mu`, and width `sigma`.
"""
_gauss(x, A, mu, sigma) = A * exp(-0.5 * ((x - mu) / sigma)^2)

"""
    _multi_gauss(x, p, n)

Evaluate a sum of `n` Gaussians plus a constant baseline.

Parameter vector layout:
    p = [baseline, A1, mu1, sigma1, A2, mu2, sigma2, ...]
"""
function _multi_gauss(x::AbstractVector, p::AbstractVector, n::Int)
    result = fill(p[1], length(x))
    for i in 1:n
        A, mu, sigma = p[2 + 3*(i-1)], p[3 + 3*(i-1)], p[4 + 3*(i-1)]
        result .+= _gauss.(x, A, mu, sigma)
    end
    return result
end

"""
    _auto_p0(x, y, n)

Estimate initial parameters for `n` Gaussians by iteratively finding the
largest residual peak and subtracting it.
"""
function _auto_p0(x::AbstractVector, y::AbstractVector, n::Int)
    baseline_est = minimum(y)
    p0 = [baseline_est]
    residual = copy(y) .- baseline_est
    for _ in 1:n
        idx     = argmax(residual)
        A_est   = residual[idx]
        mu_est  = x[idx]
        sig_est = (maximum(x) - minimum(x)) / (4 * n)
        append!(p0, [A_est, mu_est, sig_est])
        residual .-= _gauss.(x, A_est, mu_est, sig_est)
        residual   = max.(residual, 0.0)
    end
    return p0
end

"""
    _default_bounds(x, y, n)

Return default lower and upper parameter bounds for `n` Gaussians:
- baseline  : [0,  max(y)]
- amplitude : [0,  2*max(y)]
- centre    : [min(x), max(x)]
- sigma     : [1,  (max(x)-min(x))/2]
"""
function _default_bounds(x::AbstractVector, y::AbstractVector, n::Int)
    xmin, xmax = minimum(x), maximum(x)
    ymax       = maximum(y)
    sigma_max = (xmax - xmin) / (1.5 * n)
    lower = vcat([0.0],    repeat([0.0,    xmin, 1.0      ], n))
    upper = vcat([ymax],   repeat([2*ymax, xmax, sigma_max], n))
    return lower, upper
end

"""
    _safe_errors(fit)

Extract parameter uncertainties from the Jacobian-based covariance matrix.
Falls back to a pseudo-inverse (SVD) if the standard `estimate_covar` fails
due to a singular or near-singular matrix (e.g. LAPACKException).

Returns a vector of 1-sigma uncertainties, or a vector of `NaN` if both
methods fail.
"""
function _safe_errors(fit)
    # Primary: standard LsqFit covariance
    try
        cov = estimate_covar(fit)
        vars = diag(cov)
        if all(isfinite, vars) && all(>=(0), vars)
            return sqrt.(vars)
        end
    catch
    end

    # Fallback: pseudo-inverse of (J'J) via SVD
    try
        J   = fit.jacobian
        JtJ = J' * J
        F   = svd(JtJ)
        # Zero out singular values below threshold to regularise
        threshold = maximum(F.S) * length(F.S) * eps(eltype(F.S))
        S_inv = [s > threshold ? 1/s : 0.0 for s in F.S]
        cov_pinv = F.V * Diagonal(S_inv) * F.Vt
        vars = diag(cov_pinv)
        if all(isfinite, vars) && all(>=(0), vars)
            @warn "estimate_covar failed; uncertainties computed via SVD pseudo-inverse"
            return sqrt.(vars)
        end
    catch
    end

    # Last resort: return NaN so the fit result is still usable
    @warn "Could not estimate parameter uncertainties; returning NaN"
    return fill(NaN, length(coef(fit)))
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    fit_gaussians(x, y, n; p0, lower, upper) -> NamedTuple

Fit a model of `n` Gaussians plus a constant baseline to data `(x, y)`.

# Arguments
- `x`     : phase bins or longitude values
- `y`     : intensity values
- `n`     : number of Gaussian components
- `p0`    : initial parameter vector `[baseline, A1, mu1, s1, ...]`;
            estimated automatically when `nothing`
- `lower` : lower bounds on parameters (estimated automatically when `nothing`)
- `upper` : upper bounds on parameters (estimated automatically when `nothing`)

# Returns a `NamedTuple` with fields
| Field        | Description                                      |
|:-------------|:-------------------------------------------------|
| `params`     | best-fit parameter vector                        |
| `errors`     | 1-sigma uncertainties (NaN if estimation failed) |
| `baseline`   | fitted baseline value                            |
| `components` | Vector of NamedTuples `(A, mu, sigma, A_err, mu_err, sigma_err)` |
| `yfit`       | model evaluated at `x`                           |
| `residuals`  | `y .- yfit`                                      |
| `rms`        | root-mean-square of residuals                    |
| `aic`        | Akaike Information Criterion                     |
| `bic`        | Bayesian Information Criterion                   |
| `converged`  | `true` if the fit succeeded                      |
"""
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
        errors = _safe_errors(fit)   # robust — never throws

        yfit      = model(x, params)
        residuals = y .- yfit
        rms       = sqrt(sum(residuals.^2) / length(residuals))
        k         = 1 + 3 * n
        N         = length(x)
        aic       = N * log(rms^2) + 2 * k
        bic       = N * log(rms^2) + k * log(N)

        components = [(
            A         = params[2 + 3*(i-1)],
            mu        = params[3 + 3*(i-1)],
            sigma     = params[4 + 3*(i-1)],
            A_err     = errors[2 + 3*(i-1)],
            mu_err    = errors[3 + 3*(i-1)],
            sigma_err = errors[4 + 3*(i-1)],
        ) for i in 1:n]

        return (
            params     = params,
            errors     = errors,
            baseline   = params[1],
            components = components,
            yfit       = yfit,
            residuals  = residuals,
            rms        = rms,
            aic        = aic,
            bic        = bic,
            converged  = true,
        )

    catch e
        @warn "Gaussian fit failed: $e"
        return (
            params=nothing, errors=nothing, baseline=nothing,
            components=nothing, yfit=nothing, residuals=nothing,
            rms=Inf, aic=Inf, bic=Inf, converged=false,
        )
    end
end

"""
    best_ngaussians(x, y; max_n=4) -> (best_fit, best_n, aics)

Try fitting 1 to `max_n` Gaussians and return the result with the lowest AIC.

# Returns
- `best_fit` : result of `fit_gaussians` for the winning number of components
- `best_n`   : integer, optimal number of Gaussians
- `aics`     : vector of AIC values for N = 1..max_n
"""
function best_ngaussians(x::AbstractVector, y::AbstractVector; max_n::Int = 4)
    results = [fit_gaussians(x, y, n) for n in 1:max_n]
    aics    = [r.aic for r in results]
    best_n  = argmin(aics)
    return results[best_n], best_n, aics
end

"""
    component_offsets(fit_high, fit_low; nbin=1024) -> Vector

Compute the phase offset between matched Gaussian components of two profiles
observed at different radio frequencies.

# Arguments
- `fit_high` : result of `fit_gaussians` for the high-frequency profile
- `fit_low`  : result of `fit_gaussians` for the low-frequency profile
- `nbin`     : number of phase bins per period (default 1024)

# Returns a Vector of NamedTuples with fields:
`component, mu_high, mu_low, offset_bins, offset_err, offset_deg`
"""
function component_offsets(fit_high, fit_low; nbin::Int = 1024)

    nh = length(fit_high.components)
    nl = length(fit_low.components)
    @assert nh == nl "Number of components must match ($nh vs $nl)"

    # Sort both by mu so components at similar positions are matched together
    ch = sort(fit_high.components, by = c -> c.mu)
    cl = sort(fit_low.components,  by = c -> c.mu)

    return [(
        component      = i,
        mu_high        = ch[i].mu,
        mu_low         = cl[i].mu,
        offset_bins    = ch[i].mu - cl[i].mu,
        offset_err     = sqrt(ch[i].mu_err^2 + cl[i].mu_err^2),
        offset_deg     = (ch[i].mu - cl[i].mu) / nbin * 360.0,
        offset_deg_err = sqrt(ch[i].mu_err^2 + cl[i].mu_err^2) / nbin * 360.0,
        longitude      = (ch[i].mu + cl[i].mu) / 2.0 / nbin * 360.0,
        longitude_bin  = (ch[i].mu + cl[i].mu) / 2.0,
    ) for i in 1:nh]
end

"""
    print_fit_summary(fit, n; label="", nbin=1024)

Print a human-readable summary of a `fit_gaussians` result.
"""
function print_fit_summary(fit, n::Int;
                            label::String = "",
                            nbin::Int = 1024)
    isempty(label) || println("=== $label ===")
    @printf("  N components : %d\n",   n)
    @printf("  RMS          : %.5f\n", fit.rms)
    @printf("  AIC          : %.1f\n", fit.aic)
    @printf("  BIC          : %.1f\n", fit.bic)
    @printf("  baseline     : %.4f\n", fit.baseline)
    println()
    for (i, c) in enumerate(fit.components)
        deg = c.mu / nbin * 360.0
        @printf("  G%-2d  A     = %.4f ± %.4f\n", i, c.A, c.A_err)
        @printf("       mu    = %.2f ± %.2f bins  (%.2f°)\n", c.mu, c.mu_err, deg)
        @printf("       sigma = %.2f ± %.2f bins\n", c.sigma, c.sigma_err)
    end
end

end # module GaussianFit
