module GausPeakSelect

# ---------------------------------------------------------------------------
#  IRAF fxcor-style interactive Gaussian peak fitting
#
#  Workflow (mirrors fxcor 'm' key):
#    1. Plot the profile
#    2. User clicks on peaks to mark them (like pressing 'm' in fxcor)
#    3. A local Gaussian is fit independently to a window around each mark
#    4. Results are composed into the same NamedTuple format as GaussianFit
#
#  Uses PyPlot.ginput() for interactive point selection.
# ---------------------------------------------------------------------------

using LsqFit
using LinearAlgebra: diag, svd, Diagonal
using Printf
using PyPlot

export select_and_fit, detect_peaks_auto, fit_local_gaussian,
       compose_fit_result, component_offsets, print_fit_summary

# ---------------------------------------------------------------------------
# Internal helpers (shared with GaussianFit)
# ---------------------------------------------------------------------------

"""
    _gauss(x, A, mu, sigma)

Evaluate a single Gaussian with amplitude `A`, centre `mu`, and width `sigma`.
"""
_gauss(x, A, mu, sigma) = A * exp(-0.5 * ((x - mu) / sigma)^2)

"""
    _multi_gauss_eval(x, baseline, components)

Evaluate the full model: baseline + sum of Gaussian components.
Each component is a NamedTuple with fields `A`, `mu`, `sigma`.
"""
function _multi_gauss_eval(x::AbstractVector, baseline::Real, components)
    result = fill(Float64(baseline), length(x))
    for c in components
        result .+= _gauss.(x, c.A, c.mu, c.sigma)
    end
    return result
end

"""
    _safe_errors(fit)

Extract parameter uncertainties from the Jacobian-based covariance matrix.
Falls back to a pseudo-inverse (SVD) if the standard `estimate_covar` fails.
Returns a vector of 1-sigma uncertainties, or a vector of `NaN` if both fail.
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

    @warn "Could not estimate parameter uncertainties; returning NaN"
    return fill(NaN, length(coef(fit)))
end

# ---------------------------------------------------------------------------
# Phase 1 — Baseline estimation
# ---------------------------------------------------------------------------

"""
    estimate_baseline(y; quantile_frac=0.25)

Estimate baseline as the median of the lowest `quantile_frac` fraction of `y`.
More robust than `minimum(y)` against noise spikes.
"""
function estimate_baseline(y::AbstractVector; quantile_frac::Float64 = 0.25)
    sorted = sort(y)
    n_low  = max(1, round(Int, length(sorted) * quantile_frac))
    return Float64(median(sorted[1:n_low]))
end

# Need Statistics for median — inline a simple version to avoid dep
function median(v::AbstractVector)
    sv = sort(v)
    n  = length(sv)
    if isodd(n)
        return sv[div(n, 2) + 1]
    else
        return (sv[div(n, 2)] + sv[div(n, 2) + 1]) / 2.0
    end
end

# ---------------------------------------------------------------------------
# Phase 1b — Automatic peak detection (optional alternative to interactive)
# ---------------------------------------------------------------------------

"""
    detect_peaks_auto(x, y; threshold=3.0, max_peaks=4)

Automatically detect peaks above `threshold × σ_noise`.
Returns vector of `(idx, x_pos, amplitude)` sorted by descending amplitude.
"""
function detect_peaks_auto(x::AbstractVector, y::AbstractVector;
                           threshold::Float64 = 3.0,
                           max_peaks::Int     = 4)
    baseline = estimate_baseline(y)
    residual = y .- baseline

    # Estimate noise from the lower half of residual (off-pulse)
    sorted_res = sort(abs.(residual))
    n_noise    = max(3, round(Int, length(sorted_res) * 0.5))
    sigma_noise = sqrt(sum(sorted_res[1:n_noise] .^ 2) / n_noise)
    sigma_noise = max(sigma_noise, 1e-10)  # guard against zero

    cutoff = threshold * sigma_noise

    # Find local maxima above cutoff
    peaks = Tuple{Int, Float64, Float64}[]   # (idx, x_pos, amplitude)
    for i in 2:(length(y) - 1)
        if residual[i] > residual[i-1] && residual[i] > residual[i+1] && residual[i] > cutoff
            push!(peaks, (i, Float64(x[i]), residual[i]))
        end
    end

    # Sort by amplitude (descending) and take top max_peaks
    sort!(peaks, by = p -> -p[3])
    return peaks[1:min(length(peaks), max_peaks)]
end

# ---------------------------------------------------------------------------
# Phase 2 — Interactive peak selection (the 'm' key from fxcor)
# ---------------------------------------------------------------------------

"""
    _select_peaks_interactive(x, y; figax=nothing)

Display the profile and let the user click on peaks to mark them,
exactly like pressing 'm' in IRAF fxcor.

**Controls:**
- Left-click  → mark a peak at cursor x-position
- Right-click or middle-click → finish selection (or press Enter)

Returns vector of `(idx, x_pos, amplitude)` for the clicked positions,
snapped to the nearest local maximum.
"""
function _select_peaks_interactive(x::AbstractVector, y::AbstractVector;
                                   title_str::String = "Click on peaks to mark (right-click / Enter to finish)")
    fig, ax = subplots(figsize=(10, 4))
    ax.plot(x, y, color="steelblue", lw=1.2)
    ax.set_xlabel("Bin")
    ax.set_ylabel("Intensity")
    ax.set_title(title_str)
    ax.minorticks_on()
    tight_layout()
    draw()
    show(block=false)

    println()
    println("  ┌─────────────────────────────────────────────────────┐")
    println("  │  PEAK SELECTION  (fxcor 'm' mode)                  │")
    println("  │                                                     │")
    println("  │  Left-click  → mark a peak                         │")
    println("  │  Right-click / Enter → done selecting               │")
    println("  └─────────────────────────────────────────────────────┘")

    # Collect clicks — ginput with n=-1 means "until right-click / Enter"
    coords = ginput(-1, timeout=-1)

    close(fig)

    if isempty(coords)
        @warn "No peaks selected."
        return Tuple{Int, Float64, Float64}[]
    end

    # Snap each click to nearest local maximum within ±snap_radius bins
    snap_radius = max(5, round(Int, length(x) * 0.02))  # 2% of total bins or 5
    peaks = Tuple{Int, Float64, Float64}[]

    for (cx, _) in coords
        # Find closest index to the clicked x
        dists = abs.(collect(Float64.(x)) .- cx)
        center_idx = argmin(dists)

        # Search for local max within snap window
        lo = max(1, center_idx - snap_radius)
        hi = min(length(y), center_idx + snap_radius)
        local_max_idx = lo - 1 + argmax(y[lo:hi])

        push!(peaks, (local_max_idx, Float64(x[local_max_idx]), Float64(y[local_max_idx])))
    end

    # Remove duplicates (same index selected twice)
    unique!(p -> p[1], peaks)

    # Sort by x position (left to right)
    sort!(peaks, by = p -> p[2])

    println("  Marked $(length(peaks)) peak(s):")
    for (i, (idx, xp, amp)) in enumerate(peaks)
        @printf("    Peak %d: bin=%.1f  amplitude=%.4f\n", i, xp, amp)
    end

    return peaks
end

# ---------------------------------------------------------------------------
# Phase 3 — Local Gaussian fit to a window around each peak
# ---------------------------------------------------------------------------

"""
    fit_local_gaussian(x, y, peak_idx; window_fraction=0.15, min_halfwidth=5)

Fit a single Gaussian + local baseline to a window around `peak_idx`.

The window extends from the peak outward until the signal drops below
`window_fraction × peak_amplitude`, or at least `min_halfwidth` bins on each side.

Returns a NamedTuple `(A, mu, sigma, A_err, mu_err, sigma_err, baseline_local, converged)`.
"""
function fit_local_gaussian(x::AbstractVector, y::AbstractVector, peak_idx::Int;
                            window_fraction::Float64 = 0.15,
                            min_halfwidth::Int       = 5)

    baseline_est = estimate_baseline(y)
    peak_val     = y[peak_idx] - baseline_est
    threshold    = window_fraction * peak_val + baseline_est

    # Expand window left
    left = peak_idx
    while left > 1 && y[left - 1] > threshold
        left -= 1
    end
    left = min(left, max(1, peak_idx - min_halfwidth))

    # Expand window right
    right = peak_idx
    while right < length(y) && y[right + 1] > threshold
        right += 1
    end
    right = max(right, min(length(y), peak_idx + min_halfwidth))

    # Extract window data
    x_win = Float64.(x[left:right])
    y_win = Float64.(y[left:right])

    # Initial guesses: [baseline_local, A, mu, sigma]
    A_est     = peak_val
    mu_est    = Float64(x[peak_idx])
    # Estimate sigma from FWHM
    half      = A_est / 2.0 + baseline_est
    left_half = findlast(i -> y_win[i] < half, 1:argmax(y_win))
    right_rel = findfirst(i -> y_win[i] < half, argmax(y_win):length(y_win))
    right_half = isnothing(right_rel) ? nothing : argmax(y_win) - 1 + right_rel
    if !isnothing(left_half) && !isnothing(right_half)
        fwhm     = x_win[right_half] - x_win[left_half]
        sig_est  = max(fwhm / 2.355, 1.0)
    else
        sig_est  = (x_win[end] - x_win[1]) / 4.0
    end
    sig_est = max(sig_est, 0.5)

    p0 = [baseline_est, A_est, mu_est, sig_est]

    # Bounds
    xmin_w, xmax_w = x_win[1], x_win[end]
    ymax_w         = maximum(y_win)
    lower = [0.0, 0.0,       xmin_w, 0.5]
    upper = [ymax_w, 2*ymax_w, xmax_w, (xmax_w - xmin_w)]

    # Ensure p0 is within bounds
    p0 = clamp.(p0, lower, upper)

    model(xdata, p) = p[1] .+ _gauss.(xdata, p[2], p[3], p[4])

    try
        fit    = curve_fit(model, x_win, y_win, p0; lower=lower, upper=upper)
        params = coef(fit)
        errors = _safe_errors(fit)

        return (
            A              = params[2],
            mu             = params[3],
            sigma          = params[4],
            A_err          = errors[2],
            mu_err         = errors[3],
            sigma_err      = errors[4],
            baseline_local = params[1],
            window         = (left, right),
            converged      = true,
        )
    catch e
        @warn "Local Gaussian fit failed at peak_idx=$peak_idx: $e"
        return (
            A=Float64(peak_val), mu=Float64(x[peak_idx]), sigma=sig_est,
            A_err=NaN, mu_err=NaN, sigma_err=NaN,
            baseline_local=baseline_est, window=(left, right),
            converged=false,
        )
    end
end

# ---------------------------------------------------------------------------
# Phase 4 — Compose results into the standard fit format
# ---------------------------------------------------------------------------

"""
    compose_fit_result(x, y, local_fits; refine=true)

Take the independent local Gaussian fits and compose them into the same
NamedTuple format returned by `GaussianFit.fit_gaussians`.

If `refine=true` (default), uses the local fits as initial guesses for a
global multi-Gaussian refinement (like the current global approach, but with
much better starting parameters). If `refine=false`, the model is just the
sum of independent local fits.
"""
function compose_fit_result(x::AbstractVector, y::AbstractVector, local_fits;
                            refine::Bool = true)

    n = length(local_fits)
    baseline_est = estimate_baseline(y)

    if refine && n > 0
        # Build p0 from local fits for global refinement
        p0 = [baseline_est]
        for lf in local_fits
            append!(p0, [lf.A, lf.mu, lf.sigma])
        end

        # Build bounds
        xmin, xmax = minimum(x), maximum(x)
        ymax       = maximum(y)
        sigma_max  = (xmax - xmin) / max(1.5 * n, 1.0)
        lower = vcat([0.0],    repeat([0.0,    xmin, 0.5      ], n))
        upper = vcat([ymax],   repeat([2*ymax, xmax, sigma_max], n))

        # Clamp p0 within bounds
        p0 = clamp.(p0, lower, upper)

        model(xdata, p) = begin
            result = fill(p[1], length(xdata))
            for i in 1:n
                A, mu, sigma = p[2 + 3*(i-1)], p[3 + 3*(i-1)], p[4 + 3*(i-1)]
                result .+= _gauss.(xdata, A, mu, sigma)
            end
            return result
        end

        try
            fit    = curve_fit(model, Float64.(x), Float64.(y), p0; lower=lower, upper=upper)
            params = coef(fit)
            errors = _safe_errors(fit)

            yfit      = model(Float64.(x), params)
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
            @warn "Global refinement failed, falling back to local fits: $e"
            # Fall through to non-refined composition below
        end
    end

    # Non-refined path: compose from independent local fits
    components = [(
        A         = lf.A,
        mu        = lf.mu,
        sigma     = lf.sigma,
        A_err     = lf.A_err,
        mu_err    = lf.mu_err,
        sigma_err = lf.sigma_err,
    ) for lf in local_fits]

    yfit      = _multi_gauss_eval(Float64.(x), baseline_est, components)
    residuals = y .- yfit
    rms       = sqrt(sum(residuals.^2) / length(residuals))
    k         = 1 + 3 * n
    N         = length(x)
    aic       = N * log(rms^2) + 2 * k
    bic       = N * log(rms^2) + k * log(N)

    params = vcat([baseline_est], [[c.A, c.mu, c.sigma] for c in components]...)
    errors = vcat([NaN], [[c.A_err, c.mu_err, c.sigma_err] for c in components]...)

    return (
        params     = params,
        errors     = errors,
        baseline   = baseline_est,
        components = components,
        yfit       = yfit,
        residuals  = residuals,
        rms        = rms,
        aic        = aic,
        bic        = bic,
        converged  = all(lf.converged for lf in local_fits),
    )
end

# ---------------------------------------------------------------------------
# Main public API — the full fxcor 'm' workflow
# ---------------------------------------------------------------------------

"""
    select_and_fit(x, y; refine=true, threshold=3.0, mode=:interactive,
                   title_str="", nbin=1024) -> NamedTuple

Full fxcor-style Gaussian fitting pipeline:

1. Display the profile
2. User clicks on peaks to mark them (`:interactive` mode, like fxcor 'm')
   — or peaks are auto-detected (`:auto` mode)
3. Fit a local Gaussian to a window around each marked peak
4. Optionally refine with a global multi-Gaussian fit

# Arguments
- `x`          : phase bins or longitude values
- `y`          : intensity values
- `refine`     : if `true`, do a global refinement after local fits (default)
- `threshold`  : σ threshold for auto peak detection (used when `mode=:auto`)
- `mode`       : `:interactive` (click to select) or `:auto` (auto-detect)
- `title_str`  : title for the interactive plot
- `nbin`       : number of phase bins per period (for print_fit_summary)

# Returns
Same NamedTuple as `GaussianFit.fit_gaussians`:
`(params, errors, baseline, components, yfit, residuals, rms, aic, bic, converged)`
"""
function select_and_fit(x::AbstractVector, y::AbstractVector;
                        refine::Bool      = true,
                        threshold::Float64 = 3.0,
                        mode::Symbol       = :interactive,
                        title_str::String  = "Click on peaks to mark (right-click / Enter to finish)",
                        max_peaks::Int     = 8)

    # Phase 1: Peak detection / selection
    if mode == :interactive
        peaks = _select_peaks_interactive(x, y; title_str=title_str)
    elseif mode == :auto
        peaks = detect_peaks_auto(x, y; threshold=threshold, max_peaks=max_peaks)
    else
        error("Unknown mode: $mode. Use :interactive or :auto.")
    end

    if isempty(peaks)
        @warn "No peaks detected/selected. Returning empty fit."
        return (
            params=nothing, errors=nothing, baseline=nothing,
            components=nothing, yfit=nothing, residuals=nothing,
            rms=Inf, aic=Inf, bic=Inf, converged=false,
        )
    end

    # Phase 2: Independent local Gaussian fit for each peak
    println("\n  Fitting local Gaussians...")
    local_fits = []
    for (i, (idx, xp, amp)) in enumerate(peaks)
        lf = fit_local_gaussian(x, y, idx)
        push!(local_fits, lf)
        status = lf.converged ? "✓" : "✗"
        @printf("    Peak %d [%s]: mu=%.2f  A=%.4f  sigma=%.2f  (window=%d:%d)\n",
                i, status, lf.mu, lf.A, lf.sigma, lf.window[1], lf.window[2])
    end

    # Phase 3: Compose into standard result (with optional global refinement)
    result = compose_fit_result(x, y, local_fits; refine=refine)

    return result
end

# ---------------------------------------------------------------------------
# component_offsets — same logic as GaussianFit, reproduced here
# ---------------------------------------------------------------------------

"""
    component_offsets(fit_high, fit_low; nbin=1024) -> Vector

Compute the phase offset between matched Gaussian components of two profiles
observed at different radio frequencies.
Components are matched by sorting on mu (left-to-right).
"""
function component_offsets(fit_high, fit_low; nbin::Int = 1024)

    nh = length(fit_high.components)
    nl = length(fit_low.components)
    @assert nh == nl "Number of components must match ($nh vs $nl)"

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

# ---------------------------------------------------------------------------
# print_fit_summary — same as GaussianFit
# ---------------------------------------------------------------------------

"""
    print_fit_summary(fit; label="", nbin=1024)

Print a human-readable summary of a fit result.
"""
function print_fit_summary(fit; label::String = "", nbin::Int = 1024)
    if !fit.converged || isnothing(fit.components)
        println("  Fit did not converge.")
        return
    end
    n = length(fit.components)
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

end # module GausPeakSelect
