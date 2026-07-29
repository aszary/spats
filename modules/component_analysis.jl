module ComponentAnalysis

using Statistics
using Printf
using PyPlot
PyPlot.matplotlib.use("Tkagg")

include("GaussianFit.jl")
include("modules/component_analysis.jl")
using .GaussianFit

# ---------------------------------------------------------------------------
# Ask user: simple or detailed analysis mode
# ---------------------------------------------------------------------------

function ask_analysis_mode()
    println("──────────────────────────────────────────")
    println("  Select analysis mode:")
    println("  1 - Simple   (N-Gaussian on full onpulse = reference)")
    println("  2 - Detailed (same Gauss μ, windows + tracking)")
    println("──────────────────────────────────────────")
    print("Choice [Enter=1]: ")
    inp = strip(readline())
    return isempty(inp) || inp == "1" ? :simple : :detailed
end

# ---------------------------------------------------------------------------
# Plot profile + fit with component markers (PyPlot window)
# ---------------------------------------------------------------------------

function _plot_components(x, y, fit, components, title_str)
    figure(figsize=(8, 5))

    low_colors = ["#2196F3", "#4CAF50", "#9C27B0", "#E65100", "#FF5722"]

    plot(x, y, color="steelblue", lw=1.5, alpha=0.8, label="Data")
    if fit.converged
        plot(x, fit.yfit, "--", color="black", lw=1.8, label="Total fit")
        for (j, c) in enumerate(components)
            col = low_colors[mod1(j, length(low_colors))]
            gauss_y = fit.baseline .+ GaussianFit._gauss.(x, c.A, c.mu, c.sigma)
            plot(x, gauss_y, color=col, lw=1.8, alpha=0.9, label="G$j (mu=$(round(c.mu,digits=1)))")
            axvline(c.mu, color=col, lw=0.8, ls=":")
        end
    end
    legend(fontsize=8)
    xlabel("Bin")
    ylabel("Intensity")
    title(title_str)
    tight_layout()
    show()
end

# ---------------------------------------------------------------------------
# Interactive window definition using PyPlot + text commands
# Returns list of (win_st, win_end, ref_center) or nothing
# ---------------------------------------------------------------------------

function define_windows(nl, nh, p)
    bin_st  = Int(p["bin_st"])
    bin_end = Int(p["bin_end"])
    n_bins  = size(nl, 1)

    mean_low  = vec(mean(nl, dims=1))
    mean_high = vec(mean(nh, dims=1))
    mean_both = (mean_low .+ mean_high) ./ 2

    x = Float64.(bin_st:bin_end)

    while true
        # --- choose reference ---
        println("\n──────────────────────────────────────────")
        println("  Choose reference profile for component windows:")
        println("  Enter = mean profile")
        @printf("  1-%d   = specific p3fold bin\n", n_bins)
        println("──────────────────────────────────────────")
        print("Choice: ")
        inp = strip(readline())

        local y
        local title_ref
        if isempty(inp)
            y = mean_both[bin_st:bin_end]
            title_ref = "Mean profile"
        else
            b = tryparse(Int, inp)
            if isnothing(b) || b < 1 || b > n_bins
                println("  Invalid bin number.")
                continue
            end
            y = ((nl[b, bin_st:bin_end] .+ nh[b, bin_st:bin_end]) ./ 2)
            title_ref = "p3fold bin $b"
        end

        # --- initial fit ---
        print("Number of components [default=2]: ")
        n_inp = strip(readline())
        n = isempty(n_inp) ? 2 : parse(Int, n_inp)

        fit = GaussianFit.fit_gaussians(x, y, n)
        components = fit.converged ? sort(collect(fit.components), by=c->c.mu) : []

        if !fit.converged
            println("  Fit did not converge.")
            close("all")
            continue
        end

        _plot_components(x, y, fit, components, "$title_ref — component fit")

        # --- print table ---
        println()
        println("  ┌─────┬──────────────┬──────────────┬──────────────┐")
        println("  │  #  │  center(bin) │    sigma     │   amplitude  │")
        println("  ├─────┼──────────────┼──────────────┼──────────────┤")
        for (i, c) in enumerate(components)
            @printf("  │ %3d │  %9.2f   │  %9.2f   │  %9.4f   │\n",
                    i, c.mu, c.sigma, c.A)
        end
        println("  └─────┴──────────────┴──────────────┴──────────────┘")
        @printf("  RMS=%.5f   AIC=%.1f\n", fit.rms, fit.aic)

        # --- commands ---
        println()
        println("  COMMANDS:")
        println("    y                  – accept windows")
        println("    b                  – pick different bin/profile")
        println("    move N center sigma – move component N")
        println("    add  center sigma  – add component")
        println("    del  N             – delete component N")
        println("    refit              – refit Gaussians to current positions")
        println("    show               – redraw plot")

        cs = [(c.mu, c.sigma) for c in components]

        while true
            print("  cmd> ")
            cmd = strip(readline())
            parts = split(cmd)
            isempty(parts) && continue

            if parts[1] == "y"
                close("all")
                sort!(components, by=c->c.mu)
                n_c = length(components)
                mus = [c.mu for c in components]
                windows = Vector{Tuple{Int,Int,Float64}}(undef, n_c)
                # precompute inner midpoints between adjacent components
                mids = [(mus[i] + mus[i+1]) / 2.0 for i in 1:n_c-1]
                for i in 1:n_c
                    # outer: ±2.5σ so Gaussian wings are not truncated (better μ)
                    # inner: midpoint between adjacent centres (no overlap)
                    ws = i == 1 ?
                         round(Int, components[i].mu - 2.5*components[i].sigma) :
                         round(Int, mids[i-1])
                    we = i == n_c ?
                         round(Int, components[i].mu + 2.5*components[i].sigma) :
                         round(Int, mids[i])
                    windows[i] = (max(bin_st, ws), min(bin_end, we), mus[i])
                end
                println("\n  Windows accepted (midpoint boundaries, outer ±2.5σ):")
                for (i, (ws, we, wc)) in enumerate(windows)
                    @printf("    G%d: bins %d–%d  (center=%.1f)\n", i, ws, we, wc)
                end
                return windows

            elseif parts[1] == "b"
                close("all")
                break

            elseif parts[1] == "show"
                close("all")
                fit2 = _fit_with_positions(x, y, cs)
                if fit2.converged
                    fit = fit2
                    components = sort(collect(fit.components), by=c->c.mu)
                end
                _plot_components(x, y, fit, components, "$title_ref — component fit")

            elseif parts[1] == "refit"
                fit2 = _fit_with_positions(x, y, cs)
                if fit2.converged
                    fit = fit2
                    components = sort(collect(fit.components), by=c->c.mu)
                    cs = [(c.mu, c.sigma) for c in components]
                    close("all")
                    _plot_components(x, y, fit, components, "$title_ref — component fit")
                    println()
                    println("  ┌─────┬──────────────┬──────────────┬──────────────┐")
                    println("  │  #  │  center(bin) │    sigma     │   amplitude  │")
                    println("  ├─────┼──────────────┼──────────────┼──────────────┤")
                    for (i, c) in enumerate(components)
                        @printf("  │ %3d │  %9.2f   │  %9.2f   │  %9.4f   │\n",
                                i, c.mu, c.sigma, c.A)
                    end
                    println("  └─────┴──────────────┴──────────────┴──────────────┘")
                else
                    println("  Fit failed.")
                end

            elseif parts[1] == "move" && length(parts) == 4
                idx = tryparse(Int, parts[2])
                mu  = tryparse(Float64, parts[3])
                sig = tryparse(Float64, parts[4])
                if isnothing(idx) || isnothing(mu) || isnothing(sig) || idx < 1 || idx > length(cs)
                    println("  Usage: move N center sigma")
                else
                    cs[idx] = (mu, sig)
                    fit2 = _fit_with_positions(x, y, cs)
                    if fit2.converged
                        fit = fit2
                        components = sort(collect(fit.components), by=c->c.mu)
                        cs = [(c.mu, c.sigma) for c in components]
                    end
                    close("all")
                    _plot_components(x, y, fit, components, "$title_ref — component fit")
                end

            elseif parts[1] == "add" && length(parts) == 3
                mu  = tryparse(Float64, parts[2])
                sig = tryparse(Float64, parts[3])
                if isnothing(mu) || isnothing(sig)
                    println("  Usage: add center sigma")
                else
                    push!(cs, (mu, sig))
                    fit2 = _fit_with_positions(x, y, cs)
                    if fit2.converged
                        fit = fit2
                        components = sort(collect(fit.components), by=c->c.mu)
                        cs = [(c.mu, c.sigma) for c in components]
                    end
                    close("all")
                    _plot_components(x, y, fit, components, "$title_ref — component fit")
                end

            elseif parts[1] == "del" && length(parts) == 2
                idx = tryparse(Int, parts[2])
                if isnothing(idx) || idx < 1 || idx > length(cs)
                    println("  Usage: del N")
                elseif length(cs) == 1
                    println("  Cannot delete last component.")
                else
                    deleteat!(cs, idx)
                    fit2 = _fit_with_positions(x, y, cs)
                    if fit2.converged
                        fit = fit2
                        components = sort(collect(fit.components), by=c->c.mu)
                        cs = [(c.mu, c.sigma) for c in components]
                    end
                    close("all")
                    _plot_components(x, y, fit, components, "$title_ref — component fit")
                end

            else
                println("  Unknown command. Type 'y' to accept, 'b' to go back.")
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Fit Gaussians with given initial positions
# ---------------------------------------------------------------------------

function _fit_with_positions(x, y, centers_sigmas)
    n = length(centers_sigmas)
    baseline_est = minimum(y)
    ymax = maximum(y)
    p0 = [baseline_est]
    for (mu, sigma) in centers_sigmas
        idx = argmin(abs.(x .- mu))
        A_est = max(y[idx] - baseline_est, 0.01 * ymax)
        append!(p0, [A_est, mu, sigma])
    end
    lo, hi = GaussianFit._default_bounds(x, y, n)
    p0 = clamp.(p0, lo, hi)
    return GaussianFit.fit_gaussians(x, y, n; p0=p0, lower=lo, upper=hi)
end

# ---------------------------------------------------------------------------
# Measure component centre in a window.
# Primary: single-Gaussian μ (same estimator as Simple for isolated components).
# Fallback: intensity-weighted centroid when the Gaussian fit is unusable.
# ---------------------------------------------------------------------------

function _offpulse_noise(profile, bin_st, bin_end)
    n = length(profile)
    lo = max(1, Int(bin_st))
    hi = min(n, Int(bin_end))
    off = Float64[]
    lo > 1 && append!(off, profile[1:lo-1])
    hi < n && append!(off, profile[hi+1:n])
    if length(off) < 4
        on = profile[lo:hi]
        n_edge = max(2, round(Int, 0.15 * length(on)))
        off = vcat(on[1:n_edge], on[end-n_edge+1:end])
    end
    noise = std(off)
    baseline = median(off)
    return baseline, (noise > 0 ? noise : std(profile))
end

function _centroid_center(x, y, baseline, noise)
    ys = y .- baseline
    peak = maximum(ys)
    thresh = max(0.2 * peak, 2.0 * noise)
    mask = ys .>= thresh
    if sum(mask) < 2
        mask = ys .> 0
    end
    sum(mask) < 2 && return (NaN, NaN)
    xm = x[mask]
    ym = ys[mask]
    s = sum(ym)
    s <= 0 && return (NaN, NaN)
    cen = sum(ym .* xm) / s
    cerr = noise * sqrt(sum((xm .- cen).^2)) / s
    return (cen, cerr)
end

"""
    _measure_center_windowed(x, y, ref_mu; profile, bin_st, bin_end, min_snr)

Single-Gaussian (or centroid) centre inside one window — fallback path.
"""
function _measure_center_windowed(x, y, ref_mu; profile, bin_st, bin_end, min_snr=3.0)
    baseline, noise = _offpulse_noise(profile, bin_st, bin_end)
    peak = maximum(y) - baseline
    snr = noise > 0 ? peak / noise : 0.0

    ymax = maximum(y)
    peak_idx = argmax(y)
    mu0 = Float64(x[peak_idx])
    if abs(mu0 - ref_mu) > 0.35 * (last(x) - first(x))
        mu0 = Float64(ref_mu)
        peak_idx = argmin(abs.(x .- mu0))
    end
    A0 = max(y[peak_idx] - baseline, 0.01 * (ymax - baseline + eps()))
    sig0 = max((last(x) - first(x)) / 6.0, 1.5)
    p0 = [max(baseline, 0.0), A0, mu0, sig0]
    lo, hi = GaussianFit._default_bounds(x, y, 1)
    p0 = clamp.(p0, lo, hi)

    fit = GaussianFit.fit_gaussians(x, y, 1; p0=p0, lower=lo, upper=hi)
    if fit.converged
        c = fit.components[1]
        xmin, xmax = first(x), last(x)
        if c.mu >= xmin && c.mu <= xmax && c.A > 0 && isfinite(c.mu)
            mu_err = isfinite(c.mu_err) && c.mu_err > 0 ? c.mu_err :
                     noise * sig0 / max(c.A, eps())
            return (center=c.mu, err=mu_err, method=:gauss, fit=fit,
                    snr=snr, keep=(snr >= min_snr))
        end
    end

    cen, cerr = _centroid_center(x, y, baseline, noise)
    return (center=cen, err=cerr, method=:centroid,
            fit=(fit.converged ? fit : nothing),
            snr=snr, keep=(snr >= min_snr && isfinite(cen)))
end

"""
    _measure_centers_global(profile, windows, bin_st, bin_end; min_snr)

Fit N Gaussians on the full on-pulse, seeded by window centres (Detailed path —
not the same auto-p0 as Simple). Assign each fitted component to the nearest
window. For clear separated components this usually agrees with Simple within
errors, while still being a window-guided method.
"""
function _measure_centers_global(profile, windows, bin_st, bin_end; min_snr=3.0)
    n_comp = length(windows)
    x = Float64.(bin_st:bin_end)
    y = profile[bin_st:bin_end]
    baseline, noise = _offpulse_noise(profile, bin_st, bin_end)

    # Seed from window reference centres / half-widths
    p0 = [max(baseline, 0.0)]
    for (ws, we, wref) in windows
        idx = clamp(round(Int, wref) - bin_st + 1, 1, length(y))
        A0 = max(y[idx] - baseline, 0.01 * (maximum(y) - baseline + eps()))
        sig0 = max((we - ws) / 4.0, 1.5)
        append!(p0, [A0, Float64(wref), sig0])
    end
    lo, hi = GaussianFit._default_bounds(x, y, n_comp)
    p0 = clamp.(p0, lo, hi)

    fit = GaussianFit.fit_gaussians(x, y, n_comp; p0=p0, lower=lo, upper=hi)
    fit.converged || return nothing

    comps = sort(collect(fit.components), by=c -> c.mu)
    length(comps) == n_comp || return nothing

    # Greedy nearest-window assignment (windows already sorted by μ)
    used = falses(n_comp)
    assigned = Vector{Any}(undef, n_comp)
    fill!(assigned, nothing)
    for (ci, (ws, we, wref)) in enumerate(windows)
        best_j = 0
        best_d = Inf
        for (j, c) in enumerate(comps)
            used[j] && continue
            d = abs(c.mu - wref)
            in_win = (c.mu >= ws - 2) && (c.mu <= we + 2)
            if in_win && d < best_d
                best_d = d
                best_j = j
            end
        end
        best_j == 0 && return nothing
        used[best_j] = true
        c = comps[best_j]
        peak = c.A
        snr = noise > 0 ? peak / noise : 0.0
        mu_err = isfinite(c.mu_err) && c.mu_err > 0 ? c.mu_err :
                 noise * c.sigma / max(c.A, eps())
        local_fit = (
            converged = true,
            baseline = fit.baseline,
            components = [c],
            yfit = fit.baseline .+ GaussianFit._gauss.(Float64.(ws:we), c.A, c.mu, c.sigma),
        )
        assigned[ci] = (center=c.mu, err=mu_err, method=:gauss_global,
                        fit=local_fit, snr=snr, keep=(snr >= min_snr),
                        x=Float64.(ws:we), y=profile[ws:we])
    end
    any(isnothing, assigned) && return nothing
    return assigned
end

# ---------------------------------------------------------------------------
# Main detailed analysis
# ---------------------------------------------------------------------------

function analyse_components(nl, nh, p, outdir; min_snr=3.0)
    n_bins  = size(nl, 1)
    n_phase = size(nl, 2)
    bin_st  = Int(p["bin_st"])
    bin_end = Int(p["bin_end"])

    windows = define_windows(nl, nh, p)
    isnothing(windows) && return

    n_comp = length(windows)

    freq_low  = round(Int, get(p, "freq_low",  1023.0))
    freq_high = round(Int, get(p, "freq_high", 1523.0))

    low_colors  = ["#2196F3", "#0D47A1", "#64B5F6", "#1565C0", "#4CAF50", "#9C27B0", "#00897B", "#5D4037"]
    high_colors = ["#FF6F00", "#E65100", "#FFCA28", "#F57F17", "#FF7043", "#D84315", "#FFAB00", "#BF360C"]

    centers = fill(NaN, n_bins, n_comp, 2)
    c_errs  = fill(NaN, n_bins, n_comp, 2)
    keep_freq = fill(false, n_bins, n_comp, 2)

    offset_data = Dict{Int, NamedTuple{(:lon,:off,:err), Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}}()

    println("\n══════════════════════════════════════════")
    println("  DETAILED PER-COMPONENT ANALYSIS")
    println("  Primary: window-seeded N-Gaussian on onpulse (gauss_global)")
    println("  Fallback: 1-Gaussian / centroid per window")
    println("══════════════════════════════════════════")

    for bin in 1:n_bins
        @printf("\n══ p3fold bin %d / %d ══\n", bin, n_bins)

        fits_low  = Any[nothing for _ in 1:n_comp]
        fits_high = Any[nothing for _ in 1:n_comp]

        for (fi, (data, freq, fits_list)) in enumerate(zip(
                (nl, nh),
                (freq_low, freq_high),
                (fits_low, fits_high)))

            @printf("=== %d MHz ===\n", freq)
            profile = data[bin, :]

            global_ms = _measure_centers_global(profile, windows, bin_st, bin_end;
                                                min_snr=min_snr)

            for (ci, (ws, we, wref)) in enumerate(windows)
                if ws >= we
                    @printf("  G%d: empty window — skipped\n", ci)
                    continue
                end

                local m
                if !isnothing(global_ms)
                    m = global_ms[ci]
                else
                    x = Float64.(ws:we)
                    y = profile[ws:we]
                    m = _measure_center_windowed(x, y, wref;
                                                 profile=profile, bin_st=bin_st,
                                                 bin_end=bin_end, min_snr=min_snr)
                end

                if !isfinite(m.center)
                    @printf("  G%d: no valid centre\n", ci)
                    continue
                end

                centers[bin, ci, fi] = m.center
                c_errs[bin,  ci, fi] = m.err
                keep_freq[bin, ci, fi] = m.keep
                if hasproperty(m, :x)
                    fits_list[ci] = (x=m.x, y=m.y, fit=m.fit, n=1, center=m.center)
                elseif !isnothing(m.fit)
                    fits_list[ci] = (x=Float64.(ws:we), y=profile[ws:we],
                                     fit=m.fit, n=1, center=m.center)
                end

                snr_tag = m.snr < min_snr ? @sprintf(" [SNR=%.1f, low]", m.snr) : ""
                meth = String(m.method)
                @printf("  G%d: %s=%.2f bins (%.2f°)  err=%.3f bins%s\n",
                        ci, meth, m.center, m.center/n_phase*360.0, m.err, snr_tag)
            end
        end

        println("  ── Offset high − low ──")
        pending = NamedTuple{(:ci, :lon, :off, :err), Tuple{Int,Float64,Float64,Float64}}[]
        for ci in 1:n_comp
            lo_c = centers[bin, ci, 1]; hi_c = centers[bin, ci, 2]
            el   = c_errs[bin,  ci, 1]; eh   = c_errs[bin,  ci, 2]
            (isnan(lo_c) || isnan(hi_c)) && continue
            diff     = hi_c - lo_c
            diff_deg = diff / n_phase * 360.0
            lon_deg  = (lo_c + hi_c) / 2.0 / n_phase * 360.0
            err      = (isnan(el) || isnan(eh) || !(el > 0) || !(eh > 0)) ? NaN : sqrt(el^2 + eh^2)
            err_deg  = isnan(err) ? NaN : err / n_phase * 360.0
            keep = keep_freq[bin, ci, 1] && keep_freq[bin, ci, 2] && isfinite(err_deg) && err_deg > 0
            keep_tag = keep ? "" : " [excluded from mean]"
            @printf("  G%d: lon=%.2f°  %+.3f ± %.3f bins  (%+.3f° ± %.3f°)%s\n",
                    ci, lon_deg, diff, isnan(err) ? 0.0 : err,
                    diff_deg, isnan(err_deg) ? 0.0 : err_deg, keep_tag)

            if keep
                push!(pending, (ci=ci, lon=lon_deg, off=diff_deg, err=err_deg))
            end
        end

        # --- PyPlot per-bin ---
        x_data = Float64.(bin_st:bin_end)
        y_low  = nl[bin, bin_st:bin_end]
        y_high = nh[bin, bin_st:bin_end]

        figure(figsize=(8, 6))
        suptitle("DETAILED ANALYSIS — p3fold bin $bin / $n_bins", fontsize=11, fontweight="bold")

        subplot(2, 1, 1)
        plot(x_data, y_low,  color="steelblue",  lw=1.2, alpha=0.7, label="Low $(freq_low) MHz")
        plot(x_data, y_high, color="darkorange", lw=1.2, alpha=0.7, label="High $(freq_high) MHz")

        for (ci, (ws, we, wref)) in enumerate(windows)
            axvspan(ws, we, alpha=0.08, color=low_colors[mod1(ci, length(low_colors))])
            axvline(wref, color=low_colors[mod1(ci, length(low_colors))], lw=0.8, ls=":")
        end
        x_plot = collect(x_data)
        for (ci, fd) in enumerate(fits_low)
            col = low_colors[mod1(ci, length(low_colors))]
            if !isnothing(fd)
                plot(fd.x, fd.fit.yfit, "--", color=col, lw=1.8, label="Low G$ci fit")
                for c in fd.fit.components
                    g = fd.fit.baseline .+ GaussianFit._gauss.(x_plot, c.A, c.mu, c.sigma)
                    plot(x_plot, g, color=col, lw=0.8, alpha=0.5)
                end
            end
            # Always draw measured centre (matches stored result)
            cen_lo = centers[bin, ci, 1]
            isnan(cen_lo) || axvline(cen_lo, color=col, lw=1.4, ls="-",
                                     label="Low G$ci center")
        end
        for (ci, fd) in enumerate(fits_high)
            col = high_colors[mod1(ci, length(high_colors))]
            if !isnothing(fd)
                plot(fd.x, fd.fit.yfit, "--", color=col, lw=1.8, label="High G$ci fit")
                for c in fd.fit.components
                    g = fd.fit.baseline .+ GaussianFit._gauss.(x_plot, c.A, c.mu, c.sigma)
                    plot(x_plot, g, color=col, lw=0.8, alpha=0.5)
                end
            end
            cen_hi = centers[bin, ci, 2]
            isnan(cen_hi) || axvline(cen_hi, color=col, lw=1.4, ls="-",
                                     label="High G$ci center")
        end
        legend(fontsize=7)
        xlim(bin_st, bin_end)
        ylabel("Intensity")

        subplot(2, 1, 2)
        for ci in 1:n_comp
            lo_c = centers[bin, ci, 1]; hi_c = centers[bin, ci, 2]
            isnan(lo_c) || scatter([lo_c], [1.0],
                color=low_colors[mod1(ci, length(low_colors))],
                marker="o", s=80, label="Low G$ci", zorder=3)
            isnan(hi_c) || scatter([hi_c], [2.0],
                color=high_colors[mod1(ci, length(high_colors))],
                marker="s", s=80, label="High G$ci", zorder=3)
        end
        yticks([1.0, 2.0], ["Low", "High"])
        xlim(bin_st, bin_end)
        xlabel("Center position (bin)")
        legend(fontsize=7)
        tight_layout()
        show()

        println("Enter = accept bin,  s = skip bin,  q = quit")
        user_input = lowercase(strip(readline(stdin; keep=false)))
        close("all")

        if user_input == "q"
            println("Exiting analysis (current bin not saved).")
            # clear current bin from arrays so it won't appear on drift plot / file
            centers[bin, :, :] .= NaN
            c_errs[bin, :, :]  .= NaN
            break
        elseif user_input == "s" || user_input == "n" || user_input == "skip"
            println("  Bin $bin skipped — not included in mean / final plots.")
            centers[bin, :, :] .= NaN
            c_errs[bin, :, :]  .= NaN
            continue
        else
            for o in pending
                if !haskey(offset_data, o.ci)
                    offset_data[o.ci] = (lon=Float64[], off=Float64[], err=Float64[])
                end
                push!(offset_data[o.ci].lon, o.lon)
                push!(offset_data[o.ci].off, o.off)
                push!(offset_data[o.ci].err, o.err)
            end
        end
    end

    # --- weighted mean summary ---
    println("\n══════════════════════════════════════════")
    println("=== Weighted mean offsets (high − low) ===")
    for ci in sort(collect(keys(offset_data)))
        d = offset_data[ci]
        if isempty(d.err) || all(e -> !(isfinite(e) && e > 0), d.err)
            continue
        end
        w        = 1.0 ./ d.err.^2
        n        = length(d.off)
        mu       = sum(w .* d.off) / sum(w)
        sig_int  = 1.0 / sqrt(sum(w))
        chi2     = sum(w .* (d.off .- mu).^2)
        dof      = n - 1
        chi2_red = dof > 0 ? chi2/dof : NaN
        sig_ext  = sig_int * sqrt(max(1.0, chi2_red))
        @printf("G%d: offset = %+.4f° ± %.4f°  (n=%d, χ²/dof = %.2f, σ_ext = %.4f°)\n",
                ci, mu, sig_int, n, chi2_red, sig_ext)
    end
    println()

    # --- final offset plot ---
    if !isempty(offset_data)
        colors = ["#2196F3", "#E65100", "#4CAF50", "#9C27B0"]
        figure(figsize=(8, 5))
        for ci in sort(collect(keys(offset_data)))
            d = offset_data[ci]
            errorbar(d.lon, d.off, yerr=d.err,
                     fmt="o", capsize=4, lw=1.2,
                     color=colors[mod1(ci, length(colors))],
                     label="G$ci")
        end
        axhline(0.0, color="gray", lw=0.8, ls="--")
        minorticks_on()
        xlabel("Longitude (°)")
        ylabel("Offset high−low (°)")
        title("Detailed analysis — component offsets")
        legend()
        tight_layout()
        show()
        println("Press Enter to close final plot.")
        readline(stdin; keep=false)
        close("all")
    end

    # --- center vs p3fold bin ---
    bins_axis = 1:n_bins
    figure(figsize=(8, 5))
    for ci in 1:n_comp
        lo_cen = [centers[b, ci, 1] for b in 1:n_bins]
        hi_cen = [centers[b, ci, 2] for b in 1:n_bins]
        lo_err = [c_errs[b,  ci, 1] for b in 1:n_bins]
        hi_err = [c_errs[b,  ci, 2] for b in 1:n_bins]
        errorbar(collect(bins_axis), lo_cen, yerr=lo_err,
                 fmt="o-", capsize=4, lw=1.2,
                 color=low_colors[mod1(ci, length(low_colors))],
                 label="G$ci low")
        errorbar(collect(bins_axis), hi_cen, yerr=hi_err,
                 fmt="s--", capsize=4, lw=1.2,
                 color=high_colors[mod1(ci, length(high_colors))],
                 label="G$ci high")
    end
    minorticks_on()
    xlabel("p3fold bin")
    ylabel("Center position (bin)")
    title("Component center vs p3fold bin — low and high")
    legend(fontsize=8)
    tight_layout()
    show()
    println("Press Enter to close drift plot.")
    readline(stdin; keep=false)
    close("all")

    outfile = joinpath(outdir, "component_offsets.txt")
    open(outfile, "w") do f
        cols = join(["G$(i)_low  G$(i)_low_err  G$(i)_high  G$(i)_high_err"
                     for i in 1:n_comp], "  ")
        println(f, "# bin  $cols")
        for bin in 1:n_bins
            row = string(bin)
            for ci in 1:n_comp, fi in 1:2
                o = isnan(centers[bin,ci,fi]) ? "NaN" : @sprintf("%.4f", centers[bin,ci,fi])
                e = isnan(c_errs[bin,ci,fi])  ? "NaN" : @sprintf("%.4f", c_errs[bin,ci,fi])
                row *= "  $o  $e"
            end
            println(f, row)
        end
    end
    println("Saved to $outfile")
    return centers
end

end  # module ComponentAnalysis
