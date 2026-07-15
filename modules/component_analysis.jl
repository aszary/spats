module ComponentAnalysis

using Statistics
using Printf
using PyPlot
PyPlot.matplotlib.use("Tkagg")

include("GaussianFit.jl")
using .GaussianFit

# ---------------------------------------------------------------------------
# Ask user: simple or detailed analysis mode
# ---------------------------------------------------------------------------

function ask_analysis_mode()
    println("──────────────────────────────────────────")
    println("  Select analysis mode:")
    println("  1 - Simple   (standard Gaussian fitting)")
    println("  2 - Detailed (per-component precise tracking)")
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
                windows = [(max(bin_st, round(Int, c.mu - 2.5*c.sigma)),
                            min(bin_end, round(Int, c.mu + 2.5*c.sigma)),
                            c.mu)
                           for c in components]
                println("\n  Windows accepted:")
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
# Main detailed analysis
# ---------------------------------------------------------------------------

function analyse_components(nl, nh, p, outdir; max_gauss_per_window=5, min_snr=3.0)
    n_bins  = size(nl, 1)
    n_phase = size(nl, 2)

    windows = define_windows(nl, nh, p)
    isnothing(windows) && return

    n_comp = length(windows)

    freq_low  = round(Int, get(p, "freq_low",  1023.0))
    freq_high = round(Int, get(p, "freq_high", 1523.0))

    low_colors  = ["#2196F3", "#0D47A1", "#64B5F6", "#1565C0"]
    high_colors = ["#FF6F00", "#E65100", "#FFCA28", "#F57F17"]

    centers = fill(NaN, n_bins, n_comp, 2)
    c_errs  = fill(NaN, n_bins, n_comp, 2)

    offset_data = Dict{Int, NamedTuple{(:lon,:off,:err), Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}}()

    println("\n══════════════════════════════════════════")
    println("  DETAILED PER-COMPONENT ANALYSIS")
    println("  (precise center tracking per window)")
    println("══════════════════════════════════════════")

    for bin in 1:n_bins
        @printf("\n══ p3fold bin %d / %d ══\n", bin, n_bins)

        # collect per-window fits for this bin
        fits_low  = []
        fits_high = []

        for (fi, (data, freq, fits_list)) in enumerate(zip(
                (nl, nh),
                (freq_low, freq_high),
                (fits_low, fits_high)))

            @printf("=== %d MHz ===\n", freq)
            profile = data[bin, :]

            for (ci, (ws, we, _)) in enumerate(windows)
                ws >= we && continue
                x = Float64.(ws:we)
                y = profile[ws:we]

                noise   = std(y)
                snr     = noise > 0 ? maximum(y)/noise : 0.0
                snr_low = snr < min_snr

                best_fit, best_n, _ = GaussianFit.best_ngaussians(x, y; max_n=max_gauss_per_window)
                push!(fits_list, best_fit.converged ? (x=x, y=y, fit=best_fit, n=best_n) : nothing)

                if !best_fit.converged
                    @printf("  G%d: fit failed\n", ci)
                    continue
                end

                dominant  = argmax([c.A for c in best_fit.components])
                dom_c     = best_fit.components[dominant]
                # weight by A^2: dominant Gaussian dominates, others contribute a little
                w2        = [c.A^2 for c in best_fit.components]
                cen       = sum(w2[i] * best_fit.components[i].mu for i in eachindex(w2)) / sum(w2)
                cerr      = isnan(dom_c.mu_err) ? best_fit.rms : dom_c.mu_err

                centers[bin, ci, fi] = cen
                c_errs[bin,  ci, fi] = cerr

                snr_tag = snr_low ? @sprintf(" [SNR=%.1f, low]", snr) : ""
                @printf("  G%d: n_gauss=%d  center=%.2f bins (%.2f°)%s\n",
                        ci, best_n, cen, cen/n_phase*360.0, snr_tag)
            end

            # fill NaN errors using median of valid errors for this freq/bin
            valid_errs = filter(!isnan, [c_errs[bin, ci, fi] for ci in 1:n_comp])
            if !isempty(valid_errs)
                med_err = median(valid_errs)
                for ci in 1:n_comp
                    if isnan(c_errs[bin, ci, fi]) && !isnan(centers[bin, ci, fi])
                        c_errs[bin, ci, fi] = med_err
                    end
                end
            end
        end

        # offset = high - low per component (like standard analysis)
        println("  ── Offset high − low ──")
        for ci in 1:n_comp
            lo_c = centers[bin, ci, 1]; hi_c = centers[bin, ci, 2]
            el   = c_errs[bin,  ci, 1]; eh   = c_errs[bin,  ci, 2]
            (isnan(lo_c) || isnan(hi_c)) && continue
            diff     = hi_c - lo_c
            diff_deg = diff / n_phase * 360.0
            lon_deg  = (lo_c + hi_c) / 2.0 / n_phase * 360.0
            err      = (isnan(el) || isnan(eh)) ? NaN : sqrt(el^2 + eh^2)
            err_deg  = isnan(err) ? NaN : err / n_phase * 360.0
            @printf("  G%d: lon=%.2f°  %+.3f ± %.3f bins  (%+.3f° ± %.3f°)\n",
                    ci, lon_deg, diff, isnan(err) ? 0.0 : err,
                    diff_deg, isnan(err_deg) ? 0.0 : err_deg)

            if !haskey(offset_data, ci)
                offset_data[ci] = (lon=Float64[], off=Float64[], err=Float64[])
            end
            push!(offset_data[ci].lon, lon_deg)
            push!(offset_data[ci].off, diff_deg)
            push!(offset_data[ci].err, isnan(err_deg) ? 0.0 : err_deg)
        end

        # --- PyPlot per-bin window (like normal analysis) ---
        x_data = Float64.(Int(p["bin_st"]):Int(p["bin_end"]))
        y_low  = nl[bin, Int(p["bin_st"]):Int(p["bin_end"])]
        y_high = nh[bin, Int(p["bin_st"]):Int(p["bin_end"])]

        figure(figsize=(8, 6))
        suptitle("DETAILED ANALYSIS — p3fold bin $bin / $n_bins", fontsize=11, fontweight="bold")

        subplot(2, 1, 1)
        plot(x_data, y_low,  color="steelblue",  lw=1.2, alpha=0.7, label="Low $(freq_low) MHz")
        plot(x_data, y_high, color="darkorange", lw=1.2, alpha=0.7, label="High $(freq_high) MHz")

        # draw per-window fits
        for (ci, (ws, we, wref)) in enumerate(windows)
            axvspan(ws, we, alpha=0.08, color=low_colors[mod1(ci, length(low_colors))])
            axvline(wref, color=low_colors[mod1(ci, length(low_colors))], lw=0.8, ls=":")
        end
        for (ci, fd) in enumerate(fits_low)
            isnothing(fd) && continue
            col = low_colors[mod1(ci, length(low_colors))]
            plot(fd.x, fd.fit.yfit, "--", color=col, lw=1.8, label="Low G$ci fit")
            for c in fd.fit.components
                g = fd.fit.baseline .+ GaussianFit._gauss.(fd.x, c.A, c.mu, c.sigma)
                plot(fd.x, g, color=col, lw=0.8, alpha=0.5)
            end
            cen_lo = centers[bin, ci, 1]
            isnan(cen_lo) || axvline(cen_lo, color=col, lw=1.2, ls="-",
                                     label="Low G$ci center")
        end
        for (ci, fd) in enumerate(fits_high)
            isnothing(fd) && continue
            col = high_colors[mod1(ci, length(high_colors))]
            plot(fd.x, fd.fit.yfit, "--", color=col, lw=1.8, label="High G$ci fit")
            for c in fd.fit.components
                g = fd.fit.baseline .+ GaussianFit._gauss.(fd.x, c.A, c.mu, c.sigma)
                plot(fd.x, g, color=col, lw=0.8, alpha=0.5)
            end
            cen_hi = centers[bin, ci, 2]
            isnan(cen_hi) || axvline(cen_hi, color=col, lw=1.2, ls="-",
                                     label="High G$ci center")
        end
        legend(fontsize=7)
        xlim(p["bin_st"], p["bin_end"])
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
        xlim(p["bin_st"], p["bin_end"])
        xlabel("Center position (bin)")
        legend(fontsize=7)
        tight_layout()
        show()

        println("Press Enter for next bin, 'q' to quit.")
        user_input = readline(stdin; keep=false)
        close("all")
        lowercase(strip(user_input)) == "q" && break
    end

    # --- weighted mean summary ---
    println("\n══════════════════════════════════════════")
    println("=== Weighted mean offsets (high − low) ===")
    for ci in sort(collect(keys(offset_data)))
        d = offset_data[ci]
        isempty(d.err) || all(d.err .== 0.0) && continue
        w        = 1.0 ./ max.(d.err, 1e-10).^2
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

    # --- final offset plot (like normal analysis) ---
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

    # --- additional plot: center position vs p3fold bin, low and high separately ---
    bins_axis = 1:n_bins
    low_colors  = ["#2196F3", "#0D47A1", "#64B5F6", "#1565C0"]
    high_colors = ["#FF6F00", "#E65100", "#FFCA28", "#F57F17"]
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

    # --- save ---
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
