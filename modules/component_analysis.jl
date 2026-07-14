module ComponentAnalysis

using Statistics
using Printf

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
# ASCII plot of profile + fit with component markers
# ---------------------------------------------------------------------------

function _ascii_plot(x, y, yfit, components; width=70, height=14)
    ymax = maximum(y)
    ymin = minimum(y)
    yrange = ymax - ymin
    yrange < 1e-10 && (yrange = 1.0)

    x_step = max(1, div(length(x), width))
    ys_d = y[1:x_step:end]
    ys_f = yfit[1:x_step:end]
    cols  = length(ys_d)

    # build marker line showing component positions
    marker_line = fill(' ', cols)
    for (i, c) in enumerate(components)
        col = round(Int, (c.mu - x[1]) / x_step) + 1
        col = clamp(col, 1, cols)
        marker_line[col] = string(i)[1]
    end

    println("  " * String(marker_line) * "  ← component numbers")
    for row in height:-1:1
        thresh = ymin + yrange * (row - 1) / height
        line = ""
        for col in 1:cols
            d = ys_d[col]; f = ys_f[col]
            if d >= thresh && f >= thresh;  line *= "█"
            elseif d >= thresh;             line *= "▒"
            elseif f >= thresh;             line *= "░"
            else;                           line *= " "
            end
        end
        if row == height
            @printf("%.4f │%s\n", ymax, line)
        elseif row == 1
            @printf("%.4f │%s\n", ymin, line)
        else
            @printf("       │%s\n", line)
        end
    end
    println("       └" * "─"^cols)
    println("  █=data+fit  ▒=data only  ░=fit only")
end

# ---------------------------------------------------------------------------
# Print current component table
# ---------------------------------------------------------------------------

function _print_components(components)
    println("  ┌─────┬──────────────┬──────────────┬──────────────┐")
    println("  │  #  │  center(bin) │    sigma     │   amplitude  │")
    println("  ├─────┼──────────────┼──────────────┼──────────────┤")
    for (i, c) in enumerate(components)
        @printf("  │ %3d │  %9.2f   │  %9.2f   │  %9.4f   │\n",
                i, c.mu, c.sigma, c.A)
    end
    println("  └─────┴──────────────┴──────────────┴──────────────┘")
end

# ---------------------------------------------------------------------------
# Print help for manual mode
# ---------------------------------------------------------------------------

function _print_manual_help()
    println("  ┌─────────────────────────────────────────────────────────┐")
    println("  │  COMMANDS (manual component editor)                      │")
    println("  ├─────────────────────────────────────────────────────────┤")
    println("  │  move N center sigma  – move component N to new center  │")
    println("  │                         and set its sigma               │")
    println("  │  add center sigma     – add new component               │")
    println("  │  del N                – delete component N              │")
    println("  │  show                 – redraw plot                     │")
    println("  │  fit                  – re-fit Gaussians to current pos │")
    println("  │  y                    – accept and continue             │")
    println("  │  b                    – go back, pick different bin     │")
    println("  └─────────────────────────────────────────────────────────┘")
end

# ---------------------------------------------------------------------------
# Fit Gaussians with given initial positions (manual override)
# ---------------------------------------------------------------------------

function _fit_with_positions(x, y, centers_sigmas)
    n = length(centers_sigmas)
    baseline_est = minimum(y)
    p0 = [baseline_est]
    ymax = maximum(y)
    for (mu, sigma) in centers_sigmas
        # estimate amplitude from data near center
        idx = argmin(abs.(x .- mu))
        A_est = max(y[idx] - baseline_est, 0.01 * ymax)
        append!(p0, [A_est, mu, sigma])
    end
    lo, hi = GaussianFit._default_bounds(x, y, n)
    p0 = clamp.(p0, lo, hi)
    return GaussianFit.fit_gaussians(x, y, n; p0=p0, lower=lo, upper=hi)
end

# ---------------------------------------------------------------------------
# Interactive window definition: pick bin or manual mode
# Returns list of (win_st, win_end, ref_center) tuples, or nothing
# ---------------------------------------------------------------------------

function define_windows(nl, nh, p)
    bin_st    = Int(p["bin_st"])
    bin_end   = Int(p["bin_end"])
    n_bins    = size(nl, 1)
    mean_low  = vec(mean(nl, dims=1))
    mean_high = vec(mean(nh, dims=1))
    mean_both = (mean_low .+ mean_high) ./ 2

    x = Float64.(bin_st:bin_end)

    current_bin = 0   # 0 = mean profile

    while true
        # --- choose reference ---
        println("\n──────────────────────────────────────────")
        println("  Choose reference profile for component windows:")
        println("  Enter = mean profile")
        @printf("  1-%d   = specific p3fold bin\n", n_bins)
        println("──────────────────────────────────────────")
        print("Choice: ")
        inp = strip(readline())

        if isempty(inp)
            current_bin = 0
            y = mean_both[bin_st:bin_end]
            println("  Using mean profile.")
        else
            b = tryparse(Int, inp)
            if isnothing(b) || b < 1 || b > n_bins
                println("  Invalid bin number.")
                continue
            end
            current_bin = b
            y_low  = nl[b, bin_st:bin_end]
            y_high = nh[b, bin_st:bin_end]
            y = (y_low .+ y_high) ./ 2
            @printf("  Using bin %d.\n", b)
        end

        # --- initial fit ---
        print("Number of components [default=2]: ")
        n_inp = strip(readline())
        n = isempty(n_inp) ? 2 : parse(Int, n_inp)

        fit = GaussianFit.fit_gaussians(x, y, n)
        if !fit.converged
            println("  Fit did not converge. Try different bin or N.")
            continue
        end

        components = collect(fit.components)
        sort!(components, by = c -> c.mu)

        # --- review / manual edit loop ---
        while true
            println()
            _ascii_plot(x, y, fit.yfit, components)
            println()
            _print_components(components)
            @printf("  RMS=%.5f   AIC=%.1f\n", fit.rms, fit.aic)
            println()
            println("  [y] Accept   [m] Manual edit   [b] Pick different bin/profile")
            print("  Choice: ")
            ans = strip(readline())

            if ans == "y"
                # build windows
                sort!(components, by = c -> c.mu)
                windows = [(max(bin_st, round(Int, c.mu - 2.5*c.sigma)),
                            min(bin_end, round(Int, c.mu + 2.5*c.sigma)),
                            c.mu)
                           for c in components]
                println("\n  Windows accepted:")
                for (i, (ws, we, wc)) in enumerate(windows)
                    @printf("    G%d: bins %d–%d  (center=%.1f)\n", i, ws, we, wc)
                end
                return windows

            elseif ans == "b"
                break   # back to bin selection

            elseif ans == "m"
                _print_manual_help()
                # manual edit sub-loop
                cs = [(c.mu, c.sigma) for c in components]
                while true
                    print("  cmd> ")
                    cmd = strip(readline())
                    parts = split(cmd)
                    isempty(parts) && continue

                    if parts[1] == "y"
                        fit2 = _fit_with_positions(x, y, cs)
                        if fit2.converged
                            fit = fit2
                            components = collect(fit.components)
                        end
                        break

                    elseif parts[1] == "b"
                        break

                    elseif parts[1] == "show"
                        fit2 = _fit_with_positions(x, y, cs)
                        if fit2.converged
                            fit = fit2
                            components = collect(fit.components)
                        end
                        _ascii_plot(x, y, fit.yfit, components)
                        _print_components(components)

                    elseif parts[1] == "fit"
                        fit2 = _fit_with_positions(x, y, cs)
                        if fit2.converged
                            fit = fit2
                            components = collect(fit.components)
                            cs = [(c.mu, c.sigma) for c in components]
                            _ascii_plot(x, y, fit.yfit, components)
                            _print_components(components)
                        else
                            println("  Fit failed, try adjusting positions.")
                        end

                    elseif parts[1] == "move" && length(parts) == 4
                        idx = tryparse(Int, parts[2])
                        mu  = tryparse(Float64, parts[3])
                        sig = tryparse(Float64, parts[4])
                        if isnothing(idx) || isnothing(mu) || isnothing(sig) ||
                           idx < 1 || idx > length(cs)
                            println("  Usage: move N center sigma")
                        else
                            cs[idx] = (mu, sig)
                            fit2 = _fit_with_positions(x, y, cs)
                            if fit2.converged
                                fit = fit2
                                components = collect(fit.components)
                                sort!(components, by = c -> c.mu)
                                cs = [(c.mu, c.sigma) for c in components]
                            end
                            _ascii_plot(x, y, fit.yfit, components)
                            _print_components(components)
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
                                components = collect(fit.components)
                                sort!(components, by = c -> c.mu)
                                cs = [(c.mu, c.sigma) for c in components]
                            end
                            _ascii_plot(x, y, fit.yfit, components)
                            _print_components(components)
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
                                components = collect(fit.components)
                                sort!(components, by = c -> c.mu)
                                cs = [(c.mu, c.sigma) for c in components]
                            end
                            _ascii_plot(x, y, fit.yfit, components)
                            _print_components(components)
                        end

                    else
                        _print_manual_help()
                    end
                end  # manual sub-loop

                # after manual: re-enter review loop with updated fit
                fit2 = _fit_with_positions(x, y, cs)
                if fit2.converged
                    fit = fit2
                    components = collect(fit.components)
                    sort!(components, by = c -> c.mu)
                end
            end
        end  # review loop
    end  # bin selection loop
end

# ---------------------------------------------------------------------------
# Per-bin detailed analysis
# ---------------------------------------------------------------------------

function analyse_components(nl, nh, p, outdir; max_gauss_per_window=5, min_snr=3.0)
    n_bins  = size(nl, 1)
    n_phase = size(nl, 2)

    windows = define_windows(nl, nh, p)
    isnothing(windows) && return

    n_comp = length(windows)

    freq_low  = round(Int, get(p, "freq_low",  1023.0))
    freq_high = round(Int, get(p, "freq_high", 1523.0))

    centers = fill(NaN, n_bins, n_comp, 2)
    c_errs  = fill(NaN, n_bins, n_comp, 2)

    println("\n══════════════════════════════════════════")
    println("  Detailed per-bin analysis")
    println("══════════════════════════════════════════")

    for bin in 1:n_bins
        @printf("\n══ p3fold bin %d / %d ══\n", bin, n_bins)

        for (fi, (data, freq)) in enumerate(zip((nl, nh), (freq_low, freq_high)))
            @printf("=== %d MHz ===\n", freq)
            profile = data[bin, :]

            for (ci, (ws, we, wref)) in enumerate(windows)
                ws >= we && continue
                x = Float64.(ws:we)
                y = profile[ws:we]

                noise = std(y)
                snr   = noise > 0 ? maximum(y)/noise : 0.0
                if snr < min_snr
                    @printf("  G%d: SNR=%.1f — too low, skipped\n", ci, snr)
                    continue
                end

                best_fit, best_n, _ = GaussianFit.best_ngaussians(x, y; max_n=max_gauss_per_window)
                if !best_fit.converged
                    @printf("  G%d: fit failed\n", ci)
                    continue
                end

                total_amp = sum(c.A for c in best_fit.components)
                cen  = total_amp > 0 ?
                       sum(c.A * c.mu     for c in best_fit.components) / total_amp : NaN
                cerr = total_amp > 0 ?
                       sqrt(sum((c.A * c.mu_err)^2 for c in best_fit.components)) / total_amp : NaN

                centers[bin, ci, fi] = cen
                c_errs[bin,  ci, fi] = cerr

                off_bins = cen - wref
                off_deg  = off_bins / n_phase * 360.0
                @printf("  G%d: n_gauss=%d  center=%.2f  offset=%+.3f bins (%+.3f°)\n",
                        ci, best_n, cen, off_bins, off_deg)
            end
        end

        # frequency comparison
        any_valid = false
        for ci in 1:n_comp
            lo = centers[bin, ci, 1]; hi = centers[bin, ci, 2]
            isnan(lo) || isnan(hi) && continue
            any_valid = true
        end
        if any_valid
            println("  ── Δ(high − low) ──")
            for ci in 1:n_comp
                lo = centers[bin, ci, 1]; hi = centers[bin, ci, 2]
                el = c_errs[bin, ci, 1]; eh = c_errs[bin, ci, 2]
                (isnan(lo) || isnan(hi)) && continue
                diff     = hi - lo
                diff_deg = diff / n_phase * 360.0
                err      = sqrt(el^2 + eh^2)
                err_deg  = err / n_phase * 360.0
                @printf("  G%d: %+.3f ± %.3f bins  (%+.3f° ± %.3f°)\n",
                        ci, diff, err, diff_deg, err_deg)
            end
        end

        print("  Press Enter for next bin, 'q' to quit: ")
        strip(readline()) == "q" && break
    end

    # --- weighted mean summary ---
    println("\n══════════════════════════════════════════")
    println("=== Weighted mean offsets (high − low) ===")
    for ci in 1:n_comp
        diffs = [centers[b,ci,2] - centers[b,ci,1] for b in 1:n_bins]
        errs  = [sqrt(c_errs[b,ci,1]^2 + c_errs[b,ci,2]^2) for b in 1:n_bins]
        valid = .!isnan.(diffs) .& .!isnan.(errs) .& (errs .> 0)
        if !any(valid)
            println("  G$ci: no valid measurements")
            continue
        end
        w     = 1.0 ./ errs[valid].^2
        wm    = sum(w .* diffs[valid]) / sum(w)
        we    = 1.0 / sqrt(sum(w))
        wm_d  = wm / n_phase * 360.0
        we_d  = we / n_phase * 360.0
        @printf("  G%d: offset = %+.4f ± %.4f bins  (%+.4f° ± %.4f°)  n=%d\n",
                ci, wm, we, wm_d, we_d, sum(valid))
    end

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
    println("\nSaved to $outfile")
    return centers
end

end  # module ComponentAnalysis
