module ComponentAnalysis

using Statistics

include("GaussianFit.jl")
using .GaussianFit

# ---------------------------------------------------------------------------
# Ask user: simple or detailed analysis mode
# ---------------------------------------------------------------------------

function ask_analysis_mode()
    println("──────────────────────────────────────────")
    println("  Select analysis mode:")
    println("  1 - Simple   (current Gaussian fitting)")
    println("  2 - Detailed (per-component offset tracking)")
    println("──────────────────────────────────────────")
    print("Choice [Enter=1]: ")
    inp = strip(readline())
    return isempty(inp) || inp == "1" ? :simple : :detailed
end

# ---------------------------------------------------------------------------
# Fit N Gaussians to mean profile, let user accept or retry
# ---------------------------------------------------------------------------

function detect_components(mean_profile, p)
    bin_st = Int(p["bin_st"])
    bin_end = Int(p["bin_end"])
    x = Float64.(bin_st:bin_end)
    y = mean_profile[bin_st:bin_end]

    while true
        print("Number of components to fit [default=2]: ")
        inp = strip(readline())
        n = isempty(inp) ? 2 : parse(Int, inp)

        fit = GaussianFit.fit_gaussians(x, y, n)

        if !fit.converged
            println("  Fit did not converge, try a different N.")
            continue
        end

        println("──────────────────────────────────────────")
        println("  Detected components (bins):")
        for (i, c) in enumerate(fit.components)
            @printf("  %d: center=%.1f ± %.1f   sigma=%.1f   amp=%.4f\n",
                    i, c.mu, c.mu_err, c.sigma, c.A)
        end
        @printf("  RMS=%.5f   AIC=%.1f\n", fit.rms, fit.aic)
        println("──────────────────────────────────────────")
        print("Accept? [y=yes, n=retry with different N, q=quit]: ")
        ans = strip(readline())
        ans == "q" && return nothing
        ans == "y" && return fit.components
    end
end

# ---------------------------------------------------------------------------
# Measure center of signal within a window for one p3fold bin profile
# Uses multi-Gaussian fit; falls back to centroid if fit is poor
# ---------------------------------------------------------------------------

function measure_center(profile, win_st, win_end; n_gauss=2, min_snr=3.0)
    x = Float64.(win_st:win_end)
    y = profile[win_st:win_end]

    noise = std(y)
    snr   = noise > 0 ? maximum(y) / noise : 0.0
    snr < min_snr && return NaN

    fit = GaussianFit.fit_gaussians(x, y, n_gauss)

    if fit.converged
        # weighted center by amplitude across all Gaussians in window
        total_amp = sum(c.A for c in fit.components)
        if total_amp > 0
            return sum(c.A * c.mu for c in fit.components) / total_amp
        end
    end

    # fallback: intensity-weighted centroid
    total = sum(y)
    return total > 0 ? sum(x .* y) / total : NaN
end

# ---------------------------------------------------------------------------
# Main detailed analysis: per-component offset across p3fold bins
# ---------------------------------------------------------------------------

"""
    analyse_components(nl, nh, p, outdir; n_gauss_per_window=2, min_snr=3.0)

Detailed per-component drift analysis as alternative to simple Gaussian fitting.

1. Fits N Gaussians to the mean p3fold profile (user selects N and confirms).
2. Defines a window around each component (center ± 2.5σ).
3. For each p3fold bin × each window: fits `n_gauss_per_window` Gaussians
   and extracts the amplitude-weighted center.
4. Reports offset = center - reference_center per component per bin.
5. Saves results to `component_offsets.txt` in outdir.

Arguments:
- nl, nh  : normalised p3fold matrices (n_bins × n_phase), low and high freq
- p       : params dict (needs bin_st, bin_end)
- outdir  : output directory for saved results
"""
function analyse_components(nl, nh, p, outdir; n_gauss_per_window=2, min_snr=3.0)
    bin_st = Int(p["bin_st"])
    bin_end = Int(p["bin_end"])
    n_bins  = size(nl, 1)

    # mean profiles
    mean_low  = vec(mean(nl, dims=1))
    mean_high = vec(mean(nh, dims=1))
    mean_both = (mean_low .+ mean_high) ./ 2

    # detect components from combined mean profile
    println("\n--- Detecting components from mean profile ---")
    components = detect_components(mean_both, p)
    isnothing(components) && return

    n_comp = length(components)

    # windows: center ± 2.5σ, clamped to on-pulse region
    windows = [(max(bin_st, round(Int, c.mu - 2.5 * c.sigma)),
                min(bin_end, round(Int, c.mu + 2.5 * c.sigma)))
               for c in components]

    ref_centers = [c.mu for c in components]

    println("\n--- Measuring per-bin offsets ---")

    # offsets[bin, comp, freq]:  freq 1=low, 2=high
    offsets = fill(NaN, n_bins, n_comp, 2)

    for (fi, data) in enumerate((nl, nh))
        freq_label = fi == 1 ? "low" : "high"
        for bin in 1:n_bins
            for (ci, win) in enumerate(windows)
                win[1] >= win[2] && continue
                center = measure_center(data[bin, :], win[1], win[2];
                                        n_gauss=n_gauss_per_window, min_snr=min_snr)
                offsets[bin, ci, fi] = isnan(center) ? NaN : center - ref_centers[ci]
            end
        end
        println("  $freq_label freq done")
    end

    # print summary table
    println("\n  bin  |  " * join(["comp$i low  comp$i high" for i in 1:n_comp], "  |  "))
    for bin in 1:n_bins
        row = "$bin  |"
        for ci in 1:n_comp
            lo = isnan(offsets[bin, ci, 1]) ? "  NaN " : @sprintf("  %+.2f", offsets[bin, ci, 1])
            hi = isnan(offsets[bin, ci, 2]) ? "  NaN " : @sprintf("  %+.2f", offsets[bin, ci, 2])
            row *= lo * hi * "  |"
        end
        println(row)
    end

    # save to file
    outfile = joinpath(outdir, "component_offsets.txt")
    open(outfile, "w") do f
        header = "# bin  " * join(["comp$(i)_low  comp$(i)_high" for i in 1:n_comp], "  ")
        println(f, header)
        for bin in 1:n_bins
            row = string(bin)
            for ci in 1:n_comp
                for fi in 1:2
                    row *= "  " * (isnan(offsets[bin, ci, fi]) ? "NaN" : @sprintf("%.4f", offsets[bin, ci, fi]))
                end
            end
            println(f, row)
        end
    end
    println("\nOffsets saved to $outfile")

    return offsets
end

end  # module ComponentAnalysis
