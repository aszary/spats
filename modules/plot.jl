module Plot
    using LsqFit
    using FFTW
    using JLD2
    using Statistics
    using StatsBase
    using PyPlot
    using PyCall
    using Printf
    @pyimport matplotlib.patches as patch
    
    PyPlot.matplotlib.use("Tkagg") 
    #PyPlot.matplotlib.use("qt5agg") # DOES NOT WORK!! 
    using Peaks
    using Glob
    using SmoothingSplines
    using DataFrames, GLM
    using Distributions

    #using RCall
    using HypothesisTests
    using LinearAlgebra
    using Distances

    include("tools.jl")
    #include("pyrmodule.jl")
    include("functions.jl")
    include("GaussianFit.jl")


    function average(data, outdir; start=1, number=100, bin_st=nothing, bin_end=nothing, name_mod="0", show_=false)
        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1,bin_st:bin_end]
        average = Tools.average_profile(da)

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 2.362205), frameon=true)  # 8cm x 6 cm
        subplots_adjust(left=0.16, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)

        minorticks_on()
        plot(longitude, average, c="grey")
        #yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        xlabel("longitude \$(^\\circ)\$")
        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/average_$name_mod.pdf")
        savefig("$outdir/average_$name_mod.pdf")
        if show_ == true
            show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
        #clf()
    end


    function single(data, outdir; start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME", show_=false, repeat_num=nothing)
        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end

        da = data[start:start+number-1,bin_st:bin_end]

        if repeat_num !== nothing
            # repeat data
            da = repeat(da, repeat_num)

        end

        average = Tools.average_profile(da)
        intensity, pulses = Tools.intensity_pulses(da)
        intensity .-= minimum(intensity)
        intensity ./= maximum(intensity)

        pulses .+= start - 1  # julia

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))
        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 4.33071), frameon=true)  # 8cm x 11cm
        subplots_adjust(left=0.16, bottom=0.09, right=0.99, top=0.99, wspace=0., hspace=0.)

        subplot2grid((5, 3), (0, 0), rowspan=4)
        minorticks_on()
        plot(intensity, pulses, c="grey")
        ylim(pulses[1]-0.5, pulses[end]+0.5)
        xticks([0.5, 1.0])
        xlim(1.1, -0.1)
        xlabel("intensity")
        ylabel("Pulse number")

        subplot2grid((5, 3), (0, 1), rowspan=4, colspan=2)
        imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da))
        #axvline(x=563, lw=2)
        tick_params(labelleft=false, labelbottom=false)

        subplot2grid((5, 3), (4, 1), colspan=2)
        minorticks_on()
        plot(longitude, average, c="grey")
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        xlabel("longitude \$(^\\circ)\$")
        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_single.pdf")
        savefig("$outdir/$(name_mod)_single.pdf")
        println("$outdir/$(name_mod)_single.png")
        savefig("$outdir/$(name_mod)_single.png")
        if show_ == true
            PyPlot.show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
        #clf()
    end


    function lrfs_obsolete(data, outdir; start=1, number=nothing, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0", change_fftphase=true, show_=false)

        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1,bin_st:bin_end]
        average = Tools.average_profile(da)
        lrfs, intensity, freq, peak = Tools.lrfs(da)
        println("\tpeak freq $(freq[peak]) ")
        println("\tpeak P3 $(1/freq[peak])")
        phase_ = rad2deg.(angle.(view(lrfs, peak, :)))  # fft phase variation  # view used! lrfs[peak, :] -> copies data
        # skip freq = 0 and normalize intensity to 1
        inten = intensity[2:end]
        inten .-= minimum(inten)
        inten ./= maximum(inten)
        fre = freq[2:end]
        pars, errs = Tools.fit_gaussian(fre, inten; μ=freq[peak-1])  # skip zero freq
        fr = pars[2]
        frer = pars[3]  # errs[2] # yeap

        println("\tFrequency (gaussian fit): $fr, P3: $(1/fr)")
        println("\tFrequency error (gaussian fit): $frer, P3 error: $(1/fr - 1/(fr +frer))")  # TODO err ok?

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        # fixing phase continuity, estimating phase change
        dl = longitude[2]-longitude[1]
        dphase = zeros(length(phase_)-1)
        println("Longitude resolution = ", dl, " (deg.)")
        if change_fftphase == true
            for i in 1:length(phase_)-1
                global dp
                changed = true
                while changed
                    changed = false
                    dp = abs(phase_[i+1] - phase_[i])
                    dpp = abs(phase_[i+1]+360. - phase_[i])
                    dpm = abs(phase_[i+1]-360. - phase_[i])
                    if dpp < dp
                        phase_[i+1] += 360
                        changed = true
                    end
                    if dpm < dp
                        phase_[i+1] -= 360
                        changed = true
                    end
                end
                dphase[i] = dp
                #println("$i $(longitude[i])  $dp")
            end
        end
        # longitude vs. FFT phase analysis

        #=
        lin, incs = Tools.find_inclination(longitude, phase_)
        for (i,inc) in enumerate(incs)
            println("($i) Δ FFT phase / Δ longitude ", round(inc[2], digits=2))
            println("($i) sigma Δ FFT phase / Δ longitude ", round(inc[3], digits=2))
        end
        =#

        lin, inc = Tools.find_inclination2(longitude, phase_)
        println(" Δ FFT phase / Δ longitude ", round(inc[1][2], digits=2))
        println("sigma Δ FFT phase / Δ longitude ", round(inc[1][3], digits=2))

        rc("font", size=7.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 4.33071))  # 8cm x 11cm
        subplots_adjust(left=0.17, bottom=0.08, right=0.90, top=0.92, wspace=0., hspace=0.)

        ax = subplot2grid((5, 3), (0, 1), colspan=2)
        minorticks_on()
        #plot(longitude, phase_, c="grey")
        scatter(longitude, phase_, marker=".", c="grey", s=3.)
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        for l in lin
            plot(l[1], l[2])
        end
        #ylim(-500, 0)
        #xlim(-1, 2)
        xlabel("longitude \$(^\\circ)\$")
        ylabel("FFT phase \$(^\\circ)\$")
        if change_fftphase == true
            ax2 = ax.twinx()
            minorticks_on()
            plot(longitude[7:end-9], dphase[6:end-9])
            ax2.xaxis.set_label_position("top")
        end
        tick_params(labeltop=true, labelbottom=false, which="both", bottom=false, top=true)

        subplot2grid((5, 3), (1, 0), rowspan=3)
        minorticks_on()
        plot(inten, fre, c="grey")
        plot(Tools.gauss(fre, pars), fre, c="red", ls=":", lw=0.3)
        ylim(freq[1], freq[end])
        xticks([0.5, 1.0])
        xlim(1.1, -0.1)
        xlabel("intensity")
        ylabel("frequency \$(1/P)\$")

        subplot2grid((5, 3), (1, 1), rowspan=3, colspan=2)
        imshow(abs.(lrfs), origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(abs.(lrfs)))
        #println(size(lrfs))
        #return

        tick_params(labelleft=false, labelbottom=false)

        subplot2grid((5, 3), (4, 1), colspan=2)
        minorticks_on()
        plot(longitude, average, c="grey")
        axvline(x=0., ls=":", c="black")
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        xlabel("longitude \$(^\\circ)\$")
        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_lrfs_obsolete.pdf")
        savefig("$outdir/$(name_mod)_lrfs_obsolete.pdf")
        savefig("$outdir/$(name_mod)_lrf_obsolete.svg")
        if show_ == true
            show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    function p3fold(data, outdir; start=1, number=nothing, repeat_num=4, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0", show_=false)

        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end

        da = data[start:start+number-1,bin_st:bin_end]
        average = Tools.average_profile(da)
        intensity, pulses = Tools.intensity_pulses(da)
        intensity .-= minimum(intensity)
        intensity ./= maximum(intensity)

        pulses .+= start

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        # repeat data
        da = repeat(da, repeat_num)
        intensity = repeat(intensity, repeat_num)
        pulses = collect(1:length(intensity))
        le = length(pulses)
        ticks = [floor(Int, le /4), floor(Int, le /2), floor(Int, le *3 / 4)]
        fracs = [repeat_num / 4, repeat_num/2, repeat_num * 3 / 4]
        ti = ["$(fracs[i])\$P_3\$" for i in 1:3]

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)


        figure(figsize=(3.14961, 4.33071))  # 8cm x 11cm
        subplots_adjust(left=0.16, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)

        subplot2grid((5, 3), (0, 0), rowspan=4)
        minorticks_on()
        plot(intensity, pulses, c="grey")
        ylim(pulses[1], pulses[end])
        xlim(1.1, -0.1)
        xlabel("intensity")
        xticks([0.5, 1.0])
        yticks(ticks, ti)
        #ylabel("Pulse number")

        subplot2grid((5, 3), (0, 1), rowspan=4, colspan=2)
        imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da))
        tick_params(labelleft=false, labelbottom=false)

        subplot2grid((5, 3), (4, 1), colspan=2)
        minorticks_on()
        plot(longitude, average, c="grey")
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        xlabel("longitude \$(^\\circ)\$")
        #tick_params(labeltop=false, labelbottom=true)
        savefig("$outdir/$(name_mod)_p3fold.pdf")
        println("$outdir/$(name_mod)_p3fold.pdf")
        savefig("$outdir/$(name_mod)_p3fold.png")
        println("$outdir/$(name_mod)_p3fold.png")
        if show_ == true
            show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    """
    Diagnostic figure for `P3FoldViterbi.fold`: the Viterbi-refined p3-fold
    next to the constant-P3 (naive, `Tools.p3fold`) fold of the same data,
    with a panel underneath showing the per-pulse instantaneous P3 implied
    by the Viterbi phase path (`p3_per_pulse`) against the nominal P3 used
    as input (dashed line).
    """
    function p3fold_compare(folded_viterbi, folded_const, p3_per_pulse, p3_nominal, outdir;
                             bin_st=nothing, bin_end=nothing, cmap="viridis", darkness=1.0,
                             repeat_num=4, name_mod="0", show_=false, label="Viterbi refine",
                             p3_per_pulse_err=nothing, intensity=nothing, p3_ylim=nothing)

        _, bins = size(folded_viterbi)
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end

        dv = repeat(folded_viterbi[:, bin_st:bin_end], repeat_num)
        dc = repeat(folded_const[:, bin_st:bin_end], repeat_num)

        # each panel still scaled to its own max (the two folds aren't
        # necessarily on the same absolute intensity scale, e.g. psrsalsa's
        # pfold output vs Tools.p3fold's raw unnormalized sum — a shared
        # vmax washed `norefine` out entirely). Same convention as
        # `single`/`lrfs_obsolete`/`p3fold`/`twodfs`/`lrfs`: vmin is left at
        # its automatic (true min) default, only vmax is pulled down by
        # `darkness`.

        le = size(dv, 1)
        ticks = [floor(Int, le / 4), floor(Int, le / 2), floor(Int, le * 3 / 4)]
        fracs = [repeat_num / 4, repeat_num / 2, repeat_num * 3 / 4]
        ti = ["$(fracs[i])\$P_3\$" for i in 1:3]

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        # intensity panel is a thin strip, not a full extra row — keeps the figure
        # close to its original height instead of growing proportionally per panel
        nrows = intensity === nothing ? 3 : 8
        img_rowspan = intensity === nothing ? 2 : 4
        p3_rowspan  = intensity === nothing ? 1 : 3
        fig_height = intensity === nothing ? 6.29921 : 7.3
        figure(figsize=(6.29921, fig_height))

        subplots_adjust(left=0.1, bottom=0.07, right=0.99, top=0.96, wspace=0.15, hspace=0.35)

        subplot2grid((nrows, 2), (0, 0), rowspan=img_rowspan)
        imshow(dv, origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(dv))
        yticks(ticks, ti)
        title(label)

        subplot2grid((nrows, 2), (0, 1), rowspan=img_rowspan)
        imshow(dc, origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(dc))
        tick_params(labelleft=false)
        title("constant \$P_3\$")

        pulses_x = 1:length(p3_per_pulse)

        if intensity !== nothing
            subplot2grid((nrows, 2), (img_rowspan, 0), rowspan=1, colspan=2)
            minorticks_on()
            plot(pulses_x, intensity, c="grey", lw=0.6)
            xlim(1, length(p3_per_pulse))
            tick_params(labelbottom=false)
            ylabel("intensity")
        end

        subplot2grid((nrows, 2), (nrows - p3_rowspan, 0), rowspan=p3_rowspan, colspan=2)
        minorticks_on()
        if p3_per_pulse_err !== nothing
            fill_between(pulses_x, p3_per_pulse .- p3_per_pulse_err, p3_per_pulse .+ p3_per_pulse_err,
                         color="grey", alpha=0.3, linewidth=0)
        end
        plot(pulses_x, p3_per_pulse, c="grey", lw=0.8)
        axhline(p3_nominal, c="red", ls="--", lw=0.8)
        xlim(1, length(p3_per_pulse))
        if p3_ylim !== nothing
            ylim(p3_ylim)
        end
        xlabel("pulse number")
        ylabel("\$P_3\$ (P)")

        savefig("$outdir/$(name_mod)_p3fold_compare.pdf")
        println("$outdir/$(name_mod)_p3fold_compare.pdf")
        savefig("$outdir/$(name_mod)_p3fold_compare.png")
        println("$outdir/$(name_mod)_p3fold_compare.png")
        if show_ == true
            show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    """
    Plot for `P3FoldViterbi.scan_lowpass_cutoff` — helps pick `lowpass_cutoff`
    for `coherent_fold`/`p3fold_coherent` per pulsar (see that function's
    docstring for the reasoning).

    Top panel: consistency (independent-subband agreement, higher = more
    trustworthy) vs cutoff, log-x. Bottom panel: p3_std (track noise level)
    vs cutoff, same x-axis, for context on the noise/responsiveness trade-off.
    A vertical dashed line at `chosen` (if given) marks the cutoff you decide
    to use in `p3fold_coherent`.
    """
    function lowpass_cutoff_scan(results, outdir; name_mod="pulsar", chosen=nothing,
                                 chosen_metrics=nothing, refine_trace=nothing,
                                 r_threshold=nothing, show_=false)
        cutoffs     = [r.cutoff for r in results]
        consistency = [r.consistency for r in results]
        p3_std      = [r.p3_std for r in results]
        signal_std  = [hasproperty(r, :signal_std) ? r.signal_std : NaN for r in results]

        rc("font", size=8.5, family="sans-serif")
        rc("axes", linewidth=0.7, edgecolor="#444444")
        fig = figure(figsize=(7, 8))
        fig.patch.set_facecolor("white")

        subplot(3, 1, 1)
        grid(true, which="both", alpha=0.2, lw=0.4)
        minorticks_on()
        semilogx(cutoffs, consistency, "o-", color="steelblue", lw=1.3, ms=5,
                 mec="white", mew=0.4, zorder=3, label="siatka (zgrubny skan)")
        for (x, y) in zip(cutoffs, consistency)
            annotate(@sprintf("%.3f", y), (x, y), textcoords="offset points",
                     xytext=(0, 6), ha="center", fontsize=6.5, color="#333333")
        end
        # below this line the combined P3(t) track is more noise than signal
        # (reliability R < 0.5) — such cutoffs don't compete on signal_std
        if r_threshold !== nothing
            axhline(r_threshold, color="gray", ls=":", lw=1.0, zorder=2,
                    label=@sprintf("próg wiarygodności (R=0.5, r=%.3f)", r_threshold))
            legend(fontsize=6.5, loc="best", framealpha=0.9)
        end
        chosen !== nothing && axvspan(chosen*0.92, chosen*1.08, color="red", alpha=0.07, zorder=1)
        chosen !== nothing && axvline(chosen, color="red", ls="--", lw=1.0, alpha=0.6, zorder=2)
        if chosen_metrics !== nothing
            scatter([chosen_metrics.cutoff], [chosen_metrics.consistency],
                    marker="*", s=220, color="red", edgecolor="darkred", linewidth=0.6,
                    zorder=5, label="wybrany (dopracowany)")
            annotate(@sprintf("%.3f", chosen_metrics.consistency),
                     (chosen_metrics.cutoff, chosen_metrics.consistency),
                     textcoords="offset points", xytext=(0, -13), ha="center",
                     fontsize=7, color="red", fontweight="bold")
            legend(fontsize=6.5, loc="best", framealpha=0.9)
        end
        ylabel("subband consistency (r)", fontsize=8.5)
        title("lowpass_cutoff scan" * (isempty(name_mod) ? "" : "  ($name_mod)"),
              fontsize=10, fontweight="bold")
        tick_params(labelbottom=false)

        subplot(3, 1, 2)
        grid(true, which="both", alpha=0.2, lw=0.4)
        minorticks_on()
        loglog(cutoffs, p3_std, "o-", color="darkorange", lw=1.3, ms=5,
               mec="white", mew=0.4, zorder=3)
        for (x, y) in zip(cutoffs, p3_std)
            annotate(@sprintf("%.2g", y), (x, y), textcoords="offset points",
                     xytext=(0, 6), ha="center", fontsize=6.5, color="#333333")
        end
        chosen !== nothing && axvspan(chosen*0.92, chosen*1.08, color="red", alpha=0.07, zorder=1)
        chosen !== nothing && axvline(chosen, color="red", ls="--", lw=1.0, alpha=0.6, zorder=2)
        if chosen_metrics !== nothing
            scatter([chosen_metrics.cutoff], [chosen_metrics.p3_std],
                    marker="*", s=220, color="red", edgecolor="darkred", linewidth=0.6, zorder=5)
        end
        ylabel("std(P3 per pulse)", fontsize=8.5)
        tick_params(labelbottom=false)

        subplot(3, 1, 3)
        grid(true, which="both", alpha=0.2, lw=0.4)
        minorticks_on()
        semilogx(cutoffs, signal_std, "o-", color="seagreen", lw=1.3, ms=5,
                 mec="white", mew=0.4, zorder=3)
        for (x, y) in zip(cutoffs, signal_std)
            annotate(@sprintf("%.2g", y), (x, y), textcoords="offset points",
                     xytext=(0, 6), ha="center", fontsize=6.5, color="#333333")
        end
        # golden-section refinement points: show how the objective actually
        # behaves in the fine neighbourhood the chosen star came from,
        # instead of the star appearing disconnected from the coarse grid
        if refine_trace !== nothing && !isempty(refine_trace)
            tx = [t[1] for t in refine_trace]
            ty = [t[2] for t in refine_trace]
            scatter(tx, ty, marker="o", s=18, color="orange", edgecolor="#7a4a00",
                    linewidth=0.4, alpha=0.85, zorder=4, label="punkty dopracowania")
        end
        chosen !== nothing && axvspan(chosen*0.92, chosen*1.08, color="red", alpha=0.07, zorder=1)
        chosen !== nothing && axvline(chosen, color="red", ls="--", lw=1.0, alpha=0.6, zorder=2)
        if chosen_metrics !== nothing
            scatter([chosen_metrics.cutoff], [chosen_metrics.signal_std],
                    marker="*", s=220, color="red", edgecolor="darkred", linewidth=0.6, zorder=5)
            annotate(@sprintf("%.2g", chosen_metrics.signal_std),
                     (chosen_metrics.cutoff, chosen_metrics.signal_std),
                     textcoords="offset points", xytext=(0, -13), ha="center",
                     fontsize=7, color="red", fontweight="bold")
        end
        if refine_trace !== nothing && !isempty(refine_trace)
            legend(fontsize=6.5, loc="best", framealpha=0.9)
        end
        xlabel("lowpass_cutoff (cycles/pulse)", fontsize=8.5)
        ylabel("signal_std\n(p3_std·√consistency)", fontsize=8.5)

        tight_layout()
        savefig("$outdir/$(name_mod)_lowpass_cutoff_scan.pdf")
        println("$outdir/$(name_mod)_lowpass_cutoff_scan.pdf")
        savefig("$outdir/$(name_mod)_lowpass_cutoff_scan.png")
        println("$outdir/$(name_mod)_lowpass_cutoff_scan.png")
        if show_ == true
            show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    """
    For each cutoff scanned by `P3FoldViterbi.scan_lowpass_cutoff`, overlays
    the actual `p3_per_pulse(pulse)` track it produced — not just the
    summary std shown in `lowpass_cutoff_scan`, but the shape of the
    wobble itself, so you can see directly how much smoothing each cutoff
    applies and whether looser cutoffs are adding real fine structure or
    just noise wiggle.
    """
    function p3_tracks_vs_cutoff(results, outdir; name_mod="pulsar", highlight_cutoff=nothing,
                                 p3_nominal=nothing, show_=false)
        palette = ["#08519c", "#3182bd", "#6baed6", "#74c476",
                   "#fd8d3c", "#e6550d", "#a50f15", "#756bb1"]

        rc("font", size=8.5, family="sans-serif")
        rc("axes", linewidth=0.7, edgecolor="#444444")
        fig = figure(figsize=(8.5, 5.5))
        fig.patch.set_facecolor("white")
        grid(true, alpha=0.2, lw=0.4)
        minorticks_on()

        # explicit baseline so it's unambiguous what "no real variability"
        # looks like on this plot — previously some pulsars happened to have
        # a near-flat track sitting at the nominal value (which could be
        # mistaken for a deliberate "no cutoff" reference curve) while others
        # didn't, with no labelled ground truth either way
        if p3_nominal !== nothing
            axhline(p3_nominal, color="black", ls=":", lw=1.2, alpha=0.6, zorder=1,
                    label=@sprintf("nominalne P3 = %.4g", p3_nominal))
        end

        # draw the non-highlighted tracks first, muted, so the highlighted
        # one (if any) reads clearly on top instead of fighting for attention
        highlight_r = nothing
        for (i, r) in enumerate(results)
            hasproperty(r, :p3_per_pulse) || continue
            is_hl = highlight_cutoff !== nothing && isapprox(r.cutoff, highlight_cutoff; rtol=1e-6)
            if is_hl
                highlight_r = r
                continue
            end
            col = palette[mod1(i, length(palette))]
            plot(1:length(r.p3_per_pulse), r.p3_per_pulse, color=col, lw=1.0, alpha=0.5,
                 zorder=2, label=@sprintf("%.3g", r.cutoff))
        end
        if highlight_r !== nothing
            plot(1:length(highlight_r.p3_per_pulse), highlight_r.p3_per_pulse,
                 color="red", lw=2.6, alpha=0.95, zorder=5,
                 label=@sprintf("%.3g  ★ WYBRANY", highlight_r.cutoff))
        end
        xlabel("pulse number", fontsize=8.5)
        ylabel("\$P_3\$ per pulse (P)", fontsize=8.5)
        subtitle = highlight_r !== nothing ?
            @sprintf("wybrany cutoff = %.4g (cykle/impuls) — czerwona linia", highlight_r.cutoff) :
            ""
        title("P3 variability vs lowpass_cutoff" * (isempty(name_mod) ? "" : "  ($name_mod)") *
              (isempty(subtitle) ? "" : "\n$subtitle"),
              fontsize=10, fontweight="bold")
        if highlight_r !== nothing
            ax = gca()
            ax.text(0.015, 0.97, @sprintf("cutoff = %.4g", highlight_r.cutoff),
                    transform=ax.transAxes, fontsize=9, fontweight="bold", color="red",
                    va="top", ha="left",
                    bbox=Dict("boxstyle"=>"round,pad=0.35", "facecolor"=>"white",
                              "edgecolor"=>"red", "alpha"=>0.9))
        end
        legend(title="cutoff", fontsize=6.5, ncol=2, loc="best", framealpha=0.9)
        tight_layout()
        savefig("$outdir/$(name_mod)_p3_tracks_vs_cutoff.pdf")
        println("$outdir/$(name_mod)_p3_tracks_vs_cutoff.pdf")
        savefig("$outdir/$(name_mod)_p3_tracks_vs_cutoff.png")
        println("$outdir/$(name_mod)_p3_tracks_vs_cutoff.png")
        if show_ == true
            show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    """
    Diagnostic for "why do different cutoffs give such different P3(t)
    results?" — a 2D map with cutoff on the y-axis (log scale) and pulse
    number on the x-axis, colour = P3 value. A feature (bright/dark stripe)
    that runs vertically across *many* cutoff rows at the same pulse range
    is present regardless of how much smoothing is applied — strong evidence
    it's real. A feature confined to one or two rows (cutoffs) is a filter
    artifact specific to that amount of smoothing, not a genuine signal.

    `results` should be the same `valid` (non-rejected) list passed to
    `p3_tracks_vs_cutoff` — each entry needs `cutoff` and `p3_per_pulse`.
    """
    function p3_cutoff_heatmap(results, outdir; p3_nominal=nothing, name_mod="pulsar", show_=false)
        sorted = sort(collect(filter(r -> hasproperty(r, :p3_per_pulse), results)), by = r -> r.cutoff)
        isempty(sorted) && return nothing
        N = length(sorted[1].p3_per_pulse)
        M = zeros(Float64, length(sorted), N)
        for (i, r) in enumerate(sorted)
            M[i, :] = r.p3_per_pulse
        end
        cutoffs = [r.cutoff for r in sorted]

        rc("font", size=8.5, family="sans-serif")
        fig = figure(figsize=(8.5, 5))
        fig.patch.set_facecolor("white")

        # NaN = zapped pulses (gaps) — colour range from finite values only
        Mf = filter(isfinite, M)
        center = p3_nominal !== nothing ? p3_nominal : (isempty(Mf) ? 0.0 : median(Mf))
        spread = isempty(Mf) ? 1.0 : maximum(abs.(Mf .- center))
        spread = spread > 0 ? spread : 1.0

        im = imshow(M, aspect="auto", origin="lower", cmap="RdBu_r",
                    vmin=center - spread, vmax=center + spread,
                    extent=(1, N, 1, length(cutoffs)))
        cb = colorbar(im)
        cb.set_label("\$P_3\$ per pulse (P)" * (p3_nominal !== nothing ?
                     @sprintf("  (biały ≈ nominalne %.3g)", p3_nominal) : ""), fontsize=8)

        yticks(1:2:length(cutoffs), [@sprintf("%.3g", c) for c in cutoffs[1:2:end]], fontsize=6.5)
        xlabel("pulse number", fontsize=8.5)
        ylabel("lowpass_cutoff (cycles/pulse)", fontsize=8.5)
        title("Mapa P3: cutoff × czas" * (isempty(name_mod) ? "" : "  ($name_mod)") *
              "\npionowy pasek przez wiele wierszy = prawdziwy sygnał; pojedynczy wiersz = artefakt filtra",
              fontsize=9, fontweight="bold")
        tight_layout()
        savefig("$outdir/$(name_mod)_p3_cutoff_heatmap.pdf")
        println("$outdir/$(name_mod)_p3_cutoff_heatmap.pdf")
        savefig("$outdir/$(name_mod)_p3_cutoff_heatmap.png")
        println("$outdir/$(name_mod)_p3_cutoff_heatmap.png")
        if show_ == true
            show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    """
    Diagnostic for "is this feature on P3(t) real, or a filter artifact?" —
    overlays the `n_groups` independent longitude-subband P3(t) tracks
    (thin, muted colours) against the combined/main track (bold black).
    A feature (peak, dip) that shows up in most/all independent subbands is
    reproducible — likely real wobble, since each subband carries its own
    independent noise. A feature visible in only one or two subbands is more
    likely noise or a filter transient local to that subband's data, not
    genuine P3 variation.

    `group_p3` is `coherent_fold_jackknife`'s `group_p3` field (n_groups × N
    matrix); `main_p3` its `p3_per_pulse`.
    """
    function subband_consistency(group_p3, main_p3, outdir; p3_nominal=nothing,
                                 name_mod="pulsar", show_=false)
        n_groups, N = size(group_p3)
        palette = ["#08519c", "#3182bd", "#6baed6", "#74c476",
                   "#fd8d3c", "#e6550d", "#a50f15", "#756bb1"]

        rc("font", size=8.5, family="sans-serif")
        rc("axes", linewidth=0.7, edgecolor="#444444")
        fig = figure(figsize=(8.5, 5.5))
        fig.patch.set_facecolor("white")
        grid(true, alpha=0.2, lw=0.4)
        minorticks_on()

        if p3_nominal !== nothing
            axhline(p3_nominal, color="gray", ls=":", lw=1.0, alpha=0.6, zorder=1)
        end
        for g in 1:n_groups
            col = palette[mod1(g, length(palette))]
            plot(1:N, @view(group_p3[g, :]), color=col, lw=0.8, alpha=0.5,
                 zorder=2, label="podpasmo $g")
        end
        plot(1:length(main_p3), main_p3, color="black", lw=2.0, alpha=0.9,
             zorder=5, label="połączony wynik")

        xlabel("pulse number", fontsize=8.5)
        ylabel("\$P_3\$ per pulse (P)", fontsize=8.5)
        title("Spójność między niezależnymi podpasmami" *
              (isempty(name_mod) ? "" : "  ($name_mod)") *
              "\ncienkie linie = $(n_groups) niezależnych fragmentów danych, czarna = wynik połączony",
              fontsize=9.5, fontweight="bold")
        legend(fontsize=6, ncol=3, loc="best", framealpha=0.9)
        tight_layout()
        savefig("$outdir/$(name_mod)_subband_consistency.pdf")
        println("$outdir/$(name_mod)_subband_consistency.pdf")
        savefig("$outdir/$(name_mod)_subband_consistency.png")
        println("$outdir/$(name_mod)_subband_consistency.png")
        if show_ == true
            show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    # S/N of the average profile from the given (good) pulses
    function profile_snr(data, rows, bin_st, bin_end)
        isempty(rows) && return NaN
        prof = vec(mean(data[rows, :], dims=1))
        nb = length(prof)
        off = vcat(1:bin_st-1, bin_end+1:nb)
        length(off) < 3 && return NaN
        s = std(prof[off])
        s > 0 || return NaN
        return sum(prof[bin_st:bin_end] .- mean(prof[off])) / (s * sqrt(bin_end - bin_st + 1))
    end


    close_all() = close("all")


    """
    Open a saved PNG in an external image viewer (eog / display / feh /
    xdg-open, whichever exists) — independent of matplotlib's GUI backend,
    which doesn't always show windows over VNC. Returns the viewer
    process (or `nothing` if none could be started).
    """
    function open_png(path)
        if endswith(lowercase(path), ".pdf")
            for viewer in ("xpdf", "evince", "okular", "xdg-open")
                isnothing(Sys.which(viewer)) && continue
                cmd = viewer == "xpdf" ? `xpdf -z page $path` : `$viewer $path`
                return run(pipeline(cmd; stdout=devnull, stderr=devnull); wait=false)
            end
            println("Brak przeglądarki PDF — otwórz ręcznie: $path")
            return nothing
        end
        for viewer in ("eog", "display", "feh", "xdg-open")
            if !isnothing(Sys.which(viewer))
                return run(pipeline(`$viewer $path`; stdout=devnull, stderr=devnull); wait=false)
            end
        end
        # PDF saved next to the PNG -> xpdf (available on kopernik)
        pdf = replace(path, r"\.png$" => ".pdf")
        if pdf != path && isfile(pdf) && !isnothing(Sys.which("xpdf"))
            return run(pipeline(`xpdf -z page $pdf`; stdout=devnull, stderr=devnull); wait=false)
        end
        # fallback: separate python process (same python as PyPlot) showing
        # the PNG in a TkAgg window — independent of the Julia event loop
        py = PyCall.python
        if !isempty(py) && isfile(py)
            code = "import matplotlib; matplotlib.use('TkAgg'); import matplotlib.pyplot as plt; " *
                   "import matplotlib.image as mi, sys; im = mi.imread(sys.argv[1]); " *
                   "import tkinter; r = tkinter.Tk(); sw, sh = r.winfo_screenwidth(), r.winfo_screenheight(); " *
                   "r.destroy(); h, w = im.shape[:2]; s = min(1.0, 0.95 * sw / w, 0.85 * sh / h); " *
                   "f = plt.figure(figsize=(s * w / 100, s * h / 100), dpi=100); " *
                   "f.add_axes([0, 0, 1, 1]).imshow(im); plt.axis('off'); " *
                   "f.canvas.manager.set_window_title(sys.argv[1]); plt.show()"
            return run(pipeline(`$py -c $code $path`; stdout=devnull, stderr=devnull); wait=false)
        end
        println("Brak przeglądarki obrazów (eog/display/feh/xdg-open) — otwórz ręcznie: $path")
        return nothing
    end

    """
    Show a saved PNG (`open_png`), wait for Enter (if `wait_`), then close
    the viewer window so windows don't pile up. Returns the viewer process
    when not waiting (caller closes it with `close_png`).
    """
    function show_png(path; prompt="Enter = następny pulsar: ", wait_=true)
        proc = open_png(path)
        wait_ || return proc
        print(prompt)
        readline(stdin; keep=false)
        close_png(proc)
        return nothing
    end

    close_png(proc) = (isnothing(proc) || !process_running(proc) || kill(proc); nothing)


    """
    What the automatic zapper cuts and why (`Data.auto_zap_core` result `a`).

    Left: the data as they are. Right: the same data, every cut pulse
    covered by a coloured band across the full width, colour = reason:
      blue    null / no signal (stretch with too little signal)
      red     noisy stretch (stretch with too much noise)
      orange  single junk pulse (noise jump, stripe, spike)
      grey    margin around a cut region / short gap filled in
      black   empty pulse
    A thin strip right of the right panel repeats the reasons, so short
    cuts are visible even when zoomed out. Title: counts and S/N change.
    """
    function autozap_view(data, bin_st, bin_end, a, outdir; name_mod="pulsar", title_="", show_=false, wait_=true,
                          snr_pazi=nothing)
        N, nb = size(data)
        w = bin_end - bin_st + 1
        lo, hi = max(1, bin_st - w ÷ 2), min(nb, bin_end + w ÷ 2)
        win = data[:, lo:hi]
        fin = filter(isfinite, vec(win))
        vmin, vmax = isempty(fin) ? (0.0, 1.0) : (quantile(fin, 0.01), quantile(fin, 0.995))
        ext = (lo - 0.5, hi + 0.5, -0.5, N - 0.5)
        cols = Dict(1 => "#000000", 2 => "#ff7f00", 3 => "#1f78b4", 4 => "#e31a1c", 5 => "#9e9e9e", 6 => "#6a3d9a")
        names_ = Dict(1 => "pusty", 2 => "pojedynczy śmieć", 3 => "null / brak sygnału",
                      4 => "zaszumiony odcinek", 5 => "margines / domknięcie", 6 => "wygaszona emisja")

        rc("font", size=7.5, family="sans-serif")
        fig = figure(figsize=(9.2, 6.6), dpi=100)
        fig.patch.set_facecolor("white")

        ax1 = subplot2grid((1, 21), (0, 0), colspan=10)
        imshow(win, origin="lower", aspect="auto", cmap="afmhot", vmin=vmin, vmax=vmax,
               interpolation="nearest", extent=ext)
        title("PRZED (oryginał)", fontsize=9, fontweight="bold")
        xlabel("bin"); ylabel("pulse number (od 0)")

        ax2 = subplot2grid((1, 21), (0, 10), colspan=10, sharey=ax1)
        imshow(win, origin="lower", aspect="auto", cmap="afmhot", vmin=vmin, vmax=vmax,
               interpolation="nearest", extent=ext)
        # contiguous runs of the same reason -> one band each
        i = 1
        while i <= N
            r = a.reason[i]
            if r == 0
                i += 1
                continue
            end
            j = i
            while j < N && a.reason[j+1] == r; j += 1; end
            axhspan(i - 1.5, j - 0.5, color=cols[r], alpha=0.55, lw=0)
            i = j + 1
        end
        title("AUTOMAT: co wycina", fontsize=9, fontweight="bold")
        xlabel("bin")
        tick_params(labelleft=false)

        ax3 = subplot2grid((1, 21), (0, 20), sharey=ax1)
        strip = fill(NaN, N, 1)
        for i in 1:N; a.reason[i] > 0 && (strip[i, 1] = a.reason[i]); end
        cmap = pyimport("matplotlib.colors").ListedColormap([cols[k] for k in 1:6])
        imshow(strip, origin="lower", aspect="auto", cmap=cmap, vmin=0.5, vmax=6.5,
               interpolation="nearest", extent=(0, 1, -0.5, N - 0.5))
        axis("off")

        present = sort(unique(filter(>(0), a.reason)))
        handles = [patch.Patch(color=cols[k],
                   label="$(names_[k]) ($(count(==(k), a.reason)))") for k in present]
        isempty(handles) || fig.legend(handles=handles, loc="lower center", ncol=min(6, length(handles)),
                                       fontsize=7.5, frameon=false)
        n = length(a.bad)
        # given S/N pair (standard profile S/N, `Data.profile_snr`), else the automat's own measure
        sb, sa = isnothing(snr_pazi) ? (a.snr_before, a.snr_after) : snr_pazi
        lab = isnothing(snr_pazi) ? "S/N (miara automatu)" : "S/N profilu"
        suptitle(@sprintf("%s\nwycięte %d z %d (%.1f%%)   %s %.1f → %.1f", title_, n, N,
                          100 * n / N, lab, sb, sa), fontsize=9, fontweight="bold")
        tight_layout(rect=(0, isempty(handles) ? 0 : 0.05, 1, 0.92))
        png = "$outdir/$(name_mod)_autozap.pdf"   # PDF only (xpdf on kopernik)
        savefig(png, dpi=110)
        println(png)
        close()
        # external viewer instead of plt.show(): plt windows don't open here
        return show_ ? show_png(png; wait_=wait_) : nothing
    end


    """
    Final comparison of two zapping methods (`compare_zap_methods`): rows of
    `T` have psr, N, n_mine, n_maciej, none, mine, maciej (standard profile
    S/N on the same data and window).

    Left: S/N gain (vs no zaps) of my method against Maciej's — above the
    diagonal = mine better. Middle: per-pulsar difference of the gains,
    sorted. Right: fraction of pulses cut by each method.
    """
    function zap_methods_comparison(T, outdir; show_=false, tie=1.0)
        gm = [100 * (t.mine - t.none) / abs(t.none) for t in T]
        gj = [100 * (t.maciej - t.none) / abs(t.none) for t in T]
        fm = [100 * t.n_mine / t.N for t in T]
        fj = [100 * t.n_maciej / t.N for t in T]
        d = gm .- gj
        col = [x >= tie ? "tab:blue" : x <= -tie ? "tab:red" : "0.5" for x in d]

        fig, axs = subplots(1, 3, figsize=(11, 4.2))
        ax = axs[1]
        sc = [t.none >= 3 for t in T]    # S/N bez wyciec ~0 (sam RFI) daje procenty bez sensu - poza wykresem
        lo, hi = min(minimum(gm[sc]), minimum(gj[sc]), 0.0) - 2, max(maximum(gm[sc]), maximum(gj[sc]), 0.0) + 2
        ax.plot([lo, hi], [lo, hi], color="k", lw=0.8, ls="--")
        ax.scatter(gj[sc], gm[sc], c=col[sc], s=18, alpha=0.8)
        ax.axhline(0, color="0.8", lw=0.6); ax.axvline(0, color="0.8", lw=0.6)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_xlabel("zysk S/N — metoda ręczna [%]"); ax.set_ylabel("zysk S/N — automat [%]")
        ax.set_title("nad przekątną = automat lepszy", fontsize=9)

        ax = axs[2]
        o = sortperm(d)
        o = o[sc[o]]
        ax.bar(1:length(o), d[o], color=col[o], width=1.0)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlabel("pulsary (posortowane)"); ax.set_ylabel("zysk automat − ręcznie [punkty %]")
        nb, nw = count(d .>= tie), count(d .<= -tie)
        ax.set_title("automat lepszy: $nb   ręcznie lepiej: $nw   remis: $(length(d) - nb - nw)", fontsize=9)
        if length(T) <= 40
            ax.set_xticks(1:length(o)); ax.set_xticklabels([T[i].psr for i in o], rotation=90, fontsize=6)
        end

        ax = axs[3]
        m = max(maximum(fm), maximum(fj), 1.0) * 1.05
        ax.plot([0, m], [0, m], color="k", lw=0.8, ls="--")
        ax.scatter(fj, fm, c=col, s=18, alpha=0.8)
        ax.set_xlim(0, m); ax.set_ylim(0, m)
        ax.set_xlabel("wycięte — metoda ręczna [% impulsów]"); ax.set_ylabel("wycięte — automat [% impulsów]")
        ax.set_title("ile impulsów wycina każda metoda", fontsize=9)

        fig.suptitle("Wycinanie: automat vs metoda ręczna — standardowe S/N profilu " *
                     "($(length(T)) pulsarów; niebieskie = automat lepszy, czerwone = ręcznie lepiej, szare = remis <$(tie)%)",
                     fontsize=10, fontweight="bold")
        tight_layout(rect=(0, 0, 1, 0.94))
        png = joinpath(outdir, "zap_comparison.pdf")
        savefig(png, dpi=120)
        close()
        show_ && show_png(png)
        return png
    end


    """
    Hand zaps vs automatic zaps for one pulsar (`tune_autozap`).

    Left: single pulses around the on-pulse window; red bars on the left
    edge = your zaps (pazi), blue bars on the right edge = automatic.
    Right: smoothed signal-per-noise efficiency of each stretch (what the
    automatic decision uses) with the cut level frac·E dashed; your zapped
    pulses shaded red, automatic ones blue.
    """
    function zap_compare(data, bin_st, bin_end, a, user_rows, outdir; title_="", frac=1.0,
                         name_mod="pulsar", show_=false, wait_=true, labels=("Ty (pazi)", "automat"))
        N, nb = size(data)
        w = bin_end - bin_st + 1
        lo, hi = max(1, bin_st - w ÷ 2), min(nb, bin_end + w ÷ 2)
        win = data[:, lo:hi]
        fin = filter(isfinite, vec(win))
        vmin, vmax = isempty(fin) ? (0.0, 1.0) : (quantile(fin, 0.01), quantile(fin, 0.995))

        rc("font", size=7.5, family="sans-serif")
        fig = figure(figsize=(9.2, 6.6), dpi=100)
        fig.patch.set_facecolor("white")

        ax1 = subplot2grid((1, 3), (0, 0), colspan=2)
        imshow(win, origin="lower", aspect="auto", cmap="afmhot", vmin=vmin, vmax=vmax,
               interpolation="nearest", extent=(lo - 0.5, hi + 0.5, -0.5, N - 0.5))
        dx = (hi - lo) * 0.012
        isempty(user_rows) || scatter(fill(lo + dx, length(user_rows)), user_rows .- 1, marker="s",
                                      s=6, color="#e41a1c", zorder=5, label="$(labels[1]) ($(length(user_rows)))")
        isempty(a.bad) || scatter(fill(hi - dx, length(a.bad)), a.bad .- 1, marker="s",
                                  s=6, color="#377eb8", zorder=5, label="$(labels[2]) ($(length(a.bad)))")
        legend(fontsize=7, loc="upper center", ncol=2, framealpha=0.9)
        xlabel("bin"); ylabel("pulse number (od 0)")

        subplot2grid((1, 3), (0, 2), sharey=ax1)
        y = 0:N-1
        u = falses(N); u[filter(i -> 1 <= i <= N, user_rows)] .= true
        au = falses(N); au[a.bad] .= true
        for (mask, col) in ((u, "#e41a1c"), (au, "#377eb8"))
            i = 1
            while i <= N
                if mask[i]
                    j = i
                    while j < N && mask[j+1]; j += 1; end
                    axhspan(i - 1.5, j - 0.5, color=col, alpha=0.18, lw=0)
                    i = j + 1
                else
                    i += 1
                end
            end
        end
        plot(a.eff, y, color="black", lw=0.7)
        isfinite(a.E) && axvline(frac * a.E, color="gray", ls="--", lw=1.0)
        axvline(0, color="gray", ls=":", lw=0.6)
        xlabel("sygnał / szum odcinka\n(przerywana = próg cięcia)")
        tick_params(labelleft=false)

        suptitle(title_, fontsize=9, fontweight="bold")
        tight_layout(rect=(0, 0, 1, 0.94))
        png = "$outdir/$(name_mod)_zap_compare.pdf"
        savefig(png, dpi=110)
        println(png)
        close()
        # external viewer instead of plt.show(): plt windows don't open here
        return show_ ? show_png(png; wait_=wait_) : nothing
    end



    """
    pazi-like before/after view of the zapping: per column (before /
    automat / Maciej if `mrows` given) all pulses in place with the zapped
    ones left EMPTY (white), colour scale from the 1-99.5% percentiles of the
    remaining pulses only (so removed dark nulls / bright junk change the
    contrast), and the mean profile in S/N units (σ of the
    off-pulse) with the "before" profile dashed. Saved as
    `<name_mod>_contrast.pdf`.
    """
    function zap_contrast(data, bin_st, bin_end, rows, mrows, outdir; title_="", name_mod="pulsar",
                          off=nothing)
        N, nb = size(data)
        w = bin_end - bin_st + 1
        lo, hi = max(1, bin_st - w ÷ 2), min(nb, bin_end + w ÷ 2)
        isnothing(off) && (off = vcat(1:bin_st-1, bin_end+1:nb))
        empty = [all(iszero, @view data[i, :]) for i in 1:N]
        cols = Any[("PRZED", Int[]), ("AUTOMAT", rows)]
        isnothing(mrows) || push!(cols, ("MACIEJ", mrows))
        nc = length(cols)
        keeps, profs, snrs = BitVector[], Vector{Float64}[], Float64[]
        for (_, r) in cols
            keep = .!empty; keep[filter(i -> 1 <= i <= N, r)] .= false
            P = count(keep) > 0 ? vec(sum(@view(data[keep, :]), dims=1)) : zeros(nb)
            μ, σ = mean(P[off]), std(P[off])
            Pn = σ > 0 ? (P .- μ) ./ σ : P .* 0
            push!(keeps, keep); push!(profs, Pn); push!(snrs, sum(Pn[bin_st:bin_end]) / sqrt(w))
        end
        pmax = maximum(maximum(p[lo:hi]) for p in profs)
        pmin = minimum(minimum(p[lo:hi]) for p in profs)

        rc("font", size=7.5, family="sans-serif")
        fig = figure(figsize=(9.2, 6.6), dpi=100)
        fig.patch.set_facecolor("white")
        cmap = pyimport("copy").copy(PyPlot.get_cmap("afmhot"))
        cmap.set_bad("white")
        ax0 = nothing
        for (k, (name, _)) in enumerate(cols)
            win = data[:, lo:hi]
            win[.!keeps[k], :] .= NaN
            n = count(keeps[k])
            fin = vec(data[keeps[k], lo:hi])
            vmin, vmax = isempty(fin) ? (0.0, 1.0) : (quantile(fin, 0.01), quantile(fin, 0.995))
            ax = subplot2grid((4, nc), (0, k - 1), rowspan=3, sharey=ax0)
            isnothing(ax0) && (ax0 = ax)
            imshow(win, origin="lower", aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation="nearest", extent=(lo - 0.5, hi + 0.5, -0.5, N - 0.5))
            title("$name
zostało $n/$(count(.!empty))   S/N $(round(snrs[k], digits=1))",
                  fontsize=8.5, fontweight="bold")
            k == 1 ? ylabel("pulse number (od 0)") : tick_params(labelleft=false)
            tick_params(labelbottom=false)
            subplot2grid((4, nc), (3, k - 1))
            k > 1 && plot(lo:hi, profs[1][lo:hi], color="gray", lw=0.7, ls="--", label="przed")
            plot(lo:hi, profs[k][lo:hi], color="black", lw=0.9, label=name)
            xlim(lo, hi); ylim(pmin - 0.05 * (pmax - pmin), pmax * 1.05)
            xlabel("bin")
            k == 1 ? ylabel("profil [σ off-pulse]") : (tick_params(labelleft=false); legend(fontsize=6.5, loc="upper right"))
        end
        suptitle("$title_   — dane PO wycięciu (białe = wycięte, skala kolorów 1–99.5% pozostałych impulsów)",
                 fontsize=8.5, fontweight="bold")
        tight_layout(rect=(0, 0, 1, 0.95))
        pdf = "$outdir/$(name_mod)_contrast.pdf"
        savefig(pdf, dpi=110)
        println(pdf)
        close()
        return pdf
    end

    """
    Review window for bad-pulse zapping (`Data.review_zaps`).

    Left column = data BEFORE zapping, right = AFTER (current selection):
      top     single pulses around the on-pulse window (same colour scale on
              both sides, so RFI stands out on the left); zapped pulses are
              marked red on the left edge of the "before" panel
      middle  average profile, with its S/N in the title
      bottom  the two detection statistics (robust z of off-pulse noise and
              of the largest off-pulse spike), threshold dotted, selected
              pulses red. On-pulse emission is never used, so nulls stay.
    Pulse numbers are 0-based, the same numbering as `p["zaps"]`.

    With `show_=true` the window is left open (the caller asks the question
    and then calls `close_all`); otherwise it is only saved.
    """
    function zap_review(data, r, rows, bin_st, bin_end, outdir; name_mod="pulsar", show_=false)
        N, nb = size(data)
        w = bin_end - bin_st + 1
        lo, hi = max(1, bin_st - w ÷ 2), min(nb, bin_end + w ÷ 2)
        good_before = findall(i -> any(!iszero, @view data[i, :]), 1:N)
        good_after = setdiff(good_before, rows)
        after = copy(data)
        after[rows, :] .= 0

        win = data[good_before, lo:hi]
        vmin, vmax = isempty(win) ? (0.0, 1.0) : (quantile(vec(win), 0.01), quantile(vec(win), 0.995))
        snr_b = profile_snr(data, good_before, bin_st, bin_end)
        snr_a = profile_snr(data, good_after, bin_st, bin_end)
        # S/N of the sum of pulses from per-pulse signal and noise (what the
        # zapping decision optimises; the profile S/N above is often limited
        # by baseline systematics rather than noise)
        rad(rws) = isempty(rws) ? NaN : sum(r.sig[rws]) / sqrt(sum(r.noise_var[rws]))
        rad_b, rad_a = rad(good_before), rad(good_after)

        rc("font", size=7.5, family="sans-serif")
        rc("axes", linewidth=0.7, edgecolor="#444444")
        # fits a ~950 px wide VNC screen including the toolbar (zoom/pan)
        fig = figure(figsize=(9.2, 6.6), dpi=100)
        fig.patch.set_facecolor("white")
        ext = (lo - 0.5, hi + 0.5, -0.5, N - 0.5)

        for (col, (d, lab)) in enumerate(((data, "PRZED zapowaniem"), (after, "PO zapowaniu")))
            subplot2grid((7, 2), (0, col - 1), rowspan=3)
            imshow(d[:, lo:hi], origin="lower", aspect="auto", cmap="viridis",
                   vmin=vmin, vmax=vmax, interpolation="nearest", extent=ext)
            if col == 1 && !isempty(rows)
                scatter(fill(lo + 0.5, length(rows)), rows .- 1, marker="_", s=40,
                        color="#d62728", zorder=5)
            end
            title(lab, fontsize=10, fontweight="bold")
            ylabel("pulse number (od 0)")
            col == 2 && tick_params(labelleft=false)

            subplot2grid((7, 2), (3, col - 1))
            prof_rows = col == 1 ? good_before : good_after
            isempty(prof_rows) || plot(lo:hi, vec(mean(data[prof_rows, lo:hi], dims=1)),
                                       color="black", lw=0.9)
            axvspan(bin_st, bin_end, color="#3182bd", alpha=0.08)
            title(@sprintf("profil średni, S/N = %.1f", col == 1 ? snr_b : snr_a), fontsize=8.5)
            xlabel("bin", fontsize=8)
        end

        x = 0:N-1
        stats = ((r.z_rms, "skok szumu", true),
                 (r.z_struct, "linia bazowa", false),
                 (r.z_peak, "pojedynczy skok", false))
        for (k, (y, lab, twosided)) in enumerate(stats)
            subplot2grid((7, 2), (4 + k - 1, 0), colspan=2)
            grid(true, alpha=0.2, lw=0.4)
            plot(x, y, color="#3182bd", lw=0.6, zorder=2)
            isempty(r.kept_bright) || scatter(r.kept_bright .- 1, y[r.kept_bright], s=10,
                                              color="#ff9f1c", zorder=3)
            isempty(rows) || scatter(rows .- 1, y[rows], s=10, color="#d62728", zorder=4)
            axhline(r.nsigma, color="gray", ls=":", lw=1.0)
            twosided && axhline(-r.nsigma, color="gray", ls=":", lw=1.0)
            xlim(-0.5, N - 0.5)
            ylabel(lab * " (z)", fontsize=7.5)
            k < length(stats) ? tick_params(labelbottom=false) : xlabel("pulse number (od 0, jak w zaps)")
        end

        # zapping must not lower the S/N — if it does, good pulses are being cut
        warn = (snr_a < snr_b || rad_a < rad_b) ? "\n⚠ S/N SPADŁO — odrzucasz dobre impulsy?" : ""
        suptitle(@sprintf("Zapowanie: %d z %d (%.1f%%)   S/N sumy impulsów %.1f → %.1f   S/N profilu %.1f → %.1f",
                          length(rows), N, 100 * length(rows) / N, rad_b, rad_a, snr_b, snr_a) *
                 "\nczerwone = wycinane, pomarańczowe = podejrzane, ale zostawione (dużo sygnału)" * warn,
                 fontsize=9, fontweight="bold", color=(isempty(warn) ? "black" : "#b30000"))
        tight_layout(rect=(0, 0, 1, isempty(warn) ? 0.93 : 0.91))
        savefig("$outdir/$(name_mod)_bad_pulses.png", dpi=120)
        println("$outdir/$(name_mod)_bad_pulses.png")
        if show_ == true
            show()
        else
            close()
        end
    end


   function twodfs(data, outdir, params; cmap="viridis", darkness=0.3, name_mod="PSR_NAME", show_=false, average=nothing)

        p = params

        num, bins = size(data)

        st = p["bin_st"]
        en = p["bin_end"]
        nbin = p["nbin"]
        #signal_width = bin * (en - st + 1) / bins # en - st + 1 # TODO make proper range 
        #signal_width = (en - st + 1) / bins # en - st + 1 # TODO make proper range 
        signal_width = (en - st + 1) 

        # main panel data
        da = data 

        # side panels data
        sum_left = Functions.normalize(sum(da, dims=2)[:, 1])
        yleft = range(0, 0.5, length=num) # check this # should be fine
        sum_bottom = Functions.normalize(sum(da, dims=1)[1, :])
        xbottom = range(-nbin / 2, nbin/2, length=bins) 

        # plot style
        rc("font", size=7.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 4.33071))  # 8cm x 11cm
        subplots_adjust(left=0.17, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)

        # left panel
        ax_left = subplot2grid((4,3), (0,0), rowspan=3)
        minorticks_on()
        plot(sum_left, yleft, color="grey")
        ylim(yleft[1], yleft[end])
        xticks([0.5, 1.0])
        xlim(1.1, -0.1)
        ylabel(raw"Fluctuation frequency $(P/P_3)$")

        # Main panel - 2DFS
        ax_main = subplot2grid((4,3), (0,1), rowspan=3, colspan=2)
        im = imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto", extent=[xbottom[1], xbottom[end], yleft[1], yleft[end]], vmax=darkness*maximum(da))
        tick_params(left=false, labelleft=false)
        xlim(-signal_width,signal_width)
        #xlim(st-512, en - 512)
        #xlabel(raw"Fluctuation frequency $(P/P_2)$")

        # Bottom panel
        ax_bottom = subplot2grid((4,3), (3,1), colspan=2)
        minorticks_on()
        plot(xbottom, sum_bottom, color="grey")
        plot(-xbottom, sum_bottom, color="orange", ls="--")
        yticks([0., 0.5])
        xlim(-signal_width,signal_width)
        #xlim(st-512, en - 512)
        if !isnothing(average)
            ave = average[st:en]
            xave = collect(range(-signal_width/2, signal_width/2, length=length(ave))) # this is probably wrong..
            plot(xave, ave, color="red", lw=0.3, ls=":")
        end
        xlabel(raw"Fluctuation frequency $(P/P_2)$")

        savepath = "$outdir/$(name_mod)_2dfs.pdf"
        println(savepath)
        savefig(savepath)
        savefig(replace(savepath, "pdf"=>"png"))

        if show_ == true
            show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    function lrfs_plot_czarek_to_be_removed_keep_it_unchanged(data, outdir; start=1, number=nothing, bin_st=nothing, bin_end=nothing,
                       cmap="viridis", darkness=0.5, name_mod="0", show_=false)

        num, cols = size(data)
        if number === nothing
            number = num - start + 1
        end

        if bin_st === nothing
            bin_st = 1
        end
        if bin_end === nothing
            bin_end = (cols - 3) ÷ 2  # zakładam: 3 kolumny nagłówkowe + (Re,Im)*bin_count
        end

        # Rekonstrukcja zespolonych danych
        bins = bin_end - bin_st + 1
        lrfs_complex = ComplexF64[
            data[i, 4 + 2*(b-1)] + im * data[i, 5 + 2*(b-1)] for i in 1:num, b in bin_st:bin_end
        ]

        # Intensywność i częstotliwości (zakładam, że oś 1 to puls, oś 2 to longitude)
        intensity = sum(abs.(lrfs_complex), dims=2)[:]
        freq = collect(0:length(intensity)-1) ./ length(intensity)

        # Znajdź peak i P3
        peak = argmax(intensity)
        println("\tpeak freq $(freq[peak])")
        println("\tpeak P3 $(1/freq[peak])")

        # Faza dla linii peak
        phase_ = rad2deg.(angle.(view(lrfs_complex, peak, :)))

        # Ustalanie długości podłużnej
        db = bins
        dl = 360. * db / (cols ÷ 2)  # zakładam, że połowa kolumn to biny
        longitude = collect(range(-dl/2, dl/2, length=db))

        # Unwrap fazy (korekta ciągłości)
        for i in 1:length(phase_)-1
            while abs(phase_[i+1] - phase_[i]) > 180
                if phase_[i+1] > phase_[i]
                    phase_[i+1] -= 360
                else
                    phase_[i+1] += 360
                end
            end
        end

        # Styl plotu
        rc("font", size=7)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(6, 6))  # 8cm x 11cm
        subplots_adjust(left=0.17, bottom=0.08, right=0.90, top=0.92, wspace=0.0, hspace=0.0)

        # Górny panel: faza FFT
        ax_phase = subplot2grid((5, 3), (0, 1), colspan=2)
        scatter(longitude, phase_, marker=".", c="grey", s=3)
        ax_phase.xaxis.set_label_position("top")
        ax_phase.xaxis.set_ticks_position("top")
        xlabel("Phase \$(^\\circ)\$")
        #ylabel("FFT phase \$(^\\circ)\$")

        # Pobierz aktualny zakres osi Y
        ylims = ax_phase.get_ylim()
        ymin, ymax = ylims

        # Ustal pozycje ticków równomiernie od ymin do ymax
        nticks = 9  # bo mamy 9 wartości od -180 do 180 co 45 stopni
        ytick_positions = range(ymin, ymax, length=nticks)

        # Ustaw ticki na tych pozycjach, ale etykiety -180, -135, ..., 180
        ytick_labels = string.(-180:45:180)

        ax_phase.set_yticks(collect(ytick_positions))
        ax_phase.set_yticklabels(ytick_labels)

        tick_params(labeltop=true, labelbottom=false, which="both", bottom=false, top=true)



        # --- LEWY PANEL: intensywność vs częstotliwość ---
        pulses = collect(start:start+number-1)
        average_y = sum(abs.(lrfs_complex[start:start+number-1, :]), dims=2)[:,1]

        ax_intensity = subplot2grid((5, 3), (1, 0), rowspan=3)
        plot(average_y, pulses, color="grey")
        ylim(pulses[1], pulses[end])
        xticks([])
        ylabel("Fluctuation frequency (P/P₃)")

        # Tick marks Y (0.0 do 0.5 co 0.1)
        ytick_values = 0.0:0.1:0.5
        ytick_positions = pulses[1] .+ ytick_values .* (pulses[end] - pulses[1]) / 0.5
        yticks(ytick_positions, string.(ytick_values))

        # --- GŁÓWNY PANEL: mapa LRFS ---
        ax_lrfs = subplot2grid((5, 3), (1, 1), rowspan=3, colspan=2)
        imshow(abs.(lrfs_complex), origin="lower", cmap=cmap, interpolation="none", aspect="auto",
            vmax=darkness * maximum(abs.(lrfs_complex)))
        tick_params(labelleft=false, labelbottom=false)

        # --- DOLNY PANEL: średni profil ---
        ax_profile = subplot2grid((5, 3), (4, 1), colspan=2)
        average_profile = mean(real.(lrfs_complex), dims=1)[:]
        plot(longitude, average_profile, c="grey")
        axvline(x=0.0, ls=":", c="black")
        yticks([])
        xlim(longitude[1], longitude[end])
        xlabel("longitude \$(^\\circ)\$")

        # Zapis
        savefig("$outdir/$(name_mod)_lrfs.pdf")
        savefig("$outdir/$(name_mod)_lrfs.svg")

        if show_
            show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end

        close()
    end


    ##########################################################################
    ##########################################################################
    ##########################################################################


    """
        lrfs(data, outdir, params; pulses=nothing, nsig=3.0, ...)

    LRFS plot from PSRSALSA `pspec -lrfs` output (rows = fluctuation frequency 0..0.5 cycles/P,
    pol 1/2 = real/imag). Phase panel (2026-11-24 fix): phase of the LRFS at the peak frequency,
    shown only where |L| > `nsig` σ of the off-pulse LRFS at that frequency, unwrapped along
    longitude, with errors σ/|L|. With `pulses` (the debased pulse stack, pulses × bins) the phase is
    computed coherently from ALL pulses at the peak frequency (PSRSALSA averages blocks of nfft
    pulses, so its complex phase can be smeared); sign convention as PhaseDrift: kernel e^(+2πi f n),
    rising phase = drift towards later longitudes.
    """
    function lrfs(data, outdir, params; cmap="viridis", darkness=0.5, name_mod="PSR_NAME",
              show_=false, pulses=nothing, nsig=3.0)

        # === Parameters and shape ===
        p = params
        bin_st = p["bin_st"]
        bin_end = p["bin_end"]
        num, bins = size(data)

        #@assert bin_st >= 1 && bin_end <= bins "bin_st and bin_end are out of bounds"

        signal_width = bin_end - bin_st + 1

        # === Extract complex LRFS ===
        # Take the same bin range for both real and imaginary parts to form a complex matrix (num x signal_width)
        real_part = data[:, bin_st:bin_end, 1]  # real part
        imag_part = data[:, bin_st:bin_end, 2]  # imaginary part
        lrfs_complex = ComplexF64.(real_part .+ imag_part .* im)  # complex LRFS matrix

        # === Amplitude (main panel) ===
        amp = abs.(lrfs_complex)

        # === Mean profile (bottom panel) ===
        profile = if isnothing(pulses)
            mean(real.(lrfs_complex), dims=1)[1, :]
        else
            live = _live_rows(pulses)
            vec(mean(Float64.(pulses[live, bin_st:bin_end]), dims=1))
        end

        # === Power spectrum (left panel) ===
        spectrum = sum(amp, dims=2)[:, 1]

        # === Frequency axis ===
        # rows run from 0 to 0.5 cycles/P (nfft/2+1 rows); the old (0:num-1)/num ran to ~1 and the
        # printed P3 came out ~2x too small
        freq = collect(0:(num-1)) .* (0.5 / max(num - 1, 1))
        peak = argmax(spectrum[2:end]) + 1  # skip DC component at freq=0
        # znane P3 (params): szczyt szukany przy 1/P3 (±10%) - powolne zmiany / nulle nie przejmuja szczytu
        p3_ = get(p, "p3", nothing)
        if p3_ isa Number && p3_ > 0
            w = findall(f -> abs(f - 1 / p3_) <= 0.1 / p3_, freq)
            isempty(w) || (peak = w[argmax(spectrum[w])])
        end
        println("\tPeak freq: $(freq[peak]), P3 = $(1/freq[peak])")

        # === Phase at peak frequency ===
        Lrow = if isnothing(pulses)
            ComplexF64.(data[peak, :, 1] .+ im .* data[peak, :, 2])
        else
            _coherent_lrfs_row(pulses, freq[peak])
        end
        ph = _lrfs_phase_track(Lrow, bin_st, bin_end; nsig=nsig)

        # === Longitude axis ===
        # The longitude axis has the same length as signal_width and is centered around zero
        dl = 360 * signal_width / bins
        longitude = range(-dl/2, dl/2, length=signal_width)
        println("Longitude resolution = $(longitude[2] - longitude[1]) deg")



        # === Plotting settings ===
        rc("font", size=7)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(6, 7))  # figure size in inches (~15.24 cm × 17.78 cm)
        subplots_adjust(left=0.17, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)

        # --- LEFT PANEL: intensity vs frequency ---
        ax_left = subplot2grid((5, 3), (1, 0), rowspan=3)
        plot(spectrum, freq, color="grey")
        ylim(freq[1], freq[end])
        xticks([])
        ylabel("Fluctuation frequency (P/P₃)")

        # --- Custom Y ticks (0.0 to 0.5 every 0.1) ---
        ytick_values = 0.0:0.1:0.5
        ytick_positions = freq[1] .+ ytick_values .* (freq[end] - freq[1]) / 0.5
        yticks(ytick_positions, string.(ytick_values))

        # --- MAIN PANEL: LRFS map ---
        ax_main = subplot2grid((5, 3), (1, 1), rowspan=3, colspan=2)
        imshow(amp, origin="lower", cmap=cmap, interpolation="none", aspect="auto",
            extent=[longitude[1], longitude[end], freq[1], freq[end]],
            vmax=darkness * maximum(amp))
        tick_params(labelleft=false, labelbottom=false)

        # --- TOP PANEL: phase (only significant bins, unwrapped, with errors) ---
        ax_phase = subplot2grid((5, 3), (0, 1), colspan=2)
        if any(ph.sig)
            s = ph.sig
            errorbar(longitude[s], ph.phase[s], yerr=ph.err[s], fmt=".", ms=2, c="black", elinewidth=0.4, capsize=0)
            lo, hi = extrema(ph.phase[s])
            ylim(lo - 30, hi + 30)
        else
            text(0.5, 0.5, "no significant phase (|L| < $(nsig)σ)", transform=ax_phase.transAxes, ha="center", va="center")
        end
        xlim(longitude[1], longitude[end])
        xticks([])
        ylabel("Phase (°)", labelpad=10)

        # --- BOTTOM PANEL: mean profile ---
        ax_bottom = subplot2grid((5, 3), (4, 1), colspan=2)
        plot(longitude, profile, color="grey")
        axvline(x=0.0, ls=":", c="black")
        yticks([])
        xlim(longitude[1], longitude[end])
        xlabel("Longitude (°)")

        # === Save figures ===
        println("$outdir/$(name_mod)_lrfs.pdf")
        savefig("$outdir/$(name_mod)_lrfs.pdf")
        savefig("$outdir/$(name_mod)_lrfs.svg")
        savefig("$outdir/$(name_mod)_lrfs.png")

        if show_
            show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end

    # TODO TODO TODO test and remove/clean below

    function analyse_single(low, high, mid, p; pulse_start=1)
        pulses, bins = size(low)
        for i in pulse_start:pulses

            figure(figsize=(6, 7))  # figure size in inches (~15.24 cm × 17.78 cm)
            plot(low[i, p["bin_st"]:p["bin_end"]], label="low", lw=2, alpha=0.77)
            plot(high[i, p["bin_st"]:p["bin_end"]], label="high", lw=2, alpha=0.77)
            #plot(mid[i, p["bin_st"]:p["bin_end"]], label="mid", lw=2, alpha=0.77)

            legend()

            show()
            println("Press Enter for next pulse, 'q' to quit.")
            user_input = readline(stdin; keep=false)
            close("all")
            
            if lowercase(strip(user_input)) == "q"
                println("Exiting analysis.")
                break
            end

        end

    end


    function analyse_p3folds(low, high, p)
        pulses, bins = size(low)
        for i in 1:pulses


            figure(figsize=(6, 7))  # figure size in inches (~15.24 cm × 17.78 cm)
            plot(low[i, p["bin_st"]:p["bin_end"]])
            plot(high[i, p["bin_st"]:p["bin_end"]])

            show()
            println("Press Enter for next pulse, 'q' to quit.")
            user_input = readline(stdin; keep=false)
            close("all")
            
            if lowercase(strip(user_input)) == "q"
                println("Exiting analysis.")
                break
            end

        end


    end

    function analyse_p3folds2(low, high, p)
        pulses, bins = size(low)
        
        # Define model: sum of two Gaussian functions
        @. model(x, p) = p[1] * exp(-((x - p[2])^2) / (2 * p[3]^2)) + 
                        p[4] * exp(-((x - p[5])^2) / (2 * p[6]^2))
        
        for i in 1:pulses
            # Extract data in the range bin_st:bin_end
            x_data = p["bin_st"]:p["bin_end"]
            y_low = low[i, p["bin_st"]:p["bin_end"]]
            y_high = high[i, p["bin_st"]:p["bin_end"]]
            
            figure(figsize=(6, 7))
            
            # Subplot for low
            subplot(2, 1, 1)
            plot(x_data, y_low, "b.", label="Low data")
            
            # Initial parameters for fitting [amp1, mean1, sigma1, amp2, mean2, sigma2]
            # Find two maxima as starting points
            max_idx = argmax(y_low)
            p0_low = [maximum(y_low), x_data[max_idx], 10.0, 
                    maximum(y_low)*0.5, x_data[max_idx]+20, 10.0]
            
            try
                fit_low = curve_fit(model, collect(x_data), y_low, p0_low)
                plot(x_data, model(collect(x_data), fit_low.param), "r-", label="Fit", linewidth=2)
                title("Low - Pulse $i")
                legend()
            catch e
                println("Fit failed for low data, pulse $i: $e")
                title("Low - Pulse $i (fit failed)")
            end
            
            # Subplot for high
            subplot(2, 1, 2)
            plot(x_data, y_high, "b.", label="High data")
            
            max_idx = argmax(y_high)
            p0_high = [maximum(y_high), x_data[max_idx], 10.0,
                    maximum(y_high)*0.5, x_data[max_idx]+20, 10.0]
            
            try
                fit_high = curve_fit(model, collect(x_data), y_high, p0_high)
                plot(x_data, model(collect(x_data), fit_high.param), "r-", label="Fit", linewidth=2)
                title("High - Pulse $i")
                legend()
                
                # Optional: print fitting parameters
                println("High fit params: A1=$(fit_high.param[1]), μ1=$(fit_high.param[2]), σ1=$(fit_high.param[3])")
                println("                A2=$(fit_high.param[4]), μ2=$(fit_high.param[5]), σ2=$(fit_high.param[6])")
            catch e
                println("Fit failed for high data, pulse $i: $e")
                title("High - Pulse $i (fit failed)")
            end
            
            tight_layout()
            show()
            
            println("Press Enter for next pulse, 'q' to quit.")
            user_input = readline(stdin; keep=false)
            close("all")
            
            if lowercase(strip(user_input)) == "q"
                println("Exiting analysis.")
                break
            end
        end
    end    

    """
    Weighted mean offset and mean longitude per component, plus the component
    separation that follows from them. Returns nothing when there is nothing
    to summarise.

    `offset_data` maps a component index to the per-profile longitudes (deg),
    offsets (deg) and offset errors (deg) collected by the analyse_* routines.
    The longitude GaussianFit.component_offsets returns is the mean over the
    two frequencies, (mu_high + mu_low) / 2, so its error is half the offset
    error — the same weights therefore serve both means and only sigma_int
    halves.

    The separation is taken between the outermost components, matching what
    _read_offsets builds out of input/offsets.csv, so that joining the two by
    pulsar name gives the fractional narrowing delta_sep / sep. Errors quoted
    for the separation are the scatter-inflated sigma_ext, not sigma_int: the
    profile-to-profile spread is what actually limits these numbers.

    With a non-empty `psr` the same row is written to `outfile` (one line per
    pulsar, an existing row for that pulsar is replaced).
    """
    function _offset_summary(offset_data; psr="",
                             outfile=normpath(joinpath(@__DIR__, "..", "input", "separations.csv")))

        comps = sort(collect(keys(offset_data)))
        filter!(c -> !isempty(offset_data[c].err) && !all(offset_data[c].err .== 0.0), comps)
        isempty(comps) && return nothing

        # per component: weighted means of the offset and of the longitude
        stat = Dict{Int,NamedTuple}()
        for c in comps
            d = offset_data[c]
            w = 1.0 ./ (d.err .^ 2)
            n = length(d.off)
            dof = n - 1
            off_mu = sum(w .* d.off) / sum(w)
            off_int = 1.0 / sqrt(sum(w))
            off_chi2 = dof > 0 ? sum(w .* (d.off .- off_mu) .^ 2) / dof : NaN
            lon_mu = sum(w .* d.lon) / sum(w)
            lon_int = off_int / 2                     # lon_err = off_err / 2
            lon_chi2 = dof > 0 ? sum(4 .* w .* (d.lon .- lon_mu) .^ 2) / dof : NaN
            stat[c] = (n=n,
                       off=off_mu, off_int=off_int, off_chi2=off_chi2,
                       off_ext=off_int * sqrt(max(1.0, off_chi2)),
                       lon=lon_mu, lon_int=lon_int, lon_chi2=lon_chi2,
                       lon_ext=lon_int * sqrt(max(1.0, lon_chi2)))
        end

        println("\n=== Weighted mean offsets ===")
        for c in comps
            s = stat[c]
            println(@sprintf("G%d: offset = %+.4f° ± %.4f°  (n=%d, χ²/dof = %.2f, σ_ext = %.4f°)",
                             c, s.off, s.off_int, s.n, s.off_chi2, s.off_ext))
        end

        println("\n=== Component longitudes and separation ===")
        for c in comps
            s = stat[c]
            println(@sprintf("G%d: longitude = %.4f° ± %.4f°  (χ²/dof = %.2f, σ_ext = %.4f°)",
                             c, s.lon, s.lon_int, s.lon_chi2, s.lon_ext))
        end

        sep, sep_err = NaN, NaN
        dsep, dsep_err = NaN, NaN
        frac, frac_err = NaN, NaN
        if length(comps) > 1
            a, b = stat[comps[1]], stat[comps[end]]
            sep       = b.lon - a.lon
            sep_err   = sqrt(b.lon_ext^2 + a.lon_ext^2)
            dsep      = b.off - a.off
            dsep_err  = sqrt(b.off_ext^2 + a.off_ext^2)
            # frac = dsep / sep, propagated so that dsep ~ 0 stays finite
            frac      = dsep / sep
            frac_err  = sqrt((dsep_err / sep)^2 + (dsep * sep_err / sep^2)^2)
            println(@sprintf("separation G%d-G%d = %.4f° ± %.4f°  (mid-band)",
                             comps[1], comps[end], sep, sep_err))
            println(@sprintf("Δ separation       = %+.4f° ± %.4f°", dsep, dsep_err))
            println(@sprintf("fractional ΔW/W    = %+.4f ± %.4f", frac, frac_err))
        else
            println("single component — no separation")
        end
        println()

        isempty(psr) && return (comps=comps, stat=stat, sep=sep, sep_err=sep_err,
                                dsep=dsep, dsep_err=dsep_err, frac=frac, frac_err=frac_err)

        # one row per pulsar, the longitudes in the same four slots offsets.csv
        # uses so that the two files join on the name
        num(x) = isnan(x) ? "" : @sprintf("%.4f", x)
        slots = String[]
        for c in 1:4
            i = findfirst(==(c), comps)
            push!(slots, isnothing(i) ? "" : num(stat[c].lon),
                         isnothing(i) ? "" : num(stat[c].lon_ext))
        end
        row = join([psr; string(length(comps)); slots; num(sep); num(sep_err)], ",")

        header = "Nazwa,ncomp,lon 1,lon 1 err,lon 2,lon 2 err,lon 3,lon 3 err,lon 4,lon 4 err,sep,sep err"
        kept = String[]
        if isfile(outfile)
            for (i, line) in enumerate(eachline(outfile))
                i == 1 && continue
                s = strip(line)
                (isempty(s) || strip(first(split(s, ','))) == psr) && continue
                push!(kept, s)
            end
        end
        open(outfile, "w") do f
            println(f, header)
            for l in sort(push!(kept, row), by = l -> String(strip(first(split(l, ',')))))
                println(f, l)
            end
        end
        println("$psr written to $outfile")

        return (comps=comps, stat=stat, sep=sep, sep_err=sep_err,
                dsep=dsep, dsep_err=dsep_err, frac=frac, frac_err=frac_err,
                row=row, header=header)
    end


    function analyse_p3folds3(low, high, p, n_comp)
        pulses, bins = size(low)

        bin_st  = something(get(p, "bin_st",  nothing), 1)
        bin_end = something(get(p, "bin_end", nothing), bins)

        # Collected offsets per pulse: Dict(component => (longitudes, offsets, errors))
        offset_data = Dict{Int, NamedTuple{(:lon, :off, :err), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}}()

        for i in 1:pulses
            # Extract data in the range bin_st:bin_end
            x_data = bin_st:bin_end
            y_low = low[i, bin_st:bin_end]
            y_high = high[i, bin_st:bin_end]

            #=
            # Tutaj fity nie działają...
            x_data = collect(1:bins)
            y_low = low[i, :]
            y_high = high[i, :]
            =#


            fit_l = GaussianFit.fit_gaussians(x_data, y_low, n_comp)
            #println(fit_l)
            GaussianFit.print_fit_summary(fit_l, n_comp; label="1023 MHz", nbin=1024)
            fit_h = GaussianFit.fit_gaussians(x_data, y_high, n_comp)
            GaussianFit.print_fit_summary(fit_h, n_comp; label="1523 MHz", nbin=1024)


            # Component offsets
            for o in GaussianFit.component_offsets(fit_h, fit_l; nbin=1024)
                @printf("G%d: lon=%.2f° bin=%.1f  %+.3f ± %.3f bins = %+.3f° ± %.3f°\n",
                    o.component, o.longitude, o.longitude_bin, o.offset_bins, o.offset_err, o.offset_deg, o.offset_deg_err)
                if !haskey(offset_data, o.component)
                    offset_data[o.component] = (lon=Float64[], off=Float64[], err=Float64[])
                end
                push!(offset_data[o.component].lon, o.longitude)
                push!(offset_data[o.component].off, o.offset_deg)
                push!(offset_data[o.component].err, o.offset_deg_err)
            end




            figure(figsize=(6, 7))

            # Subplot for low
            low_colors  = ["#2196F3", "#0D47A1", "#64B5F6", "#1565C0"]
            high_colors = ["#FF6F00", "#E65100", "#FFCA28", "#F57F17"]

            subplot(2, 1, 1)
            plot(x_data, y_low,  color="steelblue",  lw=1.2, alpha=0.7, label="Low data")
            plot(x_data, y_high, color="darkorange",  lw=1.2, alpha=0.7, label="High data")
            if fit_l.converged
                plot(x_data, fit_l.yfit, "--", color="blue", lw=2.0, label="Fit low")
                for (j, c) in enumerate(fit_l.components)
                    plot(x_data, fit_l.baseline .+ GaussianFit._gauss.(x_data, c.A, c.mu, c.sigma),
                         color=low_colors[mod1(j, length(low_colors))], lw=1.8, alpha=0.9, label="Low G$j")
                end
            end
            if fit_h.converged
                plot(x_data, fit_h.yfit, "--", color="darkorange", lw=2.0, label="Fit high")
                for (j, c) in enumerate(fit_h.components)
                    plot(x_data, fit_h.baseline .+ GaussianFit._gauss.(x_data, c.A, c.mu, c.sigma),
                         color=high_colors[mod1(j, length(high_colors))], lw=1.8, alpha=0.9, label="High G$j")
                end
            end
            legend()
            xlim(bin_st, bin_end)

            # Bottom panel: mu vs amplitude for each component
            subplot(2, 1, 2)
            if fit_l.converged
                for (j, c) in enumerate(fit_l.components)
                    scatter([c.mu], [fit_l.baseline + c.A], color=low_colors[mod1(j, length(low_colors))],
                            marker="o", s=60, label="Low G$j", zorder=3)
                end
            end
            if fit_h.converged
                for (j, c) in enumerate(fit_h.components)
                    scatter([c.mu], [fit_h.baseline + c.A], color=high_colors[mod1(j, length(high_colors))],
                            marker="s", s=60, label="High G$j", zorder=3)
                end
            end
            xlabel("μ (bin)")
            ylabel("Amplituda")
            title("Komponenty: μ vs amplituda")
            legend()
            xlim(bin_st, bin_end)

            tight_layout()
            show()
            
            println("Press Enter for next pulse, 'q' to quit.")
            user_input = readline(stdin; keep=false)
            close("all")
            
            if lowercase(strip(user_input)) == "q"
                println("Exiting analysis.")
                break
            end




        end

        # Weighted mean offset per component
        if !isempty(offset_data)
            println("\n=== Weighted mean offsets ===")
            for comp in sort(collect(keys(offset_data)))
                d = offset_data[comp]
                if isempty(d.err) || all(d.err .== 0.0)
                    continue
                end
                w      = 1.0 ./ (d.err .^ 2)
                n      = length(d.off)
                mu     = sum(w .* d.off) / sum(w)
                sigma_int = 1.0 / sqrt(sum(w))
                chi2   = sum(w .* (d.off .- mu) .^ 2)
                dof    = n - 1
                chi2_red = dof > 0 ? chi2 / dof : NaN
                sigma_ext = sigma_int * sqrt(max(1.0, chi2_red))
                println(@sprintf("G%d: offset = %+.4f° ± %.4f°  (n=%d, χ²/dof = %.2f, σ_ext = %.4f°)",
                    comp, mu, sigma_int, n, chi2_red, sigma_ext))
            end
            println()
        end

        # Final plot: longitude vs. offset for each component
        if !isempty(offset_data)
            figure(figsize=(8, 5))
            colors = ["#2196F3", "#E65100", "#4CAF50", "#9C27B0"]
            for comp in sort(collect(keys(offset_data)))
                d = offset_data[comp]
                errorbar(d.lon, d.off, yerr=d.err,
                         fmt="o", capsize=4, lw=1.2,
                         color=colors[mod1(comp, length(colors))],
                         label="G$comp")
            end
            #ylim(-0.5, 1.5)
            axhline(0.0, color="gray", lw=0.8, ls="--")
            minorticks_on()
            xlabel("Longitude (°)")
            ylabel("Offset (°)")
            title("Longitude vs. offset")
            legend()
            tight_layout()
            show()
        end

    end


    """
    As analyse_p3folds3, with 's' added to skip a point. Passing `psr` makes
    the component longitudes and the separation land in `separations` as well
    (see _offset_summary); without it they are only printed.
    """
    function analyse_p3folds4(low, high, p, n_comp; psr="",
                              separations=normpath(joinpath(@__DIR__, "..", "input", "separations.csv")),
                              save_summary_dir=nothing)
        pulses, bins = size(low)

        bin_st  = something(get(p, "bin_st",  nothing), 1)
        bin_end = something(get(p, "bin_end", nothing), bins)

        # Collected offsets per pulse: Dict(component => (longitudes, offsets, errors))
        offset_data = Dict{Int, NamedTuple{(:lon, :off, :err), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}}()

        for i in 1:pulses
            # Extract data in the range bin_st:bin_end
            x_data = bin_st:bin_end
            y_low = low[i, bin_st:bin_end]
            y_high = high[i, bin_st:bin_end]

            #=
            # Tutaj fity nie działają...
            x_data = collect(1:bins)
            y_low = low[i, :]
            y_high = high[i, :]
            =#

            fit_l = GaussianFit.fit_gaussians(x_data, y_low, n_comp)
            GaussianFit.print_fit_summary(fit_l, n_comp; label="1023 MHz", nbin=1024)
            fit_h = GaussianFit.fit_gaussians(x_data, y_high, n_comp)
            GaussianFit.print_fit_summary(fit_h, n_comp; label="1523 MHz", nbin=1024)

            # Plot setup for visual evaluation
            figure(figsize=(6, 7))

            # Subplot for low
            low_colors  = ["#2196F3", "#0D47A1", "#64B5F6", "#1565C0"]
            high_colors = ["#FF6F00", "#E65100", "#FFCA28", "#F57F17"]

            subplot(2, 1, 1)
            plot(x_data, y_low,  color="steelblue",  lw=1.2, alpha=0.7, label="Low data")
            plot(x_data, y_high, color="darkorange",  lw=1.2, alpha=0.7, label="High data")
            if fit_l.converged
                plot(x_data, fit_l.yfit, "--", color="blue", lw=2.0, label="Fit low")
                for (j, c) in enumerate(fit_l.components)
                    plot(x_data, fit_l.baseline .+ GaussianFit._gauss.(x_data, c.A, c.mu, c.sigma),
                         color=low_colors[mod1(j, length(low_colors))], lw=1.8, alpha=0.9, label="Low G$j")
                end
            end
            if fit_h.converged
                plot(x_data, fit_h.yfit, "--", color="darkorange", lw=2.0, label="Fit high")
                for (j, c) in enumerate(fit_h.components)
                    plot(x_data, fit_h.baseline .+ GaussianFit._gauss.(x_data, c.A, c.mu, c.sigma),
                         color=high_colors[mod1(j, length(high_colors))], lw=1.8, alpha=0.9, label="High G$j")
                end
            end
            legend()
            xlim(bin_st, bin_end)

            # Bottom panel: mu vs amplitude for each component
            subplot(2, 1, 2)
            if fit_l.converged
                for (j, c) in enumerate(fit_l.components)
                    scatter([c.mu], [fit_l.baseline + c.A], color=low_colors[mod1(j, length(low_colors))],
                            marker="o", s=60, label="Low G$j", zorder=3)
                end
            end
            if fit_h.converged
                for (j, c) in enumerate(fit_h.components)
                    scatter([c.mu], [fit_h.baseline + c.A], color=high_colors[mod1(j, length(high_colors))],
                            marker="s", s=60, label="High G$j", zorder=3)
                end
            end
            xlabel("μ (bin)")
            ylabel("Amplituda")
            title("Komponenty: μ vs amplituda")
            legend()
            xlim(bin_st, bin_end)

            tight_layout()
            show()
            
            # User evaluation prompt
            println("analyse_p3folds4: Press Enter for next pulse, 's' to skip this point, 'q' to quit.")
            user_input = lowercase(strip(readline(stdin; keep=false)))
            close("all")
            
            if user_input == "q"
                println("Exiting analysis loop early.")
                break
            elseif user_input == "s"
                println("Skipped pulse $i from downstream analysis.")
                continue # Move immediately to loop step i+1 without processing offsets
            end

            # Component offsets (only processed/saved if not skipped)
            for o in GaussianFit.component_offsets(fit_h, fit_l; nbin=1024)
                @printf("G%d: lon=%.2f° bin=%.1f  %+.3f ± %.3f bins = %+.3f° ± %.3f°\n",
                    o.component, o.longitude, o.longitude_bin, o.offset_bins, o.offset_err, o.offset_deg, o.offset_deg_err)
                if !haskey(offset_data, o.component)
                    offset_data[o.component] = (lon=Float64[], off=Float64[], err=Float64[])
                end
                push!(offset_data[o.component].lon, o.longitude)
                push!(offset_data[o.component].off, o.offset_deg)
                push!(offset_data[o.component].err, o.offset_deg_err)
            end
        end

        # Weighted mean offset and longitude per component, plus the separation
        summary_result = _offset_summary(offset_data; psr=psr, outfile=separations)

        # Final plot: longitude vs. offset for each component
        if !isempty(offset_data)
            figure(figsize=(8, 5))
            colors = ["#2196F3", "#E65100", "#4CAF50", "#9C27B0"]
            for comp in sort(collect(keys(offset_data)))
                d = offset_data[comp]
                errorbar(d.lon, d.off, yerr=d.err,
                         fmt="o", capsize=4, lw=1.2,
                         color=colors[mod1(comp, length(colors))],
                         label="G$comp")
            end
            #ylim(-0.5, 1.5)
            axhline(0.0, color="gray", lw=0.8, ls="--")
            minorticks_on()
            xlabel("Longitude (°)")
            ylabel("Offset (°)")
            title("Longitude vs. offset" * (isempty(psr) ? "" : "  ($psr)"))
            legend()
            tight_layout()
            if !isnothing(save_summary_dir)
                mkpath(save_summary_dir)
                savefig(joinpath(save_summary_dir, "$(psr)_offset_summary.png"), dpi=120)
            end
            show()
        end

        return summary_result
    end







    """
    For each averaged profile (row in nl/nh) fit n_comp Gaussians to the
    low- and high-frequency profiles and collect the frequency-dependent
    offset per component.  Interactive: press Enter to advance, 'q' to quit.

    # Arguments
    - `nl`:         matrix (n_profiles × bins) of averaged low-freq profiles
    - `nh`:         matrix (n_profiles × bins) of averaged high-freq profiles
    - `p`:          params dict (must contain bin_st, bin_end)
    - `n_comp`:     number of Gaussian components
    - `npulse`:    single pulses per averaged profile (used for x-axis labels)
    - `n_pulses`:  total number of single pulses (used for last-block label)
    """
    function analyse_average_offset(nl, nh, p, n_comp; npulse=150, n_pulses=nothing, psr="",
                                    separations=normpath(joinpath(@__DIR__, "..", "input", "separations.csv")))
        n_profiles = size(nl, 1)
        bins = size(nl, 2)

        bin_st  = something(get(p, "bin_st",  nothing), 1)
        bin_end = something(get(p, "bin_end", nothing), bins)

        offset_data = Dict{Int, NamedTuple{(:idx, :lon, :off, :err),
                                           Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}}}()

        low_colors  = ["#2196F3", "#0D47A1", "#64B5F6", "#1565C0"]
        high_colors = ["#FF6F00", "#E65100", "#FFCA28", "#F57F17"]

        for i in 1:n_profiles
            x_data = bin_st:bin_end
            y_low  = nl[i, bin_st:bin_end]
            y_high = nh[i, bin_st:bin_end]

            fit_l = GaussianFit.fit_gaussians(x_data, y_low,  n_comp)
            fit_h = GaussianFit.fit_gaussians(x_data, y_high, n_comp)

            if fit_l.converged && fit_h.converged
                GaussianFit.print_fit_summary(fit_l, n_comp; label="low  (profile $i)", nbin=1024)
                GaussianFit.print_fit_summary(fit_h, n_comp; label="high (profile $i)", nbin=1024)

                pulse_center = (i - 0.5) * npulse

                for o in GaussianFit.component_offsets(fit_h, fit_l; nbin=1024)
                    @printf("G%d: lon=%.2f° bin=%.1f  %+.3f ± %.3f bins = %+.3f° ± %.3f°\n",
                        o.component, o.longitude, o.longitude_bin,
                        o.offset_bins, o.offset_err, o.offset_deg, o.offset_deg_err)
                    if !haskey(offset_data, o.component)
                        offset_data[o.component] = (idx=Float64[], lon=Float64[],
                                                    off=Float64[], err=Float64[])
                    end
                    push!(offset_data[o.component].idx, pulse_center)
                    push!(offset_data[o.component].lon, o.longitude)
                    push!(offset_data[o.component].off, o.offset_deg)
                    push!(offset_data[o.component].err, o.offset_deg_err)
                end
            else
                println("Profile $i: fit did not converge (low=$(fit_l.converged), high=$(fit_h.converged))")
            end

            figure(figsize=(6, 7))

            subplot(2, 1, 1)
            plot(x_data, y_low,  color="steelblue",  lw=1.2, alpha=0.7, label="Low")
            plot(x_data, y_high, color="darkorange",  lw=1.2, alpha=0.7, label="High")
            if fit_l.converged
                plot(x_data, fit_l.yfit, "--", color="blue", lw=2.0, label="Fit low")
                for (j, c) in enumerate(fit_l.components)
                    plot(x_data, fit_l.baseline .+ GaussianFit._gauss.(x_data, c.A, c.mu, c.sigma),
                         color=low_colors[mod1(j, length(low_colors))], lw=1.8, alpha=0.9, label="Low G$j")
                end
            end
            if fit_h.converged
                plot(x_data, fit_h.yfit, "--", color="darkorange", lw=2.0, label="Fit high")
                for (j, c) in enumerate(fit_h.components)
                    plot(x_data, fit_h.baseline .+ GaussianFit._gauss.(x_data, c.A, c.mu, c.sigma),
                         color=high_colors[mod1(j, length(high_colors))], lw=1.8, alpha=0.9, label="High G$j")
                end
            end
            st_pulse = (i - 1) * npulse + 1
            en_pulse = isnothing(n_pulses) ? i * npulse : min(i * npulse, n_pulses)
            title("Profile $i  (pulses $(st_pulse)-$(en_pulse))")
            legend(fontsize=7)
            xlim(bin_st, bin_end)

            subplot(2, 1, 2)
            if fit_l.converged
                for (j, c) in enumerate(fit_l.components)
                    scatter([c.mu], [fit_l.baseline + c.A],
                            color=low_colors[mod1(j, length(low_colors))],
                            marker="o", s=60, label="Low G$j", zorder=3)
                end
            end
            if fit_h.converged
                for (j, c) in enumerate(fit_h.components)
                    scatter([c.mu], [fit_h.baseline + c.A],
                            color=high_colors[mod1(j, length(high_colors))],
                            marker="s", s=60, label="High G$j", zorder=3)
                end
            end
            xlabel("μ (bin)")
            ylabel("Amplitude")
            title("Components: μ vs amplitude")
            legend(fontsize=7)
            xlim(bin_st, bin_end)

            tight_layout()
            show()

            println("Press Enter for next profile ($i/$n_profiles), 'q' to quit.")
            user_input = readline(stdin; keep=false)
            close("all")

            if lowercase(strip(user_input)) == "q"
                println("Exiting analysis.")
                break
            end
        end

        # Weighted mean offset and longitude per component, plus the separation
        _offset_summary(offset_data; psr=psr, outfile=separations)

        # Final summary: pulse number vs offset per component
        if !isempty(offset_data)
            figure(figsize=(8, 5))
            colors = ["#2196F3", "#E65100", "#4CAF50", "#9C27B0"]
            for comp in sort(collect(keys(offset_data)))
                d = offset_data[comp]
                errorbar(d.idx, d.off, yerr=d.err,
                         fmt="o", capsize=4, lw=1.2,
                         color=colors[mod1(comp, length(colors))],
                         label="G$comp")
            end
            axhline(0.0, color="gray", lw=0.8, ls="--")
            minorticks_on()
            xlabel("Pulse number (center of chunk)")
            ylabel("Offset (°)")
            title("Offset vs pulse number")
            legend()
            tight_layout()
            show()
        end

    end


    function position_angle(lon_l, pa_l, pa_err_l, I_l, Lin_l, V_l,
                            lon_h, pa_h, pa_err_h, I_h, Lin_h, V_h,
                            outdir; show_=true,
                            lon_rvm_l=nothing, pa_rvm_l=nothing, pa_rvm_l_ortho=nothing,
                            phi0_l=nothing,
                            lon_rvm_h=nothing, pa_rvm_h=nothing, pa_rvm_h_ortho=nothing,
                            phi0_h=nothing)
        figure(figsize=(5.354337, 6.8))
        subplots_adjust(left=0.18, bottom=0.10, right=0.99, top=0.99, hspace=0.05)

        # Top: PA with error bars
        subplot2grid((3, 1), (0, 0))
        errorbar(lon_l, pa_l, yerr=pa_err_l, fmt=".", ms=2, elinewidth=0.6,
                 c="tab:blue",   label="low",  zorder=3)
        errorbar(lon_h, pa_h, yerr=pa_err_h, fmt=".", ms=2, elinewidth=0.6,
                 c="tab:orange", label="high", zorder=3)

        # Best-fitting RVM (orange) and its 90° orthogonal mode
        if !isnothing(lon_rvm_l) && !isnothing(pa_rvm_l)
            plot(lon_rvm_l, pa_rvm_l,       c="darkorange", lw=1.2, zorder=2)
            plot(lon_rvm_l, pa_rvm_l_ortho, c="darkorange", lw=1.2, zorder=2)
        end
        if !isnothing(lon_rvm_h) && !isnothing(pa_rvm_h)
            plot(lon_rvm_h, pa_rvm_h,       c="darkorange", lw=1.2, ls="--", zorder=2)
            plot(lon_rvm_h, pa_rvm_h_ortho, c="darkorange", lw=1.2, ls="--", zorder=2)
        end

        # Inflection point — blue vertical line spanning the axes
        if !isnothing(phi0_l)
            axvline(phi0_l, c="tab:blue", lw=1.5, zorder=4)
        end
        if !isnothing(phi0_h)
            axvline(phi0_h, c="tab:blue", lw=1.5, ls="--", zorder=4)
        end

        ylabel("PA [deg]")
        tick_params(labelbottom=false)
        minorticks_on()

        # Bottom: Stokes profiles (normalized to peak I)
        subplot2grid((3, 1), (1, 0), rowspan=2)
        plot(lon_l, I_l   ./ maximum(I_l), c="black", lw=0.8)
        plot(lon_l, Lin_l ./ maximum(I_l), c="red",   lw=0.8)
        plot(lon_l, V_l   ./ maximum(I_l), c="blue",  lw=0.8)
        plot(lon_h, I_h   ./ maximum(I_h), c="black", lw=0.8, ls="--")
        plot(lon_h, Lin_h ./ maximum(I_h), c="red",   lw=0.8, ls="--")
        plot(lon_h, V_h   ./ maximum(I_h), c="blue",  lw=0.8, ls="--")
        xlabel("Longitude [deg]")
        ylabel("Flux Density [norm]")
        minorticks_on()

        savefig(joinpath(outdir, "position_angle.pdf"), dpi=200)
        if show_ show() end
    end


    # Auto-scale contour levels across the full range of the supplied maps,
    # snapped to a round 1/2/5·10^n step so labels read as e.g. 1000, 2000, 3000 km.
    function _snap_step(raw)
        raw <= 0 && return 0.0
        pw = 10.0 ^ floor(log10(raw))
        mant = raw / pw
        mant < 1.5 && return 1.0 * pw
        mant < 3.0 && return 2.0 * pw
        mant < 7.0 && return 5.0 * pw
        return 10.0 * pw
    end

    function _auto_height_levels(maps...)
        finite_h = filter(isfinite, vcat((vec(m) for m in maps)...))
        isempty(finite_h) && return Float64[]
        h_lo = minimum(finite_h)
        h_hi = maximum(finite_h)
        dh = _snap_step((h_hi - h_lo) / 8)
        if dh > 0
            return collect(ceil(h_lo / dh) * dh : dh : floor(h_hi / dh) * dh)
        else
            return Float64[]
        end
    end

    function geometry(chi2_map_l, chi2_map_h, alphas_deg, betas_deg, outdir;
                      show_=false, name_mod="PSR_NAME",
                      delta_chi2_max=nothing,
                      ndof_l=nothing, ndof_h=nothing,
                      chi2_red_max=10.0,
                      height_contours=nothing,
                      P_sec=nothing, W_deg_l=nothing, W_deg_h=nothing,
                      cmap="viridis")

        chi2_min_l = minimum(chi2_map_l[isfinite.(chi2_map_l)])
        chi2_min_h = minimum(chi2_map_h[isfinite.(chi2_map_h)])
        chi2_max_l = maximum(chi2_map_l[isfinite.(chi2_map_l)])
        chi2_max_h = maximum(chi2_map_h[isfinite.(chi2_map_h)])
        println("chi2_min: low=$(round(chi2_min_l, digits=2))  high=$(round(chi2_min_h, digits=2))")
        println("chi2_max: low=$(round(chi2_max_l, digits=2))  high=$(round(chi2_max_h, digits=2))")

        # Two display modes:
        #  (a) If ndof_l and ndof_h are given → plot reduced χ² with cutoff
        #      chi2_red_max (Johnston+ 2023 Fig. 1 convention: "χ² > 10" cutoff).
        #  (b) Else → plot Δχ² with cutoff delta_chi2_max (default 11.83 = 3σ
        #      for 2 free parameters).
        use_chi2_red = !isnothing(ndof_l) && !isnothing(ndof_h)
        if use_chi2_red
            dchi2_l = chi2_map_l ./ ndof_l
            dchi2_h = chi2_map_h ./ ndof_h
            dmax = chi2_red_max
            chi2_label = raw"$\chi^2_{\mathrm{red}}$"
            println("mode=chi2_red  ndof: low=$ndof_l  high=$ndof_h   cutoff=$(round(dmax,digits=2))")
            println("chi2_red_min: low=$(round(chi2_min_l/ndof_l, digits=3))  high=$(round(chi2_min_h/ndof_h, digits=3))")
        else
            dchi2_l = chi2_map_l .- chi2_min_l
            dchi2_h = chi2_map_h .- chi2_min_h
            dmax = isnothing(delta_chi2_max) ? 11.83 : delta_chi2_max
            chi2_label = raw"$\Delta\chi^2$"
            println("mode=delta_chi2   delta_chi2_max = $(round(dmax, digits=2))")
        end

        # Emission-height maps h(α,β) assuming a filled beam — one per band.
        # ρ from Gil, Gronkowski & Rudnicki (1984) / Johnston+ 2023 Eq. (3):
        #   cos ρ = cos α cos(α+β) + sin α sin(α+β) cos(W10/2)
        # Emission height (dipolar last-open field line; Johnston+ 2023 Eq. 2):
        #   h = 2 c P ρ² / (9π)   with ρ in radians
        function _height_map(W_deg)
            hmap = fill(NaN, length(alphas_deg), length(betas_deg))
            (isnothing(P_sec) || isnothing(W_deg)) && return hmap
            c_km_s = 2.99792458e5
            cosW2 = cos(deg2rad(W_deg) / 2)
            for i in eachindex(alphas_deg)
                a = deg2rad(alphas_deg[i])
                for j in eachindex(betas_deg)
                    b = deg2rad(betas_deg[j])
                    cr = cos(a) * cos(a + b) + sin(a) * sin(a + b) * cosW2
                    cr = clamp(cr, -1.0, 1.0)
                    rho = acos(cr)
                    hmap[i, j] = 2 * c_km_s * P_sec * rho^2 / (9 * pi)
                end
            end
            return hmap
        end

        if isnothing(P_sec) || (isnothing(W_deg_l) && isnothing(W_deg_h))
            @warn "Plot.geometry: P_sec and/or W10 not provided — height contours will not be drawn."
        end
        height_map_l = _height_map(W_deg_l)
        height_map_h = _height_map(W_deg_h)

        height_contours = isnothing(height_contours) ?
            _auto_height_levels(height_map_l, height_map_h) : height_contours
        println("height_contours [km] = $(round.(height_contours, digits=1))")

        # Mask cells outside the acceptable region (shown as white)
        dchi2_l_plot = Float64.(dchi2_l); dchi2_l_plot[dchi2_l .> dmax] .= NaN
        dchi2_h_plot = Float64.(dchi2_h); dchi2_h_plot[dchi2_h .> dmax] .= NaN

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        fig = figure(figsize=(7.0, 4.2))
        subplots_adjust(left=0.09, bottom=0.11, right=0.98, top=0.78, wspace=0.32)

        # chi2_map is [alpha_idx, beta_idx]; imshow wants [row=beta, col=alpha] → transpose
        extent = [alphas_deg[1], alphas_deg[end], betas_deg[1], betas_deg[end]]

        # --- Panel 1: low frequency ---
        ax1 = subplot(1, 2, 1)
        im1 = imshow(dchi2_l_plot',
                     origin="lower", extent=extent, aspect="auto",
                     cmap=cmap, vmin=0., vmax=dmax, interpolation="nearest")
        ax1.clabel(ax1.contour(alphas_deg, betas_deg, height_map_l',
                               levels=height_contours, colors="purple", linewidths=0.7),
                   fmt="%g", fontsize=6)
        ax1.set_xlabel("alpha [deg]")
        ax1.set_ylabel("beta [deg]")
        ax1.set_title("low frequency", pad=3)
        ax1.minorticks_on()

        # --- Panel 2: high frequency ---
        ax2 = subplot(1, 2, 2)
        imshow(dchi2_h_plot',
               origin="lower", extent=extent, aspect="auto",
               cmap=cmap, vmin=0., vmax=dmax, interpolation="nearest")
        ax2.clabel(ax2.contour(alphas_deg, betas_deg, height_map_h',
                               levels=height_contours, colors="purple", linewidths=0.7),
                   fmt="%g", fontsize=6)
        ax2.set_xlabel("alpha [deg]")
        ax2.set_ylabel("beta [deg]")
        ax2.set_title("high frequency", pad=3)
        ax2.minorticks_on()

        # Shared horizontal colorbar at the top
        cbar_ax = fig.add_axes([0.09, 0.86, 0.89, 0.05])
        cbar = fig.colorbar(im1, cax=cbar_ax, orientation="horizontal")
        cbar.set_label(chi2_label, labelpad=2)
        cbar_ax.xaxis.set_label_position("top")
        cbar_ax.xaxis.set_ticks_position("top")

        savepath = joinpath(outdir, "$(name_mod)_geometry.pdf")
        println(savepath)
        savefig(savepath)
        savefig(replace(savepath, "pdf" => "png"))

        if show_
            PyPlot.show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    """
        phase_drift(result, outdir, nbin; name_mod, show_)

    3-panel diagnostic for the phase-drift vs amplitude-modulation test.

    Panel 1: |L(φ)| – amplitude of complex LRFS slice at f3 vs longitude.
    Panel 2: ψ(φ) = arg(L(φ)) – phase vs longitude, with fitted slope line.
             Error bars (if result.p3_error > 0) come from the spread of
             ψ(φ) across the P3 values consistent with the quoted error.
    Panel 3: amplitude-null histogram with measured slope and ±1σ marked.

    Arguments:
      result  – NamedTuple from PhaseDrift.drift_test
      outdir  – output directory
      nbin    – total profile bins (for longitude axis scale)
    """
    function phase_drift(result, outdir, nbin::Int;
                         name_mod="pulsar", show_=false)
        on_bins = result.on_bins
        L_on    = result.L_on
        null    = result.null

        npts    = last(on_bins) - first(on_bins) + 1
        lon     = collect(range(-180.0 * npts / nbin, 180.0 * npts / nbin, length=npts))
        phi     = collect(on_bins)
        psi_deg = rad2deg.(angle.(L_on))
        null_deg     = rad2deg.(null)
        null_std_deg = std(null_deg)
        slope_deg    = rad2deg(result.slope)

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.5, 5.5))
        subplots_adjust(left=0.18, bottom=0.08, right=0.97, top=0.93, hspace=0.50)

        # Panel 1: Amplitude
        subplot(3, 1, 1)
        plot(lon, abs.(L_on), c="black", lw=0.8)
        xlabel("longitude (\$^\\circ\$)")
        ylabel("\$|L|\$")
        minorticks_on()

        # Panel 2: Phase
        subplot(3, 1, 2)
        psi_err_deg = rad2deg.(result.p3err_sigma)
        if any(psi_err_deg .> 0)
            errorbar(lon, psi_deg, yerr=psi_err_deg, fmt=".", ms=2, c="black",
                     ecolor="grey", elinewidth=0.4, capsize=0, zorder=3)
        else
            plot(lon, psi_deg, ".", ms=2, c="black", zorder=3)
        end
        plot(lon, rad2deg.(result.slope .* (phi .- mean(phi))), "-", c="red", lw=1.0,
             label=@sprintf("%.1f\$\\sigma\$", result.significance))
        xlabel("longitude (\$^\\circ\$)")
        ylabel("\$\\psi\$ (\$^\\circ\$)")
        ylim(-200, 200)
        minorticks_on()
        legend(fontsize=6, loc="upper right")

        # Panel 3: Null distribution
        subplot(3, 1, 3)
        hist(null_deg, bins=60, density=true, color="grey", alpha=0.7)
        axvline(x=slope_deg,      color="red",       lw=1.2, label="\$s_{\\mathrm{obs}}\$")
        axvline(x=null_std_deg,   color="steelblue", lw=0.8, ls="--")
        axvline(x=-null_std_deg,  color="steelblue", lw=0.8, ls="--", label="\$\\pm1\\sigma\$")
        xlabel("slope (\$^\\circ\$/bin)")
        ylabel("density")
        minorticks_on()
        legend(fontsize=6, loc="upper right")

        if result.p3_error > 0
            suptitle(@sprintf("P3 = %.1f \$\\pm\$ %.1f P0  |  %.1f\$\\sigma\$",
                              result.p3, result.p3_error, result.significance), fontsize=7)
        else
            suptitle(@sprintf("P3 = %.1f P0  |  %.1f\$\\sigma\$",
                              result.p3, result.significance), fontsize=7)
        end

        savepath = joinpath(outdir, "$(name_mod)_phase_drift.pdf")
        savefig(savepath)
        savefig(replace(savepath, ".pdf" => ".png"))
        println(savepath)

        if show_
            PyPlot.show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    """
    Phase-drift diagnostic, stability variant — same top two panels as
    `phase_drift` (|L| and ψ vs longitude), but the bottom panel shows the
    *local* phase gradient dψ/dφ instead of the null-slope histogram.

    That panel answers a question the single coherent slope cannot: does the
    phase change systematically, or does it jump? Pure amplitude modulation
    gives a piecewise-constant ψ with a 180° step at every node of |L|, so
    dψ/dφ sits at zero except for isolated spikes; a genuine drift gives a
    flat dψ/dφ at the global slope. Increments whose bins fall below
    snr_min·σ_off are drawn in light grey and excluded from χ² — the phase
    there is essentially random.

    Writes `<name_mod>_phase_stability.pdf/.png`, so it does not overwrite
    the `phase_drift` output.

    Arguments:
      result  – NamedTuple from PhaseDrift.drift_test
      outdir  – output directory
      nbin    – total profile bins (for longitude axis scale)
    """
    function phase_stability(result, outdir, nbin::Int;
                             name_mod="pulsar", show_=false)
        on_bins = result.on_bins
        L_on    = result.L_on

        npts    = last(on_bins) - first(on_bins) + 1
        lon     = collect(range(-180.0 * npts / nbin, 180.0 * npts / nbin, length=npts))
        phi     = collect(on_bins)
        psi_deg = rad2deg.(angle.(L_on))
        null_std_deg = std(rad2deg.(result.null))
        slope_deg    = rad2deg(result.slope)

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.5, 5.5))
        subplots_adjust(left=0.18, bottom=0.08, right=0.97, top=0.93, hspace=0.50)

        # Panel 1: Amplitude
        subplot(3, 1, 1)
        plot(lon, abs.(L_on), c="black", lw=0.8)
        xlabel("longitude (\$^\\circ\$)")
        ylabel("\$|L|\$")
        minorticks_on()

        # Panel 2: Phase
        subplot(3, 1, 2)
        psi_err_deg = rad2deg.(result.p3err_sigma)
        if any(psi_err_deg .> 0)
            errorbar(lon, psi_deg, yerr=psi_err_deg, fmt=".", ms=2, c="black",
                     ecolor="grey", elinewidth=0.4, capsize=0, zorder=3)
        else
            plot(lon, psi_deg, ".", ms=2, c="black", zorder=3)
        end
        plot(lon, rad2deg.(result.slope .* (phi .- mean(phi))), "-", c="red", lw=1.0,
             label=@sprintf("%.1f\$\\sigma\$", result.significance))
        xlabel("longitude (\$^\\circ\$)")
        ylabel("\$\\psi\$ (\$^\\circ\$)")
        ylim(-200, 200)
        minorticks_on()
        legend(fontsize=6, loc="upper right")

        # Panel 3: Local phase gradient — is the drift systematic or a jump?
        subplot(3, 1, 3)
        lon_mid  = 0.5 .* (lon[1:end-1] .+ lon[2:end])
        dpsi_deg = rad2deg.(result.dpsi)
        derr_deg = rad2deg.(result.dpsi_sigma)
        good     = result.dpsi_good
        bad      = .!good
        axhline(y=0.0, color="grey", lw=0.4, ls=":")
        # low-S/N increments: phase is essentially random there, no error bars
        if any(bad)
            plot(lon_mid[bad], dpsi_deg[bad], ".", ms=2, c="lightgrey", zorder=2)
        end
        if any(good)
            errorbar(lon_mid[good], dpsi_deg[good], yerr=derr_deg[good], fmt=".",
                     ms=2.5, c="black", ecolor="grey", elinewidth=0.4, capsize=0,
                     zorder=3)
        end
        chi2_lab = isnan(result.chi2_red) ? "\$\\chi^2_r\$ = n/a" :
            @sprintf("\$\\chi^2_r\$ = %.1f (%d)", result.chi2_red, result.chi2_dof)
        axhline(y=slope_deg, color="red", lw=1.0, label=chi2_lab)
        axhspan(slope_deg - null_std_deg, slope_deg + null_std_deg,
                color="red", alpha=0.10, lw=0)
        xlabel("longitude (\$^\\circ\$)")
        ylabel("\$d\\psi/d\\varphi\$ (\$^\\circ\$/bin)")
        ylim(-190, 190)
        minorticks_on()
        legend(fontsize=6, loc="upper right")

        if result.p3_error > 0
            suptitle(@sprintf("P3 = %.1f \$\\pm\$ %.1f P0  |  %.1f\$\\sigma\$",
                              result.p3, result.p3_error, result.significance), fontsize=7)
        else
            suptitle(@sprintf("P3 = %.1f P0  |  %.1f\$\\sigma\$",
                              result.p3, result.significance), fontsize=7)
        end

        savepath = joinpath(outdir, "$(name_mod)_phase_stability.pdf")
        savefig(savepath)
        savefig(replace(savepath, ".pdf" => ".png"))
        println(savepath)

        if show_
            PyPlot.show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    """
    Sliding-window drift diagnostic (`PhaseDrift.drift_test_sliding`) — the
    reversal-tolerant counterpart of `phase_drift`.

    Panel 1: |L_b(φ)|/σ_b map, longitude vs pulse number — where in time and
    longitude the f3 modulation power lives.
    Panel 2: per-window drift slope arg g_b vs pulse number. Windows below the
    per-window significance threshold are light grey (their phase is noise);
    detected windows are black with 1σ error bars. For a drifter that reverses
    direction this slope(t) curve flips sign between episodes — the direct
    measurement the global `phase_drift` slope averages away. Structure
    narrower than the window length is smoothing artefact, not signal.
    Panel 3: null distribution of S = Σ_b |Im g_b| (flat-phase surrogates with
    shared noise) against the observed S.

    Writes `<name_mod>_phase_drift_sliding.pdf/.png`, so `phase_drift` and
    `phase_stability` outputs are kept.

    Arguments:
      result  – NamedTuple from PhaseDrift.drift_test_sliding
      outdir  – output directory
      nbin    – total profile bins (for longitude axis scale)
    """
    function phase_drift_sliding(result, outdir, nbin::Int;
                                 name_mod="pulsar", show_=false)
        npts      = length(result.on_bins)
        lon       = collect(range(-180.0 * npts / nbin, 180.0 * npts / nbin, length=npts))
        centers   = result.centers
        slope_deg = rad2deg.(result.slope)
        err_deg   = rad2deg.(result.slope_err)
        det       = result.detected
        nd        = .!det

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.5, 5.5))
        subplots_adjust(left=0.18, bottom=0.08, right=0.97, top=0.93, hspace=0.50)

        # Panel 1: modulation power map
        subplot(3, 1, 1)
        imshow(permutedims(result.Lmap), origin="lower", aspect="auto",
               extent=(centers[1], centers[end], lon[1], lon[end]),
               cmap="viridis", interpolation="nearest")
        xlabel("pulse number")
        ylabel("longitude (\$^\\circ\$)")

        # Panel 2: time-resolved drift slope — sign flips are drift reversals
        subplot(3, 1, 2)
        axhline(y=0.0, color="grey", lw=0.4, ls=":")
        if any(nd)
            plot(centers[nd], slope_deg[nd], ".", ms=1.5, c="lightgrey", zorder=2)
        end
        if any(det)
            errorbar(centers[det], slope_deg[det], yerr=err_deg[det], fmt=".",
                     ms=2, c="black", ecolor="grey", elinewidth=0.3, capsize=0,
                     zorder=3)
            lim = min(190.0, 1.5 * maximum(abs.(slope_deg[det]) .+ err_deg[det]))
            ylim(-lim, lim)
        else
            ylim(-190, 190)
        end
        xlim(centers[1], centers[end])
        xlabel("pulse number")
        ylabel("slope (\$^\\circ\$/bin)")
        minorticks_on()

        # Panel 3: incoherent drift statistic vs its null
        subplot(3, 1, 3)
        hist(result.S_null, bins=60, density=true, color="grey", alpha=0.7)
        axvline(x=result.S, color="red", lw=1.2,
                label=@sprintf("\$S_{\\mathrm{obs}}\$ (%.1f\$\\sigma\$)",
                               result.significance))
        xlabel("\$S = \\sum_b |\\mathrm{Im}\\,g_b|\$")
        ylabel("density")
        minorticks_on()
        legend(fontsize=6, loc="upper right")

        suptitle(@sprintf("P3 = %.1f P0  |  W = %d  |  %.1f\$\\sigma\$",
                          result.p3, result.window, result.significance), fontsize=7)

        savepath = joinpath(outdir, "$(name_mod)_phase_drift_sliding.pdf")
        savefig(savepath)
        savefig(replace(savepath, ".pdf" => ".png"))
        println(savepath)

        if show_
            PyPlot.show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    """
    Windowed phase-stability diagnostic (`PhaseDrift.drift_test_profile`) —
    the longitude-resolved counterpart of `phase_stability`, for pulsars whose
    global LRFS bin is emptied by P3 wobble.

    Panel 1: on-pulse S/N profile, mean over windows of |L_b(φ)|/σ_b (pure
    noise sits at ~1.25) — says which longitudes carry modulation at all, and
    therefore which increments below are worth reading.
    Panel 2: ψ(φ) rebuilt by integrating the measured increments, against the
    fitted constant gradient. This is the panel that answers the question by
    eye: a drift follows the straight line, a phase step between components
    shows up as a staircase.
    Panel 3: dψ/dφ(φ) with 1σ surrogate errors and the fitted gradient. Bins
    whose phase error exceeds `sigma_max` are drawn light grey and are not
    worth interpreting.

    Writes `<name_mod>_phase_stability_windowed.pdf/.png`, leaving the
    `phase_stability` output intact.

    Arguments:
      result    – NamedTuple from PhaseDrift.drift_test_profile
      outdir    – output directory
      nbin      – total profile bins (for longitude axis scale)
      sigma_max – phase error [rad] above which an increment is greyed out
    """
    function phase_stability_windowed(result, outdir, nbin::Int;
                                      sigma_max=0.5, name_mod="pulsar", show_=false)
        npts    = length(result.on_bins)
        lon     = collect(range(-180.0 * npts / nbin, 180.0 * npts / nbin, length=npts))
        lon_mid = 0.5 .* (lon[1:end-1] .+ lon[2:end])
        x       = collect(1:npts-1)

        dpsi_deg = rad2deg.(result.dpsi)
        derr_deg = rad2deg.(result.dpsi_sigma)
        slope_deg = rad2deg(result.slope)
        good = result.dpsi_sigma .<= sigma_max
        bad  = .!good

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.5, 5.5))
        subplots_adjust(left=0.18, bottom=0.08, right=0.97, top=0.93, hspace=0.50)

        # Panel 1: where is there modulation to measure?
        subplot(3, 1, 1)
        plot(lon, result.amp, c="black", lw=0.8)
        axhline(y=1.25, color="grey", lw=0.4, ls=":")
        xlabel("longitude (\$^\\circ\$)")
        ylabel("\$\\langle|L_b|\\rangle/\\sigma_b\$")
        minorticks_on()

        # Panel 2: integrated phase — straight line = drift, staircase = jump.
        # Referenced to the reliable bins only: cumsum integrates the noisy
        # edge increments too, and letting those set the offset and the y-range
        # would squash the part that carries the answer.
        subplot(3, 1, 2)
        ref     = any(good) ? good : trues(length(x))
        psi_deg = rad2deg.(result.psi_cum)
        psi_deg = psi_deg .- mean(psi_deg[ref])
        fit_deg = slope_deg .* (x .- mean(x[ref]))
        psi_rel = copy(psi_deg)
        psi_rel[.!ref] .= NaN
        plot(lon_mid, psi_deg, "-", c="lightgrey", lw=0.6, zorder=2)
        plot(lon_mid, psi_rel, "-", c="black", lw=0.9, zorder=3,
             label="\$\\int d\\psi\$")
        plot(lon_mid, fit_deg, "-", c="red", lw=1.0, zorder=4,
             label=@sprintf("%.2f\$^\\circ\$/bin", slope_deg))
        lo_, hi_ = extrema(vcat(psi_deg[ref], fit_deg[ref]))
        pad = 0.18 * max(hi_ - lo_, 1.0)
        ylim(lo_ - pad, hi_ + pad)
        xlabel("longitude (\$^\\circ\$)")
        ylabel("\$\\psi\$ (\$^\\circ\$)")
        minorticks_on()
        legend(fontsize=6, loc="upper left")

        # Panel 3: the increments themselves and the constant-gradient fit
        subplot(3, 1, 3)
        axhline(y=0.0, color="grey", lw=0.4, ls=":")
        if any(bad)
            plot(lon_mid[bad], dpsi_deg[bad], ".", ms=2, c="lightgrey", zorder=2)
        end
        if any(good)
            errorbar(lon_mid[good], dpsi_deg[good], yerr=derr_deg[good], fmt=".",
                     ms=2.5, c="black", ecolor="grey", elinewidth=0.4, capsize=0,
                     zorder=3)
        end
        ptxt = result.p_value == 0 ? "p<1/n" : @sprintf("p=%.3g", result.p_value)
        axhline(y=slope_deg, color="red", lw=1.0,
                label=@sprintf("\$\\chi^2_r\$ = %.1f (%s)", result.chi2_red, ptxt))
        xlabel("longitude (\$^\\circ\$)")
        ylabel("\$d\\psi/d\\varphi\$ (\$^\\circ\$/bin)")
        ylim(-190, 190)
        minorticks_on()
        legend(fontsize=6, loc="upper right")

        suptitle(@sprintf("P3 = %.1f P0  |  W = %d  |  \$\\chi^2_r\$ = %.1f",
                          result.p3, result.window, result.chi2_red), fontsize=7)

        savepath = joinpath(outdir, "$(name_mod)_phase_stability_windowed.pdf")
        savefig(savepath)
        savefig(replace(savepath, ".pdf" => ".png"))
        println(savepath)

        if show_
            PyPlot.show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    """
    Turn one psrcat record (parameter => value) into a (P, Pdot) entry.

    A large fraction of the ATNF catalogue lists the spin frequency instead of
    the period, hence
        P = 1 / F0        Pdot = -F1 / F0^2
    Records without both P and Pdot defined are dropped.
    """
    function _psrcat_record!(names, periods, pdots, name, rec)
        f0 = get(rec, "F0", NaN)
        p0 = get(rec, "P0", NaN)
        # whichever of the two is missing, derive it from the other
        if !isfinite(f0) && isfinite(p0) && p0 > 0
            f0 = 1 / p0
        elseif !isfinite(p0) && isfinite(f0) && f0 > 0
            p0 = 1 / f0
        end
        p1 = get(rec, "P1", NaN)
        if !isfinite(p1)
            f1 = get(rec, "F1", NaN)
            isfinite(f1) && isfinite(f0) && f0 > 0 && (p1 = -f1 / f0^2)
        end
        if isfinite(p0) && isfinite(p1) && p0 > 0
            push!(names, name)
            push!(periods, p0)
            push!(pdots, p1)
        end
    end


    """
    Read (name, P, Pdot) from an ATNF psrcat database file (psrcat tarball).

    The file is a sequence of records separated by lines starting with "@-";
    within a record each line is "PARAM value [error] [reference]" and lines
    starting with "#" are comments. Only the first occurrence of a parameter is
    used (some quantities are listed several times, once per reference).
    """
    function read_psrcat(filename)
        isfile(filename) || error("psrcat database not found: $filename")

        names = String[]
        periods = Float64[]
        pdots = Float64[]

        rec = Dict{String,Float64}()
        jname = ""
        bname = ""

        for line in eachline(filename)
            if startswith(line, "@-")
                # the J name identifies the pulsar; PSRB is only a fallback, and
                # it may well come first in a record
                _psrcat_record!(names, periods, pdots, isempty(jname) ? bname : jname, rec)
                empty!(rec)
                jname = ""
                bname = ""
                continue
            end
            (isempty(strip(line)) || startswith(line, "#")) && continue
            parts = split(line)
            length(parts) < 2 && continue
            key, val = parts[1], parts[2]
            if key == "PSRJ"
                isempty(jname) && (jname = val)
            elseif key == "PSRB"
                isempty(bname) && (bname = val)
            elseif key in ("P0", "P1", "F0", "F1") && !haskey(rec, key)
                v = tryparse(Float64, val)
                isnothing(v) || (rec[key] = v)
            end
        end
        # last record (no trailing @-)
        _psrcat_record!(names, periods, pdots, isempty(jname) ? bname : jname, rec)

        return names, periods, pdots
    end


    """
    Visible part of the power law Pdot = A * P^m inside the (xlims, ylims) box.
    Returns (xs, ys), or nothing when the line does not cross the box at all.
    """
    function _ppdot_segment(A, m, xlims, ylims)
        xa = (ylims[1] / A)^(1 / m)
        xb = (ylims[2] / A)^(1 / m)
        xlo = max(xlims[1], min(xa, xb))
        xhi = min(xlims[2], max(xa, xb))
        xlo < xhi || return nothing
        xs = 10 .^ range(log10(xlo), log10(xhi), length=32)
        return xs, A .* xs .^ m
    end


    """
    On-screen slope (degrees) of a power law with exponent m in a log-log plot,
    used to rotate the line labels. Depends on the axes aspect ratio, hence the
    axes position within the figure is needed.
    """
    function _ppdot_angle(ax, m, xlims, ylims, figsize)
        bbox = ax.get_position()
        w_in = figsize[1] * bbox.width
        h_in = figsize[2] * bbox.height
        return atand(m * (h_in / log10(ylims[2] / ylims[1])) /
                         (w_in / log10(xlims[2] / xlims[1])))
    end


    # Pulsars renamed in the catalogue: a J name encodes the position, so it has
    # to follow whenever the position is corrected. Our input lists still carry
    # the old names, hence the translation on the way in.
    #   J1402-5124 -> J1402-5021 (= B1359-51, Molonglo MP 1359). Found in the
    #   first Molonglo survey (Vaughan & Large 1970) with a declination wrong by
    #   about a degree; re-localised to 14:02:56.0, -50:21:43 by Lower et al.
    #   (2020, MNRAS 494, 228), which forced the new J name.
    const PSR_RENAMED = Dict("J1402-5124" => "J1402-5021")


    """
    Read a list of pulsars to mark: one PSRJ name per line, blank lines and
    lines starting with "#" ignored, anything after the name on a line is
    dropped (so files with extra columns work as well). Names listed in
    PSR_RENAMED are translated to the ones the catalogue uses now.
    """
    function _read_pulsar_list(filename)
        selected = Set{String}()
        if !isfile(filename)
            @warn "pulsar list not found: $filename"
            return selected
        end
        for line in eachline(filename)
            s = strip(line)
            (isempty(s) || startswith(s, "#")) && continue
            name = String(first(split(s)))
            if haskey(PSR_RENAMED, name)
                println("  renamed pulsar: $name -> $(PSR_RENAMED[name])")
                name = PSR_RENAMED[name]
            end
            push!(selected, name)
        end
        return selected
    end


    """
    Read the component offsets measured between two observing frequencies
    (input/offsets.csv), one pulsar per line:

        name, offset 1, offset 1 err, ..., offset 4, offset 4 err, grade, comment

    Offset i is the shift of the i-th component, mu(1523 MHz) - mu(1023 MHz) in
    degrees of longitude, with the components ordered by longitude (see
    GaussianFit.component_offsets). For a pulsar with two or more components the
    quantity returned is the change of the component separation between the two
    frequencies,

        value = offset_last - offset_first = separation(1523) - separation(1023)

    which is free of any shift common to the whole profile: a residual DM error
    moves all components by the same amount and cancels in the difference
    (Hassall et al. 2012). Negative means the profile is narrower at the higher
    frequency, i.e. ordinary radius-to-frequency mapping. A pulsar with a single
    component gives no such handle, so its raw shift is kept instead and flagged
    with ncomp = 1 — the two are different quantities and must not be mixed.

    Returns Dict(name => (value, err, ncomp, grade)), value and err in degrees,
    grade the 1-10 quality mark of the fit.
    """
    function _read_offsets(filename)
        out = Dict{String,NamedTuple{(:value, :err, :ncomp, :grade),
                                     Tuple{Float64,Float64,Int,Int}}}()
        if !isfile(filename)
            @warn "offset file not found: $filename"
            return out
        end
        for (nline, line) in enumerate(eachline(filename))
            nline == 1 && continue  # header
            s = strip(line)
            (isempty(s) || startswith(s, "#")) && continue
            f = split(s, ',')
            length(f) < 10 && continue
            name = String(strip(f[1]))
            if haskey(PSR_RENAMED, name)
                println("  renamed pulsar: $name -> $(PSR_RENAMED[name])")
                name = PSR_RENAMED[name]
            end
            off, err = Float64[], Float64[]
            for c in 1:4
                v = tryparse(Float64, strip(f[2c]))
                e = tryparse(Float64, strip(f[2c+1]))
                (isnothing(v) || isnothing(e)) && continue
                push!(off, v)
                push!(err, e)
            end
            isempty(off) && continue
            grade = something(tryparse(Int, strip(f[10])), 0)
            out[name] = length(off) > 1 ?
                (value = off[end] - off[1],
                 err = sqrt(err[end]^2 + err[1]^2), ncomp = length(off), grade = grade) :
                (value = off[1], err = err[1], ncomp = 1, grade = grade)
        end
        return out
    end


    """ Label a power of ten as 10^n, anything else as a plain number. """
    function _pow10_label(v)
        e = log10(v)
        isapprox(e, round(e), atol=1e-9) && return "\$10^{$(Int(round(e)))}\$"
        return @sprintf("%.3g", v)
    end


    """
    P-Pdot diagram of the ATNF catalogue, shared by ppdot and ppdot_offsets —
    call one of those instead.

    Pulsars with Pdot <= 0 (mostly globular-cluster pulsars, where acceleration
    in the cluster potential dominates the observed spin-down) are skipped —
    they cannot be shown on a logarithmic axis.

    When `offsets` points to a file (see _read_offsets), the pulsars with a
    measured frequency offset are drawn on top of the population, colour-coded
    by that offset. Measurements graded below `min_grade` are dropped and those
    below `offset_nsigma` are drawn as open symbols. The symlog colour scale
    runs to `offset_vmax`, by default the `offset_vmax_quantile` quantile of
    |offset| so that a few outliers do not flatten the rest of the sample;
    whatever falls outside saturates and the colorbar gets its arrowheads.
    With `offsets = nothing` nothing of that is drawn and the plain diagram
    comes out.
    """
    function _ppdot(outdir; catalogue=normpath(joinpath(@__DIR__, "..", "input", "psrcat.db")),
                   name_mod="psrcat", show_=true,
                   b_lines=[1e10, 1e12, 1e14],     # surface magnetic field [G]
                   age_lines=[1e3, 1e6, 1e9],      # characteristic age [yr]
                   edot_lines=[1e30, 1e33, 1e36],  # spin-down luminosity [erg/s]
                   death_bp2=0.17e12,              # death line B/P^2 [G s^-2]
                   inertia=1e45,                   # moment of inertia [g cm^2]
                   highlight=normpath(joinpath(@__DIR__, "..", "input", "pulsars.txt")),
                   highlight_label="selected",
                   offsets=nothing,                # component offsets, see ppdot_offsets
                   min_grade=6,                    # lowest acceptable quality mark (1-10)
                   offset_nsigma=3.0,              # below this the symbol is left open
                   offset_linthresh=0.1,           # linear core of the symlog colour scale [deg]
                   offset_vmax=nothing,            # colour range [deg]; nothing -> from offset_vmax_quantile
                   offset_vmax_quantile=0.90,      # quantile of |offset| that sets the colour range
                   offset_cmap="RdBu_r",
                   plims=(1e-3, 2e2), pdotlims=(1e-22, 1e-8))

        names, periods, pdots = read_psrcat(catalogue)
        println("psrcat: $(length(names)) pulsars with P and Pdot ($catalogue)")

        keep = pdots .> 0
        println("psrcat: $(count(.!keep)) pulsars skipped (Pdot <= 0)")
        nam = names[keep]
        p = periods[keep]
        pd = pdots[keep]

        # pulsars to mark on top of the whole population
        selected = isnothing(highlight) ? Set{String}() : _read_pulsar_list(highlight)
        marked = [n in selected for n in nam]
        if !isempty(selected)
            println("psrcat: $(count(marked)) of $(length(selected)) pulsars from $(basename(highlight)) marked")
            absent = sort(collect(setdiff(selected, Set(nam[marked]))))
            if !isempty(absent)
                shown = join(absent[1:min(10, end)], ", ")
                println("  missing (no P/Pdot in the catalogue, or Pdot <= 0): $(length(absent)) — $shown$(length(absent) > 10 ? ", ..." : "")")
            end
        end

        rc("font", size=11.)
        rc("axes", linewidth=0.7)
        rc("lines", linewidth=0.7)

        figsize = (7.08661, 6.29921)  # 18cm x 16cm
        figure(figsize=figsize, frameon=true)
        subplots_adjust(left=0.13, bottom=0.10, right=0.98, top=0.98)

        ax = gca()
        xscale("log")
        yscale("log")
        xlim(plims)
        ylim(pdotlims)

        # Reference lines, all power laws Pdot = A P^m (straight in log-log):
        #   B     = 3.2e19 sqrt(P Pdot) G     ->  Pdot = (B / 3.2e19)^2 / P
        #   tau_c = P / (2 Pdot)              ->  Pdot = P / (2 tau_c)
        #   Edot  = 4 pi^2 I Pdot / P^3       ->  Pdot = Edot P^3 / (4 pi^2 I)
        # The death line follows the constant polar-cap potential drop
        # B / P^2 = const (Bhattacharya et al. 1992); pulsars below and to the
        # right of it are not expected to emit in the radio band.
        #   B / P^2 = death_bp2               ->  Pdot = (death_bp2 / 3.2e19)^2 P^3
        c_b, c_age, c_edot, c_death = "tab:blue", "tab:green", "tab:orange", "tab:red"
        yr_s = 3.15576e7

        function line!(A, m, txt, frac, color, ls; va="bottom")
            seg = _ppdot_segment(A, m, plims, pdotlims)
            isnothing(seg) && return
            xs, ys = seg
            ax.plot(xs, ys, ls=ls, lw=0.9, c=color, alpha=0.8, zorder=1)
            isnothing(txt) && return
            xt = 10^(log10(xs[1]) + frac * log10(xs[end] / xs[1]))
            ax.text(xt, A * xt^m, txt, fontsize=11, color=color,
                    ha="center", va=va, rotation_mode="anchor", zorder=4,
                    rotation=_ppdot_angle(ax, m, plims, pdotlims, figsize))
        end

        # the label fractions keep each family near a different edge of the box;
        # the B labels hang below their lines to stay clear of the upper frame
        for b in b_lines
            line!((b / 3.2e19)^2, -1, _pow10_label(b), 0.10, c_b, "--"; va="top")
        end
        for tau in age_lines
            line!(1 / (2 * tau * yr_s), 1, _pow10_label(tau), 0.95, c_age, "-.")
        end
        for edot in edot_lines
            line!(edot / (4 * pi^2 * inertia), 3, _pow10_label(edot), 0.85, c_edot, ":")
        end
        line!((death_bp2 / 3.2e19)^2, 3, nothing, 0.5, c_death, "-")

        plot(p, pd, ".", ms=2.8, c="black", mec="none", zorder=3)
        c_sel = "magenta"
        any(marked) && plot(p[marked], pd[marked], ".", ms=4.0, c=c_sel, mec="none", zorder=4)

        # Frequency offsets on top of the population. The colour carries the
        # change of the component separation between 1023 and 1523 MHz (circles,
        # >= 2 components) or the raw shift of the only component (diamonds);
        # the scale is symlog because the offsets span three decades. Open
        # symbols mark measurements below offset_nsigma.
        off = isnothing(offsets) ? Dict{String,Any}() : _read_offsets(offsets)
        idx = Int[]
        if !isempty(off)
            in_cat = Set(nam)
            absent = sort([k for k in keys(off) if !(k in in_cat)])
            faint = sort([k for (k, v) in off if v.grade < min_grade])
            idx = [i for i in eachindex(nam)
                   if haskey(off, nam[i]) && off[nam[i]].grade >= min_grade]
            println("offsets: $(length(idx)) of $(length(off)) pulsars from $(basename(offsets)) plotted")
            isempty(faint) || println("  below grade $min_grade: $(join(faint, ", "))")
            isempty(absent) || println("  missing (no P/Pdot in the catalogue, or Pdot <= 0): $(join(absent, ", "))")
        end

        if !isempty(idx)
            vals = [off[nam[i]].value for i in idx]
            sig = [abs(off[nam[i]].value) / off[nam[i]].err for i in idx]
            multi = [off[nam[i]].ncomp > 1 for i in idx]

            # The colour range is set by a quantile of |offset| rather than by
            # the maximum: three quarters of the sample sit below 1 deg, while
            # the largest offsets reach ~12 deg (J1733-3716, flagged in the
            # input as needing a redo), so the maximum spends a good part of
            # the scale on a handful of points. Measurements beyond +-vmax
            # saturate, which the arrowheads on the colorbar mark.
            vmax = isnothing(offset_vmax) ?
                ceil(quantile(abs.(vals), offset_vmax_quantile) * 10) / 10 :
                Float64(offset_vmax)
            vmax = max(vmax, 10 * offset_linthresh)  # at least one decade of colour
            nsat = count(>(vmax), abs.(vals))
            println(@sprintf("offsets: colour range +-%.4g deg, %d of %d measurements saturated",
                             vmax, nsat, length(vals)))
            cmap = PyPlot.matplotlib.pyplot.get_cmap(offset_cmap)
            cnorm = PyPlot.matplotlib.colors.SymLogNorm(linthresh=offset_linthresh,
                                                        vmin=-vmax, vmax=vmax, base=10)

            function offset_points!(sel, marker, ms, filled)
                any(sel) || return
                xs, ys, vs = p[idx[sel]], pd[idx[sel]], vals[sel]
                if filled
                    ax.scatter(xs, ys, s=ms, marker=marker, c=vs, cmap=cmap, norm=cnorm,
                               edgecolors="black", linewidths=0.4, zorder=5)
                else
                    ax.scatter(xs, ys, s=ms, marker=marker, facecolors="none",
                               edgecolors=cmap(cnorm(vs)), linewidths=0.9, zorder=5)
                end
            end
            strong = sig .>= offset_nsigma
            for filled in (true, false)
                keep_ = filled ? strong : .!strong
                offset_points!(multi .& keep_, "o", 26.0, filled)
                offset_points!(.!multi .& keep_, "D", 18.0, filled)
            end

            # the colorbar goes inside the axes: an external one would resize
            # them and invalidate the label rotations set by _ppdot_angle above
            cax = ax.inset_axes([0.04, 0.86, 0.38, 0.022])
            sm = PyPlot.matplotlib.cm.ScalarMappable(norm=cnorm, cmap=cmap)
            sm.set_array([])
            # Ticks: the two ends plus a 1-2-5 ladder between the linear core
            # and vmax, thinned greedily from the core outwards. The spacing
            # has to be judged in norm positions rather than in degrees — on a
            # symlog bar the last decade crowds the end as soon as vmax comes
            # near it (with linthresh 0.1 and vmax 2.2 the 1 deg tick sits only
            # 0.07 of the bar width away from 2.2), and plain decades then
            # leave the bar with no labels at all between its ends.
            inner = Float64[]
            for d in sort([m * 10.0^k for k in -3:2 for m in (1.0, 2.0, 5.0)])
                (d >= 0.999 * offset_linthresh && d < vmax) || continue
                cnorm(vmax) - cnorm(d) > 0.09 || continue
                isempty(inner) || cnorm(d) - cnorm(inner[end]) > 0.09 || continue
                push!(inner, d)
            end
            ticks = [-vmax; -reverse(inner); inner; vmax]
            cb = colorbar(sm, cax=cax, orientation="horizontal", ticks=ticks,
                          extend=(nsat > 0 ? "both" : "neither"))
            cb.ax.set_xticklabels([@sprintf("%g", t) for t in ticks], fontsize=6)
            cb.ax.tick_params(length=2, pad=1)
            cb.outline.set_linewidth(0.5)
            cb.set_label("1523 \$-\$ 1023 MHz offset (\$^\\circ\$)", fontsize=7, labelpad=2)
        end

        xlabel("\$P\$ (s)")
        ylabel("\$\\dot{P}\$ (s s\$^{-1}\$)")
        minorticks_on()

        L2D = PyPlot.matplotlib.lines.Line2D
        handles = Any[]
        any(marked) && push!(handles, L2D([], [], c=c_sel, ls="none", marker=".",
                                          ms=5.0, mec="none", label=highlight_label))
        if !isempty(idx)
            push!(handles, L2D([], [], mfc="0.6", ls="none", marker="o", ms=5.0,
                               mec="black", mew=0.4, label="\$\\Delta\$ separation (\$\\geq\$2 comp.)"))
            push!(handles, L2D([], [], mfc="0.6", ls="none", marker="D", ms=4.0,
                               mec="black", mew=0.4, label="shift (1 comp.)"))
            push!(handles, L2D([], [], mfc="none", ls="none", marker="o", ms=5.0,
                               mec="0.4", mew=0.9,
                               label=@sprintf("\$< %g\\sigma\$", offset_nsigma)))
        end
        push!(handles, L2D([], [], c=c_b, ls="--", lw=0.9, label="\$B\$ (G)"))
        push!(handles, L2D([], [], c=c_age, ls="-.", lw=0.9, label="\$\\tau_c\$ (yr)"))
        push!(handles, L2D([], [], c=c_edot, ls=":", lw=0.9, label="\$\\dot{E}\$ (erg s\$^{-1}\$)"))
        push!(handles, L2D([], [], c=c_death, ls="-", lw=0.9, label="death line"))
        legend(handles=handles, fontsize=8, loc="lower right", framealpha=0.9,
               borderpad=0.4, handlelength=2.5, labelspacing=0.35)

        savepath = joinpath(outdir, "ppdot_$(name_mod).pdf")
        savefig(savepath)
        savefig(replace(savepath, ".pdf" => ".png"))
        println(savepath)

        if show_
            PyPlot.show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        close()
    end


    """
    P-Pdot diagram of the ATNF catalogue, with the pulsars listed in `highlight`
    marked. Keywords are those of _ppdot, except `offsets` and everything that
    goes with it — use ppdot_offsets for that. Writes ppdot_<name_mod>.pdf/png.
    """
    ppdot(outdir; kwargs...) = _ppdot(outdir; kwargs..., offsets=nothing)


    """
    P-Pdot diagram with the component offsets measured between two observing
    frequencies (see _read_offsets) drawn on top of the population: colour gives
    the offset on a symlog scale, circles the change of the component separation
    (>= 2 components), diamonds the shift of a single component, open symbols
    the measurements below `offset_nsigma`. Pulsars graded below `min_grade` are
    left out, and the colour range comes from `offset_vmax` /
    `offset_vmax_quantile`. The remaining keywords are those of _ppdot.
    """
    function ppdot_offsets(outdir;
                           offsets=normpath(joinpath(@__DIR__, "..", "input", "offsets.csv")),
                           name_mod="offsets", kwargs...)
        _ppdot(outdir; kwargs..., offsets=offsets, name_mod=name_mod)
    end


    # ------------------------------------------------------------------ LRFS: faza i sliding LRFS
    "Rows that hold data (zapped pulses are all-zero rows)."
    _live_rows(X) = [any(!iszero, @view X[i, :]) for i in 1:size(X, 1)]

    "Pulse stack with the per-bin mean removed; zapped (all-zero) rows set to 0 after demeaning."
    function _demeaned(X)
        good = _live_rows(X)
        Y = Float64.(X)
        any(good) || return Y, good
        m = mean(Y[good, :], dims=1)
        Y .-= m
        Y[.!good, :] .= 0.0
        return Y, good
    end

    "Coherent complex LRFS at frequency f [cycles/P] from all pulses, kernel e^(+2πi f n) (as PhaseDrift)."
    function _coherent_lrfs_row(pulses, f)
        Y, _ = _demeaned(pulses)
        carrier = [exp(2im * π * f * (n - 1)) for n in 1:size(Y, 1)]
        return vec(sum(Y .* carrier, dims=1))
    end

    """
    Phase track of one complex LRFS row L (all bins): on-pulse phase in degrees, unwrapped along longitude
    over significant bins (|L| > nsig·σ, σ = std of Re/Im of the off-pulse bins), errors σ/|L| [deg].
    """
    function _lrfs_phase_track(L, bin_st, bin_end; nsig=3.0)
        nb = length(L)
        off = vcat(1:bin_st-1, bin_end+1:nb)
        σ = isempty(off) ? NaN : std(vcat(real.(L[off]), imag.(L[off])))
        Lon = L[bin_st:bin_end]
        amp = abs.(Lon)
        sig = isfinite(σ) ? amp .> nsig * σ : falses(length(Lon))
        ph = rad2deg.(angle.(Lon))
        prev = NaN
        for j in eachindex(ph)
            sig[j] || continue
            if isfinite(prev)
                while ph[j] - prev > 180; ph[j] -= 360; end
                while ph[j] - prev < -180; ph[j] += 360; end
            end
            prev = ph[j]
        end
        err = rad2deg.(σ ./ max.(amp, eps()))
        return (phase=ph, err=err, sig=sig, amp=amp, sigma=σ)
    end

    "Off-pulse bins (as many as on-pulse) farthest from the pulse window, on the phase circle."
    function _far_offpulse(nb, bin_st, bin_end)
        c = (bin_st + bin_end) / 2
        dist(j) = min(abs(j - c), nb - abs(j - c))
        cand = setdiff(1:nb, max(1, bin_st - 20):min(nb, bin_end + 20))
        n = min(bin_end - bin_st + 1, length(cand))
        return sort(sort(cand, by=j -> -dist(j))[1:n])
    end

    """
        sliding_lrfs_compute(data, bin_st, bin_end; window=256, step=1)

    Sliding LRFS as the S2DFS of Serylak, Stappers & Weltevrede (2009, A&A 506, 865): a window of
    `window` pulses moved by `step` pulses (1 = every pulse); in each window the LRFS (FFT along the
    pulses of every on-pulse bin, window mean removed) is collapsed over longitude (Σ|F|²). The same
    for an equal number of off-pulse bins gives the noise level, so `snr` = on/off power per frequency.
    Returns (snr, freq, centres, fpeak, peak_snr): rows = windows (centre pulse, 1-based), columns =
    fluctuation frequency k/window, k = 1..window/2.
    """
    function sliding_lrfs_compute(data, bin_st, bin_end; window=256, step=1)
        N, nb = size(data)
        N < window && error("za mało impulsów ($N) na okno $window")
        X, _ = _demeaned(data)
        on = collect(bin_st:bin_end); off = _far_offpulse(nb, bin_st, bin_end)
        starts = collect(1:step:N-window+1)
        nf = window ÷ 2
        S = zeros(length(starts), nf); So = zeros(length(starts), nf)
        Bon = zeros(window, length(on)); Boff = zeros(window, length(off))
        Pon = plan_rfft(Bon, 1); Poff = plan_rfft(Boff, 1)
        for (b, s) in enumerate(starts)
            r = s:s+window-1
            Bon .= @view X[r, on];  Bon .-= mean(Bon, dims=1)
            Boff .= @view X[r, off]; Boff .-= mean(Boff, dims=1)
            Fo = Pon * Bon; Ff = Poff * Boff
            @inbounds for k in 1:nf
                S[b, k] = sum(abs2, @view Fo[k + 1, :])
                So[b, k] = sum(abs2, @view Ff[k + 1, :])
            end
        end
        noise = vec(mean(So, dims=2)) .* length(on) ./ max(length(off), 1)
        snr = S ./ max.(noise, eps())
        freq = collect(1:nf) ./ window
        kmin = 2                                      # pomijamy najnizsza czestotliwosc (powolne zmiany)
        fpeak = [freq[kmin - 1 + argmax(@view snr[b, kmin:end])] for b in eachindex(starts)]
        peak_snr = [maximum(@view snr[b, kmin:end]) for b in eachindex(starts)]
        centres = starts .+ (window - 1) / 2
        return (snr=snr, freq=freq, centres=centres, fpeak=fpeak, peak_snr=peak_snr, window=window, step=step)
    end

    """
        sliding_lrfs_phase(data, f3, bin_st, bin_end; window=256, step=1, nsig=3.0)

    Phase of the LRFS at f3 in every sliding window (kernel e^(+2πi f3 n) with the ABSOLUTE pulse number, so
    a steady drift gives the same phase pattern in every window; a P3 change tilts it, a drift reversal
    flips the slope). Cumulative sums: O(N·bins). Insignificant bins (|L| < nsig·σ_off) are NaN.
    """
    function sliding_lrfs_phase(data, f3, bin_st, bin_end; window=256, step=1, nsig=3.0)
        N, nb = size(data)
        X, _ = _demeaned(data)
        on = collect(bin_st:bin_end); off = _far_offpulse(nb, bin_st, bin_end)
        carrier = [exp(2im * π * f3 * (n - 1)) for n in 1:N]
        Con = cumsum(X[:, on] .* carrier, dims=1); Coff = cumsum(X[:, off] .* carrier, dims=1)
        starts = collect(1:step:N-window+1)
        ph = fill(NaN, length(starts), length(on))
        for (b, s) in enumerate(starts)
            e = s + window - 1
            Lon = s == 1 ? Con[e, :] : Con[e, :] .- Con[s-1, :]
            Lof = s == 1 ? Coff[e, :] : Coff[e, :] .- Coff[s-1, :]
            σ = std(vcat(real.(Lof), imag.(Lof)))
            ok = abs.(Lon) .> nsig * σ
            ph[b, ok] .= rad2deg.(angle.(Lon[ok]))
        end
        return (phase=ph, centres=starts .+ (window - 1) / 2)
    end

    """
Podzial okna impulsu na skladowe: minimum profilu (wygladzonego) miedzy dwoma najwyzszymi maksimami (> 20%
szczytu, co najmniej 10 binow od siebie). Faza przez przerwe miedzy skladowymi jest niejednoznaczna o 360°
(Basu & Mitra 2018; TPA XX), wiec nachylenia liczymy osobno w kazdej skladowej (jak ich tabela 2).
"""
    function _profile_components(prof)
        ps = [mean(prof[max(1, j - 2):min(length(prof), j + 2)]) for j in eachindex(prof)]
        mx = [j for j in 2:length(ps)-1 if ps[j] >= ps[j-1] && ps[j] >= ps[j+1] && ps[j] > 0.2maximum(ps)]
        length(mx) < 2 && return [1:length(prof)]
        top = sort(mx, by=j -> -ps[j])
        q = findfirst(j -> abs(j - top[1]) >= 10, top[2:end])
        isnothing(q) && return [1:length(prof)]
        lo, hi = minmax(top[1], top[1 + q])
        cut = lo - 1 + argmin(ps[lo:hi])
        return [1:cut, cut+1:length(prof)]
    end

    """
Mapa fazy (okna x biny, stopnie, NaN = nieistotne) rozwinieta: w kazdym oknie 0 w binie `jref` (szczyt profilu),
potem kazdy punkt w obrocie najblizszym wspolnej sciezki (srednia kolowa po oknach, rozwinieta od szczytu profilu).
Okna bez fazy w szczycie profilu - puste.
"""
    function _unwrap_phase_map(P, jref)
        U = fill(NaN, size(P))
        for b in axes(P, 1)
            isfinite(P[b, jref]) || continue
            for j in axes(P, 2)
                isfinite(P[b, j]) && (U[b, j] = mod(P[b, j] - P[b, jref] + 180, 360) - 180)
            end
        end
        ref = [let v = filter(isfinite, U[:, j]); length(v) < 0.2 * size(U, 1) ? NaN : rad2deg(angle(sum(cis.(deg2rad.(v))))) end
               for j in axes(U, 2)]
        for dir in (1, -1)
            prev = 0.0; j = jref + dir
            while 1 <= j <= length(ref)
                if isfinite(ref[j]); ref[j] -= 360 * round((ref[j] - prev) / 360); prev = ref[j]; end
                j += dir
            end
        end
        isfinite(ref[jref]) || (ref[jref] = 0.0)
        for j in axes(U, 2)
            isfinite(ref[j]) || continue
            for b in axes(U, 1)
                isfinite(U[b, j]) && (U[b, j] += 360 * round((ref[j] - U[b, j]) / 360))
            end
        end
        return U
    end

    """
        sliding_lrfs_tracks(data, bin_st, bin_end; window=256, step=1, nsig=3.0)

    Method of Basu & Mitra (2018, MNRAS 475, 5098), with a step of `step` pulses (they used 50): the peak
    frequency fp and its width δf = FWHM/2.355 come from the time-averaged collapsed LRFS; then in every
    window and at every on-pulse longitude the largest LRFS amplitude within fp ± 3δf is taken, kept only
    if > nsig × the off-pulse rms, with its phase (kernel e^(+2πi f n): rising phase = drift to later
    longitudes). The window start is arbitrary w.r.t. P3, so each window's phases are set to 0 at the
    profile peak (as in the paper) and unwrapped outwards from it. Per window a straight line is fitted
    to the significant phases (slope in deg/deg) - its spread is the temporal variation of the drift.
    """
    function sliding_lrfs_tracks(data, bin_st, bin_end; window=256, step=1, nsig=3.0, p3=nothing, per_bin=false,
                                 common_unwrap=true)
        # V2 (jak Basu & Mitra 2018): per_bin=true - w kazdym binie najsilniejsza czestotliwosc w fp +- 3df, faza 0
        # w szczycie profilu; obrot fazy najblizszy sredniej (publikacja tego nie opisuje, ale ich zwarte chmury
        # punktow tego wymagaja; rozwijanie osobno w kazdym oknie dawalo rozrzut 124-152° zamiast 15-43° - test
        # syntetyczny). common_unwrap=false zostaje jako opcja.
        # V1 (nasza): jedna czestotliwosc na okno + rozwiniecie wzgledem wspolnej sciezki (bez kopii o 360°).
        # per_bin=true: jak w publikacji - w kazdym binie dlugosci osobno najsilniejsza czestotliwosc w fp +- 3df
        # (Basu & Mitra 2018); false: jedna czestotliwosc na okno (bez skokow fazy miedzy binami)
        N, nb = size(data)
        X, good = _demeaned(data)
        on = collect(bin_st:bin_end); off = _far_offpulse(nb, bin_st, bin_end)
        nf = window ÷ 2
        freq = collect(1:nf) ./ window
        starts = collect(1:step:N-window+1)
        Bon = zeros(window, length(on)); Boff = zeros(window, length(off))
        Pon = plan_rfft(Bon, 1); Poff = plan_rfft(Boff, 1)
        # 1. srednie widmo (zlozone po dlugosci) -> fp, FWHM
        S = zeros(nf); So = zeros(nf)
        for s in starts
            r = s:s+window-1
            Bon .= @view X[r, on]; Bon .-= mean(Bon, dims=1)
            Boff .= @view X[r, off]; Boff .-= mean(Boff, dims=1)
            Fo = Pon * Bon; Ff = Poff * Boff
            for k in 1:nf
                S[k] += sum(abs2, @view Fo[k + 1, :]); So[k] += sum(abs2, @view Ff[k + 1, :])
            end
        end
        spec = S ./ max.(So .* length(on) ./ max(length(off), 1), eps())
        kp = 1 + argmax(@view spec[2:end])
        # znane P3: szczyt przy 1/P3 (±10%, co najmniej ±1 bin) - inaczej nulle / powolne zmiany jasnosci
        # wygrywaja najnizsza czestotliwoscia (J0034, J1727, J0304: "P3 = 128")
        if p3 isa Number && p3 > 0
            w = [k for k in 1:nf if abs(freq[k] - 1 / p3) <= max(0.1 / p3, 1 / window)]
            isempty(w) || (kp = w[argmax(spec[w])])
        end
        base = median(spec); half = base + (spec[kp] - base) / 2
        lo = kp; while lo > 1 && spec[lo - 1] > half; lo -= 1; end
        hi = kp; while hi < nf && spec[hi + 1] > half; hi += 1; end
        # FWHM z liniowa interpolacja przejsc przez polowe wysokosci (inaczej najmniej 1 bin = 1/window)
        xl = lo > 1 ? (lo - 1) + (half - spec[lo - 1]) / max(spec[lo] - spec[lo - 1], eps()) : Float64(lo)
        xh = hi < nf ? hi + (spec[hi] - half) / max(spec[hi] - spec[hi + 1], eps()) : Float64(hi)
        fwhm = (xh - xl) / window
        δf = fwhm / 2.355
        band = [k for k in 1:nf if abs(freq[k] - freq[kp]) <= 3δf + 0.5 / window]
        # S jak w Basu+2016: widmo znormalizowane do 1 w szczycie, S = (szczyt - tlo) / FWHM
        Sfac = (1 - base / spec[kp]) / max(fwhm, eps())
        # 2. amplituda i faza w kazdym oknie przy fp
        prof = vec(mean(Float64.(data[good, on]), dims=1))
        jref = argmax(prof)
        A = fill(NaN, length(starts), length(on)); Φ = fill(NaN, length(starts), length(on))
        slope = fill(NaN, length(starts))
        dlon = 360 / nb
        for (b, s) in enumerate(starts)
            r = s:s+window-1
            Bon .= @view X[r, on]; Bon .-= mean(Bon, dims=1)
            Boff .= @view X[r, off]; Boff .-= mean(Boff, dims=1)
            Fo = Pon * Bon; Ff = Poff * Boff
            rms = std(vcat(real.(vec(Ff[band .+ 1, :])), imag.(vec(Ff[band .+ 1, :]))))
            # jedna czestotliwosc na okno (najsilniejsza w pasmie po calym oknie impulsu): rozne k w sasiednich
            # binach daja skoki fazy (sasiednie biny DFT maja przesunieta faze) - test syntetyczny
            k = band[argmax([sum(abs2, @view Fo[k + 1, :]) for k in band])]
            amp = zeros(length(on)); ph = zeros(length(on))
            for j in eachindex(on)
                kj = per_bin ? band[argmax([abs(Fo[kk + 1, j]) for kk in band])] : k
                c = conj(Fo[kj + 1, j])                      # conj: kernel e^(+2πi f n)
                amp[j] = abs(c); ph[j] = rad2deg(angle(c))
            end
            sig = amp .> nsig * rms
            sig[jref] || continue                            # brak odniesienia w tym oknie
            ph .-= ph[jref]
            for dir in (1, -1)                               # rozwiniecie od szczytu profilu na zewnatrz
                prev = 0.0
                j = jref + dir
                while 1 <= j <= length(on)
                    if sig[j]
                        while ph[j] - prev > 180; ph[j] -= 360; end
                        while ph[j] - prev < -180; ph[j] += 360; end
                        prev = ph[j]
                    end
                    j += dir
                end
            end
            A[b, sig] .= amp[sig] ./ rms; Φ[b, sig] .= ph[sig]
            js = findall(sig)
            if length(js) >= 3
                x = (js .- jref) .* dlon; y = ph[js]
                slope[b] = sum((x .- mean(x)) .* (y .- mean(y))) / max(sum((x .- mean(x)) .^ 2), eps())
            end
        end
        # Rozwiniecie wzgledem WSPOLNEJ sciezki: kazde okno rozwijane osobno gubilo w slabym miejscu miedzy
        # skladowymi pelny obrot, wiec czesc okien lądowala w rownoleglych "kopiach", a srednia wypadala miedzy
        # nimi (J1555). Sciezka odniesienia = srednia kolowa po oknach w kazdym binie, rozwinieta od szczytu
        # profilu; kazdy punkt przesuwany o wielokrotnosc 360° najblizej tej sciezki; nachylenia od nowa.
        if common_unwrap && any(isfinite, Φ)
            ref = [let v = filter(isfinite, Φ[:, j]); length(v) < 0.2 * size(Φ, 1) ? NaN : rad2deg(angle(sum(cis.(deg2rad.(v))))) end
                   for j in axes(Φ, 2)]
            for dir in (1, -1)
                prev = isfinite(ref[jref]) ? ref[jref] : 0.0
                j = jref + dir
                while 1 <= j <= length(ref)
                    if isfinite(ref[j])
                        ref[j] -= 360 * round((ref[j] - prev) / 360); prev = ref[j]
                    end
                    j += dir
                end
            end
            for j in axes(Φ, 2)
                isfinite(ref[j]) || continue
                for b in axes(Φ, 1)
                    isfinite(Φ[b, j]) && (Φ[b, j] += 360 * round((ref[j] - Φ[b, j]) / 360))
                end
            end
            for b in axes(Φ, 1)
                js = findall(isfinite, @view Φ[b, :])
                length(js) >= 3 || (slope[b] = NaN; continue)
                x = (js .- jref) .* dlon; y = Φ[b, js]
                slope[b] = sum((x .- mean(x)) .* (y .- mean(y))) / max(sum((x .- mean(x)) .^ 2), eps())
            end
        end
        # nachylenia osobno w kazdej skladowej (przez przerwe faza niejednoznaczna o 360°)
        comps = _profile_components(prof)
        cslope = fill(NaN, length(starts), length(comps))
        for (c, r) in enumerate(comps), b in axes(Φ, 1)
            js = [j for j in r if isfinite(Φ[b, j])]
            length(js) >= 3 || continue
            x = js .* dlon; y = Φ[b, js]
            cslope[b, c] = sum((x .- mean(x)) .* (y .- mean(y))) / max(sum((x .- mean(x)) .^ 2), eps())
        end
        return (amp=A, phase=Φ, slope=slope, cslope=cslope, comps=comps, fp=freq[kp], fwhm=fwhm, S=Sfac, jref=jref, spec=spec, freq=freq,
                starts=starts, window=window, step=step)
    end

    # Styl wykresow odtwarzanych z publikacji (V2, V3): ramka z kreskami do srodka na wszystkich bokach, czcionka
    # szeryfowa jak w MNRAS. PUB_CMAP_BASU: mapa kolorow lewego rysunku Basu & Mitra (do ustawienia wg publikacji).
    const PUB_CMAP_BASU = Ref("afmhot")          # jak w Basu & Mitra (czarny-czerwony-zolty-bialy)
    function _pub_style()
        rc("font", family="serif", size=8)
        rc("mathtext", fontset="dejavuserif")
        rc("axes", linewidth=0.7)
        rc("xtick", direction="in", top=true); rc("ytick", direction="in", right=true)
        rc("xtick.minor", visible=true); rc("ytick.minor", visible=true)
        rc("lines", linewidth=0.8)
    end
    _pub_style_off() = PyPlot.matplotlib.rcdefaults()

    """
Plot like Basu & Mitra (2018) Fig. 2 right: amplitude (top) and phase (middle) at fp vs longitude for all
windows (dots) + mean (black), profile (bottom). Longitude 0 = profile peak, as in the paper. Phase in the
paper's convention (= minus ours): rising phase = subpulses move towards EARLIER longitudes ("negative
drifting"). Several peaks (`ts`) are drawn in different colours, as their Fig. 4 for J1555-3134.
"""
    function sliding_lrfs_tracks_plot(ts, labels, prof, lon, outfile; title_="", maxpts=150, sign=-1, pub=true, note="")
        # wyglad jak Basu & Mitra (2018) rys. 4 (dla V1 i V2): chmury punktow okien w kolorach pikow (czerwony, zielony),
        # srednie czarne (pik 1: pelne kola, pik 2: puste kwadraty), opis pikow kolorowym tekstem, profil pod spodem.
        # sign = -1: konwencja Basu & Mitra (rosnaca faza = dryf ku wczesniejszym dlugosciom); +1: nasza
        _pub_style()
        o = sortperm([t.fp for t in ts])                         # pik 1 = nizsza czestotliwosc (jak w publikacji)
        ts = ts[o]
        cols = ["red", "limegreen", "royalblue"]
        mk = [(m="o", fc="black"), (m="s", fc="none"), (m="^", fc="none")]
        amax = maximum([let v = filter(isfinite, vec(mean(replace(t.amp, NaN => 0.0), dims=1))); isempty(v) ? 1.0 : maximum(v) end for t in ts])
        figure(figsize=(4.2, 6.2))
        subplots_adjust(left=0.18, bottom=0.08, right=0.96, top=0.93, hspace=0.0)
        suptitle(title_, fontsize=7.5)
        axs = []
        for (panel, key, lab) in ((0, :amp, "FFT Amp."), (2, :phase, "FFT Phase (deg)"))
            ax = subplot2grid((5, 1), (panel, 0), rowspan=2); push!(axs, ax)
            for (i, t) in enumerate(ts)
                M = key == :amp ? t.amp ./ amax : sign .* t.phase
                rows = findall(b -> any(isfinite, @view M[b, :]), axes(M, 1))
                # srednia i punkty okien tylko w binach, gdzie faza jest mierzalna w >= polowie okien (jak w publikacji:
                # pusto miedzy skladowymi i poza impulsem); rysowane najwyzej `maxpts` okien (srednia ze wszystkich)
                m = [let v = filter(isfinite, M[rows, j]); length(v) < 0.5 * length(rows) ? NaN : mean(v) end for j in axes(M, 2)]
                keep = isfinite.(m)
                sel = length(rows) > maxpts ? rows[round.(Int, range(1, length(rows), length=maxpts))] : rows
                for b in sel
                    v = M[b, :]; v[.!keep] .= NaN
                    scatter(lon, v, s=1.5, c=cols[mod1(i, 3)], alpha=0.5, lw=0)
                end
                k = mk[mod1(i, 3)]
                scatter(lon, m, s=k.m == "o" ? 7 : 12, marker=k.m, facecolors=k.fc, edgecolors="black", lw=0.6, zorder=3)
            end
            for t in ts[1:1], r in t.comps[1:end-1]
                axvline(lon[last(r)] + (lon[2] - lon[1]) / 2, c="grey", ls=":", lw=0.6)     # granica skladowych
            end
            ylabel(lab); xlim(lon[1], lon[end]); tick_params(labelbottom=false)
            if key == :phase
                for (i, t) in enumerate(ts)
                    text(0.04, 0.93 - 0.08(i - 1), @sprintf("Peak %d = %.3f cy/P", i, t.fp), transform=ax.transAxes,
                         color=cols[mod1(i, 3)], fontsize=8, va="top")
                end
                isempty(note) || text(0.96, 0.04, note, transform=ax.transAxes, fontsize=6, ha="right", va="bottom", color="grey")
            end
        end
        subplot2grid((5, 1), (4, 0))
        plot(lon, prof ./ maximum(prof), c="black"); xlim(lon[1], lon[end]); ylim(-0.05, 1.1)
        xlabel("Pulse Longitude (deg)")
        savefig(outfile); println(outfile)
        close("all"); _pub_style_off()
    end

    # V2 tracks tak jak w wersji 8b46134 (zaakceptowany wyglad): male, blade punkty wszystkich okien, srednie czarne/
    # granatowe, opisy osi i tytul z liczbami, konwencja Basu & Mitra.
    function sliding_lrfs_tracks_plot_v2(ts, labels, prof, lon, outfile; title_="", maxpts=3000, sign=-1)
        # sign = -1: konwencja Basu & Mitra (rosnaca faza = dryf ku wczesniejszym dlugosciom); +1: nasza
        cols = ["red", "royalblue", "darkgreen"]; mcols = ["black", "navy", "olive"]
        rc("font", size=7); rc("axes", linewidth=0.5)
        figure(figsize=(4.5, 7))
        subplots_adjust(left=0.16, bottom=0.07, right=0.97, top=0.92, hspace=0.0)
        suptitle(title_, fontsize=6.5)
        for (panel, key, lab) in ((0, :amp, "amplitude at fp (× rms)"), (2, :phase, sign < 0 ? "phase at fp (°), 0 at profile peak\n(Basu & Mitra convention)" :
                                                                                 "phase at fp (°), 0 at profile peak\n(ours: rising = later longitudes)"))
            subplot2grid((5, 1), (panel, 0), rowspan=2)
            for (i, t) in enumerate(ts)
                M = key == :amp ? t.amp : sign .* t.phase
                rows = findall(b -> any(isfinite, @view M[b, :]), axes(M, 1))
                sel = length(rows) > maxpts ? rows[round.(Int, range(1, length(rows), length=maxpts))] : rows
                for b in sel; scatter(lon, M[b, :], s=0.5, c=cols[mod1(i, 3)], alpha=0.12, lw=0); end
                # srednia tylko tam, gdzie bin jest istotny w >= polowie okien (inaczej pojedyncze szumowe
                # przekroczenia 3 sigma daja czarne punkty poza oknem impulsu - J0820, J1034)
                m = [let v = filter(isfinite, M[rows, j]); length(v) < 0.5 * length(rows) ? NaN : mean(v) end for j in axes(M, 2)]
                scatter(lon, m, s=3, c=mcols[mod1(i, 3)], label=labels[i])
            end
            for r in ts[1].comps[1:end-1]
                axvline(lon[last(r)] + (lon[2] - lon[1]) / 2, c="grey", ls="--", lw=0.6)    # granica skladowych
            end
            ylabel(lab); xlim(lon[1], lon[end]); tick_params(labelbottom=false)
            panel == 0 && length(ts) > 1 && legend(fontsize=6, markerscale=2)
        end
        subplot2grid((5, 1), (4, 0))
        plot(lon, prof ./ maximum(prof), c="grey"); xlim(lon[1], lon[end]); xlabel("longitude (°), 0 = profile peak"); yticks([])
        savefig(outfile); println(outfile)
        close("all")
    end

    "Second peak of the mean spectrum (Basu & Mitra 2018, J1555-3134): local maximum away from fp (> 3 FWHM,
    > 4 bins), not at the lowest frequencies (nulls), at least `frac` of the main peak height above the floor."
    function _second_peak(spec, freq, kp, fwhm, window; frac=0.3)
        base = median(spec); nf = length(spec)
        sep = max(4, ceil(Int, 3 * fwhm * window))
        cand = [k for k in 3:nf-1 if abs(k - kp) > sep && spec[k] >= spec[k-1] && spec[k] >= spec[k+1] &&
                spec[k] - base >= frac * (spec[kp] - base)]
        isempty(cand) && return nothing
        return freq[cand[argmax(spec[cand])]]
    end

    """
        sliding_lrfs(data, outdir, params; window=256, step=1, name_mod="pulsar", show_=false)

    Like Basu & Mitra (2018, MNRAS 475, 5098) Figs. 2-5, with the window moved by `step` pulses (1 =
    continuous; they used 50, which they note makes horizontal stripes). Page 1 (left figure of the paper):
    LRFS averaged over the pulse window vs window START pulse, normalised so that the time-averaged peak = 1;
    extra panels: pulse energy (nulls) on the left, phase at fp per window on the right, mean spectrum and
    profile at the bottom. Page 2 (right figure): amplitude and phase at fp for every window vs longitude
    (0 = profile peak) - for each peak if the spectrum has a second one. Files: <name_mod>_sliding_lrfs.pdf,
    _lrfs_tracks.pdf, _sliding_lrfs.csv (window, strongest frequency), _lrfs_tracks.txt (slopes).
    """
    function sliding_lrfs(data, outdir, params; window=256, step=1, name_mod="pulsar", show_=false, nsig=3.0, style=:nasza,
                          nulls=Int[], track_sign=nothing)
        # track_sign: konwencja fazy na wykresie tracks i w nachyleniach (-1 Basu & Mitra, +1 nasza = jak mapa fazy);
        # domyslnie -1 dla :basu, +1 dla :nasza; wersja koncowa: +1
        # nulls: numery impulsow (od 0) uznanych za nulle - czerwone kreski na panelu jasnosci (wersja :nasza)
        # style = :nasza (moc wzgledem szumu, srodek okna, nasza konwencja fazy) albo :basu (Basu & Mitra 2018:
        # sredni szczyt = 1, pierwszy impuls okna, ich konwencja fazy, + strona z amplituda i faza przy fp)
        basu = style == :basu
        tag = basu ? "v2_basu_mitra" : "v1_nasza"
        p = params
        bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
        N, nb = size(data)
        p3 = get(p, "p3", nothing)
        r = sliding_lrfs_compute(data, bst, bend; window=window, step=step)
        t1 = sliding_lrfs_tracks(data, bst, bend; window=window, step=step, nsig=nsig, p3=p3, per_bin=basu, common_unwrap=!basu)
        f2 = _second_peak(t1.spec, t1.freq, round(Int, t1.fp * window), t1.fwhm, window)
        t2 = isnothing(f2) ? nothing : sliding_lrfs_tracks(data, bst, bend; window=window, step=step, nsig=nsig, p3=1 / f2,
                                                           per_bin=basu, common_unwrap=!basu)
        q = sliding_lrfs_phase(data, t1.fp, bst, bend; window=window, step=step, nsig=nsig)
        mkpath(outdir)
        open(joinpath(outdir, "$(name_mod)_$(tag)_sliding.csv"), "w") do io
            println(io, "window_start,f_peak,P3_peak,peak_snr")
            for b in eachindex(r.centres)
                @printf(io, "%d,%.6f,%.4f,%.3f\n", round(Int, r.centres[b] - (window - 1) / 2) - 1, r.fpeak[b], 1 / r.fpeak[b], r.peak_snr[b])
            end
        end
        _, good = _demeaned(data)
        E = vec(sum(Float64.(data[:, bst:bend]), dims=2)); E[.!good] .= NaN
        prof = vec(mean(Float64.(data[good, bst:bend]), dims=1))
        dl = 360 * (bend - bst + 1) / nb
        lon = collect(range(-dl / 2, dl / 2, length=bend - bst + 1))
        lon0 = lon .- lon[t1.jref]                               # 0 = szczyt profilu (jak w publikacji)
        st0 = r.centres .- (window - 1) / 2 .- 1                 # pierwszy impuls okna (od 0)
        avg = vec(mean(r.snr, dims=1))
        kp = clamp(round(Int, t1.fp * window), 1, length(avg))
        img = basu ? r.snr ./ max(avg[kp], eps()) : r.snr        # basu: sredni szczyt = 1; nasza: moc / szum
        yv = basu ? st0 : r.centres .- 1
        c1, c2 = yv[1], yv[end]

        if basu
            # jak Basu & Mitra (2018) rys. 2-5 lewy: lewy panel - sila widma w oknie (niebieska, os odwrocona),
            # mapa w kolorach "hot", dolny panel - srednie widmo (czerwone, max = 1), opis w tytule
            _pub_style()
            figure(figsize=(4.6, 5.4))
            subplots_adjust(left=0.15, bottom=0.1, right=0.97, top=0.9, wspace=0.0, hspace=0.0)
            # sila dryfu w oknie: wysokosc najsilniejszego szczytu widma (moc/szum), znormalizowana do 1
            pw = copy(r.peak_snr); pw ./= max(maximum(pw), eps())
            suptitle(@sprintf("PSR %s   nfft=%d   Peak(LRFS)=%.3f", name_mod, window, t1.fp), fontsize=8)
            ax0 = subplot2grid((4, 4), (0, 0), rowspan=3)
            plot(pw, yv, c="blue", lw=0.8); ylim(c1, c2); xlim(1.02, max(0.0, minimum(pw) - 0.05))
            ax0.xaxis.set_ticks_position("top"); ax0.xaxis.set_label_position("top"); xlabel("peak (norm.)", fontsize=7)
            ylabel("start period no.")
            subplot2grid((4, 4), (0, 1), rowspan=3, colspan=3, sharey=ax0)
            imshow(img, origin="lower", aspect="auto", cmap=PUB_CMAP_BASU[], interpolation="nearest",
                   extent=[0, r.freq[end], c1, c2], vmin=0, vmax=max(1.0, quantile(vec(img), 0.995)))
            tick_params(labelleft=false, labelbottom=false, labeltop=true)
            subplot2grid((4, 4), (3, 1), colspan=3)
            plot(r.freq, avg ./ max(maximum(avg), eps()), c="red", lw=0.8); xlim(0, r.freq[end])
            xlabel("frequency (1/P)")
            f = joinpath(outdir, "$(name_mod)_$(tag)_sliding.pdf")
            savefig(f); println(f)
            close("all"); _pub_style_off()
        else
            rc("font", size=7); rc("axes", linewidth=0.5); rc("lines", linewidth=0.5)
            figure(figsize=(7.5, 7))
            subplots_adjust(left=0.08, bottom=0.08, right=0.97, top=0.9, wspace=0.0, hspace=0.0)
            suptitle("$(name_mod): sliding LRFS, okno $window impulsów, krok $step" * (basu ? " (jak Basu & Mitra 2018)" : " (nasza wersja)"), fontsize=8)
            ax0 = subplot2grid((5, 7), (0, 0), rowspan=4)
            plot(E, 0:N-1, c="grey", lw=0.3)
            if !isempty(nulls)                                     # nulle: czerwone kreski przy lewej krawedzi
                xe = minimum(filter(isfinite, E))
                scatter(fill(xe, length(nulls)), nulls, marker="_", c="red", s=12, lw=0.8)
            end
            ylim(c1, c2); xticks([]); ylabel(basu ? "start period of the window" : "pulse number (window centre)")
            subplot2grid((5, 7), (0, 1), rowspan=4, colspan=3, sharey=ax0)
            imshow(img, origin="lower", aspect="auto", cmap="viridis", interpolation="nearest",
                   extent=[0, r.freq[end], c1, c2], vmin=0, vmax=basu ? max(1.5, quantile(vec(img), 0.995)) : quantile(vec(img), 0.995))
            tick_params(labelleft=false, labelbottom=false)
            axph = subplot2grid((5, 7), (0, 4), rowspan=4, colspan=3, sharey=ax0)
            # faza ROZWINIETA: w kazdym oknie 0 w szczycie profilu, bez zawijania co 360° (wspolna sciezka jak w tracks);
            # skala niebieski (-) / czerwony (+), poziomice co 180° - ich gestosc = nachylenie fazy (P2 = 360°/nachylenie)
            U = _unwrap_phase_map(q.phase, t1.jref)
            vm = let v = filter(isfinite, U); isempty(v) ? 180.0 : max(180.0, quantile(abs.(v), 0.99)) end
            imu = imshow(U, origin="lower", aspect="auto", cmap="coolwarm", interpolation="nearest",   # NaN = puste (biale)
                         extent=[lon0[1], lon0[end], c1, c2], vmin=-vm, vmax=vm)
            if count(isfinite, U) > 10
                lv = collect(-ceil(vm / 180) * 180:180:ceil(vm / 180) * 180)
                try
                    contour(lon0, yv, U, levels=lv, colors="k", linewidths=0.3, alpha=0.6)
                catch
                end
            end
            tick_params(labelleft=false, labelbottom=false)
            title(@sprintf("faza przy fp = %.4f (P3 = %.2f P)\n0 w szczycie profilu, rozwinięta", t1.fp, 1 / t1.fp), fontsize=6.5)
            cax = axph.inset_axes([0.88, 0.62, 0.04, 0.33])
            cb = colorbar(imu, cax=cax); cb.ax.tick_params(labelsize=5); cb.set_label("faza (°)", fontsize=5)
            subplot2grid((5, 7), (4, 1), colspan=3)
            plot(r.freq, basu ? avg ./ max(avg[kp], eps()) : avg, c="grey")
            isnothing(f2) || axvline(f2, c="royalblue", ls=":", lw=0.6)
            axvline(t1.fp, c="red", ls=":", lw=0.6)
            xlim(0, r.freq[end]); xlabel("fluctuation frequency (cycles/P)"); ylabel(basu ? "mean" : "mean (S/N)")
            subplot2grid((5, 7), (4, 4), colspan=3)
            plot(lon0, prof ./ maximum(prof), c="grey"); xlim(lon0[1], lon0[end]); xlabel("longitude (°), 0 = profile peak"); yticks([])
            f = joinpath(outdir, "$(name_mod)_$(tag)_sliding.pdf")
            savefig(f); println(f)
            show_ && (show(); println("Enter - dalej"); readline(stdin; keep=false))
        end
        close("all")

        sg = isnothing(track_sign) ? (basu ? -1 : 1) : track_sign
        function _info(t)
            sl = filter(isfinite, t.slope)
            cs = join([let v = filter(isfinite, t.cslope[:, c]); length(v) < 3 ? "skł.$c -" : @sprintf("skł.%d %.1f ± %.1f", c, sg * mean(v), std(v)) end
                       for c in axes(t.cslope, 2)], ", ")
            @sprintf("fp = %.4f (P3 = %.2f P), FWHM = %.4f, S = %.1f, okien z fazą %d z %d, nachylenie fazy w składowych: %s °/° (%s)",
                     t.fp, 1 / t.fp, t.fwhm, t.S, length(sl), length(t.starts), cs,
                     sg < 0 ? "konwencja Basu & Mitra: + = dryf ku wcześniejszym długościom" : "nasza: + = dryf ku późniejszym długościom")
        end
        ts = isnothing(t2) ? [t1] : [t1, t2]
        labs = [@sprintf("P3 = %.2f P", 1 / t.fp) for t in ts]
        info = "params P3 = " * ((p3 isa Number && p3 > 0) ? @sprintf("%.2f", p3) : "brak") * " | " *
               join(["szczyt $i: " * _info(t) for (i, t) in enumerate(ts)], " | ")
        println("  ", info)
        basu && sliding_lrfs_tracks_plot_v2([t1; isnothing(t2) ? [] : [t2]], labs, prof, lon0,
                                            joinpath(outdir, "$(name_mod)_$(tag)_tracks.pdf"); sign=sg,
                                            title_="$(name_mod): okno $window, krok $step (jak Basu & Mitra 2018, rys. 2)
" *
                                                   join([_info(t) for t in ts], "
"))
        basu || sliding_lrfs_tracks_plot(ts, labs, prof, lon0, joinpath(outdir, "$(name_mod)_$(tag)_tracks.pdf"); sign=1,
                                 note=basu ? "" : "rising phase = drift to later longitudes",
                                 title_="PSR $(name_mod)   nfft=$window   Nperiod=$N" * (basu ? "" : "   (krok $step)"))
        open(joinpath(outdir, "$(name_mod)_$(tag)_tracks.txt"), "w") do io
            println(io, info)
            println(io, sg < 0 ? "peak,window_start,slope_deg_per_deg_BasuMitra" : "peak,window_start,slope_deg_per_deg_ours")
            for (i, t) in enumerate(ts), b in eachindex(t.starts)
                @printf(io, "%d,%d,%.5f\n", i, t.starts[b] - 1, sg * t.slope[b])
            end
        end
        return (spectrum=r, tracks=ts, info=info)
    end

    """
        lrfs_phase_blocks(data, bin_st, bin_end, p3s; block=64, snr_min=7.0, frac_min=60/64, nboot=200, nsig=3.0)

    Subpulse phase and modulation amplitude as in TPA XX (Hsu+2026, MNRAS 553, stag1935; Edwards &
    Stappers 2002/2003, Weltevrede 2016): independent blocks of `block` pulses; only blocks where at least
    `frac_min` of the pulses have S/N >= `snr_min` (TPA XX: 60 of 64 with S/N >= 7). For each feature
    (P3 in `p3s`) the FFT bin nearest to block/P3 (± 1 bin, strongest in the mean spectrum) is used.
    Phases of different blocks are aligned with one constant offset per block (iteratively, to the
    amplitude-weighted mean track, as Edwards & Stappers 2003); mean track = circular mean over blocks,
    errors from a bootstrap over blocks. Modulation amplitude = 2|F|/block, so a pure sinusoid has the
    amplitude of the mean profile. Sign: kernel e^(+2πi f n), rising phase = drift to LATER longitudes
    (the TPA XX convention: negative gradient = negative drift). With three features the beat residual
    φ1 − (φ2 + φ3) per block uses the raw (not aligned) phases (TPA XX Fig. 4 bottom).
    """
    function lrfs_phase_blocks(data, bin_st, bin_end, p3s; block=64, snr_min=7.0, frac_min=60 / 64, nboot=200, nsig=3.0)
        N, nb = size(data)
        X, good = _demeaned(data)
        on = collect(bin_st:bin_end); off = _far_offpulse(nb, bin_st, bin_end)
        # S/N impulsu: energia w oknie / (sigma poza impulsem * sqrt(liczba binow))
        σb = [std(@view data[i, off]) for i in 1:N]
        snr1 = [good[i] && σb[i] > 0 ? (sum(@view data[i, on]) - length(on) * mean(@view data[i, off])) / (σb[i] * sqrt(length(on))) : 0.0 for i in 1:N]
        starts = [s for s in 1:block:N-block+1 if count(>=(snr_min), @view snr1[s:s+block-1]) >= frac_min * block]
        nbk = length(starts)
        F = Array{ComplexF64}(undef, nbk, block ÷ 2 + 1, length(on)); Fo = Array{ComplexF64}(undef, nbk, block ÷ 2 + 1, length(off))
        for (b, s) in enumerate(starts)
            Bn = X[s:s+block-1, on]; Bn .-= mean(Bn, dims=1)
            Bf = X[s:s+block-1, off]; Bf .-= mean(Bf, dims=1)
            F[b, :, :] = conj.(rfft(Bn, 1)); Fo[b, :, :] = conj.(rfft(Bf, 1))      # conj: kernel e^(+2πi f n)
        end
        prof = vec(mean(Float64.(data[good, on]), dims=1))
        jref = argmax(prof)
        spec = nbk > 0 ? vec(sum(abs2.(F), dims=(1, 3)))[2:end] : Float64[]
        feats = []
        for p3 in p3s
            nbk == 0 && break
            k0 = round(Int, block / p3)
            ks = [k for k in k0-1:k0+1 if 1 <= k <= block ÷ 2]
            k = ks[argmax(spec[ks])]
            raw = fill(NaN, nbk, length(on)); amp = zeros(nbk, length(on))
            for b in 1:nbk
                σ = std(vcat(real.(Fo[b, k + 1, :]), imag.(Fo[b, k + 1, :])))
                c = F[b, k + 1, :]
                amp[b, :] = 2 .* abs.(c) ./ block
                ok = abs.(c) .> nsig * σ
                raw[b, ok] .= angle.(c[ok])
            end
            # wyrownanie: jedna stala na blok, iteracyjnie do sredniej wazonej amplituda
            off_b = zeros(nbk)
            T = zeros(length(on))
            for it in 1:5
                T = [let z = sum((amp[b, j] * cis(raw[b, j] - off_b[b]) for b in 1:nbk if isfinite(raw[b, j])); init=0.0im); angle(z) end
                     for j in eachindex(on)]
                for b in 1:nbk
                    z = sum((amp[b, j] * cis(raw[b, j] - T[j]) for j in eachindex(on) if isfinite(raw[b, j])); init=0.0im)
                    off_b[b] = angle(z)                      # przesuniecie bloku wzgledem sredniej sciezki
                end
            end
            al = raw .- off_b
            # srednia kolowa + bootstrap po blokach; odniesienie 0 w szczycie profilu
            cm(rows, j) = (v = [al[b, j] for b in rows if isfinite(al[b, j])]; isempty(v) ? NaN : angle(sum(cis.(v))))
            mean_ph = [cm(1:nbk, j) for j in eachindex(on)]
            ref = isfinite(mean_ph[jref]) ? mean_ph[jref] : 0.0
            boots = [[cm(rand(1:nbk, nbk), j) for j in eachindex(on)] for _ in 1:nboot]
            err = [let d = [angle(cis(bt[j] - mean_ph[j])) for bt in boots if isfinite(bt[j])]; length(d) > 2 ? rad2deg(std(d)) : NaN end
                   for j in eachindex(on)]
            # rozwiniecie sredniej od szczytu profilu
            mdeg = rad2deg.(angle.(cis.(mean_ph .- ref)))
            for dir in (1, -1)
                prev = 0.0; j = jref + dir
                while 1 <= j <= length(on)
                    if isfinite(mdeg[j])
                        while mdeg[j] - prev > 180; mdeg[j] -= 360; end
                        while mdeg[j] - prev < -180; mdeg[j] += 360; end
                        prev = mdeg[j]
                    end
                    j += dir
                end
            end
            aldeg = rad2deg.(angle.(cis.(al .- ref)))                  # do histogramu (-180..180)
            push!(feats, (p3=p3, k=k, f=k / block, raw=raw, aligned=aldeg, mean=mdeg, err=err,
                          amp=vec(mean(amp, dims=1))))
        end
        # reszta dudnien tylko gdy biny spelniaja k1 = k2 + k3 (inaczej faza surowa nie znosi sie miedzy blokami)
        beat = nothing; beat_ok = false
        if length(feats) == 3
            beat_ok = feats[1].k == feats[2].k + feats[3].k
            r = feats[1].raw .- (feats[2].raw .+ feats[3].raw)
            beat = beat_ok ? rad2deg.(angle.(cis.(r))) : nothing
        end
        return (feats=feats, beat=beat, beat_ok=beat_ok, nblocks=nbk, block=block, prof=prof, jref=jref, spec=spec)
    end

    """
        lrfs_tpa_plot(data, bin_st, bin_end, outfile; nfft=200, title_="")

    Power LRFS like TPA XX (Hsu+2026) Fig. 2a-b: continuous blocks of `nfft` pulses (they used 200), power
    averaged over blocks, the DC row (1/P3 = 0) set to zero, normalised so that the maximum = 1, logarithmic
    grey scale; profile on top. Longitude 0 = profile peak.
    """
    function lrfs_tpa_plot(data, bin_st, bin_end, outfile; nfft=200, title_="", pub=true)
        pub && _pub_style()
        N, nb = size(data)
        X, good = _demeaned(data)
        on = bin_st:bin_end
        nblk = N ÷ nfft
        nblk >= 1 || return nothing
        P = zeros(nfft ÷ 2 + 1, length(on))
        for b in 1:nblk
            B = X[(b-1)*nfft+1:b*nfft, on]; B .-= mean(B, dims=1)
            P .+= abs2.(rfft(B, 1))
        end
        P[1, :] .= 0.0
        P ./= maximum(P)
        prof = vec(mean(Float64.(data[good, on]), dims=1))
        lon = (collect(on) .- 1) .* (360 / nb)                  # dlugosc bezwzgledna w stopniach (jak w TPA XX)
        f = collect(0:nfft÷2) ./ nfft
        pub || (rc("font", size=7); rc("axes", linewidth=0.5))
        figure(figsize=(4.0, 5.5))
        subplots_adjust(left=0.16, bottom=0.08, right=0.97, top=0.9, hspace=0.0)
        suptitle(title_, fontsize=6.5)
        subplot2grid((4, 1), (0, 0))
        plot(lon, prof ./ maximum(prof), c="black"); xlim(lon[1], lon[end]); tick_params(labelbottom=false); ylabel("I (norm.)")
        subplot2grid((4, 1), (1, 0), rowspan=3)
        imshow(log10.(max.(P, 1e-4)), origin="lower", aspect="auto", cmap="Greys",
               extent=[lon[1], lon[end], 0, f[end]], vmin=-3, vmax=0)
        xlabel("pulse longitude (°)"); ylabel("1/P₃ (cpp)")
        savefig(outfile); println(outfile)
        close("all"); pub && _pub_style_off()
        return (nblocks=nblk,)
    end

    """
        tdfs_tpa_plot(data, bin_st, bin_end, outfile; nfft=200, title_="")

    Magnitude 2DFS like TPA XX (Hsu+2026) Fig. 3 for the whole on-pulse window: continuous blocks of
    `nfft` pulses, power averaged over blocks, magnitude = sqrt(power), DC row (1/P3 = 0) set to zero,
    normalised to max = 1, grey scale capped at 50%; horizontally FLIPPED w.r.t. the mathematical 2D FFT
    (their convention, see H25). Side panels: horizontally (left) and vertically (bottom) integrated power.
    1/P2 in cycles per period: column m of an n-bin window -> m·nbin/n cpp.
    """
    function tdfs_tpa_plot(data, bin_st, bin_end, outfile; nfft=200, title_="", pub=true)
        pub && _pub_style()
        N, nb = size(data)
        X, _ = _demeaned(data)
        on = bin_st:bin_end; n = length(on)
        nblk = N ÷ nfft
        nblk >= 1 || return nothing
        P = zeros(nfft, n)
        for b in 1:nblk
            B = X[(b-1)*nfft+1:b*nfft, on]; B .-= mean(B, dims=1)
            P .+= abs2.(fft(B))
        end
        P = fftshift(P, 2)[1:nfft÷2+1, :]                        # 1/P3 od 0 do 0.5, 1/P2 od ujemnych do dodatnich
        P[1, :] .= 0.0
        P = reverse(P, dims=2)                                    # konwencja TPA XX: odbicie w poziomie
        M = sqrt.(P); M ./= maximum(M)
        m = collect(-(n ÷ 2):(n - 1 - n ÷ 2))
        x = reverse(m) .* (nb / n)                                # 1/P2 (cpp) po odbiciu
        f = collect(0:nfft÷2) ./ nfft
        pub || (rc("font", size=7); rc("axes", linewidth=0.5))
        figure(figsize=(4.5, 4.5))
        subplots_adjust(left=0.16, bottom=0.12, right=0.97, top=0.9, wspace=0.0, hspace=0.0)
        suptitle(title_, fontsize=6.5)
        subplot2grid((4, 4), (0, 0), rowspan=3)
        plot(vec(sum(P, dims=2)), f, c="black"); ylim(0, f[end]); xticks([]); ylabel("1/P₃ (cpp)")
        gca().invert_xaxis()
        subplot2grid((4, 4), (0, 1), rowspan=3, colspan=3)
        imshow(M, origin="lower", aspect="auto", cmap="Greys", extent=[x[end], x[1], 0, f[end]], vmin=0, vmax=0.5)
        tick_params(labelleft=false, labelbottom=false)
        subplot2grid((4, 4), (3, 1), colspan=3)
        plot(x, vec(sum(P, dims=1)), c="black"); xlim(x[end], x[1]); yticks([]); xlabel("1/P₂ (cpp)")
        savefig(outfile); println(outfile)
        close("all"); pub && _pub_style_off()
        return (nblocks=nblk,)
    end

    "Plot like TPA XX Fig. 4: profile + modulation amplitudes, phase histograms per feature with block-averaged
    phases and bootstrap errors, beat residual (if 3 features). Longitude 0 = profile peak."
    function lrfs_phase_blocks_plot(res, lon, outfile; title_="", xlab="longitude (°), 0 = profile peak", pub=false)
        pub && _pub_style()
        nf = length(res.feats); nrow = 1 + nf + (isnothing(res.beat) ? 0 : 1)
        cols = ["royalblue", "green", "red"]; cmaps = ["Blues", "Greens", "Reds"]
        pub || (rc("font", size=7); rc("axes", linewidth=0.5))
        figure(figsize=(4.5, 1.6 * nrow + 0.6))
        subplots_adjust(left=0.15, bottom=0.06, right=0.97, top=0.93, hspace=0.0)
        suptitle(title_, fontsize=6.5)
        subplot(nrow, 1, 1)
        plot(lon, res.prof, c="black", lw=0.8, label="profile")
        for (i, f) in enumerate(res.feats)
            plot(lon, f.amp, c=cols[mod1(i, 3)], ls=["--", "-.", ":"][mod1(i, 3)], lw=0.8, label=@sprintf("P3 = %.2f", 1 / f.f))
        end
        legend(fontsize=5.5); xlim(lon[1], lon[end]); tick_params(labelbottom=false); yticks([])
        edges = range(-180, 180, length=37)
        for (i, f) in enumerate(res.feats)
            subplot(nrow, 1, 1 + i)
            H = zeros(length(edges) - 1, length(lon))
            for j in eachindex(lon), b in axes(f.aligned, 1)
                v = f.aligned[b, j]; isfinite(v) || continue
                H[clamp(searchsortedlast(edges, v), 1, length(edges) - 1), j] += 1
            end
            H ./= max(size(f.aligned, 1), 1)
            imshow(H, origin="lower", aspect="auto", cmap=cmaps[mod1(i, 3)], extent=[lon[1], lon[end], -180, 180])
            # srednia tylko w binach istotnych w >= polowie blokow (bez pojedynczych punktow szumu poza impulsem)
            ok = isfinite.(f.mean) .& [mean(isfinite, @view f.aligned[:, j]) >= 0.5 for j in eachindex(lon)]
            m = mod.(f.mean .+ 180, 360) .- 180
            errorbar(lon[ok], m[ok], yerr=f.err[ok], fmt=".", ms=2, c="black", elinewidth=0.4, capsize=0)
            ylim(-180, 180); yticks([-180, 0, 180]); xlim(lon[1], lon[end]); tick_params(labelbottom=false)
            ylabel(@sprintf("phase (°)\nP3 = %.2f", 1 / f.f))
        end
        if !isnothing(res.beat)
            subplot(nrow, 1, nrow)
            H = zeros(length(edges) - 1, length(lon))
            for j in eachindex(lon), b in axes(res.beat, 1)
                v = res.beat[b, j]; isfinite(v) || continue
                H[clamp(searchsortedlast(edges, v), 1, length(edges) - 1), j] += 1
            end
            imshow(H ./ max(size(res.beat, 1), 1), origin="lower", aspect="auto", cmap=pub ? PyPlot.matplotlib.colors.LinearSegmentedColormap.from_list("brown", ["white", "#7b3f00"]) : "Oranges",
                   extent=[lon[1], lon[end], -180, 180])   # TPA XX: brazowa mapa, biale = 0
            axhline(0, c="grey", ls="--", lw=0.5); ylim(-180, 180); yticks([-180, 0, 180]); xlim(lon[1], lon[end])
            ylabel("φ1 − (φ2 + φ3) (°)")
        end
        xlabel(xlab)
        savefig(outfile); println(outfile)
        close("all"); pub && _pub_style_off()
    end


end  # module Plot
