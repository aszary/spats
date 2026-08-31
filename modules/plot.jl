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


    function lrfs(data, outdir, params; cmap="viridis", darkness=0.5, name_mod="PSR_NAME",
              show_=false)
        
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
        profile = mean(real.(lrfs_complex), dims=1)[1, :]

        # === Power spectrum (left panel) ===
        spectrum = sum(amp, dims=2)[:, 1]

        # === Frequency axis ===
        freq = collect(0:(num-1)) ./ num
        peak = argmax(spectrum[2:end]) + 1  # skip DC component at freq=0
        println("\tPeak freq: $(freq[peak]), P3 = $(1/freq[peak])")

        # === Phase at peak frequency ===
        phase_ = rad2deg.(angle.(lrfs_complex[peak, :]))

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

        # --- TOP PANEL: phase ---
        ax_phase = subplot2grid((5, 3), (0, 1), colspan=2)
        #ax_phase = axes([0.17, 0.94, 0.82, 0.05])
        scatter(longitude, phase_, s=2, c="grey")
        xticks([])
        yticks([-180, -90, 0, 90, 180])
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

    function analyse_p3folds3(low, high, p, n_comp)
        pulses, bins = size(low)

        # Collected offsets per pulse: Dict(component => (longitudes, offsets, errors))
        offset_data = Dict{Int, NamedTuple{(:lon, :off, :err), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}}()

        for i in 1:pulses
            # Extract data in the range bin_st:bin_end
            x_data = p["bin_st"]:p["bin_end"]
            y_low = low[i, p["bin_st"]:p["bin_end"]]
            y_high = high[i, p["bin_st"]:p["bin_end"]]

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
            xlim(p["bin_st"], p["bin_end"])

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
            xlim(p["bin_st"], p["bin_end"])

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


    "added s to skip point"
    function analyse_p3folds4(low, high, p, n_comp)
        pulses, bins = size(low)

        # Collected offsets per pulse: Dict(component => (longitudes, offsets, errors))
        offset_data = Dict{Int, NamedTuple{(:lon, :off, :err), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}}()

        for i in 1:pulses
            # Extract data in the range bin_st:bin_end
            x_data = p["bin_st"]:p["bin_end"]
            y_low = low[i, p["bin_st"]:p["bin_end"]]
            y_high = high[i, p["bin_st"]:p["bin_end"]]

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
            xlim(p["bin_st"], p["bin_end"])

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
            xlim(p["bin_st"], p["bin_end"])

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

        # Weighted mean offset per component
        if !isempty(offset_data)
            println("\n=== Weighted mean offsets ===")
            for comp in sort(collect(keys(offset_data)))
                d = offset_data[comp]
                if isempty(d.err) || all(d.err .== 0.0)
                    continue
                end
                w       = 1.0 ./ (d.err .^ 2)
                n       = length(d.off)
                mu      = sum(w .* d.off) / sum(w)
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
    function analyse_average_offset(nl, nh, p, n_comp; npulse=150, n_pulses=nothing)
        n_profiles = size(nl, 1)

        offset_data = Dict{Int, NamedTuple{(:idx, :off, :err),
                                           Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}}()

        low_colors  = ["#2196F3", "#0D47A1", "#64B5F6", "#1565C0"]
        high_colors = ["#FF6F00", "#E65100", "#FFCA28", "#F57F17"]

        for i in 1:n_profiles
            x_data = p["bin_st"]:p["bin_end"]
            y_low  = nl[i, p["bin_st"]:p["bin_end"]]
            y_high = nh[i, p["bin_st"]:p["bin_end"]]

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
                        offset_data[o.component] = (idx=Float64[], off=Float64[], err=Float64[])
                    end
                    push!(offset_data[o.component].idx, pulse_center)
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
            xlim(p["bin_st"], p["bin_end"])

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
            xlim(p["bin_st"], p["bin_end"])

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

        # Weighted mean offset per component
        if !isempty(offset_data)
            println("\n=== Weighted mean offsets ===")
            for comp in sort(collect(keys(offset_data)))
                d = offset_data[comp]
                if isempty(d.err) || all(d.err .== 0.0)
                    continue
                end
                w         = 1.0 ./ (d.err .^ 2)
                n         = length(d.off)
                mu        = sum(w .* d.off) / sum(w)
                sigma_int = 1.0 / sqrt(sum(w))
                chi2      = sum(w .* (d.off .- mu) .^ 2)
                dof       = n - 1
                chi2_red  = dof > 0 ? chi2 / dof : NaN
                sigma_ext = sigma_int * sqrt(max(1.0, chi2_red))
                println(@sprintf("G%d: offset = %+.4f° ± %.4f°  (n=%d, χ²/dof = %.2f, σ_ext = %.4f°)",
                    comp, mu, sigma_int, n, chi2_red, sigma_ext))
            end
            println()
        end

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
        name = ""

        for line in eachline(filename)
            if startswith(line, "@-")
                _psrcat_record!(names, periods, pdots, name, rec)
                empty!(rec)
                name = ""
                continue
            end
            (isempty(strip(line)) || startswith(line, "#")) && continue
            parts = split(line)
            length(parts) < 2 && continue
            key, val = parts[1], parts[2]
            if key == "PSRJ" || (key == "PSRB" && name == "")
                name == "" && (name = val)
            elseif key in ("P0", "P1", "F0", "F1") && !haskey(rec, key)
                v = tryparse(Float64, val)
                isnothing(v) || (rec[key] = v)
            end
        end
        _psrcat_record!(names, periods, pdots, name, rec)  # last record (no trailing @-)

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


    """ Label a power of ten as 10^n, anything else as a plain number. """
    function _pow10_label(v)
        e = log10(v)
        isapprox(e, round(e), atol=1e-9) && return "\$10^{$(Int(round(e)))}\$"
        return @sprintf("%.3g", v)
    end


    """
    P-Pdot diagram of the ATNF catalogue.

    Pulsars with Pdot <= 0 (mostly globular-cluster pulsars, where acceleration
    in the cluster potential dominates the observed spin-down) are skipped —
    they cannot be shown on a logarithmic axis.
    """
    function ppdot(outdir; catalogue=normpath(joinpath(@__DIR__, "..", "input", "psrcat.db")),
                   name_mod="psrcat", show_=true,
                   b_lines=[1e10, 1e12, 1e14],     # surface magnetic field [G]
                   age_lines=[1e3, 1e6, 1e9],      # characteristic age [yr]
                   edot_lines=[1e30, 1e33, 1e36],  # spin-down luminosity [erg/s]
                   death_bp2=0.17e12,              # death line B/P^2 [G s^-2]
                   inertia=1e45,                   # moment of inertia [g cm^2]
                   plims=(1e-3, 2e2), pdotlims=(1e-22, 1e-8))

        names, periods, pdots = read_psrcat(catalogue)
        println("psrcat: $(length(names)) pulsars with P and Pdot ($catalogue)")

        keep = pdots .> 0
        println("psrcat: $(count(.!keep)) pulsars skipped (Pdot <= 0)")
        p = periods[keep]
        pd = pdots[keep]

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figsize = (3.54331, 3.14961)  # 9cm x 8cm
        figure(figsize=figsize, frameon=true)
        subplots_adjust(left=0.19, bottom=0.14, right=0.97, top=0.97)

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
            ax.plot(xs, ys, ls=ls, lw=0.6, c=color, alpha=0.8, zorder=1)
            isnothing(txt) && return
            xt = 10^(log10(xs[1]) + frac * log10(xs[end] / xs[1]))
            ax.text(xt, A * xt^m, txt, fontsize=5, color=color,
                    ha="center", va=va, rotation_mode="anchor", zorder=4,
                    rotation=_ppdot_angle(ax, m, plims, pdotlims, figsize))
        end

        # the label fractions keep each family near a different edge of the box;
        # the B labels hang below their lines to stay clear of the upper frame
        for b in b_lines
            line!((b / 3.2e19)^2, -1, _pow10_label(b), 0.10, c_b, "--"; va="top")
        end
        for tau in age_lines
            line!(1 / (2 * tau * yr_s), 1, _pow10_label(tau), 0.85, c_age, "-.")
        end
        for edot in edot_lines
            line!(edot / (4 * pi^2 * inertia), 3, _pow10_label(edot), 0.85, c_edot, ":")
        end
        line!((death_bp2 / 3.2e19)^2, 3, nothing, 0.5, c_death, "-")

        plot(p, pd, ".", ms=1.3, c="black", mec="none", zorder=3)

        xlabel("\$P\$ (s)")
        ylabel("\$\\dot{P}\$ (s s\$^{-1}\$)")
        minorticks_on()

        L2D = PyPlot.matplotlib.lines.Line2D
        legend(handles=[L2D([], [], c=c_b, ls="--", lw=0.6, label="\$B\$ (G)"),
                        L2D([], [], c=c_age, ls="-.", lw=0.6, label="\$\\tau_c\$ (yr)"),
                        L2D([], [], c=c_edot, ls=":", lw=0.6, label="\$\\dot{E}\$ (erg s\$^{-1}\$)"),
                        L2D([], [], c=c_death, ls="-", lw=0.6, label="death line")],
               fontsize=5, loc="lower right", framealpha=0.9, borderpad=0.4,
               handlelength=2.5, labelspacing=0.35)

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


end  # module Plot
