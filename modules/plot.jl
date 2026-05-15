module Plot
    using LsqFit
    using FFTW
    using JLD2
    using Statistics
    using StatsBase
    if haskey(ENV, "DISPLAY") && !isempty(get(ENV, "DISPLAY", ""))
        ENV["MPLBACKEND"] = "TkAgg"
    else
        ENV["MPLBACKEND"] = "Agg"
    end
    using PyPlot
    using PyCall
    using Printf
    @pyimport matplotlib.patches as patch
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
    include("profile_metrics.jl")
    using .ProfileMetrics


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
            close()
            
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
            close()
            
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
            close()
            
            if lowercase(strip(user_input)) == "q"
                println("Exiting analysis.")
                break
            end
        end
    end    



    # Snap a raw step to the nearest 1/2/5·10^n
    function _snap_step(raw)
        raw <= 0 && return 0.0
        pw   = 10.0 ^ floor(log10(raw))
        mant = raw / pw
        mant < 1.5 && return 1.0 * pw
        mant < 3.0 && return 2.0 * pw
        mant < 7.0 && return 5.0 * pw
        return 10.0 * pw
    end

    # Auto height contour levels spanning all supplied maps, with nice 1/2/5·10^n steps
    function _auto_height_levels(maps...)
        finite_h = filter(isfinite, vcat((vec(m) for m in maps)...))
        isempty(finite_h) && return Float64[]
        h_lo = minimum(finite_h)
        h_hi = maximum(finite_h)
        dh = _snap_step((h_hi - h_lo) / 8)
        dh > 0 || return Float64[]
        return collect(ceil(h_lo / dh) * dh : dh : floor(h_hi / dh) * dh)
    end

    """
    rvm(data4, outdir, name_mod; bin_st, bin_end,
        snr_threshold=3.5, linpol_threshold=0.8, period=nothing,
        alpha=nothing, beta=nothing, show_=false)

    4-panel RVM + emission-height figure (Johnston et al. 2023/2024 style):
      panel 1 : averaged Stokes I profile
      panel 2 : PPA scatter (grey=all, blue=filtered) + red RVM curve
      panel 3 : residuals of the RVM fit
      panel 4 : emission height per bin (Blaskiewicz + dipole formulas)
                shown only when period is provided

    data4    : 3D array (pulses × bins × 4), I=1 Q=2 U=3 V=4
    outdir   : output directory
    name_mod : filename prefix
    bin_st, bin_end : on-pulse region (1-indexed)
    period   : pulsar period in seconds  (needed for height panel)
    alpha, beta : geometry angles in degrees (needed for dipole height)
    """
    function rvm(data4, outdir, name_mod; bin_st, bin_end,
                 snr_threshold=5.0, period=nothing, show_=false, n_freq=1)

        n_full = size(data4, 2)
        n_on   = bin_end - bin_st + 1
        lon_full = collect(range(-180.0, 180.0, length=n_full))
        dl_on    = 360.0 * n_on / n_full
        lon_on   = collect(range(-dl_on/2.0, dl_on/2.0, length=n_on))

        # --- Two-pass RVM fit for one frequency band --------------------------
        # Pass 1: raw errors → chi2_red.  Pass 2: errors × max(1,√χ²ᵣ) so the
        # displayed chi2_red_min ≈ 1 (Johnston+ 2023 convention).
        function _fit_band(data, label)
            println("[rvm-$label]")
            lon_f, pa, pa_err, mask = Tools.filter_ppa(data, bin_st, bin_end;
                                                        snr_threshold=snr_threshold)
            if sum(mask) < 5
                @warn "[rvm] $label: only $(sum(mask)) valid PA bins — skipping"
                return nothing, nothing, nothing, nothing, lon_f, pa, pa_err, mask
            end
            # pass 1
            al, be, _, p0 = Tools.chi2_grid(lon_f, pa, pa_err, mask)
            res1 = Tools.fit_rvm(lon_f, pa, pa_err, mask; p0=p0)
            scale = max(1.0, sqrt(res1.chi2_red))
            println("  χ²ᵣ(pass-1)=$(round(res1.chi2_red,digits=2))  σ×$(round(scale,digits=2))")
            pa_err2 = [mask[i] ? pa_err[i]*scale : pa_err[i] for i in eachindex(pa_err)]
            # pass 2
            al2, be2, chi2_map, p1 = Tools.chi2_grid(lon_f, pa, pa_err2, mask)
            result  = Tools.fit_rvm(lon_f, pa, pa_err2, mask; p0=p1)
            return result, al2, be2, chi2_map, lon_f, pa, pa_err2, mask
        end

        # Helper: NaN-safe RVM curve over lon_on
        function _rvm_curve(res)
            pa = Tools.rvm_model(lon_on, [res.PA0, res.alpha, res.zeta, res.phi0])
            pa = mod.(pa .+ 90.0, 180.0) .- 90.0
            for i in Iterators.drop(eachindex(pa), 1)
                abs(pa[i] - pa[i-1]) > 90.0 && (pa[i-1] = NaN)
            end
            return pa
        end

        # Helper: averaged Stokes I / L / V over all bins
        function _stokes(d)
            I = vec(mean(d[:, :, 1], dims=1))
            L = sqrt.(vec(mean(d[:, :, 2], dims=1)).^2 .+
                      vec(mean(d[:, :, 3], dims=1)).^2)
            V = vec(mean(d[:, :, 4], dims=1))
            return I, L, V
        end

        # Helper: chi² display (subtract min, clip above dmax → NaN)
        function _chi2_disp(cm, dmax)
            fin = filter(isfinite, vec(cm))
            isempty(fin) && return fill(NaN, size(cm)...)
            cmin  = minimum(fin)
            delta = map(v -> isnan(v) ? NaN : max(v - cmin, 0.0), cm)
            return map(v -> isnan(v) ? NaN : (v > dmax ? NaN : v), delta)
        end

        # --- Fit per frequency band -------------------------------------------
        if n_freq == 2
            println("[rvm] Two-frequency mode: separate fits for LOW and HIGH")
            n_each    = div(size(data4, 1), 2)
            data_low  = data4[1:n_each, :, :]
            data_high = data4[n_each+1:end, :, :]

            res_l, al_l, be_l, cm_l, lon_l, pa_l, pe_l, mk_l = _fit_band(data_low,  "low")
            res_h, al_h, be_h, cm_h, lon_h, pa_h, pe_h, mk_h = _fit_band(data_high, "high")

            I_low,  L_low,  V_low  = _stokes(data_low)
            I_high, L_high, V_high = _stokes(data_high)
            I_max = max(maximum(abs.(I_low)), maximum(abs.(I_high)))
            I_ref = I_low

            # Primary result for h_blask, W10 (prefer low-freq; fall back to high)
            result = isnothing(res_l) ? res_h : res_l
            alphas = isnothing(al_l)  ? al_h  : al_l
            betas  = isnothing(be_l)  ? be_h  : be_l
            cm_primary = isnothing(cm_l) ? cm_h : cm_l

            crv_l = isnothing(res_l) ? nothing : _rvm_curve(res_l)
            crv_h = isnothing(res_h) ? nothing : _rvm_curve(res_h)
        else
            println("[rvm] Single-frequency mode")
            result, alphas, betas, cm_primary, lon_l, pa_l, pe_l, mk_l =
                _fit_band(data4, "1freq")

            I_low,  L_low,  V_low  = _stokes(data4)
            I_high = L_high = V_high = nothing
            I_max  = maximum(abs.(I_low))
            I_ref  = I_low

            cm_l = cm_primary; cm_h = nothing
            al_l = alphas;     be_l = betas
            al_h = be_h = nothing
            lon_h = pa_h = pe_h = mk_h = nothing
            res_l = result;    res_h = nothing
            crv_l = isnothing(result) ? nothing : _rvm_curve(result)
            crv_h = nothing
        end

        if isnothing(result)
            @warn "[rvm] No valid fit for $name_mod, skipping plot"
            return (h_blask=NaN, phi_center=NaN, W10=NaN, alpha=NaN, zeta=NaN,
                    phi0=NaN, chi2_red=NaN, alphas=Float64[], betas=Float64[],
                    chi2_map=Matrix{Float64}(undef,0,0))
        end

        phi_c   = Tools.pulse_center_deg(lon_full, I_ref)
        W10     = Tools.profile_width_deg(lon_full, I_ref)
        h_blask = period !== nothing ?
                  Tools.emission_height_blaskiewicz(result.phi0, phi_c, period) : NaN

        function _hcont(W, al, be)
            (period === nothing || W <= 0.0) && return nothing
            return Float64[Tools.height_from_rho_rankin(
                               Tools.rho_from_width(a, b, W), period)
                           for a in al, b in be]
        end
        hc_l = _hcont(W10, alphas, betas)
        hc_h = (n_freq == 2 && !isnothing(I_high)) ?
               _hcont(Tools.profile_width_deg(lon_full, I_high), al_h, be_h) : nothing

        # ---- shared Δχ² display scale ----------------------------------------
        all_c2 = filter(isfinite, vcat(
            vec(cm_l),
            isnothing(cm_h) ? Float64[] : vec(cm_h)))
        c2_range = isempty(all_c2) ? 1.0 : maximum(all_c2) - minimum(all_c2)
        dmax = max(_snap_step(0.02 * c2_range), 1.0)
        println("[rvm] shared Δχ² dmax=$(round(dmax,digits=2))")

        disp_l = _chi2_disp(cm_l, dmax)
        disp_h = isnothing(cm_h) ? nothing : _chi2_disp(cm_h, dmax)

        # ---- figure layout ---------------------------------------------------
        rc("font", family="sans-serif", size=9.)
        rc("axes", linewidth=0.7)
        rc("lines", linewidth=1.2)
        rc("xtick", direction="in", top=true)
        rc("ytick", direction="in", right=true)

        lon_margin = 0.30 * dl_on
        x_lo = lon_on[1]  - lon_margin
        x_hi = lon_on[end] + lon_margin

        ext_l = !isnothing(alphas) ? [alphas[1], alphas[end], betas[1], betas[end]] : [0,180,-20,20]
        ext_h = (n_freq == 2 && !isnothing(al_h)) ?
                [al_h[1], al_h[end], be_h[1], be_h[end]] : ext_l

        if n_freq == 2
            fig = figure(figsize=(12.5, 6.5))
            ax_pa  = fig.add_axes([0.06, 0.65, 0.40, 0.28])
            ax_st  = fig.add_axes([0.06, 0.09, 0.40, 0.51])
            ax_ml  = fig.add_axes([0.53, 0.09, 0.21, 0.74])   # low  chi² map
            ax_mh  = fig.add_axes([0.76, 0.09, 0.21, 0.74])   # high chi² map
            ax_cb  = fig.add_axes([0.53, 0.86, 0.44, 0.042])  # shared colorbar
        else
            fig = figure(figsize=(9.5, 6.5))
            ax_pa  = fig.add_axes([0.09, 0.65, 0.45, 0.28])
            ax_st  = fig.add_axes([0.09, 0.09, 0.45, 0.51])
            ax_ml  = fig.add_axes([0.62, 0.09, 0.36, 0.74])
            ax_mh  = nothing
            ax_cb  = fig.add_axes([0.62, 0.86, 0.36, 0.042])
        end

        # ---- PA panel --------------------------------------------------------
        # Low-freq data points (blue circles)
        pa_l_ok = .!isnan.(pa_l)
        ax_pa.errorbar(lon_l[pa_l_ok], pa_l[pa_l_ok], yerr=pe_l[pa_l_ok],
                       fmt="o", ms=3.5, color="tab:blue", alpha=0.88,
                       elinewidth=0.6, capsize=1.2, zorder=3,
                       label=n_freq==2 ? "low freq" : "data")
        # High-freq data points (orange squares)
        if n_freq == 2 && !isnothing(pa_h)
            pa_h_ok = .!isnan.(pa_h)
            ax_pa.errorbar(lon_h[pa_h_ok], pa_h[pa_h_ok], yerr=pe_h[pa_h_ok],
                           fmt="s", ms=3.0, color="tab:orange", alpha=0.85,
                           elinewidth=0.6, capsize=1.2, zorder=3, label="high freq")
        end
        # RVM curves
        if !isnothing(crv_l)
            ax_pa.plot(lon_on, crv_l, color="tab:blue",
                       lw=2.0, zorder=4, label=n_freq==2 ? "RVM low" : "RVM")
        end
        if !isnothing(crv_h)
            ax_pa.plot(lon_on, crv_h, color="tab:orange",
                       lw=1.8, ls="--", zorder=4, label="RVM high")
        end
        # φ₀ vertical lines
        if !isnothing(res_l)
            ax_pa.axvline(res_l.phi0, color="tab:blue",   lw=1.0, ls="--", alpha=0.80)
        end
        if !isnothing(res_h)
            ax_pa.axvline(res_h.phi0, color="tab:orange", lw=1.0, ls=":",  alpha=0.80)
        end

        # Annotation box
        if n_freq == 2
            txt_l = isnothing(res_l) ? "LOW:  —" :
                @sprintf("LOW:  α=%.1f° ζ=%.1f° φ₀=%.1f° χ²ᵣ=%.2f",
                         res_l.alpha, res_l.zeta, res_l.phi0, res_l.chi2_red)
            txt_h = isnothing(res_h) ? "HIGH: —" :
                @sprintf("HIGH: α=%.1f° ζ=%.1f° φ₀=%.1f° χ²ᵣ=%.2f",
                         res_h.alpha, res_h.zeta, res_h.phi0, res_h.chi2_red)
            txt_str = txt_l * "\n" * txt_h
        else
            txt_str = @sprintf("α=%.1f°  ζ=%.1f°  φ₀=%.1f°\nχ²ᵣ=%.2f  rms=%.1f°",
                               result.alpha, result.zeta, result.phi0,
                               result.chi2_red, result.rms_deg)
        end
        ax_pa.text(0.98, 0.97, txt_str,
                   transform=ax_pa."transAxes", fontsize=7.0,
                   va="top", ha="right", family="monospace",
                   bbox=Dict("boxstyle"=>"round,pad=0.3", "fc"=>"#fff8e7",
                             "ec"=>"#ccaa55", "alpha"=>0.88))
        ax_pa.set_xlim(x_lo, x_hi)
        ax_pa.set_ylim(-95, 95)
        ax_pa.set_ylabel("PA [deg]", fontsize=9)
        ax_pa.set_yticks([-90, -45, 0, 45, 90])
        ax_pa.tick_params(labelbottom=false)
        ax_pa.minorticks_on()
        ax_pa.set_title(name_mod, fontsize=10, fontweight="bold", pad=4)
        if n_freq == 2
            ax_pa.legend(fontsize=6.5, loc="upper left", framealpha=0.7,
                         ncol=2, handlelength=1.2, columnspacing=0.6)
        end

        # ---- Stokes profile panel --------------------------------------------
        ax_st.plot(lon_full, I_low ./ I_max, color="black",   lw=1.3, label="I",  zorder=3)
        ax_st.plot(lon_full, L_low ./ I_max, color="#cc3333", lw=1.3, label="L",  zorder=3)
        ax_st.plot(lon_full, V_low ./ I_max, color="#2255bb", lw=1.3, label="V",  zorder=3)
        if n_freq == 2 && !isnothing(I_high)
            ax_st.plot(lon_full, I_high ./ I_max, color="black",   lw=1.0, ls="--", alpha=0.6, zorder=2)
            ax_st.plot(lon_full, L_high ./ I_max, color="#cc3333", lw=1.0, ls="--", alpha=0.6, zorder=2)
            ax_st.plot(lon_full, V_high ./ I_max, color="#2255bb", lw=1.0, ls="--", alpha=0.6, zorder=2)
        end
        ax_st.axhline(0.0, color="grey", lw=0.5, ls=":", zorder=1)
        # φ₀ and pulse-centre reference lines
        if !isnothing(res_l)
            ax_st.axvline(res_l.phi0, color="tab:blue",   lw=1.0, ls="--", alpha=0.8,
                          label=@sprintf("φ₀=%.1f°", res_l.phi0))
        end
        if !isnothing(res_h)
            ax_st.axvline(res_h.phi0, color="tab:orange", lw=1.0, ls=":",  alpha=0.8,
                          label=@sprintf("φ₀ₕ=%.1f°", res_h.phi0))
        end
        ax_st.axvline(phi_c, color="#888888", lw=0.8, ls=":", alpha=0.65,
                      label=@sprintf("ctr=%.1f°", phi_c))
        ax_st.set_xlim(x_lo, x_hi)
        ax_st.set_ylim(-0.45, 1.18)
        ax_st.set_ylabel("Normalised flux", fontsize=9)
        ax_st.set_xlabel("Pulse longitude [deg]", fontsize=9)
        ax_st.legend(fontsize=6.5, loc="upper right", framealpha=0.7,
                     ncol=3, handlelength=1.2, columnspacing=0.6)
        ax_st.minorticks_on()
        if !isnan(h_blask)
            ax_st.text(0.02, 0.04,
                       @sprintf("h_Blask = %.0f km\nΔφ = %.1f°", h_blask, result.phi0 - phi_c),
                       transform=ax_st."transAxes", fontsize=8, va="bottom", ha="left",
                       family="monospace",
                       bbox=Dict("boxstyle"=>"round,pad=0.35", "fc"=>"#f0f8ff",
                                 "ec"=>"#5599cc", "alpha"=>0.88))
        end

        # ---- chi² map helper -------------------------------------------------
        function _draw_map(ax, disp, al, be, ext, hc, res, title_str; cb_ref=false)
            im = ax.imshow(disp', origin="lower", aspect="auto", extent=ext,
                           cmap="viridis", vmin=0.0, vmax=dmax, interpolation="nearest")
            # 1σ / 2σ contours on Δχ²
            fin = filter(isfinite, vec(disp))
            if !isempty(fin)
                try
                    ax.contour(al, be,
                               map(v -> isnan(v) ? 0.0 : v, disp)',
                               levels=[2.30, 6.17],
                               colors=["white", "#88ccff"],
                               linewidths=[0.9, 0.7],
                               linestyles=["solid", "dashed"])
                catch; end
            end
            # Best-fit star
            if !isnothing(res)
                ax.plot([res.alpha], [res.zeta - res.alpha],
                        marker="*", ms=8, color="white", mec="#ffdd00", mew=0.6, zorder=10)
            end
            # Height contours
            if !isnothing(hc)
                try
                    h_levs = _auto_height_levels(hc)
                    if !isempty(h_levs)
                        ax.clabel(ax.contour(al, be, hc', levels=h_levs,
                                             colors="orchid", linewidths=0.8),
                                  fmt="%g km", fontsize=5.5, inline=true)
                    end
                catch; end
            end
            ax.set_xlabel("α [deg]", fontsize=8)
            ax.set_title(title_str, fontsize=8, pad=3)
            ax.minorticks_on()
            return im
        end

        im_ref = _draw_map(ax_ml, disp_l, alphas, betas, ext_l, hc_l, res_l,
                           n_freq==2 ? "low freq  α–β" : "α–β confidence map")
        ax_ml.set_ylabel("β [deg]", fontsize=8)
        ax_ml.yaxis.set_label_position(n_freq==2 ? "left" : "right")
        ax_ml.yaxis.set_ticks_position(n_freq==2 ? "left" : "right")

        if n_freq == 2 && !isnothing(ax_mh) && !isnothing(disp_h)
            _draw_map(ax_mh, disp_h, al_h, be_h, ext_h, hc_h, res_h, "high freq  α–β")
            ax_mh.set_ylabel("β [deg]", fontsize=8)
            ax_mh.yaxis.set_label_position("right")
            ax_mh.yaxis.set_ticks_position("right")
            ax_mh.tick_params(labelleft=false)
        end

        cb = fig.colorbar(im_ref, cax=ax_cb, orientation="horizontal")
        cb.set_label(raw"$\Delta\chi^2_\mathrm{r}$", labelpad=2, fontsize=8)
        ax_cb.xaxis.set_label_position("top")
        ax_cb.xaxis.set_ticks_position("top")

        fig.savefig("$outdir/$(name_mod)_rvm.pdf", bbox_inches="tight", dpi=150)
        fig.savefig("$outdir/$(name_mod)_rvm.png", bbox_inches="tight", dpi=150)
        println("$outdir/$(name_mod)_rvm.pdf")

        if show_
            PyPlot.show()
            println("Press Enter to close.")
            readline(stdin; keep=false)
        end
        PyPlot.close()

        return merge(result, (h_blask=h_blask, phi_center=phi_c, W10=W10,
                              alphas=alphas, betas=betas, chi2_map=cm_primary))
    end


    # -------------------------------------------------------------------------
    # Two-panel geometry plot: low-freq and high-freq chi² maps side by side
    # -------------------------------------------------------------------------
    function geometry(chi2_map_l, chi2_map_h, alphas_deg, betas_deg, outdir;
                      show_=false, name_mod="PSR_NAME",
                      delta_chi2_max=nothing,
                      P_sec=nothing, W_deg_l=nothing, W_deg_h=nothing,
                      cmap="viridis")

        chi2_min_l = minimum(chi2_map_l[isfinite.(chi2_map_l)])
        chi2_min_h = minimum(chi2_map_h[isfinite.(chi2_map_h)])
        chi2_max_l = maximum(chi2_map_l[isfinite.(chi2_map_l)])
        chi2_max_h = maximum(chi2_map_h[isfinite.(chi2_map_h)])
        println("chi2_min: low=$(round(chi2_min_l, digits=2))  high=$(round(chi2_min_h, digits=2))")

        dchi2_l = chi2_map_l .- chi2_min_l
        dchi2_h = chi2_map_h .- chi2_min_h
        range_l = chi2_max_l - chi2_min_l
        range_h = chi2_max_h - chi2_min_h

        dmax = isnothing(delta_chi2_max) ? 0.02 * max(range_l, range_h) : delta_chi2_max
        println("delta_chi2_max = $(round(dmax, digits=2))")

        # Height map: cos ρ = cos α cos(α+β) + sin α sin(α+β) cos(W10/2)
        #             h = 2 c P ρ² / (9π)  [km]
        function _height_map(W_deg)
            hmap = fill(NaN, length(alphas_deg), length(betas_deg))
            (isnothing(P_sec) || isnothing(W_deg)) && return hmap
            c_km_s = 2.99792458e5
            cosW2  = cos(deg2rad(W_deg) / 2)
            for i in eachindex(alphas_deg)
                a = deg2rad(alphas_deg[i])
                for j in eachindex(betas_deg)
                    b  = deg2rad(betas_deg[j])
                    cr = clamp(cos(a)*cos(a+b) + sin(a)*sin(a+b)*cosW2, -1.0, 1.0)
                    hmap[i, j] = 2 * c_km_s * P_sec * acos(cr)^2 / (9π)
                end
            end
            return hmap
        end

        hmap_l = _height_map(W_deg_l)
        hmap_h = _height_map(W_deg_h)
        h_levs = _auto_height_levels(hmap_l, hmap_h)
        println("height_contours [km] = $(round.(h_levs, digits=1))")

        # Mask cells above dmax as NaN → renders white
        dchi2_l_plot = Float64.(dchi2_l); dchi2_l_plot[dchi2_l .> dmax] .= NaN
        dchi2_h_plot = Float64.(dchi2_h); dchi2_h_plot[dchi2_h .> dmax] .= NaN

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        fig = figure(figsize=(7.0, 4.2))
        subplots_adjust(left=0.09, bottom=0.11, right=0.98, top=0.78, wspace=0.32)

        extent = [alphas_deg[1], alphas_deg[end], betas_deg[1], betas_deg[end]]

        ax1 = subplot(1, 2, 1)
        im1 = imshow(dchi2_l_plot', origin="lower", extent=extent, aspect="auto",
                     cmap=cmap, vmin=0., vmax=dmax, interpolation="nearest")
        if !isempty(h_levs)
            ax1.clabel(ax1.contour(alphas_deg, betas_deg, hmap_l',
                                   levels=h_levs, colors="purple", linewidths=0.7),
                       fmt="%g", fontsize=6)
        end
        # Keep confidence contours (Czarek6 addition)
        try
            ax1.contour(alphas_deg, betas_deg, dchi2_l',
                        levels=[2.30, 6.17], colors=["cyan", "royalblue"],
                        linewidths=[1.0, 0.8], linestyles=["solid", "dashed"])
        catch; end
        ax1.set_xlabel("alpha [deg]")
        ax1.set_ylabel("beta [deg]")
        ax1.set_title("low frequency", pad=3)
        ax1.minorticks_on()

        ax2 = subplot(1, 2, 2)
        imshow(dchi2_h_plot', origin="lower", extent=extent, aspect="auto",
               cmap=cmap, vmin=0., vmax=dmax, interpolation="nearest")
        if !isempty(h_levs)
            ax2.clabel(ax2.contour(alphas_deg, betas_deg, hmap_h',
                                   levels=h_levs, colors="purple", linewidths=0.7),
                       fmt="%g", fontsize=6)
        end
        try
            ax2.contour(alphas_deg, betas_deg, dchi2_h',
                        levels=[2.30, 6.17], colors=["cyan", "royalblue"],
                        linewidths=[1.0, 0.8], linestyles=["solid", "dashed"])
        catch; end
        ax2.set_xlabel("alpha [deg]")
        ax2.set_ylabel("beta [deg]")
        ax2.set_title("high frequency", pad=3)
        ax2.minorticks_on()

        # Shared horizontal colorbar at top (master style)
        cbar_ax = fig.add_axes([0.09, 0.86, 0.89, 0.05])
        cbar = fig.colorbar(im1, cax=cbar_ax, orientation="horizontal")
        cbar.set_label(raw"$\Delta\chi^2$", labelpad=2)
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
        PyPlot.close()
    end


    # -------------------------------------------------------------------------
    # Two-frequency PA + Stokes profile plot (master style)
    # -------------------------------------------------------------------------
    function position_angle(lon_l, pa_l, pa_err_l, I_l, Lin_l, V_l,
                            lon_h, pa_h, pa_err_h, I_h, Lin_h, V_h,
                            outdir; show_=false, name_mod="PSR_NAME",
                            lon_rvm_l=nothing, pa_rvm_l=nothing,
                            pa_rvm_l_ortho=nothing, phi0_l=nothing,
                            lon_rvm_h=nothing, pa_rvm_h=nothing,
                            pa_rvm_h_ortho=nothing, phi0_h=nothing)

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(5.354337, 6.8))
        subplots_adjust(left=0.18, bottom=0.10, right=0.99, top=0.99, hspace=0.05)

        # Top panel: PA with error bars
        subplot2grid((3, 1), (0, 0))
        errorbar(lon_l, pa_l, yerr=pa_err_l, fmt=".", ms=2, elinewidth=0.6,
                 c="tab:blue",   label="low",  zorder=3)
        errorbar(lon_h, pa_h, yerr=pa_err_h, fmt=".", ms=2, elinewidth=0.6,
                 c="tab:orange", label="high", zorder=3)
        if !isnothing(lon_rvm_l) && !isnothing(pa_rvm_l)
            plot(lon_rvm_l, pa_rvm_l,       c="darkorange", lw=1.2, zorder=2)
            plot(lon_rvm_l, pa_rvm_l_ortho, c="darkorange", lw=1.2, zorder=2)
        end
        if !isnothing(lon_rvm_h) && !isnothing(pa_rvm_h)
            plot(lon_rvm_h, pa_rvm_h,       c="darkorange", lw=1.2, ls="--", zorder=2)
            plot(lon_rvm_h, pa_rvm_h_ortho, c="darkorange", lw=1.2, ls="--", zorder=2)
        end
        !isnothing(phi0_l) && axvline(phi0_l, c="tab:blue",   lw=1.5,        zorder=4)
        !isnothing(phi0_h) && axvline(phi0_h, c="tab:orange", lw=1.5, ls="--", zorder=4)
        ylabel("PA [deg]")
        ylim(-95, 95)
        yticks([-90, -45, 0, 45, 90])
        tick_params(labelbottom=false)
        minorticks_on()
        legend(fontsize=6, loc="upper right")

        # Bottom panel: Stokes I / L / V normalized to peak I
        subplot2grid((3, 1), (1, 0), rowspan=2)
        I_max_l = maximum(abs.(filter(isfinite, I_l)))
        I_max_h = maximum(abs.(filter(isfinite, I_h)))
        plot(lon_l, I_l   ./ I_max_l, c="black", lw=0.8,        label="I")
        plot(lon_l, Lin_l ./ I_max_l, c="red",   lw=0.8,        label="L")
        plot(lon_l, V_l   ./ I_max_l, c="royalblue", lw=0.8,    label="V")
        plot(lon_h, I_h   ./ I_max_h, c="black", lw=0.8, ls="--")
        plot(lon_h, Lin_h ./ I_max_h, c="red",   lw=0.8, ls="--")
        plot(lon_h, V_h   ./ I_max_h, c="royalblue", lw=0.8, ls="--")
        axhline(0.0, color="grey", lw=0.4, ls=":")
        xlabel("Longitude [deg]")
        ylabel("Flux Density [norm]")
        legend(fontsize=6, loc="upper right")
        minorticks_on()

        savepath = joinpath(outdir, "$(name_mod)_position_angle.pdf")
        println(savepath)
        savefig(savepath)
        savefig(replace(savepath, "pdf" => "png"))

        if show_
            PyPlot.show()
            println("Press Enter to close the figure.")
            readline(stdin; keep=false)
        end
        PyPlot.close()
    end


end  # module Plot
