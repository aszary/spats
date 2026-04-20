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
                 snr_threshold=3.5, linpol_threshold=0.8,
                 period=nothing, show_=false, n_freq=1)

        n_pulses_total = size(data4, 1)
        
        # Split data into low and high frequency halves if requested
        if n_freq == 2
            n_each = div(n_pulses_total, 2)
            data_low  = data4[1:n_each, :, :]
            data_high = data4[n_each+1:end, :, :]
        else
            data_low = data4
            data_high = nothing
        end

        # Create bin index arrays for full profile and on-pulse region
        n_full = size(data4, 2)       # total bins in full profile
        n_on   = bin_end - bin_st + 1  # bins in on-pulse region
        
        # Full profile longitude: centered at 0
        dl_full = 360.0  # full rotation in longitude
        lon_full = collect(range(-dl_full/2.0, dl_full/2.0, length=n_full))
        
        # On-pulse region longitude: subset
        dl_on   = 360.0 * n_on / n_full
        lon_on  = collect(range(-dl_on/2.0, dl_on/2.0, length=n_on))

        # **CHANGE: Fit RVM on COMBINED data for better SNR**
        # But still show LOW and HIGH separately on plots
        
        if n_freq == 2 && data_high !== nothing
            # Two-frequency mode: compute separate profiles but fit RVM on combined data
            println("[rvm] Two-frequency mode: showing LOW and HIGH profiles separately")
            
            # === Fit on COMBINED data for better SNR ===
            lon_f, pa_avg, pa_err, mask = 
                Tools.filter_ppa(data4, bin_st, bin_end;
                    snr_threshold=snr_threshold, linpol_threshold=linpol_threshold)
            
            alphas, betas, chi2_map, best_p =
                Tools.chi2_grid(lon_f, pa_avg, pa_err, mask)

            # Sigma-clip: remove PA points that deviate > 3σ from initial model
            model_0    = Tools.rvm_model(lon_f, best_p)
            res_0      = (pa_avg .- model_0) .- 180.0 .* round.((pa_avg .- model_0) ./ 180.0)
            sigma_0    = length(res_0[mask]) > 0 ? std(res_0[mask]) : 90.0
            mask_fit   = copy(mask)
            for i in eachindex(mask)
                mask[i] && abs(res_0[i]) > 3.0 * sigma_0 && (mask_fit[i] = false)
            end
            println("[rvm] σ-clip: $(sum(mask)) → $(sum(mask_fit)) PA bins kept")

            result = Tools.fit_rvm(lon_f, pa_avg, pa_err, mask_fit; p0=best_p)
            
            # Compute RVM curves: on-pulse region (for fitting) + full profile (for display)
            pa_curve_on = Tools.rvm_model(lon_on, 
                            [result.PA0, result.alpha, result.zeta, result.phi0])
            pa_curve_on = mod.(pa_curve_on .+ 90.0, 180.0) .- 90.0
            
            # RVM model on full profile for display
            pa_curve_full = Tools.rvm_model(lon_full,
                              [result.PA0, result.alpha, result.zeta, result.phi0])
            pa_curve_full = mod.(pa_curve_full .+ 90.0, 180.0) .- 90.0
            
            pa_curve_low = pa_curve_on   # For backward compat with on-pulse region
            pa_curve_high = nothing  # only one curve from combined fit
            
            # Combined (all-pulses) profile for display
            I_low_full  = vec(mean(data4[:, :, 1], dims=1))
            Q_mean_full = vec(mean(data4[:, :, 2], dims=1))
            U_mean_full = vec(mean(data4[:, :, 3], dims=1))
            L_low_full  = sqrt.(Q_mean_full.^2 .+ U_mean_full.^2)
            V_low_full  = vec(mean(data4[:, :, 4], dims=1))

            I_high_full = nothing
            L_high_full = nothing
            V_high_full = nothing

            I_max = maximum(abs.(I_low_full))
            I_on  = I_low_full

        else
            # Single-frequency mode
            println("[rvm] Fitting RVM for single frequency")
            
            lon_f, pa_avg, pa_err, mask = Tools.filter_ppa(data_low, bin_st, bin_end;
                snr_threshold=snr_threshold, linpol_threshold=linpol_threshold)

            alphas, betas, chi2_map, best_p = Tools.chi2_grid(lon_f, pa_avg, pa_err, mask)

            # Sigma-clip: remove PA points > 3σ from initial model
            model_0    = Tools.rvm_model(lon_f, best_p)
            res_0      = (pa_avg .- model_0) .- 180.0 .* round.((pa_avg .- model_0) ./ 180.0)
            sigma_0    = length(res_0[mask]) > 0 ? std(res_0[mask]) : 90.0
            mask_fit   = copy(mask)
            for i in eachindex(mask)
                mask[i] && abs(res_0[i]) > 3.0 * sigma_0 && (mask_fit[i] = false)
            end
            println("[rvm] σ-clip: $(sum(mask)) → $(sum(mask_fit)) PA bins kept")

            result = Tools.fit_rvm(lon_f, pa_avg, pa_err, mask_fit; p0=best_p)
            
            # Compute RVM curves on both on-pulse and full profile
            pa_curve_on = Tools.rvm_model(lon_on,
                            [result.PA0, result.alpha, result.zeta, result.phi0])
            pa_curve_on = mod.(pa_curve_on .+ 90.0, 180.0) .- 90.0
            
            pa_curve_full = Tools.rvm_model(lon_full,
                             [result.PA0, result.alpha, result.zeta, result.phi0])
            pa_curve_full = mod.(pa_curve_full .+ 90.0, 180.0) .- 90.0
            
            pa_curve_low = pa_curve_on    # For backward compat
            pa_curve_high = nothing
            
            # Full profile for display
            I_low_full = vec(mean(data_low[:, :, 1], dims=1))
            L_low_full = sqrt.(vec(mean(data_low[:, :, 2], dims=1)).^2 .+
                              vec(mean(data_low[:, :, 3], dims=1)).^2)
            V_low_full = vec(mean(data_low[:, :, 4], dims=1))
            
            I_high_full = nothing
            L_high_full = nothing
            V_high_full = nothing
            
            I_max = maximum(abs.(I_low_full))
            I_on = I_low_full
        end

        phi_c   = Tools.pulse_center_deg(lon_full, I_on)
        W10     = Tools.profile_width_deg(lon_full, I_on)
        h_blask = period !== nothing ?
                  Tools.emission_height_blaskiewicz(result.phi0, phi_c, period) : NaN

        h_contours = (period !== nothing && W10 > 0.0) ?
            Float64[Tools.height_from_rho_rankin(
                        Tools.rho_from_width(a, b, W10), period)
                    for a in alphas, b in betas] : nothing

        # ---------- figure -----------------------------------------------
        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.8)

        fig = figure(figsize=(7.5, 4.8))
        ax_pa  = fig.add_axes([0.09, 0.62, 0.37, 0.35])
        ax_st  = fig.add_axes([0.09, 0.11, 0.37, 0.44])
        ax_map = fig.add_axes([0.57, 0.11, 0.31, 0.86])
        ax_cb  = fig.add_axes([0.89, 0.11, 0.02, 0.86])

        # -- PA panel (left top) - on-pulse region only --
        # Identify sigma-clipped points (passed SNR/L-I but removed by sigma-clip)
        clipped = mask .& .!mask_fit

        # Failed SNR/L-I filter → light grey
        ax_pa.scatter(lon_f[.!mask], pa_avg[.!mask],
                      s=2, color="lightgrey", alpha=0.3, zorder=1)
        # Sigma-clipped outliers → open grey circles
        if any(clipped)
            ax_pa.errorbar(lon_f[clipped], pa_avg[clipped], yerr=pa_err[clipped],
                           fmt="o", ms=2.5, mfc="none", color="grey", alpha=0.5,
                           elinewidth=0.4, capsize=0.5, zorder=2)
        end
        # Points used in fit → black
        ax_pa.errorbar(lon_f[mask_fit], pa_avg[mask_fit], yerr=pa_err[mask_fit],
                       fmt="o", ms=2.5, color="black", alpha=0.9,
                       elinewidth=0.5, capsize=0.5, zorder=3)
        # RVM curve
        ax_pa.plot(lon_on, pa_curve_low, color="darkorange", lw=1.5, zorder=4)
        ax_pa.axvline(result.phi0, color="royalblue", lw=0.8, ls="--", alpha=0.7)

        txt_str = @sprintf("α=%.1f°  ζ=%.1f°  φ₀=%.1f°\nχ²ᵣ=%.2f  rms=%.1f°",
                           result.alpha, result.zeta, result.phi0,
                           result.chi2_red, result.rms_deg)
        
        ax_pa.set_xlim(lon_on[1], lon_on[end])
        ax_pa.set_ylim(-95, 95)
        ax_pa.set_ylabel("PA [deg]")
        ax_pa.set_yticks([-90, -45, 0, 45, 90])
        ax_pa.tick_params(labelbottom=false)
        ax_pa.minorticks_on()
        ax_pa.text(0.97, 0.97, txt_str,
                   transform=ax_pa."transAxes", fontsize=6,
                   va="top", ha="right",
                   bbox=Dict("boxstyle"=>"round", "fc"=>"wheat", "alpha"=>0.7))
        if I_high_full !== nothing && n_freq == 2
            # Show both frequencies separately with full profile
            ax_st.plot(lon_full, I_low_full ./ I_max, color="black",     lw=0.8, label="I low", alpha=0.7)
            ax_st.plot(lon_full, L_low_full ./ I_max, color="red",       lw=0.8, label="L low", alpha=0.7)
            ax_st.plot(lon_full, V_low_full ./ I_max, color="royalblue", lw=0.8, label="V low", alpha=0.7)
            ax_st.plot(lon_full, I_high_full ./ I_max, color="black",     lw=0.8, label="I high", ls="--", alpha=0.7)
            ax_st.plot(lon_full, L_high_full ./ I_max, color="red",       lw=0.8, label="L high", ls="--", alpha=0.7)
            ax_st.plot(lon_full, V_high_full ./ I_max, color="royalblue", lw=0.8, label="V high", ls="--", alpha=0.7)
        else
            # Show combined profile (full)
            ax_st.plot(lon_full, I_low_full ./ I_max, color="black",     lw=1.0, label="I")
            ax_st.plot(lon_full, L_low_full ./ I_max, color="red",       lw=1.0, label="L")
            ax_st.plot(lon_full, V_low_full ./ I_max, color="royalblue", lw=1.0, label="V")
        end
        ax_st.axhline(0.0, color="grey", lw=0.4, ls=":")
        # Show on-pulse region plus 30% margin on each side so pulse edges are visible
        lon_margin = 0.3 * (lon_on[end] - lon_on[1])
        ax_st.set_xlim(lon_on[1] - lon_margin, lon_on[end] + lon_margin)
        ax_st.set_ylim(-0.5, 1.15)
        ax_st.set_ylabel("Flux Density [norm.]")
        ax_st.set_xlabel("Longitude [deg]")
        ax_st.legend(fontsize=6, loc="upper right", framealpha=0.6, ncol=3)
        ax_st.minorticks_on()
        if !isnan(h_blask)
            ax_st.text(0.03, 0.03,
                       @sprintf("h_Blask=%.0f km  Delphi=%.1f deg",
                                h_blask, result.phi0 - phi_c),
                       transform=ax_st."transAxes", fontsize=6, va="bottom",
                       bbox=Dict("boxstyle"=>"round", "fc"=>"lightyellow", "alpha"=>0.7))
        end

        # -- chi²(alpha, beta) map (right panel): alpha on x, beta on y --
        # Δχ²ᵣ = χ²ᵣ − min(χ²ᵣ), auto-scaled so the valley is visible
        chi2_finite = filter(isfinite, vec(chi2_map))
        chi2_min    = isempty(chi2_finite) ? 0.0 : minimum(chi2_finite)
        chi2_delta  = map(v -> isnan(v) ? NaN : max(v - chi2_min, 0.0), chi2_map)

        best = argmin(map(v -> isnan(v) ? Inf : v, chi2_map))
        a_best = alphas[best[1]]
        b_best = betas[best[2]]

        # Auto-scale vmax: 2% of the finite range, clamped to [5, 50]
        chi2_d_finite = filter(isfinite, vec(chi2_delta))
        dmax = isempty(chi2_d_finite) ? 10.0 :
               clamp(0.02 * maximum(chi2_d_finite), 5.0, 50.0)

        # Mask values above dmax as NaN so they render as white (out-of-range color)
        chi2_display = map(v -> isnan(v) ? NaN : (v > dmax ? NaN : v), chi2_delta)

        # Zoom display around the minimum (±60° in α, ±10° in β) so the valley fills the panel
        a_lo = max(alphas[1],   a_best - 80.0)
        a_hi = min(alphas[end], a_best + 80.0)
        b_lo = max(betas[1],    b_best - 10.0)
        b_hi = min(betas[end],  b_best + 10.0)

        im = ax_map.imshow(chi2_display',
                           origin="lower", aspect="auto",
                           extent=[alphas[1], alphas[end], betas[1], betas[end]],
                           cmap="viridis_r", vmin=0.0, vmax=dmax,
                           interpolation="nearest")
        cb = fig.colorbar(im, cax=ax_cb)
        cb.set_label("Δχ²ᵣ")

        ax_map.scatter([a_best], [b_best],
                       marker="*", s=100, color="red", zorder=6,
                       label=@sprintf("α=%.0f° β=%.0f°", a_best, b_best))

        # Confidence-region contours: Δχ²=2.30 (1σ), 6.17 (2σ) for 2 free params
        try
            ax_map.contour(alphas, betas, chi2_delta',
                           levels=[2.30, 6.17],
                           colors=["cyan", "royalblue"],
                           linewidths=[1.0, 0.8], linestyles=["solid", "dashed"],
                           zorder=5)
        catch
        end

        ax_map.legend(fontsize=6, loc="upper right", framealpha=0.6)

        if h_contours !== nothing
            try
                h_levs = _auto_height_levels(h_contours)
                if !isempty(h_levs)
                    ax_map.clabel(
                        ax_map.contour(alphas, betas, h_contours',
                                       levels=h_levs,
                                       colors="purple", linewidths=0.7),
                        fmt="%g", fontsize=5, inline=true)
                end
            catch
            end
        end

        ax_map.set_xlim(a_lo, a_hi)
        ax_map.set_ylim(b_lo, b_hi)
        ax_map.set_xlabel("α [deg]")
        ax_map.set_ylabel("β [deg]")
        ax_map.minorticks_on()

        savefig("$outdir/$(name_mod)_rvm.pdf")
        savefig("$outdir/$(name_mod)_rvm.png")
        println("$outdir/$(name_mod)_rvm.pdf")

        if show_
            PyPlot.show()
            println("Press Enter to close.")
            readline(stdin; keep=false)
        end
        PyPlot.close()

        return merge(result, (h_blask=h_blask, phi_center=phi_c, W10=W10,
                              alphas=alphas, betas=betas, chi2_map=chi2_map))
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
