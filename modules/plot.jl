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

    function analyse_single(low, high, mid, p)
        pulses, bins = size(low)
        for i in 1:pulses

            figure(figsize=(6, 7))  # figure size in inches (~15.24 cm × 17.78 cm)
            plot(low[i, p["bin_st"]:p["bin_end"]], label="low")
            plot(high[i, p["bin_st"]:p["bin_end"]], label="high")
            plot(mid[i, p["bin_st"]:p["bin_end"]], label="mid")

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



end  # module Plot
