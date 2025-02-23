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
    
    PyPlot.matplotlib.use("Tkagg") # DOES NOT WORK on ozStar! had to set backend in matplotlib by hand
    #PyPlot.matplotlib.use("qt5agg")
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


    function average(data, outdir; start=1, number=100, bin_st=nothing, bin_end=nothing, name_mod="0")
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
        savefig("$outdir/average_$name_mod.pdf")
        close()
        #clf()
    end


    function average2(data, data2, outdir; start=1, number=100, bin_st=nothing, bin_end=nothing, name_mod="0")
        num, bins = size(data)
        num2, bins = size(data2)
        if number == nothing
            number1 = num - start  # missing one?
            number2 = num2 - start  # missing one?
        else
            number1 = number
            number2 = number
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number1-1, bin_st:bin_end]
        da2 = data2[start:start+number2-1, bin_st:bin_end]
        average = Tools.average_profile(da)
        average2 = Tools.average_profile(da2)

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 2.362205), frameon=true)  # 8cm x 6 cm
        subplots_adjust(left=0.14, bottom=0.14, right=0.99, top=0.99, wspace=0., hspace=0.)

        minorticks_on()
        plot(longitude, average, c="C1", label="low")
        plot(longitude, average2, c="C2", label="high")
        #yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        xlabel("longitude \$(^\\circ)\$")
        legend()
        #tick_params(labeltop=false, labelbottom=true)
        savefig("$outdir/$(name_mod)_average2.pdf")
        close()
        #clf()
    end



    function averageX(datas, outdir; start=1, number=100, bin_st=nothing, bin_end=nothing, name_mod="0")
        nums = []
        bins = []
        for data in datas
            nu, bi = size(data)
            push!(nums, nu)
            push!(bins, bi)
        end
        numbers = []
        if number == nothing
            for num in nums
                push!(numbers, num-start)  # missing one?
            end
        else
            for i in 1:length(nums)
                push!(numbers, number)
            end
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins[1] end
        das = []
        for (i,data) in enumerate(datas)
            da = data[start:start+numbers[i]-1, bin_st:bin_end]
            push!(das, da)
        end
        avs = []
        for da in das
            push!(avs, Tools.average_profile(da))
        end

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins[1]
        longitude = collect(range(-dl/2., dl/2., length=db))

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 2.362205), frameon=true)  # 8cm x 6 cm
        subplots_adjust(left=0.14, bottom=0.14, right=0.99, top=0.99, wspace=0., hspace=0.)

        minorticks_on()
        for i in 1:length(avs)
            plot(longitude, avs[i], c="C$i", label="Obs. num. $i")
        end
        #yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        xlabel("longitude \$(^\\circ)\$")
        legend()
        #tick_params(labeltop=false, labelbottom=true)
        savefig("$outdir/$(name_mod)_averageX.pdf")
        close()
        #clf()
    end




    function single0(data, outdir; start=1, number=100, bin_st=nothing, bin_end=nothing, norm=2.0)

        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end

        figure(figsize=(3.14961, 4.33071))  # 8cm x 11cm
        subplots_adjust(left=0.24, bottom=0.13, right=0.99, top=0.99, wspace=0., hspace=0.)
        minorticks_on()
        for i in start:start+number-1
            da = data[i,:] .* norm .+ i
            da = da[bin_st:bin_end]
            plot(da, c="grey", lw=0.3)
        end
        #PyPlot.xlim(bin_st, bin_end)
        ylim(start, start+number)
        xlabel("bin number")
        ylabel("Pulse number")
        filename = "$outdir/single0.pdf"
        println(filename)
        savefig(filename)
        close()
    end


    function single(data, outdir; start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME", show_=false)
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
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
        #clf()
    end


    function singles(datas, outdir; start=1, number=nothing, cmap="viridis", bin_st=[nothing, nothing], bin_end=[nothing, nothing], darkness=[0.5, 0.7], name_mod="PSR_NAME", show_=false)
        num, bins1 = size(datas[1])
        num2, bins2 = size(datas[2])
        if number == nothing
            number = num - start  # missing one?
        end

        noise_start = 100

        if bin_st[1] == nothing bin_st[1] = 1 end
        if bin_st[2] == nothing bin_st[2] = 1 end
        if bin_end[1] == nothing bin_end[1] = bins end
        if bin_end[2] == nothing bin_end[2] = bins end
        da1 = datas[1][start:start+number-1,bin_st[1]:bin_end[1]]
        da1full = datas[1][start:start+number-1,1:end]
        average1 = Tools.average_profile(da1)
        intensity1, pulses1 = Tools.intensity_pulses(da1)
        intensity1 .-= minimum(intensity1)
        intensity1 ./= maximum(intensity1)
        std1 = std(da1full[noise_start:noise_start + (bin_end[1]-bin_st[1])])
        println("std1 ", std1)
        average1 .-= mean(da1full[noise_start:noise_start + (bin_end[1]-bin_st[1])]) # minimum(average1) # not 0
        #average1 .-= minimum(average1)
        average1 ./= maximum(average1)

        da2 = datas[2][start:start+number-1,bin_st[2]:bin_end[2]]
        da2full = datas[2][start:start+number-1,1:end]
        average2 = Tools.average_profile(da2)
        intensity2, pulses2 = Tools.intensity_pulses(da2)
        intensity2 .-= minimum(intensity2)
        intensity2 ./= maximum(intensity2)
        std2 = std(da2full[noise_start:noise_start + (bin_end[2]-bin_st[2])])
        average2 .-= mean(da2full[noise_start:noise_start+(bin_end[2]-bin_st[2])]) # minimum(average2) # not 0
        #average2 .-= minimum(average2)
        average2 ./= maximum(average2)

        pulses1 .+= start - 1  # julia
        pulses2 .+= start - 1  # julia

        # a bad W_10 estimate with some profile smoothing (average)
        bi = []
        bi2 = []
        for i in 1:length(average1)-1
            if (average1[i] + average1[i+1])/2 > 0.1 +std1
                push!(bi, i)
            end
            if (average1[i] + average1[i+1])/2 > 0.1 - std1
                push!(bi2, i)
            end
        end
        println("W_10 Obs. min. ", (maximum(bi) - minimum(bi)) / bins1 * 360)
        println("W_10 Obs. max. ", (maximum(bi2) - minimum(bi2)) / bins1 * 360)
        #println(maximum(average1), " ", minimum(average1))
        #println(maximum(bi), " ", minimum(bi))

        bi = []
        bi2 = []
        for i in 1:length(average2)-1
            if (average2[i] + average2[i+1]) / 2 > 0.1 + std1
                push!(bi, i)
            end
            if (average2[i] + average2[i+1]) / 2 > 0.1 - std1
                push!(bi2, i)
            end
        end
        println("W_10 Mod. min. ", (maximum(bi) - minimum(bi)) / bins2 * 360)
        println("W_10 Mod. max. ", (maximum(bi2) - minimum(bi2)) / bins2 * 360)


        # Pulse longitude
        db1 = (bin_end[1] + 1) - bin_st[1]  # yes +1
        dl1 = 360. * db1 / bins1
        longitude1 = collect(range(-dl1/2., dl1/2., length=db1))

        db2 = (bin_end[2] + 1) - bin_st[2]  # yes +1
        dl2 = 360. * db2 / bins2
        longitude2 = collect(range(-dl2/2., dl2/2., length=db2))


        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 4.33071), frameon=true)  # 8cm x 11cm
        subplots_adjust(left=0.13, bottom=0.09, right=0.99, top=0.99, wspace=0., hspace=0.)

        subplot2grid((5, 2), (0, 0), rowspan=4)
        imshow(da1, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness[1]*maximum(da1), vmin=0.0)
        minorticks_on()
        tick_params(labelleft=true, labelbottom=false)

        subplot2grid((5, 2), (0, 1), rowspan=4)
        imshow(da2, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness[2]*maximum(da2), vmin=0.0)
        tick_params(labelleft=false, labelbottom=false)

        subplot2grid((5, 2), (4, 0))
        minorticks_on()
        plot(longitude1, average1, c="grey")
        yticks([0.0, 0.5])
        xlim(longitude1[1], longitude1[end])
        xlabel("longitude \$(^\\circ)\$")

        subplot2grid((5, 2), (4, 1))
        minorticks_on()
        plot(longitude2, average2, c="grey")
        yticks([0.0, 0.5])
        xlim(longitude2[1], longitude2[end])
        tick_params(labelleft=false)
        xlabel("longitude \$(^\\circ)\$")


        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_singles.pdf")
        savefig("$outdir/$(name_mod)_singles.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
        #clf()
    end



    function single_J1750(data, outdir; start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME", show_=false, panel="a")
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

        pulses .+= start - 1  # julia

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        rc("font", size=7.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        #figure(figsize=(2.362205, 3.248031), frameon=true)  # 6cm x 8.25cm
        figure(figsize=(3.14961, 4.33071))  # 8cm x 11cm
        if number > 999
            subplots_adjust(left=0.235, bottom=0.115, right=0.99, top=0.99, wspace=0., hspace=0.)
        else
            subplots_adjust(left=0.21, bottom=0.115, right=0.99, top=0.99, wspace=0., hspace=0.)
        end
        #figtext(0.01, 0.95, "$panel)", size=10)
        subplot2grid((5, 3), (0, 0), rowspan=4)
        minorticks_on()
        plot(intensity, pulses, c="grey")
        ylim(pulses[1]-0.5, pulses[end]+0.5)
        xticks([0.5, 1.0])
        xlim(1.1, -0.1)
        xlabel("intensity", labelpad=7)
        ylabel("Pulse number")

        subplot2grid((5, 3), (0, 1), rowspan=4, colspan=2)
        imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da), vmin=0)
        #axvline(x=563, lw=2)
        tick_params(labelleft=false, labelbottom=false)

        subplot2grid((5, 3), (4, 1), colspan=2)
        minorticks_on()
        plot(longitude, average, c="grey")
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        xlabel("longitude (\$^{\\circ}\$)")
        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_single.pdf")
        savefig("$outdir/$(name_mod)_single.pdf")
        savefig("$outdir/$(name_mod)_single.svg")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
        #clf()
    end



    function lrfs(data, outdir; start=1, number=nothing, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0", change_fftphase=true)

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
        println("$outdir/$(name_mod)_lrfs.pdf")
        savefig("$outdir/$(name_mod)_lrfs.pdf")
        savefig("$outdir/$(name_mod)_lrfs.svg")
        close()
    end


    function lrfses(datas, outdir; start=1, number=nothing, cmap="viridis", bin_st=[nothing, nothing], bin_end=[nothing, nothing], darkness=0.5, name_mod="0")

        num, bins1 = size(datas[1])
        num, bins2 = size(datas[2])
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st[1] == nothing bin_st[1] = 1 end
        if bin_st[2] == nothing bin_st[2] = 1 end
        if bin_end[1] == nothing bin_end[1] = bins end
        if bin_end[2] == nothing bin_end[2] = bins end

        da1 = datas[1][start:start+number-1,bin_st[1]:bin_end[1]]
        average1 = Tools.average_profile(da1)
        lrfs1, intensity1, freq1, peak1 = Tools.lrfs(da1)
        println("\tpeak freq [1] $(freq1[peak1]) ")
        println("\tpeak P3 [1] $(1/freq1[peak1])")
        phase1_ = rad2deg.(angle.(view(lrfs1, peak1, :)))  # fft phase variation  # view used! lrfs[peak, :] -> copies data

        da2 = datas[2][start:start+number-1,bin_st[2]:bin_end[2]]
        average2 = Tools.average_profile(da2)
        lrfs2, intensity2, freq2, peak2 = Tools.lrfs(da2)
        println("\tpeak freq [2] $(freq2[peak2]) ")
        println("\tpeak P3 [2] $(1/freq2[peak2])")
        phase2_ = rad2deg.(angle.(view(lrfs2, peak2, :)))  # fft phase variation  # view used! lrfs[peak, :] -> copies data
        # skip freq = 0 and normalize intensity to 1
        inten1 = intensity1[2:end]
        inten1 .-= minimum(inten1)
        inten1 ./= maximum(inten1)
        fre1 = freq1[2:end]

        inten2 = intensity2[2:end]
        inten2 .-= minimum(inten2)
        inten2 ./= maximum(inten2)
        fre2 = freq2[2:end]

        pars1, errs1 = Tools.fit_gaussian(fre1, inten1; μ=freq1[peak1-1])  # skip zero freq
        fr1 = pars1[2]
        frer1 = pars1[3]  # errs[2] # yeap
        println("\tFrequency [1](gaussian fit): $fr1, P3: $(1/fr1)")
        println("\tFrequency error [1] (gaussian fit): $frer1, P3 error: $(1/fr1 - 1/(fr1 +frer1))")

        pars2, errs2 = Tools.fit_gaussian(fre2, inten2; μ=freq2[peak2-1])  # skip zero freq
        fr2 = pars2[2]
        frer2 = pars2[3]  # errs[2] # yeap
        println("\tFrequency [2](gaussian fit): $fr2, P3: $(1/fr2)")
        println("\tFrequency error [2] (gaussian fit): $frer2, P3 error: $(1/fr2 - 1/(fr2 +frer2))")

        # Pulse longitude
        db1 = (bin_end[1] + 1) - bin_st[1]  # yes +1
        dl1 = 360. * db1 / bins1
        longitude1 = collect(range(-dl1/2., dl1/2., length=db1))

        db2 = (bin_end[2] + 1) - bin_st[2]  # yes +1
        dl2 = 360. * db2 / bins2
        longitude2 = collect(range(-dl2 / 2., dl2 / 2., length=db2))

        # longitude vs. FFT phase analysis
        #=
        lin, incs = Tools.find_inclination(longitude, phase_)
        for (i,inc) in enumerate(incs)
            println("($i) Δ FFT phase / Δ longitude ", round(inc[2], digits=2))
            println("($i) sigma Δ FFT phase / Δ longitude ", round(inc[3], digits=2))
        end
        =#

        lin1, inc1 = Tools.find_inclination2(longitude1, phase1_)
        println(" Δ FFT phase / Δ longitude [1] ", round(inc1[1][2], digits=2))
        println("sigma Δ FFT phase / Δ longitude [1] ", round(inc1[1][3], digits=2))

        lin2, inc2 = Tools.find_inclination2(longitude2, phase2_)
        println(" Δ FFT phase / Δ longitude [2] ", round(inc2[1][2], digits=2))
        println("sigma Δ FFT phase / Δ longitude [2] ", round(inc2[1][3], digits=2))


        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961*2, 4.33071))  # 8cm x 11cm
        subplots_adjust(left=0.08, bottom=0.1, right=0.92, top=0.9, wspace=0., hspace=0.)

        # phase vs longitude [1]
        ax = subplot2grid((4, 6), (0, 1), colspan=2)
        minorticks_on()
        scatter(longitude1, phase1_, marker=".", c="grey", s=3.)
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        for l in lin1
            plot(l[1], l[2], c="blue")
        end
        text(-27, -1000, string("\$\\frac{\\Delta {\\rm FFT\\, ph.}}{\\Delta {\\rm lon.}} =", round(inc1[1][2], digits=1), " \\pm ", round(inc1[1][3], digits=1), "\$"))
        ylim(-1200, 200)
        xlabel("longitude \$(^\\circ)\$")
        ylabel("FFT phase \$(^\\circ)\$")
        tick_params(labeltop=true, labelbottom=false, which="both", bottom=false, top=true)

        # frequency [1]
        ax = subplot2grid((4, 6), (1, 0), rowspan=3)
        minorticks_on()
        plot(inten1, fre1, c="grey")
        #plot(Tools.gauss(fre1, pars1), fre1, c="red", ls=":", lw=0.3)
        axhline(y=freq1[peak1], c="red", ls="--")
        ylim(freq1[2], freq1[end])
        xticks([0.5, 1.0])
        xlim(1.1, -0.1)
        xlabel("intensity")
        ylabel("frequency \$(1/P)\$")

        # lrfs [1]
        ax = subplot2grid((4, 6), (1, 1), rowspan=3, colspan=2)
        minorticks_on()
        imshow(abs.(lrfs1[2:end,1:end]), origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness[1]*maximum(abs.(lrfs1[2:end,1:end])), extent=[longitude1[1], longitude1[end], freq1[1], freq1[end]])
        xlabel("longitude \$(^\\circ)\$")
        tick_params(labelleft=false, labelbottom=true)

        # phase vs longitude [2]
        off = 100
        ax = subplot2grid((4, 6), (0, 3), colspan=2)
        minorticks_on()
        scatter(longitude2, phase2_.-off, marker=".", c="grey", s=3.)
        for l in lin2
            plot(l[1], l[2].-off, c="blue")
        end
        text(-20, -1000, string("\$\\frac{\\Delta {\\rm FFT\\, ph.}}{\\Delta {\\rm lon.}} =", round(inc2[1][2], digits=1), " \\pm ", round(inc2[1][3], digits=1), "\$"))
        ylim(-1200, 200)
        xlabel("longitude \$(^\\circ)\$")
        ylabel("FFT phase \$(^\\circ)\$")
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        ax.yaxis.set_label_position("right")
        tick_params(labeltop=true, labelbottom=false, which="both", bottom=false, top=true)
        tick_params(labelright=true, labelleft=false, which="both", right=true, left=false)

        # frequency [2]
        ax = subplot2grid((4, 6), (1, 5), rowspan=3)
        minorticks_on()
        plot(inten2, fre2, c="grey")
        #plot(Tools.gauss(fre2, pars2), fre2, c="red", ls=":", lw=0.3)
        axhline(y=freq2[peak2], c="red", ls="--")
        ylim(freq2[2], freq2[end])
        xticks([0.5, 1.0])
        xlim(-0.1, 1.1)
        xlabel("intensity")
        ylabel("frequency \$(1/P)\$")
        ax.yaxis.set_label_position("right")
        tick_params(labelbottom=true, which="both", right=true, labelleft=false, labelright=true)

        # lrfs [2]
        ax = subplot2grid((4, 6), (1, 3), rowspan=3, colspan=2)
        minorticks_on()
        imshow(abs.(lrfs2[2:end,1:end]), origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness[2]*maximum(abs.(lrfs2[2:end,1:end])), extent=[longitude1[1], longitude1[end], freq1[1], freq1[end]])
        xlabel("longitude \$(^\\circ)\$")
        tick_params(labelleft=false, labelbottom=true)
        savefig("$outdir/$(name_mod)_lrfses.pdf")
        close()
    end



    function p3_evolution(data, outdir; start=1, end_=nothing, step=10, number=256, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0", verbose=false)
        num, bins = size(data)
        if end_ == nothing end_ = num end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        intensity_ = []
        p3_ = []
        p3_err_ = []
        start_period = []
        for i in start:step:end_-number
            global p3
            global p3err
            global frequency
            da = data[i:i+number-1,bin_st:bin_end]
            lrfs, intensity, freq, peak = Tools.lrfs(da)
            # skip freq = 0 and normalize intensity to 1
            inten = intensity[2:end]
            inten ./= maximum(inten)
            fre = freq[2:end]
            try
                pars, errs = Tools.fit_gaussian(fre, inten; μ=freq[peak-1])  # skip zero freq
                f = pars[2]
                #fer = abs(pars[3])  # nope too big
                fer = abs(errs[2])  # yeap
                p3 = 1 / f
                p3err = maximum([1 / f - 1 / (f +fer), 1 / (f - fer) - 1 / f])
                if verbose == true println("\tP3 = $p3, P3 error = $p3err") end
            catch exc
                p3 = 0. #nothing
                p3err = 0. #nothing
                if verbose == true println("\t[WARNING! P3 = 0, P3 error = 0]") end
            end
            push!(intensity_, inten)
            push!(p3_, p3)
            push!(p3_err_, p3err)
            push!(start_period, i)
            frequency = fre
        end

        da = data[start:start+end_-1,bin_st:bin_end]
        average = Tools.average_profile(da)
        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        # converting intensity TODO why why why?
        x, = size(intensity_)
        y, = size(intensity_[1])
        intens = zeros((x,y))
        for i in 1:x
            for j in 1:y
                intens[i,j] = intensity_[i][j]
            end
        end

        #println(typeof(intensity_))
        #println(typeof(intens))
        #return

        #left, skip = Tools.intensity_pulses(intens)
        bottom, skip = Tools.intensity_pulses(transpose(intens))

        #println(x, " ", y)
        #println(size(intensity_))

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 4.33071))  # 8cm x 11cm
        subplots_adjust(left=0.16, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)

        subplot2grid((5, 3), (0, 0), rowspan=4)
        minorticks_on()
        errorbar(p3_, start_period, xerr=p3_err_, color="none", lw=1., marker="_", mec="grey", ecolor="grey", capsize=0, mfc="grey", ms=1.0)
        ylim(start_period[1], start_period[end])
        xlabel("\$P_3\$")
        ylabel("start period num.")

        subplot2grid((5, 3), (0, 1), rowspan=4, colspan=2)
        imshow(intens, origin="lower", cmap=cmap, interpolation="none", aspect="auto") #,  vmax=darkness*maximum(intens))
        tick_params(labelleft=false, labelbottom=false)

        subplot2grid((5, 3), (4, 1), colspan=2)
        minorticks_on()
        plot(frequency, bottom, c="grey")
        xlabel("frequency\$(1/P)\$")
        #tick_params(labeltop=false, labelbottom=true)
        xlim(frequency[1], frequency[end])
        yticks([])
        savefig("$outdir/$(name_mod)_p3_evolution.pdf")
        close()
    end


    function p3_evolution_J1750(data, outdir; panel="a", start=1, end_=nothing, step=10, number=128, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="1", verbose=false)
        num, bins = size(data)
        if end_ == nothing end_ = num end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        intensity_ = []
        p3_ = []
        p3_err_ = []
        start_period = []
        for i in start:step:end_-number
            global p3
            global p3err
            global frequency
            da = data[i:i+number-1,bin_st:bin_end]
            lrfs, intensity, freq, peak = Tools.lrfs(da)
            # skip freq = 0 and normalize intensity to 1
            inten = intensity[2:end]
            inten ./= maximum(inten)
            fre = freq[2:end]
            #inten[1] = inten[2]
            try
                pars, errs = Tools.fit_gaussian_J1750(fre, inten, i) #; μ=freq[peak-1])  # skip zero freq
                f = pars[2]
                p3 = 1 / f
                #fer = abs(pars[3])   # nope # too big?
                fer = abs(errs[2])  # yeap
                p3err = maximum([1 / f - 1 / (f +fer), 1 / (f - fer) - 1 / f])
                # try other approach to estimate error
                if p3err > 30
                    fer = abs(pars[3])
                    p3err2 = maximum([1 / f - 1 / (f +fer), 1 / (f - fer) - 1 / f])
                    if p3err2 < p3err
                        p3err = p3err2
                    end
                    #println("new $p3err $p3err2")
                end
                if verbose == true println("\t$i P3 = $p3, P3 error = $p3err") end
            catch exc
                if typeof(exc) == DivideError
                    break  # good job
                end
                p3 = 0 # nothing
                p3err = 0 # nothing
                if verbose == true println("\t[WARNING! P3 = 0, P3 error = 0]") end
            end
            push!(intensity_, inten)
            if ((p3 > 20) && (p3 <  200)) && p3err < 30
                push!(p3_, p3)
                push!(p3_err_, p3err)
                push!(start_period, i)
            end
            frequency = fre
        end
        if verbose == true println("P3 mean: $(mean(p3_))  P3 std: $(std(p3_))") end

        da = data[start:start+end_-1,bin_st:bin_end]
        average = Tools.average_profile(da)
        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        # converting intensity TODO why why why?
        x, = size(intensity_)
        y, = size(intensity_[1])
        intens = zeros((x,y))
        for i in 1:x
            for j in 1:y
                intens[i,j] = intensity_[i][j]
            end
        end

        #println(typeof(intensity_))
        #println(typeof(intens))
        #return

        #left, skip = Tools.intensity_pulses(intens)
        bottom, skip = Tools.intensity_pulses(transpose(intens))

        #println(x, " ", y)
        #println(size(intensity_))

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(2.362205, 3.248031))  # 6cm x 8.25cm
        if num > 999
            subplots_adjust(left=0.235, bottom=0.115, right=0.99, top=0.99, wspace=0., hspace=0.)
        else
            subplots_adjust(left=0.21, bottom=0.115, right=0.99, top=0.99, wspace=0., hspace=0.)
        end
        figtext(0.01, 0.95, "$panel)", size=10)

        subplot2grid((5, 3), (0, 0), rowspan=4)
        minorticks_on()
        locator_params(nbins=5)
        errorbar(p3_, start_period, xerr=p3_err_, color="none", lw=0.3, marker="_", mec="grey", ecolor="grey", capsize=0, mfc="grey", ms=0.7)
        #ylim(start_period[1], start_period[end])
        ylim(start, start+end_-number)
        #xlim(31, 55)
        if (panel != "a") && (panel != "d")
            xticks([30, 80])
        end
        xlabel("\$P_3\$")
        ylabel("start period num.")

        subplot2grid((5, 3), (0, 1), rowspan=4, colspan=2)
        minorticks_on()
        imshow(intens, origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(intens))
        tick_params(labelleft=false, labelbottom=false)
        println(size(intens))

        subplot2grid((5, 3), (4, 1), colspan=2)
        minorticks_on()
        plot(frequency, bottom, c="grey")
        xlabel("frequency \$(1/P)\$")
        #tick_params(labeltop=false, labelbottom=true)
        xlim(frequency[1]/2., frequency[end]) # nice trick
        yticks([])
        println("$outdir/$(name_mod)_p3_evolution.pdf")
        savefig("$outdir/$(name_mod)_p3_evolution.pdf")
        close()

    end



    function p3fold(data, outdir; start=1, number=nothing, repeat_num=4, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0")
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
        close()
    end


    function p3fold_two(data1, data2, outdir; start=1, number=nothing, repeat_num=4, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0", shift=7)
        num, bins = size(data1)
        if number == nothing
            number = num - start + 1  # missing one? # test it!
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data1[start:start+number-1,bin_st:bin_end]
        da2 = data2[start:start+number-1,bin_st:bin_end]
        average = Tools.average_profile(da)
        average2 = Tools.average_profile(da2)
        intensity, pulses = Tools.intensity_pulses(da)
        intensity2, pulses2 = Tools.intensity_pulses(da2)
        intensity .-= minimum(intensity)
        intensity ./= maximum(intensity)
        intensity2 .-= minimum(intensity2)
        intensity2 ./= maximum(intensity2)

        pulses .+= start

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        #longitude = collect(range(-dl/2., dl/2., length=db)) # zero in the middle
        longitude = collect(range(bin_st/bins*360, bin_end/bins*360, length=db))# raw longitude

        # repeat data
        da = repeat(da, repeat_num)
        da2 = repeat(da2, repeat_num)
        # shift the data
        one = da2[1+shift:end, :]
        two = da2[1:shift, :]
        da2 = vcat(one, two)
        intensity = repeat(intensity, repeat_num)
        intensity2 = repeat(intensity2, repeat_num)
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
        #left
        subplot2grid((5, 2), (0, 0), rowspan=4)
        minorticks_on()
        imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da))
        tick_params(labelleft=true, labelbottom=false)
        yticks(ticks, ti)
        # right
        subplot2grid((5, 2), (0, 1), rowspan=4)
        minorticks_on()
        imshow(da2, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=0.7*darkness*maximum(da2))
        tick_params(labelleft=false, labelbottom=false)
        yticks(ticks, ti)
        # bottom left
        subplot2grid((5, 2), (4, 0))
        minorticks_on()
        plot(longitude, average, c="grey")
        #xticks([-40, 0, 40])
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        ylabel("intensity (a. u.)")
        # bottom right
        subplot2grid((5, 2), (4, 1), colspan=3)
        minorticks_on()
        tick_params(labelleft=false)
        plot(longitude, average2, c="grey")
        #xticks([-40, 0, 40])
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        figtext(0.45, 0.01, "longitude (deg.)", size=8)
        #tick_params(labeltop=false, labelbottom=true)
        savefig("$outdir/$(name_mod)_p3fold_two.pdf")
        println("$outdir/$(name_mod)_p3fold_two.pdf")
        close()
    end


    function p3fold_twotracks(data1, data2, data_dir, outdir; start=1, number=nothing, repeat_num=4, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0", shift=14, show_=false)
        num, bins = size(data1)
        if number == nothing
            number = num - start + 1 # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data1[start:start+number-1,bin_st:bin_end]
        da2 = data2[start:start+number-1,bin_st:bin_end]
        average = Tools.average_profile(da)
        average2 = Tools.average_profile(da2)
        intensity, pulses = Tools.intensity_pulses(da)
        intensity2, pulses2 = Tools.intensity_pulses(da2)
        intensity .-= minimum(intensity)
        intensity ./= maximum(intensity)
        intensity2 .-= minimum(intensity2)
        intensity2 ./= maximum(intensity2)

        pulses .+= start

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        #longitude = collect(range(-dl/2., dl/2., length=db))
        longitude = collect(range(bin_st/bins*360, bin_end/bins*360, length=db))# raw longitude

        # repeat data
        da = repeat(da, repeat_num)
        da2 = repeat(da2, repeat_num)
        # shift the data
        one = da2[1+shift:end, :]
        two = da2[1:shift, :]
        da2 = vcat(one, two)
        intensity = repeat(intensity, repeat_num)
        intensity2 = repeat(intensity2, repeat_num)
        pulses = collect(1:length(intensity))
        le = length(pulses)
        ticks = [floor(Int, le /4), floor(Int, le /2), floor(Int, le *3 / 4)]
        fracs = [repeat_num / 4, repeat_num/2, repeat_num * 3 / 4]
        ti = ["$(fracs[i])\$P_3\$" for i in 1:3]

        tracks = []
        for i in 1:2
            push!(tracks, [])
            files = Glob.glob("track_*.jld2", "$(data_dir)/$i")  # no tracks? change here
            println(files)
            # load tracks
            for file in files
                @load file track
                tr = Tools.remove_duplicates(track)
                # add repeat_num data
                tr2 = [[], []]
                #println(size(tr))
                for j in 1:size(tr)[1]
                    x = tr[j,1]
                    if i == 1
                        y = tr[j,2]
                    else
                        y = tr[j,2] - shift
                    end
                    #println(num)
                    if (y >= 1)
                        push!(tr2[1], x)
                        push!(tr2[2], y)
                    end
                    for k in 1:repeat_num
                        if y + k*num < num * repeat_num
                            push!(tr2[1], x)
                            push!(tr2[2], y + k*num)
                        end
                    end
                end
                # sort
                sp = sortperm(tr2[2])
                tr2[1] = collect(view(tr2[1], sp))
                tr2[2] = collect(view(tr2[2], sp))
                push!(tracks[end], tr2)
            end
        end

        p2s = [[], []]

        for i in 1:2  # two sessions
            # two tracks? ok
            first = tracks[i][1]
            second = tracks[i][2]
            for (ii,y) in enumerate(first[2])
                for (jj,y2) in enumerate(second[2])
                    if y2-y == 0
                        dbin = abs(first[1][ii]-second[1][jj])
                        p2 = dbin / bins * 360
                        push!(p2s[i], p2)
                        #println("$i $y $p2")
                    end

                end
            end
        end
        println(mean(p2s[1]))
        println(std(p2s[1]))
        println(mean(p2s[2]))
        println(std(p2s[2]))

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 4.33071))  # 8cm x 11cm
        subplots_adjust(left=0.16, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)

        #left
        subplot2grid((5, 2), (0, 0), rowspan=4)
        minorticks_on()
        imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da), extent=[bin_st, bin_end, 1, num*repeat_num])
        tick_params(labelleft=true, labelbottom=false)
        yticks(ticks, ti)
        for tr in tracks[1]
            scatter(tr[1], tr[2], marker="o", c="none", ec="black", lw=0.3, s=2, alpha=0.7)
        end
        ylim([1, num*repeat_num])

        # right
        subplot2grid((5, 2), (0, 1), rowspan=4)
        minorticks_on()
        imshow(da2, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=0.7*darkness*maximum(da2), extent=[bin_st, bin_end, 1, num*repeat_num])
        tick_params(labelleft=false, labelbottom=false)
        yticks(ticks, ti)
        for tr in tracks[2]
            scatter(tr[1], tr[2], marker="o", c="none", ec="black", lw=0.3, s=2, alpha=0.7)
        end
        ylim([1, num*repeat_num])

        # bottom left
        subplot2grid((5, 2), (4, 0))
        minorticks_on()
        plot(longitude, average, c="grey")
        #xticks([-40, 0, 40])
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        ylabel("intensity (a. u.)")

        # bottom right
        subplot2grid((5, 2), (4, 1), colspan=3)
        minorticks_on()
        tick_params(labelleft=false)
        plot(longitude, average2, c="grey")
        #xticks([-40, 0, 40])
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])

        figtext(0.45, 0.01, "longitude (deg.)", size=8)

        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_p3fold_twotracks.pdf")
        savefig("$outdir/$(name_mod)_p3fold_twotracks.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
    end



    function offset(data, data2, outdir; start=1, number=nothing, repeat_num=4, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0")
        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1,bin_st:bin_end]
        da2 = data2[start:start+number-1,bin_st:bin_end]

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        # repeat data
        da = repeat(da, repeat_num)
        da2 = repeat(da2, repeat_num)
        #=
        le = length(pulses)
        ticks = [floor(Int, le /4), floor(Int, le /2), floor(Int, le *3 / 4)]
        fracs = [repeat_num / 4, repeat_num/2, repeat_num * 3 / 4]
        ti = ["$(fracs[i])\$P_3\$" for i in 1:3]
        =#

        #peak_data = Tools.find_peaks(da)
        peak_data = Tools.find_peaks3(da)
        #peak_data2 = Tools.find_peaks(da2)
        peak_data2 = Tools.find_peaks3(da2)
        #println("$val, $ind")
        #=
        for peak in peak_data
            for val in view(peak, 2:length(peak))
                scatter(val, peak[1], marker="x", c="black")
            end
        end

        for peak in peak_data2
            for val in view(peak, 2:length(peak))
                scatter(val, peak[1], marker="x", c="red")
            end
        end
        savefig("/home/szary/work/B0320/points.pdf")
        close()
        =#
        #return

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)


        figure(figsize=(3.14961, 4.33071))  # 8cm x 11cm
        #subplots_adjust(left=0.16, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)
        subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.99, wspace=0., hspace=0.)

        subplot2grid((5, 4), (0, 0), rowspan=4, colspan=2)
        minorticks_on()
        imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da), extent=[1, bin_end-bin_st, 1, num])
        ex = [1, bin_end-bin_st, 1, num]
        for peak in peak_data[2:end-1]
            for val in view(peak, 2:length(peak))
                println("$val $(peak[1])")
                scatter(val, peak[1], marker="x", c="green", s=2.5)
            end
        end
        grid(color="white", lw=0.4)
        grid(color="white", which="minor", lw=0.1, ls=":")
        tick_params(labelleft=false, labelbottom=false)
        #xlabel("122 MHz (BW: 9.77 MHz)")
        xlabel("1102 MHz (BW: MHz)")

        subplot2grid((5, 4), (0, 2), rowspan=4, colspan=2)
        minorticks_on()
        imshow(da2, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da2), extent=[1, bin_end-bin_st, 1, num])
        for peak in peak_data2[2:end-1]
            for val in view(peak, 2:length(peak))
                println("$val $(peak[1])")
                scatter(val, peak[1], marker="x", c="green", s=2.5)
            end
        end
        grid(color="white", lw=0.4)
        grid(color="white", which="minor", lw=0.1, ls=":")
        tick_params(labelleft=false, labelbottom=false)
        #xlabel("174 MHz (BW:29.5 MHz)")
        xlabel("1462 MHz (BW: MHz)")

        #tick_params(labeltop=false, labelbottom=true)
        savefig("$outdir/$(name_mod)_offset.pdf")
        close()
    end


    function offset_points(data, data2, outdir; start=1, number=nothing, repeat_num=3, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0")
        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1,bin_st:bin_end]
        da2 = data2[start:start+number-1,bin_st:bin_end]

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        # repeat data
        da = repeat(da, repeat_num)
        da2 = repeat(da2, repeat_num)

        #peak_data = Tools.find_peaks(da)
        peak_data = Tools.find_peaks3(da)
        #peak_data2 = Tools.find_peaks(da2)
        peak_data2 = Tools.find_peaks3(da2)

        off = []
        off_x = []
        for i in 1:size(peak_data)[1]
            for j in 2:size(peak_data[i])[1]
                mi = 1e50
                x1 = peak_data[i][j]
                y1 = peak_data[i][1]
                for k in 1:size(peak_data2)[1]
                    for l in 2:size(peak_data2[k])[1]
                        #println("$i $j $k $l")
                        x2 = peak_data2[k][l]
                        y2 = peak_data2[k][1]
                        of = sqrt((x2-x1)^2 + (y2-y1)^2)
                        if of < mi
                            mi = of
                        end
                    end
                end
                push!(off, mi)
                push!(off_x, x1)
            end
        end

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 4.33071))  # 8cm x 11cm
        #subplots_adjust(left=0.16, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)
        subplots_adjust(left=0.15, bottom=0.11, right=0.89, top=0.89, wspace=0., hspace=0.)

        subplot2grid((2, 1), (0, 0))
        minorticks_on()
        tick_params(labelleft=true, labelbottom=false, labeltop=true, which="both", top=true)
        for peak in peak_data
            for val in view(peak, 2:length(peak))
                scatter(val, peak[1], marker="x", c="black", s=2.5)
            end
        end
        for peak in peak_data2
            for val in view(peak, 2:length(peak))
                scatter(val, peak[1], marker="x", c="red", s=2.5)
            end
        end
        ylabel("y-bin (P_3)" )

        subplot2grid((2, 1), (1, 0))
        minorticks_on()
        scatter(off_x, off, marker="x", c="black", s=4)
        xlabel("x-bin (longitude)" )
        ylabel("offest (dbin)" )
        savefig("$outdir/$(name_mod)_offset_points.pdf")
        close()
    end


    function offset_points3(data, data2, data3, outdir; start=1, number=nothing, repeat_num=2, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0")
        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1,bin_st:bin_end]
        da2 = data2[start:start+number-1,bin_st:bin_end]
        da3 = data3[start:start+number-1,bin_st:bin_end]

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        # repeat data
        da = repeat(da, repeat_num)
        da2 = repeat(da2, repeat_num)
        da3 = repeat(da3, repeat_num)

        peak_data = Tools.find_peaks(da)
        peak_data2 = Tools.find_peaks(da2)
        peak_data3 = Tools.find_peaks(da3)

        off = []
        off_x = []
        off2 = []
        for i in 1:size(peak_data)[1]
            for j in 2:size(peak_data[i])[1]
                mi = 1e50
                x1 = peak_data[i][j]
                y1 = peak_data[i][1]
                for k in 1:size(peak_data2)[1]
                    for l in 2:size(peak_data2[k])[1]
                        #println("$i $j $k $l")
                        x2 = peak_data2[k][l]
                        y2 = peak_data2[k][1]
                        of = sqrt((x2-x1)^2 + (y2-y1)^2)
                        if of < mi
                            mi = of
                        end
                    end
                end
                mi2 = 1e50
                for k in 1:size(peak_data3)[1]
                    for l in 2:size(peak_data3[k])[1]
                        #println("$i $j $k $l")
                        x3 = peak_data3[k][l]
                        y3 = peak_data3[k][1]
                        of = sqrt((x3-x1)^2 + (y3-y1)^2)
                        if of < mi2
                            mi2 = of
                        end
                    end
                end
                push!(off, mi)
                push!(off2, mi2)
                push!(off_x, x1)
            end
        end

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 4.33071))  # 8cm x 11cm
        #subplots_adjust(left=0.16, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)
        subplots_adjust(left=0.15, bottom=0.11, right=0.89, top=0.89, wspace=0., hspace=0.)

        subplot2grid((2, 1), (0, 0))
        minorticks_on()
        tick_params(labelleft=true, labelbottom=false, labeltop=true, which="both", top=true)
        for peak in peak_data
            for val in view(peak, 2:length(peak))
                scatter(val, peak[1], marker="x", c="red", s=2.5)
            end
        end
        for peak in peak_data2
            for val in view(peak, 2:length(peak))
                scatter(val, peak[1], marker="x", c="green", s=2.5)
            end
        end
        for peak in peak_data3
            for val in view(peak, 2:length(peak))
                scatter(val, peak[1], marker="x", c="blue", s=2.5)
            end
        end
        ylabel("y-bin (P_3)" )

        subplot2grid((2, 1), (1, 0))
        minorticks_on()
        scatter(off_x, off, marker="x", c="green", s=4)
        scatter(off_x, off2, marker="x", c="blue", s=4)
        ylabel("offest (dbin)" )
        xlabel("x-bin (longitude)" )
        savefig("$outdir/offset_points3_$name_mod.pdf")
        close()
    end



    function crosscorplot(data, data2, outdir; start=1, number=nothing, repeat_num=2, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0")
        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1,bin_st:bin_end]
        da2 = data2[start:start+number-1,bin_st:bin_end]

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        # repeat data
        da = repeat(da, repeat_num)
        da2 = repeat(da2, repeat_num)

        cc = crosscor(da, da2, demean=true)

        println(size(cc))
        #return


        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 31*4.33071))  # 8cm x 11cm
        #subplots_adjust(left=0.16, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)
        subplots_adjust(left=0.15, bottom=0.11, right=0.89, top=0.89, wspace=0., hspace=0.)

        for i in 1:31
            subplot2grid((31, 1), (i-1, 0))
            minorticks_on()
            imshow(cc[i, :, :], origin="lower", cmap=cmap, interpolation="none", aspect="auto") #,  vmax=darkness*maximum(da))
        end
        #tick_params(labelleft=true, labelbottom=false, labeltop=true, which="both", top=true)

        savefig("$outdir/crosscorr_$name_mod.pdf")
        close()
    end



    function offset_subtract(data, data2, outdir; start=1, number=nothing, repeat_num=1, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0")
        num, bins = size(data)
        if number == nothing
            number = num - start + 1  # missing one? yes! +1
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1, bin_st:bin_end]
        da2 = data2[start:start+number-1, bin_st:bin_end]

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        da3 = Array{Float64}(undef, number, bin_end - bin_st + 1)
        for i in 1:number
            one = da[i,:]
            mi = minimum(one)
            ma = maximum(one)
            #println("$i $mi $ma")
            for j in 1:(bin_end - bin_st + 1)
                da[i,j] = (da[i,j]- mi) / (ma - mi)
            end

            two = da2[i,:]
            mi = minimum(two)
            ma = maximum(two)
            for j in 1:(bin_end - bin_st +1)
                da2[i,j] = (da2[i,j]- mi) / (ma - mi)
            end
            for j in 1:(bin_end - bin_st +1)
                da3[i,j] = da[i,j] - da2[i,j]
            end
        end
        # repeat data
        da = repeat(da, repeat_num)
        da2 = repeat(da2, repeat_num)
        da3 = repeat(da3, repeat_num)


        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(6.299213, 8.267717))  # 16cm x 21cm
        #subplots_adjust(left=0.16, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)
        subplots_adjust(left=0.01, bottom=0.05, right=0.99, top=0.99, wspace=0., hspace=0.)

        subplot2grid((1, 3), (0, 0))
        minorticks_on()
        imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da), extent=[1, bin_end-bin_st, 1, number*repeat_num])
        #xlabel("122 MHz (BW: 9.77 MHz)")
        xlabel("1102 MHz")

        subplot2grid((1, 3), (0, 1))
        minorticks_on()
        imshow(da2, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da2), extent=[1, bin_end-bin_st, 1, number*repeat_num])
        #xlabel("174 MHz (BW:29.5 MHz)")
        xlabel("1462 MHz")

        subplot2grid((1, 3), (0, 2))
        minorticks_on()
        #imshow(da3, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da3), extent=[1, bin_end-bin_st, 1, num])
        imshow(da3, origin="lower", cmap=cmap, aspect="auto",  vmax=darkness*maximum(da3),  extent=[1, bin_end-bin_st, 1, number*repeat_num])
        #xlabel("122 MHz minus 174 MHz")
        xlabel("1102 MHz minus 1462 MHz")

        savefig("$outdir/$(name_mod)_offset_subtract.pdf")
        show()
        readline(stdin; keep=false)
        close()
    end


    function subpulses(data, outdir, peaks; start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME")
        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1, bin_st:bin_end]
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
        subplots_adjust(left=0.16, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)

        subplot2grid((5, 3), (0, 0), rowspan=4)
        minorticks_on()
        plot(intensity, pulses, c="grey")
        ylim(pulses[1]-0.5, pulses[end]+0.5)
        xticks([0.5, 1.0])
        xlim(1.1, -0.1)
        xlabel("intensity")
        ylabel("Pulse number")
        extent = [bin_st, bin_end, start-0.5, start+number-1+0.5]
        subplot2grid((5, 3), (0, 1), rowspan=4, colspan=2)
        imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da), extent=extent)
        xlim([extent[1], extent[2]])
        ylim([extent[3], extent[4]])
        #println([bin_st, bin_end, start-0.5, start+number+0.5])
        for peak in peaks
            if (peak[1] >= start) && (peak[1] <= start + number)
                for x in peak[2]
                    scatter(x, peak[1], marker="x", c="red", s=12.5, lw=1)
                    #println("$(peak[1]) $x")
                end
            end
        end
        tick_params(labelleft=false, labelbottom=false)

        subplot2grid((5, 3), (4, 1), colspan=2)
        minorticks_on()
        plot(longitude, average, c="grey")
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        xlabel("longitude \$(^\\circ)\$")
        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_subpulses.pdf")
        savefig("$outdir/$(name_mod)_subpulses.pdf")
        #show()
        st = readline(stdin; keep=false)
        close()
        #clf()
    end


    function tracks(data, outdir, tracks, peaks; start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME")
        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1, bin_st:bin_end]
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
        subplots_adjust(left=0.16, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)

        subplot2grid((5, 3), (0, 0), rowspan=4)
        minorticks_on()
        plot(intensity, pulses, c="grey")
        ylim(pulses[1]-0.5, pulses[end]+0.5)
        xticks([0.5, 1.0])
        xlim(1.1, -0.1)
        xlabel("intensity")
        ylabel("Pulse number")
        extent = [bin_st, bin_end, start-0.5, start+number-1+0.5]
        subplot2grid((5, 3), (0, 1), rowspan=4, colspan=2)
        imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da), extent=extent)
        xlim([extent[1], extent[2]])
        ylim([extent[3], extent[4]])
        #println([bin_st, bin_end, start-0.5, start+number+0.5])
        println(length(tracks))
        for track in tracks
            if (track[1][1] >= start) && (track[1][1] <= start + number)
                plot(track[2], track[1], lw=3)
            end
        end
        for peak in peaks
            if (peak[1] >= start) && (peak[1] <= start + number)
                for x in peak[2]
                    scatter(x, peak[1], marker="x", c="red", s=9.5, lw=1)
                    #println("$(peak[1]) $x")
                end
            end
        end
        tick_params(labelleft=false, labelbottom=false)

        subplot2grid((5, 3), (4, 1), colspan=2)
        minorticks_on()
        plot(longitude, average, c="grey")
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        xlabel("longitude \$(^\\circ)\$")
        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_tracks.pdf")
        savefig("$outdir/$(name_mod)_tracks.pdf")
        #show()
        st = readline(stdin; keep=false)
        close()
        #clf()
    end

    "start : start character ext: file type"
    function getlastnumber(files, start="_", ext=".jld2")
        numbers = []
        for file in files
            num = parse(Int, replace(split(split(file, "/")[end], start)[end], ext=>"")) # wow
            push!(numbers, num)
        end
        nrs = sort(numbers)
        return nrs[end]
    end


    function group_tracks(data, outdir, peaks; start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME")
        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1, bin_st:bin_end]
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

        picks = []
        add = true
        function onkey(event)
            if event.key == "a"
                if add == true
                    add = false
                else
                    add = true
                end
                println("Adding points: ", add)
            elseif event.key == "w"
                files = Glob.glob("track_*.jld2", outdir)
                files = sort(files)
                if length(files) > 0
                   lastindex = getlastnumber(files)
                else
                    lastindex = 0
                end
                filename = joinpath(outdir, "track_$(lastindex+1).jld2")
                if length(picks) > 0
                    # sort and convert to 2D
                    picks = sort(picks, by=x->x[2])  # nice and easy
                    sz = size(picks)
                    track = Array{Float64}(undef, sz[1], 2)
                    for i in 1:sz[1]
                        track[i, 1] = picks[i][1] # IMPORTANT do not use convert_tracks anymore! extent changed
                        track[i, 2] = picks[i][2]
                    end
                    #println(track[:,1])
                    @save filename track
                    # @load "example.jld2" hello foo
                    println("File $filename saved.")
                    #println(picks)
                    picks = []
                else
                    println("Nothing to save")
                end
            elseif event.key == "d"
                println("Droping all picks...")
                picks = []
            end

        end

        function onclick(event)
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            for i in ind
                pi = (xdata[i+1], ydata[i+1])
                if add == true
                    push!(picks, pi)
                    thisline.set_color("green")
                    thisline.set_alpha(0.0) # transparent!
                else
                    #pop!(picks, pi)  # does not work
                    filter!(p->p!=pi, picks) # works, but I still not get it
                    thisline.set_color("red")
                    thisline.set_alpha(0.3)
                end
            end
            #println(picks)
        end


        tracks = []
        # load tracks
        files = Glob.glob("track_*.jld2", outdir)
        for file in files
            @load file track
            push!(tracks, track)
        end


        f = figure(figsize=(3.14961, 4.33071), frameon=true)  # 8cm x 11cm
        f.canvas.mpl_connect("pick_event", onclick)
        f.canvas.mpl_connect("key_press_event", onkey)
        subplots_adjust(left=0.16, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)

        subplot2grid((5, 3), (0, 0), rowspan=4)
        minorticks_on()
        plot(intensity, pulses, c="grey")
        ylim(pulses[1]-0.5, pulses[end]+0.5)
        xticks([0.5, 1.0])
        xlim(1.1, -0.1)
        xlabel("intensity")
        ylabel("Pulse number")
        extent = [bin_st/bins*360, bin_end/bins*360, start-0.5, start+number-1+0.5]
        subplot2grid((5, 3), (0, 1), rowspan=4, colspan=2)
        #imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da), extent=extent)
        imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(da), vmin=0, extent=extent, alpha=1.0)
        xlim([extent[1], extent[2]])
        ylim([extent[3], extent[4]])
        #println([bin_st, bin_end, start-0.5, start+number+0.5])
        for peak in peaks
            if (peak[1] >= start) && (peak[1] <= start + number)
                for x in peak[2]
                    plot(x/bins*360, peak[1], marker="x", markersize=3, markeredgewidth=2.2, c="red", fillstyle="full", mfc="red", lw=0, picker=10, alpha=0.5)
                    #scatter(x, peak[1], marker="x", c="red", s=9.5, lw=1, picker=2)
                    #println("$(peak[1]) $x")
                end
            end
        end
        for track in tracks
            plot(track[:,1], track[:,2], marker="x", markersize=2.5, c="green", lw=0, alpha=1.0)
        end
        tick_params(labelleft=false, labelbottom=false)

        subplot2grid((5, 3), (4, 1), colspan=2)
        minorticks_on()
        plot(longitude, average, c="grey")
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        xlabel("longitude \$(^\\circ)\$")
        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_tracks_grouping.pdf")
        savefig("$outdir/$(name_mod)_tracks_grouping.pdf")
        #show()
        st = readline(stdin; keep=false)
        close()
        #clf()
    end


    function show_tracks(data, outdir, peaks; start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME")
        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1, bin_st:bin_end]
        average = Tools.average_profile(da)
        intensity, pulses = Tools.intensity_pulses(da)
        intensity .-= minimum(intensity)
        intensity ./= maximum(intensity)

        pulses .+= start - 1  # julia

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        tracks = []
        # load tracks
        files = Glob.glob("track_*.jld2", outdir)
        for file in files
            @load file track
            push!(tracks, track)
        end

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        f = figure(figsize=(3.14961, 4.33071), frameon=true)  # 8cm x 11cm
        subplots_adjust(left=0.19, bottom=0.09, right=0.99, top=0.99, wspace=0., hspace=0.)

        subplot2grid((5, 3), (0, 0), rowspan=4)
        minorticks_on()
        plot(intensity, pulses, c="grey")
        ylim(pulses[1]-0.5, pulses[end]+0.5)
        xticks([0.5, 1.0])
        xlim(1.1, -0.1)
        xlabel("intensity")
        ylabel("Pulse number")
        extent = [bin_st, bin_end, start-0.5, start+number-1+0.5]
        subplot2grid((5, 3), (0, 1), rowspan=4, colspan=2)
        imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da), extent=extent)
        xlim([extent[1], extent[2]])
        ylim([extent[3], extent[4]])
        for track in tracks
            plot(track[:,1], track[:,2], marker="x", markersize=1.5, c="red", lw=0, alpha=1.5)
        end
        #=
        for peak in peaks
            if (peak[1] >= start) && (peak[1] <= start + number)
                for x in peak[2]
                    plot(x, peak[1], marker="x", markersize=2.5, c="red", lw=1, picker=10)
                end
            end
        end
        =#
        tick_params(labelleft=false, labelbottom=false)

        subplot2grid((5, 3), (4, 1), colspan=2)
        minorticks_on()
        plot(longitude, average, c="grey")
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        xlabel("longitude \$(^\\circ)\$")
        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_tracks.pdf")
        savefig("$outdir/$(name_mod)_tracks.pdf")
        #show()
        st = readline(stdin; keep=false)
        close()
        #clf()
    end


    function tracks_analysis(outdir; win=100, start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME")

        tracks = []
        # load tracks
        files = Glob.glob("track_*.jld2", outdir)
        for file in files
            @load file track
            push!(tracks, track)
        end
        lines = []
        inclines = []
        for track in tracks#[1:2]
            ll, inc = Tools.analyse_track(track[:,2], track[:,1]; win=win)
            push!(lines, ll)
            push!(inclines, inc)
        end

        spl = fit(SmoothingSpline, tracks[1][:,2], tracks[1][:,1], 12500.0)
        ysp = SmoothingSplines.predict(spl)

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961*2, 2.362205*2), frameon=true)  # 16cm x 12 cm
        subplots_adjust(left=0.13, bottom=0.08, right=0.99, top=0.90, wspace=0., hspace=0.)

        ax = subplot2grid((2, 1), (0, 0))
        ax.xaxis.set_label_position("top")
        tick_params(labeltop=true, labelbottom=false, which="both", bottom=false, top=true)
        minorticks_on()
        xlabel("Pulse number")
        ylabel("Longitude \$(^\\circ)\$")
        for track in tracks#[1:2]
            plot(track[:,2], track[:,1]) #, marker="x", markersize=2.5, lw=1)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=2.5, lw=0)
        end
        #plot(tracks[1][:,2], ysp, lw=2, alpha=0.4, marker="x")
        for line in lines
            for ll in line
                plot(ll[1],  ll[2], lw=0.3)
            end
        end
        xl = xlim()
        subplot2grid((2, 1), (1, 0))
        minorticks_on()
        ylabel("Drift rate \$(^\\circ / P)\$")
        xlabel("Pulse number")
        for i in 1:length(inclines)
            for j in 1:length(inclines[i])
                errorbar(inclines[i][j][1], inclines[i][j][2], yerr=inclines[i][j][3], color="none", lw=1., marker="_", mec="grey", ecolor="grey", capsize=0, mfc="grey", ms=1.0)
            end
            #scatter(inclines[i][1], inclines[i][2])
        end
        axhline(y=0, lw=1, ls="--")
        xlim(xl)
        println("$outdir/$(name_mod)_tracks_analysis.pdf")
        savefig("$outdir/$(name_mod)_tracks_analysis.pdf")
        show()
        st = readline(stdin; keep=false)
        close()
        #clf()
    end


    function tracks_analysis2(outdir; win=100, start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME")

        tracks = []
        # load tracks
        files = Glob.glob("track_*.jld2", outdir)
        for file in files
            @load file track
            push!(tracks, track)
        end
        lines = []
        inclines = []
        drift_rate = [[], []]
        for track in tracks#[1:3]
            ll, inc = Tools.analyse_track(track[:,2], track[:,1]; win=win)
            push!(lines, ll)
            push!(inclines, inc)
        end
        for inc in inclines
            #println(length(inc))
            for pu in inc
                push!(drift_rate[1], pu[1])
                push!(drift_rate[2], pu[2])
            end
        end
        sp = sortperm(drift_rate[1])
        drift_rate[1] = view(drift_rate[1], sp)
        drift_rate[2] = view(drift_rate[2], sp)
        dr1 = convert(Array{Float64,1}, drift_rate[1])
        dr2 = convert(Array{Float64,1}, drift_rate[2])

        spl = fit(SmoothingSpline, dr1, dr2, 100.0)
        #spl = fit(SmoothingSpline, tracks[1][:,2], tracks[1][:,1], 12500.0)
        ysp = SmoothingSplines.predict(spl)
        # add fft here
        half = floor(Int, length(ysp) / 2) # one side frequency range?
        ff = fft(ysp)[1:half]
        freq = Tools.fftfreq(length(ysp))

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961*2, 2.362205*2), frameon=true)  # 16cm x 12 cm
        subplots_adjust(left=0.13, bottom=0.08, right=0.99, top=0.90, wspace=0., hspace=0.)

        ax = subplot2grid((3, 1), (0, 0))
        ax.xaxis.set_label_position("top")
        tick_params(labeltop=true, labelbottom=false, which="both", bottom=false, top=true)
        minorticks_on()
        xlabel("Pulse number")
        ylabel("Longitude \$(^\\circ)\$")
        for track in tracks#[1:3]
            plot(track[:,2], track[:,1]) #, marker="x", markersize=2.5, lw=1)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=2.5, lw=0)
        end
        for line in lines
            for ll in line
                plot(ll[1],  ll[2], lw=0.3)
            end
        end
        xl = xlim()
        subplot2grid((3, 1), (1, 0))
        tick_params(labeltop=true, labelbottom=true, which="both", bottom=true, top=true)
        ylabel("Drift rate \$(^\\circ / P)\$")
        minorticks_on()
        for i in 1:length(inclines)
            for j in 1:length(inclines[i])
                errorbar(inclines[i][j][1], inclines[i][j][2], yerr=inclines[i][j][3], color="none", lw=0.5, marker="_", mec="black", ecolor="black", capsize=0, mfc="grey", ms=0.5)
            end
        end
        axhline(y=0, lw=1, ls="--")
        plot(dr1, ysp, lw=4, alpha=0.5, c="grey")
        xlim(xl)
        subplot2grid((3, 1), (2, 0))
        minorticks_on()
        #xlabel("Pulse number")
        #xlim(xl)
        #plot(dr1, dr2, lw=0.3)
        #plot(dr1, ysp, lw=2, alpha=0.7)
        plot(freq, ff)

        println("$outdir/$(name_mod)_tracks_analysis2.pdf")
        savefig("$outdir/$(name_mod)_tracks_analysis2.pdf")
        show()
        st = readline(stdin; keep=false)
        close()
        #clf()
    end


    function tracks_analysis3(outdir; lambda=100.0, start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME")

        tracks = []
        # load tracks
        files = Glob.glob("track_*.jld2", outdir)
        for file in files
            @load file track
            tr = Tools.remove_duplicates(track)
            push!(tracks, tr)
        end

        lines = []
        inclines = []
        drift_rate = [[], []]
        for track in tracks#[1:3]
            #ll, inc = Tools.analyse_track_sl(track[:,3], track[:,1]; lambda=lambda)
            ll, inc = Tools.analyse_track_simple(track[:,2], track[:,1]; lambda=lambda)
            push!(lines, ll)
            push!(inclines, inc)
        end
        for inc in inclines
            for i in 1:length(inc[1])
                push!(drift_rate[1], inc[1][i])
                push!(drift_rate[2], inc[2][i])
            end
        end
        sp = sortperm(drift_rate[1])
        drift_rate[1] = view(drift_rate[1], sp)
        drift_rate[2] = view(drift_rate[2], sp)
        dr1 = convert(Array{Float64,1}, drift_rate[1])
        dr2 = convert(Array{Float64,1}, drift_rate[2])

        spl = fit(SmoothingSpline, dr1, dr2, 100.0)
        ysp = SmoothingSplines.predict(spl)
        # add fft here
        half = floor(Int, length(ysp) / 2) # one side frequency range?
        ff = fft(ysp)[1:half]
        freq = Tools.fftfreq(length(ysp))

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961*2, 2.362205*2), frameon=true)  # 16cm x 12 cm
        subplots_adjust(left=0.13, bottom=0.08, right=0.99, top=0.90, wspace=0., hspace=0.)

        ax = subplot2grid((3, 1), (0, 0))
        ax.xaxis.set_label_position("top")
        tick_params(labeltop=true, labelbottom=false, which="both", bottom=false, top=true)
        minorticks_on()
        xlabel("Pulse number")
        ylabel("Longitude \$(^\\circ)\$")
        for track in tracks#[1:3]
            #plot(track[:,2], track[:,1]) #, marker="x", markersize=2.5, lw=1)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=2.5, lw=0)
        end
        for line in lines
            plot(line[1], line[2], lw=2, alpha=0.5)
        end
        xl = xlim()

        subplot2grid((3, 1), (1, 0))
        tick_params(labeltop=true, labelbottom=true, which="both", bottom=true, top=true)
        ylabel("Drift rate \$(^\\circ / P)\$")
        minorticks_on()
        for inc in inclines
            plot(inc[1], inc[2], lw=1, c="black")
        end
        axhline(y=0, lw=1, ls="--")
        plot(dr1, ysp, lw=4, alpha=0.5, c="grey")
        xlim(xl)

        subplot2grid((3, 1), (2, 0))
        minorticks_on()
        #xlabel("Pulse number")
        #xlim(xl)
        #plot(dr1, dr2, lw=0.3)
        #plot(dr1, ysp, lw=2, alpha=0.7)
        plot(freq, ff, lw=1)

        println("$outdir/$(name_mod)_tracks_analysis3.pdf")
        savefig("$outdir/$(name_mod)_tracks_analysis3.pdf")
        show()
        st = readline(stdin; keep=false)
        close()
        #clf()
    end


    function tracks_analysis4(outdir; start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME")

        tracks = []
        # load tracks
        files = Glob.glob("track_*.jld2", outdir)
        for file in files
            @load file track
            tr = Tools.remove_duplicates(track)
            push!(tracks, tr)
        end

        lines = []
        inclines = []
        drift_rate = [[], []]
        for track in tracks#[1:3]
            #ll, inc = Tools.analyse_track_sl(track[:,3], track[:,1]; lambda=lambda)
            ll, inc = PyRModule.cubic_spline(track[:,2], track[:,1])
            push!(lines, ll)
            push!(inclines, inc)
        end
        for inc in inclines
            for i in 1:length(inc[1])
                push!(drift_rate[1], inc[1][i])
                push!(drift_rate[2], inc[2][i])
            end
        end
        sp = sortperm(drift_rate[1])
        drift_rate[1] = view(drift_rate[1], sp)
        drift_rate[2] = view(drift_rate[2], sp)
        dr1 = convert(Array{Float64,1}, drift_rate[1])
        dr2 = convert(Array{Float64,1}, drift_rate[2])

        spl = fit(SmoothingSpline, dr1, dr2, 100.0)
        ysp = SmoothingSplines.predict(spl)
        # add fft here
        half = floor(Int, length(ysp) / 2) # one side frequency range?
        ff = fft(ysp)[1:half]
        freq = Tools.fftfreq(length(ysp))

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961*2, 2.362205*2), frameon=true)  # 16cm x 12 cm
        subplots_adjust(left=0.13, bottom=0.08, right=0.99, top=0.90, wspace=0., hspace=0.)

        ax = subplot2grid((3, 1), (0, 0))
        ax.xaxis.set_label_position("top")
        tick_params(labeltop=true, labelbottom=false, which="both", bottom=false, top=true)
        minorticks_on()
        xlabel("Pulse number")
        ylabel("Longitude \$(^\\circ)\$")
        for track in tracks#[1:3]
            #plot(track[:,2], track[:,1]) #, marker="x", markersize=2.5, lw=1)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=2.5, lw=0)
        end
        for line in lines
            plot(line[1], line[2], lw=2, alpha=0.5)
        end
        xl = xlim()

        subplot2grid((3, 1), (1, 0))
        tick_params(labeltop=true, labelbottom=true, which="both", bottom=true, top=true)
        ylabel("Drift rate \$(^\\circ / P)\$")
        minorticks_on()
        for inc in inclines
            plot(inc[1], inc[2], lw=1, c="black")
        end
        axhline(y=0, lw=1, ls="--")
        plot(dr1, ysp, lw=4, alpha=0.5, c="grey")
        xlim(xl)

        subplot2grid((3, 1), (2, 0))
        minorticks_on()
        #xlabel("Pulse number")
        #xlim(xl)
        #plot(dr1, dr2, lw=0.3)
        #plot(dr1, ysp, lw=2, alpha=0.7)
        plot(freq, ff, lw=1)

        println("$outdir/$(name_mod)_tracks_analysis4.pdf")
        savefig("$outdir/$(name_mod)_tracks_analysis4.pdf")
        show()
        st = readline(stdin; keep=false)
        close()
        #clf()
    end


    function driftrate_J1750(outdir; lambda=100, name_mod="123456", show_=false)

        tracks1, lines1, inclines1, dr1, ysp1, dof1 = Tools.get_driftrate("$outdir/tracks1", lambda)
        tracks2, lines2, inclines2, dr2, ysp2, dof2 = Tools.get_driftrate("$outdir/tracks2", lambda)
        tracks3, lines3, inclines3, dr3, ysp3, dof3 = Tools.get_driftrate("$outdir/tracks3", lambda)
        tracks4, lines4, inclines4, dr4, ysp4, dof4 = Tools.get_driftrate("$outdir/tracks4", lambda)
        tracks5, lines5, inclines5, dr5, ysp5, dof5 = Tools.get_driftrate("$outdir/tracks5", lambda)
        tracks6, lines6, inclines6, dr6, ysp6, dof6 = Tools.get_driftrate("$outdir/tracks6", lambda)

        # check all tracks in second session
        for i in 1:length(tracks2)
            y_ = tracks2[i][:,1]
            y_fit = lines2[i][2]
            chisq = Tools.chisquare_reduced(y_, y_fit, dof2[i])
            rsq = Tools.rsquared(y_, y_fit)
            println("Track: $i Chi^2_ν: ", chisq, " R^2: ", rsq)
        end

        y0 = 401 # bins!
        y1 = 599

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(7.086614, 6.299213))  # 18 cm x 16 cm
        subplots_adjust(left=0.08, bottom=0.07, right=0.95, top=0.92, wspace=0.0, hspace=0.0)
        figtext(0.09, 0.9, "a)", size=10)
        figtext(0.239, 0.9, "b)", size=10)
        figtext(0.74, 0.9, "c)", size=10)
        figtext(0.09, 0.44, "d)", size=10)
        figtext(0.43, 0.44, "e)", size=10)
        figtext(0.74, 0.44, "f)", size=10)

        # first session
        ax = subplot2grid((34, 1800), (0, 0), colspan=294, rowspan=8)  # row column
        minorticks_on()
        for track in tracks1#[1:3]
            #plot(track[:,2], track[:,1]) #, marker="x", markersize=2.5, lw=1)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=2.5, lw=0)
        end
        for line in lines1
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        xl = xlim()
        ylim([y0, y1])
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        ylabel("Longitude (\$ ^{\\circ}\$)") # TODO no! bins!
        # first session
        subplot2grid((34, 1800), (8, 0), colspan=294, rowspan=8)  # row column
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr1, ysp1, lw=3, alpha=1.0, c="grey")
        for inc in inclines1
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        xlim(xl)
        ylabel("Drift rate \$(^\\circ / P)\$")
        ylim([-2.3, 2.3])
        tick_params(axis="x", which="both", direction="out", labelbottom=false)

        # second session
        ax = subplot2grid((34, 1800), (0, 310), colspan=1030, rowspan=8)  # row column
        minorticks_on()
        for track in tracks2#[1:3]
            #plot(track[:,2], track[:,1]) #, marker="x", markersize=2.5, lw=1)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=2.5, lw=0)
        end
        for line in lines2
            plot(line[1], line[2], lw=2, alpha=0.5, c="C2")
        end
        xl = xlim()
        ylim([y0, y1])
        xlabel("Pulse number")
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=true)
        # second session
        subplot2grid((34, 1800), (8, 310), colspan=1030, rowspan=8)  # row column
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr2, ysp2, lw=4, alpha=1.0, c="grey")
        for inc in inclines2
            plot(inc[1], inc[2], lw=0.5, alpha=0.5, c="black")
        end
        xlim(xl)
        ylim([-2.3, 2.3])
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=true)
        tick_params(axis="x", which="both", direction="out", labelbottom=false)

        # third session
        ax = subplot2grid((34, 1800), (0, 1355), colspan=445, rowspan=8)  # row column
        minorticks_on()
        for track in tracks3#[1:3]
            #plot(track[:,2], track[:,1]) #, marker="x", markersize=2.5, lw=1)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=2.5, lw=0)
        end
        for line in lines3
            plot(line[1], line[2], lw=2, alpha=0.5, c="C2")
        end
        xl = xlim()
        ylim([y0, y1])
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        ax.yaxis.set_label_position("right")
        ax.yaxis.set_ticks_position("right")
        # third session
        ax = subplot2grid((34, 1800), (8, 1355), colspan=445, rowspan=8)  # row column
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr3, ysp3, lw=4, alpha=1.0, c="grey")
        for inc in inclines3
            plot(inc[1], inc[2], lw=0.5, alpha=0.5, c="black")
        end
        xlim(xl)
        ylim([-2.3, 2.3])
        ax.yaxis.set_label_position("right")
        ax.yaxis.set_ticks_position("right")
        tick_params(axis="x", which="both", direction="out", labelbottom=false)

        # fourth session
        ax = subplot2grid((34, 1800), (18, 0), colspan=449, rowspan=8)  # row column
        minorticks_on()
        for track in tracks4#[1:3]
            #plot(track[:,2], track[:,1]) #, marker="x", markersize=2.5, lw=1)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=2.5, lw=0)
        end
        for line in lines4
            plot(line[1], line[2], lw=2, alpha=0.5, c="C2")
        end
        xl = xlim()
        ylim([y0, y1])
        ylabel("Longitude (\$ ^{\\circ}\$)")
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        # fourth session
        subplot2grid((34, 1800), (26, 0), colspan=449, rowspan=8)  # row column
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr4, ysp4, lw=4, alpha=1.0, c="grey")
        for inc in inclines4
            plot(inc[1], inc[2], lw=0.5, alpha=0.5, c="black")
        end
        ylabel("Drift rate \$(^\\circ / P)\$")
        xlim(xl)
        ylim([-2.3, 2.3])

        # fifth session
        ax = subplot2grid((34, 1800), (18, 700), colspan=446, rowspan=8)  # row column
        minorticks_on()
        for track in tracks5#[1:3]
            #plot(track[:,2], track[:,1]) #, marker="x", markersize=2.5, lw=1)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=2.5, lw=0)
        end
        for line in lines5
            plot(line[1], line[2], lw=2, alpha=0.5, c="C2")
        end
        xl = xlim()
        ylim([y0, y1])
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        # fifth session
        subplot2grid((34, 1800), (26, 700), colspan=446, rowspan=8)  # row column
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr5, ysp5, lw=4, alpha=1.0, c="grey")
        for inc in inclines5
            plot(inc[1], inc[2], lw=0.5, alpha=0.5, c="black")
        end
        xlim(xl)
        ylim([-2.3, 2.3])
        xlabel("Pulse number")

        # sixth session
        ax = subplot2grid((34, 1800), (18, 1350), colspan=446, rowspan=8)  # row column
        minorticks_on()
        for track in tracks6#[1:3]
            #plot(track[:,2], track[:,1]) #, marker="x", markersize=2.5, lw=1)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=2.5, lw=0)
        end
        for line in lines6
            plot(line[1], line[2], lw=2, alpha=0.5, c="C2")
        end
        xl = xlim()
        ylim([y0, y1])
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        # sixth session
        subplot2grid((34, 1800), (26, 1350), colspan=446, rowspan=8)  # row column
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr6, ysp6, lw=4, alpha=1.0, c="grey")
        for inc in inclines6
            plot(inc[1], inc[2], lw=0.5, alpha=0.5, c="black")
        end
        xlim(xl)
        ylim([-2.3, 2.3])

        println("$outdir/$(name_mod)_driftrate.pdf")
        savefig("$outdir/$(name_mod)_driftrate.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
        #clf()
    end

    function get_length(tracks)
        mi = 1e50
        ma = -1e50
        for track in tracks
            ex = extrema(track[:,2])
            if ex[1] < mi
                mi = ex[1]
            end
            if ex[2] > ma
                ma = ex[2]
            end
        end
        return ma - mi + 1 # +1?!
    end


    function driftrate_J1750_2(outdir; spar=0.6, lambda=1000, name_mod="1234567", show_=false)

        #=
        tracks1, lines1, inclines1, dr1, ysp1, dof1 = Tools.get_driftrate2("$outdir/tracks/1", spar)
        tracks2, lines2, inclines2, dr2, ysp2, dof2 = Tools.get_driftrate2("$outdir/tracks/2", spar)
        tracks3, lines3, inclines3, dr3, ysp3, dof3 = Tools.get_driftrate2("$outdir/tracks/3", spar)
        tracks4, lines4, inclines4, dr4, ysp4, dof4 = Tools.get_driftrate2("$outdir/tracks/4", spar)
        tracks5, lines5, inclines5, dr5, ysp5, dof5 = Tools.get_driftrate2("$outdir/tracks/5", spar)
        tracks6, lines6, inclines6, dr6, ysp6, dof6 = Tools.get_driftrate2("$outdir/tracks/6", spar)
        tracks7, lines7, inclines7, dr7, ysp7, dof7 = Tools.get_driftrate2("$outdir/tracks/7", spar)
        =#

        tracks1, lines1, inclines1, dr1, ysp1, dof1 = Tools.get_driftrate("$outdir/tracks/1", lambda)
        tracks2, lines2, inclines2, dr2, ysp2, dof2 = Tools.get_driftrate("$outdir/tracks/2", lambda)
        tracks3, lines3, inclines3, dr3, ysp3, dof3 = Tools.get_driftrate("$outdir/tracks/3", lambda)
        tracks4, lines4, inclines4, dr4, ysp4, dof4 = Tools.get_driftrate("$outdir/tracks/4", lambda)
        tracks5, lines5, inclines5, dr5, ysp5, dof5 = Tools.get_driftrate("$outdir/tracks/5", lambda)
        tracks6, lines6, inclines6, dr6, ysp6, dof6 = Tools.get_driftrate("$outdir/tracks/6", lambda)
        tracks7, lines7, inclines7, dr7, ysp7, dof7 = Tools.get_driftrate("$outdir/tracks/7", lambda)

        # check all tracks in second session
        for i in 1:length(tracks2)
            y_ = tracks2[i][:,1]
            y_fit = lines2[i][2]
            chisq = Tools.chisquare_reduced(y_, y_fit, dof2[i])
            rsq = Tools.rsquared(y_, y_fit)
            println("Track: $i Chi^2_ν: ", chisq, " R^2: ", rsq)
        end

        #=
        println(get_length(tracks1))
        println(get_length(tracks2))
        println(get_length(tracks3))
        println(get_length(tracks4))
        println(get_length(tracks5))
        println(get_length(tracks6))
        println(get_length(tracks7))
        =#

        a = get_length(tracks1) + get_length(tracks2) + get_length(tracks3) # 1900 is fine
        b = get_length(tracks4) + get_length(tracks5) + get_length(tracks6) + get_length(tracks7) # 1900 is fine
        println("$a $((1900-a)/2) ")
        println("$b $((1900-b)/3)")

        #yl = (141, 210)
        yl = (210, 131)
        #yl = (135, 215)
        #yl = (125, 225)
        yl2 = (-0.9, 0.9)

        # CHECKING periodicity here
        #ys_ = convert(Array{Float64,1}, ysp2)
        #Tools.periodicity(ys_)
        #return

        rc("font", size=7.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        #figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
        #figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
        figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
        subplots_adjust(left=0.08, bottom=0.07, right=0.95, top=0.92, wspace=0.0, hspace=0.0)
        figtext(0.09, 0.89, "a)", size=10)
        figtext(0.25, 0.89, "b)", size=10)
        figtext(0.75, 0.89, "c)", size=10)
        figtext(0.09, 0.43, "d)", size=10)
        figtext(0.31, 0.43, "e)", size=10)
        figtext(0.54, 0.43, "f)", size=10)
        figtext(0.76, 0.43, "g)", size=10)

        # first session
        ax = subplot2grid((43, 1900), (0, 0), colspan=295, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        minorticks_on()
        ylabel("Longitude (\$ ^{\\circ}\$)")
        for track in tracks1#[1:3]
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.0, lw=0, markeredgewidth=0.2)
        end
        for line in lines1
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        xl = xlim()
        ylim(yl)
        # first session
        subplot2grid((43, 1900), (10, 0), colspan=295, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr1, ysp1, lw=3, alpha=1.0, c="grey")
        for inc in inclines1
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        xlim(xl)
        ylabel("Drift rate \$(^\\circ / P)\$")
        ylim(yl2)

        # second session
        ax = subplot2grid((43, 1900), (0, 359), colspan=1031, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        xlabel("Pulse number")
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        for track in tracks2
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.5, lw=0, markeredgewidth=0.3)
        end
        for line in lines2
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        xl = xlim()
        ylim(yl)
        # second session
        subplot2grid((43, 1900), (10, 359), colspan=1031, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr2, ysp2, lw=3, alpha=1.0, c="grey")
        for inc in inclines2
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        xlim(xl)
        ylim(yl2)
        # save driftrate for modeling
        #=
        f = open("$outdir/$(name_mod)_driftrate.txt", "w")
        for i in 1:length(dr2)
            write(f, "$(dr2[i]) $(ysp2[i])\n")
        end
        close(f)
        =#

        # third session
        ax = subplot2grid((43, 1900), (0, 1454), colspan=446, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true, labelright=true)
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        #ax.yaxis.set_label_position("right")
        #ax.yaxis.set_ticks_position("right")
        minorticks_on()
        for track in tracks3
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.5, lw=0, markeredgewidth=0.3)
        end
        for line in lines3
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        xl = xlim()
        ylim(yl)
        # third session
        ax = subplot2grid((43, 1900), (10, 1454), colspan=446, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true, labelright=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr3, ysp3, lw=3, alpha=1.0, c="grey")
        for inc in inclines3
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        xlim(xl)
        ylim(yl2)
        #ax.yaxis.set_label_position("right")
        #ax.yaxis.set_ticks_position("right")

        # fourth session
        ax = subplot2grid((43, 1900), (23, 0), colspan=450, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        minorticks_on()
        for track in tracks4
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.5, lw=0, markeredgewidth=0.3)
        end
        for line in lines4
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        xl = xlim()
        ylim(yl)
        ylabel("Longitude (\$ ^{\\circ}\$)")
        # fourth session
        subplot2grid((43, 1900), (33, 0), colspan=450, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr4, ysp4, lw=3, alpha=1.0, c="grey")
        for inc in inclines4
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        xlim(xl)
        ylim(yl2)
        ylabel("Drift rate \$(^\\circ / P)\$")
        locator_params(nbins=3)

        # fifth session
        ax = subplot2grid((43, 1900), (23, 487), colspan=447, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        for track in tracks5
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.5, lw=0, markeredgewidth=0.3)
        end
        for line in lines5
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        xl = xlim()
        ylim(yl)
        # fifth session
        subplot2grid((43, 1900), (33, 487), colspan=447, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr5, ysp5, lw=3, alpha=1.0, c="grey")
        for inc in inclines5
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        xlim(xl)
        ylim(yl2)
        #xlabel("Pulse number")

        # sixth session
        ax = subplot2grid((43, 1900), (23, 971), colspan=447, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        for track in tracks6
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.5, lw=0, markeredgewidth=0.3)
        end
        for line in lines6
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        xl = xlim()
        ylim(yl)
        # sixth session
        subplot2grid((43, 1900), (33, 971), colspan=447, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr6, ysp6, lw=3, alpha=1.0, c="grey")
        for inc in inclines6
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        xlim(xl)
        ylim(yl2)
        figtext(0.5, 0.01, "Pulse number", size=7, ha="center")

        # seventh session
        ax = subplot2grid((43, 1900), (23, 1454), colspan=443, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        #tick_params(axis="x", which="both", direction="in", labeltop=false, labelbottom=true, top=false)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true, labelright=true)
        minorticks_on()
        for track in tracks7
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.5, lw=0, markeredgewidth=0.3)
        end
        for line in lines7
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        xl = xlim()
        ylim(yl)
        # seventh session
        subplot2grid((43, 1900), (33, 1454), colspan=443, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true, labelright=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr7, ysp7, lw=3, alpha=1.0, c="grey")
        for inc in inclines7
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        xlim(xl)
        ylim(yl2)

        println("$outdir/$(name_mod)_driftrate_2.pdf")
        savefig("$outdir/$(name_mod)_driftrate_2.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
        #clf()

    end


    function get_slopes_breakpoints(tracks, lambda; npsi=nothing, fixedpsi=nothing, preview=false)
        if npsi == nothing
            npsi = zeros(length(tracks))
        end

        if fixedpsi == nothing
            fixedpsi = []
            for i in 1:length(tracks)
                push!(fixedpsi, nothing)
            end
        end

        slopes_= []
        eslopes_ = []
        bps_ = []
        ebps_ = []
        xs_ = []
        exs_ = []
        x_ = []
        y_ = []

        # fit broken lines to get slopes and breakpoints
        for (i, track) in enumerate(tracks)
            println(i, " ", npsi[i])
            x = track[:, 2] # pulse number
            y = track[:, 1] # longitude
            x, y2, ysl, slopes, eslopes, bp, ebp, xl, exl = PyRModule.segmented(x, y, npsi[i]; lambda=lambda, fixedpsi=fixedpsi[i], preview=preview)
            slopes_ = vcat(slopes_, slopes)
            eslopes_ = vcat(eslopes_, eslopes)
            bps_ = vcat(bps_, bp)
            ebps_ = vcat(ebps_, ebp)
            xs_ = vcat(xs_, xl)
            exs_ = vcat(exs_, exl)
            push!(x_, x)
            push!(y_, y2)
        end
        return slopes_, eslopes_, bps_, ebps_, xs_, exs_, x_, y_

    end

    function driftrate_J1750_3(outdir; spar=0.6, lambda=1000, name_mod="1234567", bin=1024, show_=false)

        tracks1, lines1, inclines1, dr1, ysp1, dof1 = Tools.get_driftrate("$outdir/tracks/1", lambda)
        tracks2, lines2, inclines2, dr2, ysp2, dof2 = Tools.get_driftrate("$outdir/tracks/2", lambda)
        tracks3, lines3, inclines3, dr3, ysp3, dof3 = Tools.get_driftrate("$outdir/tracks/3", lambda)
        tracks4, lines4, inclines4, dr4, ysp4, dof4 = Tools.get_driftrate("$outdir/tracks/4", lambda)
        tracks5, lines5, inclines5, dr5, ysp5, dof5 = Tools.get_driftrate("$outdir/tracks/5", lambda)
        tracks6, lines6, inclines6, dr6, ysp6, dof6 = Tools.get_driftrate("$outdir/tracks/6", lambda)
        tracks7, lines7, inclines7, dr7, ysp7, dof7 = Tools.get_driftrate("$outdir/tracks/7", lambda)

        slopes1_, eslopes1_, bps1_, ebps1_, xs1_, exs1_, x1_, y1_ = get_slopes_breakpoints(tracks1, lambda, preview=false)
        slopes2_, eslopes2_, bps2_, ebps2_, xs2_, exs2_, x2_, y2_ = get_slopes_breakpoints(tracks2, lambda; npsi=[1, 4, 1, 0, 1, 5, 7, 3, 2, 4, 1, 0, 6, 12, 3, 1], fixedpsi=[nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, [496, 510, 530], nothing], preview=false)
        # 13 - 4,6, 14 - 9,12
        # 7 - [889, 895, 909]
        slopes3_, eslopes3_, bps3_, ebps3_, xs3_, exs3_, x3_, y3_ = get_slopes_breakpoints(tracks3, lambda; npsi=[0, 2, 4, 1, 1, 2, 2, 2])
        slopes4_, eslopes4_, bps4_, ebps4_, xs4_, exs4_, x4_, y4_ = get_slopes_breakpoints(tracks4, lambda; npsi=[0, 0, 0, 2, 2, 0, 1, 3, 4, 2])
        slopes5_, eslopes5_, bps5_, ebps5_, xs5_, exs5_, x5_, y5_ = get_slopes_breakpoints(tracks5, lambda; npsi=[0, 0, 1, 1, 3, 2, 3, 1, 0])
        slopes6_, eslopes6_, bps6_, ebps6_, xs6_, exs6_, x6_, y6_ = get_slopes_breakpoints(tracks6, lambda; npsi=[3, 2, 1, 3, 6, 8, 0])
        slopes7_, eslopes7_, bps7_, ebps7_, xs7_, exs7_, x7_, y7_ = get_slopes_breakpoints(tracks7, lambda; npsi=[2, 4, 2, 5, 1, 1, 1, 2, 1])

        a = get_length(tracks1) + get_length(tracks2) + get_length(tracks3) # 1900 is fine
        b = get_length(tracks4) + get_length(tracks5) + get_length(tracks6) + get_length(tracks7) # 1900 is fine
        println("$a $((1900-a)/2) ")
        println("$b $((1900-b)/3)")

        #yl = (141, 210)
        yl = (210, 131)
        #yl = (135, 215)
        #yl = (125, 225)

        yl2 = [-0.9, 0.9]

        rc("font", size=7.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        #figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
        #figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
        figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
        subplots_adjust(left=0.08, bottom=0.07, right=0.95, top=0.92, wspace=0.0, hspace=0.0)
        figtext(0.09, 0.89, "a)", size=10)
        figtext(0.25, 0.89, "b)", size=10)
        figtext(0.75, 0.89, "c)", size=10)
        figtext(0.09, 0.43, "d)", size=10)
        figtext(0.31, 0.43, "e)", size=10)
        figtext(0.54, 0.43, "f)", size=10)
        figtext(0.76, 0.43, "g)", size=10)

        # first session
        ax = subplot2grid((43, 1900), (0, 0), colspan=295, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        minorticks_on()
        ylabel("Longitude (\$ ^{\\circ}\$)")
        for track in tracks1#[1:3]
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.0, lw=0, markeredgewidth=0.2)
        end
        for i in 1:length(x1_)
            plot(x1_[i], y1_[i], lw=1.5, alpha=0.7, c="C2")
        end
        #for bp in bps1_
        #    axvline(x=bp)
        #end
        #=
        for line in lines1
            #plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        =#
        xl = xlim()
        ylim(yl)
        # first session
        subplot2grid((43, 1900), (10, 0), colspan=295, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr1, ysp1, lw=1, alpha=1.0, c="grey")
        #=
        for inc in inclines1
            plot(inc[1], inc[2], lw=0.2, c="black", alpha=1.0)
        end
        =#
        errorbar(xs1_, slopes1_, yerr=eslopes1_, xerr=exs1_, color="none", lw=0.5, marker="_", mec="blue", ecolor="blue", capsize=0, mfc="blue", ms=0.0, zorder=999)
        xlim(xl)
        ylabel("Drift rate \$(^\\circ / P)\$")
        ylim(yl2)

        # second session
        ax = subplot2grid((43, 1900), (0, 359), colspan=1031, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        xlabel("Pulse number")
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        for (i, track) in enumerate(tracks2)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.5, lw=0, markeredgewidth=0.3)
            #text(track[:,2][1], track[:,1][1], "$i")
        end
        for i in 1:length(x2_)
            plot(x2_[i], y2_[i], lw=1.5, alpha=0.7, c="C2")
        end
        #for bp in bps2_
        #    axvline(x=bp)
        #end
        #=
        for line in lines2
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        =#
        xl = xlim()
        ylim(yl)
        # second session
        subplot2grid((43, 1900), (10, 359), colspan=1031, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr2, ysp2, lw=1, alpha=1.0, c="grey")
        #=
        for inc in inclines2
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        =#
        # TODO what is wrong here?
        println(length(xs2_))
        println(length(exs2_))
        println(length(slopes2_))
        println(length(eslopes2_))
        errorbar(xs2_, slopes2_[1:length(xs2_)], yerr=eslopes2_[1:length(xs2_)], xerr=exs2_, color="none", lw=0.5, marker="_", mec="blue", ecolor="blue", capsize=0, mfc="blue", ms=0.0, zorder=999)
        xlim(xl)
        ylim(yl2)

        # third session
        ax = subplot2grid((43, 1900), (0, 1454), colspan=446, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true, labelright=true)
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        #ax.yaxis.set_label_position("right")
        #ax.yaxis.set_ticks_position("right")
        minorticks_on()
        for (i, track) in enumerate(tracks3)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.5, lw=0, markeredgewidth=0.3)
            #text(track[:,2][1], track[:,1][1], "$i")
        end
        for i in 1:length(x3_)
            plot(x3_[i], y3_[i], lw=1.5, alpha=0.7, c="C2")
        end
        #=
        for line in lines3
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        =#
        xl = xlim()
        ylim(yl)
        # third session
        ax = subplot2grid((43, 1900), (10, 1454), colspan=446, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true, labelright=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr3, ysp3, lw=1, alpha=1.0, c="grey")
        #=
        for inc in inclines3
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        =#
        errorbar(xs3_, slopes3_, yerr=eslopes3_, xerr=exs3_, color="none", lw=0.5, marker="_", mec="blue", ecolor="blue", capsize=0, mfc="blue", ms=0.0, zorder=999)
        xlim(xl)
        ylim(yl2)
        #ax.yaxis.set_label_position("right")
        #ax.yaxis.set_ticks_position("right")

        # fourth session
        ax = subplot2grid((43, 1900), (23, 0), colspan=450, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        minorticks_on()
        for (i, track) in enumerate(tracks4)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.5, lw=0, markeredgewidth=0.3)
            #text(track[:,2][1], track[:,1][1], "$i")
        end
        for i in 1:length(x4_)
            plot(x4_[i], y4_[i], lw=1.5, alpha=0.7, c="C2")
        end
        #=
        for line in lines4
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        =#
        xl = xlim()
        ylim(yl)
        ylabel("Longitude (\$ ^{\\circ}\$)")
        # fourth session
        subplot2grid((43, 1900), (33, 0), colspan=450, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr4, ysp4, lw=1, alpha=1.0, c="grey")
        #=
        for inc in inclines4
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        =#
        errorbar(xs4_, slopes4_, yerr=eslopes4_, xerr=exs4_, color="none", lw=0.5, marker="_", mec="blue", ecolor="blue", capsize=0, mfc="blue", ms=0.0, zorder=999)
        xlim(xl)
        ylim(yl2)
        ylabel("Drift rate \$(^\\circ / P)\$")
        locator_params(nbins=3)

        # fifth session
        ax = subplot2grid((43, 1900), (23, 487), colspan=447, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        for (i, track) in enumerate(tracks5)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.5, lw=0, markeredgewidth=0.3)
            #text(track[:,2][1], track[:,1][1], "$i")
        end
        for i in 1:length(x5_)
            plot(x5_[i], y5_[i], lw=1.5, alpha=0.7, c="C2")
        end
        #=
        for line in lines5
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        =#
        xl = xlim()
        ylim(yl)
        # fifth session
        subplot2grid((43, 1900), (33, 487), colspan=447, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr5, ysp5, lw=1, alpha=1.0, c="grey")
        #=
        for inc in inclines5
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        =#
        errorbar(xs5_, slopes5_, yerr=eslopes5_, xerr=exs5_, color="none", lw=0.5, marker="_", mec="blue", ecolor="blue", capsize=0, mfc="blue", ms=0.0, zorder=999)
        xlim(xl)
        ylim(yl2)
        #xlabel("Pulse number")

        # sixth session
        ax = subplot2grid((43, 1900), (23, 971), colspan=447, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        for (i, track) in enumerate(tracks6)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.5, lw=0, markeredgewidth=0.3)
            #text(track[:,2][1], track[:,1][1], "$i")
        end
        for i in 1:length(x6_)
            plot(x6_[i], y6_[i], lw=1.5, alpha=0.7, c="C2")
        end
        #=
        for line in lines6
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        =#
        xl = xlim()
        ylim(yl)
        # sixth session
        subplot2grid((43, 1900), (33, 971), colspan=447, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr6, ysp6, lw=1, alpha=1.0, c="grey")
        #=
        for inc in inclines6
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        =#
        errorbar(xs6_, slopes6_, yerr=eslopes6_, xerr=exs6_, color="none", lw=0.5, marker="_", mec="blue", ecolor="blue", capsize=0, mfc="blue", ms=0.0, zorder=999)
        xlim(xl)
        ylim(yl2)
        figtext(0.5, 0.01, "Pulse number", size=7, ha="center")

        # seventh session
        ax = subplot2grid((43, 1900), (23, 1454), colspan=443, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        #tick_params(axis="x", which="both", direction="in", labeltop=false, labelbottom=true, top=false)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true, labelright=true)
        minorticks_on()
        for (i, track) in enumerate(tracks7)
            plot(track[:,2], track[:,1], marker="x", color="red", markersize=1.5, lw=0, markeredgewidth=0.3)
            #text(track[:,2][1], track[:,1][1], "$i")
        end
        for i in 1:length(x7_)
            plot(x7_[i], y7_[i], lw=1.5, alpha=0.7, c="C2")
        end
        #=
        for line in lines7
            plot(line[1], line[2], lw=1.5, alpha=0.7, c="C2")
        end
        =#
        xl = xlim()
        ylim(yl)
        # seventh session
        subplot2grid((43, 1900), (33, 1454), colspan=443, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true, labelright=true)
        minorticks_on()
        axhline(y=0, lw=1, ls="--")
        plot(dr7, ysp7, lw=1, alpha=1.0, c="grey")
        #=
        for inc in inclines7
            plot(inc[1], inc[2], lw=0.4, c="black", alpha=0.6)
        end
        =#
        errorbar(xs7_, slopes7_, yerr=eslopes7_, xerr=exs7_, color="none", lw=0.5, marker="_", mec="blue", ecolor="blue", capsize=0, mfc="blue", ms=0.0, zorder=999)
        xlim(xl)
        ylim(yl2)

        println("$outdir/$(name_mod)_driftrate_3.pdf")
        savefig("$outdir/$(name_mod)_driftrate_3.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
        #clf()

    end


    function driftrate_analysis_J1750(outdir; lambda=100, name_mod="123456", show_=false)

        nd, pd = Tools.driftrate_analysis_J1750(outdir, lambda)

        println("Negative mean: ", mean(nd))
        println("Negative median: ", median(nd))
        println("Negative std: ", std(nd))

        println("Positive mean: ", mean(pd))
        println("Positive median: ", median(pd))
        println("Positive std: ", std(pd))

        println("Negative sum: ", sum(nd))
        println("Positive sum: ", sum(pd))

        #figure(figsize=(7.086614, 6.299213))  # 18 cm x 16 cm
        hist(nd, alpha=0.7)
        hist(pd, alpha=0.7)

        println("$outdir/$(name_mod)_ratehist.pdf")
        savefig("$outdir/$(name_mod)_ratehist.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
    end


    function driftrate_analysis_J1750_2(outdir; lambda=100, name_mod="1234567", show_=false, datas=nothing)

        nd, pd, null_pulses = Tools.driftrate_analysis_J1750_2(outdir, lambda)

        println()
        println("(in p.) OLD Negative mean: ", mean(nd))
        println("Negative median: ", median(nd))
        println("Negative std: ", std(nd))
        println("(in p.) OLD Neg. Standard error of the mean: ", std(nd) / sqrt(length(nd)))
        println()
        println("(in p.) OLD Positive mean: ", mean(pd))
        println("Positive median: ", median(pd))
        println("Positive std: ", std(pd))
        println("(in p.) OLD Pos. Standard error of the mean: ", std(pd) / sqrt(length(pd)))
        println()
        println("Negative sum: ", sum(nd))
        println("Positive sum: ", sum(pd))
        println("(in p.) OLD Positive coverage: ",  sum(pd)/ (sum(pd) + sum(nd)))

        # not very useful
        #=
        if datas != nothing
            for i in 1:size(datas)[1]
                for j in 1:size(null_pulses[i])[1]
                    println(null_pulses[i][j])
                    single(datas[i], outdir; darkness=1.0, start=trunc(Int, null_pulses[i][j])-2, number=5, show_=true)
                    q = readline(stdin; keep=false)
                    if q == "q"
                        return
                    end
                end
            end

        end
        =#

        #figure(figsize=(7.086614, 6.299213))  # 18 cm x 16 cm
        hist(nd, alpha=0.7)
        hist(pd, alpha=0.7)

        println("$outdir/$(name_mod)_ratehist.pdf")
        savefig("$outdir/$(name_mod)_ratehist.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
    end


    function driftdirection_analysis_J1750_4(outdir; lambda=200, name_mod="1234567", show_=false)

        tracks = []
        lines = []
        inclines = []
        pulses = []
        ysps = []

        for i in 1:7
            tracks1, lines1, inclines1, dr1, ysp1 = Tools.get_driftrate("$outdir/tracks/$i", lambda)
            push!(tracks, tracks1)
            push!(lines, lines1)
            push!(inclines, inclines1)
            push!(pulses, dr1)
            push!(ysps, ysp1)
        end

        breaks_longitude = []


        for i in 1:length(tracks)
            negative = []
            positive = []
            for (j, ysp) in enumerate(ysps[i])
                if signbit(ysp)
                    push!(negative, pulses[i][j])
                else
                    push!(positive, pulses[i][j])
                end
            end

            #if i == 2
            for j in 1:size(positive)[1]-1
                dpulse = positive[j+1] - positive[j]
                if dpulse > 2
                    break_pulse = trunc(positive[j])
                    #println("Break pulse: $break_pulse")
                    #ff = findfirst(x->x==positive[j], pulses[i])
                    for k in 1:size(tracks[i])[1]
                        #println("$k", size(tracks[i][k]))
                        #println("$k ", pulses[i])
                        for l in 1:size(tracks[i][k])[1]
                            pulse_number = tracks[i][k][l,2]
                            lon = tracks[i][k][l,1]
                            #println(pulse_number, " ", break_pulse)
                            if pulse_number == break_pulse
                                #println(lon)
                                push!(breaks_longitude, lon)
                            end
                        end
                    end
                end
                #end
                #println(positive)
                #println(negative)
            end
        end


        mi = minimum(breaks_longitude)
        ma = maximum(breaks_longitude)
        println(length(breaks_longitude))


        rc("font", size=7.5)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 1.946563744), frameon=true)  # 8cm x 4.94427191 cm (golden)
        subplots_adjust(left=0.14, bottom=0.20, right=0.99, top=0.99, wspace=0., hspace=0.)
        minorticks_on()

        res = hist(breaks_longitude, bins=12, color="tab:blue", edgecolor="black", lw=0.5)
        #xlim([mi, ma])
        xlim([147, 207])

        #y_ = res[1]
        y_ = vcat([0,0], res[1], [0, 0])
        x_ = []
        for i in 1:size(res[2])[1]-1
            push!(x_, res[2][i] + (res[2][i+1]-res[2][i]) / 2)
        end
        x_ = vcat([x_[1]-(x_[2]-x_[1])], x_, [x_[end]+(x_[end]-x_[end-1])])
        x_ = vcat([x_[1]-(x_[2]-x_[1])], x_, [x_[end]+(x_[end]-x_[end-1])])

        p, err = Tools.fit_twogaussians(x_, y_, 5.0, 8.0, 162.0, 180.0, 10.0, 10.0; baselevel=0.0)
        #p, err = Tools.fit_gaussian(x_, y_; a=8.0, μ=180.0, σ=5.0, baselevel=0.0)
        println(p)
        println(err)
        xx = collect(x_[1]:0.1:x_[end])
        ga = Tools.twogauss(xx, p)
        plot(xx, ga, lw=1, c="tab:red")
        xlabel("longitude (deg.)")
        ylabel("Number of direction changes")

        savefig("$outdir/$(name_mod)_driftdirection_4.pdf")
        println("$outdir/$(name_mod)_driftdirection_4.pdf")

        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
    end



    function p2_analysis_J1750(outdir; lambda=200, name_mod="1234567", show_=false)

        tracks = []
        lines = []
        inclines = []
        pulses = []
        ysps = []

        for i in 1:7
            tracks1, lines1, inclines1, dr1, ysp1 = Tools.get_driftrate("$outdir/tracks/$i", lambda)
            push!(tracks, tracks1)
            push!(lines, lines1)
            push!(inclines, inclines1)
            push!(pulses, dr1)
            push!(ysps, ysp1)
        end

        lons = []
        lonsp = []
        lonsn = []
        # i - session
        added_pulses = []
        for i in 1:length(tracks)
            push!(added_pulses, [])
            for j in 1:length(tracks[i])
                for k in 1:length(tracks[i])
                    for l in 1:size(tracks[i][j])[1]
                        for m in 1:size(tracks[i][k])[1]
                            dpulse = abs(tracks[i][j][l, 2] - tracks[i][k][m, 2])
                            pulse = trunc(tracks[i][j][l, 2])
                            if (dpulse == 0) && ~(pulse in added_pulses[end])
                                dlon = abs(tracks[i][j][l, 1] - tracks[i][k][m, 1])
                                if (dlon != 0) && (dlon < 30)
                                    push!(lons, dlon)
                                    push!(added_pulses[end], pulse)
                                    # check drift rate
                                    for n in 1:size(pulses[i])[1]
                                        p = trunc(pulses[i][n])
                                        if pulse == p
                                            if ysps[i][n] > 0
                                                push!(lonsp, dlon)
                                            else
                                                push!(lonsn, dlon)
                                            end
                                            break
                                        end
                                    end
                                    #println(dlon)
                                end
                                #println(j, " ", k)
                            end

                        end
                    end
                end
            end
        end

        #=
        # Nope use gaussian fitting instead
        println(round(median(lons), digits=2), " +- ", round(std(lons), digits=2))
        println("\t pos:", round(median(lonsp), digits=2), " +- ", round(std(lonsp), digits=2))
        if length(lonsn) > 0
            println("\t neg:", round(median(lonsn), digits=2), " +- ", round(std(lonsn), digits=2))
        end=#

        # save data
        #=
        f = open("all.txt", "w")
        write(f, "# P2\n")
        for (i,l) in enumerate(lons)
            write(f, "$l\n")
        end
        close(f)

        f = open("pos.txt", "w")
        write(f, "# P2\n")
        for (i,l) in enumerate(lonsp)
            write(f, "$l\n")
        end
        close(f)

        f = open("neg.txt", "w")
        write(f, "# P2\n")
        for (i,l) in enumerate(lonsn)
            write(f, "$l\n")
        end
        close(f)
        =#


        rc("font", size=7.5)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 1.946563744), frameon=true)  # 8cm x 4.94427191 cm (golden)
        subplots_adjust(left=0.15, bottom=0.20, right=0.99, top=0.99, wspace=0., hspace=0.)
        minorticks_on()
        bins = 18
        res = hist(lons, color="black", bins=bins, lw=0.5, alpha=0.5, edgecolor="black", label="All pulses", ls="--")
        res2 = hist(lonsp, color="tab:red", bins=res[2], lw=0.5, alpha=0.5, edgecolor="black", label="Positive drift", ls="--")
        res3 = hist(lonsn, color="tab:blue", bins=res[2], lw=0.5, alpha=0.5, edgecolor="black", label="Negative drift", ls="--")

        y_ = res[1]
        #y_ = vcat([0,0], res[1], [0, 0])
        x_ = []
        for i in 1:size(res[2])[1]-1
            push!(x_, res[2][i] + (res[2][i+1]-res[2][i]) / 2)
        end
        p, err = Tools.fit_gaussian(x_, y_; a=1700.0, μ=20.0, σ=3.0, baselevel=0.0)
        println(round(p[2], digits=2), " +- ", round(err[2], digits=2))
        xx = collect(x_[1]:0.1:x_[end])
        ga = Tools.gauss(xx, p)

        plot(xx, ga, lw=1, c="black")

        y_ = res2[1]
        #y_ = vcat([0,0], res[1], [0, 0])
        x_ = []
        for i in 1:size(res2[2])[1]-1
            push!(x_, res2[2][i] + (res2[2][i+1]-res2[2][i]) / 2)
        end
        p, err = Tools.fit_gaussian(x_, y_; a=1700.0, μ=20.0, σ=3.0, baselevel=0.0)
        println("pos. :", round(p[2], digits=2), " +- ", round(err[2], digits=2))
        xx = collect(x_[1]:0.1:x_[end])
        ga = Tools.gauss(xx, p)

        plot(xx, ga, lw=1, c="tab:red")

        y_ = res3[1]
        #y_ = vcat([0,0], res[1], [0, 0])
        x_ = []
        for i in 1:size(res3[2])[1]-1
            push!(x_, res3[2][i] + (res3[2][i+1]-res3[2][i]) / 2)
        end
        p, err = Tools.fit_gaussian(x_, y_; a=1000.0, μ=20.0, σ=3.0, baselevel=0.0)
        println("neg. :", round(p[2], digits=2), " +- ", round(err[2], digits=2))
        xx = collect(x_[1]:0.1:x_[end])
        ga = Tools.gauss(xx, p)

        plot(xx, ga, lw=1, c="tab:blue")
        legend(prop=Dict("size"=> 6.7))

        xlabel("\$P_{2}\$ (deg.)")
        ylabel("Number of pulses")

        savefig("$outdir/$(name_mod)_p2_analysis.pdf")
        println("$outdir/$(name_mod)_p2_analysis.pdf")

        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
    end


    function driftdirection_analysis_J1750_3(outdir; lambda=200, name_mod="1234567", show_=false)

        tracks = []
        lines = []
        inclines = []
        pulses = []
        ysps = []

        for i in 1:7
            tracks1, lines1, inclines1, dr1, ysp1 = Tools.get_driftrate("$outdir/tracks/$i", lambda)
            push!(tracks, tracks1)
            push!(lines, lines1)
            push!(inclines, inclines1)
            push!(pulses, dr1)
            push!(ysps, ysp1)
        end

        npsis = [nothing,
            [1, 4, 1, 0, 1, 5, 7, 3, 2, 4, 1, 0, 6, 12, 3, 1],
            [0, 2, 4, 1, 1, 2, 2, 2],
            [0, 0, 0, 2, 2, 0, 1, 3, 4, 2],
            [0, 0, 1, 1, 3, 2, 3, 1, 0],
            [3, 2, 1, 3, 6, 8, 0],
            [2, 4, 2, 5, 1, 1, 1, 2, 1]
        ]
        fixedpsis = [nothing,
            [nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, [496, 510, 530], nothing],
            nothing, nothing, nothing, nothing, nothing
        ]
        slopes, eslopes, bps, ebps, xs, exs = [], [], [], [], [], []

        for i in 1:7
            slopesi_, eslopesi_, bpsi_, ebpsi_, xsi_, exsi_, xi_, yi_ = get_slopes_breakpoints(tracks[i], lambda; npsi=npsis[i], fixedpsi=fixedpsis[i], preview=false)
            push!(slopes, slopesi_)
            push!(eslopes, eslopesi_)
            push!(bps, bpsi_)
            push!(ebps, ebpsi_)
            push!(xs, xsi_)
            push!(exs, exsi_)
        end

        println("\n\n")

        ndur = []
        pdur = []
        ndrate = []
        pdrate = []
        for i in 1:7
            (nd,ndr, pd, pdr) = Tools.driftdirection_analysis_J1750_3(slopes, eslopes, bps, ebps, xs, exs, i)
            pdur = vcat(pdur, pd)
            ndur = vcat(ndur, nd)
            ndrate = vcat(ndrate, ndr)
            pdrate = vcat(pdrate, pdr)
        end
        println("Positive:")
        #println(pdur)
        println("Duration: ", mean(pdur), " +/- ", std(pdur))
        println("Drift rate: ", mean(pdrate), " +/- ", std(pdrate), "max. ", maximum(pdrate))

        println("Negative:")
        #println(ndur)
        println("Duration: ", mean(ndur), " +/- ", std(ndur))
        println("Drift rate: ", mean(ndrate), " +/- ", std(ndrate), "min. ", minimum(ndrate))

    end


    """ calculates average profiels in four ranges (too much, check _2) """
    function average_J1750(datas, outdir; lambda=1000.0, bin_st=nothing, bin_end=nothing, name_mod="0", show_=false)

        nums = []
        bins = []
        for data in datas
            nu, bi = size(data)
            push!(nums, nu)
            push!(bins, bi)
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins[1] end
        das = [] # single pulse data
        for data in datas
            da = data[:, bin_st:bin_end]
            push!(das, da)
        end

        # drift rate data
        pulses = []
        ysps = []  # drift rate

        for i in 1:length(das)
            tracks1, lines1, inclines1, dr1, ysp1 = Tools.get_driftrate("$outdir/tracks/$i", lambda)
            (ysp, dr) = Tools.remove_duplicates_ysp(ysp1, dr1)
            push!(pulses, dr)
            push!(ysps, ysp)
        end

        ranges = [(1, 0), (0, -1), (1, 0.455), (0.455, 0), (0, -0.365), (-0.365, -1)]
        profiles = []
        averages = []

        for r in ranges
            push!(profiles, [])
            for i in 1:length(das)
                for (j, driftrate) in enumerate(ysps[i])
                    if (driftrate < r[1]) && (driftrate > r[2])
                        push!(profiles[end], das[i][pulses[i][j], :] )
                    end
                end
            end
        end

        ##=
        # make profiles equal
        le = 270
        for i in 1:length(profiles)
            indx = []
            num = 0
            while num < le
                ind = rand(1:length(profiles[i]))
                if ~(ind in indx)
                    push!(indx, ind)
                    num += 1
                end
            end
            pr = [profiles[i][x] for x in indx]
            profiles[i] = pr
        end
        #

        # converting profiles TODO why why why?
        for ii in 1:length(profiles)
            println(ii)
            x, = size(profiles[ii])
            y, = size(profiles[ii][1]) # bins
            profs = zeros((x,y))
            for i in 1:x
                for j in 1:y
                    profs[i,j] = profiles[ii][i][j]
                end
            end
            profiles[ii] = profs
        end

        nums = []
        for profile in profiles
            push!(nums, size(profile)[1])
            push!(averages, Tools.average_profile(profile; norm=true))
            #println(size(profile)[1])
        end

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins[1]
        longitude = collect(range(-dl/2., dl/2., length=db))

        (p1, err1) = Tools.fit_twogaussians(longitude, averages[1], 0.5, 0.9, -13.0, 3.0, 5.0, 5.0)
        (p2, err2) = Tools.fit_twogaussians(longitude, averages[2], 0.5, 0.9, -13.0, 3.0, 5.0, 5.0)
        (p3, err3) = Tools.fit_twogaussians(longitude, averages[3], 0.5, 0.9, -13.0, 3.0, 5.0, 5.0)
        (p4, err4) = Tools.fit_twogaussians(longitude, averages[4], 0.5, 0.9, -13.0, 3.0, 5.0, 5.0)
        ga = Tools.twogauss(longitude, p4)
        (p5, err5) = Tools.fit_twogaussians(longitude, averages[5], 0.6, 0.9, -13.0, 3.0, 5.0, 5.0)
        (p6, err6) = Tools.fit_gaussian(longitude, averages[6], a=1.0, μ=-4.0, σ=15.0)
        println("1: ", p1)
        println("3: ", p3)
        println("4: ", p4)
        #println("err1: ", err1)
        println("2: ", p2)
        println("5: ", p5)
        println("6: ", p6)
        #println(err2)

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        labels = ["\$ \\qquad \\;\\; {\\rm D} > 0\$", "\$\\qquad \\;\\; {\\rm D} < 0 \$", "\$0.91 > {\\rm D} > \\quad \\; 0.5 \$", "\$\\quad 0.5\\;\\; > {\\rm D} > \\quad \\; 0 \$", "\$ \\quad 0 \\;\\;  > {\\rm D} > \\;-0.5 \$", "\$\\!-0.5 \\;\\; > {\\rm D} > --0.73 \$"]


        figure(figsize=(3.14961, 1.946563744), frameon=true)  # 8cm x 4.94427191 cm (golden)
        subplots_adjust(left=0.17, bottom=0.19, right=0.99, top=0.99, wspace=0., hspace=0.)

        minorticks_on()
        #plot(longitude, ga * sqrt(nums[4]), c="green", lw=1.3, zorder=1200)
        plot(longitude, averages[1] * sqrt(nums[1]), c="black", label=labels[1], lw=0.3, zorder=200)
        plot(longitude, averages[3] * sqrt(nums[3]), c="C1", label=labels[3], lw=0.7, alpha=0.7, ls=(0, (1, 1)))
        plot(longitude, averages[4] * sqrt(nums[4]), c="C2", label=labels[4], lw=0.7, alpha=0.7, ls=(0, (5, 1)))
        plot(longitude, averages[2] * sqrt(nums[2]), c="black", label=labels[2], lw=0.3, zorder=201)
        plot(longitude, averages[5] * sqrt(nums[5]), c="C3", label=labels[5], lw=0.7,alpha=0.7, ls=(0, (1, 1)))
        plot(longitude, averages[6] * sqrt(nums[6]), c="C4", label=labels[6], lw=0.7, alpha=0.7, ls=(0, (5, 1)))
        #for i in 1:length(averages)
        #    plot(longitude, averages[i] * sqrt(nums[i]), c="C$i", label=labels[i])
            #plot(longitude, averages[i], c="C$i", label=labels[i])
        #end
        #yticks([0.0, 0.5])
        #xlim(longitude[1], longitude[end])
        xlabel("longitude (deg.)")
        ylabel("\$\\sqrt{\\rm Number \\; of \\; pulses}\$")
        legend(prop=Dict("size"=> 5.1))
        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_averages.pdf")
        savefig("$outdir/$(name_mod)_averages.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()

    end


    """ calculates average profiels in two ranges """
    function average_J1750_2(datas, outdir; lambda=1000.0, bin_st=nothing, bin_end=nothing, name_mod="0", show_=false)

        nums = []
        bins = []
        for data in datas
            nu, bi = size(data)
            push!(nums, nu)
            push!(bins, bi)
        end
        # take full but plot only part
        bin_st0 = 1
        bin_end0 = bins[1]
        das = [] # single pulse data
        for data in datas
            da = data[:, bin_st0:bin_end0]
            push!(das, da)
        end

        # drift rate data
        pulses = []
        ysps = []  # drift rate

        for i in 1:length(das)
            tracks1, lines1, inclines1, dr1, ysp1 = Tools.get_driftrate("$outdir/tracks/$i", lambda)
            (ysp, dr) = Tools.remove_duplicates_ysp(ysp1, dr1)
            push!(pulses, dr)
            push!(ysps, ysp)
        end

        ranges = [(1, 0), (0, -1), (2, -2)]
        profiles = []
        averages = []

        for r in ranges
            push!(profiles, [])
            for i in 1:length(das)
                for (j, driftrate) in enumerate(ysps[i])
                    if (driftrate < r[1]) && (driftrate > r[2])
                        push!(profiles[end], das[i][pulses[i][j], :] )
                    end
                end
            end
        end

        #println(length(profiles[3]))
        #return

        #=
        # make profiles from equal number of pulses
        le = 730
        for i in 1:length(profiles)
            indx = []
            num = 0
            while num < le
                ind = rand(1:length(profiles[i]))
                if ~(ind in indx)
                    push!(indx, ind)
                    num += 1
                end
            end
            pr = [profiles[i][x] for x in indx]
            profiles[i] = pr
        end
        =#


        # converting profiles TODO why why why?
        for ii in 1:length(profiles)
            println(ii)
            x, = size(profiles[ii])
            y, = size(profiles[ii][1]) # bins
            profs = zeros((x,y))
            for i in 1:x
                for j in 1:y
                    profs[i,j] = profiles[ii][i][j]
                end
            end
            profiles[ii] = profs
        end

        nums = []
        for profile in profiles
            push!(nums, size(profile)[1])
            push!(averages, Tools.average_profile(profile; norm=false))
            #println(size(profile)[1])
        end

        # normalize profiles here
        for i in 1:3
            me = mean(averages[i][700:1000])
            averages[i] .-= me
        end
        ma = maximum(averages[3])
        averages[3] ./= ma
        ma1 = maximum(averages[1])
        ma2 = maximum(averages[2])
        m1 = mean(averages[1][350:650])
        m2 = mean(averages[2][350:650])
        m3 = mean(averages[3][350:650])
        #averages[1] ./= ma1
        #averages[2] ./= ma2
        averages[1] .*= m3 / m1
        averages[2] .*= m3 / m2

        # Pulse longitude around 0
        #db = (bin_end + 1) - bin_st  # yes +1
        #dl = 360. * db / bins[1]
        #longitude = collect(range(-dl/2., dl/2., length=db))

        # Pulse longitude (raw)
        longitude = collect(range(bin_st/bins[1] * 360, bin_end/bins[1] * 360, length=bin_end - bin_st + 1))

        #=
        (p1, err1) = Tools.fit_twogaussians(longitude, averages[1], 0.5, 0.9, -13.0, 3.0, 5.0, 5.0)
        (p2, err2) = Tools.fit_twogaussians(longitude, averages[2], 0.5, 0.9, -13.0, 3.0, 5.0, 5.0)
        (p3, err3) = Tools.fit_twogaussians(longitude, averages[3], 0.5, 0.9, -13.0, 3.0, 5.0, 5.0)
        (p4, err4) = Tools.fit_twogaussians(longitude, averages[4], 0.5, 0.9, -13.0, 3.0, 5.0, 5.0)
        ga = Tools.twogauss(longitude, p4)
        (p5, err5) = Tools.fit_twogaussians(longitude, averages[5], 0.6, 0.9, -13.0, 3.0, 5.0, 5.0)
        (p6, err6) = Tools.fit_gaussian(longitude, averages[6], a=1.0, μ=-4.0, σ=15.0)
        println("1: ", p1)
        println("3: ", p3)
        println("4: ", p4)
        #println("err1: ", err1)
        println("2: ", p2)
        println("5: ", p5)
        println("6: ", p6)
        #println(err2)
        =#
        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        labels = ["\$ \\qquad \\;\\; {\\rm D} > 0\$ (N=$(size(profiles[1])[1]))", "\$\\qquad \\;\\; {\\rm D} < 0 \$ (N=$(size(profiles[2])[1]))", "\$0.91 > {\\rm D} > -0.73 \$ (N=$(size(profiles[3])[1]))", "\$\\quad 0.5\\;\\; > {\\rm D} > \\quad \\; 0 \$", "\$ \\quad 0 \\;\\;  > {\\rm D} > \\;-0.5 \$", "\$\\!-0.5 \\;\\; > {\\rm D} > --0.73 \$"]


        figure(figsize=(3.14961, 1.946563744), frameon=true)  # 8cm x 4.94427191 cm (golden)
        subplots_adjust(left=0.17, bottom=0.19, right=0.99, top=0.99, wspace=0., hspace=0.)

        minorticks_on()
        plot(longitude, averages[1][bin_st:bin_end], c="C1", label=labels[1], lw=0.3, zorder=200)
        plot(longitude, averages[2][bin_st:bin_end], c="C2", label=labels[2], lw=0.3, zorder=201)
        plot(longitude, averages[3][bin_st:bin_end], c="black", label=labels[3], lw=0.5, zorder=199, alpha=0.9)
        #plot(longitude, ga * sqrt(nums[4]), c="green", lw=1.3, zorder=1200)
        #plot(longitude, averages[1] * sqrt(nums[1]), c="C1", label=labels[1], lw=0.3, zorder=200)
        #plot(longitude, averages[3] * sqrt(nums[3]), c="C1", label=labels[3], lw=0.7, alpha=0.7, ls=(0, (1, 1)))
        #plot(longitude, averages[4] * sqrt(nums[4]), c="C2", label=labels[4], lw=0.7, alpha=0.7, ls=(0, (5, 1)))
        #plot(longitude, averages[2] * sqrt(nums[2]), c="C2", label=labels[2], lw=0.3, zorder=201)
        #plot(longitude, averages[5] * sqrt(nums[5]), c="C3", label=labels[5], lw=0.7,alpha=0.7, ls=(0, (1, 1)))
        #plot(longitude, averages[6] * sqrt(nums[6]), c="C4", label=labels[6], lw=0.7, alpha=0.7, ls=(0, (5, 1)))
        #for i in 1:length(averages)
        #    plot(longitude, averages[i] * sqrt(nums[i]), c="C$i", label=labels[i])
            #plot(longitude, averages[i], c="C$i", label=labels[i])
        #end
        yticks(collect(0.0:0.2:1))
        #xlim(longitude[1], longitude[end])
        xlabel("longitude (deg.)")
        #ylabel("\$\\sqrt{\\rm Number \\; of \\; pulses}\$")
        #ylabel("Intensity (arbitrary units)")
        ylabel("Intensity (a. u.)")
        legend(prop=Dict("size"=> 4.1))
        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_averages2.pdf")
        savefig("$outdir/$(name_mod)_averages2.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()

    end


    """ calculates average profiels in two ranges - different normalization"""
    function average_J1750_3(datas, outdir; lambda=1000.0, bin_st=nothing, bin_end=nothing, name_mod="0", show_=false)

        nums = []
        bins = []
        for data in datas
            nu, bi = size(data)
            push!(nums, nu)
            push!(bins, bi)
        end
        # take full but plot only part
        bin_st0 = 1
        bin_end0 = bins[1]
        das = [] # single pulse data
        for data in datas
            da = data[:, bin_st0:bin_end0]
            push!(das, da)
        end

        # drift rate data
        pulses = []
        ysps = []  # drift rate

        for i in 1:length(das)
            tracks1, lines1, inclines1, dr1, ysp1 = Tools.get_driftrate("$outdir/tracks/$i", lambda)
            (ysp, dr) = Tools.remove_duplicates_ysp(ysp1, dr1)
            push!(pulses, dr)
            push!(ysps, ysp)
        end

        ranges = [(2, 0), (0, -2), (2, -2)]
        profiles = []
        averages = []

        for r in ranges
            push!(profiles, [])
            for i in 1:length(das)
                for (j, driftrate) in enumerate(ysps[i])
                    if (driftrate < r[1]) && (driftrate > r[2])
                        push!(profiles[end], das[i][pulses[i][j], :] )
                    end
                end
            end
        end

        # converting profiles TODO why why why?
        for ii in 1:length(profiles)
            println(ii)
            x, = size(profiles[ii])
            y, = size(profiles[ii][1]) # bins
            profs = zeros((x,y))
            for i in 1:x
                for j in 1:y
                    profs[i,j] = profiles[ii][i][j]
                end
            end
            profiles[ii] = profs
        end

        nums = []
        for profile in profiles
            push!(nums, size(profile)[1])
            push!(averages, Tools.average_profile(profile; norm=false))
            #println(size(profile)[1])
        end

        # normalize profiles here
        # OLD
        #=
        for i in 1:3
            me = mean(averages[i][700:1000])
            averages[i] .-= me
        end
        ma = maximum(averages[3])
        averages[3] ./= ma
        ma1 = maximum(averages[1])
        ma2 = maximum(averages[2])
        m1 = mean(averages[1][350:650])
        m2 = mean(averages[2][350:650])
        m3 = mean(averages[3][350:650])
        #averages[1] ./= ma1
        #averages[2] ./= ma2
        averages[1] .*= m3 / m1
        averages[2] .*= m3 / m2
        =#
        # NEW
        # around zero
        for i in 1:3
            me = mean(averages[i][700:1000])
            averages[i] .-= me
            averages[i] ./= size(profiles[i])[1] # multiply based on number of pulses
        end

        # Pulse longitude (raw)
        longitude = collect(range(bin_st/bins[1] * 360, bin_end/bins[1] * 360, length=bin_end - bin_st + 1))

        stds = []

        ma = maximum(averages[3]) # normalize to profile for all single pulses
        for i in 1:3
            averages[i] ./= ma
            #println("mean ", mean(averages[i][350:650]))
            #me = mean(averages[i][700:1000])
            st = std(averages[i][700:1000])
            push!(stds, st)
            #rms = sqrt(mean(averages[i][700:1000].^2))
            #rms2 = norm(averages[i][700:1000]) / sqrt(length(averages[i][700:1000]))
            #println("me ", me)
            #println("st ", st)
            #println("rms ", rms)
            #println("rms2 ", rms2)
            #auc = PyRModule.auc(longitude, averages[i][700:1000])
            #auc2 = PyRModule.auc(longitude, averages[i][bin_st:bin_end])
            #println("auc ", auc)
            #println("auc2 ", auc2)
        end

        # normalize to one
        no = mean(averages[3][bin_st:bin_end])
        for i in 1:3
            av = averages[i][bin_st:bin_end]
            println("Mean: ", mean(av) / no)
            println("std prof.: ", std(av) / no)
            println("SE: ", std(av) / sqrt(bin_end-bin_st) /no)
            println("std baseline: ", stds[i] / no)
        end
        println("CONSIDER ABOVE!")

        # random area calculations
        repeat = 10000
        areas = []

        for i in 1:3
            push!(areas, [])
            for j in 1:repeat
                # randomize signal based on std/rms
                av = averages[i][bin_st:bin_end]
                for k in 1:size(av)[1]
                    d = Normal(av[k], stds[i])
                    av[k] = rand(d)
                end
                auc = PyRModule.auc(longitude, av)
                push!(areas[end], auc)
            end
        end
        # normalize areas
        meme = mean(areas[3])
        areas ./= meme
        for i in 1:3
            println("Area = ", mean(areas[i]), " ± ", std(areas[i]), "underestimated errors!")
        end

        # simple error calculation (overestimate?)
        areas = []
        for i in 1:3
            push!(areas, [])
            # randomize signal based on std/rms
            av = averages[i][bin_st:bin_end]
            auc1 = PyRModule.auc(longitude, av.-stds[i])
            auc2 = PyRModule.auc(longitude, av)
            auc3 = PyRModule.auc(longitude, av.+stds[i])
            push!(areas[end], auc1)
            push!(areas[end], auc2)
            push!(areas[end], auc3)
        end

        areas ./= areas[3][2]

        for i in 1:3
            println("Area Simple... (min., mean, max.) = ", areas[i])
        end

        #auc = PyRModule.auc_r(longitude, averages[1][bin_st:bin_end])
        #auc = PyRModule.auc_r(longitude, averages[2][bin_st:bin_end])
        #auc = PyRModule.auc_r(longitude, averages[3][bin_st:bin_end])
        #return

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        labels = ["\$ \\qquad \\;\\; {\\rm D} > 0\$ (N=$(size(profiles[1])[1]))", "\$\\qquad \\;\\; {\\rm D} < 0 \$ (N=$(size(profiles[2])[1]))", "\$0.83 > {\\rm D} > -0.80 \$ (N=$(size(profiles[3])[1]))", "\$\\quad 0.5\\;\\; > {\\rm D} > \\quad \\; 0 \$", "\$ \\quad 0 \\;\\;  > {\\rm D} > \\;-0.5 \$", "\$\\!-0.5 \\;\\; > {\\rm D} > --0.73 \$"]


        figure(figsize=(3.14961, 1.946563744), frameon=true)  # 8cm x 4.94427191 cm (golden)
        subplots_adjust(left=0.17, bottom=0.19, right=0.99, top=0.99, wspace=0., hspace=0.)

        minorticks_on()
        plot(longitude, averages[1][bin_st:bin_end], c="C1", label=labels[1], lw=0.3, zorder=200)
        plot(longitude, averages[2][bin_st:bin_end], c="C2", label=labels[2], lw=0.3, zorder=201)
        plot(longitude, averages[3][bin_st:bin_end], c="black", label=labels[3], lw=0.5, zorder=199, alpha=0.9)
        #yticks(collect(0.0:0.2:1)) # TODO
        #xlim(longitude[1], longitude[end])
        xlabel("longitude (deg.)")
        #ylabel("\$\\sqrt{\\rm Number \\; of \\; pulses}\$")
        #ylabel("Intensity (arbitrary units)")
        ylabel("Intensity (a. u.)")
        legend(prop=Dict("size"=> 4.1))
        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_averages3.pdf")
        savefig("$outdir/$(name_mod)_averages3.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()

    end


    """ calculates average profiels in two ranges - for different ranges tests"""
    function average_J1750_4(datas, outdir; lambda=1000.0, bin_st=nothing, bin_end=nothing, name_mod="0", show_=false)

        nums = []
        bins = []
        for data in datas
            nu, bi = size(data)
            push!(nums, nu)
            push!(bins, bi)
        end
        # take full but plot only part
        bin_st0 = 1
        bin_end0 = bins[1]
        das = [] # single pulse data
        for data in datas
            da = data[:, bin_st0:bin_end0]
            push!(das, da)
        end

        # drift rate data
        pulses = []
        ysps = []  # drift rate

        for i in 1:length(das)
            tracks1, lines1, inclines1, dr1, ysp1 = Tools.get_driftrate("$outdir/tracks/$i", lambda)
            (ysp, dr) = Tools.remove_duplicates_ysp(ysp1, dr1)
            push!(pulses, dr)
            push!(ysps, ysp)
        end

        ranges = [(2, 0.2), (-0.0, -1), (2, -2)]
        profiles = []
        averages = []

        for r in ranges
            push!(profiles, [])
            for i in 1:length(das)
                for (j, driftrate) in enumerate(ysps[i])
                    if (driftrate < r[1]) && (driftrate > r[2])
                        push!(profiles[end], das[i][pulses[i][j], :] )
                    end
                end
            end
        end

        # converting profiles TODO why why why?
        for ii in 1:length(profiles)
            println(ii)
            x, = size(profiles[ii])
            y, = size(profiles[ii][1]) # bins
            profs = zeros((x,y))
            for i in 1:x
                for j in 1:y
                    profs[i,j] = profiles[ii][i][j]
                end
            end
            profiles[ii] = profs
        end

        nums = []
        for profile in profiles
            push!(nums, size(profile)[1])
            push!(averages, Tools.average_profile(profile; norm=false))
            #println(size(profile)[1])
        end

        # around zero
        for i in 1:3
            me = mean(averages[i][700:1000])
            averages[i] .-= me
            averages[i] ./= size(profiles[i])[1] # multiply based on number of pulses
        end

        # Pulse longitude (raw)
        longitude = collect(range(bin_st/bins[1] * 360, bin_end/bins[1] * 360, length=bin_end - bin_st + 1))

        stds = []
        ma = maximum(averages[3])
        for i in 1:3
            averages[i] ./= ma
            st = std(averages[i][700:1000])
            push!(stds, st)
        end

        # random area calculations
        repeat = 10000
        areas = []

        for i in 1:3
            push!(areas, [])
            for j in 1:repeat
                # randomize signal based on std/rms
                av = averages[i][bin_st:bin_end]
                for k in 1:size(av)[1]
                    d = Normal(av[k], stds[i])
                    av[k] = rand(d)
                end
                auc = PyRModule.auc(longitude, av)
                push!(areas[end], auc)
            end
        end
        # normalize areas
        meme = mean(areas[3])
        areas ./= meme
        for i in 1:3
            println("Area = ", mean(areas[i]), " ± ", std(areas[i]))
        end

        # simple error calculation (overestimate?)
        areas = []
        for i in 1:3
            push!(areas, [])
            # randomize signal based on std/rms
            av = averages[i][bin_st:bin_end]
            auc1 = PyRModule.auc(longitude, av.-stds[i])
            auc2 = PyRModule.auc(longitude, av)
            auc3 = PyRModule.auc(longitude, av.+stds[i])
            push!(areas[end], auc1)
            push!(areas[end], auc2)
            push!(areas[end], auc3)
        end

        areas ./= areas[3][2]

        for i in 1:3
            println("Area Simple... (min., mean, max.) = ", areas[i])
        end

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        labels = ["\$ \\qquad \\;\\; {\\rm D} > 0\$ (N=$(size(profiles[1])[1]))", "\$\\qquad \\;\\; {\\rm D} < 0 \$ (N=$(size(profiles[2])[1]))", "\$0.91 > {\\rm D} > -0.73 \$ (N=$(size(profiles[3])[1]))", "\$\\quad 0.5\\;\\; > {\\rm D} > \\quad \\; 0 \$", "\$ \\quad 0 \\;\\;  > {\\rm D} > \\;-0.5 \$", "\$\\!-0.5 \\;\\; > {\\rm D} > --0.73 \$"]


        figure(figsize=(3.14961, 1.946563744), frameon=true)  # 8cm x 4.94427191 cm (golden)
        subplots_adjust(left=0.17, bottom=0.19, right=0.99, top=0.99, wspace=0., hspace=0.)

        minorticks_on()
        plot(longitude, averages[1][bin_st:bin_end], c="C1", label=labels[1], lw=0.3, zorder=200)
        plot(longitude, averages[2][bin_st:bin_end], c="C2", label=labels[2], lw=0.3, zorder=201)
        plot(longitude, averages[3][bin_st:bin_end], c="black", label=labels[3], lw=0.5, zorder=199, alpha=0.9)
        #yticks(collect(0.0:0.2:1)) # TODO
        #xlim(longitude[1], longitude[end])
        xlabel("longitude (deg.)")
        #ylabel("\$\\sqrt{\\rm Number \\; of \\; pulses}\$")
        #ylabel("Intensity (arbitrary units)")
        ylabel("Intensity (a. u.)")
        legend(prop=Dict("size"=> 4.1))
        #tick_params(labeltop=false, labelbottom=true)
        println("$outdir/$(name_mod)_averages4.pdf")
        savefig("$outdir/$(name_mod)_averages4.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()

    end



    """ checks average profiels stability """
    function average_J1750_stability(datas, outdir; lambda=1000.0, bin_st=nothing, bin_end=nothing, name_mod="0", show_=false)

        nums = []
        bins = []
        for data in datas
            nu, bi = size(data)
            push!(nums, nu)
            push!(bins, bi)
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins[1] end
        das = [] # single pulse data
        for data in datas
            da = data[:, bin_st:bin_end]
            push!(das, da)
        end

        # drift rate data
        pulses = []
        ysps = []  # drift rate

        for i in 1:length(das)
            tracks1, lines1, inclines1, dr1, ysp1 = Tools.get_driftrate("$outdir/tracks/$i", lambda)
            (ysp, dr) = Tools.remove_duplicates_ysp(ysp1, dr1)
            push!(pulses, dr)
            push!(ysps, ysp)
        end

        ranges = [(1, 0), (0, -1), (2, -2)]
        profiles = []
        averages = []

        for r in ranges
            push!(profiles, [])
            for i in 1:length(das)
                for (j, driftrate) in enumerate(ysps[i])
                    if (driftrate < r[1]) && (driftrate > r[2])
                        push!(profiles[end], das[i][pulses[i][j], :] )
                        #println(pulses[i][j])
                    end
                end
            end
        end

        # converting profiles TODO why why why?
        for ii in 1:length(profiles)
            x, = size(profiles[ii])
            y, = size(profiles[ii][1]) # bins
            profs = zeros((x,y))
            for i in 1:x
                for j in 1:y
                    profs[i,j] = profiles[ii][i][j]
                end
            end
            profiles[ii] = profs
        end

        binoff_st = 10
        binoff_end = 310
        binon_st = 350
        binon_end = 650

        nums = []
        for profile in profiles
            push!(nums, size(profile)[1])
            #push!(averages, Tools.average_profile(profile; norm=true))
            push!(averages, Tools.average_profile(profile; norm=[binoff_st, binoff_end]))
            #println(size(profile)[1])
        end

        # Pulse longitude (raw)
        longitude = collect(range(bin_st/bins[1] * 360, bin_end/bins[1] * 360, length=bin_end - bin_st + 1))

        bin_num = bins[1]
        # calculate stability for the whole observation
        onemrhos, onemrhose, npulses = Tools.analyse_profile_stability(profiles[3], averages[3], bin_num, binoff_st, binoff_end, binon_st, binon_end, longitude)
        onemrhos1, onemrhose1, npulses1 = Tools.analyse_profile_stability(profiles[1], averages[1], bin_num, binoff_st, binoff_end, binon_st, binon_end, longitude)
        onemrhos2, onemrhose2, npulses2 = Tools.analyse_profile_stability(profiles[2], averages[2], bin_num, binoff_st, binoff_end, binon_st, binon_end, longitude)


        lines = []
        stabilization_time = 1
        stabilization_time2 = 1
        thresh = 0.0005

        # positive
        @. f(x, p) = p[4] * x^3 + p[1] * x ^ 2 + p[2] * x + p[3]
        #p0 = [-0.5, 0, 1]
        p0 = [-0.23506367429667388, 0.12580584169631978, -0.14562974681462082, -0.1330915339905905]
        onemrhos_ = convert(Array{Float64,1}, onemrhos1[2:end]) # skipping first
        npulses_ = convert(Array{Float64,1}, npulses1[2:end])
        x_ = log10.(npulses_)
        y_ = log10.(onemrhos_)

        fit = curve_fit(f, x_, y_, p0)
        p = coef(fit)
        err = stderror(fit)
        x_ = collect(0:0.01:4.3)
        line = [10 .^ x_, 10 .^ f(x_, p)]
        push!(lines, line)
        #println(p)
        #println(err)
        for i in 1:size(line[1])[1]
            #println(line[2][i], " $i")
            if line[2][i] < thresh
                stabilization_time = line[1][i]
                println("Stabilization time: $stabilization_time")
                break
            end
        end

        # negative
        @. f(x, p) = p[4] * x^3 + p[1] * x ^ 2 + p[2] * x + p[3]
        #p0 = [-0.5, 0, 1]
        p0 = [-0.23506367429667388, 0.12580584169631978, -0.14562974681462082, -0.1330915339905905]
        onemrhos_ = convert(Array{Float64,1}, onemrhos2[1:end]) # skipping first
        npulses_ = convert(Array{Float64,1}, npulses2[1:end])
        x_ = log10.(npulses_)
        y_ = log10.(onemrhos_)

        fit = curve_fit(f, x_, y_, p0)
        p = coef(fit)
        err = stderror(fit)
        x_ = collect(0:0.01:4.3)
        line2 = [10 .^ x_, 10 .^ f(x_, p)]
        push!(lines, line2)
        #println(p)
        #println(err)
        for i in 1:size(line2[1])[1]
            #println(line[2][i], " $i")
            if line2[2][i] < thresh
                stabilization_time2 = line2[1][i]
                println("Stabilization time 2: $stabilization_time2")
                break
            end
        end

        #=
        # line fitting (not so good)
        @. f(x, p) = p[2] * x + p[1]
        # why why why?
        onemrhos_ = convert(Array{Float64,1}, onemrhos[9:end]) # skipping some
        npulses_ = convert(Array{Float64,1}, npulses[9:end])
        x_ = log10.(npulses_)
        y_ = log10.(onemrhos_)
        data = DataFrame(X=x_, Y=y_)
        ols = lm(@formula(Y ~ X), data)
        x_ = collect(2.2:0.1:4.3)
        line = [10 .^ x_, 10 .^ f(x_, coef(ols))]
        push!(lines, line)
        for i in 1:size(line[1])[1]
            println(line[2][i], " $i")
        onemrhos_ = convert(Array{Float64,1}, onemrhos[9:end]) # skipping some
        npulses_ = convert(Array{Float64,1}, npulses[9:end])
        x_ = log10.(npulses_)
        y_ = log10.(onemrhos_)
        data = DataFrame(X=x_, Y=y_)
        ols = lm(@formula(Y ~ X), data)
        x_ = collect(2.2:0.1:4.3)
        line = [10 .^ x_, 10 .^ f(x_, coef(ols))]
        push!(lines, line)
        for i in 1:size(line[1])[1]
            println(line[2][i], " $i")
            if line[2][i] < 0.0005
                stabilization_time = line[1][i]
                println("Stabilization time: $stabilization_time")
                break
            end
        end
        =#

        labels = ["\$ \\qquad \\;\\; {\\rm D} > 0\$ (N=$(size(profiles[1])[1]))", "\$\\qquad \\;\\; {\\rm D} < 0 \$ (N=$(size(profiles[2])[1]))", "\$0.91 > {\\rm D} > -0.73 \$ (N=$(size(profiles[3])[1]))", "\$\\quad 0.5\\;\\; > {\\rm D} > \\quad \\; 0 \$", "\$ \\quad 0 \\;\\;  > {\\rm D} > \\;-0.5 \$", "\$\\!-0.5 \\;\\; > {\\rm D} > --0.73 \$"]

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 1.946563744), frameon=true)  # 8cm x 4.94427191 cm (golden)
        subplots_adjust(left=0.19, bottom=0.21, right=0.99, top=0.99, wspace=0., hspace=0.)
        minorticks_on()

        errorbar(npulses, onemrhos, yerr=onemrhose, color="black", lw=1.5, alpha=0.79, label=labels[3])
        errorbar(npulses1, onemrhos1, yerr=onemrhose1, color="tab:orange", lw=1.5, alpha=0.79, label=labels[1])
        errorbar(npulses2, onemrhos2, yerr=onemrhose2, color="tab:green", lw=1.5, alpha=0.79, label=labels[2])
        plot(lines[1][1], lines[1][2], lw=0.7, color="tab:orange", ls="--")
        plot(lines[2][1], lines[2][2], lw=0.7, color="tab:green", ls="--")
        #axvline(x = stabilization_time)
        axhline(y = thresh, ls="--", color="black")
        #semilogx()

        legend(prop=Dict("size"=> 4.1))
        loglog()
        xlim([0.7, 9e3])
        ylim([9e-5, 1.4])
        xlabel("Number of pulses")
        ylabel("\$1-\\rho\$")
        println("$outdir/$(name_mod)_averages_stability.pdf")
        savefig("$outdir/$(name_mod)_averages_stability.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()

    end





    """
    get data for driftdirection_J1750 plot
    """
    function getdata_driftdirection(ii, jj, selected_tracks, fitted_lines)
        x_ = []
        y_ = []
        ysl_ = []
        y2_ = []
        slopes = []
        eslopes = []
        xlin = []
        exlin = []
        for i in ii
            for j in jj
                push!(x_, selected_tracks[i][j][1])
                push!(y_, selected_tracks[i][j][2])
                #println("ysl (fl)", size(fitted_lines[i][j][3][1]))
                push!(y2_, fitted_lines[i][j][2][1]) # WHY [1]?
                push!(ysl_, fitted_lines[i][j][3][1]) # WHY [1]?
                push!(slopes, fitted_lines[i][j][4][1]) # WHY [1]?
                push!(eslopes, fitted_lines[i][j][5][1]) # WHY [1]?
                push!(xlin, fitted_lines[i][j][8][1]) # WHY [1]?
                push!(exlin, fitted_lines[i][j][9][1]) # WHY [1]?
            end
        end
        bps = []
        ebps = []
        # add only once
        for j in jj
            push!(bps, fitted_lines[ii[1]][j][6][1]) # WHY [1]?
            push!(ebps, fitted_lines[ii[1]][j][7][1]) # WHY [1]?
        end

        return x_, y_, ysl_, y2_, slopes, eslopes, xlin, exlin, bps, ebps
    end


    function driftdirection_J1750(datas, outdir; lambda=1000.0, bin_st=nothing, bin_end=nothing, name_mod="0", show_=false)

        nums = []
        bins = []
        for data in datas
            nu, bi = size(data)
            push!(nums, nu)
            push!(bins, bi)
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins[1] end
        das = [] # single pulse data
        for data in datas
            da = data[:, bin_st:bin_end]
            push!(das, da)
        end

        tracks = []
        inclines = []
        dr_x = []
        dr_y = []
        for i in 1:6
            tracks1, lines1, inclines1, dr1, ysp1, dof1 = Tools.get_driftrate("$outdir/tracks$i", lambda)
            push!(tracks, tracks1)
            push!(inclines, inclines1)
            push!(dr_x, dr1)
            push!(dr_y, ysp1)
        end

        ranges = Dict(2=>(510, 589), 3=>(50, 150), 3=>(375, 445), 4=>(100, 150), 5=>(250, 350), 6=>(200, 365))
        #ranges = Dict(2=>(510, 589), 3=>(50, 150), 3=>(375, 445), 4=>(85, 150), 5=>(250, 350), 6=>(200, 365))
        #ranges = Dict(2=>(500, 590))

        selected_tracks = [[] for i in 1:6] # all six slots (first will be empty)
        #println(tracks[1][1][:, 2])  # pulse number
        for (k, v) in ranges
            println("Session: ", k)
            for i in 1:length(tracks[k])
                create_new = true
                for (ii, pulse) in enumerate(tracks[k][i][:, 2])
                    if (pulse > v[1]) && (pulse < v[2])
                        if create_new
                            push!(selected_tracks[k], [[], []])
                            create_new = false
                            println("\tNew track session:$k track:$i pulse:$pulse loc:", tracks[k][i][ii, 1])
                        end
                        push!(selected_tracks[k][end][1], pulse) # x - pulse number
                        push!(selected_tracks[k][end][2], tracks[k][i][ii, 1]) # location
                        #println("\t $v $pulse ", tracks[k][i][ii, 1])
                    end
                end
            end
        end

        #npsi = [3, 3, 2, 2, 1, 2, 1, 1, 2, 2, 1, 4, 1, 4]
        npsi = [3, 3, 2, 2, 1, 2, 1, 1, 2, 2, 1, 4, 1, 4]
        nn = 1
        fitted_lines = [[] for i in 1:6] # all six slots (first will be empty)
        for i in 1:length(selected_tracks)
            for j in 1:length(selected_tracks[i])
                push!(fitted_lines[i], [[], [], [], [], [], [], [], [], []])
                println("$i $j nn:$nn")
                x = convert(Array{Float64,1}, selected_tracks[i][j][1])
                y = convert(Array{Float64,1}, selected_tracks[i][j][2])
                push!(fitted_lines[i][end][1], x)
                println("nn $nn")
                x, y2, ysl, slopes, eslopes, bp, ebp, xl, exl = PyRModule.segmented(x, y, npsi[nn]; lambda=lambda)
                nn += 1
                push!(fitted_lines[i][end][2], y2)
                push!(fitted_lines[i][end][3], ysl)  # for tracks
                push!(fitted_lines[i][end][4], slopes)
                push!(fitted_lines[i][end][5], eslopes)
                push!(fitted_lines[i][end][6], bp)
                push!(fitted_lines[i][end][7], ebp)
                push!(fitted_lines[i][end][8], xl) # xlin later on
                push!(fitted_lines[i][end][9], exl) # exlin later on
            end
        end

        # NOPE
        #=
        p0 = [510, 530, 550, 570, 500, 1, 500, -1, 500, 1, 500, -1]
        println(PyModule.least_sq(x, y, p0; num=1, show_=true))
        @. model(x, p) = p[1] + p[2] * x
        fi = curve_fit(model, x, y, p0)
        #	dof(fit): degrees of freedom
        #	coef(fit): best fit parameters
        #	fit.resid: residuals = vector of residuals
        #	fit.jacobian: estimated Jacobian at solution
        p = coef(fi)
        err = stderror(fi)
        =#

        iis = [[2, 2], [3, 3], [4, 4], [5, 5, 5], [6, 6]]  # obs. session
        jjs = [[1, 2], [1, 2], [1, 2], [1, 3, 4], [2, 4]]  # tracks

        #=
        ddx = 0
        for wh in 1:length(iis)
            x_, y_, ysl_, y2_, slopes, eslopes, xlin, exlin, bps, ebps = getdata_driftdirection(iis[wh], jjs[wh], selected_tracks, fitted_lines)

            mx = 0
            for i in 1:length(x_)
                dx = maximum(x_[i]) - minimum(x_[i])
                if dx > mx
                    mx = dx
                end
            end
            ddx += mx
            println(wh, " dx: ", mx)

            #println(minimum(x_), " ", maximum(x_))
        end
        println(" ddx: ", ddx)
        println(" ddx/2: ", ddx/2)
        return
        =#

        #println(y2_)
        #println(size(y2_))

        #name_mod = "$(iis[1][1])"
        name_mod = "23456"

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(7.086614, 6.299213))  # 18 cm x 16 cm

        subplots_adjust(left=0.08, bottom=0.07, right=0.95, top=0.92, wspace=0.0, hspace=0.0)
        figtext(0.09, 0.9, "a)", size=10)
        figtext(0.5, 0.9, "b)", size=10)
        figtext(0.84, 0.9, "c)", size=10)
        figtext(0.09, 0.44, "d)", size=10)
        figtext(0.43, 0.44, "e)", size=10)
        figtext(0.74, 0.44, "f)", size=10)

        driftdirection_J1750_subplot(1, 0, 0, 77, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        driftdirection_J1750_subplot(2, 0, 120, 68, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        driftdirection_J1750_subplot(3, 0, 210, 48, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        driftdirection_J1750_subplot(4, 35, 0, 98, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        driftdirection_J1750_subplot(5, 35, 100, 163, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)

        savefig("$outdir/$(name_mod)_driftdirection.pdf")
        println("$outdir/$(name_mod)_driftdirection.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
    end


    function driftdirection_J1750_subplot(which, row, col, colspan, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)

        x_, y_, ysl_, y2_, slopes, eslopes, xlin, exlin, bps, ebps = getdata_driftdirection(iis[which], jjs[which], selected_tracks, fitted_lines)

        ax = subplot2grid((65, 270), (row, col), colspan=colspan, rowspan=10)  # row column
        minorticks_on()
        for i in 1:length(x_)
            #scatter(x_[i], y_[i])
            plot(x_[i], y_[i], marker="x", color="red", markersize=2.5, lw=0)
            plot(x_[i], ysl_[i], lw=1.5, alpha=0.7, c="C2")
            plot(x_[i], y2_[i], lw=1.0, alpha=0.7, c="C1")
        end
        yl = ylim()
        #println(bps)
        for i in 1:length(bps)
            for j in 1:length(bps[i])
                if i == 1
                    color = "red"
                else
                    color = "blue"
                end
                axvline(x=bps[i][j], ls="--", c=color, lw=0.3)
                rectangle = patch.Rectangle((bps[i][j]-ebps[i][j], yl[1]), width=2*ebps[i][j], height=yl[2]-yl[1], fc=color, alpha=0.5)
                #println((yl[1], bps[i][j]-ebps[i][j])," ",  2*ebps[i][j], " ", yl[2]-yl[1])
                ax.add_patch(rectangle)
            end
        end
        xl = xlim()

        ax = subplot2grid((65, 270), (row+10, col), colspan=colspan, rowspan=10)  # row column
        minorticks_on()
        for i in 1:length(slopes)
            errorbar(xlin[i], slopes[i], yerr=eslopes[i], xerr=exlin[i], color="none", lw=0.5, marker="_", mec="grey", ecolor="grey", capsize=0, mfc="grey", ms=1.0)
        end
        plot(dr_x[iis[which][1]], dr_y[iis[which][1]])  # smoothed drift rate
        xlim(xl)
        ylim([-2.3, 2.3])

        ax = subplot2grid((65, 270), (row+20, col), colspan=colspan, rowspan=10)  # row column
        minorticks_on()
        for i in 1:length(x_)
            residuals1 = Tools.residuals(y_[i], ysl_[i])
            residuals2 = Tools.residuals(y_[i], y2_[i])
            scatter(x_[i], residuals1, marker="+", color="C2", s=4, lw=0.5, alpha=0.5) # green
            scatter(x_[i], residuals2, marker="o", color="C1", fc="none", s=4, lw=0.5, alpha=0.5) # yellow
            chi1 = Tools.chisquare(y_[i], ysl_[i])
            chi2 = Tools.chisquare(y_[i], y2_[i])
            rs = Tools.rsquared(y_[i], ysl_[i])
            rs2 = Tools.rsquared(y_[i], y2_[i])
            println("chi1 ", chi1)
            println("chi2 ", chi2)
            println("rs ", rs)
            println("rs2 ", rs2)
            println("AAA")
        end
    end


    function driftdirection_J1750_2(datas, outdir; lambda=1000.0, bin_st=nothing, bin_end=nothing, name_mod="0", show_=false)

        nums = []
        bins = []
        for data in datas
            nu, bi = size(data)
            push!(nums, nu)
            push!(bins, bi)
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins[1] end
        das = [] # single pulse data
        for data in datas
            da = data[:, bin_st:bin_end]
            push!(das, da)
        end

        tracks = []
        inclines = []
        dr_x = []
        dr_y = []
        for i in 1:7
            tracks1, lines1, inclines1, dr1, ysp1, dof1 = Tools.get_driftrate("$outdir/tracks/$i", lambda)
            push!(tracks, tracks1)
            push!(inclines, inclines1)
            push!(dr_x, dr1)
            push!(dr_y, ysp1)
        end

        ranges = [[2, (300, 350)], [2, (375, 420)], [2, (485, 595)], [2, (955, 1040)], [3, (50, 350)], [4, (1, 450)]]#, [4, (100, 150)], [5, (250, 350)], [6, (200, 365)]]

        selected_tracks = [[] for i in 1:7] # all seven slots (first will be empty)
        #println(tracks[1][1][:, 2])  # pulse number

        for ra in ranges
            k = ra[1]
            v = ra[2]
            println("Session: ", k)
            for i in 1:length(tracks[k])
                create_new = true
                for (ii, pulse) in enumerate(tracks[k][i][:, 2])
                    if (pulse > v[1]) && (pulse < v[2])
                        if create_new
                            push!(selected_tracks[k], [[], []])
                            create_new = false
                            println("\tNew track session:$k track:$i pulse:$pulse loc:", tracks[k][i][ii, 1])
                        end
                        push!(selected_tracks[k][end][1], pulse) # x - pulse number
                        push!(selected_tracks[k][end][2], tracks[k][i][ii, 1]) # location
                        #println("\t $v $pulse ", tracks[k][i][ii, 1])
                    end
                end
            end
        end


        #npsi = [3, 3, 2, 2, 1, 2, 1, 1, 2, 2, 1, 4, 1, 4]
        npsi = [1, 1,   1, 1, 1,   0,4, 4,   2, 2,   2, 4, 4, 1,    2, 2, 2, 2,   0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0] # 0 is for skipped track
        # bad fit for 2 4 3 1 (4th)
        nn = 1
        fitted_lines = [[] for i in 1:7] # all seven slots (first will be empty)
        for i in 1:length(selected_tracks)
            for j in 1:length(selected_tracks[i])
                push!(fitted_lines[i], [[], [], [], [], [], [], [], [], []])
                println("$i $j nn:$nn")
                x = convert(Array{Float64,1}, selected_tracks[i][j][1])
                y = convert(Array{Float64,1}, selected_tracks[i][j][2])
                push!(fitted_lines[i][end][1], x)
                println("nn $nn")
                println("npsi $(npsi[nn])")
                x, y2, ysl, slopes, eslopes, bp, ebp, xl, exl = PyRModule.segmented(x, y, npsi[nn]; lambda=lambda, preview=false)
                nn += 1
                push!(fitted_lines[i][end][2], y2)
                push!(fitted_lines[i][end][3], ysl)  # for tracks
                push!(fitted_lines[i][end][4], slopes)
                push!(fitted_lines[i][end][5], eslopes)
                push!(fitted_lines[i][end][6], bp)
                push!(fitted_lines[i][end][7], ebp)
                push!(fitted_lines[i][end][8], xl) # xlin later on
                push!(fitted_lines[i][end][9], exl) # exlin later on
            end
        end

        iis = [[2, 2], [2, 2, 2], [2, 2], [2, 2], [3, 3, 3, 3], [4, 4, 4],    [5, 5], [6, 6]]  # obs. session why iis?
        jjs = [[1, 2], [3, 4, 5], [7, 8], [9, 10], [1, 2, 3, 4], [1, 2, 3],    [1, 2], [1, 3, 4], [2, 4]]  # tracks

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(7.086614, 6.299213))  # 18 cm x 16 cm

        subplots_adjust(left=0.08, bottom=0.07, right=0.95, top=0.92, wspace=0.0, hspace=0.0)
        figtext(0.09, 0.9, "a)", size=10)
        figtext(0.5, 0.9, "b)", size=10)
        figtext(0.84, 0.9, "c)", size=10)
        figtext(0.09, 0.44, "d)", size=10)
        figtext(0.43, 0.44, "e)", size=10)
        figtext(0.74, 0.44, "f)", size=10)

        driftdirection_J1750_subplot(1, 0, 0, 77, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        driftdirection_J1750_subplot(2, 0, 120, 68, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        driftdirection_J1750_subplot(3, 0, 210, 48, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        driftdirection_J1750_subplot(4, 35, 0, 98, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        driftdirection_J1750_subplot(5, 35, 100, 163, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        driftdirection_J1750_subplot(6, 35, 100, 163, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)

        savefig("$outdir/$(name_mod)_driftdirection_2.pdf")
        println("$outdir/$(name_mod)_driftdirection_2.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
    end



    function driftdirection_J1750_subplot3(which, row, col, colspan, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=600)

        x_, y_, ysl_, y2_, slopes, eslopes, xlin, exlin, bps, ebps = getdata_driftdirection(iis[which], jjs[which], selected_tracks, fitted_lines)

        colors = ["tab:orange", "tab:green", "tab:purple", "tab:brown"]
        tracks_num = length(jjs[which])

        ax = subplot2grid((65, colnum), (row, col), colspan=colspan, rowspan=20)  # row column
        #tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        if which in (2, 3, 4, 5, 7, 8, 9)
            tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        end
        if which < 6
            ax.xaxis.set_label_position("top")
            ax.xaxis.set_ticks_position("top")
        else
            tick_params(axis="x", which="both", direction="out", labelbottom=false)
        end
        minorticks_on()
        for i in 1:length(x_)
            #scatter(x_[i], y_[i])
            plot(x_[i], y_[i], marker="x", color="tab:red", markersize=2.5, lw=0)
            # plot(x_[i], ysl_[i], lw=1.3, alpha=0.7, c="tab:blue") # not visible
            #plot(x_[i], y2_[i], lw=2.0, alpha=0.7, c=colors[i]) # why?
        end
        for i in 1:tracks_num #length(x_)
                plot(x_[i], y2_[i], lw=1.7, alpha=0.7, c=colors[i])
        end
        yl = (195, 146)
        #yl = ylim() # TODO just for now
        for i in 1:length(bps)
            for j in 1:length(bps[i])
                if ebps[i][j] < 10 # do not plot wrong estimates
                    axvline(x=bps[i][j], ls="--", c=colors[i], lw=0.3)
                    rectangle = patch.Rectangle((bps[i][j]-ebps[i][j], yl[1]), width=2*ebps[i][j], height=yl[2]-yl[1], fc=colors[i], alpha=0.1)
                    #println("ERR ", ebps[i][j])
                    #println((yl[1], bps[i][j]-ebps[i][j])," ",  2*ebps[i][j], " ", yl[2]-yl[1])
                    ax.add_patch(rectangle)
                end
            end
        end
        ylim(yl)
        xl = xlim()
        if (which == 1) || (which == 6)
            ylabel("Longitude (\$ ^{\\circ}\$)")
        end

        ax = subplot2grid((65, colnum), (row+20, col), colspan=colspan, rowspan=10)  # row column
        if which in (2, 3, 4, 5, 7, 8, 9)
            tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        end
        if which < 6
            tick_params(axis="x", which="both", direction="out", labelbottom=false)
        end
        minorticks_on()
        #println(size(slopes))
        for i in 1:tracks_num # length(slopes)
            errorbar(xlin[i], slopes[i], yerr=eslopes[i], xerr=exlin[i], color="none", lw=0.5, marker="_", mec=colors[i], ecolor=colors[i], capsize=0, mfc=colors[i], ms=1.0)
        end
        plot(dr_x[iis[which][1]], dr_y[iis[which][1]], color="tab:blue")  # smoothed drift rate
        xlim(xl)
        ylim([-0.999, 0.999])
        if (which == 1) || (which == 6)
            ylabel("Drift rate \$(^\\circ / P)\$")
        end

        # report on fitting, plotting turned off
        #ax = subplot2grid((65, 270), (row+20, col), colspan=colspan, rowspan=10)  # row column
        #minorticks_on()
        rs_ = []
        rs2_ = []
        for i in 1:length(x_)
            #residuals1 = Tools.residuals(y_[i], ysl_[i])
            #residuals2 = Tools.residuals(y_[i], y2_[i])
            #scatter(x_[i], residuals1, marker="+", color="C2", s=4, lw=0.5, alpha=0.5) # green
            #scatter(x_[i], residuals2, marker="o", color="C1", fc="none", s=4, lw=0.5, alpha=0.5) # yellow
            chi1 = Tools.chisquare(y_[i], ysl_[i])
            chi2 = Tools.chisquare(y_[i], y2_[i])
            rs = Tools.rsquared(y_[i], ysl_[i])
            rs2 = Tools.rsquared(y_[i], y2_[i])
            #println("chi1 ", chi1)
            #println("chi2 ", chi2)
            push!(rs_, rs)
            push!(rs2_, rs2)
            println("Panel: $which")
            println("\tSpline: r^2 ", rs)
            println("\tLinear reg.: r^2 ", rs2)
        end
        return rs_, rs2_
    end


    function driftdirection_J1750_3(datas, outdir; lambda=1000.0, bin_st=nothing, bin_end=nothing, name_mod="0", show_=false)

        nums = []
        bins = []
        for data in datas
            nu, bi = size(data)
            push!(nums, nu)
            push!(bins, bi)
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins[1] end
        das = [] # single pulse data
        for data in datas
            da = data[:, bin_st:bin_end]
            push!(das, da)
        end

        tracks = []
        inclines = []
        dr_x = []
        dr_y = []
        for i in 1:7
            tracks1, lines1, inclines1, dr1, ysp1, dof1 = Tools.get_driftrate("$outdir/tracks/$i", lambda)
            push!(tracks, tracks1)
            push!(inclines, inclines1)
            push!(dr_x, dr1)
            push!(dr_y, ysp1)
        end

        ranges = [[2, (477, 600)], [2, (955, 1050)], [3, (50, 170)], [3, (380, 470)], [4, (85, 144)], [5, (240, 350)], [6, (200, 400)], [7, (0, 80)], [7, (125, 275)]]
        # rejected [2, (825, 922)] (npsi, 0, 2, 2)

        sum = 0
        i = 0
        for r in ranges
            i += 1
            pulses = r[2][2] - r[2][1]
            sum += pulses
            println("$i $pulses")
        end
        println("Sum: $sum")


        selected_tracks = [[] for i in 1:7] # all seven slots (first will be empty)
        #println(tracks[1][1][:, 2])  # pulse number

        for ra in ranges
            k = ra[1]
            v = ra[2]
            println("Session: ", k)
            for i in 1:length(tracks[k])
                create_new = true
                for (ii, pulse) in enumerate(tracks[k][i][:, 2])
                    if (pulse > v[1]) && (pulse < v[2])
                        if create_new
                            push!(selected_tracks[k], [[], []])
                            create_new = false
                            println("\tNew track session:$k track:$i pulse:$pulse loc:", tracks[k][i][ii, 1])
                        end
                        push!(selected_tracks[k][end][1], pulse) # x - pulse number
                        push!(selected_tracks[k][end][2], tracks[k][i][ii, 1]) # location
                        #println("\t $v $pulse ", tracks[k][i][ii, 1])
                    end
                end
            end
        end

        npsi = [0, 6, 6,   3, 3,   2, 2, 1,   2, 2,   2, 2,   2, 2, 2,    2, 6, 7,   2, 2,    0, 2, 2, 1, 0 ] # 0 is for skipped track
        nn = 1
        fitted_lines = [[] for i in 1:7] # all seven slots (first will be empty)
        for i in 1:length(selected_tracks)
            println("\nSession: $i")
            for j in 1:length(selected_tracks[i])
                push!(fitted_lines[i], [[], [], [], [], [], [], [], [], []])
                println("\tTrack: $j nn:$nn npsi $(npsi[nn])")
                x = convert(Array{Float64,1}, selected_tracks[i][j][1])
                y = convert(Array{Float64,1}, selected_tracks[i][j][2])
                push!(fitted_lines[i][end][1], x)
                x, y2, ysl, slopes, eslopes, bp, ebp, xl, exl = PyRModule.segmented(x, y, npsi[nn]; lambda=lambda, preview=false)
                nn += 1
                push!(fitted_lines[i][end][2], y2)
                push!(fitted_lines[i][end][3], ysl)  # for tracks
                push!(fitted_lines[i][end][4], slopes)
                push!(fitted_lines[i][end][5], eslopes)
                push!(fitted_lines[i][end][6], bp)
                push!(fitted_lines[i][end][7], ebp)
                push!(fitted_lines[i][end][8], xl) # xlin later on
                push!(fitted_lines[i][end][9], exl) # exlin later on
            end
        end


        iis = [[2], [2], [3], [3], [4], [5], [6], [7], [7]]  # obs. session why iis?
        jjs = [[2, 3], [4, 5],   [1, 2, 3], [4, 5],   [1, 2],   [1, 2, 3],   [1, 2, 3],    [1, 2],   [4, 5, 6]]  # tracks

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
        subplots_adjust(left=0.08, bottom=0.09, right=0.99, top=0.92, wspace=0.0, hspace=0.0)
        figtext(0.09, 0.89, "a)", size=10)
        figtext(0.32, 0.89, "b)", size=10)
        figtext(0.51, 0.89, "c)", size=10)
        figtext(0.73, 0.89, "d)", size=10)
        figtext(0.905, 0.89, "e)", size=10)
        figtext(0.085, 0.425, "f)", size=10)
        figtext(0.55, 0.43, "g)", size=10)
        figtext(0.615, 0.43, "h)", size=10)
        figtext(0.97, 0.43, "i)", size=10)

        (r1,r11) = driftdirection_J1750_subplot3(1, 0, 0, 123, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        (r2, r22) = driftdirection_J1750_subplot3(2, 0, 151, 95, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        (r3, r33) = driftdirection_J1750_subplot3(3, 0, 274, 120, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        (r4, r44) = driftdirection_J1750_subplot3(4, 0, 422, 90, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        (r5, r55) = driftdirection_J1750_subplot3(5, 0, 540, 59, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        (r6, r66) = driftdirection_J1750_subplot3(6, 35, 0, 110, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        (r7, r77) = driftdirection_J1750_subplot3(7, 35, 130, 200, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        (r8, r88) = driftdirection_J1750_subplot3(8, 35, 350, 80, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        (r9, r99) = driftdirection_J1750_subplot3(9, 35, 450, 150, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y)
        figtext(0.5, 0.97, "Pulse number", size=7, ha="center")
        figtext(0.5, 0.01, "Pulse number", size=7, ha="center")

        savefig("$outdir/$(name_mod)_driftdirection_3.pdf")
        println("$outdir/$(name_mod)_driftdirection_3.pdf")

        r = vcat(r1, r2, r3, r4, r5, r6, r7, r8, r9)
        rr = vcat(r11, r22, r33, r44, r55, r66, r77, r88, r99)
        println("Spline R^2: ", mean(r), " ", std(r))
        println("Linear reg. R^2: ", mean(rr), " ", std(rr))
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
    end


    function driftdirection_J1750_3b(datas, outdir; lambda=1000.0, bin_st=nothing, bin_end=nothing, name_mod="0", show_=false)

        nums = []
        bins = []
        for data in datas
            nu, bi = size(data)
            push!(nums, nu)
            push!(bins, bi)
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins[1] end
        das = [] # single pulse data
        for data in datas
            da = data[:, bin_st:bin_end]
            push!(das, da)
        end

        tracks = []
        inclines = []
        dr_x = []
        dr_y = []
        for i in 1:7
            tracks1, lines1, inclines1, dr1, ysp1, dof1 = Tools.get_driftrate("$outdir/tracks/$i", lambda)
            push!(tracks, tracks1)
            push!(inclines, inclines1)
            push!(dr_x, dr1)
            push!(dr_y, ysp1)
        end

        #ranges = [[2, (477, 600)], [2, (955, 1050)], [3, (50, 170)], [3, (380, 470)], [4, (85, 144)], [5, (240, 350)], [6, (200, 400)], [7, (0, 80)], [7, (125, 275)]]
        ranges = [[2, (0, 70)], [2, (500, 610)], [2, (630, 700)], [2, (960, 1020)], [4, (85, 140)], [5, (295, 355)], [6, (200, 350)], [7, (0, 90)], [7, (190, 270)]]
        # rejected [2, (825, 922)] (npsi, 0, 2, 2)

        sum = 0
        i = 0
        for r in ranges
            i += 1
            pulses = r[2][2] - r[2][1]
            sum += pulses
            println("$i $pulses")
        end
        println("Sum: $sum")

        selected_tracks = [[] for i in 1:7] # all seven slots (first will be empty)
        #println(tracks[1][1][:, 2])  # pulse number

        for ra in ranges
            k = ra[1]
            v = ra[2]
            println("Session: ", k)
            for i in 1:length(tracks[k])
                create_new = true
                for (ii, pulse) in enumerate(tracks[k][i][:, 2])
                    if (pulse > v[1]) && (pulse < v[2])
                        if create_new
                            push!(selected_tracks[k], [[], []])
                            create_new = false
                            println("\tNew track session:$k track:$i pulse:$pulse loc:", tracks[k][i][ii, 1])
                        end
                        push!(selected_tracks[k][end][1], pulse) # x - pulse number
                        push!(selected_tracks[k][end][2], tracks[k][i][ii, 1]) # location
                        #println("\t $v $pulse ", tracks[k][i][ii, 1])
                    end
                end
            end
        end

        #npsi = [0, 6, 6,   3, 3,   2, 2, 1,   2, 2,   2, 2,   2, 2, 2,    2, 6, 7,   2, 2,    0, 2, 2, 1, 0 ] # 0 is for skipped track
        npsi = [0, 0, 2, 2,   3, 6, 0,   3, 3,   2, 2,   0, 2, 2,   1, 1, 0,   0, 4, 1, 4,   2, 2, 0,   1, 1, 1, 0, 0, 0] # 0 is for skipped track
        nn = 1
        fitted_lines = [[] for i in 1:7] # all seven slots (first will be empty)
        for i in 1:length(selected_tracks)
            println("\nSession: $i")
            for j in 1:length(selected_tracks[i])
                push!(fitted_lines[i], [[], [], [], [], [], [], [], [], []])
                println("\tTrack: $j nn:$nn npsi $(npsi[nn])")
                x = convert(Array{Float64,1}, selected_tracks[i][j][1])
                y = convert(Array{Float64,1}, selected_tracks[i][j][2])
                push!(fitted_lines[i][end][1], x)
                x, y2, ysl, slopes, eslopes, bp, ebp, xl, exl = PyRModule.segmented(x, y, npsi[nn]; lambda=lambda, preview=false)
                nn += 1
                push!(fitted_lines[i][end][2], y2)
                push!(fitted_lines[i][end][3], ysl)  # for tracks
                push!(fitted_lines[i][end][4], slopes)
                push!(fitted_lines[i][end][5], eslopes)
                push!(fitted_lines[i][end][6], bp)
                push!(fitted_lines[i][end][7], ebp)
                push!(fitted_lines[i][end][8], xl) # xlin later on
                push!(fitted_lines[i][end][9], exl) # exlin later on
            end
        end


        #iis = [[2], [2], [3], [3], [4], [5], [6], [7], [7]]  # obs. session why iis?
        #jjs = [[2, 3], [4, 5],   [1, 2, 3], [4, 5],   [1, 2],   [1, 2, 3],   [1, 2, 3],    [1, 2],   [4, 5, 6]]  # tracks
        iis = [[2], [2], [2], [2], [4], [5], [6], [7], [7]]  # obs. session why iis?
        jjs = [[3, 4], [5, 6], [8, 9], [10,11], [2, 3], [1, 2], [2, 3, 4], [1, 2], [4, 5, 6]]  # tracks


        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
        subplots_adjust(left=0.08, bottom=0.09, right=0.99, top=0.92, wspace=0.0, hspace=0.0)
        figtext(0.09, 0.89, "a)", size=10)
        figtext(0.27, 0.89, "b)", size=10)
        figtext(0.53, 0.89, "c)", size=10)
        figtext(0.71, 0.89, "d)", size=10)
        figtext(0.865, 0.89, "e)", size=10)
        figtext(0.09, 0.425, "f)", size=10)
        figtext(0.24, 0.43, "g)", size=10)
        figtext(0.60, 0.43, "h)", size=10)
        figtext(0.82, 0.43, "i)", size=10)
        colnum = 800
        (r1,r11) = driftdirection_J1750_subplot3(1, 0, 0, 140, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r2, r22) = driftdirection_J1750_subplot3(2, 0, 159, 210, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r3, r33) = driftdirection_J1750_subplot3(3, 0, 388, 140, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r4, r44) = driftdirection_J1750_subplot3(4, 0, 547, 120, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r5, r55) = driftdirection_J1750_subplot3(5, 0, 686, 110, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r6, r66) = driftdirection_J1750_subplot3(6, 35, 0, 120, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r7, r77) = driftdirection_J1750_subplot3(7, 35, 133, 300, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r8, r88) = driftdirection_J1750_subplot3(8, 35, 447, 180, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r9, r99) = driftdirection_J1750_subplot3(9, 35, 641, 160, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        figtext(0.5, 0.97, "Pulse number", size=7, ha="center")
        figtext(0.5, 0.01, "Pulse number", size=7, ha="center")

        savefig("$outdir/$(name_mod)_driftdirection_3b.pdf")
        println("$outdir/$(name_mod)_driftdirection_3b.pdf")

        r = vcat(r1, r2, r3, r4, r5, r6, r7, r8, r9)
        rr = vcat(r11, r22, r33, r44, r55, r66, r77, r88, r99)
        println("Spline R^2: ", mean(r), " ", std(r))
        println("Linear reg. R^2: ", mean(rr), " ", std(rr))
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
    end


    function driftdirection_J1750_4(datas, outdir; lambda=1000.0, bin_st=nothing, bin_end=nothing, name_mod="0", show_=false)

        nums = []
        bins = []
        for data in datas
            nu, bi = size(data)
            push!(nums, nu)
            push!(bins, bi)
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins[1] end
        das = [] # single pulse data
        for data in datas
            da = data[:, bin_st:bin_end]
            push!(das, da)
        end

        tracks = []
        inclines = []
        dr_x = []
        dr_y = []
        for i in 1:7
            tracks1, lines1, inclines1, dr1, ysp1, dof1 = Tools.get_driftrate("$outdir/tracks/$i", lambda)
            push!(tracks, tracks1)
            push!(inclines, inclines1)
            push!(dr_x, dr1)
            push!(dr_y, ysp1)
        end

        #ranges = [[2, (0, 70)], [2, (500, 610)], [2, (630, 700)], [2, (960, 1020)], [4, (85, 140)], [5, (295, 355)], [6, (200, 350)], [7, (0, 90)], [7, (190, 270)]]
        #ranges = [[2, (0, 30)], [2, (290, 350)], [2, (500, 600)], [2, (630, 700)], [2, (960, 1020)], [4, (85, 140)], [5, (295, 355)], [6, (200, 350)], [7, (0, 90)]]
        ranges = [[2, (290, 350)], [2, (500, 600)], [2, (630, 700)], [2, (960, 1020)], [4, (85, 140)], [5, (295, 355)], [6, (200, 350)], [7, (0, 90)], [7, (190, 270)]]

        # rejected [2, (825, 922)] (npsi, 0, 2, 2)

        sum = 0
        i = 0
        for r in ranges
            i += 1
            pulses = r[2][2] - r[2][1]
            sum += pulses
            println("$i $pulses")
        end
        println("Sum: $sum")

        selected_tracks = [[] for i in 1:7] # all seven slots (first will be empty)
        #println(tracks[1][1][:, 2])  # pulse number

        for ra in ranges
            k = ra[1]
            v = ra[2]
            println("Session: ", k)
            for i in 1:length(tracks[k])
                create_new = true
                for (ii, pulse) in enumerate(tracks[k][i][:, 2])
                    if (pulse > v[1]) && (pulse < v[2])
                        if create_new
                            push!(selected_tracks[k], [[], []])
                            create_new = false
                            println("\tNew track session:$k track:$i pulse:$pulse loc:", tracks[k][i][ii, 1])
                        end
                        push!(selected_tracks[k][end][1], pulse) # x - pulse number
                        push!(selected_tracks[k][end][2], tracks[k][i][ii, 1]) # location
                        #println("\t $v $pulse ", tracks[k][i][ii, 1])
                    end
                end
            end
        end

        #npsi = [0, 0, 2, 2,   3, 6, 0,   3, 3,   2, 2,   0, 2, 2,   1, 1, 0,   0, 4, 1, 4,   2, 2, 0,   1, 1, 1, 0, 0, 0] # 0 is for skipped track
        #npsi = [1, 1,   1, 1,   3, 3,   3, 3,   2, 2,   0, 2, 2,   1, 1, 0,   0, 4, 4, 0,   2, 2, 0 ] # 0 is for skipped track
        npsi = [1, 1,   3, 3,   3, 3,   2, 2,   0, 2, 2,   1, 1, 0,   0, 4, 4, 0,   2, 2, 0,   1, 1] # 0 is for skipped track

        nn = 1
        fitted_lines = [[] for i in 1:7] # all seven slots (first will be empty)
        for i in 1:length(selected_tracks)
            println("\nSession: $i")
            for j in 1:length(selected_tracks[i])
                push!(fitted_lines[i], [[], [], [], [], [], [], [], [], []])
                println("\tTrack: $j nn:$nn npsi $(npsi[nn])")
                x = convert(Array{Float64,1}, selected_tracks[i][j][1])
                y = convert(Array{Float64,1}, selected_tracks[i][j][2])
                push!(fitted_lines[i][end][1], x)
                x, y2, ysl, slopes, eslopes, bp, ebp, xl, exl = PyRModule.segmented(x, y, npsi[nn]; lambda=lambda, preview=false)
                nn += 1
                push!(fitted_lines[i][end][2], y2)
                push!(fitted_lines[i][end][3], ysl)  # for tracks
                push!(fitted_lines[i][end][4], slopes)
                push!(fitted_lines[i][end][5], eslopes)
                push!(fitted_lines[i][end][6], bp)
                push!(fitted_lines[i][end][7], ebp)
                push!(fitted_lines[i][end][8], xl) # xlin later on
                push!(fitted_lines[i][end][9], exl) # exlin later on
            end
        end

        #iis = [[2], [2], [2], [2], [4], [5], [6], [7], [7]]  # obs. session why iis?
        #jjs = [[3, 4], [5, 6], [8, 9], [10,11], [2, 3], [1, 2], [2, 3, 4], [1, 2], [4, 5, 6]]  # tracks
        #iis = [[2], [2], [2], [2], [2], [4], [5], [6], [7]]  # obs. session why iis?
        #jjs = [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [2, 3], [1, 2], [2, 3], [1, 2]]  # tracks
        iis = [[2], [2], [2], [2], [4], [5], [6], [7], [7]]  # obs. session why iis?
        jjs = [[1, 2], [3, 4], [5, 6], [7, 8], [2, 3], [1, 2], [2, 3], [1, 2], [4, 5]]  # tracks


        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
        subplots_adjust(left=0.08, bottom=0.09, right=0.99, top=0.92, wspace=0.0, hspace=0.0)
        figtext(0.09, 0.89, "a)", size=10)
        figtext(0.25, 0.89, "b)", size=10)
        figtext(0.51, 0.89, "c)", size=10)
        figtext(0.70, 0.89, "d)", size=10)
        figtext(0.865, 0.89, "e)", size=10)
        figtext(0.09, 0.425, "f)", size=10)
        figtext(0.24, 0.43, "g)", size=10)
        figtext(0.60, 0.43, "h)", size=10)
        figtext(0.81, 0.43, "i)", size=10)
        colnum = 400
        (r1,r11) = driftdirection_J1750_subplot3(1, 0, 0, 60, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r2, r22) = driftdirection_J1750_subplot3(2, 0, 73, 100, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r3, r33) = driftdirection_J1750_subplot3(3, 0, 186, 70, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r4, r44) = driftdirection_J1750_subplot3(4, 0, 269, 60, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r5, r55) = driftdirection_J1750_subplot3(5, 0, 342, 55, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r6, r66) = driftdirection_J1750_subplot3(6, 35, 0, 60, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r7, r77) = driftdirection_J1750_subplot3(7, 35, 66, 150, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r8, r88) = driftdirection_J1750_subplot3(8, 35, 222, 90, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        (r9, r99) = driftdirection_J1750_subplot3(9, 35, 318, 80, selected_tracks, fitted_lines, iis, jjs, dr_x, dr_y; colnum=colnum)
        figtext(0.5, 0.97, "Pulse number", size=7, ha="center")
        figtext(0.5, 0.01, "Pulse number", size=7, ha="center")

        savefig("$outdir/$(name_mod)_driftdirection_4.pdf")
        println("$outdir/$(name_mod)_driftdirection_4.pdf")

        r = vcat(r1, r2, r3, r4, r5, r6, r7, r8, r9)
        rr = vcat(r11, r22, r33, r44, r55, r66, r77, r88, r99)
        println("Spline R^2: ", mean(r), " ", std(r))
        println("Linear reg. R^2: ", mean(rr), " ", std(rr))
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
    end



    function lrfs_average_J1750(slices, pulses, outdir; start=1, end_=nothing, step=10, number=128, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="1", verbose=false)
        num, bins = size(slices[1])
        if end_ == nothing end_ = num end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end

        LRFS = nothing
        for (i, sl) in enumerate(slices)
            da = sl[:, bin_st:bin_end]
            lrfs, intensity, freq, peak = Tools.lrfs(da)
            if LRFS == nothing
                LRFS = lrfs
            else
                LRFS += lrfs
            end
            #println(i)
            println(size(lrfs))
            println(size(intensity))
            println(size(freq))


        end


            imshow(abs.(lrfs), origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(abs.(lrfs)))

            savefig("$outdir/$(name_mod)_lrfs_average_J1750.pdf")
            close()


        return

        average = Tools.average_profile(da)
        lrfs, intensity, freq, peak = Tools.lrfs(da)
        phase_ = rad2deg.(angle.(view(lrfs, peak, :)))  # fft phase variation  # view used! lrfs[peak, :] -> copies data
        # skip freq = 0 and normalize intensity to 1
        inten = intensity[2:end]
        inten .-= minimum(inten)
        inten ./= maximum(inten)
        fre = freq[2:end]
        pars, errs = Tools.fit_gaussian(fre, inten; μ=freq[peak-1])  # skip zero freq
        fr = pars[2]
        frer = pars[3]  # errs[2] # yeap

        println("\tFrequency: $fr, P3: $(1/fr)")
        println("\tFrequency error: $frer, P3 error: $(1/fr - 1/(fr +frer))")  # TODO err ok?

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))


        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 4.33071))  # 8cm x 11cm
        subplots_adjust(left=0.13, bottom=0.08, right=0.90, top=0.92, wspace=0., hspace=0.)

        ax = subplot2grid((5, 3), (0, 1), colspan=2)
        minorticks_on()
        #plot(longitude, phase_, c="grey")
        scatter(longitude, phase_, marker=".", c="grey", s=3.)
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
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
        tick_params(labelleft=false, labelbottom=false)

        subplot2grid((5, 3), (4, 1), colspan=2)
        minorticks_on()
        plot(longitude, average, c="grey")
        axvline(x=0., ls=":", c="black")
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        xlabel("longitude \$(^\\circ)\$")
        #tick_params(labeltop=false, labelbottom=true)
        savefig("$outdir/$(name_mod)_lrfs_average_J1750.pdf")
        close()

    end


    function singlepulses_J1750(datas, outdir; name_mod="123456", darkness=1.0, cmap="viridis", show_=false)

        pulses = 0
        extents = [] # x - npulses, y - bins
        intens = []
        for data in datas
            pulse_num = size(data)[1]
            pulses += pulse_num
            println(pulse_num)

            intensity, pulsess = Tools.intensity_pulses(data)
            intensity .-= minimum(intensity)
            intensity ./= maximum(intensity)
            push!(extents, [1, pulse_num, 0.0, 360.0])
            push!(intens, [pulsess, intensity])
        end
        println(pulses)
        println(pulses/2)
        a = size(datas[1])[1] + size(datas[2])[1] + size(datas[3])[1]
        b = size(datas[4])[1] + size(datas[5])[1] + size(datas[6])[1] + size(datas[7])[1]

        println("$a $((1900-a)/2) ")
        println("$b $((1900-b)/3)")

        #return

        # bin range
        y1 = 600 /1024 *360.0
        y2 = 360 / 1024 * 360.0 # to avoid fix_ticsks!
        #cmap = "binary"

        #=
        # HERE
        fl = Base.Iterators.flatten(datas[2])
        PyPlot.hist(collect(fl), bins=100)
        show()
        #println(size(vcat(transpose(datas[1]))))
        readline(stdin; keep=false)
        #return
        =#

        rc("font", size=7.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        #figure(figsize=(7.086614, 6.299213))  # 18 cm x 16 cm
        figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
        subplots_adjust(left=0.08, bottom=0.07, right=0.95, top=0.92, wspace=0.0, hspace=0.0)
        figtext(0.09, 0.89, "a)", size=10, c="white")
        figtext(0.25, 0.89, "b)", size=10, c="white")
        figtext(0.76, 0.89, "c)", size=10, c="white")
        figtext(0.09, 0.44, "d)", size=10, c="white")
        figtext(0.32, 0.44, "e)", size=10, c="white")
        figtext(0.54, 0.44, "f)", size=10, c="white")
        figtext(0.76, 0.44, "g)", size=10, c="white")

        # first session
        ax = subplot2grid((32, 1900), (0, 0), colspan=295, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        minorticks_on()
        #imshow(transpose(datas[1]), origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(datas[1]), extent=extents[1])
        imshow(transpose(datas[1]), origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(datas[1]), vmin=0, extent=extents[1], alpha=1.0)
        println("max ", maximum(datas[1]))
        println("min ", minimum(datas[1]))
        ylim(y1, y2)
        ylabel("Longitude (\$ ^{\\circ}\$)")
        # first session
        subplot2grid((32, 1900), (10, 0), colspan=295, rowspan=5)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        minorticks_on()
        plot(intens[1][1], intens[1][2], c="grey")
        xlim(intens[1][1][1]-0.5, intens[1][1][end]+0.5)
        ylim(-0.2, 1.2)
        #xticks([0.5, 1.0])
        ylabel("Intensity (a. u.)")

        # second session
        ax = subplot2grid((32, 1900), (0, 359), colspan=1031, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        #imshow(transpose(datas[2]), origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(datas[2]), extent=extents[2])
        imshow(transpose(datas[2]), origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(datas[2]), vmin=0.0, extent=extents[2], alpha=1.0)
        ylim(y1, y2)
        xlabel("Pulse number")
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        # second session
        subplot2grid((32, 1900), (10, 359), colspan=1031, rowspan=5)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true)
        minorticks_on()
        plot(intens[2][1], intens[2][2], c="grey")
        xlim(intens[2][1][1]-0.5, intens[2][1][end]+0.5)
        ylim(-0.2, 1.2)

        # third session
        ax = subplot2grid((32, 1900), (0, 1454), colspan=446, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true, labelright=true)
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        #ax.yaxis.set_label_position("right")
        #ax.yaxis.set_ticks_position("right")
        minorticks_on()
        #imshow(transpose(datas[3]), origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(datas[3]), extent=extents[3])
        imshow(transpose(datas[3]), origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(datas[3]), vmin=0.0, extent=extents[3], alpha=1.0)
        ylim(y1, y2)
        # third session
        ax = subplot2grid((32, 1900), (10, 1454), colspan=446, rowspan=5)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true, labelright=true)
        minorticks_on()
        #ax.yaxis.set_label_position("right")
        #ax.yaxis.set_ticks_position("right")
        plot(intens[3][1], intens[3][2], c="grey")
        xlim(intens[3][1][1]-0.5, intens[3][1][end]+0.5)
        ylim(-0.2, 1.2)

        # fourth session
        ax = subplot2grid((32, 1900), (17, 0), colspan=450, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        minorticks_on()
        #imshow(transpose(datas[4]), origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(datas[4]), extent=extents[4])
        imshow(transpose(datas[4]), origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(datas[4]), vmin=0.0, extent=extents[4], alpha=1.0)
        ylim(y1, y2)
        ylabel("Longitude (\$ ^{\\circ}\$)")
        # fourth session
        subplot2grid((32, 1900), (27, 0), colspan=450, rowspan=5)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        minorticks_on()
        plot(intens[4][1], intens[4][2], c="grey")
        xlim(intens[4][1][1]-0.5, intens[4][1][end]+0.5)
        locator_params(nbins=2)
        ylim(-0.2, 1.2)
        ylabel("Intensity (a. u.)")

        # fifth session
        ax = subplot2grid((32, 1900), (17, 487), colspan=447, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        #imshow(transpose(datas[5]), origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(datas[5]), extent=extents[5])
        imshow(transpose(datas[5]), origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(datas[5]), vmin=0.0, extent=extents[5], alpha=1.0)
        ylim(y1, y2)
        # fifth session
        subplot2grid((32, 1900), (27, 487), colspan=447, rowspan=5)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        plot(intens[5][1], intens[5][2], c="grey")
        xlim(intens[5][1][1]-0.5, intens[5][1][end]+0.5)
        ylim(-0.2, 1.2)
        #xlabel("Pulse number")

        # sixth session
        ax = subplot2grid((32, 1900), (17, 971), colspan=447, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        #imshow(transpose(datas[6]), origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(datas[6]), extent=extents[6])
        imshow(transpose(datas[6]), origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(datas[6]), vmin=0.0, extent=extents[6], alpha=1.0)
        ylim(y1, y2)
        # sixth session
        subplot2grid((32, 1900), (27, 971), colspan=447, rowspan=5)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true)
        minorticks_on()
        plot(intens[6][1], intens[6][2], c="grey")
        xlim(intens[6][1][1]-0.5, intens[6][1][end]+0.5)
        ylim(-0.2, 1.2)
        #xlabel("Pulse number")
        figtext(0.5, 0.01, "Pulse number", size=7, ha="center")


        # seventh session
        ax = subplot2grid((32, 1900), (17, 1454), colspan=446, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, labelbottom=false, top=true)
        #tick_params(axis="x", which="both", direction="in", labeltop=false, labelbottom=true, top=false)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true, labelright=true)
        minorticks_on()
        #imshow(transpose(datas[7]), origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(datas[7]), extent=extents[7])
        imshow(transpose(datas[7]), origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(datas[7]), vmin=0.0, extent=extents[7], alpha=1.0)
        ylim(y1, y2)
        # seventh session
        subplot2grid((32, 1900), (27, 1454), colspan=446, rowspan=5)  # row column
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true, left=true, labelright=true)
        minorticks_on()
        plot(intens[7][1], intens[7][2], c="grey")
        xlim(intens[7][1][1]-0.5, intens[7][1][end]+0.5)
        ylim(-0.2, 1.2)

        println("$outdir/$(name_mod)_singlepulses.pdf")
        savefig("$outdir/$(name_mod)_singlepulses.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
        #clf()
    end


    function modeledpulses_J1750(datas, outdir; name_mod="123456", darkness=[1.0, 1.0], cmap="viridis", show_=false)

        pulses = 0
        extents = [] # x - npulses, y - bins
        intens = []
        for data in datas
            pulse_num = size(data)[1]
            pulses += pulse_num
            println(pulse_num)
            intensity, pulsess = Tools.intensity_pulses(data)
            intensity .-= minimum(intensity)
            intensity ./= maximum(intensity)
            push!(extents, [1, pulse_num, 0.0, 360.0])
            push!(intens, [pulsess, intensity])
        end

        # read p3
        f = open("/home/szary/work/MeerTime/J1750/1234567_p3.txt")
        pul = []
        p3s = []
        p3obs = []
        p3max = 150
        for (i, line) in enumerate(readlines(f))
            res = split(line)
            push!(pul, parse(Int, res[1]))
            push!(p3s, parse(Float64, res[2]))
            #if (p3s[end] > 0.99) &&  (p3s[end] < 1.01)
            push!(p3obs, Functions.p3_obs(p3s[end], 1))
            if p3obs[end] > p3max
                p3obs[end] = p3max
            end
            if p3obs[end] < -p3max
                p3obs[end] = -p3max
            end
        end
        close(f)

        mi = minimum(p3s)
        ma = maximum(p3s)
        println("P3 min: ", mi, " max ", ma)


        rc("font", size=7.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        #figure(figsize=(7.086614, 6.299213))  # 18 cm x 16 cm
        figure(figsize=(6.299213, 3.8931373288047713))  # 16 cm x # golden ratio
        subplots_adjust(left=0.09, bottom=0.1, right=0.91, top=0.92, wspace=0.0, hspace=0.0)

        # first session
        ax = subplot2grid((3, 1), (0, 0))
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=false, left=true)
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        minorticks_on()
        imshow(transpose(datas[1]), origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness[1]*maximum(datas[1]), vmin=-0.03,  extent=extents[1])
        ylim(210, 131)
        ylabel("Longitude (\$ ^{\\circ}\$)")

        ax = subplot2grid((3, 1), (1, 0))
        imshow(transpose(datas[2]), origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness[2]*maximum(datas[2]), extent=extents[1])
        #ylim(165, 85)
        ylim(210-15, 131-15)
        ylabel("Longitude (\$ ^{\\circ}\$)")

        ax = subplot2grid((3, 1), (2, 0))
        minorticks_on()
        col1 = "C1" # "tab:blue"
        col2 = "C2" # "tab:red"
        plot(pul, p3s, lw=1.5, color=col1)
        ax.tick_params(axis="y", labelcolor=col1)
        xlim(pul[1], pul[end])
        xlabel("pulse number")
        ylabel("\$P_3\$", color=col1)
        ax2 = ax.twinx()
        minorticks_on()
        plot(pul, p3obs, color=col2, alpha=0.7, ls="--")
        ax2.tick_params(axis="y", labelcolor=col2)
        ylabel("\$P_3^{\\rm obs}\$", color=col2)


        println("$outdir/$(name_mod)_modeledpulses.pdf")
        savefig("$outdir/$(name_mod)_modeledpulses.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
        #clf()
    end


    function p3evolutions_J1750(datas, outdir, step, number; bin_st=nothing, bin_end=nothing, name_mod="1234567", darkness=1.0, cmap="viridis", show_=false)

        none, bins = size(datas[1])
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end

        pulses = 0

        extents = [] # x - npulses, y - bins
        intensities = []
        p3s = []
        p3errs = []
        freqs = []
        inten_freqs = []

        for data in datas
        #for i in 1:3 #datas
        #    data = datas[i]
            (intens, p3, p3err, freq) = Tools.p3evolution(data, step, number, bin_st, bin_end)
            push!(intensities, intens)
            push!(p3s, p3)
            push!(p3errs, p3err)
            push!(freqs, freq)
            inten_freq, skip = Tools.intensity_pulses(transpose(intens))
            push!(inten_freqs, inten_freq)
            println(size(intens))
            pulses += size(intens)[1]
        end
        #println(minimum(intensities[1]))
        #println(maximum(intensities[1]))

        #bottom, skip = Tools.intensity_pulses(transpose(intens))

        # bin range why? remove
        #y1 = 360 / 1024 * 360.0
        #y2 = 650 /1024 *360.0

        stable_indexes = [[1, 167]], [[660, 773], [852, 903]], [#= 3 =#], [ #= 4 =#],  [[1, 88]], [[8, 73]], [[20, 103], [228, 318]]

        change_indx = [[#=1=#],  [#=2=#[131, 165], [463, 534]], [#=3=#], [#=4=#[64, 158]], [#=5=#], [#=6=#[74, 115], [135, 156]], [#=7=#]]

        dps = []
        println("Stable P3:")
        for (i, rngs) in enumerate(stable_indexes)
            for rng in rngs
                p3 = mean(p3s[i][rng[1]:rng[2]])
                err = std(p3s[i][rng[1]:rng[2]])
                push!(dps, rng[2]-rng[1] )
                @printf "%d P3: %.1f Err: %.1f %d %d dpulse = %d\n" i p3 err rng[1] rng[2] dps[end]
            end
        end
        println("\t Mean length: $(mean(dps)) std:$(std(dps))")

        dps = []
        println("Increasing/decreesing P3:")
        for (i, rngs) in enumerate(change_indx)
            for rng in rngs
                p3_1 = p3s[i][rng[1]]
                p3_2 = p3s[i][rng[2]]
                ep3_1 = p3errs[i][rng[1]]
                ep3_2 = p3errs[i][rng[2]]
                push!(dps, rng[2]-rng[1] )
                @printf "%d P3[1]: %.1f +/- %.1f f P3[2]: %.1f +/- %.1f %d %d dpulse=%d\n" i p3_1 ep3_1 p3_2 ep3_2 rng[1] rng[2] dps[end]
            end
        end
        println("\t Mean length: $(mean(dps)) std:$(std(dps))")
        # CHECKING periodicity here
        #=
        for set in 1:size(p3s, 1)
            p3s_ = convert(Array{Float64,1}, p3s[set])
            Tools.periodicity(p3s_)
        end
        =#

        p3_min = 20
        p3_max = 95

        rc("font", size=7.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        #figure(figsize=(7.086614, 6.299213))  # 18 cm x 16 cm
        figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
        subplots_adjust(left=0.08, bottom=0.07, right=0.95, top=0.92, wspace=0.0, hspace=0.0)
        figtext(0.083, 0.89, "a)", size=10)
        figtext(0.22, 0.89, "b)", size=10)
        figtext(0.75, 0.89, "c)", size=10)
        figtext(0.083, 0.44, "d)", size=10)
        figtext(0.31, 0.44, "e)", size=10)
        figtext(0.535, 0.44, "f)", size=10)
        figtext(0.755, 0.44, "g)", size=10)
        figtext(0.5, 0.01, "start period num.", size=7, ha="center")

        # first session
        subplot2grid((32, 1600), (0, 0), colspan=60, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, top=false, bottom=false, labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        minorticks_on()
        plot(inten_freqs[1], freqs[1], c="grey")
        ylabel("frequency (1/P)")
        xl = xlim()
        xlim(xl[2], xl[1])
        ylim(freqs[1][1]-(freqs[1][2]-freqs[1][1])/2.0, freqs[1][end]+(freqs[1][end]- freqs[1][end-1])/2.0) # WOW
        # first session
        ax = subplot2grid((32, 1600), (0, 60), colspan=167, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=true, left=true)
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        minorticks_on()
        imshow(transpose(intensities[1]), origin="lower", cmap=cmap, interpolation="none", aspect="auto")
        #imshow(transpose(intensities[1]), origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness*maximum(intensities[1]), vmin=darkness*minimum(intensities[1]))
        #imshow(intensities[1], origin="lower", cmap=cmap, interpolation="none", aspect="auto", , vmax=darkness*maximum(intensities[1]))
        #ylim(y1, y2)
        # first session
        subplot2grid((32, 1600), (10, 60), colspan=167, rowspan=5)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true, left=true)
        minorticks_on()
        errorbar(collect(1:length(p3s[1])), p3s[1], yerr=p3errs[1], color="none", lw=0.1, marker="_", mec="grey", ecolor="grey", capsize=0, mfc="grey", ms=0.7)
        for rng in stable_indexes[1]
            errorbar(collect(rng[1]:rng[2]), p3s[1][rng[1]:rng[2]], yerr=p3errs[1][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:blue", ecolor="tab:blue", capsize=0, mfc="tab:blue", ms=0.7)
        end
        for rng in change_indx[1]
            errorbar(collect(rng[1]:rng[2]), p3s[1][rng[1]:rng[2]], yerr=p3errs[1][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:red", ecolor="tab:red", capsize=0, mfc="tab:red", ms=0.7)
        end
        xlim(1, length(p3s[1]))
        ylim(p3_min, p3_max)
        #xticks([0.5, 1.0])
        ylabel("\$P_{3}\$")

        # second session
        subplot2grid((32, 1600), (0, 243), colspan=60, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, top=false, bottom=false, labelbottom=false)
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=false, left=true)
        minorticks_on()
        plot(inten_freqs[2], freqs[2], c="grey")
        xl = xlim()
        xlim(xl[2], xl[1])
        ylim(freqs[2][1]-(freqs[2][2]-freqs[2][1])/2.0, freqs[2][end]+(freqs[2][end]- freqs[2][end-1])/2.0) # WOW
        # second session
        ax = subplot2grid((32, 1600), (0, 303), colspan=903, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=true, left=true)
        minorticks_on()
        imshow(transpose(intensities[2]), origin="lower", cmap=cmap, interpolation="none", aspect="auto")
        xlabel("start period num.")
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        # second session
        subplot2grid((32, 1600), (10, 303), colspan=903, rowspan=5)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true)
        minorticks_on()
        errorbar(collect(1:length(p3s[2])), p3s[2], yerr=p3errs[2], color="none", lw=0.1, marker="_", mec="grey", ecolor="grey", capsize=0, mfc="grey", ms=0.7)
        for rng in stable_indexes[2]
            errorbar(collect(rng[1]:rng[2]), p3s[2][rng[1]:rng[2]], yerr=p3errs[2][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:blue", ecolor="tab:blue", capsize=0, mfc="tab:blue", ms=0.7)
        end
        for rng in change_indx[2]
            errorbar(collect(rng[1]:rng[2]), p3s[2][rng[1]:rng[2]], yerr=p3errs[2][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:red", ecolor="tab:red", capsize=0, mfc="tab:red", ms=0.7)
        end
        xlim(1, length(p3s[2]))
        ylim(p3_min, p3_max)

        # third session
        subplot2grid((32, 1600), (0, 1222), colspan=60, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, top=false, bottom=false, labelbottom=false)
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=false, left=true)
        minorticks_on()
        plot(inten_freqs[3], freqs[3], c="grey")
        xl = xlim()
        xlim(xl[2], xl[1])
        ylim(freqs[3][1]-(freqs[3][2]-freqs[3][1])/2.0, freqs[3][end]+(freqs[3][end]- freqs[3][end-1])/2.0) # WOW
        # third session
        ax = subplot2grid((32, 1600), (0, 1282), colspan=318, rowspan=10)  # row column
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=true, left=true)
        minorticks_on()
        imshow(transpose(intensities[3]), origin="lower", cmap=cmap, interpolation="none", aspect="auto")
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        # third session
        subplot2grid((32, 1600), (10, 1282), colspan=318, rowspan=5)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true)
        minorticks_on()
        errorbar(collect(1:length(p3s[3])), p3s[3], yerr=p3errs[3], color="none", lw=0.1, marker="_", mec="grey", ecolor="grey", capsize=0, mfc="grey", ms=0.7)
        for rng in stable_indexes[3]
            errorbar(collect(rng[1]:rng[2]), p3s[3][rng[1]:rng[2]], yerr=p3errs[3][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:blue", ecolor="tab:blue", capsize=0, mfc="tab:blue", ms=0.7)
        end
        for rng in change_indx[3]
            errorbar(collect(rng[1]:rng[2]), p3s[3][rng[1]:rng[2]], yerr=p3errs[3][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:red", ecolor="tab:red", capsize=0, mfc="tab:red", ms=0.7)
        end
        xlim(1, length(p3s[3]))
        ylim(p3_min, p3_max)

        # fourth session
        subplot2grid((32, 1600), (17, 0), colspan=60, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, top=false, bottom=false, labelbottom=false)
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=false, left=true)
        minorticks_on()
        plot(inten_freqs[4], freqs[4], c="grey")
        xl = xlim()
        xlim(xl[2], xl[1])
        ylim(freqs[4][1]-(freqs[4][2]-freqs[4][1])/2.0, freqs[4][end]+(freqs[4][end]- freqs[4][end-1])/2.0) # WOW
        ylabel("frequency (1/P)")
        # fourth session
        ax = subplot2grid((32, 1600), (17, 60), colspan=322, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="in", top=true, labeltop=false, bottom=false, labelbottom=false)
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=true, left=true)
        minorticks_on()
        imshow(transpose(intensities[4]), origin="lower", cmap=cmap, interpolation="none", aspect="auto")
        # fourth session
        subplot2grid((32, 1600), (27, 60), colspan=322, rowspan=5)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=true)
        tick_params(axis="y", which="both", direction="out", labelleft=true, right=true)
        minorticks_on()
        errorbar(collect(1:length(p3s[4])), p3s[4], yerr=p3errs[4], color="none", lw=0.1, marker="_", mec="grey", ecolor="grey", capsize=0, mfc="grey", ms=0.7)
        for rng in stable_indexes[4]
            errorbar(collect(rng[1]:rng[2]), p3s[4][rng[1]:rng[2]], yerr=p3errs[4][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:blue", ecolor="tab:blue", capsize=0, mfc="tab:blue", ms=0.7)
        end
        for rng in change_indx[4]
            errorbar(collect(rng[1]:rng[2]), p3s[4][rng[1]:rng[2]], yerr=p3errs[4][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:red", ecolor="tab:red", capsize=0, mfc="tab:red", ms=0.7)
        end
        xlim(1, length(p3s[4]))
        ylim(p3_min, p3_max)
        ylabel("\$P_{3}\$")

        # fifth session
        subplot2grid((32, 1600), (17, 409), colspan=60, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, top=false, bottom=false, labelbottom=false)
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=false, left=true)
        minorticks_on()
        plot(inten_freqs[5], freqs[5], c="grey")
        xl = xlim()
        xlim(xl[2], xl[1])
        ylim(freqs[5][1]-(freqs[5][2]-freqs[5][1])/2.0, freqs[5][end]+(freqs[5][end]- freqs[5][end-1])/2.0) # WOW
        # fifth session
        ax = subplot2grid((32, 1600), (17, 469), colspan=319, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="in", top=true, labeltop=false, bottom=false, labelbottom=false)
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=true, left=true)
        minorticks_on()
        imshow(transpose(intensities[5]), origin="lower", cmap=cmap, interpolation="none", aspect="auto")
        # fifth session
        subplot2grid((32, 1600), (27, 469), colspan=319, rowspan=5)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=true)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true)
        minorticks_on()
        errorbar(collect(1:length(p3s[5])), p3s[5], yerr=p3errs[5], color="none", lw=0.1, marker="_", mec="grey", ecolor="grey", capsize=0, mfc="grey", ms=0.7)
        for rng in stable_indexes[5]
            errorbar(collect(rng[1]:rng[2]), p3s[5][rng[1]:rng[2]], yerr=p3errs[5][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:blue", ecolor="tab:blue", capsize=0, mfc="tab:blue", ms=0.7)
        end
        for rng in change_indx[5]
            errorbar(collect(rng[1]:rng[2]), p3s[5][rng[1]:rng[2]], yerr=p3errs[5][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:red", ecolor="tab:red", capsize=0, mfc="tab:red", ms=0.7)
        end
        xlim(1, length(p3s[5]))
        ylim(p3_min, p3_max)

        # sixth session
        subplot2grid((32, 1600), (17, 815), colspan=60, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="out", labeltop=false, top=false, bottom=false, labelbottom=false)
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=false, left=true)
        minorticks_on()
        plot(inten_freqs[6], freqs[6], c="grey")
        xl = xlim()
        xlim(xl[2], xl[1])
        ylim(freqs[6][1]-(freqs[6][2]-freqs[6][1])/2.0, freqs[6][end]+(freqs[6][end]- freqs[6][end-1])/2.0) # WOW
        # sixth session
        ax = subplot2grid((32, 1600), (17, 875), colspan=319, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="in", top=true, labeltop=false, bottom=false, labelbottom=false)
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=true, left=true)
        minorticks_on()
        imshow(transpose(intensities[6]), origin="lower", cmap=cmap, interpolation="none", aspect="auto")
        # sixth session
        subplot2grid((32, 1600), (27, 875), colspan=319, rowspan=5)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=true)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true)
        minorticks_on()
        errorbar(collect(1:length(p3s[6])), p3s[6], yerr=p3errs[6], color="none", lw=0.1, marker="_", mec="grey", ecolor="grey", capsize=0, mfc="grey", ms=0.7)
        for rng in stable_indexes[6]
            errorbar(collect(rng[1]:rng[2]), p3s[6][rng[1]:rng[2]], yerr=p3errs[6][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:blue", ecolor="tab:blue", capsize=0, mfc="tab:blue", ms=0.7)
        end
        for rng in change_indx[6]
            errorbar(collect(rng[1]:rng[2]), p3s[6][rng[1]:rng[2]], yerr=p3errs[6][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:red", ecolor="tab:red", capsize=0, mfc="tab:red", ms=0.7)
        end
        xlim(1, length(p3s[6]))
        ylim(p3_min, p3_max)

        # seventh session
        subplot2grid((32, 1600), (17, 1221), colspan=60, rowspan=10)  # row column
        minorticks_on()
        plot(inten_freqs[7], freqs[7], c="grey")
        xl = xlim()
        xlim(xl[2], xl[1])
        ylim(freqs[7][1]-(freqs[7][2]-freqs[7][1])/2.0, freqs[7][end]+(freqs[7][end]- freqs[7][end-1])/2.0) # WOW
        tick_params(axis="x", which="both", direction="out", labeltop=false, top=false, bottom=false, labelbottom=false)
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=false, left=true)
        # seventh session
        ax = subplot2grid((32, 1600), (17, 1281), colspan=318, rowspan=10)  # row column
        tick_params(axis="x", which="both", direction="in", top=true, labeltop=false, bottom=false, labelbottom=false)
        tick_params(axis="y", which="both", direction="in", labelleft=false, right=true, left=true)
        minorticks_on()
        imshow(transpose(intensities[7]), origin="lower", cmap=cmap, interpolation="none", aspect="auto")
        # seventh session
        subplot2grid((32, 1600), (27, 1281), colspan=318, rowspan=5)  # row column
        tick_params(axis="x", which="both", direction="out", labelbottom=true)
        tick_params(axis="y", which="both", direction="out", labelleft=false, right=true)
        minorticks_on()
        errorbar(collect(1:length(p3s[7])), p3s[7], yerr=p3errs[7], color="none", lw=0.1, marker="_", mec="grey", ecolor="grey", capsize=0, mfc="grey", ms=0.7)
        for rng in stable_indexes[7]
            errorbar(collect(rng[1]:rng[2]), p3s[7][rng[1]:rng[2]], yerr=p3errs[7][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:blue", ecolor="tab:blue", capsize=0, mfc="tab:blue", ms=0.7)
        end
        for rng in change_indx[7]
            errorbar(collect(rng[1]:rng[2]), p3s[7][rng[1]:rng[2]], yerr=p3errs[7][rng[1]:rng[2]], color="none", lw=0.1, marker="_", mec="tab:red", ecolor="tab:red", capsize=0, mfc="tab:red", ms=0.7)
        end
        xlim(1, length(p3s[7]))
        ylim(p3_min, p3_max)

        println("$outdir/$(name_mod)_p3evolutions.pdf")
        savefig("$outdir/$(name_mod)_p3evolutions.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()
        #clf()
    end


    function test_track_subpulses_snr(data, outdir, p2, snrfile; threshs=(7, 5), thresh2=0.6, bin_st=350, bin_end=650, name_mod="X", show_=true)

        thresh2_old = thresh2
        thresh2 = []
        pulse_num, bins = size(data)


        f = open(snrfile)
        snrs = []
        for line in readlines(f)
            push!(snrs, parse(Float64, line))
        end

        # pick random pulses with thresh tolerance tol
        tol = 0.2
        pulses = []
        for th in threshs
            pulse = nothing
            while pulse == nothing
                    pu = rand(1:pulse_num)
                    #println(snrs[pu])
                    if (snrs[pu] > th - tol) && (snrs[pu] < th + tol)
                        pulse = pu
                    end
            end
            push!(pulses, pulse)
        end

        #pulses = [544, 725]
        pulses = [38, 725]

        # detect peaks (copied from Tools.track_subpulses_snr)
        p2_bins = floor(Int, p2 / 360 * bins)
        peaks = [] # [pulse_num, [p1, p2, p3...]]
        ppeaks = [] # [p1, p2, p3...]
        σ = p2_bins / 2 / 2.35482
        kernel = Tools.gauss(collect(1:p2_bins), [1, p2_bins/2, σ, 0])
        ress = []
        ress_acorr = []
        for (kk, i) in enumerate(pulses)
            #y = view(data, i, on_st:on_end)
            y = view(data, i, :)
            (mi, ma) = extrema(y)
            #y = (y .- mi) / (ma - mi)
            y = y ./ ma # this is much much better!

            # convolution with Gaussian
            res = Tools.conv(y, kernel)
            (mi, ma) = extrema(res)
            res = (res .- mi) / (ma - mi)
            re = res[floor(Int,p2_bins/2)+1:end-floor(Int,p2_bins/2)]
            #println(floor(Int,p2_bins/2), " ", floor(Int,p2_bins/2))
            #println(length(y))
            #println(length(re))

            #println(length(y))
            #println(length(kernel))
            #println(length(res))

            # auto correlation # just for some tests..
            resa = crosscor(y, y, collect(floor(Int,-length(y)/2):floor(Int,length(y)/2)))
            (mi, ma) = extrema(resa)
            resa = (resa .- mi) / (ma - mi)
            rea = resa[floor(Int,p2_bins/2):end-floor(Int,p2_bins/2)]

            peak = Tools.peaks(re)
            # new syntax?
            ma, pa = peakprom(Maxima(), re, floor(Int, p2_bins/4))
            #ma, pa = peakprom(re, Maxima(), floor(Int, p2_bins/4))
            #ma, pa = peakprom(re, Maxima(), p2_bins/2)
            inds = sortperm(pa, rev=true)
            # get new thresh2 (maximum in off pulse region)
            new_thresh = 0.
            for ii in inds
                if ((ma[ii] < bin_st) || (ma[ii] > bin_end)) && (new_thresh < re[ma[ii]])
                    new_thresh = re[ma[ii]]
                end
            end
            push!(thresh2, maximum([new_thresh, thresh2_old]))
            ppeaks = []
            for ii in inds
                if (re[ma[ii]] >= thresh2[kk]) &&  (ma[ii] > bin_st) && (ma[ii] < bin_end)
                    push!(ppeaks, ma[ii])
                end
            end

            #println(length(ppeaks))
            comps = length(ppeaks)
            if (comps == 1)
                cen = ppeaks[end]
                x_ = collect(cen-p2_bins:cen+p2_bins)
                pars, errs = Tools.fit_gaussian(x_, y[cen-p2_bins:cen+p2_bins]; μ=ppeaks[end])  # skip zero freq
                #println(pars)
                #println(errs)
                println(cen)
                println(round(pars[2]), " ± ", round(errs[2]), "(bins) ± ", round(errs[2]/bins*360.0, digits=1), "(deg)")
                g = Tools.gauss(x_, pars)
            elseif (comps == 2)
                m1 = minimum(ppeaks)
                m2 = maximum(ppeaks)
                x_ = collect(m1-p2_bins:m2+p2_bins)
                y_ = y[m1-p2_bins:m2+p2_bins]
                pars, errs = Tools.fit_twogaussians(x_, y_, re[m1], re[m2], m1, m2,p2_bins/2, p2_bins/2)
                #println(pars)
                println(m1)
                println(round(pars[2]), " ± ", round(errs[2]), "(bins) ± ", round(errs[2]/bins*360.0, digits=1), "(deg)")
                println(m2)
                println(round(pars[5]), " ± ", round(errs[5]), "(bins) ± ", round(errs[5]/bins*360.0, digits=1), "(deg)")
                g = Tools.twogauss(x_, pars)
            else
                println("Not implemented...")
                x_ = []
                g = []
            end

            axhline(y=thresh2[end], c="black", lw=2)
            plot(y)
            plot(re)
            plot(x_, g)
            show()
            xlim((bin_st, bin_end))
            readline(stdin; keep=false)
            PyPlot.close()

            push!(ress, re)
            push!(ress_acorr, rea)
            push!(peaks, [i, ppeaks])
        end

        #return

        # Pulse longitude
        longitude = collect(range(bin_st/bins * 360, bin_end/bins*360, length=bin_end+1-bin_st))

        println(pulses)

        rc("font", size=7.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 2.362205), frameon=true)  # 8cm x 6 cm
        #figure(figsize=(3.14961, 1.946563744), frameon=true)  # 8cm x 4.94427191 cm (golden)
        subplots_adjust(left=0.17, bottom=0.15, right=0.99, top=0.99, wspace=0.0, hspace=0.0)

        subplot2grid((2, 1), (0, 0))  # row column
        tick_params(labelbottom=true)
        minorticks_on()
        y = view(data, pulses[1], bin_st:bin_end)
        (mi, ma) = extrema(y)
        y = y ./ ma # this is much much better!
        plot(longitude, y, c="black", lw=0.5)
        plot(longitude, ress[1][bin_st:bin_end], c="tab:red")
        #plot(longitude, ress_acorr[1][bin_st:bin_end], c="tab:blue") # for tests only
        axhline(y=thresh2[1], ls=":", lw=0.7, c="tab:green")
        for p in peaks[1][2]
            axvline(x=p / bins * 360, c="tab:blue", ls="--")
        end
        xlim(longitude[1], longitude[end])
        yticks([-0.5, 0.0, 0.5])
        text(205, 0.65, "S/N \$\\approx\$ 3", size=7)
        ylabel("intensity (a. u.)")

        subplot2grid((2, 1), (1, 0))  # row column
        minorticks_on()
        y = view(data, pulses[2], bin_st:bin_end)
        (mi, ma) = extrema(y)
        y = y ./ ma # this is much much better!
        plot(longitude, y, c="black", lw=0.5)
        plot(longitude, ress[2][bin_st:bin_end], c="tab:red")
        #plot(longitude, ress_acorr[2][bin_st:bin_end], c="tab:blue") # for tests only
        axhline(y=thresh2[2], ls=":", lw=0.7, c="tab:green")
        for p in peaks[2][2]
            axvline(x=p / bins * 360, c="tab:blue", ls="--")
        end
        xlim(longitude[1], longitude[end])
        yticks([-0.5, 0.0, 0.5])
        text(205, 0.77, "S/N \$\\approx\$ 5", size=7)
        ylabel("intensity (a. u.)")
        xlabel("longitude \$(^\\circ)\$")

        filename = "$outdir/$(name_mod)_track_subpulses.pdf"
        println(filename)
        savefig(filename)
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        PyPlot.close()

    end


    function test_track_subpulses_snr_new(data, outdir, p2, snrfile; pulse_snr=7, thresh=5, bin_st=350, bin_end=650, off_st=20, off_end=320, name_mod="X", show_=true)
        """ New approch based on S/N of the pulse"""
        pulse_num, bins = size(data)

        f = open(snrfile)
        snrs = []
        for line in readlines(f)
            push!(snrs, parse(Float64, line))
        end

        # pick random pulse with thresh tolerance tol
        tol = 0.1
        pulse = nothing
        while pulse == nothing
                pu = rand(1:pulse_num)
                #println(snrs[pu])
                if (snrs[pu] > pulse_snr - tol) && (snrs[pu] < pulse_snr + tol)
                    pulse = pu
                end
        end
        #pulse = rand(1:pulse_num)
        pulse = 1012 # looks good
        println("pulse S/N: ", snrs[pulse])
        # olds
        #pulses = [544, 725]
        #pulses = [38, 725]

        # detect peaks (copied from Tools.track_subpulses_snr)
        p2_bins = floor(Int, p2 / 360 * bins)
        peaks = [] # [p1, p2, p3...]
        σ = p2_bins / 2 / 2.35482
        kernel = Tools.gauss(collect(1:p2_bins), [1, p2_bins/2, σ, 0])

        #i = pulse # kk - index
        #y = view(data, i, on_st:on_end)

        y = view(data, pulse, :)
        (mi, ma) = extrema(y)
        #y = (y .- mi) / (ma - mi)
        y = y ./ ma # this is much much better!

        # convolution with Gaussian
        res = Tools.conv(y, kernel)
        (mi, ma) = extrema(res)
        res = (res .- mi) / (ma - mi)
        re = res[floor(Int,p2_bins/2)+1:end-floor(Int, p2_bins/2)]

        # auto correlation # just for some tests..
        #=
        resa = crosscor(y, y, collect(floor(Int,-length(y)/2):floor(Int,length(y)/2)))
        (mi, ma) = extrema(resa)
        resa = (resa .- mi) / (ma - mi)
        rea = resa[floor(Int,p2_bins/2):end-floor(Int,p2_bins/2)]
        =#

        peak = Tools.peaks(re)
        # new syntax?
        ma, pa = peakprom(Maxima(), re, floor(Int, p2_bins/4))
        #ma, pa = peakprom(re, Maxima(), floor(Int, p2_bins/4))
        #ma, pa = peakprom(re, Maxima(), p2_bins/2)
        inds = sortperm(pa, rev=true)

        peaks = []
        yys = []
        xyys = []
        psnrs = []
        for ii in inds
            st = floor(Int, ma[ii] - p2_bins / 2)
            en = ceil(Int, ma[ii] + p2_bins / 2)
            if (st >= bin_st) && (en <= bin_end)
                #println(y)
                yy = y[st:en]
                signal = Tools.simps(yy , collect(st:en))
                #signal2 = trapz(collect(st:en), yy) # should work
                #yy3 = y[st:en] .- minimum(y[st:en]) # nope
                #signal3 = Tools.simps(yy3 , collect(st:en))

                #nn = y[off_st:off_st+(en-st)]
                nn = y[off_st:off_end]
                noise = std(nn) * (en-st)^0.5
                #nn2 = y[off_st:off_st+(en-st)]
                #noise2 = std(nn2) * (en-st)^0.5
                #nn3 = y[off_st:off_end] .- minimum(y[off_st:off_end]) # nope
                #noise3 = std(nn3) * (en-st)^0.5

                #println("S/N (3)", signal3 / noise3) # nope
                println("S/N ", signal / noise)

                if signal / noise > thresh
                    push!(peaks, ma[ii])
                    push!(psnrs, signal/noise)
                    push!(yys, yy)
                    push!(xyys, collect(range(st/bins*360, en/bins*360, length=length(yy))))
                end
            #if (re[ma[ii]] >= 0.5) &&  (ma[ii] > bin_st) && (ma[ii] < bin_end)
            #    push!(peaks, ma[ii])
            end
        end


        # Pulse longitude
        longitude = collect(range(bin_st/bins * 360, bin_end/bins*360, length=bin_end+1-bin_st))

        # save data for classes
        y = view(data, pulse, bin_st:bin_end)
        (mi, ma) = extrema(y)
        y = y ./ ma # this is much much better!
        f = open("output/pulse.txt", "w")
        write(f, "# longitude, intensity\n")
        for i in 1:length(longitude)
            write(f, "$(longitude[i]) $(y[i])\n")
        end
        close(f)
        ko = re[bin_st:bin_end]
        f = open("output/convolution.txt", "w")
        write(f, "# longitude, intensity\n")
        for i in 1:length(longitude)
            write(f, "$(longitude[i]) $(ko[i])\n")
        end
        close(f)
        f = open("output/data.txt", "w")
        write(f, "Pulse S/N: $(round(Int, snrs[pulse]))\n")
        write(f, "Peaks: [$(peaks[1]/bins*360), $(peaks[2]/bins*360)]\n")
        write(f, "S/Ns: [$(round(psnrs[1], digits=1)), $(round(psnrs[2], digits=1))] ")
        close(f)
        for i in 1:length(yys)
            f = open("output/subpulse_$i.txt", "w")
            write(f, "# longitude, intensity")
            for j in 1:length(xyys[i])
                write(f, "$(xyys[i][j]) $(yys[i][j])\n")
                #plot(xyys[i], yys[i], c="tab:green", lw=2, alpha=0.5)
            end
            close(f)
        end


        println("pulse: ", pulse)
        println("peaks: ", peaks)

        rc("font", size=7.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        #figure(figsize=(3.14961, 2.362205), frameon=true)  # 8cm x 6 cm
        figure(figsize=(3.14961, 1.946563744), frameon=true)  # 8cm x 4.94427191 cm (golden)
        subplots_adjust(left=0.165, bottom=0.19, right=0.99, top=0.99, wspace=0.0, hspace=0.0)

        minorticks_on()
        y = view(data, pulse, bin_st:bin_end)
        (mi, ma) = extrema(y)
        y = y ./ ma # this is much much better!
        plot(longitude, y, c="black", lw=0.5)
        plot(longitude, re[bin_st:bin_end], c="tab:red")
        for i in 1:length(yys)
            plot(xyys[i], yys[i], c="tab:green", lw=2, alpha=0.5)
        end
        #plot(longitude, ress_acorr[2][bin_st:bin_end], c="tab:blue") # for tests only
        #axhline(y=thresh2[2], ls=":", lw=0.7, c="tab:green")
        for p in peaks
            axvline(x=p / bins * 360, c="tab:blue", ls="--")
        end
        xlim(longitude[1], longitude[end])
        yticks([-0.5, 0.0, 0.5])
        text(195, 0.87, "pulse S/N \$\\approx\$ $(round(Int, snrs[pulse]))", size=7)
        for i in 1:length(psnrs)
            text(peaks[i]/bins*360+2, -0.5, "S/N = $(round(psnrs[i], digits=1))", ha="center", size=6)
        end

        ylabel("intensity (a. u.)")
        xlabel("longitude \$(^\\circ)\$")

        filename = "$outdir/$(name_mod)_track_subpulses_new.pdf"
        println(filename)
        savefig(filename)
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        PyPlot.close()

    end


end  # module Plot
