module Plot
    using LsqFit
    using FFTW
    using JLD2
    using Statistics
    using StatsBase
    using PyPlot
    using PyCall
    @pyimport matplotlib.patches as patch

    #PyPlot.matplotlib.use("agg") # DOES NOT WORK on ozStar! had to set backend in matplotlib by hand
    PyPlot.matplotlib.use("qt5agg")
    using Peaks
    using Glob
    using SmoothingSplines
    using DataFrames, GLM

    using RCall

    include("tools.jl")
    include("pyrmodule.jl")


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
        savefig("$outdir/single0.pdf")
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

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(2.362205, 3.248031), frameon=true)  # 6cm x 8.25cm
        if number > 999
            subplots_adjust(left=0.235, bottom=0.115, right=0.99, top=0.99, wspace=0., hspace=0.)
        else
            subplots_adjust(left=0.21, bottom=0.115, right=0.99, top=0.99, wspace=0., hspace=0.)
        end
        figtext(0.01, 0.95, "$panel)", size=10)
        subplot2grid((5, 3), (0, 0), rowspan=4)
        minorticks_on()
        plot(intensity, pulses, c="grey")
        ylim(pulses[1]-0.5, pulses[end]+0.5)
        xticks([0.5, 1.0])
        xlim(1.1, -0.1)
        xlabel("intensity", labelpad=7)
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
        xlabel("longitude (deg.)")
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
        phase_ = rad2deg.(angle.(view(lrfs, peak, :)))  # fft phase variation  # view used! lrfs[peak, :] -> copies data
        # skip freq = 0 and normalize intensity to 1
        inten = intensity[2:end]
        inten .-= minimum(inten)
        inten ./= maximum(inten)
        fre = freq[2:end]
        pars, errs = Tools.fit_gaussian(fre, inten; μ=freq[peak-1])  # skip zero freq
        fr = pars[2]
        frer = pars[3]  # errs[2] # yeap

        println("\tFrequancy: $fr, P3: $(1/fr)")
        println("\tFrequancy error: $frer, P3 error: $(1/fr - 1/(fr +frer))")  # TODO err ok?

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
        ylabel("frequancy \$(1/P)\$")

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
        savefig("$outdir/$(name_mod)_lrfs.pdf")
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
        close()
    end


    function p3fold_two(data1, data2, outdir; start=1, number=nothing, repeat_num=4, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0", shift=11)
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
        longitude = collect(range(-dl/2., dl/2., length=db))

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
        imshow(da2, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da2))
        tick_params(labelleft=false, labelbottom=false)
        yticks(ticks, ti)
        # bottom left
        subplot2grid((5, 2), (4, 0))
        minorticks_on()
        plot(longitude, average, c="grey")
        xticks([-40, 0, 40])
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        ylabel("intensity (a. u.)")
        # bottom right
        subplot2grid((5, 2), (4, 1), colspan=3)
        minorticks_on()
        tick_params(labelleft=false)
        plot(longitude, average2, c="grey")
        xticks([-40, 0, 40])
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        figtext(0.45, 0.01, "longitude (deg.)", size=8)
        #tick_params(labeltop=false, labelbottom=true)
        savefig("$outdir/$(name_mod)_p3fold_two.pdf")
        close()
    end


    function p3fold_twotracks(data1, data2, data_dir, outdir; start=1, number=nothing, repeat_num=4, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0", shift=7)
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
        longitude = collect(range(-dl/2., dl/2., length=db))

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
            files = Glob.glob("track_*.jld2", "$(data_dir)$i")
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
            # two tracks
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
        imshow(da2, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da2), extent=[bin_st, bin_end, 1, num*repeat_num])
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
        xticks([-40, 0, 40])
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])
        ylabel("intensity (a. u.)")

        # bottom right
        subplot2grid((5, 2), (4, 1), colspan=3)
        minorticks_on()
        tick_params(labelleft=false)
        plot(longitude, average2, c="grey")
        xticks([-40, 0, 40])
        yticks([0.0, 0.5])
        xlim(longitude[1], longitude[end])

        figtext(0.45, 0.01, "longitude (deg.)", size=8)

        #tick_params(labeltop=false, labelbottom=true)
        savefig("$outdir/$(name_mod)_p3fold_twotracks.pdf")
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
                        track[i, 1] = picks[i][1]
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
                    thisline.set_alpha(0.1)
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
        extent = [bin_st, bin_end, start-0.5, start+number-1+0.5]
        subplot2grid((5, 3), (0, 1), rowspan=4, colspan=2)
        imshow(da, origin="lower", cmap=cmap, interpolation="none", aspect="auto",  vmax=darkness*maximum(da), extent=extent)
        xlim([extent[1], extent[2]])
        ylim([extent[3], extent[4]])
        #println([bin_st, bin_end, start-0.5, start+number+0.5])
        for peak in peaks
            if (peak[1] >= start) && (peak[1] <= start + number)
                for x in peak[2]
                    plot(x, peak[1], marker="x", markersize=2.5, c="red", lw=0, picker=10, alpha=0.3)
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


    function tracks_analysis3(outdir; lambda=100, start=1, number=100, cmap="viridis", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="PSR_NAME")

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



    function driftrate_J1750(outdir; lambda=100, name_mod="123456", show_=false)

        tracks1, lines1, inclines1, dr1, ysp1 = Tools.get_driftrate("$outdir/tracks1", lambda)
        tracks2, lines2, inclines2, dr2, ysp2 = Tools.get_driftrate("$outdir/tracks2", lambda)
        tracks3, lines3, inclines3, dr3, ysp3 = Tools.get_driftrate("$outdir/tracks3", lambda)
        tracks4, lines4, inclines4, dr4, ysp4 = Tools.get_driftrate("$outdir/tracks4", lambda)
        tracks5, lines5, inclines5, dr5, ysp5 = Tools.get_driftrate("$outdir/tracks5", lambda)
        tracks6, lines6, inclines6, dr6, ysp6 = Tools.get_driftrate("$outdir/tracks6", lambda)

        y0 = 401
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
        ylabel("Longitude (\$ ^{\\circ}\$)")
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

        for i in 1:6
            tracks1, lines1, inclines1, dr1, ysp1 = Tools.get_driftrate("$outdir/tracks$i", lambda)
            (ysp, dr) = Tools.remove_duplicates_ysp(ysp1, dr1)
            push!(pulses, dr)
            push!(ysps, ysp)
        end

        ranges = [(3, 0), (0, -3), (3, 1), (1, 0), (0, -1), (-1, -3)]
        profiles = []
        averages = []

        for r in ranges
            push!(profiles, [])
            for i in 1:6
                for (j, driftrate) in enumerate(ysps[i])
                    if (driftrate < r[1]) && (driftrate > r[2])
                        push!(profiles[end], das[i][pulses[i][j], :] )
                    end
                end
            end
        end

        #=
        # make profiles equal
        le = 210
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

        labels = ["\$ \\qquad \\;\\; {\\rm D} > 0\$", "\$\\qquad \\;\\; {\\rm D} < 0 \$", "\$2.15 > {\\rm D} > \\quad \\; 1 \$", "\$\\quad 1\\;\\; > {\\rm D} > \\quad \\; 0 \$", "\$ \\quad 0 \\;\\;  > {\\rm D} > \\;-1 \$", "\$\\!-1 \\;\\; > {\\rm D} > -1.88 \$"]


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
            tracks1, lines1, inclines1, dr1, ysp1 = Tools.get_driftrate("$outdir/tracks$i", lambda)
            push!(tracks, tracks1)
            push!(inclines, inclines1)
            push!(dr_x, dr1)
            push!(dr_y, ysp1)
        end

        ranges = Dict(2=>(510, 590), 3=>(50, 150), 3=>(375, 445), 4=>(100, 150), 5=>(250, 350), 6=>(200, 365))
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
                x, y2, ysl, slopes, eslopes, bp, ebp, xl, exl = PyRModule.segmented(x, y, npsi[nn])
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

        which = 4


        x_, y_, ysl_, y2_, slopes, eslopes, xlin, exlin, bps, ebps = getdata_driftdirection(iis[which], jjs[which], selected_tracks, fitted_lines)
        #println(y2_)
        #println(size(y2_))

        name_mod = "$(iis[1][1])"


        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 1.946563744), frameon=true)  # 8cm x 4.94427191 cm (golden)
        subplots_adjust(left=0.17, bottom=0.19, right=0.99, top=0.99, wspace=0., hspace=0.)

        ax = subplot2grid((2, 1), (0, 0))
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
            for j in 1:length(bps[1])
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
        ax = subplot2grid((2, 1), (1, 0))
        minorticks_on()
        for i in 1:length(slopes)
            errorbar(xlin[i], slopes[i], yerr=eslopes[i], xerr=exlin[i], color="none", lw=0.5, marker="_", mec="grey", ecolor="grey", capsize=0, mfc="grey", ms=1.0)
        end
        plot(dr_x[iis[which][1]], dr_y[iis[which][1]])  # smoothed drift rate
        xlim(xl)
        ylim([-2.3, 2.3])
        savefig("$outdir/$(name_mod)_driftdirection.pdf")
        println("$outdir/$(name_mod)_driftdirection.pdf")
        if show_ == true
            show()
            readline(stdin; keep=false)
        end
        close()

    end



end  # modul Plot
