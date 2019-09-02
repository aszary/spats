module Plot
    using PyPlot
    PyPlot.matplotlib.use("Agg")

    include("tools.jl")


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


    function single(data, outdir; start=1, number=100, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0")
        num, bins = size(data)
        if number == nothing
            number = num - start  # missing one?
        end
        if bin_st == nothing bin_st = 1 end
        if bin_end == nothing bin_end = bins end
        da = data[start:start+number-1,bin_st:bin_end]
        average = Tools.average_profile(da)
        intensity, pulses = Tools.intensity_pulses(da)
        pulses .+= start

        # Pulse longitude
        db = (bin_end + 1) - bin_st  # yes +1
        dl = 360. * db / bins
        longitude = collect(range(-dl/2., dl/2., length=db))

        rc("font", size=6.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 4.33071))  # 8cm x 11cm
        subplots_adjust(left=0.16, bottom=0.08, right=0.99, top=0.99, wspace=0., hspace=0.)

        subplot2grid((5, 3), (0, 0), rowspan=4)
        minorticks_on()
        plot(intensity, pulses, c="grey")
        xlim(1.1, -0.1)
        ylim(pulses[1], pulses[end])
        xticks([0.5, 1.0])
        xlabel("intensity")
        ylabel("Pulse number")

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
        savefig("$outdir/single_$name_mod.pdf")
        close()
    end


    function lrfs(data, outdir; start=1, number=256, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0")

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

        rc("font", size=6.)
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
        ax2 = ax.twinx()
        minorticks_on()
        plot(longitude[7:end-9], dphase[6:end-9])
        ax2.xaxis.set_label_position("top")
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
        savefig("$outdir/lrfs_$name_mod.pdf")
        close()
    end


    function p3_evolution(data, outdir; start=1, end_=nothing, step=10, number=256, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0", verbose=false)
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
                fer = abs(pars[3])  # errs[2] # yeap
                p3 = 1 / f
                p3err = maximum([1 / f - 1 / (f +fer), 1 / (f - fer) - 1 / f])
                if verbose == true println("\tP3 = $p3, P3 error = $p3err") end
            catch exc
                p3 = nothing
                p3err = nothing
                if verbose == true println("\tP3 = nothing, P3 error = nothing") end
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

        rc("font", size=6.)
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
        savefig("$outdir/p3_evolution_$name_mod.pdf")
        close()

    end

end  # modul Plot
