module Tools
    using Glob
    using JLD2
    using FFTW
    using Peaks
    using LsqFit
    using Statistics
    using SmoothingSplines
    using StatsBase
    using CubicSplines
    using Trapz
    using JSON
    using SavitzkyGolay

    using PyPlot
    using DSP
    using LinearAlgebra
    using DataFrames, GLM

    #include("pyrmodule.jl")

    function rms(data)
        s = 0.0
        for a in data
            s += a * a
        end
        return sqrt(s / length(data))
    end


    function average_profile(data; norm=true)
        pulses, bins = size(data)
        ave = zeros(bins)
        for i in 1:pulses
            for j in 1:bins
                ave[j] += data[i,j]
            end
        end
        if norm == true
            ma = maximum(ave)
            ave = ave / ma
        elseif length(norm) == 2
            me = mean(ave[norm[1]:norm[2]])
            ave .-= me
            ma = maximum(ave)
            ave = ave ./ ma
        end
        return ave
    end


    function intensity_pulses(data)
        pulse_num, bins = size(data)
        intensity = zeros(pulse_num)
        for i in 1:pulse_num
            intensity[i] = sum(data[i, :])
        end
        pulses = collect(1:pulse_num)
        mi = minimum(intensity)
        if mi < 0
            intensity .+= -mi
        end
        intensity /= maximum(intensity)
        return (intensity, pulses)
    end

    function fftfreq(n, d=1.0)
        # one side (positive) frequancy range
        if iseven(n)
            rf = range(0, convert(Int, n/2 - 1), step=1)
        else
            rf = range(0, convert(Int, (n-1)/2 - 1), step=1)  # -1 added to match size
        end
        f = collect(rf) ./ (d * n)
        return f
    end


    function lrfs(data)
        da = transpose(data)
        bins, pulse_num = size(da)
        half = floor(Int, pulse_num / 2) # one side frequency range?
        lrfs = fft(da, 2)[:, 1:half] # second dim! important!
        lrfs = transpose(lrfs)
        freq = fftfreq(pulse_num)
        intensity = zeros(half)
        ab = abs.(lrfs)
        for i in 1:half
            intensity[i] = sum(ab[i,:]) # this is important!
        end
        pk = peaks(intensity)
        return lrfs, intensity, freq, pk
    end


    function periodicity(ydata)
        si = size(ydata, 1)
        # one side (positive) frequancy range
        if iseven(si)
            rf = range(0, convert(Int, si/2 - 1), step=1)
        else
            rf = range(0, convert(Int, (si-1)/2 - 1), step=1)  # -1 added to match size
        end
        freq = collect(rf) ./ si

        half = floor(Int, si / 2)
        ff = fft(ydata)[1:half]
        pk = peaks(abs.(ff))
        #println(freq[pk])

        subplot2grid((2, 1), (0, 0))  # row column
        plot(ydata)
        axvline(x=1/freq[pk], c="red")

        subplot2grid((2, 1), (1, 0))  # row column
        plot(freq, ff)
        #plot([1])
        show()
        readline(stdin; keep=false)
        close()
        return
    end


    function peaks(intensity)
        y = intensity
        pks, vals = findmaxima(y)
        pks, proms = peakproms(pks, y)
        maxprom , ind = findmax(proms)
        return pks[ind] # index of peak with maximum prominence
    end


    """
    module Peaks changed? No Maxima()?
    """
    function peaks_old(intensity)

        #max = maxima(intensity, 30)
        # new syntax?
        ma, pa = peakprom(Maxima(), intensity) #floor(Int, length(intensity)/20))  # works for J1750
        #ma, pa = peakprom(intensity, Maxima(), 1) #floor(Int, length(intensity)/20))  # works for J1750
        #ma, pa = peakprom(intensity, Maxima(), floor(Int, length(intensity)/20))
        #println(ma)
        #println(pa)
        val , ind = findmax(pa)
        return ma[ind]
    end

    #=
    @. # NOPE!
    function gauss(x, p)
        return p[1] * exp(-0.5*((x-p[2])/p[3])^2) + p[4]
    end
    =#

    @. gauss(x, p) = p[1] * exp(-0.5*((x-p[2])/p[3])^2) + p[4]

    function fit_gaussian(xdata, ydata; a=nothing, μ=nothing, σ=nothing, baselevel=nothing)
        if baselevel == nothing
            baselevel = mean(ydata)
        end
        if a == nothing
            a = maximum(ydata) - baselevel
        end
        if μ == nothing
            ind = trunc(Int, length(xdata) / 2)
            μ = xdata[ind]
        end
        if σ == nothing
            σ = xdata[20] - xdata[1]
        end
        # how to use gauss here? done!
        #@. model(x, p) = p[1] * exp(-0.5*((x-p[2])/p[3])^2) + p[4]

        p0 = [a, μ, σ, baselevel]  # look here
        fit = curve_fit(gauss, xdata, ydata, p0)
        p = coef(fit)
        err = stderror(fit)
        #=
        figure()
        plot(xdata, ydata, c="black")
        plot(xdata, gauss(xdata, p0), c="blue")
        plot(xdata, gauss(xdata, p), c="red")
        savefig("output/test.pdf")
        close()
        =#
        return p, err
    end


    function fit_gaussian_J1750(xdata, ydata, start_pulse; a=nothing, μ=nothing, σ=nothing, baselevel=nothing)
        if baselevel == nothing
            baselevel = mean(ydata)
        end
        if a == nothing
            a = maximum(ydata) - baselevel
        end
        if μ == nothing
            val , ind = findmax(ydata) # what?
            μ = xdata[ind]
        end
        if σ == nothing
            σ = xdata[2] - xdata[1]
        end
        # how to use gauss here? done!
        #@. model(x, p) = p[1] * exp(-0.5*((x-p[2])/p[3])^2) + p[4]

        p0 = [a, μ, σ, baselevel]  # look here
        fit = curve_fit(gauss, xdata, ydata, p0)
        p = coef(fit)
        err = stderror(fit)
        if start_pulse >= 17900
            println(start_pulse)
            figure()
            plot(xdata, ydata, c="black")
            plot(xdata, gauss(xdata, p0), c="blue")
            plot(xdata, gauss(xdata, p), c="red")
            q = readline(stdin)
            if q == "q"
                throw(DivideError())
            end
            close()
        end
        return p, err
    end


    @. twogauss(x, p) = p[1] * exp(-0.5*((x-p[2])/p[3])^2)  + p[4] * exp(-0.5*((x-p[5])/p[6])^2) + p[7]

    function fit_twogaussians(xdata, ydata, a1, a2, μ1, μ2, σ1, σ2; baselevel=nothing)
        if baselevel == nothing
            baselevel = median(ydata)
        end

        p0 = [a1, μ1, σ1, a2, μ2, σ2, baselevel]  # look here
        fit = curve_fit(twogauss, xdata, ydata, p0)
        p = coef(fit)
        err = stderror(fit)
        #=
        figure()
        plot(xdata, ydata, c="black")
        plot(xdata, twogauss(xdata, p0), c="blue")
        plot(xdata, twogauss(xdata, p), c="red")
        savefig("output/test.pdf")
        close()
        #println("new")
        #readline(stdin)
        =#
        return p, err
    end


    @. threegauss(x, p) = p[1] * exp(-0.5*((x-p[2])/p[3])^2)  + p[4] * exp(-0.5*((x-p[5])/p[6])^2) + p[7] * exp(-0.5*((x-p[8])/p[9])^2) + p[10]

    function fit_threegaussians(xdata, ydata, a1, a2, a3, μ1, μ2, μ3, σ1, σ2, σ3; baselevel=nothing)
        if baselevel == nothing
            baselevel = median(ydata)
        end
        p0 = [a1, μ1, σ1, a2, μ2, σ2, a3, μ3, σ3, baselevel]  # look here
        fit = curve_fit(threegauss, xdata, ydata, p0)
        p = coef(fit)
        err = stderror(fit)
        #=
        figure()
        plot(xdata, ydata, c="black")
        plot(xdata, threegauss(xdata, p0), c="blue")
        plot(xdata, threegauss(xdata, p), c="red")
        savefig("output/test.pdf")
        close()
        #println("new")
        #readline(stdin)
        =#
        return p, err
    end


    function find_peaks(data)
        sz = size(data)
        peaks = []
        for i in 1:sz[1]
            signal = data[i, :]
            xdata = collect(0:length(signal)-1) # from 0 python plotting
            peak = Tools.peaks(signal)
            ma, pa = peakprom(signal, Maxima(), 3)
            inds = sortperm(pa, rev=true)
            if length(inds) > 2 && (signal[ma[inds[3]]]> 0.5 * signal[ma[inds[2]]])
            #if length(inds) > 2 && ( signal[ma[inds[3]]]> 0.7 * signal[ma[inds[2]]]) # lazy, but works
                a1 = signal[ma[inds[1]]]
                a2 = signal[ma[inds[2]]]
                a3 = signal[ma[inds[3]]]
                try
                    pa, errs = Tools.fit_threegaussians(xdata, signal, a1, a2, a3, xdata[ma[inds[1]]], xdata[ma[inds[2]]], xdata[ma[inds[3]]], 3, 3, 3)
                    push!(peaks, [i-1, pa[2], pa[5], pa[8]])  # python plotting!
                catch
                    pa, errs = Tools.fit_twogaussians(xdata, signal, a1, a2, xdata[ma[inds[1]]], xdata[ma[inds[2]]], 3, 3)
                    push!(peaks, [i-1, pa[2], pa[5]])  # python plotting!
                    #=
                    plot(xdata, signal)
                    axvline(x=xdata[ma[inds[1]]])
                    axvline(x=xdata[ma[inds[2]]])
                    axvline(x=xdata[ma[inds[3]]])
                    axvline(x=pa[2], c="green")
                    axvline(x=pa[5], c="green")
                    axhline(y=a1, c="red")
                    axhline(y=a2, c="red")
                    axhline(y=a3, c="red")
                    savefig("output/test.pdf")
                    close()
                    =#
                end
                #=
                plot(xdata, signal)
                axvline(x=pa[2])
                axvline(x=pa[5])
                axvline(x=pa[8])
                axhline(y=a1, c="red")
                axhline(y=a2, c="red")
                axhline(y=a3, c="red")
                savefig("output/test.pdf")
                close()
                =#
            elseif length(inds) > 1
                a1 = signal[ma[inds[1]]]
                a2 = signal[ma[inds[2]]]
                pa, errs = Tools.fit_twogaussians(xdata, signal, a1, a2, xdata[ma[inds[1]]], xdata[ma[inds[2]]], 3, 3)
                push!(peaks, [i-1, pa[2], pa[5]])  # python plotting!
                #=
                plot(xdata, signal)
                axvline(x=pa[2])
                axvline(x=pa[5])
                axhline(y=a1, c="red")
                axhline(y=a2, c="red")
                savefig("output/test.pdf")
                close()
                =#
            end
        end
        return peaks
    end


    function find_peaks2(data)
        sz = size(data)
        peaks = []
        for i in 1:sz[1]
            signal = data[i, :]
            xdata = collect(0:length(signal)-1) # from 0 python plotting
            peak = Tools.peaks(signal)
            ma, pa = peakprom(signal, Maxima(), 3)
            inds = sortperm(pa, rev=true)
            if length(inds) > 2 && (signal[ma[inds[3]]]> 0.5 * signal[ma[inds[2]]]) # lazy, but works
                push!(peaks, [i-1, xdata[ma[inds[1]]], xdata[ma[inds[2]]], xdata[ma[inds[3]]]])  # python plotting!
            elseif length(inds) > 1
                push!(peaks, [i-1, xdata[ma[inds[1]]], xdata[ma[inds[2]]]])  # python plotting!
            end
        end
        return peaks
    end


    function find_peaks3(data)
        sz = size(data)
        peaks = []
        for i in 1:sz[1]
            signal = data[i, :]
            xdata = collect(1:length(signal))
            #println("$i $inds")
            spl = fit(SmoothingSpline, map(Float64, xdata), signal, 500.0)
            ysp = SmoothingSplines.predict(spl)

            peak = Tools.peaks(ysp)
            ma, pa = peakprom(ysp, Maxima(), 3)
            inds = sortperm(pa, rev=true)

            #push!(peaks, [i, xdata[ma[inds[1]]], xdata[ma[inds[2]]], xdata[ma[inds[3]]]])
            push!(peaks, [i, xdata[ma[inds[1]]]])
            #=
                a1 = signal[ma[inds[1]]]

                figure()
                axvline(x=xdata[ma[inds[1]]], lw=1, c="black")
                plot(xdata, signal)
                plot(xdata, ysp, label="spline")
                legend()
                savefig("output/test.pdf")
                readline(stdin; keep=false)
                close()
            =#

            #=
                a1 = signal[ma[inds[1]]]
                a2 = signal[ma[inds[2]]]
                a3 = signal[ma[inds[3]]]

                println(inds[1:3])

                    plot(xdata, signal)
                    plot(xdata, ysp, label="ysp")
                    axvline(x=xdata[ma[inds[1]]], lw=3)
                    axvline(x=xdata[ma[inds[1]]])
                    axvline(x=xdata[ma[inds[2]]], lw=3)
                    axvline(x=xdata[ma[inds[3]]])
                    axvline(x=pa[2], c="green")
                    axvline(x=pa[5], c="green")
                    axhline(y=a1, c="red")
                    axhline(y=a2, c="red")
                    axhline(y=a3, c="red")
                    #show()  # does not work
                    savefig("output/test.pdf")
                    readline(stdin; keep=false)
                    close()
                =#
                #=
                plot(xdata, signal)
                axvline(x=pa[2])
                axvline(x=pa[5])
                axvline(x=pa[8])
                axhline(y=a1, c="red")
                axhline(y=a2, c="red")
                axhline(y=a3, c="red")
                savefig("output/test.pdf")
                close()
                =#
            #=
            elseif length(inds) > 1
                a1 = signal[ma[inds[1]]]
                a2 = signal[ma[inds[2]]]
                pa, errs = Tools.fit_twogaussians(xdata, signal, a1, a2, xdata[ma[inds[1]]], xdata[ma[inds[2]]], 3, 3)
                push!(peaks, [i-1, pa[2], pa[5]])  # python plotting!
                #=
                plot(xdata, signal)
                axvline(x=pa[2])
                axvline(x=pa[5])
                axhline(y=a1, c="red")
                axhline(y=a2, c="red")
                savefig("output/test.pdf")
                close()
                =#
            end
            =#
        end
        return peaks
    end


    function p2_estimate(data; on_st=450, on_end=700, off_st=100, off_end=350, thresh=3.3, win=6, template_num=nothing)
        pulses, bins = size(data)
        average = average_profile(data)
        avs = []
        for i in on_st:on_end
            #av = Array{Float64}(undef, bins)
            av = zeros(bins)
            for j in 1:pulses
                (ma, ind) = findmax(data[j, :])
                if (ind >= i - win) && (ind <= i + win)
                    #println("$(length(av))")
                    #println("$(length(data[j, on_st:on_end]))")
                    av += data[j, :]
                end
            end
            on_rms = rms(av[on_st:on_end])
            off_rms = rms(av[off_st:off_end])
            if on_rms > thresh * off_rms
                push!(avs, av)
            end
        end

        # normalise avs
        for i in 1:length(avs)
            (mi, ma) = extrema(avs[i])
            #println(mi)
            avs[i] .-= mi
            avs[i] ./= (ma - mi)
        end
        # normalise average
        (mi, ma) = extrema(average)
        average .-= mi
        average ./= (ma - mi)

        println("Avs. num: $(length(avs))")
        dbins = []
        #p2errs = []
        p2s = []
        for av in avs
            p0 = [1.0, 0.5, 460.0, 510.0, 15, 15]
            #p0 = [1.0, 0.5, 430.0, 500.0, 15, 15]
            xdata = collect(on_st:on_end)
            pa, errs = Tools.fit_twogaussians(xdata, av[on_st:on_end], p0[1], p0[2], p0[3], p0[4], p0[5], p0[6])
            dbin = abs(pa[2] - pa[5])
            p2 = dbin / bins * 360
            # nesty hack!
            if (p2 > 16) && p2 < (20) && (pa[1] > 0) && (pa[2] > 0)
                push!(p2s, p2)
            end
            #=
            ga = twogauss(xdata, pa)
            PyPlot.close()
            plot(average, c="grey", lw=0.5)
            plot(av, lw=0.3)
            plot(xdata, ga, lw=0.6, c="red")
            savefig("output/test.pdf")
            show()
            q = readline(stdin; keep=false)
            if q == "q"
                break
            end
            =#


            #p2err = win * 360 / bins * sqrt(1 + (rms(av[off_st:off_end])/maximum(av))^2)
            #push!(p2errs, p2err)
        end

        #p2 = median(p2s)
        p2 = mean(p2s)
        println(p2)
        p2err = std(p2s)
        println("$p2 $p2err")
        #=
        for (i, av) in enumerate(avs)
            println("$i")
            PyPlot.close()
            plot(average, c="black", lw=2)
            plot(av, lw=0.7)
            savefig("output/test.pdf")
            q = readline(stdin; keep=false)
            if q == "q"
                break
            end
        end
        =#
        if template_num == nothing
            return p2, nothing
        else
            return p2, avs[template_num]#[on_st+120:on_end-80] # TODO -50?
        end
    end


    function track_subpulses(data, p2; on_st=350, on_end=650, off_st=20, off_end=320, thresh=2.1, thresh2=0.8)
        pulses, bins = size(data)
        p2_bins = floor(Int, p2 / 360 * bins)
        peaks = [] # [pulse_num, [p1, p2, p3...]]
        ppeaks = [] # [p1, p2, p3...]
        σ = p2_bins / 2 / 2.35482
        kernel = gauss(collect(1:p2_bins), [1, p2_bins/2, σ, 0])
        for i in 1:pulses
            #y = view(data, i, on_st:on_end)
            y = view(data, i, :)
            (mi, ma) = extrema(y)
            y = (y .- mi) / (ma - mi)

            res = conv(y, kernel)
            (mi, ma) = extrema(res)
            res = (res .- mi) / (ma - mi)
            re = res[floor(Int,p2_bins/2):end-floor(Int,p2_bins/2)]

            on = maximum(re[on_st:on_end])
            off = rms(y[off_st:off_end])
            sigma = on / off
            #println("$i $sigma")
            if sigma > thresh
                peak = Tools.peaks(re)
                # new syntax?
                ma, pa = peakprom(Maxima(), re, floor(Int, p2_bins/4))
                #ma, pa = peakprom(re, Maxima(), floor(Int, p2_bins/4))
                #ma, pa = peakprom(re, Maxima(), p2_bins/2)
                inds = sortperm(pa, rev=true)
                ppeaks = []
                for ii in inds
                    if (re[ma[ii]] >= thresh2) &&  (ma[ii] > on_st) && (ma[ii] < on_end)
                        push!(ppeaks, ma[ii])
                    end
                end
                push!(peaks, [i, ppeaks])
                #println("$i $ppeaks")
                #=
                PyPlot.close()
                plot(y, c="black")
                plot(re, c="red")
                plot(kernel, c="blue")
                for ii in inds
                    axvline(x=ma[ii]-1, c="pink")
                end
                axvline(x=ma[inds[1]]-1, c="green") # TODO plotting starts from 0!
                axvline(x=ma[inds[2]]-1, c="magenta") # TODO plotting starts from 0!
                #axvline(x=ma[inds[3]]-1, c="brown") # TODO plotting starts from 0!
                show()
                st = readline(stdin; keep=false)
                if st == "q"
                    break
                end
                =#
            end
        end
        return peaks
    end


    function track_subpulses_snr(data, p2, snrfile; on_st=350, on_end=650, off_st=20, off_end=320, thresh=2.1, thresh2=0.8)
        """ based on adaptive threshold """

        thresh2_old = thresh2

        f = open(snrfile)
        snrs = []
        for line in readlines(f)
            push!(snrs, parse(Float64, line))
        end

        println("obs. SNR (no?): ", round(sum(snrs) / sqrt(size(snrs)[1])))
        println("mean SNR: ", mean(snrs))
        println("median SNR: ", median(snrs))

        pulses, bins = size(data)
        p2_bins = floor(Int, p2 / 360 * bins)
        peaks = [] # [pulse_num, [p1, p2, p3...]]
        ppeaks = [] # [p1, p2, p3...]
        σ = p2_bins / 2 / 2.35482
        kernel = gauss(collect(1:p2_bins), [1, p2_bins/2, σ, 0])
        detected_pulses = 0
        located_subpulses = 0
        for i in 1:pulses
            #y = view(data, i, on_st:on_end)
            y = view(data, i, :)
            (mi, ma) = extrema(y)
            #y = (y .- mi) / (ma - mi)
            y = y ./ ma # this is much much better!

            res = conv(y, kernel)
            (mi, ma) = extrema(res)
            res = (res .- mi) / (ma - mi)
            re = res[floor(Int,p2_bins/2):end-floor(Int,p2_bins/2)]

            #on = maximum(re[on_st:on_end])
            #off = rms(y[off_st:off_end])
            #sigma = on / off
            #println("$i $sigma")
            #println("$i $(snrs[i])")
            snr = snrs[i]

            if (snr > thresh) && ~isnan(re[1])  # why re is sometimes NaN?
                peak = Tools.peaks(re)
                # new syntax?
                ma, pa = peakprom(Maxima(), re, floor(Int, p2_bins/4))
                #ma, pa = peakprom(re, Maxima(), floor(Int, p2_bins/4))
                #ma, pa = peakprom(re, Maxima(), p2_bins/2)
                inds = sortperm(pa, rev=true)
                ppeaks = []
                # get new thresh2 (maximum in off pulse region)
                new_thresh = 0.
                for ii in inds
                    if ((ma[ii] < on_st) || (ma[ii] > on_end)) && (new_thresh < re[ma[ii]])
                        new_thresh = re[ma[ii]]
                    end
                end
                thresh2 = maximum([new_thresh, thresh2_old])
                #println(thresh2)
                for ii in inds
                    if (re[ma[ii]] >= thresh2) &&  (ma[ii] > on_st) && (ma[ii] < on_end)
                        push!(ppeaks, ma[ii])
                    end
                end
                if length(ppeaks) > 0
                    located_subpulses += 1
                end
                push!(peaks, [i, ppeaks])
                detected_pulses += 1
                #println("$i $ppeaks")
                #=
                PyPlot.close()
                plot(y, c="black")
                plot(re, c="red")
                plot(kernel, c="blue")
                #for ii in inds
                for p in ppeaks
                    axvline(x=p-1, c="pink")
                end
                #axvline(x=ma[inds[1]]-1, c="green") # TODO plotting starts from 0!
                #axvline(x=ma[inds[2]]-1, c="magenta") # TODO plotting starts from 0!
                #axvline(x=ma[inds[3]]-1, c="brown") # TODO plotting starts from 0!
                show()
                st = readline(stdin; keep=false)
                if st == "q"
                    break
                end
                =#
            end
        end
        println("Number of pulses: $pulses  Included pulses : $detected_pulses  Pulses with located subpulses: $located_subpulses")
        frac = located_subpulses / pulses * 100
        println("Fraction of pulses used: $(round(frac)) (before grouping (overestimated!))")
        return peaks
    end




    function track_subpulses_template(data, template; on_st=450, on_end=700, off_st=100, off_end=350, thresh=2.1, thresh2=0.8)
        pulses, bins = size(data)
        #pulses = 1
        peaks = [] # [pulse_num, [p1, p2, p3...]]
        ppeaks = [] # [p1, p2, p3...]
        kernel = template
        for i in 1:pulses
            #y = view(data, i, on_st:on_end)
            y = view(data, i, :)
            #y = data
            (mi, ma) = extrema(y)
            y = (y .- mi) / (ma - mi)

            res = conv(y, kernel)
            (mi, ma) = extrema(res)
            res = (res .- mi) / (ma - mi)
            re = res # [floor(Int,p2_bins/2):end-floor(Int,p2_bins/2)]

            peak = Tools.peaks(re)
            ma, pa = peakprom(re, Maxima(), 10)
            inds = sortperm(pa, rev=true)
            #on = maximum(re[on_st:on_end])
            #off = rms(y[off_st:off_end])
            #sigma = on / off
            #println("$i $sigma")
            #if sigma > thresh
            #end
            println("ind. $(ma[inds[1]])")

            PyPlot.close()
            plot(y, c="black")
            plot(re, c="red")
            plot(kernel, c="blue")
            axvline(x=ma[inds[1]]-1, c="green") # TODO plotting starts from 0!
            #axvline(x=ma[inds[3]]-1, c="magenta") # TODO plotting starts from 0!
            show()
            st = readline(stdin; keep=false)
            if st == "q"
                break
            end
        end
        return peaks
    end


    "based on find_peaks3"
    function track_subpulses_smoothing(data; lambda=500)
        pulses, bins = size(data)
        peaks = []
        for i in 1:pulses
            signal = data[i, :]
            xdata = collect(1:length(signal))
            #println("$i $inds")
            spl = fit(SmoothingSpline, map(Float64, xdata), signal, lambda)
            ysp = SmoothingSplines.predict(spl)

            peak = Tools.peaks(ysp)
            ma, pa = peakprom(ysp, Maxima(), 3)
            inds = sortperm(pa, rev=true)

            #push!(peaks, [i, xdata[ma[inds[1]]], xdata[ma[inds[2]]], xdata[ma[inds[3]]]])
            push!(peaks, [i, xdata[ma[inds[1]]]])
            #=
                a1 = signal[ma[inds[1]]]
                figure()
                axvline(x=xdata[ma[inds[1]]], lw=1, c="black")
                plot(xdata, signal)
                plot(xdata, ysp, label="spline")
                legend()
                savefig("output/test.pdf")
                readline(stdin; keep=false)
                if readline(stdin; keep=false) == "q"
                    break
                end
                close()
            =#
            #=
                a1 = signal[ma[inds[1]]]
                a2 = signal[ma[inds[2]]]
                a3 = signal[ma[inds[3]]]

                println(inds[1:3])

                    plot(xdata, signal)
                    plot(xdata, ysp, label="ysp")
                    axvline(x=xdata[ma[inds[1]]], lw=3)
                    axvline(x=xdata[ma[inds[1]]])
                    axvline(x=xdata[ma[inds[2]]], lw=3)
                    axvline(x=xdata[ma[inds[3]]])
                    axvline(x=pa[2], c="green")
                    axvline(x=pa[5], c="green")
                    axhline(y=a1, c="red")
                    axhline(y=a2, c="red")
                    axhline(y=a3, c="red")
                    #show()  # does not work
                    savefig("output/test.pdf")
                    readline(stdin; keep=false)
                    close()
                =#
                #=
                plot(xdata, signal)
                axvline(x=pa[2])
                axvline(x=pa[5])
                axvline(x=pa[8])
                axhline(y=a1, c="red")
                axhline(y=a2, c="red")
                axhline(y=a3, c="red")
                savefig("output/test.pdf")
                close()
                =#
            #=
            elseif length(inds) > 1
                a1 = signal[ma[inds[1]]]
                a2 = signal[ma[inds[2]]]
                pa, errs = Tools.fit_twogaussians(xdata, signal, a1, a2, xdata[ma[inds[1]]], xdata[ma[inds[2]]], 3, 3)
                push!(peaks, [i-1, pa[2], pa[5]])  # python plotting!
                #=
                plot(xdata, signal)
                axvline(x=pa[2])
                axvline(x=pa[5])
                axhline(y=a1, c="red")
                axhline(y=a2, c="red")
                savefig("output/test.pdf")
                close()
                =#
            end
            =#
        end
        return peaks
    end


    @. fourgauss(x, p) = p[1] * exp(-0.5*((x-p[2])/p[3])^2)  + p[4] * exp(-0.5*((x-p[5])/p[6])^2) + p[7] * exp(-0.5*((x-p[8])/p[9])^2) + p[10] * exp(-0.5*((x-p[11])/p[12])^2) + p[13]

    function template(data, p2; on_st=450, on_end=700, off_st=100, off_end=350, dbins=LinRange(-0.55, -0.6, 100))
        pulses, bins = size(data)
        te = Array{Float64}(undef, pulses, bins)
        sigma_max = -1e50
        dbin_best = nothing
        singlepulses = nothing
        for dbin in dbins
            for i in 1:pulses
                for j in 1:bins
                    try
                        te[i, j] = data[i, j-trunc(Int, i*dbin+0.5)]
                    catch err
                        if j-trunc(Int, i*dbin+0.5) <= 0
                            te[i, j] = data[i, end+j-trunc(Int, i*dbin+0.5)]
                        else
                            te[i, j] = data[i, j-trunc(Int, i*dbin+0.5)-bins]
                        end
                    end
                end
            end
            on = maximum(te[on_st:on_end])
            off = rms(te[off_st:off_end])
            sigma = on / off
            println("$dbin $sigma")
            if sigma > sigma_max
                sigma_max = sigma
                singlepulses = copy(te)
                dbin_best = dbin
            end
        end
        println("sigma = $sigma_max dbins = $dbin_best")
        template = average_profile(singlepulses)

        #p0 = [0.5, 510.0, 15, 1.0, 560.0, 15, 0.]
        #xdata = collect(1:length(template))
        #pa, errs = Tools.fit_twogaussians(xdata, template, p0[1], p0[2], p0[3], p0[4], p0[5], p0[6])
        #ga = twogauss(xdata, p0)

        #p0 = [a1, μ1, σ1, a2, μ2, σ2, a3, μ3, σ3, baselevel]  # look here
        #p0 = [0.5, 510.0, 15, 1.0, 560.0, 15, 0.2, 610., 15, 0, 0]
        #xdata = collect(1:length(template))
        #fit = curve_fit(threegauss, xdata, template, p0)
        #p = coef(fit)
        #err = stderror(fit)
        #ga = threegauss(xdata, p0)
        #=
        p0 = [0.5, 510.0, 15, 1.0, 560.0, 15, 0.2, 610., 15, 0.1, 450, 15, 0]
        xdata = collect(1:length(template))
        fit = curve_fit(fourgauss, xdata, template, p0)
        p = coef(fit)
        err = stderror(fit)

        ga = fourgauss(xdata, p)

        println(p)
        plot(template)
        plot(ga)
        show()
        readline(stdin; keep=false)
        close()

        template = twogauss(xdata, [1.0, p[2], 17, 1.0, p[5], 17, 0])
        =#
        return template, singlepulses
    end


    function group_tracks_obsolete(peaks, p2; bins=1024)
        # peaks [pulse_num, [peak1, peak2, ...]]
        pulses = size(peaks)[1]
        p2_bins = floor(Int, p2 / 360 * bins)
        # first points
        #tracks = [[peaks[1][1], peaks[1][2]]]
        #tracks = [[[peaks[1][1]], [p]] for p in peaks[1][2]]
        #println(tracks)
        println(p2_bins)
        tracks = []
        for i in 1:1
            for peak in peaks[i][2]#[2] # TODO
                track = [[peaks[i][1]], [peak]]
                for j in i+1:pulses
                    x1 = track[2][end]
                    y1 = track[1][end]
                    y2 = peaks[j][1]
                    dls = []
                    pe = []
                    for peak2 in peaks[j][2]
                        x2 = peak2
                        dl = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2)
                        push!(dls, dl)
                        push!(pe, peak2)
                        println("$j ($x1, $y1) ($x2, $y2) $dl")
                    end
                    dl, ind = findmin(dls)
                    if dl < p2_bins / 3
                        push!(track[1], j)
                        push!(track[2], pe[ind])
                    elseif dl > 2*p2_bins
                        break
                    end
                end
                #println("$i $track")
                push!(tracks, copy(track))
            end
        end
        return tracks
    end

    "using lines fitting"
    function analyse_track(x, y; win=10)
        @. f(x, p) = p[2] * x + p[1]
        lines = []
        inclines = []
        len = length(x)
        for i in 1:len-win
            x_ = x[i:i+win]
            y_ = y[i:i+win]
            data = DataFrame(X=x_, Y=y_)
            ols = lm(@formula(Y ~ X), data)
            line = [x_, f(x_, coef(ols))]
            push!(lines, line)
            # integer pulse number # no need?!
            #push!(inclines, [trunc(Int, mean(x_)+0.5), coef(ols)[2], stderror(ols)[2]])
            push!(inclines, [mean(x_), coef(ols)[2], stderror(ols)[2]])

        end
        return lines, inclines
    end

    """ fit line and get inclinations for data in range i:j - used by find_inclination"""
    function fit_line(x, y, i, j)
        @. f(x, p) = p[2] * x + p[1]
        lines = []
        inclines = []
        x_ = x[i:j]
        y_ = y[i:j]
        data = DataFrame(X=x_, Y=y_)
        ols = lm(@formula(Y ~ X), data)
        line = [x_, f(x_, coef(ols))]
        push!(lines, line)
        push!(inclines, [mean(x_), coef(ols)[2], stderror(ols)[2]])
        return lines, inclines
    end


    function find_inclination(x, y; )
        lines = []
        inclines = []
        len = size(x)[1]

        # find breakpoints
        break_inds = []
        for i in 1 : len-1
            dy = abs(y[i+1] - y[i])
            if dy > 200
                push!(break_inds, i)
                #println(dy, i+1)
            end
        end

        # so lazy
        if length(break_inds) == 0
            li, inc = fit_line(x, y, 1, len)
            lines = vcat(lines, li)
            inclines = vcat(inclines, inc)
        end

        if length(break_inds) == 1
            li, inc = fit_line(x, y, 1, break_inds[1])
            lines = vcat(lines, li)
            inclines = vcat(inclines, inc)

            li, inc = fit_line(x, y, break_inds[1]+1, len)
            lines = vcat(lines, li)
            inclines = vcat(inclines, inc)
        end

        if length(break_inds) > 1
            li, inc = fit_line(x, y, 1, break_inds[1])
            lines = vcat(lines, li)
            inclines = vcat(inclines, inc)
            # Not tested yet
            for j in 1:length(break_inds)-1
                li, inc = fit_line(x, y, break_inds[j], break_inds[j+1])
                lines = vcat(lines, li)
                inclines = vcat(inclines, inc)
            end
            li, inc = fit_line(x, y, break_inds[end]+1, len)
            lines = vcat(lines, li)
            inclines = vcat(inclines, inc)
        end
        return lines, inclines
    end


    """ connects phase first then makes a fit"""
    function find_inclination2(x, y; )
        len = size(x)[1]

        # continuous phase change
        for i in 1 : len-1
            dy = y[i+1] - y[i]
            while abs(dy) > 200
                if dy > 0
                    y[i+1] -= 360
                else
                    y[i+1] += 360
                end
                dy = y[i+1] - y[i]
            end
        end
        #println(y)
        line, incline = fit_line(x, y, 1, len)
        return line, incline
    end




    "using Smooth lines"
    function analyse_track_sl(x, y; lambda=10)
        spl = fit(SmoothingSpline, x, y, lambda)
        ysp = SmoothingSplines.predict(spl)

        line = [x, ysp]
        df = diff(ysp) # this one works strange (peaks)
        x2 = Array{Float64}(undef, length(df))
        for i in 1:length(x)-1
            x2[i] = (x[i+1] + x[i]) / 2.0
        end
        inclines = [x2, df]
        return line, inclines
    end


    function analyse_track_simple(x, y; lambda=10)
        spl = fit(SmoothingSpline, x, y, lambda)
        ysp = SmoothingSplines.predict(spl)

        line = [x, ysp]
        x2 = Array{Float64}(undef, length(y)-1)
        dr = Array{Float64}(undef, length(y)-1)
        for i in 1:length(x)-1
            dx = x[i+1] - x[i]
            dy = ysp[i+1] - ysp[i]
            x2[i] = (x[i+1] + x[i]) / 2.0
            dr[i] = dy / dx
        end
        inclines = [x2, dr]
        return line, inclines
    end


    function analyse_track_simple_dof(x, y, spar)
        #spl = fit(SmoothingSpline, x, y, lambda)
        #ysp = SmoothingSplines.predict(spl)
        #ysp, dof = PyRModule.smooth_spline(x, y; df_start=trunc(Int, length(x)*0.3))
        ysp, dof = PyRModule.smooth_spline(x, y, spar)
        line = [x, ysp]
        x2 = Array{Float64}(undef, length(y)-1)
        dr = Array{Float64}(undef, length(y)-1)
        for i in 1:length(x)-1
            dx = x[i+1] - x[i]
            dy = ysp[i+1] - ysp[i]
            x2[i] = (x[i+1] + x[i]) / 2.0
            dr[i] = dy / dx
        end
        inclines = [x2, dr]
        return line, inclines, dof
    end


    "Why there are duplicates in tracks?"
    function remove_duplicates(track)
        tr = [[track[1,1], track[1,2]]]
        s = size(track)
        for i in 2:s[1]
            if (track[i, 1] != tr[end][1]) &&  (track[i, 2] != tr[end][2])
                push!(tr, [track[i, 1], track[i, 2]])
            end
        end
        # there has to be a batter way
        sz = size(tr)[1]
        tr_final = Array{Float64}(undef, sz, 2)
        for i in 1:sz
            tr_final[i, 1] = tr[i][1]
            tr_final[i, 2] = tr[i][2]
        end
        #println(typeof(tr_final))
        #println(size(tr_final))
        return tr_final
    end



    "Why there are duplicates in in ysp and dr?"
    function remove_duplicates_ysp(ysp, dr)
        ysp_ = []
        dr_ = []

        for (i, pulse) in enumerate(dr)
            try
                prev = dr_[end]
                if trunc(Int64, pulse) != prev
                    push!(ysp_, ysp[i])
                    push!(dr_, trunc(Int64, pulse))
                end
            catch ex
                if typeof(ex) == BoundsError
                    push!(ysp_, ysp[i])
                    push!(dr_, trunc(Int64, pulse))
                end
            end
        end
        return ysp_, dr_
    end


    function get_driftrate(trackdir, lambda)
        tracks = []
        # load tracks
        files = Glob.glob("track_*.jld2", trackdir)
        for file in files
            @load file track
            tr = Tools.remove_duplicates(track)
            push!(tracks, tr)
        end

        lines = []
        dofs = []
        inclines = []
        drift_rate = [[], []]
        for track in tracks#[1:3]
            #ll, inc = Tools.analyse_track_sl(track[:,3], track[:,1]; lambda=lambda)
            ll, inc = Tools.analyse_track_simple(track[:,2], track[:,1]; lambda=lambda)
            #ll, inc, dof = Tools.analyse_track_simple_dof(track[:,2], track[:,1]; lambda=lambda) # do not use it yet!
            dof = 1  # do not use it yet
            push!(lines, ll)
            push!(inclines, inc)
            push!(dofs, dof)
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
        spl = SmoothingSplines.fit(SmoothingSpline, dr1, dr2, lambda)
        #spl = SmoothingSplines.fit(SmoothingSpline, dr1, dr2, 1000.0)
        ysp = SmoothingSplines.predict(spl)
        #println(length(spl.g), " ", length(ysp), " ", length(dr1), " ", spl.Xcount) # works here? does not change with lambda (100)?
        return tracks, lines, inclines, dr1, ysp, dofs
    end


    function get_driftrate2(trackdir, spar)
        tracks = []
        # load tracks
        files = Glob.glob("track_*.jld2", trackdir)
        for file in files
            @load file track
            tr = Tools.remove_duplicates(track)
            push!(tracks, tr)
        end

        lines = []
        dofs = []
        inclines = []
        drift_rate = [[], []]
        for track in tracks#[1:3]
            #ll, inc = Tools.analyse_track_sl(track[:,3], track[:,1]; lambda=lambda)
            #ll, inc = Tools.analyse_track_simple(track[:,2], track[:,1]; lambda=lambda)
            ll, inc, dof = Tools.analyse_track_simple_dof(track[:,2], track[:,1], spar)
            dof = 1  # do not use it yet
            push!(lines, ll)
            push!(inclines, inc)
            push!(dofs, dof)
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
        ysp, dof_ = PyRModule.smooth_spline(dr1, dr2, spar)
        return tracks, lines, inclines, dr1, ysp, dofs
    end



    function driftrate_analysis_J1750(outdir, lambda)

        tracks1, lines1, inclines1, dr1, ysp1 = Tools.get_driftrate("$outdir/tracks1", lambda)
        tracks2, lines2, inclines2, dr2, ysp2 = Tools.get_driftrate("$outdir/tracks2", lambda)
        tracks3, lines3, inclines3, dr3, ysp3 = Tools.get_driftrate("$outdir/tracks3", lambda)
        tracks4, lines4, inclines4, dr4, ysp4 = Tools.get_driftrate("$outdir/tracks4", lambda)
        tracks5, lines5, inclines5, dr5, ysp5 = Tools.get_driftrate("$outdir/tracks5", lambda)
        tracks6, lines6, inclines6, dr6, ysp6 = Tools.get_driftrate("$outdir/tracks6", lambda)

        tracks = []
        lines = []
        inclines = []
        pulses = []
        ysps = []

        for i in 1:6
            tracks1, lines1, inclines1, dr1, ysp1 = Tools.get_driftrate("$outdir/tracks$i", lambda)
            push!(tracks, tracks1)
            push!(lines, lines1)
            push!(inclines, inclines1)
            push!(pulses, dr1)
            push!(ysps, ysp1)
        end

        min_length = 3
        neg_duration = []
        pos_duration = []
        drn_ = []
        drp_ = []
        for i in 1:6
            negative = []
            positive = []
            push!(neg_duration,[])
            push!(pos_duration,[])
            for (j, ysp) in enumerate(ysps[i])
                if signbit(ysp)
                    push!(negative, pulses[i][j])
                    push!(drn_, ysp)
                else
                    push!(positive, pulses[i][j])
                    push!(drp_, ysp)
                end
            end
            # negative duration
            if length(negative) >= 1
                st = negative[1]
                en = negative[end]
            end
            for j in 1:length(negative)-1
                if (negative[j+1] - negative[j] > min_length)
                    push!(neg_duration[end], Dict(trunc(st)=>trunc((negative[j]-st))))
                    st = negative[j+1]
                end
            end
            if length(negative) >= 1
                push!(neg_duration[end], Dict(trunc(st)=>trunc(en-st)))
            end
            # positive duration
            if length(positive) >= 1
                stp = positive[1]
                enp = positive[end]
            end
            for j in 1:length(positive)-1
                if (positive[j+1] - positive[j] > min_length)
                    push!(pos_duration[end], Dict(trunc(stp)=>trunc((positive[j]-stp))))
                    stp = positive[j+1]
                end
            end
            if length(positive) >= 1
                push!(pos_duration[end], Dict(trunc(stp)=>trunc(enp-stp)))
            end
        end

        # negatives
        #println(neg_duration)
        nnum = 0
        ndurations = []
        for (i, neg) in enumerate(neg_duration)
            nnum += length(neg)
            for n in neg
                push!(ndurations, collect(values(n))[1])
            end
        end

        # positives
        println(pos_duration)
        pnum = 0
        pdurations = []
        for (i, pos) in enumerate(pos_duration)
            pnum += length(pos)
            for p in pos
                push!(pdurations, collect(values(p))[1])
            end
        end

        println("Negative instances: ", nnum)
        println("Longest negative: ", maximum(ndurations))
        println("Smallest negative: ", minimum(drn_))

        println("Positive instences: ", pnum)
        println("Positive longest: ", maximum(pdurations))
        println("Biggest positive: ", maximum(drp_))


        return ndurations, pdurations

    end


    function driftrate_analysis_J1750_2(outdir, lambda)

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

        min_length = 3 # minimum mode length in periods
        neg_duration = []
        pos_duration = []
        drn_ = []
        drp_ = []

        null_pulses = []
        for i in 1:length(tracks)
            negative = []
            positive = []
            push!(neg_duration,[])
            push!(pos_duration,[])
            for (j, ysp) in enumerate(ysps[i])
                if signbit(ysp)
                    push!(negative, pulses[i][j])
                    push!(drn_, ysp)
                else
                    push!(positive, pulses[i][j])
                    push!(drp_, ysp)
                end
            end
            # negative duration
            if length(negative) >= 1
                st = negative[1]
                en = negative[end]
            end
            for j in 1:length(negative)-1
                # if no continuity in pulse number
                if (negative[j+1] - negative[j] > min_length)
                    push!(neg_duration[end], Dict(trunc(st)=>trunc((negative[j]-st))))
                    #println("$i ", neg_duration[end][end])
                    st = negative[j+1]
                end
            end
            if (length(negative) >= 1) && (en-st > 0)
                push!(neg_duration[end], Dict(trunc(st)=>trunc(en-st)))
                #println("$i ", neg_duration[end][end])
            end
            # positive duration
            if length(positive) >= 1
                stp = positive[1]
                enp = positive[end]
            end
            for j in 1:length(positive)-1
                if (positive[j+1] - positive[j] > min_length)
                    push!(pos_duration[end], Dict(trunc(stp)=>trunc((positive[j]-stp))))
                    #println("$i ", pos_duration[end][end])
                    stp = positive[j+1]
                end
            end
            if (length(positive)) >= 1 && (enp-stp > 0)
                push!(pos_duration[end], Dict(trunc(stp)=>trunc(enp-stp)))
            end
            # check pulse continuity (nulls?)
            push!(null_pulses, [])
            all = vcat(negative, positive)
            sort!(all)
            nos = 0
            for j in 1:length(all)-1
                if all[j+1] - all[j] > 1
                    push!(null_pulses[end], all[j+1])
                    nos += 1
                    println("Session: $i " ," pulse: ", all[j], " dp: ", all[j+1] - all[j])
                end
            end
            println("\n Session: $i Coverage: ", (1- nos/all[end]) * 100, " Pulses: ", all[end]-nos, " Null fraction: ", nos / all[end], "\n")

        end

        # negatives # based on single pulses information (problem with continuity for low detections!)
        nnum = 0
        ndurations = []
        mimi = 1e50
        mama = -1e50
        nsession = nothing
        nsession2 = nothing
        npulse = nothing
        npulse2 = nothing
        for (i, neg) in enumerate(neg_duration)
            nnum += length(neg)
            for n in neg
                #println(n)
                val = collect(values(n))[1]
                if val < mimi
                    nsession = i
                    #println("$i ", n)
                    npulse = collect(keys(n))[1]
                    mimi = val
                end
                if val > mama
                    nsession2 = i
                    #println("$i ", n)
                    npulse2 = collect(keys(n))[1]
                    mama = val
                end
                push!(ndurations, val)
            end
        end
        #println(ndurations)

        # positives  # based on single pulses information (problem with continuity for low detections!)
        pnum = 0
        pdurations = []
        mama = -1e50
        psession = nothing
        ppulse = nothing
        for (i, pos) in enumerate(pos_duration)
            pnum += length(pos)
            for p in pos
                val = collect(values(p))[1]
                if val > mama
                    psession = i
                    ppulse = collect(keys(p))[1]
                    mama = val
                end
                push!(pdurations, val)
            end
        end

        # negative / positive - based on ysps
        fix_driftrates!(pulses, ysps)

        directions = []
        for i in 1:size(pulses)[1]
            push!(directions, [])
            for j in 1:size(pulses[i])[1]
                if ysps[i][j] > 0
                    push!(directions[end], 1)
                else
                    push!(directions[end], -1)
                end
            end
        end

        pdur = []
        ndur = []
        ptime = [] # [session, start, end]
        ntime = [] # [session, start, end]
        for i in 1:size(directions)[1]
            #push!(pdur, [])
            #push!(ndur, [])
            start = pulses[i][1]
            dir = directions[i][1]
            for j in 2:size(directions[i])[1]
                if dir != directions[i][j]
                    ln = pulses[i][j] - start
                    if dir < 0
                        push!(ndur, ln)
                        push!(ntime, [i, start, pulses[i][j]])
                    else
                        push!(pdur, ln)
                        push!(ptime, [i, start, pulses[i][j]])
                    end
                    dir = directions[i][j]
                    start = pulses[i][j]
                end
             end
            # solve all positive / negative sessions
            if (start == pulses[i][1])
                #println([i, start, pulses[i][end]])
                if directions[i][1] > 0
                    push!(ptime, [i, start, pulses[i][end]])
                    push!(pdur, pulses[i][end]-start)
                else
                    push!(ntime, [i, start, pulses[i][end]])
                    push!(ndur, pulses[i][end]-start)
                end
            end
        end

        #=
        set = 1
        plot(pulses[set], ysps[set], lw=4)
        plot(pulses[set], directions[set], lw=3)
        show()
        readline(stdin; keep=false)
        close()


        hist(ndur, alpha=0.7)
        hist(pdur, alpha=0.7)
        show()
        readline(stdin; keep=false)
        close()
        =#
        ysps_flat = Base.Iterators.flatten(ysps)
        ysps_neg = []
        ysps_pos = []
        for ys in ysps_flat
            if ys < 0
                push!(ysps_neg, ys)
            else
                push!(ysps_pos, ys)
            end
        end

        #println("(in p.) OLD Negative instances: ", nnum)
        println("Negative instances: ", length(ndur))
        #println("(in p.) OLD Longest negative: ", maximum(ndurations), " for session ", nsession2, " pulse ", npulse2)
        println("Longest negative: ", maximum(ndur), " for session ", ntime[findmax(ndur)[2]][1], " pulse ", ntime[findmax(ndur)[2]][2])
        #println("(in p.) OLD shortest negative: ", minimum(ndurations), " for session ", nsession, " pulse ", npulse)
        println("Shortest negative: ", minimum(ndur), " for session ", ntime[findmin(ndur)[2]][1], " pulse ", ntime[findmin(ndur)[2]][2])
        println("(in p.) Smallest negative: ", minimum(drn_))
        #println("Smallest negative: ", minimum(ysps_flat))
        println("(in p.) mean neg.: ", mean(drn_))
        #println("mean neg: ", mean(ysps_neg))
        println("std neg.: ", std(drn_))
        println("(in p.) Neg. Standard error of the mean: ", std(drn_) / sqrt(length(drn_)))
        println()

        #println("(in p.) OLD Positive instences: ", pnum)
        println("Positive instances: ", length(pdur))
        #println("(in p.) OLD Positive longest: ", maximum(pdurations), " for session ", psession, " pulse ", ppulse)
        println("Longest positive: ", maximum(pdur), " for session ", ptime[findmax(pdur)[2]][1], " pulse ", ptime[findmax(pdur)[2]][2])
        #println("(in p.) OLD shortest positive: ", minimum(pdurations))
        println("Shortest positive: ", minimum(pdur), " for session ", ptime[findmin(pdur)[2]][1], " pulse ", ptime[findmin(pdur)[2]][2])
        println("(in p.) Biggest positive: ", maximum(drp_))
        #println("Biggest positive: ", maximum(ysps_flat))
        println("(in p.) mean pos.: ", mean(drp_))
        #println("mean pos.: ", mean(ysps_pos))
        println("std pos.: ", std(drp_))
        println("(in p.) Pos. Standard error of the mean: ", std(drp_) / sqrt(length(drp_)))

        #return ndurations, pdurations, null_pulses
        return ndur, pdur, nothing

    end

    "fixes driftrates arrays to be continuous, changes both pulses and ysps"
    function fix_driftrates!(pulses, ysps)
        for j in 1:size(pulses)[1]
            ysps_ = [ysps[j][1]]
            pulses_ = [pulses[j][1]]

            for i in 2:(size(pulses[j])[1]-1)
                if pulses[j][i] != pulses_[end]
                    push!(pulses_, pulses[j][i])
                    push!(ysps_, ysps[j][i])
                end
            end
            splines = CubicSpline(pulses_, ysps_)
            xs = range(pulses_[1], length=trunc(Int, pulses_[end]-pulses_[1]))
            ys = splines[xs]
            pulses[j] = xs
            ysps[j] = ys
        end
    end


    function driftdirection_analysis_J1750_3(slopes, eslopes, bps, ebps, xs, exs, set)

        positive = []
        negative = []
        pdrate = []
        ndrate = []

        for i in 1:length(xs[set])
            if slopes[set][i] > 0
                for pulse in trunc(Int, xs[set][i]-exs[set][i]):1:trunc(Int,xs[set][i] + exs[set][i])
                    if ~(pulse in positive)
                        push!(positive, pulse)
                        push!(pdrate, slopes[set][i])
                    end
                end
            else
                for pulse in trunc(Int, xs[set][i]-exs[set][i]):1:trunc(Int,xs[set][i] + exs[set][i])
                    if ~(pulse in negative)
                        push!(negative, pulse)
                        push!(ndrate, slopes[set][i])
                    end
                end
            end
        end

        sort!(positive)
        sort!(negative)

        rejected = []
        pos = [] # only positive
        neg = [] # only negative
        pdr = [] # only negative
        ndr = [] # only negative

        # add pos and mark rejections
        for (i, pulse) in enumerate(positive)
            if pulse in negative
                push!(rejected, pulse)
                #println("Set: ", set, " pulse: ", pulse, " pos./neg.")
            else
                push!(pos, pulse)
                push!(pdr, pdrate[i])
            end
        end

        # add neg
        for (i, pulse) in enumerate(negative)
            if ~(pulse in rejected)
                push!(neg, pulse)
                push!(ndr, ndrate[i])
            end
        end

        # POSITIVE
        breaks = []
        lengths = []

        # shortest positive
        start_pulse = pos[1]
        added = 0
        for i in 1:(length(pos)-1)
            if pos[i+1]-pos[i] != 1
                added += 1
                push!(lengths, pos[i]-start_pulse + 1) # yes +1 here
                start_pulse = pos[i+1]
                push!(breaks, i)
            end
        end
        if added == 0 # first observing session
            push!(lengths, pos[end]-pos[1] + 1) # yes +1 here
        end

        println("\n$set")
        println("\tpositive lengths: ", lengths)
        println("\tpositive breaks: ", pos[breaks])
        println("\tPositive: ", length(positive), " Negative: ", length(negative), " Rejected: ",length(rejected))
        println("Breakpoints\n")
        for i in 1:length(bps[set])
            println("\t\t$i ", bps[set][i])
        end
        println("positive drift rate for $set session: ", mean(pdr), " +/- ", std(pdr))

        # NEGATIVE
        nbreaks = []
        nlengths = []
        # shortest positive
        if length(neg) > 0
            start_pulse = neg[1]
            for i in 1:(length(neg)-1)
                if neg[i+1]-neg[i] != 1
                    push!(nlengths, neg[i]-start_pulse + 1) # yes +1 here
                    start_pulse = neg[i+1]
                    push!(nbreaks, i)
                end
            end
        end

        println("\n$set")
        println("\tnegative lengths: ", nlengths)
        println("\tnegative breaks: ", neg[nbreaks])
        println("\tPositive: ", length(positive), " Negative: ", length(negative), " Rejected: ",length(rejected))
        println("Breakpoints\n")
        for i in 1:length(bps[set])
            println("\t\t$i ", bps[set][i])
        end
        return nlengths, ndr, lengths, pdr

    end



    function squared_error(ydata, yfit)
        return sum((yfit .- ydata) .* (yfit .- ydata))
    end


    """ https://pythonprogramming.net/how-to-program-r-squared-machine-learning-tutorial/ """
    function rsquared(ydata, yfit)
        mean_line = [mean(ydata) for y in ydata]
        squared_error_fit = squared_error(ydata, yfit)
        squared_error_mean = squared_error(ydata, mean_line)
        return 1 - (squared_error_fit / squared_error_mean)
    end


    function chisquare(ydata, yfit)
        return sum((ydata .- yfit) .* (ydata .- yfit) ./ ydata)  # there may be some dot problem?
    end


    # https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
    function chisquare_reduced(o, c, dof)
        return sum((o .- c) .^ 2 ./ var(o)) / dof
    end


    function residuals(ydata, yfit)
        return ydata .- yfit

    end

    function p3evolution(data, step, number, bin_st, bin_end; verbose=false)
        num, bins = size(data)
        intensity_ = []
        p3_ = []
        p3_err_ = []
        start_period = []
        for i in 1:step:(num-number)
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
                pars, errs = Tools.fit_gaussian_J1750(fre, inten, 12345) #; μ=freq[peak-1])   # 12345 wtf?!
                f = pars[2]
                #fer = abs(pars[3])  # nope too big
                fer = abs(errs[2])  # yeap
                p3 = 1 / f
                p3err = maximum([1 / f - 1 / (f +fer), 1 / (f - fer) - 1 / f])
                # try other approach to estimate error if too big
                if p3err > 30
                    fer = abs(pars[3])
                    p3err2 = maximum([1 / f - 1 / (f +fer), 1 / (f - fer) - 1 / f])
                    if p3err2 < p3err
                        p3err = p3err2
                    end
                    #println("new $p3err $p3err2")
                end
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

        # converting intensity why why why?
        x, = size(intensity_)
        y, = size(intensity_[1])
        intens = zeros((x,y))
        for i in 1:x
            for j in 1:y
                intens[i,j] = intensity_[i][j]
            end
        end

        return intens, p3_, p3_err_, frequency
    end


    function convert_tracks(tracks_dir; bin=1024)
        for i in 1:7
            dir_ = "$tracks_dir/$i/"
            files = Glob.glob("track_*.jld2", dir_)
            for file in files
                @load file track
                #println(track)
                for i in 1:size(track)[1]
                    track[i, 1] = track[i, 1] / bin * 360
                end
                #@save file track # disabled just to be safe
            end
        end
    end


    """ analyses profile stability using
        https://ui.adsabs.harvard.edu/abs/1975ApJ...198..661H/abstract """
    function analyse_profile_stability(pulses, average, bin_num, binoff_st, binoff_end, binon_st, binon_end, longitude)
        # pulses = profiles[3]
        #average = averages[3]

        # calculate stability for the whole observation
        pulse_num = size(pulses)[1]
        m = trunc(Int, log2(pulse_num / 2))
        #m += 1 # nope
        #avers = Array{Float64}(undef, m, bins[1])
        avers = []
        for i in 0:m
            push!(avers, [])
        end
        for i in 0:m
            mm = 2^i
            av = zeros(bin_num)
            for j in 1:pulse_num
                for k in 1:bin_num
                    av[k] += pulses[j, k]
                end
                if j % mm == 0
                    # normalizes here
                    #m1 = mean(averages[3][binoff_st:binoff_end])
                    m2 = mean(av[binoff_st:binoff_end])
                    av .-= m2
                    # same maximum
                    #maa = maximum(av)
                    #av = av ./ maa
                    # not the maximum, but mean
                    m3 = mean(average[binon_st:binon_end])
                    m4 = mean(av[binon_st:binon_end])
                    av = av .* (m3 / m4)
                    #=
                    m5 = mean(av[binon_st:binon_end])
                    m6 = mean(av[binoff_st:binoff_end])
                    println("$m3 $m4 $m5 $m6")
                    plot(longitude, average, c="tab:blue", lw=0.7)
                    plot(longitude, av, c="tab:red", lw=0.7)
                    show()
                    r = readline(stdin; keep=false)
                    close()
                    if r == "q"
                        return
                    end
                    =#
                    push!(avers[i+1], av)
                    #println("$i $j $(m+1)")
                    av = zeros(bin_num)

                end
            end
        end
        onemrhos = []
        onemrhose = []
        npulses = []
        for i in 1:size(avers)[1]
            push!(npulses, 2^(i-1))
            omrhos = []
            for j in 1:size(avers[i])[1]
                rho = cor(avers[i][j], average)
                #rho = crosscor(avers[i][j], average, [0])[1] # the same yep!
                push!(omrhos, 1 - rho)
            end
            push!(onemrhos, mean(omrhos))
            push!(onemrhose, std(omrhos))
        end
        return onemrhos, onemrhose, npulses
    end

    """ SNR calculation using psrchive """
    function generate_snr(filename)
        #println(filename)
        # get pulse number
        cmd = `psrstat -c nsubint $filename`
        res = read(pipeline(cmd), String)
        #f = findfirst("nsubint=", res)
        #nsubint = parse(res[f[end]+1:end])
        nsubint = parse(Int, split(split(res)[2], "=")[2])
        # generate file SNR file
        #cmd = `psrstat -l subint=0-$(nsubint-1) -qc snr $filename`
        cmd = `psrstat -l subint=0-$(nsubint-1) -c snr=pdmp -qc snr $filename` # use it in the future
        #run(pipeline(cmd, stdout="$filename.snr.txt"))
        res = split(read(pipeline(cmd), String), "\n")
        f = open("$filename.snr.txt", "w")
        snrs = []
        for i in 1:nsubint
            snr = parse(Float64, split(res[i], "=")[2])
            push!(snrs, snr)
            write(f, "$snr\n")
            #println("$i $snr")
        end
        close(f)
        println("SNR (old): ", round(sum(snrs) / sqrt(nsubint))) # why it is different?
        # psrstat -j FTp -c snr=pdmp -c snr filename
        #cmd = `psrstat -j FTp -c snr $filename`
        cmd = `psrstat -j FTp -c snr=pdmp -c snr $filename` # use it in the future
        #cmd = `psrstat -j FTp -c on=set -c on:range=350:650%bin -c off=set -c off:range=700:1000%bin -c snr=modular -c snr -c snr:on:found -c snr:off:found $filename`
        res = read(pipeline(cmd), String)
        snr = parse(Float64, split(res, "=")[2])
        println("SNR: ", round(snr))
    end


    function simps(y::Vector, x::Union{Vector,UnitRange})
    #function simps(y, x)
    # https://mmas.github.io/simpson-integration-julia
        n = length(y) - 1
        n % 2 == 0 || error("`y` length (number of intervals) must be odd")
        length(x) - 1 == n || error("`x` and `y` length must be equal")
        h = (x[end] - x[1])/n
        s = sum(view(y, 1:2:n) + 4view(y, 2:2:n) + view(y, 3:2:n+1))
        return h / 3 * s
    end

    
    """
    Creates folded profile
    :param data: single pulses
    :param p3: P_3 periodicity
    :return: folded profile
    """
    function p3fold(data, p3, ybins=10, start_ind=nothing)
        if isnothing(start_ind)
            start_ind = Int(floor(p3 / 2))
        end
    
        data = data[start_ind+1:end, :]
        folded_ = zeros(Float64, ybins, size(data, 2))
    
        dp3 = p3 / ybins
    
        for i in 1:size(data, 1)
            new_ind = i / dp3
            j = Int(floor(mod(new_ind, ybins)))
            folded_[j+1, :] += data[i, :]
        end
    
        return folded_
    end

    """
    Save parameters used in the analysis to a JSON file
    """
    function save_params(filename, params)
        open(filename, "w") do file
            JSON.print(file, params)
        end
    end


    """
    Read parameters used in the analysis from a JSON file
    """
    function read_params(filename)
        f = open(filename, "r")
        params = JSON.parse(f)
        close(f)
        return params
    end


    """
    Write default parameters to a JSON file
    """
    function default_params(filename)
        p = Dict(
            "nbin" => 1024,
            "_nbin" => "Number of pulse phase bins",
            "nfft" => 256,
            "_nfft" => "Set size of fft's [default=256]",
            "bin_st" => nothing,
            "bin_end" => nothing,
            "pulse_start" => 1,
            "pulse_end" => 150,
            "nsubint" => nothing,
            "p3" => -1.0,
            "p3_error" => -1.0,
            "p3_ybins" => nothing,
            "clean_threshold" => 0.5,
            "_clean_threshold" => "threshold for polarization cleaning"
        )
        f = open(filename, "w")
        JSON.print(f, p)
        close(f)
        return p
    end

    # TODO TODO TODO revive/use function

    # track_subpulses_snr
    # track_subpulses_snr2 <-- ten!
    function track_subpulses_snr2(data, p2, snrfile; on_st=350, on_end=650, off_st=20, off_end=320, thresh=3)
        """ based on S/N of the signal """

        f = open(snrfile)
        snrs = []
        for line in readlines(f)
            push!(snrs, parse(Float64, line))
        end

        println("obs. SNR (no?): ", round(sum(snrs) / sqrt(size(snrs)[1])))
        println("mean SNR: ", mean(snrs))
        println("median SNR: ", median(snrs))

        pulses, bins = size(data)
        p2_bins = floor(Int, p2 / 360 * bins)
        peaks = [] # [pulse_num, [p1, p2, p3...]]
        ppeaks = [] # [p1, p2, p3...]
        σ = p2_bins / 2 / 2.35482
        kernel = gauss(collect(1:p2_bins), [1, p2_bins/2, σ, 0])
        detected_pulses = 0
        located_subpulses = 0
        for i in 1:pulses
            #y = view(data, i, on_st:on_end)
            y = view(data, i, :)
            (mi, ma) = extrema(y)
            #y = (y .- mi) / (ma - mi)
            y = y ./ ma # this is much much better!

            res = conv(y, kernel)
            (mi, ma) = extrema(res)
            res = (res .- mi) / (ma - mi)
            re = res[floor(Int,p2_bins/2):end-floor(Int,p2_bins/2)]

            #on = maximum(re[on_st:on_end])
            #off = rms(y[off_st:off_end])
            #sigma = on / off
            #println("$i $sigma")
            #println("$i $(snrs[i])")
            snr = snrs[i]

            if (snr > thresh) && ~isnan(re[1])  # why re is sometimes NaN?
                peak = Tools.peaks(re)
                # new syntax?
                ma, pa = peakprom(Maxima(), re, floor(Int, p2_bins/4))
                #ma, pa = peakprom(re, Maxima(), floor(Int, p2_bins/4))
                #ma, pa = peakprom(re, Maxima(), p2_bins/2)
                inds = sortperm(pa, rev=true)
                ppeaks = []
                # get new thresh2 (maximum in off pulse region)
                #=
                new_thresh = 0.
                for ii in inds
                    if ((ma[ii] < on_st) || (ma[ii] > on_end)) && (new_thresh < re[ma[ii]])
                        new_thresh = re[ma[ii]]
                    end
                end
                thresh2 = maximum([new_thresh, thresh2_old])
                #println(thresh2)
                =#
                # TODO add singnal to noise calculation HERE
                #thresh2 = 0.7
                for ii in inds
                    st = floor(Int, ma[ii] - p2_bins / 2)
                    en = ceil(Int, ma[ii] + p2_bins / 2)
                    if (st >= on_st) && (en <= on_end)
                        #println(y)
                        yy = y[st:en]
                        signal = simps(yy, collect(st:en))
                        #signal2 = trapz(collect(st:en), yy) # should work
                        #println("$signal, $signal2")

                        # do not use integration for noise!
                        #nn = y[off_st:off_st+(en-st)]
                        #noise = simps(nn, collect(off_st:off_st+(en-st)))

                        #nn = y[off_st:off_st+(en-st)]
                        nn = y[off_st:off_end]
                        noise = std(nn) * (en-st)^0.5
                        #nn2 = y[off_st:off_st+(en-st)]
                        #noise2 = std(nn2) * (en-st)^0.5
                        #println("$noise $noise2")

                        #println(signal / noise)

                        if signal / noise > thresh
                            push!(ppeaks, ma[ii])
                        end

                    end
                    #if (re[ma[ii]] >= thresh2) &&  (ma[ii] > on_st) && (ma[ii] < on_end)
                    #    push!(ppeaks, ma[ii])
                    #end
                end
                if length(ppeaks) > 0
                    located_subpulses += 1
                end
                push!(peaks, [i, ppeaks])
                detected_pulses += 1
                #println("$i $ppeaks")
                #=
                PyPlot.close()
                plot(y, c="black")
                plot(re, c="red")
                plot(kernel, c="blue")
                #for ii in inds
                for p in ppeaks
                    axvline(x=p-1, c="pink")
                end
                #axvline(x=ma[inds[1]]-1, c="green") # TODO plotting starts from 0!
                #axvline(x=ma[inds[2]]-1, c="magenta") # TODO plotting starts from 0!
                #axvline(x=ma[inds[3]]-1, c="brown") # TODO plotting starts from 0!
                show()
                st = readline(stdin; keep=false)
                if st == "q"
                    break
                end
                =#
            end
        end
        println("Number of pulses: $pulses  Included pulses : $detected_pulses  Pulses with located subpulses: $located_subpulses")
        frac = located_subpulses / pulses * 100
        println("Fraction of pulses used: $(round(frac)) (before grouping (overestimated!))")
        return peaks
    end


    """ 
    Track subpulses based on S/N of the signal 
    new version 
    """
    function track_subpulses_snr3(data, data2, p2, snrfile; on_st=350, on_end=650, off_st=20, off_end=320, thresh=3)

        f = open(snrfile)
        snrs = []
        for line in readlines(f)
            push!(snrs, parse(Float64, line))
        end

        println("obs. SNR (no?): ", round(sum(snrs) / sqrt(size(snrs)[1])))
        println("mean SNR: ", mean(snrs))
        println("median SNR: ", median(snrs))

        data = integrate(data, 8)
        data2 = integrate(data, 8)

        pulses, bins = size(data)

        p2_bins = floor(Int, p2 / 360 * bins)
        if p2_bins % 2 == 0
            p2_bins += 1  # we force an odd length
        end
        σ = p2_bins / 2 / 2.35482 # why like that?
        kernel = gauss(collect(1:p2_bins), [1, p2_bins/2, σ, 0])

        pulse_st = 1
        pulse_end = 10

        pu_1 = []
        loc_1 = []
        pu_2 = []
        loc_2 = []

        for i in pulse_st:pulse_end # pulses
            y = view(data, i, on_st:on_end) # single pulse

            # normalize the data
            (mi, ma) = extrema(y)
            #y = (y .- mi) / (ma - mi)
            #y = y ./ ma # this is much much better!

            # convolution with gauss function
            res = conv(y, kernel)
            (mi, ma) = extrema(res)
            # full convolution result
            res = (res .- mi) / (ma - mi) 
            # removes boundry artifacts => re and y should have the same sizes
            half_kernel = div(p2_bins, 2)
            re = res[half_kernel+1:end-half_kernel] # this is why ind-1?

            #sg = savitzky_golay(y, 21, 3)

            # work with data2
            y2 = view(data2, i, on_st:on_end) # single pulse
            # normalize the data
            (mi, ma) = extrema(y2)
            #y2 = y2 ./ ma
            res2 = conv(y2, kernel)
            (mi, ma) = extrema(res2)
            # full convolution result
            res2 = (res2 .- mi) / (ma - mi) 
            # removes boundry artifacts => re and y should have the same sizes
            half_kernel = div(p2_bins, 2)
            re2 = res2[half_kernel+1:end-half_kernel] # this is why ind-1?
            #sg2 = savitzky_golay(y2, 21, 3)

            snr = snrs[i]

            # just maximum finding
            val, idx = findmax(y)
            println("Max. val=$val i=$idx pulse=$i")
            push!(pu_1, i)
            push!(loc_1, idx)
            val2, idx2 = findmax(y2)
            push!(pu_2, i)
            push!(loc_2, idx2)

            #=
            if snr > thresh

                pks = findmaxima(re)
                pks = peakproms(pks)
                #pks = peakwidths(pks)

                pks2 = findmaxima(re2)
                pks2 = peakproms(pks2)
                #pks2 = peakwidths(pks2)

                #println(pks)
                println(i)
                println(keys(pks))
                
                #println(pks[:proms])
                #println(pks[:widths])
                

                PyPlot.close()

                for (j,ind) in enumerate(pks[:indices])
                    if pks[:proms][j] > 0.5
                        scatter([ind-1],[pks[:heights][j]]) # why ind-1?
                        plot([ind-1, ind-1],[pks[:heights][j], pks[:heights][j]-pks[:proms][j]], c="red")
                        push!(pu_1,i)
                        push!(loc_1, ind-1)
                    end
                end

                for (j,ind) in enumerate(pks2[:indices])
                    if pks2[:proms][j] > 0.5
                        scatter([ind-1],[pks2[:heights][j]]) # why ind-1?
                        plot([ind-1, ind-1],[pks2[:heights][j], pks2[:heights][j]-pks2[:proms][j]], c="C2")
                        push!(pu_2, i)
                        push!(loc_2, ind-1)
                    end
                end


                plot(y, c="black", lw=1)
                #plot(sg.y, c="black", lw=3, alpha=0.5)
                plot(re, c="red")
                #plot(sg.y, c="blue")


                plot(y2, c="C1", lw=1)
                #plot(sg2.y, c="C1", lw=3, alpha=0.5)
                plot(re2, c="C2")
                
                #xlim(on_st, on_end) # if full
                
                show()
                st = readline(stdin; keep=false)
                #st = 1
                if st == "q"
                    break
                end

            end
            =#

        end

        PyPlot.close()
        figure()

        subplot2grid((1, 2), (0, 0))
        imshow(data[:,on_st:on_end], origin="lower", cmap="viridis", interpolation="none", aspect="auto", extent=[1, on_end-on_st + 1+1, 1, pulses+1])
        scatter(loc_1 .+0.5, pu_1.+0.5, marker="x", color="C1", s=50)


        subplot2grid((1, 2), (0, 1))
        imshow(data2[:,on_st:on_end], origin="lower", cmap="viridis", interpolation="none", aspect="auto", extent=[1, on_end-on_st + 1+1, 1, pulses+1])
        scatter(loc_2 .+ 0.5, pu_2 .+0.5, marker="x", color="C2", s=50)
        scatter([69.5], [248.5], marker="x", color="C3", s=50)
        println()

        PyPlot.show()
        st = readline(stdin; keep=false)


    end


    """
    function integrates / adds (pulses_num) single pulses
    """
    function integrate(data, pulses_num)
        pulses, bins = size(data)
        groups = div(pulses, pulses_num)
        result = zeros(groups, bins)
    
        for g in 1:groups
            start_idx = (g - 1) * pulses_num + 1
            end_idx = g * pulses_num
            result[g, :] = sum(data[start_idx:end_idx, :], dims=1)
        end
        
        return result
    end
    

end  # module Tools
