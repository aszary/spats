module Tools
    using Glob
    using JLD2
    using FFTW
    using Peaks
    using LsqFit
    using Statistics
    using SmoothingSplines

    using PyPlot
    using DSP
    using LinearAlgebra
    using DataFrames, GLM


    function rms(data)
        s = 0.0
        for a in data
            s += a * a
        end
        return sqrt(s / length(data))
    end


    function average_profile(data)
        pulses, bins = size(data)
        ave = zeros(bins)
        for i in 1:pulses
            for j in 1:bins
                ave[j] += data[i,j]
            end
        end
        ma = maximum(ave)
        ave = ave / ma
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


    function peaks(intensity)
        #max = maxima(intensity, 30)
        ma, pa = peakprom(intensity, Maxima(), 1) #floor(Int, length(intensity)/20))  # works for J1750
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


    function track_subpulses(data, p2; on_st=450, on_end=700, off_st=100, off_end=350, thresh=2.1, thresh2=0.8)
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
                ma, pa = peakprom(re, Maxima(), floor(Int, p2_bins/4))
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

        return tracks, lines, inclines, dr1, ysp
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


end  # module Tools
