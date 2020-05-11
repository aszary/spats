module Tools
    using FFTW
    using Peaks
    using LsqFit
    using Statistics
    using SmoothingSplines

    using PyPlot
    using DSP
    using LinearAlgebra


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
        if a == nothing
            a = maximum(ydata)
        end
        if μ == nothing
            ind = trunc(Int, length(xdata) / 2)
            μ = xdata[ind]
        end
        if σ == nothing
            σ = xdata[20] - xdata[1]
        end
        if baselevel == nothing
            baselevel = mean(ydata)
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


    function p2_estimate(data; on_st=450, on_end=700, off_st=100, off_end=350, thresh=3.3, win=6)
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
            p0 = [1.0, 0.5, 530.0, 600.0, 15, 15]
            #p0 = [1.0, 0.5, 430.0, 500.0, 15, 15]
            xdata = collect(on_st:on_end)
            pa, errs = Tools.fit_twogaussians(xdata, av[on_st:on_end], p0[1], p0[2], p0[3], p0[4], p0[5], p0[6])
            dbin = abs(pa[2] - pa[5])
            p2 = dbin / bins * 360
            # nesty hack!
            if (p2 > 16) && p2 < (20)
                push!(p2s, p2)
            end
            #=
            ga = twogauss(xdata, pa)
            PyPlot.close()
            plot(average, c="black", lw=2)
            plot(av, lw=0.3)
            plot(xdata, ga, lw=0.6, c="red")
            savefig("output/test.pdf")
            readline(stdin; keep=false)
            #show()
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
        PyPlot.close()
        plot(average, c="black", lw=2)
        for av in avs
            plot(av, lw=0.3)
        end
        savefig("output/test.pdf")
        #readline(stdin; keep=false)
        =#
        return p2
    end


    function track_subpulses(data, p2; on_st=450, on_end=700, off_st=100, off_end=350, thresh=2.1, thresh2=0.8)
        pulses, bins = size(data)
        p2_bins = floor(Int, p2 / 360 * bins)
        peaks = [] # [pulse_num, [p1, p2, p3...]]
        ppeaks = [] # [p1, p2, p3...]
        for i in 1:pulses
            σ = p2_bins / 2 / 2.35482
            kernel = gauss(collect(1:p2_bins), [1, p2_bins/2, σ, 0])
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
                    if (re[ma[ii]] > thresh2) &&  (ma[ii] > on_st) && (ma[ii] < on_end)
                        push!(ppeaks, ma[ii])
                    end
                end
                push!(peaks, [i, ppeaks])
                #=
                println("$i $ppeaks")
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

end  # module Tools
