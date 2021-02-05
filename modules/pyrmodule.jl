module PyRModule
    using PyCall
    using RCall

    using Peaks
    using SmoothingSplines
    using PyPlot
    using Statistics
    using DataFrames, GLM


    function __init__()
        py"""
        import numpy as np
        from scipy.optimize import leastsq

        def fun(v, x, num):
            if x < v[0]:
                print("x:{} < v0:{}".format(x, v[0]))
                return v[num] + v[num+1] * x
            elif x > v[num-1]:
                print("x:{} > v{}:{}".format(x, num-1, v[num-1]))
                return v[num+(num-1)*2] + v[num+(num-1)*2+1] * x
            for i in range(0, num-1):
                if (x > v[i]) and (x < v[i+1]):
                    print("v{}:{} < x:{} < v{}:{}".format(i, v[i], x, i+1, v[i+1]))
                    return v[num+i*2] + v[num+i*2+1] * x

        def fun_vec(v, x, num):
            y = np.zeros(x.shape)
            for i in range(len(y)):
                y[i] = fun(v, x[i], num)
            return y

        def least_sq(x, y, fun, v0, num=4, show_=False):
            ## Error function
            errfunc = lambda v, x, y, num: (fun(v, x, num) - y)
            res = leastsq(errfunc, v0, args=(np.array(x), np.array(y), num), maxfev=10000, full_output=True)
            v, conv = res[0], res[1]
            chi_sq_red2 = (res[2]['fvec'] ** 2.).sum() / (len(y)-len(v))
            res_errs = conv * chi_sq_red2
            errs = (np.absolute(res_errs[0][0])**0.5, np.absolute(res_errs[1][1])**0.5)
            if show_ is True:
                print("convolution (jacobian around the solution)", conv)
                print("chi square:", sum(pow(errfunc(v, np.array(x), np.array(y), num), 2.)))
                print("chi^2_red = ", chi_sq_red2)
                print("Parameters:", v)
                print("Errors:", errs)
            return v, errs
        """
    end


    function least_sq(x, y, v0; num=1, show_=false)
        return py"least_sq"(x, y, py"fun_vec", v0, num=num, show_=show_)
    end


    function segmented(x, y, npsi; lambda=1000.0, fixedpsi=nothing, preview=false)
        @. f(x, p) = p[2] * x + p[1]

        # splines fitting to estimate number of extrema
        sl = fit(SmoothingSpline, x, y, lambda)
        ysl = SmoothingSplines.predict(sl)

        # broken-lines fitting
        @rput x
        @rput y
        @rput npsi
        R"linmod <- lm(y ~ x)"
        #R"summary(linmod)"
        R"library(segmented)"
        #=
        if npsi == 3
            R"segmod <- segmented(linmod, seg.Z = ~x, psi=c(520, 542, 565))"
        else
            R"segmod <- segmented(linmod, seg.Z = ~x, npsi=npsi)"
        end
        =#
        # lines
        try
            if fixedpsi != nothing
                R"segmod <- segmented(linmod, seg.Z = ~x, npsi=npsi, fixed.psi=$fixedpsi)"
                #R"segmod <- segmented(linmod, seg.Z = ~x, npsi=npsi)"
            else
                R"segmod <- segmented(linmod, seg.Z = ~x, npsi=npsi)"
            end
            R"line <- broken.line(segmod)"
        catch
            R"con <- confint(linmod, level=0.68)"
            @rget con
            bp = []
            ebp = []
            for i in 1:size(con)[1]
                push!(bp, con[i, 1])
                push!(ebp, con[i, 1] - con[i, 2]) # it is fine
                #println(con[i, 2] - con[i, 2], " ", con[i, 3] - con[i, 1])
            end
            #R"lnn <- line(x, y)"
            #R"y2 <- fitted(lnn)"
            #R"co <- coef(lnn)"
            data = DataFrame(X=x, Y=y)
            ols = lm(@formula(Y ~ X), data)
            y2 = f(x, coef(ols))
            slopes = [coef(ols)[2]]
            eslopes = [stderror(ols)[2]]
            # get lines means (for slopes plotting)
            xl = []
            exl = []
            push!(xl, x[1] + (x[end] - x[1]) / 2) # first point # is it enougth?
            push!(exl, (x[end] - x[1]) / 2) # first point
            #push!(inclines, [mean(x_), coef(ols)[2], stderror(ols)[2]])
            println("BBBB")
            if preview
                plot(x, y2, lw=3)
                scatter(x, y)
                plot(x, ysl)
                show()
                readline(stdin; keep=false)
                close()
            end
            return x, y2, ysl, slopes, eslopes, bp, ebp, xl, exl
        end
        @rget line
        y2 = line[:fit]
        # slopes
        R"slo <- slope(segmod)"
        @rget slo
        slopes = []
        eslopes = []
        for i in 1:size(slo[:x])[1]
            push!(slopes, slo[:x][i,1])
            push!(eslopes, slo[:x][i,2])
        end
        # break points one sigma
        R"con <- confint(segmod, level=0.68)"
        @rget con
        bp = []
        ebp = []
        for i in 1:size(con)[1]
            push!(bp, con[i, 1])
            push!(ebp, con[i, 1] - con[i, 2]) # it is fine
            #println(con[i, 2] - con[i, 2], " ", con[i, 3] - con[i, 1])
        end
        #show(R"summary(segmod)")
        #show(R"print(segmod)")
        #show(R"plot(segmod)")
        #show(R"points(segmod)")
        #show(R"print(segmod)")
        #show(R"lines.segmented(segmod)")
        #show(R"confint(segmod, level=0.95)")
        if preview
            plot(x, y2, lw=3)
            scatter(x, y)
            plot(x, ysl)
            show()
            readline(stdin; keep=false)
            close()
        end

        # get lines means (for slopes plotting)
        xl = []
        exl = []
        push!(xl, x[1] + (bp[1] - x[1]) / 2) # first point
        push!(exl, (bp[1] - x[1]) / 2) # first point
        for i in 1:npsi - 1
            push!(xl, bp[i] + (bp[i+1] - bp[i]) / 2) # first point
            push!(exl, (bp[i+1] -bp[i]) / 2) # first point
        end
        push!(xl, bp[npsi] + (x[end] - bp[npsi]) / 2) # last point
        push!(exl, (x[end] - bp[npsi]) / 2) # last point
        #println("xl ", xl)
        #println("slopes ", slopes)
        return x, y2, ysl, slopes, eslopes, bp, ebp, xl, exl
    end


    """ overfits """
    function cubic_spline(x, y)
        interpolate =  pyimport("scipy.interpolate")
        #cs = interpolate.CubicSpline(x, y)
        spl = interpolate.UnivariateSpline(x, y)
        #spl.set_smoothing_factor(5000.0)  # the idea was not to use any smoothing factor
        line = [x, spl(x)]
        ysp = spl(x)
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


    function smooth_spline(x, y, spar; lambda=0.001, df_start=nothing)
        #println(x)
        #println(y)
        #df_start = length(x) - 50
        @rput x
        @rput y
        #@rput lambda
        #R"spl <- smooth.spline(x, y, lambda=$lambda)"
        R"spl <- smooth.spline(x, y, spar=$spar)"
        #R"spl <- smooth.spline(x, y, cv=FALSE)" #GCV
        #R"spl <- smooth.spline(x, y, df=$df_start)"
        R"yspp <- predict(spl, x)"
        @rget spl
        @rget yspp
        ysp = yspp[:y]
        dof = spl[:df]
        #println()
        #println(spl[:lambda])
        #println(ysp)
        #println("dof ", dof)
        #=
        rep = 0
        r2 = rsquared(y, ysp)
        while (r2 < 0.9) || (r2 > 0.96)
            rep += 1
            if r2 < 0.9
                dof += 1
            else
                dof -= 1
            end
            R"spl <- smooth.spline(x, y, df=$dof)"
            R"yspp <- predict(spl, x)"
            @rget spl
            @rget yspp
            ysp = yspp[:y]
            dof = spl[:df]
            r2 = rsquared(y, ysp)
            chisq_red = chisquare_reduced(y, ysp, dof)
            println("$rep Chi^2_red $dof ", chisq_red)
            println("$rep R^2 ", r2)
        end
        =#
        #=
        chisq_red = chisquare_reduced(y, ysp, dof)
        while (chisq_red > 1.2)||(chisq_red < 0.8)
            rep += 1
            if chisq_red > 1.1
                dof += 1
            else
                dof -= 1
            end
            R"spl <- smooth.spline(x, y, df=$dof)"
            R"yspp <- predict(spl, x)"
            @rget spl
            @rget yspp
            ysp = yspp[:y]
            dof = spl[:df]
            chisq_red = chisquare_reduced(y, ysp, dof)
            println("$rep Chi^2_red $dof", chisq_red)
        end
        =#
        return ysp, dof
    end


    """ lazy  """
    function chisquare_reduced(o, c, dof)
        return sum((o .- c) .^ 2 ./ var(o)) / dof
    end


    """ https://pythonprogramming.net/how-to-program-r-squared-machine-learning-tutorial/ """
    function rsquared(o, c)
        mean_line = [mean(o) for y in o]
        squared_error_fit = squared_error(o, c)
        squared_error_mean = squared_error(o, mean_line)
        return 1 - (squared_error_fit / squared_error_mean)
    end


    function squared_error(o, c)
        return sum((c .- o) .* (c .- o))
    end


end
