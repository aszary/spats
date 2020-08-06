module PyRModule
    using PyCall
    using RCall

    using Peaks
    using SmoothingSplines
    using PyPlot

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

    function segmented(x, y, npsi; lambda=1000.0, preview=false)
        # splines fitting to estimate number of extrema
        sl = fit(SmoothingSpline, x, y, lambda)
        ysl = SmoothingSplines.predict(sl)

        if preview
            scatter(x, y)
            plot(x, ysl)
            savefig("test.pdf")
            close()
        end

        # broken-lines fitting
        @rput x
        @rput y
        @rput npsi
        R"linmod <- lm(y ~ x)"
        #R"summary(linmod)"
        R"library(segmented)"
        R"segmod <- segmented(linmod, seg.Z = ~x, npsi=npsi)"
        #=
        if npsi == 3
            R"segmod <- segmented(linmod, seg.Z = ~x, psi=c(520, 542, 565))"
        else
            R"segmod <- segmented(linmod, seg.Z = ~x, npsi=npsi)"
        end
        =#
        # lines
        try
            R"line <- broken.line(segmod)"
        catch
            println("BBBB")
            return x, y, ysl, nothing, nothing, nothing, nothing, nothing, nothing
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
        for i in 1:npsi-1
            push!(xl, bp[i] + (bp[i+1] - bp[i]) / 2) # first point
            push!(exl, (bp[i+1] -bp[i]) / 2) # first point
        end
        push!(xl, bp[npsi] + (x[end] - bp[npsi]) / 2) # last point
        push!(exl, (x[end] - bp[npsi]) / 2) # last point
        #println("xl ", xl)
        #println("slopes ", slopes)
        return x, y2, ysl, slopes, eslopes, bp, ebp, xl, exl
    end


end
