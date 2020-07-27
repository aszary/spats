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

    function segmented(x, y, npsi; preview=false)
        # splines fitting to estimate number of extrema
        sl = fit(SmoothingSpline, x, y, 1000.0)
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
        try
            R"line <- broken.line(segmod)"
        catch
            println("BBBB")
            return x, y, ysl
        end
        @rget line
        #show(R"summary(segmod)")
        #show(R"print(segmod)")
        #show(R"slope(segmod)")
        #show(R"plot(segmod)")
        #show(R"points(segmod)")
        #show(R"print(segmod)")
        #show(R"lines.segmented(segmod)")
        y2 = line[:fit]
        if preview
            plot(x, y2, lw=3)
            scatter(x, y)
            plot(x, ysl)
            show()
            readline(stdin; keep=false)
            close()
        end

        return x, y2, ysl


    end


end
