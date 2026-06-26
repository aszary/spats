module Functions


    function p3_obs(p3, n) # P3 in Periods
        p3_obs = 1 / (1 / p3 - n)
        return p3_obs
    end


    function p3(p3_obs, n) # P3_obs in Periods
        p3 = 1 / (n + 1 / (p3_obs))
        return p3
    end


    function p3_edwards(p3_obs, p, n) # P3_obs in Periods
        p3_sec = p / (n + p / (p3_obs*p))
        p3 = p3_sec / p
        return p3
    end


    function p3_gupta(p3_obs, p, n) # P3_obs in Periods
        fobs = 1 / (p3_obs * p)
        fN = 0.5 / p
        k = trunc(Int64, (n+1) / 2)
        l = mod(n, 2)
        f = 2 * k * fN + (-1)^l * fobs
        # in seconds
        p3_sec = 1 / f
        p3 = p3_sec / p
        return p3
    end


    """
    Find the best number of bins for a given P3 value
    """
    function find_ybins_old(p3)
        if p3 > 8
            return ceil(Int, p3)
        else
            return ceil(Int, 2 * p3)
        end
    end


    function find_ybins(p3, n_pulses=nothing; min_ppb=30, max_ybins=10)
        # Step 1: find k (1..20) that makes k*p3 closest to an integer
        # (best rational approximation of p3 as n/k)
        best_k = 1
        best_rel_err = abs(p3 - round(p3)) / p3
        for k in 2:20
            rel_err = abs(k * p3 - round(k * p3)) / (k * p3)
            if rel_err < best_rel_err
                best_rel_err = rel_err
                best_k = k
            end
        end
        base = round(Int, best_k * p3)

        # Step 2: apply SNR cap
        ybins_max = max_ybins
        if !isnothing(n_pulses)
            ybins_max = min(ybins_max, max(1, n_pulses ÷ min_ppb))
        end

        # Step 3: largest divisor of base that fits within ybins_max
        result = 1
        for k in 1:base
            d, r = divrem(base, k)
            r == 0 && 1 <= d <= ybins_max && (result = max(result, d))
        end

        return p3 > 30 ? max(5, result) : result
    end



    """
        normalize(data)

    Normalize numeric data to range [0,1] using min-max normalization.
    Modifies input array in-place.

    # Arguments
    - `data`: Numeric array to normalize
    """
    function normalize(data)
        data .-= minimum(data)
        data ./= maximum(data)
    end

end  # module Functions
