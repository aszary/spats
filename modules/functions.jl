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


    function find_ybins(p3, n_pulses=nothing; min_ppb=50)
        ybins = floor(Int, 2 * p3)
        if !isnothing(n_pulses)
            ybins = min(ybins, floor(Int, n_pulses / min_ppb))
        end
        return max(4, ybins)
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
