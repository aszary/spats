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
        p3lo = max(1, floor(Int, p3))
        p3hi = max(1, ceil(Int, p3))

        target = if p3 <= 5
            ceil(Int, 2p3)
        elseif p3 <= 15
            ceil(Int, p3)
        else
            max(ceil(Int, p3 / 2), 5)
        end

        if !isnothing(n_pulses)
            target = min(target, max(1, n_pulses ÷ min_ppb))
        end
        target = min(target, max_ybins)

        # candidates from divisors and multiples of both floor and ceil of p3
        candidates = Set{Int}()
        for p3r in (p3lo, p3hi)
            # multiples: only small k needed, stop when exceeds target
            k = 1
            while k * p3r <= target
                push!(candidates, k * p3r)
                k += 1
            end
            # divisors: must iterate up to p3r itself (e.g. 85÷17=5 needs k=17)
            for k in 1:p3r
                d, r = divrem(p3r, k)
                r == 0 && 1 <= d <= target && push!(candidates, d)
            end
        end

        push!(candidates, target)
        return p3 > 30 ? max(5, maximum(candidates)) : maximum(candidates)

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
