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
    Count total pulses in an observation directory by parsing filenames.
    Files are expected to match pattern: ..._NNNNN-MMMMM.spCf16 or .spCF
    Returns the last pulse index + 1, or nothing if no files match.
    """
    function count_pulses(indir)
        files = filter(f -> endswith(f, ".spCf16") || endswith(f, ".spCF"), readdir(indir))
        max_pulse = -1
        for f in files
            m = match(r"_(\d+)-(\d+)\.(spCf16|spCF)$", f)
            if !isnothing(m)
                max_pulse = max(max_pulse, parse(Int, m.captures[2]))
            end
        end
        return max_pulse >= 0 ? max_pulse + 1 : nothing
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
