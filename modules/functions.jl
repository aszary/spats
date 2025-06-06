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
    function find_ybins(p3)
        max_bins = round(Int, 5 * p3)

        best_diff = Inf  # Initialize the best difference with a large number
        best_ybins = 0   # Initialize the best number of bins
        
        # Start checking from p3
        start_ybins = round(Int, p3)
        
        # Loop through possible values of ybins from start_ybins to max_bins
        for ybins in start_ybins:max_bins
            # Calculate the ratio ybins / p3
            ratio = ybins / p3
            
            # Find the nearest integer to the ratio
            nearest_int = round(ratio)
            
            # Calculate the absolute difference between the ratio and the nearest integer
            diff = abs(ratio - nearest_int)
            
            # If this difference is smaller than the current best difference, update best_ybins
            if diff < best_diff
                best_diff = diff
                best_ybins = ybins
            end
        end
        return best_ybins
    end

    
    function fftshift(x)
        n = size(x, 1)
        m = size(x, 2)
        return circshift(x, (div(n,2), div(m,2)))
    end


end  # module Functions
