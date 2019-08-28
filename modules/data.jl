module Data
    using Glob
    using FITSIO
    using ProgressMeter


    """
    Loads PSRCHIVE ASCII files
    """
    function load_ascii(regex; bins=1024)
        files = glob(regex)
        data = Array{Float64}(undef, length(files), bins)
        for file in files
            # get pulse number # wirdo
            fl = split(file, "/")[end]
            pulse = parse(Int, replace(replace(fl, "pulse_" => ""), ".txt" => ""))
            f = open(file)
            lines = readlines(f)
            for line in lines
                if ~startswith(line, "#")
                    res = split(line)
                    bin = parse(Int, res[3]) + 1
                    inte = parse(Float64, res[4])
                    data[pulse, bin] = inte
                end
            end
            close(f)
        end
        return data
    end

    """
    Loads data in FITS format
    """
    function load_fits(filename)
        f = FITS(filename)
        for hdu in f; println(typeof(hdu)); end
        he = read_header(f[5])
        println(he)
        close(f)
    end

    """
    Converts PSRFIT file to ASCII (using PSRCHIVE tools)
    """
    function convert_psrfit_ascii(infile, outdir)
        # get number of pulses
        pn = parse(Int, split(read(`psrstat -q -c nsubint $infile`, String), "=")[2])

        @showprogress 1 for i in 1:pn  # psrchive indexing
            run(pipeline(`pdv -A -F -p -K -i $(i-1) $infile`, stdout="$outdir/pulse_$i.txt", stderr="$outdir/errs.txt"))
        end
    end

end # module
