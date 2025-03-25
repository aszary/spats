module Data
    using Glob
    using FITSIO
    using ProgressMeter

    """
    Zero selected pulses

    """
    function zap!(data; ranges=nothing)
        pulses, bins = size(data)
        zaps = []
        for rn in ranges
            if length(rn) == 1
                push!(zaps, rn)
            else
                zaps = vcat(zaps, collect(rn[1]:rn[2]))
            end
            for z in zaps
                for i in 1:bins
                    data[z,i] = 0
                end
            end
        end
    end


    """
    Loads PSRCHIVE ASCII file
    """
    function load_ascii(infile)
        f = open(infile)
        lines = readlines(f)
        res = split(lines[1])
        pulses = parse(Int, res[6])
        bins = parse(Int, res[12])
        data = Array{Float64}(undef, pulses, bins)
        for i in 2: length(lines)
            res = split(lines[i])
            bin = parse(Int, res[3]) + 1
            pulse = parse(Int, res[1]) + 1
            inte = parse(Float64, res[4])
            data[pulse, bin] = inte
        end
        close(f)
        return data
    end



    """
    Save PSRCHIVE ASCII file
    """
    function save_ascii(data, outfile)
        f = open(outfile, "w")
        pulse_num = size(data)[1]
        bin_num = size(data)[2]
        header = "File: filename Src: J1750-3503 Nsub: $pulse_num Nch: 1 Npol: 1 Nbin: $bin_num RMS: 0.0\n"
        write(f, header)

        for i in 1:pulse_num
            for j in 1:bin_num
                write(f, "$(i-1) 0 $(j-1) $(data[i,j])\n")
            end
        end
        close(f)
        return
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
    function convert_psrfit_ascii(infile, outfile)
        run(pipeline(`pdv -t -F -p $infile`, stdout="$outfile", stderr="errs.txt"))
        # change -t to -A to get frequancy information
        #@showprogress 1 for i in 1:pn  # psrchive indexing
        #end
    end

    """
    Process data with PSRCHIVE
    """
    function process_psrchive(indir, outdir, files; outfile="pulsar.spCF")
        file_names = [joinpath(indir, file) for file in files]
        outfile = joinpath(outdir, outfile)
        # connecting all files
        run(pipeline(`psradd $file_names -o $outfile`, stderr="errs.txt")) # PSRCHIVE
        # debase the data
        #run(pipeline(`pmod -debase $outfile`,  stderr="errs.txt")) # PSRSALSA
        buffer = IOBuffer()
        run(pipeline(`pmod -debase $outfile`, stdout=buffer, stderr=buffer))
        seekstart(buffer)
        output = String(read(buffer))
        println(output)
        return
        debased_file = replace(outfile, ".spCF" => ".debase.gg")
        run(pipeline(`pspec -w -2dfs -lrfs -nfft 256 $debased_file`,  stderr="errs.txt"))
        run(pipeline(`pspecDetect -v $debased_file`,  stderr="errs.txt"))
        # TODO read P3 from pspecDetect output  
        #run(pipeline(`pfold -p3fold "43.5 87" -p3fold_nritt 50 -p3fold_cpb 50 -w -oformat ascii $debased_file`,  stderr="errs.txt"))
        run(pipeline(`pfold -p3fold "41 82" -w -oformat ascii $debased_file`,  stderr="errs.txt"))

    end

end # module
