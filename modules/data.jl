module Data
    using Glob
    using FITSIO
    using ProgressMeter

    include("functions.jl")

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
    Process data with PSRCHIVE and PSRSALSA
    """
    function process_psrdata(indir, outdir; outfile="pulsar.spCF", files=nothing)
        # Ensure output directory exists
        mkpath(outdir)
    
        if files === nothing
            # Find all .spCF files in the input directory
            files = filter(f -> endswith(f, ".spCF"), readdir(indir))
        end
    
        # Sort files based on pulse numbers (e.g., 00000-00255, 00256-00511, etc.)
        sort!(files, by = f -> begin
            m = match(r"_(\d+)-(\d+)\.spCF$", f)
            if isnothing(m)
                return typemax(Int)  # Files without proper format go to the end
            else
                return parse(Int, m.captures[1])  # Sort by the starting pulse number
            end
        end)
    
        if isempty(files)
            error("No .spCF files found in directory: $indir")
        end
    
        println("Processing files in order:")
        for (i, f) in enumerate(files)
            println("$i. $f")
        end
    
        file_names = [joinpath(indir, file) for file in files]
    
        # Extract pulsar name dynamically from the directory structure
        pulsar_name_match = match(r"J\d{4}\+\d{4}", indir)
        if isnothing(pulsar_name_match)
            error("Could not extract pulsar name from directory path: $indir")
        end
        pulsar_name = pulsar_name_match.match
    
        # Create pulsar-specific output directory
        pulsar_outdir = joinpath(outdir, pulsar_name)
        mkpath(pulsar_outdir)
        println("Created new directory: $pulsar_outdir")
    
        # Define file paths inside the pulsar-specific directory
        outfile = joinpath(pulsar_outdir, "$pulsar_name.spCF")
        debased_file = joinpath(pulsar_outdir, "$pulsar_name.debase.gg")
    
        # Combine all files using psradd
        run(pipeline(`psradd $file_names -o $outfile`, stderr="errs.txt"))
    
        # Debase the data
        run(pipeline(`pmod -device "/xw" -debase $outfile`, `tee pmod_output.txt`))
    
        # Read captured output
        output = read("pmod_output.txt", String)
        rm("pmod_output.txt")  # cleanup
    
        # Extract onpulse values
        m = match(r"-onpulse '(\d+) (\d+)'", output)
        if isnothing(m)
            error("Failed to extract onpulse range from pmod output")
        end
        bin_st, bin_end = parse.(Int, m.captures)
    
        println("Found onpulse range: $bin_st to $bin_end")
    
        # Move the debased file into the pulsar directory (force overwrite if needed)
        mv(outfile, debased_file; force=true)
    
        # Calculate 2dfs and lrfs
        run(pipeline(`pspec -w -2dfs -lrfs -profd "/NULL" -onpulsed "/NULL" -2dfsd "/NULL" -lrfsd "/NULL" -nfft 256 -onpulse "$bin_st $bin_end" $debased_file`, stderr="errs.txt"))
    
        # Find P3
        run(pipeline(`pspecDetect -v -device "/xw" $debased_file`, `tee pspecDetect_output.txt`))
        output = read("pspecDetect_output.txt", String)
        rm("pspecDetect_output.txt")  # cleanup
    
        # Extract P3 value from the last occurrence
        p3_matches = collect(eachmatch(r"P3\[P0\]\s*=\s*(\d+\.\d+)\s*\+-\s*(\d+\.\d+)", output))
        if isempty(p3_matches)
            error("Failed to extract P3 value from pspecDetect output")
        end
        last_match = p3_matches[end]
        p3_value = parse(Float64, last_match.captures[1])
        p3_error = parse(Float64, last_match.captures[2])
        println("Found P3 = $p3_value ± $p3_error P0")
    
        ybins = Functions.find_ybins(p3_value)
        println("Number of ybins: $ybins")
    
        # Compute pfold and save results in the pulsar-specific directory
        pfold_file = joinpath(pulsar_outdir, "$pulsar_name.debase.p3fold")
        println("Running pfold command...")
        run(pipeline(`pfold -p3fold "$p3_value $ybins" -onpulse "$bin_st $bin_end" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $debased_file -o $pfold_file`, stderr="errs.txt"))
        println("pfold command finished")
    
        # Check if the pfold file was created
        if !isfile(pfold_file)
            error("ERROR: pfold output file was not created: $pfold_file")
        end
    
        println("All files are stored in: $pulsar_outdir")
        return bin_st - 20, bin_end + 20
    end
    
end # module
