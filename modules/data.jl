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
    
        # Sort files based on the pulse range extracted from the filename (e.g., "00000-00255")
        sort!(files, by = f -> begin
            # Extract the pulse range from the filename (e.g., "2020-02-02-11:45:29_00000-00255.spCF")
            m = match(r"_(\d+-\d+)\.spCF$", f)
            if isnothing(m)
                return typemax(Int)  # Files without the expected format go to the end
            else
                return parse(Int, split(m.captures[1], "-")[1])  # Sort by the starting pulse number
            end
        end)
    
        # Check if there are any files to process
        if isempty(files)
            error("No .spCF files found in directory: $indir")
        end
    
        println("Processing files in order:")
        for (i, f) in enumerate(files)
            println("$i. $f")
        end
    
        file_names = [joinpath(indir, file) for file in files]
    
        # Extract pulsar name from the first file (pulsar name is the part before the pulse range)
        first_file = files[1]
        pulsar_match = match(r"_(\d+-\d+)\.spCF$", first_file)  # Extract the pulse range from the first file
        if pulsar_match === nothing
            error("Could not determine pulsar name from filename: $first_file")
        end
    
        pulsar_name = first_file |> x -> split(x, "_")[1]  # Get the pulsar name, which is the part before the "_"
        pulsar_outdir = joinpath(outdir, pulsar_name)  # Create a subdirectory for the pulsar
        mkpath(pulsar_outdir)  # Ensure pulsar-specific directory exists
    
        outfile_path = joinpath(pulsar_outdir, outfile)  # Define the output file path
    
        # Connect all files using PSRCHIVE
        run(pipeline(`psradd $file_names -o $outfile_path`, stderr="errs.txt"))
    
        # Debase the data using pmod tool
        run(pipeline(`pmod -device "/xw" -debase $outfile_path`, `tee pmod_output.txt`))
        
        # Read captured output
        output = read("pmod_output.txt", String)
        rm("pmod_output.txt")  # Clean up temporary files
    
        # Extract onpulse values from the output
        m = match(r"-onpulse '(\d+) (\d+)'", output)
        if !isnothing(m)
            bin_st, bin_end = parse.(Int, m.captures)
            # Ensure the onpulse region length is even
            region_length = bin_end - bin_st + 1
            if region_length % 2 != 0
                println("Warning: Onpulse region length ($region_length) is not even. Adjusting bin_end to make it even.")
                bin_end -= 1  # Adjust bin_end to make the region length even
            end
            println("Found onpulse range: $bin_st to $bin_end")
        end
    
        debased_file = replace(outfile_path, ".spCF" => ".debase.gg")  # Define the debased file path
    
        # Calculate 2dfs and lrfs using pspec
        run(pipeline(`pspec -w -2dfs -lrfs -onpulsed "/NULL" -nfft 256 -onpulse "$(bin_st) $(bin_end)" $debased_file`, stderr="errs.txt"))
    
        # Find P3 value using pspecDetect
        run(pipeline(`pspecDetect -v -device "/xw" $debased_file`, `tee pspecDetect_output.txt`))
    
        # Read captured output
        output = read("pspecDetect_output.txt", String)
        rm("pspecDetect_output.txt")  # Clean up temporary files
    
        # Extract P3 value from the output
        p3_matches = collect(eachmatch(r"P3\[P0\]\s*=\s*(\d+\.\d+)\s*\+-\s*(\d+\.\d+)", output))
        if !isempty(p3_matches)
            last_match = p3_matches[end]
            p3_value = parse(Float64, last_match.captures[1])
            p3_error = parse(Float64, last_match.captures[2])
            println("Found P3 = $p3_value ± $p3_error P0")
        end
    
        # Calculate ybins based on the P3 value
        ybins = round(Int, p3_value * 10)  # Generate ybins dynamically based on P3 value
        println("Number of ybins: $ybins")
    
        # Perform the fold operation using pfold
        run(pipeline(`pfold -p3fold "$p3_value $ybins" -onpulse "$bin_st $bin_end" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $debased_file`, stderr="errs.txt"))
    
        # Return adjusted onpulse region
        return bin_st-20, bin_end+20
    end
    

end # module
