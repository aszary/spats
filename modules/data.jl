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
    try
        # Run the pdv tool to convert the file to ASCII
        run(pipeline(`pdv -t -F -p $infile`, stdout="$outfile", stderr="errs.txt"))
    catch e
        println("Error during conversion of PSRFIT file to ASCII: $e")
        rethrow(e)
    end
    # If you want frequency info, change `-t` to `-A` in the pdv command
    # @showprogress 1 for i in 1:pn  # psrchive indexing
    # end
end

"""
Process data with PSRCHIVE and PSRSALSA
"""
function process_psrdata(indir, outdir; outfile="pulsar.spCF", files=nothing)
    # Check if the output directory exists, if not, create it
    if !isdir(outdir)
        println("Directory $outdir does not exist. Creating directory...")
        mkpath(outdir)  # Creates the output directory if it doesn't exist
        println("Directory $outdir has been created.")
    end

    if files === nothing
        # Find all .spCF files in the input directory
        files = filter(f -> endswith(f, ".spCF"), readdir(indir))
    end

    # Sort files based on pulse numbers (e.g., 00000-00255, 00256-00511, etc.)
    sort!(files, by = f -> begin
        # Extract the pulse range from the filename (e.g., "2019-12-15-03:19:04_00000-00255.spCF")
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

    # Define a unique output folder based on the pulsar and create it
    pulsar_name = "pulsar_output"  # You can change this to use a specific pulsar identifier if needed
    pulsar_folder = joinpath(outdir, pulsar_name)

    if !isdir(pulsar_folder)
        println("Creating unique folder for pulsar: $pulsar_folder")
        mkpath(pulsar_folder)
    else
        println("Directory for pulsar $pulsar_name already exists. Skipping...")
        return  # Skip further processing if the directory already exists
    end

    # Define the output file path inside the pulsar folder
    outfile = joinpath(pulsar_folder, outfile)

    file_names = [joinpath(indir, file) for file in files]

    # Connect all files into the output file
    println("Running psradd to combine files into $outfile...")
    try
        run(pipeline(`psradd $file_names -o $outfile`, stderr="errs.txt"))  # PSRCHIVE
    catch e
        println("Error during psradd execution: $e")
        rethrow(e)
    end

    # Debase the data
    println("Running pmod to debase the data...")
    try
        run(pipeline(`pmod -device "/xw" -debase $outfile`, `tee pmod_output.txt`))
    catch e
        println("Error during pmod execution: $e")
        rethrow(e)
    end

    # Read captured output from pmod
    try
        output = read("pmod_output.txt", String)
        rm("pmod_output.txt")  # Cleanup temporary file
    catch e
        println("Error reading pmod output file: $e")
        rethrow(e)
    end

    # Extract onpulse values from output
    m = match(r"-onpulse '(\d+) (\d+)'", output)
    if !isnothing(m)
        bin_st, bin_end = parse.(Int, m.captures)
        # Check if onpulse region length is even
        region_length = bin_end - bin_st + 1
        if region_length % 2 != 0
            println("Warning: Onpulse region length ($region_length) is not even. Adjusting bin_end to make it even.")
            bin_end -= 1
            println("Adjusted onpulse range: $bin_st to $bin_end")
        end
        println("Found onpulse range: $bin_st to $bin_end")
    else
        println("No onpulse range found in pmod output.")
    end

    debased_file = replace(outfile, ".spCF" => ".debase.gg")

    # Calculate 2dfs and lrfs
    println("Running pspec to calculate 2dfs and lrfs...")
    try
        run(pipeline(`pspec -w -2dfs -lrfs -profd "/NULL" -onpulsed "/NULL" -2dfsd "/NULL" -lrfsd "/NULL" -nfft 256 -onpulse "$(bin_st) $(bin_end)" $debased_file`, stderr="errs.txt"))
    catch e
        println("Error during pspec execution: $e")
        rethrow(e)
    end

    # Find P3
    println("Running pspecDetect to find P3...")
    try
        run(pipeline(`pspecDetect -v -device "/xw" $debased_file`, `tee pspecDetect_output.txt`))
    catch e
        println("Error during pspecDetect execution: $e")
        rethrow(e)
    end

    # Read captured output from pspecDetect
    try
        output = read("pspecDetect_output.txt", String)
        rm("pspecDetect_output.txt")  # Cleanup temporary file
    catch e
        println("Error reading pspecDetect output file: $e")
        rethrow(e)
    end

    # Extract P3 value from the last occurrence
    p3_matches = collect(eachmatch(r"P3\[P0\]\s*=\s*(\d+\.\d+)\s*\+-\s*(\d+\.\d+)", output))
    if !isempty(p3_matches)
        last_match = p3_matches[end]
        p3_value = parse(Float64, last_match.captures[1])
        p3_error = parse(Float64, last_match.captures[2])
        println("Found P3 = $p3_value ± $p3_error P0")
    else
        println("P3 not found in pspecDetect output.")
    end

    # Ensure that the function `Functions.find_ybins(p3_value)` is defined elsewhere in the code
    ybins = Functions.find_ybins(p3_value)
    println("Number of ybins: $ybins")

    # Perform folding on the data
    println("Running pfold for folding the data...")
    try
        run(pipeline(`pfold -p3fold "$p3_value $ybins" -onpulse "$bin_st $bin_end" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $debased_file`, stderr="errs.txt"))
    catch e
        println("Error during pfold execution: $e")
        rethrow(e)
    end

    return bin_st - 20, bin_end + 20
end

    

end # module
