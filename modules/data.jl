module Data
    using Glob
    using FITSIO
    using ProgressMeter

    include("functions.jl")
    include("tools.jl")


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
        run(pipeline(`pdv -t -F $infile`, stdout="$outfile", stderr="errs.txt"))
        # change -t to -A to get frequancy information
        #@showprogress 1 for i in 1:pn  # psrchive indexing
        #end
    end


    """
    Uses PSRCHIVE to add .spCF files 
    
    indir: input directory
    outfile: output file name (full path)
    files: list of files to add (filenames only)
    """
    function add_psrfiles(indir, outfile; files=nothing)

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
    
        file_names = [joinpath(indir, file) for file in files]

        # connecting all files
        run(pipeline(`psradd $file_names -o $outfile`, stderr="errs.txt")) # PSRCHIVE

    end


    function get_nsubint(outfile, params_file, params)
        # gets number of single pulses
        outstr = read(pipeline(`psrstat -c nsubint $outfile`, stderr="errs.txt"), String)
        nsubint = parse(Int, split(outstr, "=")[2])
        params["nsubint"] = nsubint
        println("Number of subintervals: $nsubint")
        Tools.save_params(params_file, params)
    end
 
    
    """
    Debase the data
    """
    function debase(outfile, params_file, params)

        if isnothing(params["bin_st"]) || isnothing(params["bin_end"])

            run(pipeline(`pmod -device "/xw" -iformat PSRFITS -debase $outfile`, `tee pmod_output.txt`))

            # Read captured output to get bin_st and bin_end
            output = read("pmod_output.txt", String)
            rm("pmod_output.txt")  # cleanup

            # Extract onpulse values
            m = match(r"-onpulse '(.+) (.+)'", output)
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
                params["bin_st"] = bin_st
                params["bin_end"] = bin_end
                Tools.save_params(params_file, params)
            end
        else
            run(pipeline(`pmod -onpulse "$(params["bin_st"]) $(params["bin_end"])" -device "/NULL" -iformat PSRFITS -debase $outfile`))
        end
    end


    """
    Calculate 2dfs and lrfs using PSRSALSA
    """
    function twodfs_lrfs(debased_file, params_file, p; detect=false)

        # Calculate 2dfs and lrfs
        run(pipeline(`pspec -w -2dfs -lrfs -profd "/NULL" -onpulsed "/NULL" -2dfsd "/NULL" -lrfsd "/NULL" -nfft $(p["nfft"]) -onpulse "$(p["bin_st"]) $(p["bin_end"])" $debased_file`,  stderr="errs.txt"))

        if p["p3"] == -1.0 || detect == true
            # Find P3
            run(pipeline(`pspecDetect -v -device "/xw" $debased_file`, `tee pspecDetect_output.txt`))

            # Read captured output
            output = read("pspecDetect_output.txt", String)
            rm("pspecDetect_output.txt")  # cleanup

            # Extract P3 value from the last occurrence
            p3_matches = collect(eachmatch(r"P3\[P0\]\s*=\s*(\d+\.\d+)\s*\+-\s*(\d+\.\d+)", output))
            if !isempty(p3_matches)
                last_match = p3_matches[end]
                p3_value = parse(Float64, last_match.captures[1])
                p3_error = parse(Float64, last_match.captures[2])
                println("Found P3 = $p3_value ± $p3_error P0")
            end

            ybins = Functions.find_ybins(p3_value)
            println("Number of ybins: $ybins")
            p["p3"] = p3_value
            p["p3_error"] = p3_error
            p["p3_ybins"] = ybins
            Tools.save_params(params_file, p)
        end
    end


    """
    Process data with PSRCHIVE and PSRSALSA
    """
    function process_psrdata(indir, outdir; files=nothing, outfile="pulsar.spCF", params_file="params.json")

        # full paths
        params_file = joinpath(outdir, params_file)
        outfile = joinpath(outdir, outfile)
        debased_file = replace(outfile, ".spCF" => ".debase.gg")

        # check if params_file exists if not creating default one
        if !isfile(params_file)
            println("File $params_file does not exist, creating default one.")
            p = Tools.default_params(params_file)
        else
            p = Tools.read_params(params_file)
        end

        # add all .spCF files 
        add_psrfiles(indir, outfile; files=files)

        # gets number of single pulses if needed
        if isnothing(p["nsubint"])
            get_nsubint(outfile, params_file, p)
        end

        # debase the data
        debase(outfile, params_file, p)

        # Calculate 2dfs and lrfs
        twodfs_lrfs(debased_file, params_file, p; detect=false)

        # calculate p3-folded profile
        run(pipeline(`pfold  -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $debased_file`,  stderr="errs.txt"))

        return p
    end


    function process_psrfit_files(base_dir, output_dir; name_mod=nothing)
        # Step 1: Extract base directory 
        name = basename(base_dir)

        # Step 2: Create output subdirectory
        output_subdir = joinpath(output_dir, name)
        println(output_subdir)
        if !isdir(output_subdir)
            mkpath(output_subdir)
            println("Created output directory: ", output_subdir)
        end

        # Step 3 create JSON file name
        params_file = joinpath(output_subdir, "params.json")

        #Step 3.5 Check if already processed
        if isdir(output_subdir) && isfile(joinpath(output_subdir, "params.json"))
            println("Skipping already processed catalogue: $name")
            return
        else
            if isfile(params_file)
                p = Tools.read_params(params_file)
            else
                p = Tools.default_params(params_file)
                println("New file with parameters created: $params_file")
            end
        end

        # Step 4: Find second-level catalogue
        second_catalogue = joinpath(base_dir, readdir(base_dir)[1])
        println("Second-level catalogue found: ", second_catalogue)
    
        # Step 5: Find .spCF files
        spcf_files = filter(f -> occursin("spCF", f), readdir(second_catalogue, join=true))
        Data.sort!(spcf_files)

        #step 6: combining .spCF into one 
        output_file = joinpath(output_subdir, "converted.spCF")
        file_names = [joinpath(name, file) for file in spcf_files]

        run(pipeline(`psradd $file_names -o $output_file`, stderr="errs.txt"))
        out_txt = replace(output_file ,".spCF" => ".txt")


        # Step 7: Convert spCF -> ascii
        Data.convert_psrfit_ascii(output_file, out_txt)
        rm(output_file) #cleanup .spCF file 


        # Step 8: Debase the combined ASCII file using pmod
        #debased_file = replace(out_txt, ".txt" => ".debase.gg")
        run(pipeline(`pmod -device "/xw" -debase $out_txt`, `tee pmod_output.txt`))
        # Read captured output and cleanup
        output = read("pmod_output.txt", String)
        m = match(r"-onpulse '(\d+) (\d+)'", output)
            if !isnothing(m)
                bin_st, bin_end = parse.(Int, m.captures)
                # Ensure onpulse region length is even
                region_length = bin_end - bin_st + 1
                if region_length % 2 != 0
                  println("Warning: Onpulse region length ($region_length) is not even. Adjusting bin_end to make it even.")
                   bin_end -= 1
                   println("Adjusted onpulse range: $bin_st to $bin_end")
                end
                p["bin_st"] = bin_st
                p["bin_end"] = bin_end
                println("Found onpulse range: $bin_st to $bin_end")
            end
            Tools.save_params(params_file, p)
            println("Parameters updated and saved to $params_file")
        end
        #rm("pmod_output.txt")  # Cleanup 

        return

        # Step 9: Load combined data
        combined_data = Data.load_ascii(out_txt)



           #= Tools.save_params(params_file, p)
            println("Parameters updated and saved to $params_file")
            =#
    
        # Step 10: Plot
        Plot.single(combined_data, output_subdir; darkness=0.5, bin_st= p["bin_st"], bin_end= p["bin_end"], start= p["pulse_start"], number= (p["pulse_end"] - p["pulse_start"]), name_mod=name_mod, show_=false)

        Plot.lrfs(combined_data, output_subdir; darkness=0.1, start= p["pulse_start"], bin_st= p["bin_st"], bin_end= p["bin_end"], name_mod=name_mod, change_fftphase=false, show_=false)

        #Plot.average(combined_data, output_subdir; bin_st=p["bin_st"], bin_end=p["bin_end"], number=(p["pulse_end"]-p["pulse_start"]), name_mod=name_mod, show_=false)
        
        #=
        Plot.single(combined_data, output_subdir, darkness=0.5, bin_st=bin_st, bin_end=bin_end, number=nothing, name_mod=name_mod, show_=false)
        Plot.lrfs(combined_data, output_subdir, darkness=0.1, start=1, bin_st=bin_st, bin_end=bin_end, name_mod=name_mod, change_fftphase=false, show_=false)
        Plot.average(combined_data, output_subdir,bin_st=bin_st, bin_end=bin_end, number=nothing, name_mod=name_mod, show_=false)
        =#
    end
    

    """

    """
    function process_all_data(outdir; base_root="/home/psr/data/new")
        # Get all subdirectories in base_root
        catalogues = filter(isdir, readdir(base_root, join=true))
    
        if isempty(catalogues)
            println("No data found in $base_root. Exiting...")
            return
        end
    
        for catalogue in catalogues
            name = basename(catalogue)  # Extract directory name
            println("Processing catalogue: ", name)
            process_psrfit_files(catalogue, outdir, name_mod=name)
        end
    end

end # module
