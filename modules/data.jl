module Data
    using Glob
    using FITSIO
    using ProgressMeter
    using CairoMakie, FileIO
    using PDFIO
    using Statistics
    using FFTW
    using PyPlot
    
    using LsqFit, DSP 
    using DelimitedFiles, Printf
    
    include("functions.jl")
    include("tools.jl")
    include("plot.jl")
    include("profile_metrics.jl")
    using .ProfileMetrics


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
        Loads PSRCHIVE ASCII file with all 4 polarizations.
        Returns a 3D array: (pulse, bin, pol)
    """
    function load_ascii_all(infile)
        f = open(infile)
        lines = readlines(f)
        
        header = split(lines[1])
        pulses = parse(Int, header[6])    # Nsub
        bins = parse(Int, header[12])     # Nbin
        pols = parse(Int, header[10])     # Npol

        data = Array{Float64}(undef, pulses, bins, pols)

        for i in 2:length(lines)
            res = split(lines[i])
            pulse = parse(Int, res[1]) + 1
            bin = parse(Int, res[3]) + 1
            for pol in 1:pols
                val = parse(Float64, res[3 + pol])
                data[pulse, bin, pol] = val
            end
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
    txt: convert to txt?
    """
    function add_psrfiles(indir, outfile; files=nothing, sixteen=false)

        if files === nothing
            # Find all .spCF files in the input directory
            if sixteen == false
                files = filter(f -> endswith(f, ".spCF"), readdir(indir))
            else
                files = filter(f -> endswith(f, ".spCf16"), readdir(indir))
            end
        end
        
        # Sort files based on pulse numbers (e.g., 00000-00255, 00256-00511, etc.)
        sort!(files, by = f -> begin
            # Extract the pulse range from the filename (e.g., "2019-12-15-03:19:04_00000-00255.spCF")
            m = match(r"_(\d+)-(\d+)\.spC(?:F|f16)$", f)
            if isnothing(m)
                return typemax(Int)  # Files without proper format go to the end
            else
                return parse(Int, m.captures[1])  # Sort by the starting pulse number
            end
        end)
    
        if isempty(files)
            error("No .spCF or .spCf16 files found in directory: $indir")
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

    # Arguments
    """
    function debase(infile, params_file, params)


        if isnothing(params["bin_st"]) || isnothing(params["bin_end"])

            println("pmod -device \"/xw\" -debase $infile")
            run(pipeline(`pmod -device "/xw" -debase $infile`, `tee pmod_output.txt`))
            #run(pipeline(`pmod -device "/xw" -iformat PSRFITS -debase $outfile`, `tee pmod_output.txt`))

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
            # ASCII single pulses
            convert_psrfit_ascii(replace(infile, ".spCF"=>".debase.gg"), replace(infile, ".spCF"=>".debase.txt"))

        else
            run(pipeline(`pmod -onpulse "$(params["bin_st"]) $(params["bin_end"])" -device "/NULL" -debase $infile`))
            #run(pipeline(`pmod -onpulse "$(params["bin_st"]) $(params["bin_end"])" -device "/NULL" -iformat PSRFITS -debase $outfile`))
            # ASCII single pulses
            convert_psrfit_ascii(replace(infile, ".spCF"=>".debase.gg"), replace(infile, ".spCF"=>".debase.txt"))
        end
    end


    """
    Calculate 2dfs and lrfs using PSRSALSA

    # ADD description!
    debased_file ??
    """
    function twodfs_lrfs(debased_file, params_file, p; detect=false)

        # Calculate 2dfs and lrfs
        run(pipeline(`pspec -w -oformat ASCII -2dfs -lrfs -profd "/NULL" -onpulsed "/NULL" -2dfsd "/NULL" -lrfsd "/NULL" -nfft $(p["nfft"]) -onpulse "$(p["bin_st"]) $(p["bin_end"])" $debased_file`,  stderr="errs.txt"))

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
    Process data with my old routines
    """
    function process_data_andrzej(debased_file, outdir, p)

        println("DE $debased_file")
        # all data
        da0 = Data.load_ascii(replace(debased_file, ".gg"=>".txt"))  
        Plot.single(da0, outdir; darkness=0.7, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="all", show_=true)

        folded = Tools.p3fold(da0, p["p3"],  p["p3_ybins"])
        Plot.single(folded, outdir; darkness=0.97, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="p3fold", show_=true, repeat_num=4)

        println("Polarization cleaning threshold: ", p["clean_threshold"])
        # polarisation cleaning
        da = Data.load_ascii_all(replace(debased_file, ".gg"=>".txt"))  
        da2 = clean(da; threshold=p["clean_threshold"]) 
        Plot.single(da2, outdir; darkness=0.7, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="cleaned", show_=true)

        folded = Tools.p3fold(da2, p["p3"],  p["p3_ybins"])
        Plot.single(folded, outdir; darkness=0.97, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="p3fold_clean", show_=true, repeat_num=4)

        
    end



    """
    Process data with PSRCHIVE and PSRSALSA
    """
    function process_psrdata(indir, outdir; files=nothing, outfile="pulsar.spCF", params_file="params.json")

        # full paths
        params_file = joinpath(outdir, params_file)
        outfile = joinpath(outdir, outfile)
        debased_file = replace(outfile, ".spCF" => ".debase.gg")

        # ensure directory exists
        params_dir = dirname(params_file)
        if !isdir(params_dir)
            println("Directory $params_dir does not exist, creating it.")
            mkpath(params_dir)
        end

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
        println("P3-folding with:")
        # TODO experiment here
        println("pfold -p3fold_noonpulse -p3fold \"$(p["p3"]) $(p["p3_ybins"])\" -onpulse \"$(p["bin_st"]) $(p["bin_end"])\" -onpulsed \"/NULL\" -p3foldd \"/NULL\" -w -oformat ascii $debased_file")
        run(pipeline(`pfold  -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $debased_file`,  stderr="errs.txt"))

        return p, debased_file, outdir
    end


    function plot_psrdata(data_dir, params)

        p = params

        d = load_ascii(joinpath(data_dir, "pulsar.debase.txt"))

        # Single pulses
        Plot.single(d, data_dir; darkness=0.5, bin_st= p["bin_st"], bin_end=p["bin_end"], start=p["pulse_start"], number= (p["pulse_end"] - p["pulse_start"]), name_mod="pulsar", show_=true)

        # 2DFS
        d2 = load_ascii(joinpath(data_dir, "pulsar.debase.1.2dfs"))
        Plot.twodfs(d2, data_dir, p; name_mod="pulsar", show_=true, average=Tools.average_profile(d))

        # LRFS
        d3 = load_ascii_all(joinpath(data_dir, "pulsar.debase.lrfs"))
        Plot.lrfs(d3, data_dir, p; name_mod="pulsar", show_=true)

        # P3-folded
        d5 = Data.load_ascii(joinpath(data_dir, "pulsar.debase.p3fold"))
        Plot.p3fold(d5, data_dir; start=1, bin_st=p["bin_st"]-20, bin_end=p["bin_end"]+20, name_mod="pulsar", show_=true, repeat_num=4)
 
    end


    function process_all_data(outdir; base_root="/home/psr/data/new")
        # Get all subdirectories in base_root
        dirs = filter(isdir, readdir(base_root, join=true))

        if isempty(dirs)
            println("No data found in $base_root. Exiting...")
            return
        end
        
        for dir in dirs
            name = basename(dir)  # Extract directory name -> pulsar name
            data_dir = joinpath(dir, readdir(dir)[1]) # one level in
            out_dir = joinpath(outdir, name)

            println("\n\n"*"#"^79*"\n"*"#"^79*"\n"*"#"^79)
            println("Processing pulsar: ", name)
            println(out_dir)
            
            # TODO temporatory fix
            if isfile(joinpath(out_dir,"pulsar_single.pdf"))
                println("SKIPPING pulsar $name !!!!\n\n")
                continue 
            end

            # Creates out_dir if does not exists...
            if !isdir(out_dir)
                mkpath(out_dir)
                println("Created output directory: ", out_dir)
            end

            # TODO temporatory fix
            try
                # process data
                p = process_psrdata(data_dir, out_dir)
            
                # plot data
                plot_psrdata(out_dir, p)
            catch
                println("ERROR in pulsar $name !!!!\n\n")
            end
        end
    end

    function combine_pngs(outdir)
        paths = filter(isdir, readdir(outdir, join=true))
        pulsars = basename.(paths)

        pdfs = []

        # create pdf for all pulsars
        for (i,pulsar) in enumerate(pulsars)

            fig = Figure(size = (500, 800))
            Label(fig[0, :], text = pulsar, fontsize = 20, tellwidth = false)

            ax1 = Axis(fig[1, 1])
            ax2 = Axis(fig[1, 2])
            ax3 = Axis(fig[2, 1])
            ax4 = Axis(fig[2, 2])
            ax5 = Axis(fig[3, 1])
            ax6 = Axis(fig[3, 2])


            # load images
            try
                img1 = load(joinpath(paths[i], "pulsar_single.png"))
                image!(ax1, rotr90(img1))
            catch e
                println("No pulsar_single.png for $pulsar")
            end
            try
                img2 = load(joinpath(paths[i], "pulsar_lrfs.png"))
                image!(ax2, rotr90(img2))
            catch e
                println("No pulsar_lrfs.png for $pulsar")
            end
            try
                img3 = load(joinpath(paths[i], "pulsar_2dfs.png"))
                image!(ax3, rotr90(img3))
            catch e
                println("No pulsar_2dfs.png for $pulsar")
            end
            try
                img4 = load(joinpath(paths[i], "pulsar_p3fold.png"))
                image!(ax4, rotr90(img4))
            catch e
                println("No pulsar_p3fold.png for $pulsar")
            end
            try
                img5 = load(joinpath(paths[i], "all_single.png"))
                image!(ax5, rotr90(img5))
            catch e
                println("No  all_single.png for $pulsar")
            end
            try
                img6 = load(joinpath(paths[i], "p3fold_single.png"))
                image!(ax6, rotr90(img6))
            catch e
                println("No p3fold_single.png for $pulsar")
            end


            hidedecorations!(ax1)
            hidedecorations!(ax2)
            hidedecorations!(ax3)
            hidedecorations!(ax4)
            hidedecorations!(ax5)
            hidedecorations!(ax6)


            filename = joinpath(outdir, "page_$i.pdf")
            save(filename, fig)
            push!(pdfs, basename(filename))
        end

        # combine pdfs
        cmd = vcat(["pdfunite"], pdfs, ["all.pdf"])
        write("convert.txt", join(cmd, " "))
        println(join(cmd, " "))

    end


    function combine_pngs_to_pdf(out_dir::String)

        # Collect and sort PNG files recursively
        png_paths = sort([
            joinpath(root, file)
            for (root, _, files) in walkdir(out_dir)
            for file in files
            if endswith(lowercase(file), ".png")
        ])

        isempty(png_paths) && return

        pulsar_name = basename(out_dir)
        output_pdf_dir = "/home/psr/output/"
        mkpath(output_pdf_dir)

        pages = Figure[]
        pdf_paths = String[]
        images_per_page = 4
        fig=nothing

        for (i, path) in enumerate(png_paths)
            page_index = div(i - 1, images_per_page) + 1
            position_in_page = (i - 1) % images_per_page

            if position_in_page == 0
                fig = Figure(size = (400, 400))
            end

            row = div(position_in_page, 2) + 1
            col = mod(position_in_page, 2) + 1

            ax = Axis(fig[row, col];
                xticksvisible = false, 
                yticksvisible = false,
                xgridvisible = false, 
                ygridvisible = false,
                leftspinevisible = false, 
                rightspinevisible = false,
                topspinevisible = false, 
                bottomspinevisible = false,
                xlabelvisible = false, 
                ylabelvisible = false,
                titlevisible = false,
            )

            image!(ax, load(path))

            if (i % images_per_page == 0) || (i == length(png_paths))
                push!(pages, fig)
                pdf_full_path = joinpath(output_pdf_dir, "$(pulsar_name)_page_$(page_index).pdf")
                pdf_name = "$(pulsar_name)_page_$(page_index).pdf"
                CairoMakie.save(pdf_full_path, fig)
                push!(pdf_paths, pdf_name)
            end
        end

        # Combine all PDF pages into one file using pdfunite
        combined_pdf_path = joinpath(output_pdf_dir, "combined.pdf")
        cmd_args = vcat(["pdfunite"], pdf_paths, [combined_pdf_path])
        println("Running command: ", join(cmd_args, " "))
        write("convert.txt", join(cmd_args, " "))
        #run(`$(cmd_args...)`)
        #println("Combined PDF created at: $combined_pdf_path")
    end


    """
        clean(data; threshold=0.7)

    Filter polarization data based on the ratio of linear polarization to intensity.

    # Arguments
    - `data`: 3D array with dimensions (num, bins, npol) containing polarization data
    - `data[:, :, 1]`: Stokes I intensity
    - `data[:, :, 2]`: Stokes Q parameter  
    - `data[:, :, 3]`: Stokes U parameter
    - `threshold`: filtering threshold (default 0.7) - minimum value for L/I ratio

    # Returns
    - `data_clean`: 2D array containing only intensity I for points meeting the L/I > threshold criterion

    # Description
    The function calculates linear polarization L = √(Q² + U²) for each data point and retains 
    only those points where the ratio of linear polarization to intensity exceeds the given threshold.
    """
    function clean(data; threshold=0.7)
        # four polarisation data
        num, bins, npol = size(data)
        data_clean = zeros(num, bins)
        for i in 1:num
            for j in 1:bins 
                I = data[i, j, 1]
                Q = data[i, j, 2]
                U = data[i, j, 3]
                L = sqrt(Q^2 + U^2)
                if L / I > threshold
                    data_clean[i, j] = I
                end
            end
        end
        return data_clean
    end


    """
    just some quick cleaning..
    """
    function remove_folders(dirname)
        # Get all subdirectories in dirname
        dirs = filter(isdir, readdir(dirname, join=true))
        for dir in dirs
            # Get all files in the subdirectory
            files = filter(isfile, readdir(dir, join=true))
            filename = basename(files[1])
            if length(files)==1 && filename=="params.json"
                println(basename(dir))
                rm(dir, recursive=true)
            end
        end
    end

    """
    some more cleaning
    """
    function remove_notinteresting(interesting_filename, data_dir)

        keep_dirs = Set(readlines(interesting_filename))
        all_dirs = filter(isdir, readdir(data_dir, join=true))

        for dir in all_dirs
            dir_name = basename(dir)
            if dir_name ∉ keep_dirs
                rm(dir, recursive=true, force=true)
            end
        end
        
    end


    """
    Process multifrequency data with PSRCHIVE and PSRSALSA
    """
    function process_psrdata_16(indir, outdir; files=nothing, outfile="pulsar.spCf16", params_file="params.json")

        # full paths
        params_file = joinpath(outdir, params_file)
        outfile = joinpath(outdir, outfile)

        # if no catalog
        if !isdir(outdir)
            mkdir(outdir)
        end

        # check if params_file exists if not creating default one
        if !isfile(params_file)
            println("File $params_file does not exist, creating default one.")
            p = Tools.default_params(params_file)
        else
            p = Tools.read_params(params_file)
        end

        # add all .spCf16 files 
        add_psrfiles(indir, outfile; files=files, sixteen=true)

        if haskey(p, "zaps") && !isnothing(p["zaps"])
            println("ZAPPING subints starts")
            run(pipeline(`paz -m -w "$(p["zaps"])" $outfile`)) # zero weights  # this is the way
            println("ZAPPING subints ends")
        end
        #run(pipeline(`paz -m -W "49 194" -W "532 580" -W "880 982" $outfile`)) # zero weights  # this is the way

        # divide to two frequencies
        low_filename, high_filename, mid_filename = multifrequency_split(outfile)

        # gets number of single pulses if needed
        if isnothing(p["nsubint"])
            get_nsubint(low_filename, params_file, p)
        end

        # debase the data
        debase_16(low_filename, high_filename, mid_filename, params_file, p)

        # Calculate 2dfs and lrfs => finds P3
        low_debase = replace(low_filename, ".low"=>"_low.debase.gg")
        twodfs_lrfs(low_debase, params_file, p; detect=false)
        high_debase = replace(high_filename, ".high"=>"_high.debase.gg")
        twodfs_lrfs(high_debase, params_file, p; detect=false)

        # single pulses - low
        da_lo = Data.load_ascii(replace(low_filename, ".low"=>"_low_debase.txt"))
        Plot.single(da_lo, outdir; darkness=0.7, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="low", show_=true)

        # single pulses - high
        da_hi = Data.load_ascii(replace(high_filename, ".high"=>"_high_debase.txt"))
        Plot.single(da_hi, outdir; darkness=0.7, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="high", show_=true)

        # low freq p3-fold with cleaning
        da_low = Data.load_ascii_all(replace(low_filename, ".low"=>"_low_debase.txt"))
        da_low = clean(da_low; threshold=p["clean_threshold"]) 
        folded_low = Tools.p3fold(da_low, p["p3"],  p["p3_ybins"])
        Plot.single(folded_low, outdir; darkness=0.9, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="p3fold_low", show_=true, repeat_num=4)

        # high freq p3-fold with cleaning
        da_high = Data.load_ascii_all(replace(high_filename, ".high"=>"_high_debase.txt"))
        da_high = clean(da_high; threshold=p["clean_threshold"]) 
        folded_high = Tools.p3fold(da_high, p["p3"],  p["p3_ybins"])
        Plot.single(folded_high, outdir; darkness=0.9, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="p3fold_high", show_=true, repeat_num=4)

        # PSRSALSA p3folding low freq
        run(pipeline(`pfold -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $low_debase`,  stderr="errs.txt"))
        mv(replace(low_debase, ".gg"=>".p3fold"), replace(low_debase, ".gg"=>".p3fold_refine"), force=true)
        d5 = Data.load_ascii(replace(low_debase, ".gg"=>".p3fold_refine"))
        Plot.p3fold(d5, outdir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low", show_=true, repeat_num=4)
        # norefine
        #println("pfold -p3fold_norefine -p3fold \"$(p["p3"]) $(p["p3_ybins"])\" -onpulse \"$(p["bin_st"]) $(p["bin_end"])\" -onpulsed \"/NULL\" -p3foldd \"/NULL\" -w -oformat ascii $low_debase")
        run(pipeline(`pfold -p3fold_norefine -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $low_debase`,  stderr="errs.txt"))
        mv(replace(low_debase, ".gg"=>".p3fold"), replace(low_debase, ".gg"=>".p3fold_norefine"), force=true)
        d5 = Data.load_ascii(replace(low_debase, ".gg"=>".p3fold_norefine"))
        Plot.p3fold(d5, outdir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low_norefine", show_=true, repeat_num=4)
 
        # PSRSALSA p3folding high freq
        run(pipeline(`pfold -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $high_debase`,  stderr="errs.txt"))
        #run(pipeline(`pfold -p3fold_nritt 50 -p3fold_cpb 50 -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $high_debase`,  stderr="errs.txt"))
        mv(replace(high_debase, ".gg"=>".p3fold"), replace(high_debase, ".gg"=>".p3fold_refine"), force=true)
        d5 = Data.load_ascii(replace(high_debase, ".gg"=>".p3fold_refine"))
        Plot.p3fold(d5, outdir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_high", show_=true, repeat_num=4)
        # norefine
        run(pipeline(`pfold -p3fold_norefine -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $high_debase`,  stderr="errs.txt"))
        mv(replace(high_debase, ".gg"=>".p3fold"), replace(high_debase, ".gg"=>".p3fold_norefine"), force=true)
        d5 = Data.load_ascii(replace(high_debase, ".gg"=>".p3fold_norefine"))
        Plot.p3fold(d5, outdir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_high_norefine", show_=true, repeat_num=4)

        return p, nothing, outdir

    end


    function multifrequency_split(spCf16_file)
        println(spCf16_file)

        low = replace(spCf16_file, ".spCf16"=>".low")
        low_txt = replace(spCf16_file, ".spCf16"=>"_low.txt")
        high = replace(spCf16_file, ".spCf16"=>".high")
        high_txt = replace(spCf16_file, ".spCf16"=>"_high.txt")
        mid = replace(spCf16_file, ".spCf16"=>".mid")
        mid_txt = replace(spCf16_file, ".spCf16"=>"_mid.txt")

        #run(pipeline(`paz -Z 0-7 -e high $spCf16_file`,stderr="errs.txt"))
        #run(pipeline(`paz -Z 8-15 -e low $spCf16_file`,stderr="errs.txt"))

        #run(pipeline(`paz -Z 0-12 -e high $spCf16_file`,stderr="errs.txt"))
        #run(pipeline(`paz -Z 3-15 -e low $spCf16_file`,stderr="errs.txt"))
        run(pipeline(`paz -Z 0-6 -Z 9-15 -e mid $spCf16_file`,stderr="errs.txt"))
        #run(pipeline(`paz -Z 0-0 -e mid $spCf16_file`,stderr="errs.txt"))
        
        run(pipeline(`paz -Z 0-9 -e high $spCf16_file`,stderr="errs.txt"))
        run(pipeline(`paz -Z 6-15 -e low $spCf16_file`,stderr="errs.txt"))

        # fscrunch
        run(pipeline(`pam -F -e high $high`,stderr="errs.txt"))
        run(pipeline(`pam -F -e low $low`,stderr="errs.txt"))
        run(pipeline(`pam -F -e mid $mid`,stderr="errs.txt"))

        run(pipeline(`pdv -A -F $high`, stdout=high_txt, stderr="errs.txt"))
        run(pipeline(`pdv -A -F $low`, stdout=low_txt, stderr="errs.txt"))
        run(pipeline(`pdv -A -F $mid`, stdout=mid_txt, stderr="errs.txt"))
        #println("done")
        # change -t to -A to get frequancy information
        return low, high, mid

    end


    """
    Debase the data

    # Arguments
    """
    function debase_16(low, high, mid, params_file, params)

        # finding bin_st and bin_end based on low frequency file
        if isnothing(params["bin_st"]) || isnothing(params["bin_end"])

            println("pmod -device \"/xw\" -debase $low")
            run(pipeline(`pmod -device "/xw" -debase $low`, `tee pmod_output.txt`))
            #run(pipeline(`pmod -device "/xw" -iformat PSRFITS -debase $outfile`, `tee pmod_output.txt`))

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
            # ASCII single pulses
            convert_psrfit_ascii(replace(low, ".low"=>".debase.gg"), replace(low, ".low"=>"_low_debase.txt"))

        else
            run(pipeline(`pmod -onpulse "$(params["bin_st"]) $(params["bin_end"])" -device "/NULL" -debase $low`))
            #run(pipeline(`pmod -onpulse "$(params["bin_st"]) $(params["bin_end"])" -device "/NULL" -iformat PSRFITS -debase $outfile`))
            # ASCII single pulses
            convert_psrfit_ascii(replace(low, ".low"=>".debase.gg"), replace(low, ".low"=>"_low_debase.txt"))
        end

        # changing low frequency psrfit filename
        mv(replace(low, ".low"=>".debase.gg"), replace(low, "pulsar.low"=>"pulsar_low.debase.gg"), force=true)

        # high frequency runs
        run(pipeline(`pmod -onpulse "$(params["bin_st"]) $(params["bin_end"])" -device "/NULL" -debase $high`))
        convert_psrfit_ascii(replace(high, ".high"=>".debase.gg"), replace(high, ".high"=>"_high_debase.txt"))
        # changing high frequency psrfit filename
        mv(replace(high, ".high"=>".debase.gg"), replace(high, "pulsar.high"=>"pulsar_high.debase.gg"), force=true)

        # mid frequency runs
        run(pipeline(`pmod -onpulse "$(params["bin_st"]) $(params["bin_end"])" -device "/NULL" -debase $mid`))
        convert_psrfit_ascii(replace(mid, ".mid"=>".debase.gg"), replace(mid, ".mid"=>"_mid_debase.txt"))
        # changing high frequency psrfit filename
        mv(replace(mid, ".mid"=>".debase.gg"), replace(mid, "pulsar.mid"=>"pulsar_mid.debase.gg"), force=true)

    end


    """
    Analyse p3folds

    # Arguments
    - `indir`: data input directory
    - `type`: type of p3folds "refine" or "norefine" filenames are based on that
 
    """
    function analyse_p3folds_16(indir, type)

        # parameters file 
        p = Tools.read_params(joinpath(indir, "params.json"))

        # low, high frequancy filename end
        low = joinpath(indir, "pulsar_low.debase.p3fold_" * type)
        high = joinpath(indir, "pulsar_high.debase.p3fold_" * type)

        #println(low)
        l = Data.load_ascii(low)
        Plot.p3fold(l, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low_$(type)_analyse", show_=true, repeat_num=1)

        #println(high)
        h = Data.load_ascii(high)
        Plot.p3fold(h, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_high_$(type)_analyse", show_=true, repeat_num=1)
 
        #diff = normalize_01(l) .- normalize_01(h)
        diff = normalize_per_pulse(l) .- normalize_per_pulse(h)
        Plot.p3fold(diff, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low_high_$(type)_analyse", show_=true, repeat_num=1)

        nl = normalize_per_pulse(l)
        nh = normalize_per_pulse(h)
        shifts, corr_values = find_shift_per_pulse(nl, nh)

        println("Mean shift: $(mean(shifts)) bins")
        println("Median shift: $(median(shifts)) bins")

        Plot.analyse_p3folds(nl, nh, p)
        Plot.analyse_p3folds2(nl, nh, p)

    end


    




    """
    Main function called from SpaTs.main() for robust offset analysis.
    Uses non-parametric Fourier Phase Gradient and Barycenter methods.
    
    # Arguments
    - `indir`: data input directory
    - `type`: type of p3folds ("refine" or "norefine")
    - `period`: pulsar period (not used currently, for future reference)
    - `n_comp`: number of components to detect (default: 3)
    """
    function analyse_offsets(indir, type, period; n_comp=3)
        """
        Fourier-based component offset analysis (PHASE-RESOLVED).
        
        Methodology:
        - Global Shift: Computed via Fourier Phase Gradient (Shift Theorem).
        - Component Centroids: Computed using Barycenter (Center of Mass) algorithm.
          This provides a model-independent measure of the energy balance, 
          which is physically robust for asymmetric pulsar profiles.
        
        Output axes:
        - X: Longitude in DEGREES (matching standard pulsar timing literature)
        - Y: Phase offset in DEGREES
        """
        
        # --- 1. Configuration and Data Loading ---
        p = Tools.read_params(joinpath(indir, "params.json"))
        low_path = joinpath(indir, "pulsar_low.debase.p3fold_" * type)
        high_path = joinpath(indir, "pulsar_high.debase.p3fold_" * type)
        
        l_matrix = load_ascii(low_path)
        h_matrix = load_ascii(high_path)
        
        n_phases, n_bins = size(l_matrix)
        
        println("\n" * "="^60)
        println(">>> Starting Phase-Resolved Fourier Analysis: $indir")
        println(">>> P3-folded matrix size: $n_phases phases × $n_bins bins")
        
        pulsar_name = basename(rstrip(indir, '/'))
        
        # Get on-pulse region (Region of Interest)
        bin_st  = get(p, "bin_st", 1)
        bin_end = get(p, "bin_end", n_bins)
        bin_st  = isnothing(bin_st) ? 1 : bin_st
        bin_end = isnothing(bin_end) ? n_bins : bin_end
        
        println(">>> On-Pulse region: bins $bin_st to $bin_end")
        
        # Conversion factors
        db = (bin_end + 1) - bin_st
        bins_to_deg = 360.0 / n_bins
        
        println(">>> Degrees per bin: $(round(bins_to_deg, digits=4))°")
        
        # --- 2. Integrated Profile Analysis (Global reference) ---
        # We calculate the mean profile across all P3 phases to establish 
        # a high-SNR reference for component windows.
        mean_l_roi = copy(normalize_01(vec(mean(l_matrix, dims=1))))
        mean_h_roi = copy(normalize_01(vec(mean(h_matrix, dims=1))))
        
        if bin_st > 1
            mean_l_roi[1:Int(bin_st)-1] .= 0.0
            mean_h_roi[1:Int(bin_st)-1] .= 0.0
        end
        if bin_end < n_bins
            mean_l_roi[Int(bin_end)+1:end] .= 0.0
            mean_h_roi[Int(bin_end)+1:end] .= 0.0
        end
        
        # Global offset using Fourier method
        global_offset_res = ProfileMetrics.global_fourier_shift(mean_l_roi, mean_h_roi)
        global_offset_bins = global_offset_res.shift
        @printf(">>> Global Fourier Shift (High vs Low): %.4f ± %.4f bins\n", global_offset_bins, global_offset_res.error)
        
        # Detect static component boundaries (Barycenters) on the integrated profile
        # This prevents the algorithm from jumping to noise spikes in low-SNR individual phases.
        ref_comps = ProfileMetrics.component_barycenters(mean_l_roi, mean_h_roi; threshold=0.10, n_comp=n_comp)
        n_comps = length(ref_comps)
        windows = [(c.left_bound, c.right_bound) for c in ref_comps]
        
        println(">>> Detected $n_comps components using Barycenter Method.")
        
        # --- 3. Plotting: Integrated Profile ---
        PyPlot.figure(figsize=(10, 9))
        colors = ["#2196F3", "#E65100", "#4CAF50", "#9C27B0", "#F44336", "#795548"]
        
        # Przejście z binów na fizyczne stopnie długości geograficznej (Longitude) dla poprawności naukowej
        lon_arr = (1:n_bins) .* bins_to_deg
        
        PyPlot.subplot(3, 1, 1)
        PyPlot.title("Low Frequency (Integrated Profile & Barycentric Windows)", fontsize=13, fontweight="bold")
        PyPlot.plot(lon_arr, mean_l_roi, color="#1f77b4", linewidth=2.0, label="Low Data")
        PyPlot.fill_between(lon_arr, zeros(length(lon_arr)), mean_l_roi, color="#1f77b4", alpha=0.15)
        PyPlot.ylabel("Normalized Intensity", fontsize=11)
        PyPlot.grid(true, linestyle=":", alpha=0.6)
        PyPlot.minorticks_on()
        for (k, c) in enumerate(ref_comps)
            col = colors[mod1(k, length(colors))]
            PyPlot.axvspan(c.left_bound * bins_to_deg, c.right_bound * bins_to_deg, color=col, alpha=0.15)
            PyPlot.scatter([c.com_low * bins_to_deg], [mean_l_roi[round(Int, c.com_low)]], color=col, s=40, zorder=5)
        end
        PyPlot.legend(loc="upper right", fontsize=10)
        
        PyPlot.subplot(3, 1, 2)
        PyPlot.title("High Frequency (Integrated Profile & Barycentric Windows)", fontsize=13, fontweight="bold")
        PyPlot.plot(lon_arr, mean_h_roi, color="#ff7f0e", linewidth=2.0, label="High Data")
        PyPlot.fill_between(lon_arr, zeros(length(lon_arr)), mean_h_roi, color="#ff7f0e", alpha=0.15)
        PyPlot.ylabel("Normalized Intensity", fontsize=11)
        PyPlot.grid(true, linestyle=":", alpha=0.6)
        PyPlot.minorticks_on()
        for (k, c) in enumerate(ref_comps)
            col = colors[mod1(k, length(colors))]
            PyPlot.axvspan(c.left_bound * bins_to_deg, c.right_bound * bins_to_deg, color=col, alpha=0.15)
            PyPlot.scatter([c.com_high * bins_to_deg], [mean_h_roi[round(Int, c.com_high)]], color=col, s=40, zorder=5)
        end
        PyPlot.legend(loc="upper right", fontsize=10)
        
        PyPlot.subplot(3, 1, 3)
        PyPlot.title("Profile Comparison (Low vs High)", fontsize=13, fontweight="bold")
        PyPlot.plot(lon_arr, mean_l_roi, color="#1f77b4", linewidth=2.0, label="Low Frequency")
        PyPlot.plot(lon_arr, mean_h_roi, color="#ff7f0e", linestyle="--", linewidth=2.0, label="High Frequency")
        PyPlot.xlabel("Longitude [deg]", fontsize=12, fontweight="bold")
        PyPlot.ylabel("Intensity", fontsize=11)
        PyPlot.grid(true, linestyle=":", alpha=0.6)
        PyPlot.minorticks_on()
        PyPlot.legend(loc="upper right", fontsize=10)
        for (k, c) in enumerate(ref_comps)
            col = colors[mod1(k, length(colors))]
            offset_deg = c.offset * bins_to_deg
            PyPlot.text(c.com_low * bins_to_deg, 0.5, @sprintf("Δ=%.3f°", offset_deg), color=col, fontsize=12, ha="center", va="bottom", fontweight="bold")
        end
        PyPlot.tight_layout()
        PyPlot.savefig(joinpath(indir, "$(pulsar_name)_integrated_offsets_$(type).pdf"), dpi=150)
        PyPlot.close()
        
        # Save integrated fit statistics to CSV
        out_path = joinpath(indir, "$(pulsar_name)_integrated_offsets_$(type).csv")
        csv_data = Matrix{Any}(undef, n_comps + 1, 2)
        csv_data[1, 1] = "Global_Fourier"
        csv_data[1, 2] = global_offset_bins
        for k in 1:n_comps
            csv_data[k+1, 1] = "Component_$k"
            csv_data[k+1, 2] = ref_comps[k].offset
        end
        writedlm(out_path, csv_data, ',')

        # --- 4. LRFS P3 vs Phase Analysis ---
        println("\n--- P3 vs Phase Analysis ---")
        lrfs_file = joinpath(indir, "pulsar_low.debase.lrfs")
        if !isfile(lrfs_file)
            lrfs_file = joinpath(indir, "pulsar.debase.lrfs") 
        end

        p3_values = Float64[]
        long_values = Float64[]

        if isfile(lrfs_file)
            lrfs_data = load_ascii_all(lrfs_file) 
            n_freqs, n_lrfs_bins, n_pols = size(lrfs_data)
            freqs = range(0.0, 0.5, length=n_freqs)
            
            for k in 1:n_comps
                phi_avg = (ref_comps[k].com_low + ref_comps[k].com_high) / 2.0
                b_idx = clamp(round(Int, phi_avg), 1, n_lrfs_bins)
                
                spectrum = lrfs_data[:, b_idx, 1]
                valid_range = 3:n_freqs 
                if !isempty(valid_range)
                    max_i = argmax(spectrum[valid_range]) + (valid_range[1] - 1)
                    peak_freq = freqs[max_i]
                    if peak_freq > 0.01
                        push!(p3_values, 1.0 / peak_freq)
                        push!(long_values, phi_avg)
                    end
                end
            end
            
            if !isempty(p3_values)
                PyPlot.figure(figsize=(8, 5))
                PyPlot.title("Modulation Period (P3) vs Phase", fontsize=13, fontweight="bold")
                PyPlot.scatter(long_values .* bins_to_deg, p3_values, color="#e53935", s=50, edgecolor="black", linewidth=0.8, zorder=3)
                PyPlot.xlabel("Longitude [deg]", fontsize=12, fontweight="bold")
                PyPlot.ylabel("P3 [P0 units]", fontsize=12)
                PyPlot.grid(true, linestyle=":", alpha=0.6)
                PyPlot.minorticks_on()
                PyPlot.tight_layout()
                PyPlot.savefig(joinpath(indir, "$(pulsar_name)_fourier_p3_vs_phase_$(type).pdf"))
                PyPlot.close()
            end
        end

        # --- 5. Phase-Resolved P3 Analysis ---
        println("\n--- Starting Phase-Resolved Analysis (per P3 bin) ---")
        
        # Collect measurements in flat arrays for dynamic scattered plotting
        lon_all = Float64[]
        off_all = Float64[]
        err_all = Float64[]
        snr_all = Float64[]
        phase_all = Float64[]
        
        valid_phases = 0
        
        for phase_idx in 1:n_phases
            row_l = vec(l_matrix[phase_idx, :])
            row_h = vec(h_matrix[phase_idx, :])
            
            row_l_roi = copy(normalize_01(row_l))
            row_h_roi = copy(normalize_01(row_h))
            if bin_st > 1
                row_l_roi[1:Int(bin_st)-1] .= 0.0
                row_h_roi[1:Int(bin_st)-1] .= 0.0
            end
            if bin_end < n_bins
                row_l_roi[Int(bin_end)+1:end] .= 0.0
                row_h_roi[Int(bin_end)+1:end] .= 0.0
            end
            
            # Simple SNR metric: signal peak vs standard deviation of off-pulse region
            # We take a small sample outside the ROI to represent noise
            off_pulse_l = row_l[1:max(1, bin_st-10)] 
            noise_l = isempty(off_pulse_l) ? 1.0 : std(off_pulse_l)
            snr_val = (maximum(row_l) - mean(off_pulse_l)) / (noise_l + 1e-6)
            
            # Weryfikacja: Zachowujemy słabsze sygnały (próg 5%), ale odrzucamy fałszywy szum sprawdzając SNR > 3.0
            if maximum(row_l_roi) < 0.05 || maximum(row_h_roi) < 0.05 || snr_val < 3.0
                continue
            end
            
            valid_phases += 1
            
            # Dynamicznie znajdujemy komponenty obecne w tej konkretnej fazie:
            comps = ProfileMetrics.dynamic_component_fourier_offsets(row_l_roi, row_h_roi; threshold=0.05)
            
            for res in comps
                # Odrzucamy pomiary niefizyczne oraz te z drastycznie wielkim błędem (np. > 40 stopni)
                if isnan(res.offset) || abs(res.offset) > n_bins * 0.1 || res.error > (40.0 / bins_to_deg)
                    continue
                end
                
                com_low = res.com_ref
                com_high = res.com_ref + res.offset
                lon_deg = (com_low + com_high) / 2.0 * bins_to_deg
                
                offset_deg = res.offset * bins_to_deg
                # Błąd jest teraz rygorystycznie liczony ze stromości wariancji faz Fouriera!
                error_deg = res.error * bins_to_deg
                
                push!(lon_all, lon_deg)
                push!(off_all, offset_deg)
                push!(err_all, error_deg)
                push!(snr_all, snr_val)
                push!(phase_all, Float64(phase_idx))
            end
        end
        
        println(">>> Processed $valid_phases valid phases out of $n_phases")
        
        # --- 6. Plotting: Phase-Resolved Offsets ---
        # Ulepszono wizualizację, aby była bardziej czytelna i estetyczna.
        # Dodano rygorystyczne naukowo obliczenie 95% przedziału ufności dla linii trendu
        # przy użyciu metody bootstrap, co pozwala na ocenę niepewności dopasowania.
        # Wprowadzono statystyczne ważenie punktów (Inverse Variance) oraz słupki błędów (Error Bars).
        # Zwiększono wysokość figury, aby pomieścić dwa panele (główny i zbliżenie).
        PyPlot.figure(figsize=(12, 11))
        
        if !isempty(lon_all)
            # Funkcja pomocnicza do ważonego wygładzania jądrowego (Nadaraya-Watson z Inverse-Variance Weighting)
            function kernel_smooth_weighted(x_data, y_data, err_data, x_grid, sigma)
                y_grid = Float64[]
                for x in x_grid
                    # Waga odległości (Jądro Gaussa)
                    w_kernel = exp.(-0.5 .* ((x_data .- x) ./ sigma).^2)
                    # Waga ufności statystycznej (odwrotność wariancji)
                    w_stat = 1.0 ./ (err_data.^2 .+ 1e-6)
                    weights = w_kernel .* w_stat
                    
                    weight_sum = sum(weights)
                    if weight_sum > 1e-6
                        push!(y_grid, sum(weights .* y_data) / weight_sum)
                    else
                        push!(y_grid, NaN)
                    end
                end
                return y_grid
            end

            # Sortowanie danych wejściowych jest kluczowe dla wygładzania
            perm = sortperm(lon_all)
            slon = lon_all[perm]
            soff = off_all[perm]
            serr = err_all[perm]
            
            # Definicja siatki i parametru wygładzania (sigma)
            smooth_lon = range(minimum(slon), maximum(slon), length=200)
            sigma = (maximum(slon) - minimum(slon)) / 15.0 # Szersze jądro dla gładszego trendu

            # Obliczenie głównego trendu (z rygorystycznym uwzględnieniem błędów)
            smooth_off = kernel_smooth_weighted(slon, soff, serr, smooth_lon, sigma)

            # Pętla bootstrap do estymacji przedziału ufności
            n_bootstrap = 200 # Więcej próbek dla stabilniejszych wyników
            bootstrap_curves = zeros(n_bootstrap, length(smooth_lon))
            
            for i in 1:n_bootstrap
                indices = rand(1:length(lon_all), length(lon_all))
                # Próbkowanie z powtórzeniami
                boot_lon = lon_all[indices]
                boot_off = off_all[indices]
                boot_err = err_all[indices]
                
                perm_boot = sortperm(boot_lon)
                bootstrap_curves[i, :] = kernel_smooth_weighted(boot_lon[perm_boot], boot_off[perm_boot], boot_err[perm_boot], smooth_lon, sigma)
            end

            # Obliczanie kwantyli dla 95% przedziału ufności
            lower_ci, upper_ci = Float64[], Float64[]
            for j in 1:length(smooth_lon)
                valid_vals = filter(!isnan, bootstrap_curves[:, j])
                if length(valid_vals) > 20 # Wymagamy minimum punktów do wiarygodnego kwantyla
                    push!(lower_ci, quantile(valid_vals, 0.025))
                    push!(upper_ci, quantile(valid_vals, 0.975))
                else
                    push!(lower_ci, NaN); push!(upper_ci, NaN)
                end
            end

            # --- Wykres 1: Główny rozkład z punktami (Górny Panel) ---
            ax1 = PyPlot.subplot(2, 1, 1)

            # Rysowanie rzeczywistych, formalnych błędów matematycznych (Error bars)
            PyPlot.errorbar(lon_all, off_all, yerr=err_all, fmt="none", ecolor="black", elinewidth=0.8, capsize=1.5, alpha=0.4, zorder=2)

            # Wykres punktowy (Scatter) z ulepszoną stylistyką dla lepszej widoczności
            # Rozmiar punktu skalujemy z SNR, by wyraźne sygnały były łatwiejsze do interpretacji wizualnej
            point_sizes = clamp.(snr_all .* 2.5, 20.0, 90.0)
            sc = PyPlot.scatter(lon_all, off_all, c=phase_all, cmap="viridis", alpha=0.85, s=point_sizes, edgecolor="black", linewidth=0.6, zorder=3)
            cbar = PyPlot.colorbar(sc)
            cbar.set_label("P3 Phase (Row Index)", fontsize=12, fontweight="bold")
            
            # Rysowanie przedziału ufności jako cieniowany obszar
            PyPlot.fill_between(smooth_lon, lower_ci, upper_ci, color="#d62728", alpha=0.3, label="95% Bootstrap Confidence Interval", zorder=4)
            
            # Rysowanie wygładzonego trendu
            PyPlot.plot(smooth_lon, smooth_off, color="#d62728", lw=3.2, label="Weighted Trend (Inv. Variance)", zorder=5)
        
            PyPlot.axhline(0.0, color="gray", ls="--", lw=1.5, alpha=0.8, zorder=2)
            PyPlot.minorticks_on()
            PyPlot.grid(true, which="major", linestyle="--", color="gray", alpha=0.6)
            PyPlot.grid(true, which="minor", linestyle=":", color="gray", alpha=0.3)
            
            PyPlot.ylabel("Offset [deg]", fontsize=13, fontweight="bold")
            PyPlot.title("Phase-Resolved Offset vs. Longitude with 95% Confidence Interval", fontsize=15, fontweight="bold")
            PyPlot.legend(loc="best", fontsize=11, frameon=true, framealpha=0.95, edgecolor="black")

            # --- Wykres 2: Zbliżenie na sam wygładzony trend (Dolny Panel) ---
            ax2 = PyPlot.subplot(2, 1, 2, sharex=ax1)
            
            # Rysujemy wyłącznie przedział ufności i trend
            PyPlot.fill_between(smooth_lon, lower_ci, upper_ci, color="#d62728", alpha=0.3, label="95% Bootstrap Confidence Interval", zorder=4)
            PyPlot.plot(smooth_lon, smooth_off, color="#d62728", lw=3.2, label="Weighted Trend (Inv. Variance)", zorder=5)
            
            PyPlot.axhline(0.0, color="gray", ls="--", lw=1.5, alpha=0.8, zorder=2)
            PyPlot.minorticks_on()
            PyPlot.grid(true, which="major", linestyle="--", color="gray", alpha=0.6)
            PyPlot.grid(true, which="minor", linestyle=":", color="gray", alpha=0.3)
            
            PyPlot.xlabel("Longitude [deg]", fontsize=13, fontweight="bold")
            PyPlot.ylabel("Trend Offset [deg]", fontsize=13, fontweight="bold")
            PyPlot.title("Zoom-in: Fitted Trend Variations", fontsize=14, fontweight="bold")
            
            # Automatyczny dobór ciasnego zakresu Y dla zbliżenia, bazujący na bandzie błędu
            valid_lower = filter(!isnan, lower_ci)
            valid_upper = filter(!isnan, upper_ci)
            if !isempty(valid_lower) && !isempty(valid_upper)
                y_min, y_max = minimum(valid_lower), maximum(valid_upper)
                margin = max((y_max - y_min) * 0.25, 0.05) # Minimum 0.05 stopnia marginesu
                PyPlot.ylim(y_min - margin, y_max + margin)
            end
        else
            # Ubezpieczenie na wypadek braku punktów w danej fazie
            PyPlot.title("Phase-Resolved Offset vs. Longitude (No Data)", fontsize=15, fontweight="bold")
            PyPlot.xlabel("Longitude [deg]", fontsize=13, fontweight="bold")
            PyPlot.ylabel("Offset [deg]", fontsize=13, fontweight="bold")
        end
        
        PyPlot.tight_layout(pad=1.5)
        
        out_name = "$(pulsar_name)_fourier_offsets_$(type).pdf"
        PyPlot.savefig(joinpath(indir, out_name), dpi=150)
        PyPlot.close()
        
        # --- 7. Output to CSV ---
        csv_path = joinpath(indir, "$(pulsar_name)_fourier_offsets_$(type).csv")
        open(csv_path, "w") do io
            println(io, "Longitude_deg,Offset_deg,Error_deg,Phase_SNR,P3_Phase_Idx")
            for i in 1:length(lon_all)
                println(io, "$(lon_all[i]),$(off_all[i]),$(err_all[i]),$(snr_all[i]),$(phase_all[i])")
            end
        end

        # --- 8. Component Separation (P2) Calculation ---
        # Since we use continuous tracking, static P2 is less defined here.
        # Keeping previous Barycenter-based P2 from the integrated profile:
        println("\n--- Integrated Profile Separation (P2) ---")
        for k in 1:(length(ref_comps) - 1)
            lon1 = ref_comps[k].com_low * bins_to_deg
            lon2 = ref_comps[k+1].com_low * bins_to_deg
            p2_sep = abs(lon2 - lon1)
            @printf(">>> Base P2 between C%d and C%d: %.3f°\n", k, k+1, p2_sep)
        end
        
        println("="^60 * "\n")
        return ref_comps
    end

    # Pomocnicza normalizacja wewnątrz modułu
    function normalize_02(vec)
        m, M = minimum(vec), maximum(vec)
        return (vec .- m) ./ (M - m)
    end

    


    function normalize_01(data)
        min_val = minimum(data)
        max_val = maximum(data)
        return (data .- min_val) ./ (max_val - min_val)
    end

    function normalize_per_pulse(data)
        normalized = similar(data)
        for i in 1:size(data, 1)
            pulse = data[i, :]
            min_val = minimum(pulse)
            max_val = maximum(pulse)
            normalized[i, :] = (pulse .- min_val) ./ (max_val - min_val)
        end
        return normalized
    end


    function find_shift_per_pulse(data1, data2)
        n_pulses = size(data1, 1)
        shifts = zeros(Int, n_pulses)
        correlations = zeros(n_pulses)
        
        for i in 1:n_pulses
            # Cross-correlation
            ccf = ifft(fft(data1[i, :]) .* conj.(fft(data2[i, :])))
            ccf = real.(ccf)
            
            # Znajdź maksimum
            max_idx = argmax(ccf)
            correlations[i] = ccf[max_idx]
            
            # Przelicz na przesunięcie (z uwzględnieniem wraparound)
            n_bins = length(ccf)
            shifts[i] = max_idx - 1
            if shifts[i] > n_bins ÷ 2
                shifts[i] = shifts[i] - n_bins
            end
        end
        
        return shifts, correlations
    end

    """
    Process multifrequency data with PSRCHIVE and PSRSALSA single pulse analysis
    """
    function process_psrdata_single(indir, outdir; files=nothing, outfile="pulsar.spCf16", params_file="params.json")

        # full paths
        params_file = joinpath(outdir, params_file)
        outfile = joinpath(outdir, outfile)

        # if no catalog
        if !isdir(outdir)
            mkdir(outdir)
        end

        # check if params_file exists if not creating default one
        if !isfile(params_file)
            println("File $params_file does not exist, creating default one.")
            p = Tools.default_params(params_file)
        else
            p = Tools.read_params(params_file)
        end

        # speeding up
        # add all .spCf16 files 
        add_psrfiles(indir, outfile; files=files, sixteen=true)

        # divide to two frequencies
        low_filename, high_filename, mid_filename = multifrequency_split(outfile)

        # gets number of single pulses if needed
        if isnothing(p["nsubint"])
            get_nsubint(low_filename, params_file, p)
        end

        # debase the data
        debase_16(low_filename, high_filename, mid_filename, params_file, p)

        # speeding ends here

        # comment those three if speeding up is off
        #low_filename = "/home/psr/output/J2139+2242_single/pulsar.low"
        #high_filename = "/home/psr/output/J2139+2242_single/pulsar.high"
        #mid_filename = "/home/psr/output/J2139+2242_single/pulsar.mid"

        low_debase = replace(low_filename, ".low"=>"_low.debase.gg")
        Tools.generate_snr(low_debase) 
        low_snrfile = low_debase * ".snr.txt"

        twodfs_lrfs(low_debase, params_file, p; detect=false) # just to detect p2

        # single pulses low
        da_lo = Data.load_ascii(replace(low_filename, ".low"=>"_low_debase.txt"))

        high_debase = replace(high_filename, ".high"=>"_high.debase.gg")
        Tools.generate_snr(high_debase) 
        high_snrfile = high_debase * ".snr.txt"
 
        # single pulses high 
        da_hi = Data.load_ascii(replace(high_filename, ".high"=>"_high_debase.txt"))

        # single pulses mid 
        da_mi = Data.load_ascii(replace(mid_filename, ".mid"=>"_mid_debase.txt"))

        Tools.track_subpulses_snr3(da_lo, da_hi, da_mi, 10, low_snrfile, on_st=p["bin_st"], on_end=p["bin_end"])
        
        Plot.analyse_single(da_lo, da_hi, da_mi, p; pulse_start=250)

        # TODO single pulse detection work starts here
        #Tools.track_subpulses_manual(da_mi, on_st=p["bin_st"], on_end=p["bin_end"])

        return 
    end





end # module