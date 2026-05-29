module Data
    using Glob
    using FITSIO
    using ProgressMeter
    using CairoMakie, FileIO
    using PDFIO
    using Statistics
    using FFTW
    using DelimitedFiles
    using PyCall

    include("functions.jl")
    include("tools.jl")
    include("plot.jl")


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

        for i in eachindex(lines)
            startswith(lines[i], "F") && continue
            startswith(lines[i], "M") && continue
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

        if !isnothing(p["zaps"])
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
        
        # 3 channels  
        run(pipeline(`paz -Z 0-12 -e high $spCf16_file`,stderr="errs.txt"))
        run(pipeline(`paz -Z 3-15 -e low $spCf16_file`,stderr="errs.txt"))
        run(pipeline(`paz -Z 0-6 -Z 9-15 -e mid $spCf16_file`,stderr="errs.txt")) # 2 channels?

        # 6 channels  
        #run(pipeline(`paz -Z 0-9 -e high $spCf16_file`,stderr="errs.txt"))
        #run(pipeline(`paz -Z 6-15 -e low $spCf16_file`,stderr="errs.txt"))
        #run(pipeline(`paz -Z 0-4 -Z 11-15 -e mid $spCf16_file`,stderr="errs.txt"))

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

        profile_idx = 1
        p1 = nl[profile_idx, :]
        p2 = nh[profile_idx, :]
        bins = 1:length(p1)

        writedlm(joinpath(indir,"p_low.txt"), hcat(bins, p1), ' ')
        writedlm(joinpath(indir,"p_high.txt"), hcat(bins, p2), ' ')


        println(size(nl))

        Plot.analyse_p3folds(nl, nh, p)
        Plot.analyse_p3folds2(nl, nh, p)

    end


    """
    Analyse p3folds

    # Arguments
    - `indir`: data input directory
    - `type`: type of p3folds "refine" or "norefine" filenames are based on that
 
    """
    function analyse_p3folds_16_new(indir, type; n_comp=3)

        # parameters file 
        p = Tools.read_params(joinpath(indir, "params.json"))

        # low, high frequancy filename end
        low = joinpath(indir, "pulsar_low.debase.p3fold_" * type)
        high = joinpath(indir, "pulsar_high.debase.p3fold_" * type)

        #println(low)
        l = Data.load_ascii(low)
        #Plot.p3fold(l, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low_$(type)_analyse", show_=true, repeat_num=1)

        #println(high)
        h = Data.load_ascii(high)
        #Plot.p3fold(h, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_high_$(type)_analyse", show_=true, repeat_num=1)
 
        #diff = normalize_01(l) .- normalize_01(h)
        diff = normalize_per_pulse(l) .- normalize_per_pulse(h)
        #Plot.p3fold(diff, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low_high_$(type)_analyse", show_=true, repeat_num=1)

        nl = normalize_per_pulse(l)
        nh = normalize_per_pulse(h)
        #=
        profile_idx = 1
        p1 = nl[profile_idx, :]
        p2 = nh[profile_idx, :]
        bins = 1:length(p1)

        writedlm(joinpath(indir,"p_low.txt"), hcat(bins, p1), ' ')
        writedlm(joinpath(indir,"p_high.txt"), hcat(bins, p2), ' ')
        =#
        # TODO work here...

        Plot.analyse_p3folds3(nl, nh, p, n_comp)


        #println(size(nl))


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

    function p3fold_psrdata(indir)

        # parameters file 
        p = Tools.read_params(joinpath(indir, "params.json"))

        # low, high frequancy filenames
        low = joinpath(indir, "pulsar.low")
        high = joinpath(indir, "pulsar.high")

        low_debase = replace(low, ".low"=>"_low.debase.gg")
        high_debase = replace(high, ".high"=>"_high.debase.gg")
        
        #-p3fold_norefine -p3fold_smooth 3 -p3fold_nritt 1 -p3fold_cpb 1

        
        # PSRSALSA p3folding low freq
        run(pipeline(`pfold -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $low_debase`,  stderr="errs.txt"))
        mv(replace(low_debase, ".gg"=>".p3fold"), replace(low_debase, ".gg"=>".p3fold_refine"), force=true)
        d5 = Data.load_ascii(replace(low_debase, ".gg"=>".p3fold_refine"))
        Plot.p3fold(d5, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low", show_=true, repeat_num=4)
        # norefine
        #println("pfold -p3fold_norefine -p3fold \"$(p["p3"]) $(p["p3_ybins"])\" -onpulse \"$(p["bin_st"]) $(p["bin_end"])\" -onpulsed \"/NULL\" -p3foldd \"/NULL\" -w -oformat ascii $low_debase")
        run(pipeline(`pfold -p3fold_norefine -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $low_debase`,  stderr="errs.txt"))
        mv(replace(low_debase, ".gg"=>".p3fold"), replace(low_debase, ".gg"=>".p3fold_norefine"), force=true)
        d5 = Data.load_ascii(replace(low_debase, ".gg"=>".p3fold_norefine"))
        Plot.p3fold(d5, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low_norefine", show_=true, repeat_num=4)
 
        # PSRSALSA p3folding high freq
        run(pipeline(`pfold -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $high_debase`,  stderr="errs.txt"))
        mv(replace(high_debase, ".gg"=>".p3fold"), replace(high_debase, ".gg"=>".p3fold_refine"), force=true)
        d5 = Data.load_ascii(replace(high_debase, ".gg"=>".p3fold_refine"))
        Plot.p3fold(d5, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_high", show_=true, repeat_num=4)
        # norefine
        run(pipeline(`pfold -p3fold_norefine -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $high_debase`,  stderr="errs.txt"))
        mv(replace(high_debase, ".gg"=>".p3fold"), replace(high_debase, ".gg"=>".p3fold_norefine"), force=true)
        d5 = Data.load_ascii(replace(high_debase, ".gg"=>".p3fold_norefine"))
        Plot.p3fold(d5, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_high_norefine", show_=true, repeat_num=4)

        #=
        =#




    end



    """
    Remove orthogonal polarization mode (OPM) jumps from a PA track by
    unwrapping along longitude with a 90° period. Walks bin-by-bin keeping
    a running "continuous" PA; whenever the next sample lies farther than
    `tolerance` degrees from it, shift the sample by ±90° (the OPM step)
    until it rejoins the track. The accumulated offset persists, so a whole
    OPM island is corrected end-to-end regardless of length. Finally the
    output is wrapped back to (-90°, 90°] — which also undoes any shifts
    of ±180° that correspond to the intrinsic mod-π PA wrap, leaving those
    pixels unchanged.
    """
    function deopm_pa(pa_deg; tolerance=45.0)
        pa_out = copy(pa_deg)
        n = length(pa_out)
        last_unwrapped = NaN
        offset = 0.0
        flipped = 0
        for i in 1:n
            if isnan(pa_out[i])
                continue  # keep offset/last_unwrapped across NaN gaps
            end
            pa_in = pa_out[i]
            if isnan(last_unwrapped)
                last_unwrapped = pa_in
                continue
            end
            shifted = pa_in + offset
            d = shifted - last_unwrapped
            while d > tolerance
                offset -= 90; shifted -= 90; d -= 90
            end
            while d < -tolerance
                offset += 90; shifted += 90; d += 90
            end
            pa_out[i] = mod(shifted + 90, 180) - 90
            abs(pa_out[i] - pa_in) > 0.1 && (flipped += 1)
            last_unwrapped = shifted
        end
        return pa_out, flipped
    end


    """
    Fit the Rotating Vector Model (RVM) to position angle data.

    Faithful Johnston+2023 method:
    - 85 trials α ∈ (5°,175°), 40 trials β ∈ (±20°)  [paper values]
    - scipy.optimize.least_squares for (φ₀, PA₀) at each (α,β) pair
    - OPM built into residuals: min(|resA|, |resB|) for each data point
    - Fast: vectorised numpy loop over φ₀ only; scipy for the global best
    - α normalised to < 90° in reported values (paper convention)
    - Minimum 10 valid PA points required (paper criterion)

    Returns NamedTuple with fields: alpha, beta, phi0, pa0 (all in degrees),
    chi2, chi2_red, ndof.  Returns nothing if < 10 valid PA points.
    """
    const _rvm_py_ready = Ref(false)

    function _ensure_rvm_python()
        _rvm_py_ready[] && return
        py"""
import numpy as _np
from scipy.optimize import least_squares as _ls

def _rvm_res(params, lon, pa, sigma, alpha, beta):
    phi0, pa0 = params
    zeta = alpha + beta
    ca, sa = _np.cos(alpha), _np.sin(alpha)
    cz, sz = _np.cos(zeta),  _np.sin(zeta)
    dphi = lon - phi0
    rvm  = pa0 + _np.arctan2(-sa * _np.sin(dphi),
                               sz*ca - cz*sa*_np.cos(dphi))
    diff = pa - rvm
    # OPM: for each point try both modes, take the nearer one
    resA = _np.mod(diff + _np.pi/2, _np.pi) - _np.pi/2
    resB = _np.mod(diff + _np.pi,   _np.pi) - _np.pi/2
    return _np.where(_np.abs(resA) <= _np.abs(resB), resA, resB) / sigma

def _scipy_ab(lon, pa, sig, alpha, beta, phi0_init, pa0_init):
    # scipy from grid best; also tries PA0 +/- pi/2 to overcome OPM-biased start
    best_c2 = _np.inf
    best_x  = [phi0_init, pa0_init]
    for pa0_try in [pa0_init, pa0_init + _np.pi/2, pa0_init - _np.pi/2]:
        try:
            r = _ls(_rvm_res, [phi0_init, pa0_try],
                    args=(lon, pa, sig, alpha, beta), method='lm')
            c2 = float(_np.sum(r.fun**2))
            if c2 < best_c2:
                best_c2 = c2
                best_x  = [float(r.x[0]), float(r.x[1])]
        except Exception:
            pass
    return best_c2, best_x[0], best_x[1]

def _fit_rvm_py(lon_deg, pa_deg, pa_err_deg,
                n_alpha, n_beta, n_phi0, n_scipy_top,
                alpha_min, alpha_max, beta_min, beta_max):
    lon_deg    = _np.asarray(lon_deg,    dtype=float)
    pa_deg     = _np.asarray(pa_deg,     dtype=float)
    pa_err_deg = _np.asarray(pa_err_deg, dtype=float)
    mask = _np.isfinite(pa_deg) & _np.isfinite(pa_err_deg) & (pa_err_deg > 0)
    if int(mask.sum()) < 10:
        return None

    lon = _np.deg2rad(lon_deg[mask])
    pa  = _np.deg2rad(pa_deg[mask])
    sig = _np.deg2rad(pa_err_deg[mask])
    inv_var = 1.0 / sig**2

    na  = int(n_alpha)
    nb  = int(n_beta)
    np_ = int(n_phi0)
    ntop = int(n_scipy_top)

    alphas = _np.linspace(_np.deg2rad(alpha_min), _np.deg2rad(alpha_max), na)
    betas  = _np.linspace(_np.deg2rad(beta_min),  _np.deg2rad(beta_max),  nb)
    lon_min, lon_max = float(lon.min()), float(lon.max())
    margin = (lon_max - lon_min) * 0.5
    phi0s  = _np.linspace(lon_min - margin, lon_max + margin, np_)

    # --- Vectorised grid: loop over phi0, vectorise over (alpha, beta, pts) ---
    a2 = alphas[:, None];  b2 = betas[None, :]
    z2 = a2 + b2
    sa2 = _np.sin(a2);  ca2 = _np.cos(a2)
    sz2 = _np.sin(z2);  cz2 = _np.cos(z2)

    chi2_map = _np.full((na, nb), _np.inf)
    phi0_map = _np.full((na, nb), phi0s[0])
    pa0_map  = _np.zeros((na, nb))

    iv3 = inv_var[None, None, :]
    pa3 = pa     [None, None, :]

    for phi0 in phi0s:
        dphi = lon - phi0
        sd = _np.sin(dphi);  cd = _np.cos(dphi)
        num = -sa2[:, :, None] * sd[None, None, :]
        den = (sz2[:, :, None] * ca2[:, :, None]
               - cz2[:, :, None] * sa2[:, :, None] * cd[None, None, :])
        rvm   = _np.arctan2(num, den)
        diffs = pa3 - rvm
        Sx  = _np.sum(_np.sin(2 * diffs) * iv3, axis=2)
        Cx  = _np.sum(_np.cos(2 * diffs) * iv3, axis=2)
        pa0 = _np.arctan2(Sx, Cx) / 2
        d2   = diffs - pa0[:, :, None]
        resA = _np.mod(d2 + _np.pi/2, _np.pi) - _np.pi/2
        resB = _np.mod(d2 + _np.pi,   _np.pi) - _np.pi/2
        res  = _np.where(_np.abs(resA) <= _np.abs(resB), resA, resB)
        chi2 = _np.sum(res**2 * iv3, axis=2)
        improved = chi2 < chi2_map
        chi2_map = _np.where(improved, chi2,  chi2_map)
        phi0_map = _np.where(improved, phi0,  phi0_map)
        pa0_map  = _np.where(improved, pa0,   pa0_map)

    # --- scipy refinement for top-N (alpha,beta) pairs (Johnston+2023 method) ---
    # Running scipy for top-N handles OPM-biased analytical PA0: even if the grid
    # ranks a pair sub-optimally due to OPM, scipy with ±PA0 tries fixes it.
    flat_order = _np.argsort(chi2_map, axis=None)   # ascending chi2
    best_c2   = _np.inf
    best_ia_s = int(_np.unravel_index(flat_order[0], chi2_map.shape)[0])
    best_ib_s = int(_np.unravel_index(flat_order[0], chi2_map.shape)[1])
    best_phi0 = float(phi0_map[best_ia_s, best_ib_s])
    best_pa0  = float(pa0_map [best_ia_s, best_ib_s])

    for k in range(min(ntop, na * nb)):
        ia, ib = _np.unravel_index(flat_order[k], chi2_map.shape)
        ia, ib = int(ia), int(ib)
        c2, p0, pa0_s = _scipy_ab(lon, pa, sig,
                                   float(alphas[ia]), float(betas[ib]),
                                   float(phi0_map[ia, ib]),
                                   float(pa0_map [ia, ib]))
        chi2_map[ia, ib] = c2          # update map with scipy-refined value
        if c2 < best_c2:
            best_c2   = c2
            best_ia_s = ia
            best_ib_s = ib
            best_phi0 = p0
            best_pa0  = pa0_s

    alpha_b = float(alphas[best_ia_s])
    beta_b  = float(betas [best_ib_s])
    ndof    = max(int(mask.sum()) - 4, 1)

    # Wrap PA0 to (−90°, 90°]
    pa0_deg = float(_np.mod(_np.rad2deg(best_pa0) + 90.0, 180.0) - 90.0)

    # Normalised alpha (< 90°) — used for printing/tables only (Johnston+2023 convention)
    a_norm = float(_np.rad2deg(alpha_b))
    b_norm = float(_np.rad2deg(beta_b))
    if a_norm > 90.0:
        a_norm =  180.0 - a_norm
        b_norm = -b_norm

    return {
        'alpha':      float(_np.rad2deg(alpha_b)),   # raw best-fit (for display)
        'beta':       float(_np.rad2deg(beta_b)),
        'alpha_norm': a_norm,                         # normalised (for tables)
        'beta_norm':  b_norm,
        'phi0':       float(_np.rad2deg(best_phi0)),
        'pa0':        pa0_deg,
        'chi2':       float(best_c2),
        'chi2_red':   float(best_c2 / ndof),
        'ndof':       int(ndof),
        'chi2_map':   chi2_map.T,
        'alphas':     _np.rad2deg(alphas),
        'betas':      _np.rad2deg(betas)
    }
"""
        _rvm_py_ready[] = true
    end

    function fit_rvm(lon_deg, pa_deg, pa_err_deg;
                     return_map=false,
                     n_alpha=85, n_beta=40, n_phi0=120, n_scipy_top=20,
                     alpha_range=(5.0, 175.0), beta_range=(-20.0, 20.0))

        _ensure_rvm_python()

        result_py = py"_fit_rvm_py"(
            collect(lon_deg), collect(pa_deg), collect(pa_err_deg),
            n_alpha, n_beta, n_phi0, n_scipy_top,
            alpha_range[1], alpha_range[2],
            beta_range[1],  beta_range[2]
        )

        (result_py === nothing || result_py == py"None") && return nothing

        # alpha/beta are raw best-fit values (used for display & chi2 map marker).
        # alpha_norm/beta_norm are normalised to α<90° (Johnston+2023 table convention).
        result = (alpha      = Float64(result_py["alpha"]),
                  beta       = Float64(result_py["beta"]),
                  alpha_norm = Float64(result_py["alpha_norm"]),
                  beta_norm  = Float64(result_py["beta_norm"]),
                  phi0       = Float64(result_py["phi0"]),
                  pa0        = Float64(result_py["pa0"]),
                  chi2       = Float64(result_py["chi2"]),
                  chi2_red   = Float64(result_py["chi2_red"]),
                  ndof       = Int(result_py["ndof"]))

        if return_map
            chi2_map_j = convert(Matrix{Float64}, result_py["chi2_map"])
            alphas_j   = Float64.(Array(result_py["alphas"]))
            betas_j    = Float64.(Array(result_py["betas"]))
            return result, chi2_map_j, alphas_j, betas_j
        else
            return result
        end
    end


    """
    Evaluate RVM curve over a dense longitude grid (degrees in, degrees out).
    Returns (lon_dense, pa_rvm, pa_rvm_ortho) where pa_rvm_ortho is the 90° mode.
    """
    function rvm_curve(params, lon_min_deg, lon_max_deg; npts=500)
        lon = collect(range(lon_min_deg, lon_max_deg, length=npts))
        lon_r = lon .* (π / 180)
        phi0  = params.phi0 * (π / 180)
        pa0   = params.pa0  * (π / 180)
        alpha = params.alpha * (π / 180)
        beta  = params.beta  * (π / 180)
        zeta  = alpha + beta

        dphi = lon_r .- phi0
        # Komesaroff (1970) sign convention (matches publication).
        pa_r = pa0 .+ atan.(.-sin(alpha) .* sin.(dphi),
                             sin(zeta) .* cos(alpha) .- cos(zeta) .* sin(alpha) .* cos.(dphi))
        pa_deg       = mod.(pa_r .* (180/π) .+ 90, 180) .- 90
        pa_ortho_deg = mod.(pa_deg .+ 180, 180) .- 90

        # Break line at wrap jumps so renderers don't draw vertical segments
        idx = eachindex(pa_deg)
        for i in Iterators.drop(idx, 1)
            if abs(pa_deg[i] - pa_deg[i-1]) > 90
                pa_deg[i-1] = NaN
            end
            if abs(pa_ortho_deg[i] - pa_ortho_deg[i-1]) > 90
                pa_ortho_deg[i-1] = NaN
            end
        end

        return lon, pa_deg, pa_ortho_deg
    end


    function position_angle(indir)
        # parameters file
        p = Tools.read_params(joinpath(indir, "params.json"))

        # low, high frequancy filenames
        low = joinpath(indir, "pulsar.low")
        high = joinpath(indir, "pulsar.high")

        # txt files
        lt = joinpath(indir, "pulsar_low.txt") # no dabese here needed?
        ht = joinpath(indir, "pulsar_high.txt") # no dabese here needed?

        l = Data.load_ascii_all(lt)
        h = Data.load_ascii_all(ht)

        bin_st  = p["bin_st"]
        bin_end = p["bin_end"]

        pulses_l, bins_l, _ = size(l)
        pulses_h, bins_h, _ = size(h)

        # Longitude axis (centered at 0)
        db_l = (bin_end + 1) - bin_st
        dl_l = 360.0 * db_l / bins_l
        lon_l = collect(range(-dl_l/2.0, dl_l/2.0, length=db_l))

        db_h = (bin_end + 1) - bin_st
        dl_h = 360.0 * db_h / bins_h
        lon_h = collect(range(-dl_h/2.0, dl_h/2.0, length=db_h))

        # Average Stokes profiles over on-pulse window
        I_l = vec(mean(l[:, bin_st:bin_end, 1], dims=1))
        Q_l = vec(mean(l[:, bin_st:bin_end, 2], dims=1))
        U_l = vec(mean(l[:, bin_st:bin_end, 3], dims=1))
        V_l = vec(mean(l[:, bin_st:bin_end, 4], dims=1))
        Lin_l = sqrt.(Q_l.^2 .+ U_l.^2)

        I_h = vec(mean(h[:, bin_st:bin_end, 1], dims=1))
        Q_h = vec(mean(h[:, bin_st:bin_end, 2], dims=1))
        U_h = vec(mean(h[:, bin_st:bin_end, 3], dims=1))
        V_h = vec(mean(h[:, bin_st:bin_end, 4], dims=1))
        Lin_h = sqrt.(Q_h.^2 .+ U_h.^2)

        # PA masked where L < 5σ on averaged Q, U (Johnston et al. 2023)
        off_noise_l = std(l[:, 1:bin_st-1, 1])
        off_noise_h = std(h[:, 1:bin_st-1, 1])
        sigma_avg_l = off_noise_l / sqrt(pulses_l)
        sigma_avg_h = off_noise_h / sqrt(pulses_h)
        thresh_l = 5.0 * sigma_avg_l
        thresh_h = 5.0 * sigma_avg_h

        pa_l = [Lin_l[i] > thresh_l ? 0.5 * atan(U_l[i], Q_l[i]) * (180.0/pi) : NaN for i in 1:db_l]
        pa_h = [Lin_h[i] > thresh_h ? 0.5 * atan(U_h[i], Q_h[i]) * (180.0/pi) : NaN for i in 1:db_h]

        # PA errors: σ_PA = 0.5 * σ_noise_avg / L  (in degrees)
        pa_err_l = [Lin_l[i] > thresh_l ? 0.5 * sigma_avg_l / Lin_l[i] * (180.0/pi) : NaN for i in 1:db_l]
        pa_err_h = [Lin_h[i] > thresh_h ? 0.5 * sigma_avg_h / Lin_h[i] * (180.0/pi) : NaN for i in 1:db_h]

        # Detect and undo orthogonal polarization mode jumps before RVM fit
        pa_l, flipped_l = deopm_pa(pa_l)
        pa_h, flipped_h = deopm_pa(pa_h)
        println("De-OPM: flipped $flipped_l (low) and $flipped_h (high) bins")

        # Select OPM branch (publication convention): global +90° shift
        pa_shift = 90.0
        pa_l = [isnan(x) ? x : mod(x + pa_shift + 90, 180) - 90 for x in pa_l]
        pa_h = [isnan(x) ? x : mod(x + pa_shift + 90, 180) - 90 for x in pa_h]

        # Fit RVM to low-frequency PA (more points typically)
        println("Fitting RVM (low frequency)...")
        rvm_params_l = fit_rvm(lon_l, pa_l, pa_err_l)
        if !isnothing(rvm_params_l)
            println("  α = $(round(rvm_params_l.alpha, digits=1))°, " *
                    "β = $(round(rvm_params_l.beta, digits=1))°, " *
                    "φ₀ = $(round(rvm_params_l.phi0, digits=2))°, " *
                    "PA₀ = $(round(rvm_params_l.pa0, digits=1))°, " *
                    "χ²/ndof = $(round(rvm_params_l.chi2_red, digits=2))")
            lon_rvm_l, pa_rvm_l, pa_rvm_l_ortho = rvm_curve(rvm_params_l,
                                                              minimum(lon_l), maximum(lon_l))
        else
            println("  Not enough PA points for RVM fit (low)")
            lon_rvm_l = pa_rvm_l = pa_rvm_l_ortho = nothing
        end

        println("Fitting RVM (high frequency)...")
        rvm_params_h = fit_rvm(lon_h, pa_h, pa_err_h)
        if !isnothing(rvm_params_h)
            println("  α = $(round(rvm_params_h.alpha, digits=1))°, " *
                    "β = $(round(rvm_params_h.beta, digits=1))°, " *
                    "φ₀ = $(round(rvm_params_h.phi0, digits=2))°, " *
                    "PA₀ = $(round(rvm_params_h.pa0, digits=1))°, " *
                    "χ²/ndof = $(round(rvm_params_h.chi2_red, digits=2))")
            lon_rvm_h, pa_rvm_h, pa_rvm_h_ortho = rvm_curve(rvm_params_h,
                                                              minimum(lon_h), maximum(lon_h))
        else
            println("  Not enough PA points for RVM fit (high)")
            lon_rvm_h = pa_rvm_h = pa_rvm_h_ortho = nothing
        end

        phi0_l = isnothing(rvm_params_l) ? nothing : rvm_params_l.phi0
        phi0_h = isnothing(rvm_params_h) ? nothing : rvm_params_h.phi0

        Plot.position_angle(lon_l, pa_l, pa_err_l, I_l, Lin_l, V_l,
                            lon_h, pa_h, pa_err_h, I_h, Lin_h, V_h,
                            indir; show_=true,
                            lon_rvm_l=lon_rvm_l, pa_rvm_l=pa_rvm_l,
                            pa_rvm_l_ortho=pa_rvm_l_ortho, phi0_l=phi0_l,
                            lon_rvm_h=lon_rvm_h, pa_rvm_h=pa_rvm_h,
                            pa_rvm_h_ortho=pa_rvm_h_ortho, phi0_h=phi0_h)
    end

    function geometry_analysis(indir; chi2_red_max=10.0)

        # parameters file
        p = Tools.read_params(joinpath(indir, "params.json"))

        lt = joinpath(indir, "pulsar_low.txt")
        ht = joinpath(indir, "pulsar_high.txt")

        l = load_ascii_all(lt)
        h = load_ascii_all(ht)

        bin_st  = p["bin_st"]
        bin_end = p["bin_end"]

        pulses_l, bins_l, _ = size(l)
        pulses_h, bins_h, _ = size(h)

        db_l = (bin_end + 1) - bin_st
        lon_l = collect(range(-360.0 * db_l / bins_l / 2, 360.0 * db_l / bins_l / 2, length=db_l))
        db_h = (bin_end + 1) - bin_st
        lon_h = collect(range(-360.0 * db_h / bins_h / 2, 360.0 * db_h / bins_h / 2, length=db_h))

        Q_l   = vec(mean(l[:, bin_st:bin_end, 2], dims=1))
        U_l   = vec(mean(l[:, bin_st:bin_end, 3], dims=1))
        Lin_l = sqrt.(Q_l.^2 .+ U_l.^2)
        Q_h   = vec(mean(h[:, bin_st:bin_end, 2], dims=1))
        U_h   = vec(mean(h[:, bin_st:bin_end, 3], dims=1))
        Lin_h = sqrt.(Q_h.^2 .+ U_h.^2)

        # Mean Stokes I for W10 estimation
        I_l_full = vec(mean(l[:, :, 1], dims=1))
        I_h_full = vec(mean(h[:, :, 1], dims=1))

        # sigma_avg: statistically correct noise on the mean profile, used both
        # for detection threshold (5σ on L) and for PA error (0.5·σ/L). Combined
        # with the χ²_red-display mode in Plot.geometry (cutoff 10), this gives
        # the Johnston+ 2023 Fig. 1 convention: minimum χ²_red ≈ 1 at best fit,
        # the "χ² > 10" region masked as blank.
        sigma_avg_l = std(l[:, 1:bin_st-1, 1]) / sqrt(pulses_l)
        sigma_avg_h = std(h[:, 1:bin_st-1, 1]) / sqrt(pulses_h)
        thresh_l = 5.0 * sigma_avg_l
        thresh_h = 5.0 * sigma_avg_h

        pa_l     = [Lin_l[i] > thresh_l ? 0.5 * atan(U_l[i], Q_l[i]) * (180.0/π) : NaN for i in 1:db_l]
        pa_err_l = [Lin_l[i] > thresh_l ? 0.5 * sigma_avg_l / Lin_l[i] * (180.0/π) : NaN for i in 1:db_l]
        pa_h     = [Lin_h[i] > thresh_h ? 0.5 * atan(U_h[i], Q_h[i]) * (180.0/π) : NaN for i in 1:db_h]
        pa_err_h = [Lin_h[i] > thresh_h ? 0.5 * sigma_avg_h / Lin_h[i] * (180.0/π) : NaN for i in 1:db_h]

        # Match the preprocessing used in position_angle: undo OPM jumps,
        # then select the same OPM branch (global +90° shift).
        pa_l, _ = deopm_pa(pa_l)
        pa_h, _ = deopm_pa(pa_h)
        pa_shift = 90.0
        pa_l = [isnan(x) ? x : mod(x + pa_shift + 90, 180) - 90 for x in pa_l]
        pa_h = [isnan(x) ? x : mod(x + pa_shift + 90, 180) - 90 for x in pa_h]

        # Two-pass with floored rescaling: first fit with σ_avg, then inflate
        # errors by √χ²_red (if > 1) so χ²_red_min = 1 by construction. Floor at
        # 1 prevents deflating σ when the model fits below statistical noise
        # (e.g., sparse high-freq PA points) — there the region stays wide.
        println("Fitting RVM (low frequency, pass 1)...")
        res_l_pre = fit_rvm(lon_l, pa_l, pa_err_l;
            n_alpha=171, n_beta=161, n_phi0=200)
        scale_l = max(1.0, sqrt(res_l_pre.chi2_red))
        println("  χ²_red=$(round(res_l_pre.chi2_red,digits=2))  → σ scale factor = $(round(scale_l,digits=2))")

        println("Fitting RVM chi² map (low frequency, pass 2)...")
        res_l, chi2_l, alphas_deg, betas_deg = fit_rvm(lon_l, pa_l, pa_err_l .* scale_l;
            return_map=true, n_alpha=171, n_beta=161, n_phi0=200)
        println("  best: α=$(round(res_l.alpha,digits=1))° β=$(round(res_l.beta,digits=1))° χ²/ndof=$(round(res_l.chi2_red,digits=4))")

        println("Fitting RVM (high frequency, pass 1)...")
        res_h_pre = fit_rvm(lon_h, pa_h, pa_err_h;
            n_alpha=171, n_beta=161, n_phi0=200)
        scale_h = max(1.0, sqrt(res_h_pre.chi2_red))
        println("  χ²_red=$(round(res_h_pre.chi2_red,digits=2))  → σ scale factor = $(round(scale_h,digits=2))")

        println("Fitting RVM chi² map (high frequency, pass 2)...")
        res_h, chi2_h, _, _ = fit_rvm(lon_h, pa_h, pa_err_h .* scale_h;
            return_map=true, n_alpha=171, n_beta=161, n_phi0=200)
        println("  best: α=$(round(res_h.alpha,digits=1))° β=$(round(res_h.beta,digits=1))° χ²/ndof=$(round(res_h.chi2_red,digits=4))")

        # W10 (pulse width at 10 % of peak flux), cf. Johnston et al. 2023, Eq. (3)
        nbin  = p["nbin"]
        W10_l = pulse_width_fraction(I_l_full, nbin; frac=0.1,
                                     on_st=bin_st, on_end=bin_end)
        W10_h = pulse_width_fraction(I_h_full, nbin; frac=0.1,
                                     on_st=bin_st, on_end=bin_end)
        P_sec = p["period"]
        println("P = $P_sec s   W10: low=$(round(W10_l, digits=2))°  high=$(round(W10_h, digits=2))°")

        # Diagnostic: χ²_red percentiles — helps see if a seemingly "missing"
        # region is just above the cutoff or genuinely a bad fit.
        function _stats(chi2_map, ndof)
            red = chi2_map[isfinite.(chi2_map)] ./ ndof
            sort!(red)
            n = length(red)
            q(p) = red[clamp(round(Int, p*n), 1, n)]
            (min=red[1], p50=q(0.5), p90=q(0.9), p95=q(0.95), p99=q(0.99), max=red[end])
        end
        st_l = _stats(chi2_l, res_l.ndof)
        st_h = _stats(chi2_h, res_h.ndof)
        println("χ²_red low:  min=$(round(st_l.min,digits=2))  p50=$(round(st_l.p50,digits=2))  p90=$(round(st_l.p90,digits=2))  p95=$(round(st_l.p95,digits=2))  p99=$(round(st_l.p99,digits=2))  max=$(round(st_l.max,digits=2))")
        println("χ²_red high: min=$(round(st_h.min,digits=2))  p50=$(round(st_h.p50,digits=2))  p90=$(round(st_h.p90,digits=2))  p95=$(round(st_h.p95,digits=2))  p99=$(round(st_h.p99,digits=2))  max=$(round(st_h.max,digits=2))")

        Plot.geometry(chi2_l, chi2_h, alphas_deg, betas_deg, indir;
                      show_=true, P_sec=P_sec, W_deg_l=W10_l, W_deg_h=W10_h,
                      ndof_l=res_l.ndof, ndof_h=res_h.ndof,
                      chi2_red_max=chi2_red_max)
    end

    """
    W_frac (e.g. W10): distance between the outermost longitude bins above
    `frac * peak`, following Posselt et al. 2021 / Johnston et al. 2023.
    If `on_st`/`on_end` are given, the peak and search are restricted to that
    on-pulse window; this avoids noise spikes extending the width artificially.
    """
    function pulse_width_fraction(profile, nbin; frac=0.1,
                                   on_st=nothing, on_end=nothing)
        st = isnothing(on_st)  ? 1              : on_st
        en = isnothing(on_end) ? length(profile) : on_end
        seg = @view profile[st:en]
        threshold = frac * maximum(seg)
        above = findall(>=(threshold), seg)
        isempty(above) && return NaN
        return (last(above) - first(above) + 1) * 360.0 / nbin
    end

end # module
