module Data
    using Glob
    using FITSIO
    using ProgressMeter
    using CairoMakie, FileIO
    using PDFIO
    using Statistics
    using FFTW
    
    # --- Implementacja funkcji dopasowania Gaussowskiego (inline) ---
    struct ComponentOffset
        offset_bins::Float64
    end

    function gaussian_model(x, p)
        y = zeros(eltype(x), length(x))
        n = length(p) ÷ 3
        for i in 1:n
            amp = p[3*(i-1)+1]
            mu  = p[3*(i-1)+2]
            sig = p[3*(i-1)+3]
            @. y += amp * exp(-0.5 * ((x - mu) / sig)^2)
        end
        return y
    end

    function fit_gaussians(x, y, n::Int; p0, lower=Float64[], upper=Float64[])
        if isempty(lower) || isempty(upper)
            return curve_fit((x, p) -> gaussian_model(x, p), x, y, p0)
        else
            return curve_fit((x, p) -> gaussian_model(x, p), x, y, p0; lower=lower, upper=upper)
        end
    end

    function component_offsets(fit_h, fit_l; nbin=1024, period=1.0)
        offsets = ComponentOffset[]
        n = length(fit_h.param) ÷ 3
        for i in 1:n
            mu_h = fit_h.param[3*(i-1)+2]
            mu_l = fit_l.param[3*(i-1)+2]
            push!(offsets, ComponentOffset(mu_h - mu_l))
        end
        return offsets
    end
    
    # Funkcja pomocnicza do sortowania parametrów Gaussów (wymusza kolejność L->R)
    function sort_gaussians_params(p)
        n = length(p) ÷ 3
        comps = []
        for i in 1:n
            push!(comps, p[3*(i-1)+1:3*i])
        end
        sort!(comps, by = x -> x[2]) # sortuj po pozycji mu (center)
        return vcat(comps...)
    end
    
    # Prosta funkcja do znajdowania szczytów (zastępuje findpeaks z DSP)
    function find_peaks_simple(data; threshold=0.1)
        peaks = Tuple{Int, Float64}[]
        len = length(data)
        # Wygładzanie prostą średnią ruchomą (okno 3) dla redukcji szumu przed szukaniem
        smoothed = [mean(data[max(1, i-1):min(len, i+1)]) for i in 1:len]
        
        for i in 2:len-1
            if smoothed[i] > smoothed[i-1] && smoothed[i] > smoothed[i+1] && smoothed[i] > threshold
                push!(peaks, (i, smoothed[i]))
            end
        end
        # Sortuj malejąco po wysokości
        sort!(peaks, by=x->x[2], rev=true)
        return peaks
    end
    
    # Funkcja normalizująca wektor do zakresu 0-1
    function normalize_vec(v)
        return (v .- minimum(v)) ./ (maximum(v) - minimum(v))
    end

    using LsqFit, DSP 
    using DelimitedFiles, Printf
    
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
    Główna funkcja wywoływana ze SpaTs.main()
    """
    function analyse_offsets_gauss(indir, type, period)
        # --- 1. Konfiguracja i Wczytanie danych ---
        p = Tools.read_params(joinpath(indir, "params.json"))
        low_path = joinpath(indir, "pulsar_low.debase.p3fold_" * type)
        high_path = joinpath(indir, "pulsar_high.debase.p3fold_" * type)
        
        l_matrix = load_ascii(low_path)
        h_matrix = load_ascii(high_path)

        # Obliczamy profil średni (zintegrowany) z P3 foldów - to drastycznie redukuje szum
        # Używamy dims=1 aby uśrednić po fazach P3
        mean_l = vec(mean(l_matrix, dims=1))
        mean_h = vec(mean(h_matrix, dims=1))
        
        n_bins = length(mean_l)
        bins = collect(1.0:n_bins)

        println(">>> Starting Integrated Gaussian Analysis (High SNR): ", indir)

        # Normalizacja
        norm_l = normalize_vec(mean_l)
        norm_h = normalize_vec(mean_h)

        # --- 2. Automatyczna detekcja liczby komponentów ---
        # Analizujemy profil Low Freq (zazwyczaj szerszy/wyraźniejszy)
        detected_peaks = find_peaks_simple(norm_l, threshold=0.15)
        
        # Jeśli nie wykryto, zakładamy 1 główny pik
        if isempty(detected_peaks)
            max_val, max_idx = findmax(norm_l)
            push!(detected_peaks, (max_idx, max_val))
        end
        
        # Ograniczamy do max 5 komponentów, sortujemy po pozycji (od lewej do prawej)
        n_comps = min(length(detected_peaks), 5)
        sort!(detected_peaks, by=x->x[1]) # sortuj po pozycji (bin)
        
        println(">>> Automatically detected $n_comps Gaussian component(s) at bins: ", [p[1] for p in detected_peaks])

        # Generowanie parametrów startowych (p0) i granic (bounds)
        p0_guess = Float64[]
        lower_bounds = Float64[]
        upper_bounds = Float64[]
        
        est_width = n_bins / 20.0 # Startowa szerokość ~5% profilu
        
        for (pos, amp) in detected_peaks[1:n_comps]
            # p0: [Amp, Mu, Sigma]
            push!(p0_guess, amp, Float64(pos), est_width)
            
            # Bounds
            append!(lower_bounds, [0.0, 1.0, 0.5])                         # Min: Amp>0, Pos>1, Sig>0.5
            append!(upper_bounds, [2.0, Float64(n_bins), Float64(n_bins)/4.0]) # Max: Amp<2, Pos<N, Sig<1/4 profilu
        end

        # --- 3. Dopasowanie Profilu Średniego (Integrated) ---
        println("--- Fitting Integrated Low Frequency Profile ---")
        fit_l = fit_gaussians(bins, norm_l, n_comps; p0=p0_guess, lower=lower_bounds, upper=upper_bounds)
        params_l = sort_gaussians_params(fit_l.param)

        # Wyrównanie zgrubne (CCF) dla High Freq
        println("--- Aligning High Frequency Profile (CCF) ---")
        ccf_val = real.(ifft(fft(norm_h) .* conj.(fft(norm_l))))
        global_shift = argmax(ccf_val) - 1
        if global_shift > n_bins ÷ 2; global_shift -= n_bins; end
        println(">>> Detected global profile shift (High vs Low): ", global_shift, " bins")

        # Startowe dla High = parametry z Low + globalne przesunięcie
        p0_h = copy(params_l)
        for k in 1:n_comps
            idx = 3*(k-1) + 2
            new_pos = p0_h[idx] + global_shift
            if new_pos < 1.0; new_pos += n_bins; end
            if new_pos > n_bins; new_pos -= n_bins; end
            p0_h[idx] = clamp(new_pos, 2.0, Float64(n_bins)-2.0)
        end

        println("--- Fitting Integrated High Frequency Profile ---")
        fit_h = fit_gaussians(bins, norm_h, n_comps; p0=p0_h, lower=lower_bounds, upper=upper_bounds)
        params_h = sort_gaussians_params(fit_h.param)

        # Obliczanie offsetów dla średniego profilu
        offsets = ComponentOffset[]
        for k in 1:n_comps
            mu_l = params_l[3*(k-1)+2]
            mu_h = params_h[3*(k-1)+2]
            diff = mu_h - mu_l
            diff = mod(diff + n_bins/2, n_bins) - n_bins/2
            push!(offsets, ComponentOffset(diff))
        end

        # --- Rysowanie: Profil Średni ---
        fig = Figure(size = (1000, 900))
        colors = [:red, :green, :blue, :orange, :purple] # Zapas kolorów
        safe_colors = [colors[mod(i-1, length(colors))+1] for i in 1:n_comps]

        # Panel 1: Low Freq
        ax1 = Axis(fig[1, 1], title="Low Frequency (Integrated)", xlabel="Biny", ylabel="Intensywność")
        scatter!(ax1, bins, norm_l, color=:gray, markersize=3, label="Dane")
        lines!(ax1, bins, gaussian_model(bins, params_l), color=:black, linewidth=2, label="Model")
        for k in 1:n_comps
            lines!(ax1, bins, gaussian_model(bins, params_l[3*(k-1)+1:3*k]), color=safe_colors[k], linestyle=:dash, linewidth=2, label="G$k")
        end
        axislegend(ax1)

        # Panel 2: High Freq
        ax2 = Axis(fig[2, 1], title="High Frequency (Integrated)", xlabel="Biny", ylabel="Intensywność")
        scatter!(ax2, bins, norm_h, color=:gray, markersize=3, label="Dane")
        lines!(ax2, bins, gaussian_model(bins, params_h), color=:black, linewidth=2, label="Model")
        for k in 1:n_comps
            lines!(ax2, bins, gaussian_model(bins, params_h[3*(k-1)+1:3*k]), color=safe_colors[k], linestyle=:dash, linewidth=2, label="G$k")
        end
        axislegend(ax2)

        # Panel 3: Porównanie komponentów
        ax3 = Axis(fig[3, 1], title="Przesunięcie średnie (Ciągła=Low, Przerywana=High)", xlabel="Biny")
        for k in 1:n_comps
            p_l_k = params_l[3*(k-1)+1:3*k]
            p_h_k = params_h[3*(k-1)+1:3*k]
            
            # Rysujemy znormalizowane gauss'y
            y_l = gaussian_model(bins, p_l_k) ./ maximum(gaussian_model(bins, p_l_k))
            y_h = gaussian_model(bins, p_h_k) ./ maximum(gaussian_model(bins, p_h_k))
            
            lines!(ax3, bins, y_l, color=safe_colors[k], linewidth=2, label="L G$k")
            lines!(ax3, bins, y_h, color=safe_colors[k], linestyle=:dash, linewidth=2, label="H G$k")
            off_val = offsets[k].offset_bins
            text!(ax3, params_l[3*(k-1)+2], 0.5, text=@sprintf("Δ=%.2f", off_val), color=safe_colors[k], fontsize=14, align=(:center, :bottom))
        end
        
        out_name_base = "gaussian_fit_integrated_$(type)"
        save(joinpath(indir, out_name_base * ".pdf"), fig)
        
        # Zapis statystyk średnich
        out_path = joinpath(indir, "gaussian_offsets_integrated_$(type).csv")
        csv_data = Matrix{Any}(undef, n_comps, 2)
        println("\n--- Integrated Fit Results ---")
        for k in 1:n_comps
             @printf("Component G%d Offset: %.4f bins\n", k, offsets[k].offset_bins)
            csv_data[k, 1] = k
            csv_data[k, 2] = offsets[k].offset_bins
        end
        writedlm(out_path, csv_data, ',')

        # --- 4. Analiza fazowa P3 (Phase-resolved) ---
        println("\n--- Starting Phase-Resolved Analysis (per P3 bin) ---")
        
        n_phases = size(l_matrix, 1)
        phase_offsets = fill(NaN, n_phases, n_comps)
        
        # Pętla po wierszach (fazach P3)
        for i in 1:n_phases
            row_l = vec(l_matrix[i, :])
            row_h = vec(h_matrix[i, :])
            
            # Warunek SNR: jeśli max wiersza jest mały w porównaniu do szumu średniego
            if maximum(row_l) < 0.2 * maximum(mean_l) 
                continue 
            end

            try
                n_l_row = normalize_vec(row_l)
                n_h_row = normalize_vec(row_h)

                # Fit Low: Startujemy z parametrów średnich (params_l) aby zachować tożsamość komponentów
                # Ale pozwalamy na lekkie zmiany
                f_l = fit_gaussians(bins, n_l_row, n_comps; p0=params_l, lower=lower_bounds, upper=upper_bounds)
                
                # Fit High: Startujemy z parametrów średnich (params_h)
                f_h = fit_gaussians(bins, n_h_row, n_comps; p0=params_h, lower=lower_bounds, upper=upper_bounds)
                
                # Wyciągnij parametry (zakładamy, że kolejność się zachowała dzięki dobremu p0)
                pl = f_l.param
                ph = f_h.param
                
                for k in 1:n_comps
                    mu_l_ph = pl[3*(k-1)+2]
                    mu_h_ph = ph[3*(k-1)+2]
                    
                    diff = mu_h_ph - mu_l_ph
                    diff = mod(diff + n_bins/2, n_bins) - n_bins/2
                    
                    # Odrzucamy absurdalne skoki (np. > 20% szerokości pulsu od średniej)
                    if abs(diff) < n_bins * 0.2
                        phase_offsets[i, k] = diff
                    end
                end
            catch
                # Błąd dopasowania w danej fazie - ignorujemy
            end
        end

        # --- Rysowanie: Offset vs Faza P3 ---
        fig_phase = Figure(size = (900, 600))
        ax_ph = Axis(fig_phase[1, 1], 
            title="Offset między częstotliwościami w funkcji fazy P3", 
            xlabel="Faza P3 (nr wiersza)", ylabel="Offset (biny) [High - Low]",
            xgridvisible=true, ygridvisible=true)
        
        phases = 1:n_phases
        for k in 1:n_comps
            vals = phase_offsets[:, k]
            # Filtrujemy NaN do rysowania
            mask = .!isnan.(vals)
            if any(mask)
                lines!(ax_ph, phases[mask], vals[mask], label="G$k", color=safe_colors[k], linewidth=2)
                scatter!(ax_ph, phases[mask], vals[mask], color=safe_colors[k], markersize=8)
            end
        end
        axislegend(ax_ph)
        
        save(joinpath(indir, "gaussian_offsets_phase_resolved_$(type).pdf"), fig_phase)
        println(">>> Phase-resolved plot saved: gaussian_offsets_phase_resolved_$(type).pdf")
        
        # Zapis danych fazowych do CSV
        out_ph_csv = joinpath(indir, "gaussian_offsets_phase_resolved_$(type).csv")
        # Nagłówek i dane
        open(out_ph_csv, "w") do io
            println(io, "Phase_Bin," * join(["G$(k)_Offset" for k in 1:n_comps], ","))
            for i in 1:n_phases
                line = ["$i"]
                for k in 1:n_comps
                    push!(line, @sprintf("%.4f", phase_offsets[i, k]))
                end
                println(io, join(line, ","))
            end
        end

        return offsets
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
