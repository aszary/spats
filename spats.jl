module SpaTs
    using ArgParse
    using Glob
    using JSON

    include("modules/data.jl")
    include("modules/plot.jl")
    include("modules/tools.jl")
    include("modules/functions.jl")


    function test(outdir)
        d = Data.load_ascii("input/1.txt")
        Plot.single(d, outdir; darkness=0.3, number=256, bin_st=400, bin_end=600, start=1, name_mod="1", show_=true)
        Plot.average(d, outdir; number=256, bin_st=400, bin_end=600, start=1, name_mod="1", show_=true)
        Plot.lrfs(d, outdir; darkness=0.1, start=1, name_mod="1", bin_st=500, bin_end=530, show_=true)
    end

    
    function test2(outdir)    
        Data.convert_psrfit_ascii("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00000-00255.spCF", outdir*"1.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00256-00511.spCF", outdir*"2.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00512-00767.spCF", outdir*"3.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00768-01029.spCF", outdir*"4.txt")  
        data1 = Data.load_ascii(outdir*"1.txt")
        data2 = Data.load_ascii(outdir*"2.txt")
        data3 = Data.load_ascii(outdir*"3.txt")
        data4 = Data.load_ascii(outdir*"4.txt")
        data = vcat(data1, data2, data3, data4)
        Plot.single(data, outdir; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319", show_=true)
        Plot.average(data, outdir; number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319", show_=true)
        Plot.lrfs(data, outdir; darkness=0.1, start=1, name_mod="J1319", bin_st=400, bin_end=600, show_=true)
        folded = Tools.p3fold(data, 20, 40)
        Plot.single(folded, outdir; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319_p3fold", show_=true)
    end


    function process_psrdata(indir, outdir)
        p = Data.process_psrdata(indir, outdir)
        folded = Data.load_ascii(outdir*"/pulsar.debase.p3fold")
        Plot.p3fold(folded, outdir; start=3, bin_st=p["bin_st"]-20, bin_end=p["bin_end"]+20, name_mod="test", show_=true, repeat_num=4)
    end



    function process_psrdata2(indir, outdir; outfile="pulsar.spCF", files=nothing)

        # === 1. Extract pulsar name from indir path ===
        split_path = splitpath(indir)
        pulsar_name = split_path[end-1]  # assuming indir ends with /<pulsar>/<date>/
        println("Detected pulsar: $pulsar_name")

        # === 2. Create output subdirectory ===
        pulsar_outdir = joinpath(outdir, pulsar_name)
        if !isdir(pulsar_outdir)
            println("Creating output directory: $pulsar_outdir")
            mkpath(pulsar_outdir)
        end

        # === 3. Redirect all output to the new subdirectory ===
        outdir = pulsar_outdir

        # === 4. Prepare output file path ===
        outfile = joinpath(outdir, outfile)

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

        # Combine all files
        run(pipeline(`psradd $file_names -o $outfile`, stderr="errs.txt"))  # PSRCHIVE

        # Debase the data
        run(pipeline(`pmod -device "/xw" -debase $outfile`, `tee pmod_output.txt`))
        output = read("pmod_output.txt", String)
        rm("pmod_output.txt")  # cleanup

        # Extract onpulse values
        m = match(r"-onpulse '(\d+) (\d+)'", output)
        if !isnothing(m)
            bin_st, bin_end = parse.(Int, m.captures)
            region_length = bin_end - bin_st + 1
            if region_length % 2 != 0
                println("Warning: Onpulse region length ($region_length) is not even. Adjusting bin_end to make it even.")
                bin_end -= 1
                println("Adjusted onpulse range: $bin_st to $bin_end")
            end
            println("Found onpulse range: $bin_st to $bin_end")
        end

        debased_file = replace(outfile, ".spCF" => ".debase.gg")

        # Compute 2DFS and LRFS
        run(pipeline(`pspec -w -2dfs -lrfs -profd "/NULL" -onpulsed "/NULL" -2dfsd "/NULL" -lrfsd "/NULL" -nfft 256 -onpulse "$(bin_st) $(bin_end)" $debased_file`, stderr="errs.txt"))

        # Find P3 value
        run(pipeline(`pspecDetect -v -device "/xw" $debased_file`, `tee pspecDetect_output.txt`))
        output = read("pspecDetect_output.txt", String)
        rm("pspecDetect_output.txt")  # cleanup

        # Extract the last reported P3 value
        p3_matches = collect(eachmatch(r"P3\[P0\]\s*=\s*(\d+\.\d+)\s*\+-\s*(\d+\.\d+)", output))
        if !isempty(p3_matches)
            last_match = p3_matches[end]
            p3_value = parse(Float64, last_match.captures[1])
            p3_error = parse(Float64, last_match.captures[2])
            println("Found P3 = $p3_value ± $p3_error P0")
        end

        ybins = Functions.find_ybins(p3_value)
        println("Number of ybins: $ybins")

        
        run(pipeline(`pfold -p3fold "$p3_value $ybins" -onpulse "$bin_st $bin_end" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $debased_file`, stderr="errs.txt"))

        # --- Dodajemy teraz ładowanie i rysowanie wykresu p3fold ---
        p3fold_file = joinpath(outdir, "pulsar.debase.p3fold")


        # Załaduj dane
        folded = Data.load_ascii(p3fold_file)

        # Wywołaj funkcję rysującą
        Plot.p3fold(folded, joinpath(outdir, pulsar_name);
                    start=3,
                    bin_st=bin_st - 20,
                    bin_end=bin_end + 20,
                    name_mod="auto_p3fold",
                    show_=true,
                    repeat_num=4)

        # --- Zwróć wartości potrzebne na zewnątrz ---
        return bin_st - 20, bin_end + 20
    end



    function main()
        # output directory for VPM
        vpmout = "/home/psr/output/"

        #d = Data.load_ascii("/home/psr/output/pulsar.debase.1.2dfs")
        #println(size(d))
        #Plot.twodfs_plot(d, vpmout; show_=true)

        #d = Data.load_ascii_all("/home/psr/output/pulsar.debase.lrfs")
        #println(size(d))
        #Plot.lrfs_plot(d[:,:,3], vpmout; show_=true)

        #test(vpmout)
        #test2(vpmout)
        process_psrdata2("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)

    end

end # module

SpaTs.main()

println("Bye")