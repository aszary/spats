module SpaTs
    using ArgParse
    using Glob
    using JSON

    include("modules/data.jl")
    include("modules/plot.jl")
    include("modules/tools.jl")


    function test(outdir)
        d = Data.load_ascii("input/1.txt")
        Plot.single(d, outdir; darkness=0.3, number=256, bin_st=400, bin_end=600, start=1, name_mod="1", show_=true)
        Plot.average(d, outdir; number=256, bin_st=400, bin_end=600, start=1, name_mod="1", show_=true)
        Plot.lrfs_obsolete(d, outdir; darkness=0.1, start=1, name_mod="1", bin_st=500, bin_end=530, show_=true)
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
        Plot.lrfs_obsolete(data, outdir; darkness=0.1, start=1, name_mod="J1319", bin_st=400, bin_end=600, show_=true)
        folded = Tools.p3fold(data, 20, 40)
        Plot.single(folded, outdir; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319_p3fold", show_=true)
    end


    function test3(indir, outdir)
        p = Data.process_psrdata(indir, outdir)
        Data.convert_psrfit_ascii(joinpath(outdir, "pulsar.debase.gg"), joinpath(outdir, "pulsar.debase.txt"))

        outfile="pulsar.spCF"
        outfile = joinpath(outdir, outfile)
        debased_file = replace(outfile, ".spCF" => ".debase.gg")

        p = Tools.read_params(joinpath(outdir, "params.json"))
        d4 = Data.load_ascii_all(joinpath(outdir, "pulsar.debase.txt"))
        d1 = Data.clean(d4; threshold=0.0031)
        Plot.single(d1, outdir; darkness=0.7, number=100, bin_st=p["bin_st"], bin_end=p["bin_end"], start=210, name_mod="pulsar", show_=true)
        Plot.lrfs_obsolete(d1, outdir; darkness=0.3, start=210, name_mod="pulsar", bin_st=p["bin_st"], bin_end=p["bin_end"], show_=true)
        Data.twodfs_lrfs(debased_file, outdir, p)
        lrfs_file = replace(debased_file, "gg"=>"lrfs")
        data = Data.load_ascii_all(lrfs_file)
        Plot.lrfs(data, outdir, p; darkness=0.3, name_mod="pulsar", show_=true)
        twodfs_file = replace(debased_file, "gg"=>"1.2dfs")
        data_2dfs = Data.load_ascii_all(twodfs_file)
        Plot.twodfs(data_2dfs, outdir, p; darkness=0.3, name_mod="pulsar", show_=true)
        folded = Data.load_ascii(joinpath(outdir, "pulsar.debase.p3fold"))
        Plot.p3fold(folded, outdir; start=3, bin_st=p["bin_st"]-20, bin_end=p["bin_end"]+20, name_mod="pulsar", show_=true, repeat_num=4)        

    end

    """
    For Sardinia poster
    """
    function J1539_6322_Sard(indir, outdir)
        p = Data.process_psrdata(indir, outdir)
        Data.convert_psrfit_ascii(joinpath(outdir, "pulsar.debase.gg"), joinpath(outdir, "pulsar.debase.txt"))
        outfile="pulsar.spCF"
        outfile = joinpath(outdir, outfile)
        debased_file = replace(outfile, ".spCF" => ".debase.gg")

        p = Tools.read_params(joinpath(outdir, "params.json"))
        d4 = Data.load_ascii_all(joinpath(outdir, "pulsar.debase.txt"))
        d1 = Data.clean(d4; threshold=0.001)
        Plot.single(d1, outdir; darkness=0.9, number=150, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="pulsar", show_=true)
    end

    function process_psrdata(indir, outdir)
        p, debased_file, outdir = Data.process_psrdata(indir, outdir)
        Data.plot_psrdata(outdir, p)
        Data.process_data_andrzej(debased_file, outdir, p)
    end

    function process_psrdata_16(indir, outdir)
        p, debased_file, outdir = Data.process_psrdata_16(indir, outdir)
    end

    function process_psrdata_single(indir, outdir)
        Data.process_psrdata_single(indir, outdir)
    end


    function run_all_rvm(outdirs::Vector{String}, data_dir, rvm_outdir)
        mkpath(rvm_outdir)

        # collect unique pulsar names from all output directories
        seen = Set{String}()
        entries = Tuple{String,String}[]  # (name, params_file)
        for outdir in outdirs
            isdir(outdir) || continue
            for d in readdir(outdir)
                endswith(d, "_16") || continue
                isdir(joinpath(outdir, d)) || continue
                name = d[1:end-3]
                name in seen && continue
                params_file = joinpath(outdir, d, "params.json")
                isfile(params_file) || continue
                push!(seen, name)
                push!(entries, (name, params_file))
            end
        end

        for (name, params_file) in sort(entries, by=x->x[1])
            p = Tools.read_params(params_file)
            bin_st_raw = get(p, "bin_st", nothing)
            bin_end_raw = get(p, "bin_end", nothing)
            if isnothing(bin_st_raw) || isnothing(bin_end_raw)
                @warn "[$name] missing bin_st/bin_end in params.json, skipping"
                continue
            end
            bin_st  = Int(bin_st_raw)
            bin_end = Int(bin_end_raw)

            spcf = Glob.glob("*/*_00000-00255.spCF", joinpath(data_dir, name))
            if isempty(spcf)
                @warn "[$name] no .spCF file found, skipping"
                continue
            end
            infile = first(spcf)

            outdir = joinpath(rvm_outdir, name * "_rvm") * "/"
            try
                mkpath(outdir)
            catch e
                @warn "[$name] cannot create output dir $outdir: $e, skipping"
                continue
            end

            # try to extract pulsar period from params.json or via vap
            period = let raw = get(p, "period", nothing)
                isnothing(raw) || raw == -1.0 ? nothing : Float64(raw)
            end
            if isnothing(period)
                try
                    out = readchomp(`vap -c period $infile`)
                    period = parse(Float64, split(strip(out))[end])
                    println("[$name] period from vap: $(period) s")
                catch
                end
            end

            println("\n=== RVM: $name (bin $(bin_st)-$(bin_end)) ===")
            try
                rvm_analysis(infile, outdir, name; bin_st=bin_st, bin_end=bin_end,
                             snr_threshold=3.0, linpol_threshold=0.3,
                             period=period, show_=false)
            catch e
                @warn "[$name] failed: $e"
            end
        end
    end

    function main()
        outdirs  = ["/home/psr/data/OUTPUT/czarek/", "/home/psr/data/OUTPUT/andrzej/"]
        data_dir = "/home/psr/data/new/"
        rvm_out  = "/home/psr/output/"

        run_all_rvm(outdirs, data_dir, rvm_out)








        #Data.process_all_data(vpmout)
        #Data.combine_pngs_to_pdf(vpmout)
        #Data.combine_pngs(vpmout)
        
        #Data.remove_folders(vpmout)
        #Data.remove_notinteresting("input/pulsars_interesting.txt", vpmout)
    end

    function test_rvm(infile, outdir, name_mod; bin_st, bin_end, period=nothing)
        rvm_analysis(infile, outdir, name_mod;
                     bin_st=bin_st, bin_end=bin_end,
                     period=period, show_=true)
    end

    # -------------------------------------------------------------------------
    # RVM analysis
    # Usage:
    #   SpaTs.rvm_analysis(
    #       "/home/psr/data/.../pulsar.debase.gg",   # 4-pol PSRFIT file
    #       "output/J1842-0359/",                    # output directory
    #       "J1842-0359";                            # name prefix
    #       bin_st=450, bin_end=580,                 # on-pulse region
    #       period=0.536,                            # pulse period [s] — enables height panel
    #       alpha=30.0,                              # magnetic inclination [deg] — enables dipole heights
    #       beta=5.0                                 # impact parameter [deg]
    #   )
    # -------------------------------------------------------------------------
    function rvm_analysis(infile, outdir, name_mod;
                          bin_st, bin_end,
                          snr_threshold=3.5, linpol_threshold=0.8,
                          period=nothing,
                          show_=false)

        mkpath(outdir)
        txt_file = joinpath(outdir, name_mod * "_pol.txt")
        isfile(txt_file) || Data.convert_psrfit_ascii_pol(infile, txt_file)
        data4  = Data.load_ascii_all(txt_file)
        result = Plot.rvm(data4, outdir, name_mod;
                          bin_st=bin_st, bin_end=bin_end,
                          snr_threshold=snr_threshold,
                          linpol_threshold=linpol_threshold,
                          period=period, show_=show_, n_freq=2)

        println("RVM fit:  alpha=$(round(result.alpha, digits=1)) deg  " *
                "zeta=$(round(result.zeta, digits=1)) deg  " *
                "phi0=$(round(result.phi0, digits=1)) deg  " *
                "chi2r=$(round(result.chi2_red, digits=2))")
        println("          W10=$(round(result.W10, digits=1)) deg  " *
                "phi_center=$(round(result.phi_center, digits=1)) deg")
        if !isnan(result.h_blask)
            println("  h_Blask = $(round(result.h_blask, digits=0)) km  " *
                    "(Delphi=$(round(result.phi0 - result.phi_center, digits=1)) deg)")
        end
        return result
    end

end # module

SpaTs.main()

println("Bye")