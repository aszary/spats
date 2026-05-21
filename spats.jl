module SpaTs
    using ArgParse
    using Glob
    using JSON
    using Printf
    using Statistics

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
                             snr_threshold=5.0,
                             period=period, show_=false)
            catch e
                @warn "[$name] failed: $e"
            end
        end
    end

    # -------------------------------------------------------------------------
    # Johnston+2023 batch analysis
    # Automatically finds .ar files in data_dir, converts, auto-detects
    # on-pulse region, gets period via vap, and runs rvm_analysis for all
    # pulsars that appear in figures of Johnston et al. 2023 (MNRAS 520, 4801).
    # -------------------------------------------------------------------------
    function run_johnston2023(data_dir, outdir_base; snr_threshold=3.0, show_=false)
        mkpath(outdir_base)

        # (name, figure_label) — every pulsar shown in a Johnston+2023 figure
        johnston_targets = [
            ("J1048-5832", "Fig.1 RVM"),
            ("J1722-4400", "Fig.1 RVM"),
            ("J1110-5637", "Fig.1 RVM"),
            ("J1842-0359", "Fig.1 RVM"),
            ("J0108-1431", "Fig.2 flat"),
            ("J1232-4742", "Fig.2 flat"),
            ("J1511-5414", "Fig.2 flat"),
            ("J1524-5625", "Fig.2 flat"),
            ("J1625-4048", "Fig.3 non-RVM"),
            ("J2046-0421", "Fig.3 non-RVM"),
            ("J1607-0032", "Fig.3 non-RVM"),
            ("J0846-3533", "Fig.3 non-RVM"),
            ("J0134-2937", "Fig.8 partial"),
            ("J1648-6044", "Fig.8 partial"),
            ("J1615-5444", "Fig.8 partial"),
            ("J0631+1036", "Fig.8 partial"),
            ("J0738-4042", "Fig.9 OPM"),
        ]

        results  = Dict{String,Any}()
        skipped  = String[]

        for (name, fig_label) in johnston_targets
            println("\n=== $name  [$fig_label] ===")

            # --- find .ar file (glob: any date suffix) ---
            candidates = Glob.glob(name * "*.ar", data_dir)
            if isempty(candidates)
                @warn "[$name] no .ar file found in $data_dir — skipping"
                push!(skipped, name)
                continue
            end
            infile = first(candidates)
            println("  file: $(basename(infile))")

            outdir = joinpath(outdir_base, name) * "/"
            mkpath(outdir)

            # --- period via vap (needed for h_Blask) ---
            period = nothing
            try
                raw = readchomp(`vap -c period $infile`)
                period = parse(Float64, split(strip(raw))[end])
                println("  period: $(round(period, digits=6)) s")
            catch
                @warn "[$name] vap failed — h_Blask will not be computed"
            end

            # --- convert .ar → ASCII (cached: skip if already done) ---
            txt_file = joinpath(outdir, name * "_pol.txt")
            if !isfile(txt_file)
                println("  converting .ar → $txt_file ...")
                Data.convert_psrfit_ascii_pol(infile, txt_file)
            else
                println("  using cached $txt_file")
            end

            # --- load and auto-detect on-pulse region ---
            data4 = Data.load_ascii_all(txt_file)
            nbin  = size(data4, 2)
            I_avg = vec(mean(data4[:, :, 1], dims=1))
            thr   = 0.1 * maximum(I_avg)
            on    = findall(I_avg .> thr)
            if isempty(on)
                @warn "[$name] on-pulse detection failed — skipping"
                push!(skipped, name)
                continue
            end
            margin  = max(20, round(Int, 0.05 * nbin))
            bin_st  = max(1,    first(on) - margin)
            bin_end = min(nbin, last(on)  + margin)
            println("  on-pulse: bins $bin_st – $bin_end  ($(bin_end - bin_st + 1) / $nbin bins)")

            # --- run RVM fit and plot ---
            try
                result = Plot.rvm(data4, outdir, name;
                                  bin_st=bin_st, bin_end=bin_end,
                                  snr_threshold=snr_threshold,
                                  period=period, show_=show_, n_freq=1)
                beta = round(result.zeta - result.alpha, digits=1)
                println("  α=$(round(result.alpha,digits=1))°  β=$(beta)°  " *
                        "φ₀=$(round(result.phi0,digits=1))°  " *
                        "χ²ᵣ=$(round(result.chi2_red,digits=2))")
                if !isnan(result.h_blask)
                    println("  h_Blask = $(round(result.h_blask, digits=0)) km")
                end
                results[name] = (result=result, fig=fig_label)
            catch e
                @warn "[$name] RVM failed: $e"
                push!(skipped, name)
            end
        end

        # --- summary table ---
        println("\n" * "="^72)
        println("Johnston+2023 RVM summary")
        println("="^72)
        @printf("%-14s  %-18s  %6s  %6s  %6s  %6s  %8s\n",
                "Pulsar", "Figure", "α[°]", "β[°]", "φ₀[°]", "χ²ᵣ", "h_Blask")
        println("-"^72)
        for (name, fig_label) in johnston_targets
            if haskey(results, name)
                r = results[name].result
                hb = isnan(r.h_blask) ? "    —  " : @sprintf("%6.0f km", r.h_blask)
                @printf("%-14s  %-18s  %6.1f  %6.1f  %6.1f  %6.2f  %s\n",
                        name, fig_label,
                        r.alpha, r.zeta - r.alpha, r.phi0, r.chi2_red, hb)
            else
                @printf("%-14s  %-18s  %s\n", name, fig_label, "SKIPPED (no file or fit failed)")
            end
        end
        println("="^72)
        if !isempty(skipped)
            println("Skipped: ", join(skipped, ", "))
        end

        return results
    end

    function main()
        data_dir = "/home/psr/data/posselt/ar_files/"
        rvm_out  = "/home/psr/output/posselt_rvm/"

        run_johnston2023(data_dir, rvm_out)








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
                          snr_threshold=5.0,
                          period=nothing,
                          n_freq=1,
                          show_=false)

        mkpath(outdir)
        txt_file = joinpath(outdir, name_mod * "_pol.txt")
        isfile(txt_file) || Data.convert_psrfit_ascii_pol(infile, txt_file)
        data4  = Data.load_ascii_all(txt_file)
        result = Plot.rvm(data4, outdir, name_mod;
                          bin_st=bin_st, bin_end=bin_end,
                          snr_threshold=snr_threshold,
                          period=period, show_=show_, n_freq=n_freq)

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