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


    function main()
        # output directory for VPM
        vpmout = "/home/psr/output/"

        Data.analyse_offsets(vpmout*"J1842-0359_16", "norefine")
        Data.analyse_offsets(vpmout*"J1539-6322_16", "norefine")
        Data.analyse_offsets(vpmout*"J2139+2242_16", "norefine")
        Data.analyse_offsets(vpmout*"J1232-4742_16", "norefine")
        #Data.analyse_offsets(vpmout*"J1741-0840_16", "norefine", 0.1549; n_comp=2)
        Data.analyse_offsets(vpmout*"J1512-5431_16", "norefine")

        RVM polarimetry — J1842-0359 (odkomentuj aby uruchomić)
        test_rvm("/home/psr/data/new/J1842-0359/2019-11-05-18:03:43/pulsar.debase.gg",
                 vpmout*"J1842-0359_rvm/", "J1842-0359";
                 bin_st=450, bin_end=580, period=0.536)








        #Data.process_all_data(vpmout)
        #Data.combine_pngs_to_pdf(vpmout)
        #Data.combine_pngs(vpmout)
        
        #Data.remove_folders(vpmout)
        #Data.remove_notinteresting("input/pulsars_interesting.txt", vpmout)
    end

    function test_rvm(infile, outdir, name_mod; bin_st, bin_end,
                      period=nothing, alpha=nothing, beta=nothing)
        rvm_analysis(infile, outdir, name_mod;
                     bin_st=bin_st, bin_end=bin_end,
                     period=period, alpha=alpha, beta=beta,
                     show_=true)
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
                          period=nothing, alpha=nothing, beta=nothing,
                          show_=false)

        txt_file = joinpath(outdir, name_mod * "_pol.txt")

        # export 4-pol ASCII from PSRFIT file (requires PSRCHIVE on PATH)
        Data.convert_psrfit_ascii_pol(infile, txt_file)

        # load full Stokes array  (pulses × bins × 4)
        data4 = Data.load_ascii_all(txt_file)

        # plot + fit
        result = Plot.rvm(data4, outdir, name_mod;
                          bin_st=bin_st, bin_end=bin_end,
                          snr_threshold=snr_threshold,
                          linpol_threshold=linpol_threshold,
                          period=period, alpha=alpha, beta=beta,
                          show_=show_)

        println("RVM fit done: alpha=$(round(result.alpha, digits=2)) deg  " *
                "zeta=$(round(result.zeta, digits=2)) deg  " *
                "phi0=$(round(result.phi0, digits=2)) deg  " *
                "chi2r=$(round(result.chi2_red, digits=3))")
        if period !== nothing && !isempty(result.h_blask)
            h_med = median(result.h_blask)
            println("  Blaskiewicz height (median): $(round(h_med, digits=0)) km")
        end
        return result
    end

end # module

SpaTs.main()

println("Bye")