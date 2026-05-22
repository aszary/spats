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


    """
    Automatically process all pulsars found in `dataroot`.
    For each pulsar directory, picks the lexicographically first observation
    subdirectory and runs process_psrdata_16 + analyse_p3folds_16_new.
    Errors for individual pulsars are caught so the loop continues.

    Arguments:
      dataroot  – parent directory containing JXXXX-XXXX subdirs
      vpmout    – output directory prefix (same convention as main())
      n_comp    – number of components passed to analyse_p3folds_16_new
    """
    function analyse_all(dataroot="/home/psr/data/new/", vpmout="/home/psr/output/"; n_comp=2)
        isdir(dataroot) || error("dataroot not found: $dataroot")
        psr_dirs = sort(filter(d -> isdir(joinpath(dataroot, d)), readdir(dataroot)))
        isempty(psr_dirs) && (@warn "No pulsar directories found in $dataroot"; return)

        for psr in psr_dirs
            psr_path = joinpath(dataroot, psr)
            obs_dirs = sort(filter(d -> isdir(joinpath(psr_path, d)), readdir(psr_path)))
            if isempty(obs_dirs)
                @warn "No observation subdirectory for $psr, skipping"
                continue
            end
            indir  = joinpath(psr_path, obs_dirs[1]) * "/"
            outdir = vpmout * psr * "_16"
            println("=== Processing $psr (obs: $(obs_dirs[1])) ===")
            try
                process_psrdata_16(indir, outdir)
                Data.analyse_p3folds_16_new(outdir, "norefine"; n_comp=n_comp)
            catch e
                @warn "Failed for $psr: $e"
            end
        end
    end

    function main()
        # output directory for VPM
        vpmout = "/home/psr/output/"

        # PSR J0034-0721
        #process_psrdata_16("/home/psr/data/new/J0034-0721/2019-10-18-22:29:51/", vpmout*"J0034-0721_16")
        #Data.analyse_p3folds_16(vpmout*"J0034-0721_16", "norefine")
        #process_psrdata(vpmout*"J0034-0721", vpmout*"J0034-0721") # P. nice
        #Data.analyse_p3folds_16_new(vpmout*"J0034-0721_16", "norefine"; n_comp=2)
        #Data.analyse_p3folds_16_new(vpmout*"J0034-0721_16", "refine"; n_comp=2)
        #Data.position_angle(vpmout*"J0034-0721_16")
        #Data.geometry_analysis(vpmout*"J0034-0721_16")

        # PSR J0134-2937
        #process_psrdata_16("/home/psr/data/new/J0134-2937/2019-10-18-22:46:59/", vpmout*"J0134-2937_16")
        #Data.analyse_p3folds_16_new(vpmout*"J0134-2937_16", "norefine"; n_comp=2)


        # PSR J0151-0635
        #process_psrdata_16("/home/psr/data/new/J0151-0635/2020-04-13-10:00:14/", vpmout*"J0151-0635_16")
        #Data.analyse_p3folds_16(vpmout*"J0151-0635_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J0151-0635_16", "norefine"; n_comp=2)

        # PSR J0255-5304
        #process_psrdata_16("/home/psr/data/new/J0255-5304/2019-10-18-23:29:57/", vpmout*"J0255-5304_16")
        #Data.analyse_p3folds_16(vpmout*"J0255-5304_16", "norefine")

        # PSR J0421-0345
        #process_psrdata_16("/home/psr/data/new/J0421-0345/2019-10-27-21:58:44/", vpmout*"J0421-0345_16")
        #Data.analyse_p3folds_16(vpmout*"J0421-0345_16", "norefine")
        #process_psrdata(vpmout*"J0421-0345", vpmout*"J0421-0345") #

        # PSR J0630-2834
        #process_psrdata_16("/home/psr/data/new/J0630-2834/2019-11-18-21:05:35/", vpmout*"J0630-2834_16")
        #Data.analyse_p3folds_16(vpmout*"J0630-2834_16", "norefine")
        #process_psrdata(vpmout*"J0630-2834", vpmout*"J0630-2834") # not nice

        # PSR J0820-1350
        #process_psrdata("/home/psr/data/uhf/J0820-1350/2020-07-11-10:30:12/816/single/", vpmout*"J0820-1350_uhf")
        #process_psrdata_16("/home/psr/data/new/J0820-1350/2020-01-11-01:05:56/", vpmout*"J0820-1350_16")
        #Data.analyse_p3folds_16(vpmout*"J0820-1350_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J0820-1350_16", "norefine", n_comp=2)
        #process_psrdata_single("/home/psr/data/new/J0820-1350/2020-01-11-01:05:56/", vpmout*"J0820-1350_single") # very bad no signal at least at the beginning
        #process_psrdata(vpmout*"J0820-1350", vpmout*"J0820-1350") # nice

        # PSR J0823+0159
        #process_psrdata("/home/psr/data/uhf/J0823+0159/2023-02-27-20:38:58/816/single", vpmout*"J0823+0159_uhf")
        #process_psrdata_16("/home/psr/data/new/J0823+0159/2019-11-05-23:33:48/", vpmout*"J0823+0159_16")
        #process_psrdata_single("/home/psr/data/new/J0823+0159/2019-11-05-23:33:48/", vpmout*"J0823+0159_single")
        #process_psrdata(vpmout*"J0823+0159", vpmout*"J0823+0159") #  nice single pulses, but p3folds no

        # PSR J0856-6137
        #process_psrdata_16("/home/psr/data/new/J0856-6137/2020-01-04-19:54:45/", vpmout*"J0856-6137_16")
        #Data.analyse_p3folds_16(vpmout*"J0856-6137_16", "norefine")

        # PSR J0904-7459
        #process_psrdata_16("/home/psr/data/new/J0904-7459/2019-10-19-05:41:13", vpmout*"J0904-7459_16")
        #Data.analyse_p3folds_16(vpmout*"J0904-7459_16", "norefine")

        # PSR J0934-5249
        #process_psrdata_16("/home/psr/data/new/J0934-5249/2019-11-05-23:07:57/", vpmout*"J0934-5249_16")
        #Data.analyse_p3folds_16(vpmout*"J0934-5249_16", "norefine")

        # PSR J0959-4809
        #process_psrdata_16("/home/psr/data/new/J0959-4809/2019-10-19-06:23:58/", vpmout*"J0959-4809_16")
        #Data.analyse_p3folds_16(vpmout*"J0959-4809_16", "norefine")

        # PSR J1034-3224
        #process_psrdata_16("/home/psr/data/new/J1034-3224/2019-10-19-07:31:26/", vpmout*"J1034-3224_16")
        #Data.analyse_p3folds_16(vpmout*"J1034-3224_16", "norefine")
        #process_psrdata(vpmout*"J1034-3224", vpmout*"J1034-3224") # single bad, p3folds bad, P3 stable

        # PSR J1048-5832
        #process_psrdata_16("/home/psr/data/new/J1048-5832/2020-08-29-13:00:11/", vpmout*"J1048-5832_16")
        #Data.analyse_p3folds_16(vpmout*"J1048-5832_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1048-5832_16", "norefine"; n_comp=2)
        #Data.position_angle(vpmout*"J1048-5832_16")
        #Data.geometry_analysis(vpmout*"J1048-5832_16")

        # PSR J1110-5637
        #process_psrdata("/home/psr/data/new/J1110-5637/2019-10-19-08:18:28/", vpmout*"J1110-5637")
        #process_psrdata_16("/home/psr/data/new/J1110-5637/2019-10-19-08:18:28/", vpmout*"J1110-5637_16")
        #Data.analyse_p3folds_16(vpmout*"J1110-5637_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1110-5637_16", "norefine"; n_comp=2)
        #Data.position_angle(vpmout*"J1110-5637_16")
        #Data.geometry_analysis(vpmout*"J1110-5637_16")

        # PSR J1114-6100
        #process_psrdata_16("/home/psr/data/new/J1114-6100/2019-10-19-08:30:30/", vpmout*"J1114-6100_16")
        #Data.analyse_p3folds_16(vpmout*"J1114-6100_16c", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1114-6100_16c", "norefine")

        # PSR J1133-6250
        #process_psrdata_16("/home/psr/data/new/J1133-6250/2019-11-06-00:46:43/", vpmout*"J1133-6250_16")
        #process_psrdata_single("/home/psr/data/new/J1133-6250/2019-11-06-00:46:43/", vpmout*"J1133-6250_single") # some offset? probably not
        #process_psrdata(vpmout*"J1133-6250", vpmout*"J1133-6250") # single not stable, P. p3fold (10 ybins)

        # PSR J1137-6700
        #process_psrdata_16("/home/psr/data/new/J1137-6700/2019-11-06-00:36:02/", vpmout*"J1137-6700_16")
        #Data.analyse_p3folds_16_new(vpmout*"J1137-6700_16", "norefine"; n_comp=2)
        #Data.position_angle(vpmout*"J1137-6700_16")

        # PSR J1232-4742
        #process_psrdata_16("/home/psr/data/new/J1232-4742/2020-04-11-22:21:18/", vpmout*"J1232-4742_16")
        #Data.analyse_p3folds_16(vpmout*"J1232-4742_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1232-4742_16", "norefine")

        # PSR J1312-5402
        #process_psrdata_16("/home/psr/data/new/J1312-5402/2020-03-30-19:11:06/", vpmout*"J1312-5402_16")
        #Data.analyse_p3folds_16(vpmout*"J1312-5402_16", "norefine")

        # PSR J1345-6115
        #Data.analyse_p3folds_16_new(vpmout*"J1345-6115_16", "norefine"; n_comp=2)

        # PSR J1414-6802
        #process_psrdata_16("/home/psr/data/new/J1414-6802/2020-07-24-16:17:58/", vpmout*"J1414-6802_16") # start diging there...
        #process_psrdata("/home/psr/data/new/J1414-6802/2020-07-24-16:17:58/", vpmout*"J1414-6802_test") # start diging there...
        #process_psrdata_single("/home/psr/data/new/J1414-6802/2020-07-24-16:17:58/", vpmout*"J1414-6802_single") # weak signal no offset
        #process_psrdata(vpmout*"J1414-6802", vpmout*"J1414-6802") # start diging there...

        # PSR J1512-5431
        #process_psrdata_16("/home/psr/data/new/J1512-5431/2020-04-11-23:16:42/", vpmout*"J1512-5431_16")
        #Data.analyse_p3folds_16(vpmout*"J1512-5431_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1512-5431_16", "norefine", n_comp=2)

        # PSR J1514-4834
        #process_psrdata_16("/home/psr/data/new/J1514-4834/2019-12-16-08:22:16/", vpmout*"J1514-4834_16")
        #process_psrdata(vpmout*"J1514-4834", vpmout*"J1514-4834") #

        # PSR J1524-5706
        #Data.analyse_p3folds_16(vpmout*"JJ1524-5706_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1524-5706_16", "norefine")

        # PSR J1534-4428
        #process_psrdata_16("/home/psr/data/new/J1534-4428/2020-01-03-01:03:02/", vpmout*"J1534-4428_16")
        #Data.analyse_p3folds_16(vpmout*"J1534-4428_16", "norefine")

        # PSR J1539-6322
        #J1539_6322_Sard("/home/psr/data/new/J1539-6322/2020-04-11-23:52:16/", vpmout)
        #process_psrdata_16("/home/psr/data/new/J1539-6322/2020-04-11-23:52:16/", vpmout*"J1539-6322_16")
        #Data.p3fold_psrdata(vpmout*"J1539-6322_16")
        #Data.analyse_p3folds_16(vpmout*"J1539-6322_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1539-6322_16", "norefine", n_comp=2)
        #Data.position_angle(vpmout*"J1539-6322_16")
        #process_psrdata_single("/home/psr/data/new/J1539-6322/2020-04-11-23:52:16/", vpmout*"J1539-6322_single")
        #process_psrdata(vpmout*"J1539-6322", vpmout*"J1539-6322") # single not stable, P. p3fold (10 ybins)

        # PSR J1555-3134
        #process_psrdata_16("/home/psr/data/new/J1555-3134/2019-10-31-06:34:39/", vpmout*"J1555-3134_16")
        #Data.analyse_p3folds_16(vpmout*"J1555-3134_16", "norefine")
        #Data.analyse_p3folds_16(vpmout*"J1555-3134_16", "refine")
        #process_psrdata_single("/home/psr/data/new/J1555-3134/2019-10-31-06:34:39/", vpmout*"J1555-3134_single") # no offset
        #process_psrdata(vpmout*"J1555-3134", vpmout*"J1555-3134") # single not stable, P. p3fold (10 ybins)

        # PSR J1651-5222
        #process_psrdata_16("/home/psr/data/new/J1651-5222/2019-10-31-15:25:01/", vpmout*"J1651-5222_16")
        #Data.analyse_p3folds_16(vpmout*"J1651-5222_16", "norefine")

        # PSR J1651-7642
        #Data.analyse_p3folds_16_new(vpmout*"J1651-7642_16", "norefine"; n_comp=2)
        #Data.position_angle(vpmout*"J1651-7642_16")

        # PSR J1703-1846
        #process_psrdata_16("/home/psr/data/new/J1703-1846/2019-11-05-16:26:36/", vpmout*"J1703-1846_16")
        #Data.analyse_p3folds_16(vpmout*"J1703-1846_16", "norefine")


        # PSR J1741-0840
        #process_psrdata_16("/home/psr/data/new/J1741-0840/2019-11-05-16:44:10/", vpmout*"J1741-0840_16")
        #Data.analyse_p3folds_16(vpmout*"J1741-0840_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1741-0840_16", "norefine")

        # PSR J1744-1610
        #test3("/home/psr/data/new/J1744-1610/2020-06-12-20:00:13/", vpmout)

        # PSR J1807+0756
        #process_psrdata_16("/home/psr/data/new/J1807+0756/2020-06-07-22:03:05/", vpmout*"J1807+0756_16")
        #Data.analyse_p3folds_16(vpmout*"J1807+0756_16", "refine")

        # PSR J1834-0010
        #process_psrdata(vpmout*"J1834-0010", vpmout*"J1834-0010") # P. nice

        # PSR J1842-0359
        #process_psrdata_16("/home/psr/data/new/J1842-0359/2019-11-05-18:03:43/", vpmout*"J1842-0359_16")
        #Data.analyse_p3folds_16(vpmout*"J1842-0359_16", "refine")
        #Data.analyse_p3folds_16(vpmout*"J1842-0359_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1842-0359_16", "norefine", n_comp=4)
        #process_psrdata("/home/psr/data/uhf/J1842-0359/2020-07-25-19:58:19/816/single/", vpmout*"J1842-0359_uhf")
        #process_psrdata_16("/home/psr/data/uhf/J1842-0359/2020-07-25-19:58:19/816/single/", vpmout*"J1842-0359_uhf_16")
        #Data.analyse_p3folds_16(vpmout*"J1842-0359_uhf_16", "norefine")
        #process_psrdata_single("/home/psr/data/new/J1842-0359/2019-11-05-18:03:43/", vpmout*"J1842-0359_single")
        #process_psrdata(vpmout*"J1842-0359", vpmout*"J1842-0359") #  nice single pulses, nice P. p3fold
        #Data.position_angle(vpmout*"J1842-0359_16")
        #Data.geometry_analysis(vpmout*"J1842-0359_16")

        # PSR J1857+0057
        #process_psrdata_16("/home/psr/data/new/J1857+0057/2020-02-03-06:00:48/", vpmout*"J1857+0057_16")
        #Data.analyse_p3folds_16(vpmout*"J1857+0057_16", "norefine")

        # PSR J1900-7951
        #process_psrdata_16("/home/psr/data/new/J1900-7951/2020-05-31-22:21:41/", vpmout*"J1900-7951_16")
        #Data.analyse_p3folds_16(vpmout*"J1900-7951_16", "norefine")

        # PSR J1901-0906
        #process_psrdata_16("/home/psr/data/new/J1901-0906/2019-12-16-13:57:56/", vpmout*"J1901-0906_16")    # P3=3
        #Data.analyse_p3folds_16(vpmout*"J1901-0906_16", "norefine")
        #process_psrdata_16("/home/psr/data/new/J1901-0906/2019-12-16-13:57:56/", vpmout*"J1901-0906_16b")   # P3=5
        #Data.analyse_p3folds_16(vpmout*"J1901-0906_16b", "norefine")

        # PSR J1903+2225
        #process_psrdata_16("/home/psr/data/new/J1903+2225/2020-08-20-20:53:16/", vpmout*"J1903+2225_16")
        #Data.analyse_p3folds_16_new(vpmout*"J1903+2225_16", "norefine"; n_comp=1)

        # PSR J1904-1224
        #process_psrdata_16("/home/psr/data/new/J1904-1224/2020-04-14-05:57:45", vpmout*"J1904-1224_16")
        #Data.analyse_p3folds_16(vpmout*"J1904-1224_16", "norefine")

        # PSR J1907+0731
        #process_psrdata(vpmout*"J1907+0731", vpmout*"J1907+0731") #

        # PSR J1919+0134
        #test3("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)
        #process_psrdata("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)
        #process_psrdata_16("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout*"J1919+0134_16")
        #Data.analyse_p3folds_16(vpmout*"J1919+0134_16", "norefine")

        # PSR J1921+2153
        #process_psrdata_16("/home/psr/data/new/J1921+2153/2020-08-08-19:17:21/", vpmout*"J1921+2153_16")
        #Data.analyse_p3folds_16(vpmout*"J1921+2153_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1921+2153_16", "norefine", n_comp=2)
        #process_psrdata(vpmout*"J1921+2153", vpmout*"J1921+2153") # not nice

        # PSR J2046+1540
        #process_psrdata_16("/home/psr/data/new/J2046+1540/2019-11-26-17:18:47/", vpmout*"J2046+1540_16")
        #Data.analyse_p3folds_16_new(vpmout*"J2046+1540_16", "norefine", n_comp=2)
        #Data.position_angle(vpmout*"J2046+1540_16")


        # PSR J2053-7200
        #process_psrdata_16("/home/psr/data/new/J2053-7200/2019-11-26-18:06:44/", vpmout*"J2053-7200_16")
        #Data.analyse_p3folds_16(vpmout*"J2053-7200_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J2053-7200_16", "norefine", n_comp=2)
        #Data.position_angle(vpmout*"J2053-7200_16")
        #Data.geometry_analysis(vpmout*"J2053-7200_16")

        # PSR J2139+2242
        #process_psrdata_16("/home/psr/data/new/J2139+2242/2020-09-07-20:01:41/", vpmout*"J2139+2242_16")
        #Data.analyse_p3folds_16(vpmout*"J2139+2242_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J2139+2242_16", "norefine", n_comp=2)

        #test(vpmout)
        #test2(vpmout)
        #Data.process_all_data(vpmout)
        #Data.combine_pngs_to_pdf(vpmout)
        #Data.combine_pngs(vpmout)
        #Data.remove_folders(vpmout)
        #Data.remove_notinteresting("input/pulsars_interesting.txt", vpmout)

        Tools.clean_all(vpmout)
        #analyse_all()
    end


    """
    Replicate Johnston+2023 RVM geometry analysis.

    For each pulsar subdirectory under `dataroot` containing a `.ar` file:
      1. Convert .ar → ASCII (cached)
      2. Load data, detect on-pulse, compute off-pulse noise from Q/U
      3. Compute PA, apply deopm_pa and optional OPM shift
      4. Two-pass Data.fit_rvm (scale errors by sqrt(χ²ᵣ) between passes)
      5. Save position_angle.pdf and geometry.pdf to `outdir/<PSR>`

    Parameters:
      dataroot      – parent directory with JXXXX-XXXX subdirs containing .ar files
      outdir        – root output directory (one subdir per pulsar is created)
      snr_threshold – L/σ threshold for PA inclusion (default 5)
      max_pa_err_deg – secondary ceiling on PA error (default 10°)
      pa_shift_deg  – global OPM branch rotation applied after deopm (default 90°)
    """
    function run_johnston2023(dataroot, outdir;
                              snr_threshold   = 5.0,
                              max_pa_err_deg  = 10.0,
                              pa_shift_deg    = 90.0)
        mkpath(outdir)

        # Files are flat in dataroot: JXXXX-XXXX_DATE_....ar
        # Group by pulsar name (everything before the first '_')
        all_ar = sort(Glob.glob("*.ar", dataroot))
        isempty(all_ar) && (@warn "No .ar files found in $dataroot"; return)

        psr_to_file = Dict{String,String}()
        for f in all_ar
            base = basename(f)
            psr  = split(base, "_")[1]
            haskey(psr_to_file, psr) || (psr_to_file[psr] = f)
        end

        for psr in sort(collect(keys(psr_to_file)))
            ar_file    = psr_to_file[psr]
            psr_outdir = joinpath(outdir, psr * "_rvm")
            mkpath(psr_outdir)

            println("\n=== $psr ===")
            println("  ar: $(basename(ar_file))")

            # --- convert to ASCII (cached) ---
            txt_file = joinpath(psr_outdir, "$(psr).txt")
            if !isfile(txt_file)
                try
                    Data.convert_psrfit_ascii(ar_files[1], txt_file)
                catch e
                    @warn "  convert failed for $(ar_files[1]): $e"
                    continue
                end
            end

            # --- load ---
            data4 = try
                Data.load_ascii_all(txt_file)
            catch e
                @warn "  load failed for $txt_file: $e"
                continue
            end
            n_pulses, n_bins, _ = size(data4)
            println("  n_pulses=$n_pulses  n_bins=$n_bins")

            # --- average profiles ---
            I_avg = vec(mean(data4[:, :, 1], dims=1))
            Q_avg = vec(mean(data4[:, :, 2], dims=1))
            U_avg = vec(mean(data4[:, :, 3], dims=1))
            V_avg = vec(mean(data4[:, :, 4], dims=1))
            L_avg = sqrt.(Q_avg .^ 2 .+ U_avg .^ 2)

            # --- on-pulse detection ---
            I_max = maximum(abs.(I_avg))
            on_mask = I_avg .> 0.1 * I_max

            if sum(on_mask) > 0.5 * n_bins
                # broad profile: use ±8% window around peak
                pk     = argmax(I_avg)
                half   = round(Int, 0.08 * n_bins)
                bin_st = max(1, pk - half)
                bin_end = min(n_bins, pk + half)
            else
                idxs   = findall(on_mask)
                bin_st  = first(idxs)
                bin_end = last(idxs)
            end
            println("  on-pulse bins: $bin_st:$bin_end")

            # --- off-pulse noise from Q/U ---
            off_rng = vcat(1:(bin_st - 1), (bin_end + 1):n_bins)
            if length(off_rng) < 10
                @warn "  Too few off-pulse bins for $psr — skipping"
                continue
            end
            sigma_Q   = std(vec(data4[:, off_rng, 2]))
            sigma_U   = std(vec(data4[:, off_rng, 3]))
            sigma_avg = (sigma_Q + sigma_U) / 2.0 / sqrt(n_pulses)
            thresh    = snr_threshold * sigma_avg
            println("  σ_avg=$(round(sigma_avg, sigdigits=3))  thresh=$(round(thresh, sigdigits=3))")

            # --- on-pulse slices & longitude grid ---
            on_rng  = bin_st:bin_end
            n_on    = length(on_rng)
            lon_full = collect(range(-180.0, 180.0, length=n_bins + 1)[1:end-1])
            lon_on  = lon_full[on_rng]

            I_on = I_avg[on_rng]
            L_on = L_avg[on_rng]
            V_on = V_avg[on_rng]
            Q_on = Q_avg[on_rng]
            U_on = U_avg[on_rng]

            # --- PA + error ---
            pa_raw  = 0.5 .* atan.(U_on, Q_on) .* (180.0 / π)
            pa_err  = fill(NaN, n_on)
            for i in 1:n_on
                if L_on[i] > thresh
                    e = 0.5 * sigma_avg / L_on[i] * (180.0 / π)
                    if e <= max_pa_err_deg
                        pa_err[i] = e
                    else
                        pa_raw[i] = NaN
                    end
                else
                    pa_raw[i] = NaN
                end
            end

            # --- remove OPM jumps ---
            pa_deopm, n_flipped = Data.deopm_pa(pa_raw)
            println("  deopm: $n_flipped bins flipped")

            # --- optional global OPM shift ---
            pa_plot = if pa_shift_deg != 0.0
                map(x -> isnan(x) ? NaN : mod(x + pa_shift_deg + 90.0, 180.0) - 90.0, pa_deopm)
            else
                pa_deopm
            end

            n_valid = sum(.!isnan.(pa_plot))
            println("  valid PA bins: $n_valid")
            if n_valid < 5
                @warn "  Too few valid PA points for $psr — skipping RVM"
                continue
            end

            # --- two-pass RVM fit ---
            res1 = Data.fit_rvm(lon_on, pa_plot, pa_err)
            if isnothing(res1)
                @warn "  fit_rvm pass1 returned nothing for $psr"
                continue
            end
            println("  pass1: α=$(round(res1.alpha,digits=1))°  β=$(round(res1.beta,digits=1))°  χ²ᵣ=$(round(res1.chi2_red,digits=2))")

            scale      = max(1.0, sqrt(res1.chi2_red))
            pa_err_sc  = map(x -> isnan(x) ? NaN : x * scale, pa_err)

            res2, chi2_map, alphas_deg, betas_deg =
                Data.fit_rvm(lon_on, pa_plot, pa_err_sc; return_map=true)
            if isnothing(res2)
                @warn "  fit_rvm pass2 returned nothing for $psr"
                continue
            end
            println("  pass2: α=$(round(res2.alpha,digits=1))°  β=$(round(res2.beta,digits=1))°  φ₀=$(round(res2.phi0,digits=1))°  PA₀=$(round(res2.pa0,digits=1))°  χ²ᵣ=$(round(res2.chi2_red,digits=2))")

            # --- RVM curve ---
            lon_rvm, pa_rvm, pa_rvm_ortho =
                Data.rvm_curve(res2, lon_on[1], lon_on[end])

            # --- plots ---
            # position_angle expects two bands; pass single-freq data for both
            Plot.position_angle(
                lon_on, pa_plot, pa_err,
                I_on, L_on, V_on,
                lon_on, pa_plot, pa_err,
                I_on, L_on, V_on,
                psr_outdir;
                show_        = false,
                lon_rvm_l    = lon_rvm,
                pa_rvm_l     = pa_rvm,
                pa_rvm_l_ortho = pa_rvm_ortho,
                phi0_l       = res2.phi0
            )

            Plot.geometry(
                chi2_map, chi2_map,
                alphas_deg, betas_deg,
                psr_outdir;
                show_    = false,
                name_mod = psr,
                ndof_l   = res2.ndof,
                ndof_h   = res2.ndof
            )

            PyPlot.close("all")
            println("  Saved to: $psr_outdir")
        end
    end


end # module

SpaTs.run_johnston2023("/home/psr/data/posselt/ar_files/", "/home/psr/output/")

println("Bye")