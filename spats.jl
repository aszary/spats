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

        #test(vpmout)
        #test2(vpmout)
        #test3("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)
        #test3("/home/psr/data/new/J1744-1610/2020-06-12-20:00:13/", vpmout)
        #process_psrdata("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)
        #J1539_6322_Sard("/home/psr/data/new/J1539-6322/2020-04-11-23:52:16/", vpmout)

        # uhf data
        #process_psrdata("/home/psr/data/uhf/J0820-1350/2020-07-11-10:30:12/816/single/", vpmout*"J0820-1350_uhf")
        #process_psrdata("/home/psr/data/uhf/J0823+0159/2023-02-27-20:38:58/816/single", vpmout*"J0823+0159_uhf")
        # TODO work on PSR J1842-0359
        # L-band
        #process_psrdata_single("/home/psr/data/new/J1842-0359/2019-11-05-18:03:43/", vpmout*"J1842-0359_single")
        #process_psrdata_16("/home/psr/data/new/J1842-0359/2019-11-05-18:03:43/", vpmout*"J1842-0359_16")
        #Data.analyse_p3folds_16(vpmout*"J1842-0359_16", "refine")
        # uhf
        #process_psrdata("/home/psr/data/uhf/J1842-0359/2020-07-25-19:58:19/816/single/", vpmout*"J1842-0359_uhf")
        # uhf multifrequency
        #process_psrdata_16("/home/psr/data/uhf/J1842-0359/2020-07-25-19:58:19/816/single/", vpmout*"J1842-0359_uhf_16") 
        #Data.analyse_p3folds_16(vpmout*"J1842-0359_uhf_16", "norefine")

        # work on pulsar J2139+2242
        #process_psrdata_16("/home/psr/data/new/J2139+2242/2020-09-07-20:01:41/", vpmout*"J2139+2242_16")
        #Data.analyse_p3folds_16(vpmout*"J2139+2242_16", "norefine")

        # work on pulsar J1807+0756
        #process_psrdata_16("/home/psr/data/new/J1807+0756/2020-06-07-22:03:05/", vpmout*"J1807+0756_16")
        #Data.analyse_p3folds_16(vpmout*"J1807+0756_16", "refine")


        #process_psrdata(vpmout*"J2139+2242", vpmout*"J2139+2242") # moderate
        #process_psrdata(vpmout*"J1921+2153", vpmout*"J1921+2153") # not nice
        #process_psrdata(vpmout*"J0630-2834", vpmout*"J0630-2834") # not nice
        #process_psrdata(vpmout*"J0820-1350", vpmout*"J0820-1350") # nice

        #process_psrdata(vpmout*"J1834-0010", vpmout*"J1834-0010") # P. nice
        #process_psrdata(vpmout*"J0034-0721", vpmout*"J0034-0721") # P. nice
        #process_psrdata(vpmout*"J0823+0159", vpmout*"J0823+0159") #  nice single pulses, but p3folds no
        #process_psrdata(vpmout*"J1842-0359", vpmout*"J1842-0359") #  nice single pulses, nice P. p3fold
        #process_psrdata(vpmout*"J1034-3224", vpmout*"J1034-3224") # single bad, p3folds bad, P3 stable
        #process_psrdata(vpmout*"J1133-6250", vpmout*"J1133-6250") # single not stable, P. p3fold (10 ybins)
        #process_psrdata(vpmout*"J1539-6322", vpmout*"J1539-6322") # single not stable, P. p3fold (10 ybins)
        #process_psrdata(vpmout*"J1555-3134", vpmout*"J1555-3134") # single not stable, P. p3fold (10 ybins)
        #process_psrdata(vpmout*"J1907+0731", vpmout*"J1907+0731") #
        #process_psrdata(vpmout*"J0421-0345", vpmout*"J0421-0345") #
        #process_psrdata(vpmout*"J1514-4834", vpmout*"J1514-4834") #

        # checking multifrequency offsets
        #process_psrdata_16("/home/psr/data/new/J2139+2242/2020-09-07-20:01:41/", vpmout*"J2139+2242_16")
        #process_psrdata(vpmout*"J1414-6802", vpmout*"J1414-6802") # start diging there...
        #process_psrdata("/home/psr/data/new/J1414-6802/2020-07-24-16:17:58/", vpmout*"J1414-6802_test") # start diging there...
        #process_psrdata_16("/home/psr/data/new/J0820-1350/2020-01-11-01:05:56/", vpmout*"J0820-1350_16")
        #Data.analyse_p3folds_16(vpmout*"J0820-1350_16", "norefine")
        #process_psrdata_16("/home/psr/data/new/J1414-6802/2020-07-24-16:17:58/", vpmout*"J1414-6802_16") # start diging there...
        #process_psrdata_16("/home/psr/data/new/J1539-6322/2020-04-11-23:52:16/", vpmout*"J1539-6322_16") 
        #process_psrdata_16("/home/psr/data/new/J1514-4834/2019-12-16-08:22:16/", vpmout*"J1514-4834_16")
        #process_psrdata_16("/home/psr/data/new/J0034-0721/2019-10-18-22:29:51/", vpmout*"J0034-0721_16")
        #process_psrdata_16("/home/psr/data/new/J0823+0159/2019-11-05-23:33:48/", vpmout*"J0823+0159_16")
        #process_psrdata_16("/home/psr/data/new/J1842-0359/2019-11-05-18:03:43/", vpmout*"J1842-0359_16")
        #process_psrdata_16("/home/psr/data/new/J1034-3224/2019-10-19-07:31:26/", vpmout*"J1034-3224_16")
        #process_psrdata_16("/home/psr/data/new/J1133-6250/2019-11-06-00:46:43/", vpmout*"J1133-6250_16")
        #process_psrdata_16("/home/psr/data/new/J1555-3134/2019-10-31-06:34:39/", vpmout*"J1555-3134_16")
        #Data.analyse_p3folds_16(vpmout*"J1555-3134_16", "refine")
        #process_psrdata_16("/home/psr/data/new/J1312-5402/2020-03-30-19:11:06/", vpmout*"J1312-5402_16")

        # checking single pulses offset
        #process_psrdata_single("/home/psr/data/new/J2139+2242/2020-09-07-20:01:41/", vpmout*"J2139+2242_single") # no visible offset
        #process_psrdata_single("/home/psr/data/new/J0820-1350/2020-01-11-01:05:56/", vpmout*"J0820-1350_single") # very bad no signal at least at the beginning
        #process_psrdata_single("/home/psr/data/new/J0823+0159/2019-11-05-23:33:48/", vpmout*"J0823+0159_single")
        #process_psrdata_single("/home/psr/data/new/J1842-0359/2019-11-05-18:03:43/", vpmout*"J1842-0359_single")
        #process_psrdata_single("/home/psr/data/new/J1555-3134/2019-10-31-06:34:39/", vpmout*"J1555-3134_single") # no offset
        #process_psrdata_single("/home/psr/data/new/J1539-6322/2020-04-11-23:52:16/", vpmout*"J1539-6322_single") 
        #process_psrdata_single("/home/psr/data/new/J1414-6802/2020-07-24-16:17:58/", vpmout*"J1414-6802_single") # weak signal no offset
        #process_psrdata_single("/home/psr/data/new/J1133-6250/2019-11-06-00:46:43/", vpmout*"J1133-6250_single") # some offset? probably not
        #process_psrdata_16("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout*"J1919+0134_16") # visible offset
        #Data.analyse_p3folds_16(vpmout*"J1919+0134_16", "refine")
        #process_psrdata_16("/home/psr/data/new/J1312-5402/2020-03-30-19:11:06/", vpmout*"J1312-5402_16") # visible offset
        #Data.analyse_p3folds_16(vpmout*"J1312-5402_16", "refine")
        #process_psrdata_16("/home/psr/data/new/J1404+1159/2020-03-30-02:48:21/", vpmout*"J1404+1159_16") # visible offset
        #Data.analyse_p3folds_16(vpmout*"J1404+1159_16", "refine")
        #process_psrdata_16("/home/psr/data/new/J1133-6250/2019-11-06-00:46:43/", vpmout*"J1133-6250_16") # visible offset
        #Data.analyse_p3folds_16(vpmout*"J1133-6250_16", "norefine")
        #process_psrdata_16("/home/psr/data/new/J1414-6802/2020-07-24-16:17:58/", vpmout*"J1414-6802_16") # visible offset
        #Data.analyse_p3folds_16(vpmout*"J1414-6802_16", "norefine")
        #process_psrdata_16("/home/psr/data/new/J2046+1540/2019-11-26-17:18:47/", vpmout*"J2046+1540_16") # visible offset
        #Data.analyse_p3folds_16(vpmout*"J2046+1540_16", "norefine")
        process_psrdata_16("/home/psr/data/new/J1910+0714/2020-05-24-21:50:19/", vpmout*"J1910+0714_16") # visible offset
        Data.analyse_p3folds_16(vpmout*"J1910+0714_16", "norefine")


        #Data.process_all_data(vpmout)
        #Data.combine_pngs_to_pdf(vpmout)
        #Data.combine_pngs(vpmout)
        
        #Data.remove_folders(vpmout)
        #Data.remove_notinteresting("input/pulsars_interesting.txt", vpmout)
    end

end # module

SpaTs.main()

println("Bye")