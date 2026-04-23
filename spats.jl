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

        # PSR J0034-0721
        #process_psrdata_16("/home/psr/data/new/J0034-0721/2019-10-18-22:29:51/", vpmout*"J0034-0721_16")
        #Data.analyse_p3folds_16(vpmout*"J0034-0721_16", "norefine")
        #process_psrdata(vpmout*"J0034-0721", vpmout*"J0034-0721") # P. nice

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

        # PSR J1110-5637  # TODO HERE
        #process_psrdata_16("/home/psr/data/new/J1110-5637/2019-10-19-08:18:28/", vpmout*"J1110-5637_16")
        #Data.analyse_p3folds_16(vpmout*"J1110-5637_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1110-5637_16", "norefine"; n_comp=2)
        #Data.position_angle(vpmout*"J1110-5637_16")
        Data.geometry_analysis(vpmout*"J1110-5637_16")



        # PSR J1048-5832 # TODO HERE2
        #process_psrdata_16("/home/psr/data/new/J1048-5832/2020-08-29-13:00:11/", vpmout*"J1048-5832_16")
        #Data.analyse_p3folds_16(vpmout*"J1048-5832_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1048-5832_16", "norefine"; n_comp=2)
        #Data.position_angle(vpmout*"J1048-5832_16")
        #Data.geometry_analysis(vpmout*"J1048-5832_16")





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
        #process_psrdata(vpmout*"J1921+2153", vpmout*"J1921+2153") # not nice

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
    end

end # module

SpaTs.main()

println("Bye")