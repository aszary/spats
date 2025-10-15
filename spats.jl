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
        p = Data.process_psrdata(indir, outdir)
        Data.plot_psrdata(outdir, p)
        #folded = Data.load_ascii(outdir*"/pulsar.debase.p3fold")
        #Plot.p3fold(folded, outdir; start=3, bin_st=p["bin_st"]-20, bin_end=p["bin_end"]+20, name_mod="test", show_=true, repeat_num=4)
        
    end


    function main()
        # output directory for VPM
        vpmout = "/home/psr/data/OUTPUT/maciej/"

        #test(vpmout)
        #test2(vpmout)
        #test3("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)
        #test3("/home/psr/data/new/J1744-1610/2020-06-12-20:00:13/", vpmout)
        #process_psrdata("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)
        #J1539_6322_Sard("/home/psr/data/new/J1539-6322/2020-04-11-23:52:16/", vpmout)

        #process_psrdata(vpmout*"J1921+2153", vpmout*"J1921+2153")
        process_psrdata(vpmout*"J1220-6318", vpmout*"J1220-6318")
        process_psrdata(vpmout*"J1904+1011", vpmout*"J1904+1011")
        process_psrdata(vpmout*"J1720-0212", vpmout*"J1720-0212")
        process_psrdata(vpmout*"J1655-3048", vpmout*"J1655-3048")
        process_psrdata(vpmout*"J0711+0931", vpmout*"J0711+0931")
        process_psrdata(vpmout*"J1850+0026", vpmout*"J1850+0026")
        process_psrdata(vpmout*"J1745-0129", vpmout*"J1745-0129")
        process_psrdata(vpmout*"J1512-5431", vpmout*"J1512-5431")
        process_psrdata(vpmout*"J1844-0433", vpmout*"J1844-0433")
        process_psrdata(vpmout*"J1137-6700", vpmout*"J1137-6700")
        process_psrdata(vpmout*"J1945-0040", vpmout*"J1945-0040")
        process_psrdata(vpmout*"J1809-0119", vpmout*"J1809-0119")
        process_psrdata(vpmout*"J0421-0345", vpmout*"J0421-0345")
        process_psrdata(vpmout*"J1534-4428", vpmout*"J1534-4428")
        process_psrdata(vpmout*"J1857+0057", vpmout*"J1857+0057")
        process_psrdata(vpmout*"J1750-3503", vpmout*"J1750-3503")
        process_psrdata(vpmout*"J1036-6559", vpmout*"J1036-6559")
        process_psrdata(vpmout*"J1627-5936", vpmout*"J1627-5936")
        process_psrdata(vpmout*"J1645+1012", vpmout*"J1645+1012")
        process_psrdata(vpmout*"J0812-3905", vpmout*"J0812-3905")
        process_psrdata(vpmout*"J1903+2225", vpmout*"J1903+2225")
        process_psrdata(vpmout*"J1035-6345", vpmout*"J1035-6345")
        process_psrdata(vpmout*"J1703-4442", vpmout*"J1703-4442")
        process_psrdata(vpmout*"J1547-5750", vpmout*"J1547-5750")
        process_psrdata(vpmout*"J1926-0652", vpmout*"J1926-0652")
        process_psrdata(vpmout*"J1927+0911", vpmout*"J1927+0911")
        process_psrdata(vpmout*"J1907+0731", vpmout*"J1907+0731")

        #Data.process_all_data(vpmout)
        #Data.combine_pngs_to_pdf(vpmout)
        #Data.combine_pngs(vpmout)
        
        #Data.remove_folders(vpmout)
        #Data.remove_notinteresting("input/pulsars_interesting.txt", vpmout)
    end

end # module

SpaTs.main()

println("Bye")