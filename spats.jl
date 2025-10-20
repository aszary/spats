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
        #vpmout = "/home/psr/output/"
        vpmout = "/home/psr/data/OUTPUT/maciej/"
        #test(vpmout)
        #test2(vpmout)
        #test3("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)
        #test3("/home/psr/data/new/J1744-1610/2020-06-12-20:00:13/", vpmout)
        #process_psrdata("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)
        #J1539_6322_Sard("/home/psr/data/new/J1539-6322/2020-04-11-23:52:16/", vpmout)

        #process_psrdata(vpmout*"J1921+2153", vpmout*"J1921+2153") # not nice
        #process_psrdata(vpmout*"J0630-2834", vpmout*"J0630-2834") # not nice
        #process_psrdata(vpmout*"J2139+2242", vpmout*"J2139+2242") # moderate
        #process_psrdata(vpmout*"J0820-1350", vpmout*"J0820-1350") # nice

        #process_psrdata(vpmout*"J1834-0010", vpmout*"J1834-0010") # P. nice
        #process_psrdata(vpmout*"J0034-0721", vpmout*"J0034-0721") # P. nice
        #process_psrdata(vpmout*"J0823+0159", vpmout*"J0823+0159") #  nice single pulses, but p3folds no
        #process_psrdata(vpmout*"J1842-0359", vpmout*"J1842-0359") #  nice single pulses, nice P. p3fold
        #process_psrdata(vpmout*"J1034-3224", vpmout*"J1034-3224") # single bad, p3folds bad, P3 stable
        #process_psrdata(vpmout*"J1133-6250", vpmout*"J1133-6250") # single not stable, P. p3fold (10 ybins)
        #process_psrdata(vpmout*"J1539-6322", vpmout*"J1539-6322") # single not stable, P. p3fold (10 ybins)
        #process_psrdata(vpmout*"J1555-3134", vpmout*"J1555-3134") # single not stable, P. p3fold (10 ybins)
        """
        process_psrdata(vpmout*"J1414-6802", vpmout*"J1414-6802") 
        
        process_psrdata(vpmout*"J0108-1431", vpmout*"J0108-1431") 
        process_psrdata(vpmout*"J0856-6137", vpmout*"J0856-6137") 
        process_psrdata(vpmout*"J2046+1540", vpmout*"J2046+1540")
        process_psrdata(vpmout*"J0904-7459", vpmout*"J0904-7459") 
        process_psrdata(vpmout*"J1914+0219", vpmout*"J1914+0219") 
        process_psrdata(vpmout*"J1902+0556", vpmout*"J1902+0556") 
        process_psrdata(vpmout*"J1514-4834", vpmout*"J1514-4834")
        process_psrdata(vpmout*"J2053-7200", vpmout*"J2053-7200") 
        process_psrdata(vpmout*"J1404+1159", vpmout*"J1404+1159") 
        process_psrdata(vpmout*"J1900-7951", vpmout*"J1900-7951") 
        process_psrdata(vpmout*"J0151-0635", vpmout*"J0151-0635") 
        process_psrdata(vpmout*"J1246+2253", vpmout*"J1246+2253") 
        process_psrdata(vpmout*"J0533+0402", vpmout*"J0533+0402") 
        process_psrdata(vpmout*"J1232-4742", vpmout*"J1232-4742") 
        process_psrdata(vpmout*"J1910+0714", vpmout*"J1910+0714") 
        process_psrdata(vpmout*"J1720-2933", vpmout*"J1720-2933") 
        process_psrdata(vpmout*"J1345-6115", vpmout*"J1345-6115") 
        process_psrdata(vpmout*"J2253+1516", vpmout*"J2253+1516") 
        process_psrdata(vpmout*"J1919+0134", vpmout*"J1919+0134") 
        process_psrdata(vpmout*"J0459-0210", vpmout*"J0459-0210") 
        process_psrdata(vpmout*"J0959-4809", vpmout*"J0959-4809") 
        process_psrdata(vpmout*"J1312-5402", vpmout*"J1312-5402") 
        process_psrdata(vpmout*"J1900-0933", vpmout*"J1900-0933") 
        process_psrdata(vpmout*"J1921+1948", vpmout*"J1921+1948") 
        process_psrdata(vpmout*"J1519-6106", vpmout*"J1519-6106") 
        process_psrdata(vpmout*"J1951+1123", vpmout*"J1951+1123") 
        process_psrdata(vpmout*"J1651-7642", vpmout*"J1651-7642") 
        process_psrdata(vpmout*"J1821+1715", vpmout*"J1821+1715") 
        process_psrdata(vpmout*"J0952-3839", vpmout*"J0952-3839")
        process_psrdata(vpmout*"J1638-3815", vpmout*"J1638-3815") 
        """
        #process_psrdata(vpmout*"J0421-0345", vpmout*"J0421-0345") 

       #process_psrdata(vpmout*"J1741-0840", vpmout*"J1741-0840")
        #process_psrdata(vpmout*"J1651-5222", vpmout*"J1651-5222")
        #process_psrdata(vpmout*"J1901-0906", vpmout*"J1901-0906")


    




        #=
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
        =#
        process_psrdata(vpmout*"J1944+1755", vpmout*"J1944+1755")
        process_psrdata(vpmout*"J1727-2739", vpmout*"J1727-2739")
        process_psrdata(vpmout*"J0725-1635", vpmout*"J0725-1635")
        process_psrdata(vpmout*"J1807+0756", vpmout*"J1807+0756")
        process_psrdata(vpmout*"J1746+2245", vpmout*"J1746+2245")
        process_psrdata(vpmout*"J1927+1852", vpmout*"J1927+1852")
        process_psrdata(vpmout*"J1904-1224", vpmout*"J1904-1224")
        process_psrdata(vpmout*"J1220-6318", vpmout*"J1220-6318")
        #process_psrdata(vpmout*"J0421-0345", vpmout*"J0421-0345")

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
        process_psrdata(vpmout*"J1834-0010", vpmout*"J1834-0010")
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