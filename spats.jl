module SpaTs
    using ArgParse
    using Glob
    using JSON

    include("modules/data.jl")
    include("modules/plot.jl")
    include("modules/tools.jl")


    function J1319(outdir)    
        Data.convert_psrfit_ascii("/home/psr/data/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00000-00255.spCF", outdir*"1.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00256-00511.spCF", outdir*"2.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00512-00767.spCF", outdir*"3.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00768-01029.spCF", outdir*"4.txt")  
        data1 = Data.load_ascii(outdir*"1.txt")
        data2 = Data.load_ascii(outdir*"2.txt")
        data3 = Data.load_ascii(outdir*"3.txt")
        data4 = Data.load_ascii(outdir*"4.txt")
        data = vcat(data1, data2, data3, data4)
        Plot.single(data, outdir; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319", show_=true)
        Plot.average(data, outdir; number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319", show_=true)
        Plot.lrfs(data, outdir; darkness=0.1, start=1, name_mod="J1319", bin_st=400, bin_end=600, show_=true)
        #Plot.p3fold(data, outdir; darkness=0.1, start=1, name_mod="J1319", bin_st=400, bin_end=600, show_=true)
    end

    function J1319_psrchive(indir, outdir)
        Data.process_psrchive(indir, outdir, ["2019-12-15-03:19:04_00000-00255.spCF", "2019-12-15-03:19:04_00256-00511.spCF", "2019-12-15-03:19:04_00512-00767.spCF", "2019-12-15-03:19:04_00768-01029.spCF"], "J1319-6105.spCF")
    end

    function main()
        # output directory for VPM
        vpmout = "/home/psr/output/"

        J1319(vpmout)
        #J1319_psrchive("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/", vpmout)

    end

end # module

SpaTs.main()

println("Bye")