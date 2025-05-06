module SpaTs
    using ArgParse
    using Glob
    using JSON
    using Printf

    include("modules/data.jl")
    include("modules/plot.jl")
    include("modules/tools.jl")


    function test(outdir)
        d = Data.load_ascii("input/1.txt")
        Plot.single(d, outdir; darkness=0.3, number=256, bin_st=400, bin_end=600, start=1, name_mod="1", show_=true)
        Plot.average(d, outdir; number=256, bin_st=400, bin_end=600, start=1, name_mod="1", show_=true)
        Plot.lrfs(d, outdir; darkness=0.1, start=1, name_mod="1", bin_st=500, bin_end=530, show_=true)
    end

    
    function test2(outdir)    
        Data.convert_psrfit_ascii("/home/psr/data/new/J1232-4742/2020-04-11-22:21:18/2020-04-11-22:21:18_00000-00255.spCF", outdir*"1.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/new/J1232-4742/2020-04-11-22:21:18/2020-04-11-22:21:18_00256-00511.spCF", outdir*"2.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/new/J1232-4742/2020-04-11-22:21:18/2020-04-11-22:21:18_00512-00767.spCF", outdir*"3.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/new/J1232-4742/2020-04-11-22:21:18/2020-04-11-22:21:18_00768-00961.spCF", outdir*"4.txt")  
        data1 = Data.load_ascii(outdir*"1.txt")
        data2 = Data.load_ascii(outdir*"2.txt")
        data3 = Data.load_ascii(outdir*"3.txt")
        data4 = Data.load_ascii(outdir*"4.txt")
        data = vcat(data1, data2, data3, data4)
        Plot.single(data, outdir; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1232-4742", show_=true)
        Plot.average(data, outdir; number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1232-4742", show_=true)
        Plot.lrfs(data, outdir; darkness=0.1, start=1, name_mod="J1232-4742", bin_st=400, bin_end=600, show_=true)
        folded = Tools.p3fold(data, 20, 40)
        Plot.single(folded, outdir; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1232-4742_p3fold", show_=true)
    end



    function repuls(vpmout::String, num_files::Int)
        max_rms_value = 0.0
    
        for i in 1:num_files
            infile = "$(vpmout)$(i).txt" 

            for line in readlines(infile)
                parts = split(strip(line))  
    

                if length(parts) == 7
                    fifth_col = parse(Float64, parts[5])
                    sixth_col = parse(Float64, parts[6])
                    rms_value = sqrt(fifth_col^2 + sixth_col^2)


                    if rms_value > max_rms_value
                        max_rms_value = rms_value
                    end
                end
            end
        end
    
        threshold = 0.3 * max_rms_value
    
        for i in 1:num_files
            infile = "$(vpmout)$(i).txt" 
            outfile = "$(vpmout)$(i)_zmiany.txt"
    
            open(outfile, "w") do io
                for line in readlines(infile)
                    parts = split(strip(line))  
    
                    if length(parts) == 7
                        fifth_col = parse(Float64, parts[5])
                        sixth_col = parse(Float64, parts[6])
                        rms_value = sqrt(fifth_col^2 + sixth_col^2)
    
                        if rms_value < threshold
                            parts[4] = "0.0"  
                        end
    
                        println(io, join(parts, " "))
                    else
                        println(io, line)
                    end
                end
            end
        end
    
        data1 = Data.load_ascii(vpmout * "1_zmiany.txt")
        data2 = Data.load_ascii(vpmout * "2_zmiany.txt")
        data3 = Data.load_ascii(vpmout * "3_zmiany.txt")
        data4 = Data.load_ascii(vpmout * "4_zmiany.txt")
        
        combined_data = vcat(data1, data2, data3, data4)
    
        Plot.single(combined_data, vpmout; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1232-4742", show_=true)
        #Plot.average(combined_data, vpmout; number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319", show_=true)
        Plot.lrfs(combined_data, vpmout; darkness=0.1, start=1, name_mod="J1232-4742", bin_st=400, bin_end=600, show_=true)
        #folded = Tools.p3fold(combined_data, 20, 40)
        #Plot.single(folded, vpmout; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319_p3fold", show_=true)
    end
    
    
    
    
    
    
    
    
    



    function process_psrdata(indir, outdir)
        p = Data.process_psrdata(indir, outdir)
        folded = Data.load_ascii(outdir*"/pulsar.debase.p3fold")
        Plot.p3fold(folded, outdir; start=3, bin_st=p["bin_st"]-20, bin_end=p["bin_end"]+20, name_mod="test", show_=true, repeat_num=4)
        
    end


    function main()
        # output directory for VPM
        vpmout = "/home/psr/output/"
        num_files = 4  # lub inna liczba plików, które chcesz przetworzyć
        #repuls(vpmout, num_files)
        #test(vpmout)
        test2(vpmout)
        #process_psrdata("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)

    end

end # module

SpaTs.main()

println("Bye")