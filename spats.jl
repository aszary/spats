module SpaTs
    using ArgParse

    include("modules/data.jl")
    include("modules/plot.jl")


    function mkieth()

        data = []
        push!(data, Data.load_ascii("/home/szary/work/J1651-4246/mkieth/c_1.txt"))
        push!(data, Data.load_ascii("/home/szary/work/J1651-4246/mkieth/c_101.txt"))
        push!(data, Data.load_ascii("/home/szary/work/J1651-4246/mkieth/c_201.txt"))
        #push!(data, Data.load_ascii("/home/szary/work/J1651-4246/mkieth/b_301.txt"))
        println(size(data))

        d = Array{Float64}(undef, 100 * (size(data)[1]-1) + size(data[end])[1], size(data[end])[2])

        for (i,da) in enumerate(data)
            pulses, bins = size(da)
            for j in 1:pulses
                for k in i:bins
                    d[100*(i-1)+j, k] = da[j, k]
                end
            end
        end
        println(size(d))
        Plot.single(d, args["outdir"]; darkness=0.5, number=nothing, name_mod="J1651-4246")
        Plot.lrfs(d, args["outdir"]; darkness=0.1, start=1, name_mod="J1651-4246", change_fftphase=false)

    end


    function main()
        args = parse_commandline()
        for (arg, val) in args
            println("  $arg  =>  $val")
        end


        data = Data.load_ascii("/home/szary/work/MeerTime/J1651/1700.txt")
        Plot.single(data, args["outdir"]; start=1, number=nothing, darkness=0.5, name_mod="J1651")
        Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1,  bin_st=300, bin_end=700, name_mod="J1651", change_fftphase=false)

        #data = Data.load_ascii("/fred/oz005/users/aszary/search/J1705-1906/2019-04-03T02:19:47/1444.5/grand.txt")
        #Plot.single(data, args["outdir"]; darkness=0.5, number=nothing, name_mod="J1705-1906")
        #Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1,  bin_st=400, bin_end=500, name_mod="J1705-1906", change_fftphase=false)

        # B0320
        #Data.convert_psrfit_ascii("/data/szary/B0320/add.ar.high4", "input/high4.txt")
        #Data.convert_psrfit_ascii("/data/szary/B0320/add.ar.low", "input/low.txt")
        #Data.convert_psrfit_ascii("/data/szary/B0320/add.ar.pazi", "input/full_full.txt")
        #Data.convert_psrfit_ascii("/data/szary/B0320/add.ar.mid", "input/mid.txt")
        #Data.convert_psrfit_ascii("/data/szary/B0320/core/rawvoltages/SAP0/BEAM0/B0320+39_L570031_SAP0_BEAM0.paz.fscr.pdmp.AR", "input/core/")
        #data = Data.load_ascii("input/low.txt")

        #data = Data.load_ascii("input/low.p3fold")
        #data2 = Data.load_ascii("input/mid.p3fold")
        #data3 = Data.load_ascii("input/high4.p3fold")
        #Plot.offset(data, data2, "output/"; bin_st=50, bin_end=80, name_mod="low_mid", darkness=1.0)
        #Plot.offset_points(data, data2, "output/"; bin_st=50, bin_end=80, name_mod="db_low_mid")
        #Plot.offset_points3(data, data2, data3, "output/"; bin_st=50, bin_end=80, name_mod="low_mid_high4")
        #Plot.crosscorplot(data, data2, "/home/szary/work/B0320/"; bin_st=50, bin_end=80, name_mod="low_high4")
        #Plot.offset(data, data2, args["outdir"]; bin_st=50, bin_end=80, name_mod="low_high4", darkness=1.0)

        #Plot.single0(data, args["outdir"]; number=250, bin_st=50, bin_end=100)
        #Plot.single0(data, args["outdir"])
        #Plot.single(data, args["outdir"]; darkness=0.5, number=nothing, name_mod="all")
        #Plot.single(data, args["outdir"]; start=1, number=17, bin_st=50, bin_end=80, name_mod="low_p3fold", darkness=1.0)
        #Plot.p3fold(data, args["outdir"]; start=1, bin_st=50, bin_end=80, name_mod="high4", darkness=1.0)
        #Plot.single(data, args["outdir"]; start=1, number=512, bin_st=30, bin_end=95, name_mod="p3fold")
        #Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1, number=512, bin_st=50, bin_end=80, name_mod="low")
        #Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1000, number=512, bin_st=40, bin_end=85, name_mod="2")
        #Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1250, number=512, bin_st=40, bin_end=85, name_mod="3")
        #Plot.p3_evolution(data, args["outdir"]; darkness=0.1, bin_st=40, bin_end=85, name_mod="1", number=128)
        #Plot.p3_evolution(data, args["outdir"]; darkness=0.1, bin_st=40, bin_end=85, name_mod="2", number=256)
    end

    function parse_commandline()
        s = ArgParseSettings()
        @add_arg_table s begin
            "--indir", "-i"
                help = "input directory"
                default = "input"
            "--outdir", "-o"
                help = "output directory"
                default = "output"
            "--plot", "-p"
                help = "plots to create"
                default = []
                nargs = '*'
        end
        return parse_args(s)
    end

end # module

#SpaTs.mkieth()
SpaTs.main()

println("Bye")
