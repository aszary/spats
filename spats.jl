module SpaTs
    using ArgParse

    include("modules/data.jl")
    include("modules/plot.jl")

    function main()
        args = parse_commandline()
        for (arg, val) in args
            println("  $arg  =>  $val")
        end

        #Data.convert_psrfit_ascii("/data/szary/B0320/add.ar.high", "input/high/")
        Data.convert_psrfit_ascii("/data/szary/B0320/add.ar.high4", "input/high4/")
        #Data.convert_psrfit_ascii("/data/szary/B0320/core/rawvoltages/SAP0/BEAM0/B0320+39_L570031_SAP0_BEAM0.paz.fscr.pdmp.AR", "input/core/")
        #data = Data.load_ascii("input/high2/pulse*")
        #Plot.single0(data, args["outdir"]; number=250, bin_st=50, bin_end=100)
        #Plot.single0(data, args["outdir"])
        #Plot.single(data, args["outdir"]; darkness=0.5, number=nothing, name_mod="all")
        #Plot.single(data, args["outdir"]; start=1, number=141, bin_st=30, bin_end=95, name_mod="core")
        #Plot.single(data, args["outdir"]; start=1, number=512, bin_st=30, bin_end=95, name_mod="high2")
        #Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1, number=512, bin_st=45, bin_end=85, name_mod="high2")
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

SpaTs.main()

println("Bye")
