module SpaTs
    using ArgParse

    include("modules/data.jl")

    function main()
        args = parse_commandline()
        for (arg, val) in args
            println("  $arg  =>  $val")
        end
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
