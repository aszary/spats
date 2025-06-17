module SpaTs
    using ArgParse
    using Glob
    using JSON

    include("modules/data.jl")
    include("modules/plot.jl")
    include("modules/tools.jl")


    function main()
        vpmout = "/home/psr/output/"
        #Data.process_all_data(vpmout)
        Data.plot_psrdata(vpmout)
        Data.combine_4page(vpmout)
        
    end

end # module

SpaTs.main()

println("Bye")