#!/usr/bin/env julia
# Minimal direct test
include("modules/tools.jl")
include("modules/plot.jl")
include("modules/functions.jl")
include("modules/data.jl")

# Now test RVM
Base.include(Main.Plot, "modules/plot.jl")
result = Main.Plot.rvm(data4, outdir, name_mod; 
                       bin_st=450, bin_end=580, period=0.536, n_freq=2, show_=true)
