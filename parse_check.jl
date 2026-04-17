src_plot = read("modules/plot.jl", String)
ex_plot = Meta.parse(src_plot; raise=true)
println("OK: modules/plot.jl parsed")

src_tools = read("modules/tools.jl", String)
ex_tools = Meta.parse(src_tools; raise=true)
println("OK: modules/tools.jl parsed")
