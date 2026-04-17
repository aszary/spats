#!/usr/bin/env julia

# Fast parse check for plot.jl and tools.jl
import Base.Meta.parse

files_to_check = [
    "modules/plot.jl",
    "modules/tools.jl"
]

for fname in files_to_check
    println("Checking $fname...")
    try
        code = read(fname, String)
        Base.parse(code; raise=true)
        println("  ✓ OK: $fname parsed successfully")
    catch e
        println("  ✗ ERROR in $fname:")
        println("    $(e)")
    end
end
println("\nParse check complete!")
