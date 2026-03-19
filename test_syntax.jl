#!/usr/bin/env julia

# Quick syntax check
try
    include("spats.jl")
    println("SUCCESS: Code compiled without errors!")
catch e
    println("ERROR: $(e)")
    showerror(stdout, e)
end
