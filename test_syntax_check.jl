for f in ["modules/tools.jl", "modules/plot.jl"]
    src = read(f, String)
    ok = true
    for ex in Meta.parseall(src).args
        if ex isa Expr && ex.head == :error
            println("SYNTAX ERROR in $f: ", ex.args[1]); ok = false
        end
    end
    ok && println("OK: $f")
end
