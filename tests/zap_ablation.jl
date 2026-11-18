# Test v7: czy kazdy krok automatu jest potrzebny (dane SYNTETYCZNE ze znana prawda).
# Kazdy wariant wylacza jeden krok; porownanie z pelnym automatem v6 w kazdym scenariuszu.
#   julia tests/zap_ablation.jl [plik.json]
# Krok jest zbedny na danych syntetycznych, gdy bez niego w ZADNYM scenariuszu nie spada wykrywanie
# smieci ani nulli o > 2 pp, S/N wzgledem idealu o > 0.5 pp, a falszywe ciecia nie rosna o > 0.2 pp.
include(joinpath(@__DIR__, "zap_synth.jl"))
include(joinpath(@__DIR__, "..", "tools", "v7_steps.jl"))
D.RFI_GREEDY[] = true
# V7_EXTRA="(use_jump = false, min_keep = 0)": tylko v6 i ten jeden wariant ("uproszczony")
haskey(ENV, "V7_EXTRA") && (V7_STEPS = [V7_STEPS[1], ("uproszczony", eval(Meta.parse(ENV["V7_EXTRA"])))])
seeds = 1:4
snr(d, rows) = Base.invokelatest(D.profile_snr, d, first(ON), last(ON); drop=rows)
R = Dict{String,Dict{Symbol,Vector{Float64}}}()      # wariant => scenariusz => [smieci, nulle, falszywe, vs ideal]
data = Dict(sc => [gen(sc, s) for s in seeds] for (sc, _) in scen)
for (name, st) in V7_STEPS
    R[name] = Dict{Symbol,Vector{Float64}}()
    for (sc, _) in scen
        r = zeros(4); nj = nn = 0
        for (d, junk, nulls) in data[sc]
            N = size(d, 1)
            clean = setdiff(1:N, junk, nulls)
            c = Set(Base.invokelatest(D.auto_zap, d, first(ON), last(ON); _v7_kw(st)...).bad)
            isempty(junk) || (r[1] += count(in(c), junk) / length(junk); nj += 1)
            isempty(nulls) || (r[2] += count(in(c), nulls) / length(nulls); nn += 1)
            r[3] += count(in(c), clean) / length(clean) / length(seeds)
            r[4] += 100 * snr(d, collect(c)) / snr(d, vcat(junk, nulls)) / length(seeds)
        end
        r[1] = nj > 0 ? r[1] / nj : NaN; r[2] = nn > 0 ? r[2] / nn : NaN
        R[name][sc] = r
    end
    println("policzone: ", name); flush(stdout)
end
base = R[V7_STEPS[1][1]]
println("\n=== Dane syntetyczne: zmiana względem pełnego automatu (najgorszy scenariusz) ===")
@printf("%-40s %10s %10s %10s %10s  %s\n", "wariant", "śmieci", "nulle", "fałszywe", "vs ideał", "ocena")
verdict = Dict{String,Any}()
for (name, _) in V7_STEPS
    dj = dn = df = dv = 0.0; wj = wn = wv = ""
    for (sc, _) in scen
        b, r = base[sc], R[name][sc]
        x = isnan(b[1]) ? 0.0 : 100 * (r[1] - b[1]); x < dj && (dj = x; wj = string(sc))
        x = isnan(b[2]) ? 0.0 : 100 * (r[2] - b[2]); x < dn && (dn = x; wn = string(sc))
        df = max(df, 100 * (r[3] - b[3]))
        x = r[4] - b[4]; x < dv && (dv = x; wv = string(sc))
    end
    ok = dj >= -2 && dn >= -2 && df <= 0.2 && dv >= -0.5
    verdict[name] = Dict("smieci_pp" => dj, "nulle_pp" => dn, "falszywe_pp" => df, "vs_ideal_pp" => dv, "zbedny" => ok,
                         "gdzie" => join(filter(!isempty, [wj, wn, wv]), ","))
    @printf("%-40s %+9.1f %+9.1f %+9.2f %+9.2f   %s\n", name, dj, dn, df, dv,
            name == V7_STEPS[1][1] ? "-" : ok ? "zbędny na syntetycznych" : "potrzebny (" * verdict[name]["gdzie"] * ")")
end
println("\n=== Średnio po scenariuszach: fałszywe cięcia / S/N jako % ideału ===")
for (name, _) in V7_STEPS
    v = collect(values(R[name]))
    @printf("%-40s %6.2f%%  %6.2f%%\n", name, 100 * mean(x[3] for x in v), mean(x[4] for x in v))
end
if !isempty(ARGS)
    open(ARGS[1], "w") do io
        println(io, "{", join(["\"$k\": {" * join(["\"$a\": " * (b isa Bool ? string(b) : b isa String ? "\"$b\"" : @sprintf("%.4f", b))
                                               for (a, b) in v], ", ") * "}" for (k, v) in verdict], ",\n"), "}")
    end
end
