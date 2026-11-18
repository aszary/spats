# Czułość progów automatu (v5) na danych syntetycznych z benchmarku: julia tests/threshold_sensitivity.jl
using Statistics, Random, Printf
include(joinpath(@__DIR__, "zap_synth.jl"))
D.RFI_GREEDY[] = true
D.NULL_MARGIN[] = true

snr(d, rows) = Base.invokelatest(D.profile_snr, d, first(ON), last(ON); drop=rows)
scen = [:clean, :bright, :nulls, :faded, :noisy, :stripes, :weakstripes, :badstart, :impulse, :mixed, :rfimany, :slownoise]
data = [(gen(sc, s)) for sc in scen for s in 1:3]
function score(kw)
    fc = Float64[]; vi = Float64[]; jr = Float64[]; nr = Float64[]
    for (d, junk, nulls) in data
        N = size(d, 1); clean = setdiff(1:N, junk, nulls)
        c = Set(Base.invokelatest(D.auto_zap, d, first(ON), last(ON); kw...).bad)
        push!(fc, count(in(c), clean) / length(clean))
        push!(vi, 100 * snr(d, collect(c)) / snr(d, vcat(junk, nulls)))
        isempty(junk) || push!(jr, count(in(c), junk) / length(junk))
        isempty(nulls) || push!(nr, count(in(c), nulls) / length(nulls))
    end
    @sprintf("fałszywe %.2f%%  S/N %.1f%% ideału  RFI %.0f%%  nulle %.0f%%", 100mean(fc), mean(vi), 100mean(jr), 100mean(nr))
end
println(rpad("domyślne (v5)", 34), score((;)))
for (lab, kws) in [("próg RFI rfi_k (3σ)", [(rfi_k=2.5,), (rfi_k=3.5,), (rfi_k=4.0,)]),
                   ("bezpiecznik rfi_tol (1%)", [(rfi_tol=0.005,), (rfi_tol=0.02,)]),
                   ("wygaszenia min_run_fade (30)", [(min_run_fade=20,), (min_run_fade=40,)]),
                   ("czyste nulle null_faded_frac (0.15)", [(null_faded_frac=0.10,), (null_faded_frac=0.25,)]),
                   ("wykrycie impulsu zdet (2σ)", [(zdet=1.5,), (zdet=2.5,)])]
    println(lab)
    for kw in kws
        println("   ", rpad(string(kw), 31), score(kw))
    end
end
