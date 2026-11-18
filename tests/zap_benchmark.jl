# Benchmark metod wycinania na danych SYNTETYCZNYCH ze znaną prawdą (żadne prawdziwe dane
# nie są czytane ani zmieniane). Uruchomienie:  julia tests/zap_benchmark.jl
# Kolumny: śmieci/nulle = ile wykryto, fałszywe = ile dobrych impulsów wycięto,
# vs ideał = S/N profilu względem wycięcia dokładnie prawdziwych śmieci i nulli.
# Metody: automat (z bramką / bez), TPA-like (off-pulse rms, 1.5 IQR),
# clfd-like (cechy off-pulse std/ptp/|ACF1|, Tukey q=2), wyrocznia (prawda).
include(joinpath(@__DIR__, "zap_synth.jl"))

methods = [
    ("automat", d -> Base.invokelatest(D.auto_zap, d, first(ON), last(ON)).bad),
    ("automat, nowy bezp.", d -> Base.invokelatest(D.auto_zap, d, first(ON), last(ON); rfi_greedy=true).bad),
    ("automat v6 (bez marg.)", d -> Base.invokelatest(D.auto_zap, d, first(ON), last(ON); rfi_greedy=true, null_margin=false).bad),
    ("automat bez RFI", d -> Base.invokelatest(D.auto_zap, d, first(ON), last(ON); rfi=false).bad),
    ("TPA-like (rms)", tpa_like),
    ("clfd-like (cechy)", clfd_like),
]
seeds = 1:4
snr(d, rows) = Base.invokelatest(D.profile_snr, d, first(ON), last(ON); drop=rows)

agg = Dict{String,Vector{Float64}}()
for (sc, lab) in scen
    println("\n=== $lab ===")
    @printf("%-20s %8s %8s %9s %9s %9s\n", "metoda", "śmieci", "nulle", "fałszywe", "S/N zysk", "vs ideał")
    res = Dict(m => zeros(5) for (m, _) in methods)
    for s in seeds
        d, junk, nulls = gen(sc, s)
        N = size(d, 1)
        clean = setdiff(1:N, junk, nulls)
        s0 = snr(d, Int[]); sid = snr(d, vcat(junk, nulls))
        for (m, f) in methods
            c = Set(f(d))
            jr = isempty(junk) ? NaN : count(in(c), junk) / length(junk)
            nr = isempty(nulls) ? NaN : count(in(c), nulls) / length(nulls)
            fc = count(in(c), clean) / length(clean)
            s1 = snr(d, collect(c))
            res[m] .+= [jr, nr, fc, 100(s1 / s0 - 1), 100(s1 / sid)] ./ length(seeds)
        end
    end
    for (m, _) in methods
        r = res[m]
        f(x) = isnan(x) ? "     -" : @sprintf("%5.0f%%", 100x)
        @printf("%-20s %8s %8s %8.1f%% %+8.1f%% %8.1f%%\n", m, f(r[1]), f(r[2]), 100r[3], r[4], r[5])
        push!(get!(agg, m, Float64[]), r[3], r[5])
    end
end
println("\n=== ŚREDNIO po scenariuszach ===")
for (m, _) in methods
    v = agg[m]; fc = v[1:2:end]; vi = v[2:2:end]
    @printf("%-20s  fałszywe cięcia %.1f%%   S/N jako %% ideału %.1f%%\n", m, 100mean(fc), mean(vi))
end

# RFI-only data (nulls kept = zaps2_rfi, for nulling studies) must not lose S/N
# against no cuts; old safety (rfi_ref=:cut) measured it only after the nulls
println("
=== Dane tylko po RFI (nulle zostawione) vs bez wycięć: S/N zmiana ===")
@printf("%-34s %12s %12s
", "scenariusz", "stary (:cut)", "nowy (:both)")
for (sc, lab) in scen
    sc in (:nulls, :faded, :mixed, :bright) || continue
    r = zeros(2)
    for s in seeds
        d, _, _ = gen(sc, s)
        s0 = snr(d, Int[])
        for (k, ref) in enumerate((:cut, :both))
            a = Base.invokelatest(D.auto_zap, d, first(ON), last(ON); rfi_ref=ref)
            r[k] += 100(snr(d, Base.invokelatest(D.split_zaps, a.reason).rfi) / s0 - 1) / length(seeds)
        end
    end
    @printf("%-34s %+11.2f%% %+11.2f%%
", lab, r[1], r[2])
end
