# Benchmark metod wycinania na danych SYNTETYCZNYCH ze znaną prawdą (żadne prawdziwe dane
# nie są czytane ani zmieniane). Uruchomienie:  julia tests/zap_benchmark.jl
# Kolumny: śmieci/nulle = ile wykryto, fałszywe = ile dobrych impulsów wycięto,
# vs ideał = S/N profilu względem wycięcia dokładnie prawdziwych śmieci i nulli.
# Metody: automat (z bramką / bez), TPA-like (off-pulse rms, 1.5 IQR),
# clfd-like (cechy off-pulse std/ptp/|ACF1|, Tukey q=2), wyrocznia (prawda).
using Statistics, Random, Printf
function load(file)
    src = replace(read(file, String), "\r\n" => "\n")
    m = Module()
    Core.eval(m, :(using Statistics, Random))
    for name in ("offpulse_bins", "pulse_signal_noise", "noise_jump_z", "auto_zap_core", "pulses_to_zaps",
                 "profile_snr", "profile_gate", "zap_snr", "rfi_rms", "auto_zap")
        r = match(Regex("    function " * name * "\\(.*?\\n    end\\n", "s"), src)
        isnothing(r) || Core.eval(m, Meta.parse(r.match))
    end
    m
end
const D = load(joinpath(@__DIR__, "..", "modules", "data.jl"))

const NB, ON = 1024, 456:573
prof(j) = exp(-((j - 497) / 6)^2) + 1.3 * exp(-((j - 527) / 8)^2)

# returns data, junk (RFI/noise – ALWAYS should go), nulls (no emission)
function gen(sc, seed; N=1000)
    rng = MersenneTwister(seed)
    d = randn(rng, N, NB)
    junk, nulls = Int[], Int[]
    amp = fill(1.0, N)
    base = sc in (:weakstripes,) ? 0.4 : sc in (:bright,) ? 3.0 : 1.0
    if sc == :bright                                   # J1048-like correlated strong modulation
        x = 0.0
        for i in 1:N; x = 0.85x + sqrt(1 - 0.85^2) * randn(rng); amp[i] = base * exp(1.0x); end
    else
        amp .= base .* exp.(0.5 .* randn(rng, N))
    end
    if sc in (:nulls, :mixed)
        for _ in 1:5
            s = rand(rng, 30:N-100); L = rand(rng, 20:80)
            append!(nulls, s:s+L-1)
        end
    end
    sc == :faded && append!(nulls, 300:599)
    nulls = sort(unique(nulls))
    for i in 1:N
        a = i in nulls ? (sc == :faded ? 0.08 * amp[i] : 0.0) : amp[i]
        for j in 440:600; d[i, j] += 2a * prof(j); end
    end
    if sc in (:noisy, :mixed)
        r = 700:759; d[r, :] .*= 2.5; append!(junk, r)
    end
    if sc in (:stripes, :mixed, :weakstripes)
        A = sc == :weakstripes ? 0.6 : 1.5
        for i in randperm(rng, N)[1:(sc == :weakstripes ? 30 : 15)]
            i in junk && continue
            d[i, :] .+= A .* sin.(2pi .* (1:NB) ./ rand(rng, 60:300) .+ 6rand(rng)); push!(junk, i)
        end
    end
    if sc in (:badstart, :mixed)
        d[1:20, :] .+= 3 .* randn(rng, 20, NB); append!(junk, 1:20)
    end
    if sc == :impulse
        for i in randperm(rng, N)[1:15]
            j0 = rand(rng, [50:350; 650:950]); d[i, j0:j0+4] .+= 8; push!(junk, i)
        end
    end
    if sc == :mixed                                      # correlated baseline everywhere
        for i in 1:N; d[i, :] .+= 0.5 .* sin.(2pi .* (1:NB) ./ rand(rng, 80:400) .+ 6rand(rng)); end
    end
    junk = sort(unique(junk))
    return d, junk, setdiff(nulls, junk)
end

# --- reference methods from the RFI survey (whole-pulse, off-pulse only) ---
offb() = vcat(1:first(ON)-41, last(ON)+41:NB)
function tpa_like(d)                                    # Keith+2026: off-pulse rms > Q3 + 1.5 IQR
    o = offb(); r = [std(@view d[i, o]) for i in axes(d, 1)]
    q1, q3 = quantile(r, 0.25), quantile(r, 0.75)
    findall(r .> q3 + 1.5 * (q3 - q1))
end
function clfd_like(d; q=2.0)                            # Morello+2019 profile masking, off-pulse features
    o = offb(); N = size(d, 1)
    F = zeros(N, 3)
    for i in 1:N
        x = d[i, o] .- mean(d[i, o]); v = sum(abs2, x)
        F[i, :] = [std(x), maximum(x) - minimum(x), v > 0 ? abs(sum(x[1:end-1] .* x[2:end]) / v) : 0.0]
    end
    bad = falses(N)
    for k in 1:3
        q1, q3 = quantile(F[:, k], 0.25), quantile(F[:, k], 0.75)
        bad .|= F[:, k] .> q3 + q * (q3 - q1)
    end
    findall(bad)
end

methods = [
    ("automat", d -> Base.invokelatest(D.auto_zap, d, first(ON), last(ON)).bad),
    ("automat bez RFI", d -> Base.invokelatest(D.auto_zap, d, first(ON), last(ON); rfi=false).bad),
    ("TPA-like (rms)", tpa_like),
    ("clfd-like (cechy)", clfd_like),
]
scen = [(:clean, "czysty (nic do cięcia)"), (:bright, "jasny, silna modulacja (J1048)"),
        (:nulls, "nulle"), (:faded, "wygaszony odcinek (J1110)"), (:noisy, "zaszumiony odcinek"),
        (:stripes, "paski RFI"), (:weakstripes, "słabe paski, słaby pulsar (J1750)"),
        (:badstart, "zły początek"), (:impulse, "impulsowe RFI w off-pulse"),
        (:mixed, "mieszany realistyczny")]
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
