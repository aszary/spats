# Wspolne dla testow syntetycznych (zap_benchmark.jl, zap_ablation.jl): wczytanie funkcji automatu
# z modules/data.jl i generator danych ze znana prawda.
using Statistics, Random, Printf
function load(file)
    src = replace(read(file, String), "\r\n" => "\n")
    m = Module()
    Core.eval(m, :(using Statistics, Random))
    Core.eval(m, :(RFI_GREEDY = Ref(false)))
    Core.eval(m, :(NULL_MARGIN = Ref(true)))
    Core.eval(m, :(NULL_MODE = Ref(:auto)))
    Core.eval(m, :(ZAP_STEPS = Ref{Any}(NamedTuple())))
    Core.eval(m, :(NULL_EXTRA_CUR = Ref{Any}(nothing)))
    for name in ("offpulse_bins", "pulse_signal_noise", "noise_jump_z", "auto_zap_core", "pulses_to_zaps",
                 "profile_snr", "profile_gate", "zap_snr", "rfi_rms", "rfi_rms_z", "split_zaps", "fake_nulls_to_noisy!", "auto_zap")
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
    sc == :slownoise && (d[300:700, :] .*= 1.5)         # J1110/J1057: noisier stretch, signal fine
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
    if sc == :rfimany                                    # J2136-1606: duzo RFI roznej sily + jasne impulsy
        idx = randperm(rng, N)                           # z lekko podwyzszonym szumem i dobrym sygnalem
        for i in idx[1:40]                               # (te NIE sa smieciami - ich wyciecie obniza S/N)
            d[i, :] .*= 1.25; d[i, 440:600] .+= [4amp[i] * prof(j) for j in 440:600]
        end
        for (k, i) in enumerate(idx[41:100])             # impulsowe RFI w off-pulse, od slabego do silnego
            j0 = rand(rng, [50:350; 650:950]); d[i, j0:j0+9] .+= 1.0 + 4.0 * k / 60; push!(junk, i)
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

scen = [(:clean, "czysty (nic do cięcia)"), (:bright, "jasny, silna modulacja (J1048)"),
        (:nulls, "nulle"), (:faded, "wygaszony odcinek (J1110)"), (:noisy, "zaszumiony odcinek"),
        (:stripes, "paski RFI"), (:weakstripes, "słabe paski, słaby pulsar (J1750)"),
        (:badstart, "zły początek"), (:impulse, "impulsowe RFI w off-pulse"),
        (:mixed, "mieszany realistyczny"),
        (:rfimany, "dużo RFI różnej siły + jasne zaszumione impulsy (J2136)"),
        (:slownoise, "dłuższy odcinek z wyższym szumem, dobry sygnał (J1110/J1057)")]
