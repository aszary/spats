# Test wyboru reguly wycinania nulli na danych SYNTETYCZNYCH ze znana prawda (krok 1 z 2).
# Generuje impulsy (nulle o roznych dlugosciach, slaba i silnie zmienna emisja), uruchamia automat
# (z marginesem przy nullach i bez), liczy energie jak tools/final_v6.jl i zapisuje prawde.
#   julia tests/null_rule_gen.jl <katalog>
#   python3 tests/null_rule_eval.py <katalog>      (krok 2: model energii + porownanie wariantow)
using Statistics, Random, Printf
include(joinpath(@__DIR__, "zap_synth.jl"))
D.RFI_GREEDY[] = true

out = ARGS[1]; mkpath(joinpath(out, "energies"))
# scenariusz => (N, amplituda (S/N impulsu), sigma lognormalna, NF, dlugosci nulli)
SC = [("jasny_nuller", (1000, 3.0, 0.4, 0.30, 1:60)), ("krotkie_nulle", (1000, 3.0, 0.4, 0.15, 1:3)),
      ("slaby_nuller", (1000, 0.8, 0.4, 0.40, 5:60)), ("modulowany", (1000, 2.0, 0.9, 0.00, 1:1)),
      ("modul_nuller", (1000, 2.0, 0.9, 0.25, 1:40)), ("bardzo_slaby", (1000, 0.5, 0.4, 0.30, 5:60))]
nullpart(r) = begin                     # nulle + wygaszenia + margines przy nullach (to, co automat tnie jako "null")
    m = falses(length(r))
    for i in eachindex(r)
        r[i] in (3, 6) && (m[i] = true)
        if r[i] == 5
            nb = [r[j] for j in max(1, i - 1):min(length(r), i + 1) if j != i]
            any(x -> x in (3, 6), nb) && (m[i] = true)
        end
    end
    m
end
for (name, (N, A, sig, nf, lens)) in SC, seed in 1:3
    rng = MersenneTwister(100seed + length(name))
    isnull = falses(N)
    while count(isnull) < nf * N
        s = rand(rng, 1:N); L = rand(rng, lens)
        isnull[s:min(N, s + L - 1)] .= true
    end
    amp = A .* exp.(sig .* randn(rng, N)); amp[isnull] .= 0.0
    d = randn(rng, N, NB)
    for i in 1:N, j in 440:600
        d[i, j] += 2amp[i] * prof(j)
    end
    a = Base.invokelatest(D.auto_zap, d, first(ON), last(ON); null_margin=true)
    b = Base.invokelatest(D.auto_zap, d, first(ON), last(ON); null_margin=false)
    lab = zeros(Int, N); lab[a.reason .== 3] .= 2; lab[(a.reason .== 5) .| (a.reason .== 6)] .= 3
    lab[(a.reason .== 1) .| (a.reason .== 2) .| (a.reason .== 4)] .= 1
    on = collect(ON); c = (first(ON) + last(ON)) / 2
    dist(j) = min(abs(j - c), NB - abs(j - c))
    off = setdiff(1:NB, (first(ON) - 40):(last(ON) + 40))
    W = sort(sort(off, by=j -> -dist(j))[1:length(on)])
    Wx = Set(mod1(j + k, NB) for j in W for k in -10:10)
    B = filter(j -> !(j in Wx), off)
    base = [mean(@view d[i, B]) for i in 1:N]
    Eon = [sum(@view d[i, on]) - length(on) * base[i] for i in 1:N]
    Eoff = [sum(@view d[i, W]) - length(on) * base[i] for i in 1:N]
    σb = std(vec(d[:, W] .- base))
    psr = "J$(name)_$seed"
    open(joinpath(out, "energies", psr * ".csv"), "w") do io
        println(io, "# pulsar=$psr\n# N=$N\n# n_on=$(length(on))\n# n_B=$(length(B))\n# sigma_bin=$σb\n# bin_st=$(first(ON))\n# bin_end=$(last(ON))\n# W=$(first(W))-$(last(W))\ni,E_on,E_off,label")
        for i in 1:N; println(io, "$(i-1),$(Eon[i]),$(Eoff[i]),$(lab[i])"); end
    end
    ma, mb = nullpart(a.reason), nullpart(b.reason)
    open(joinpath(out, psr * ".truth"), "w") do io      # null, amplituda/srednia, automat z marg., automat bez marg.
        for i in 1:N; println(io, "$(Int(isnull[i])),$(amp[i] / A),$(Int(ma[i])),$(Int(mb[i]))"); end
    end
    @printf("%-22s NF prawdziwe %.3f  automat: z marg. %d, bez marg. %d\n", psr, mean(isnull), count(ma), count(mb))
end
