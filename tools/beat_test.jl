# Test hipotezy "systemu dudnien" (Hsu+2026, TPA XX) dla dwoch pikow w widmie fluktuacji, np. J1555-3134:
# jeden dryf D modulowany amplitudowo M (bez nosnej) daje piki boczne f1 = D - M i f2 = D + M.
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/beat_test.jl"); beat_test_run()
# Najpierw te same testy na danych syntetycznych ze znana odpowiedzia (dudnienia / dwa niezalezne dryfy /
# przelaczanie trybow), zeby bylo widac, co kazdy test rozroznia; potem pulsary. Wynik: OUTPUT/beat_test/WYNIK.txt
#
# Testy (okna W impulsow przesuwane co 1, faza wzgledem bezwzglednego numeru impulsu):
#  T1 czy moc obu pikow rosnie i maleje RAZEM (dudnienia: r > 0; przelaczanie trybow: r < 0; niezalezne: r ~ 0);
#     istotnosc: przesuniecia cykliczne jednego szeregu (zachowuja autokorelacje okien)
#  T3 nachylenie fazy obu pikow (czysta modulacja amplitudy: takie samo)
#  T4 ksztalt amplitudy w dlugosci (czysta modulacja amplitudy: ten sam ksztalt dla obu pikow)
#  T5 slabe piki przy M i D w widmie calej obserwacji (S/N wzgledem lokalnej mediany) - pomocniczo: przy
#     modulacji bez nosnej nie ma ich (tez w syntetycznych); pik przy D = czesc nosnej nie jest wygaszona
# (Test zmian czestotliwosci i sygnalu przy 2M usuniete: na danych syntetycznych wychodzily istotne takze dla
#  niezaleznych dryfow / przelaczania trybow, wiec niczego nie rozstrzygaja.)
# T1 rozroznia wszystkie trzy przypadki (syntetyczne: dudnienia r ~ +0.9, niezalezne ~ +0.1..0.2 nieistotne,
# tryby ~ -0.95). T3: rozne nachylenia wykluczaja dudnienia, takie same ich nie dowodza.
using Statistics, Random, Printf, LinearAlgebra

"Zespolone LRFS przy czestotliwosci f w oknach W (krok 1): L[b, j], okno b = impulsy b:b+W-1."
function _bt_windowed(X, f, W)
    N = size(X, 1)
    car = [cis(2π * f * (n - 1)) for n in 1:N]
    C = cumsum(X .* car, dims=1)
    nw = N - W + 1
    L = Matrix{ComplexF64}(undef, nw, size(X, 2))
    L[1, :] = C[W, :]
    for b in 2:nw; L[b, :] = C[b + W - 1, :] .- C[b - 1, :]; end
    return L
end

_rank(x) = (p = sortperm(x); r = similar(x, Float64); r[p] = 1:length(x); r)
_spearman(a, b) = cor(_rank(a), _rank(b))

"Korelacja Spearmana + p z przesuniec cyklicznych (minimalne przesuniecie W)."
function _bt_corr(a, b, W; nperm=1000, rng=MersenneTwister(1))
    ok = isfinite.(a) .& isfinite.(b); a = a[ok]; b = b[ok]; n = length(a)
    n < 3W && return (r=NaN, p=NaN)
    r = _spearman(a, b)
    null = [_spearman(a, circshift(b, rand(rng, W:n-W))) for _ in 1:nperm]
    return (r=r, p=(1 + count(x -> abs(x) >= abs(r), null)) / (nperm + 1))
end

"Odchylenie czestotliwosci w czasie: faza wspolrzednej projekcji L na sredni wzor, pochodna po oknach."
function _bt_freqdev(L, W)
    ref = vec(sum(L, dims=1)); ref ./= max(norm(ref), eps())
    g = L * conj.(ref)
    θ = angle.(g)
    for b in 2:length(θ)                                  # rozwiniecie
        d = θ[b] - θ[b-1]; θ[b] -= 2π * round(d / 2π)
    end
    # pochodna po oknach wygladzona na dlugosci W (odchylenie czestotliwosci w cyklach/P)
    df = [ (θ[min(end, b + W ÷ 2)] - θ[max(1, b - W ÷ 2)]) / (2π * (min(length(θ), b + W ÷ 2) - max(1, b - W ÷ 2))) for b in eachindex(θ)]
    return df, abs.(g)
end

"Nachylenie fazy (deg/bin) w kazdym oknie, faza wzgledem bezwzglednego impulsu; tylko istotne biny."
function _bt_slopes(L, Loff)
    nw = size(L, 1); s = fill(NaN, nw)
    for b in 1:nw
        σ = std(vcat(real.(Loff[b, :]), imag.(Loff[b, :])))
        js = findall(abs.(L[b, :]) .> 3σ)
        length(js) >= 4 || continue
        ph = rad2deg.(angle.(L[b, js]))
        for k in 2:length(ph); d = ph[k] - ph[k-1]; ph[k] -= 360 * round(d / 360); end
        x = Float64.(js)
        s[b] = sum((x .- mean(x)) .* (ph .- mean(ph))) / sum((x .- mean(x)) .^ 2)
    end
    return s
end

"Wszystkie testy dla danych X (impulsy x biny, zera = wyciete), okna bin_st:bin_end, piki f1 < f2."
function beat_test_core(X0, bin_st, bin_end, f1, f2; W=64, nperm=1000)
    N, nb = size(X0)
    good = [any(!iszero, @view X0[i, :]) for i in 1:N]
    X = Float64.(X0); m = mean(X[good, :], dims=1); X .-= m; X[.!good, :] .= 0
    on = bin_st:bin_end
    c = (bin_st + bin_end) / 2
    off = sort(sort(setdiff(1:nb, max(1, bin_st - 20):min(nb, bin_end + 20)), by=j -> -min(abs(j - c), nb - abs(j - c)))[1:min(length(on), nb - length(on) - 40)])
    Xon = X[:, on]; Xoff = X[:, off]
    D = (f1 + f2) / 2; M = (f2 - f1) / 2
    L1 = _bt_windowed(Xon, f1, W); L2 = _bt_windowed(Xon, f2, W)
    O1 = _bt_windowed(Xoff, f1, W); O2 = _bt_windowed(Xoff, f2, W)
    n1 = mean(sum(abs2.(O1), dims=2)) * length(on) / length(off); n2 = mean(sum(abs2.(O2), dims=2)) * length(on) / length(off)
    P1 = vec(sum(abs2.(L1), dims=2)) .- n1; P2 = vec(sum(abs2.(L2), dims=2)) .- n2
    t1 = _bt_corr(P1, P2, W; nperm=nperm)
    # T1n: moc wzgledem jasnosci w oknie (kwadrat energii impulsow) - wspolne zmiany jasnosci calego pulsara
    # (scyntylacja, nulle) nie daja wtedy pozornego "razem"
    Eraw = vec(sum(Float64.(X0[:, on]), dims=2)); cE = cumsum(Eraw)
    Ew = [b == 1 ? cE[W] : cE[b + W - 1] - cE[b - 1] for b in 1:N-W+1]
    t1n = _bt_corr(P1 ./ max.(Ew, eps()) .^ 2, P2 ./ max.(Ew, eps()) .^ 2, W; nperm=nperm)
    # T3 osobno w kazdej skladowej profilu (przez przerwe miedzy skladowymi faza moze zgubic pelny obrot -
    # J1901: F1 i F3 wychodzily 3x rozne, choc w TPA XX i na naszym wykresie V3 sa podobne); podzial w minimum
    # profilu miedzy dwoma najwyzszymi maksimami
    prof = vec(mean(Float64.(X0[good, on]), dims=1))
    ps = [mean(prof[max(1, j - 2):min(end, j + 2)]) for j in eachindex(prof)]
    mx = [j for j in 2:length(ps)-1 if ps[j] >= ps[j-1] && ps[j] >= ps[j+1] && ps[j] > 0.2maximum(ps)]
    comps = [1:length(on)]
    if length(mx) >= 2
        top = sort(mx, by=j -> -ps[j])
        a = top[1]; bq = findfirst(j -> abs(j - a) >= 10, top[2:end])
        if !isnothing(bq)
            b2 = top[1 + bq]; lo, hi = minmax(a, b2)
            cut = lo - 1 + argmin(ps[lo:hi])
            comps = [1:cut, cut+1:length(on)]
        end
    end
    neff = max(1, (N - W + 1) ÷ W)                         # niezalezne okna
    function sl(L, O, r)
        v = filter(isfinite, _bt_slopes(L[:, r], O))
        length(v) < 3 ? (NaN, NaN) : (mean(v), std(v) / sqrt(neff))
    end
    T3 = [(c=k, s1=sl(L1, O1, r), s2=sl(L2, O2, r)) for (k, r) in enumerate(comps)]
    zz = [(t.s1[1] - t.s2[1]) / sqrt(t.s1[2]^2 + t.s2[2]^2) for t in T3]
    zs = isempty(filter(isfinite, zz)) ? NaN : filter(isfinite, zz)[argmax(abs.(filter(isfinite, zz)))]
    m1, e1 = T3[1].s1; m2, e2 = T3[1].s2
    snrP = mean(P1 .+ P2) / (n1 + n2)                      # sila pikow wzgledem szumu (do dopasowania syntetycznych)
    A1 = vec(mean(abs.(L1), dims=1)); A2 = vec(mean(abs.(L2), dims=1))
    a1 = A1 ./ sum(A1); a2 = A2 ./ sum(A2)
    rA = cor(a1, a2)
    dA = maximum(abs.(a1 .- a2)) / maximum(a1 .+ a2) * 2
    # T5: widmo calej obserwacji (bloki 256), S/N przy M i D wzgledem lokalnej mediany (bez pikow)
    Wl = min(256, N)
    nbk = N ÷ Wl
    spec = zeros(Wl ÷ 2)
    for k in 1:nbk
        B = Xon[(k-1)*Wl+1:k*Wl, :]; B = B .- mean(B, dims=1)
        for q in 1:Wl÷2
            spec[q] += sum(abs2, vec(sum(B .* [cis(2π * q * (n - 1) / Wl) for n in 1:Wl], dims=1)))
        end
    end
    fr = collect(1:Wl÷2) ./ Wl
    function snr_at(f)
        q = argmin(abs.(fr .- f))
        loc = [i for i in max(1, q - 12):min(length(fr), q + 12) if abs(i - q) > 2 &&
               abs(fr[i] - f1) > 2 / Wl && abs(fr[i] - f2) > 2 / Wl]
        isempty(loc) && return NaN
        base = median(spec[loc]); mad = 1.4826 * median(abs.(spec[loc] .- base))
        maximum(spec[max(1, q - 1):min(length(fr), q + 1)]) > 0 ? (maximum(spec[max(1, q - 1):min(length(fr), q + 1)]) - base) / max(mad, eps()) : NaN
    end
    return (D=D, M=M, t1=t1, t1n=t1n, slope1=(m1, e1), slope2=(m2, e2), zslope=zs, T3=T3, snrP=snrP, rA=rA, dA=dA,
            snrM=snr_at(M), snrD=snr_at(D), W=W, N=N, P1=P1, P2=P2, a1=a1, a2=a2,
            spec=spec, fr=fr)
end

beat_test_core(X0, bin_st, bin_end; f1, f2, kw...) = beat_test_core(X0, bin_st, bin_end, f1, f2; kw...)

"Syntetyczny pulsar: :beat (D z modulacja amplitudy M, bez nosnej), :indep (dwa niezalezne dryfy),
:indep_scint (to samo + wspolne wolne zmiany jasnosci calego pulsara, jak scyntylacja), :modes (przelaczanie
trybow f1 <-> f2). Dwie skladowe profilu jak J1555."
function beat_synth(kind; N=1045, f1=0.0508, f2=0.0938, P2=40.0, snr=1.0, seed=1)
    rng = MersenneTwister(seed); nb = 1024; bst, bend = 440, 600
    prof(j) = exp(-0.5 * ((j - 480) / 12)^2) + 0.8exp(-0.5 * ((j - 560) / 12)^2)
    D = (f1 + f2) / 2; M = (f2 - f1) / 2
    rw(σ) = cumsum(σ .* randn(rng, N))                            # bladzenie fazy
    slow() = (s = cumsum(randn(rng, N)); s = s .- [mean(s[max(1, i - 40):min(N, i + 40)]) for i in 1:N] .+ 0; 0.6 .+ 0.4 .* tanh.(s ./ (std(s) + eps())))
    X = randn(rng, N, nb)
    if kind == :beat
        θD = rw(0.05); θM = rw(0.03); s = slow()
        for n in 1:N, j in bst:bend
            X[n, j] += 4snr * prof(j) * (1 + 0.9s[n] * cos(2π * M * (n - 1) + θM[n]) * cos(2π * (D * (n - 1) - (j - 480) / P2) + θD[n]))
        end
    elseif kind in (:indep, :indep_scint)
        θ1 = rw(0.05); θ2 = rw(0.05); s1 = slow(); s2 = slow()
        g = kind == :indep_scint ? (0.2 .+ 1.6 .* (slow() .- 0.2) ./ 0.8) : ones(N)     # wspolna jasnosc 0.2..1.8
        for n in 1:N, j in bst:bend
            X[n, j] += 4snr * g[n] * prof(j) * (1 + 0.45s1[n] * cos(2π * (f1 * (n - 1) - (j - 480) / P2) + θ1[n]) +
                                            0.45s2[n] * cos(2π * (f2 * (n - 1) - (j - 480) / (1.6P2)) + θ2[n]))
        end
    elseif kind == :modes
        mode = 1; n = 1; md = zeros(Int, N)
        while n <= N; L = 20 + round(Int, -60log(rand(rng))); md[n:min(N, n + L - 1)] .= mode; mode = 3 - mode; n += L; end
        θ = rw(0.05)
        for n in 1:N, j in bst:bend
            f = md[n] == 1 ? f1 : f2
            X[n, j] += 4snr * prof(j) * (1 + 0.9cos(2π * (f * (n - 1) - (j - 480) / P2) + θ[n]))
        end
    end
    return X, bst, bend
end

function _bt_line(name, r)
    v = r.t1.p < 0.05 ? (r.t1.r > 0 ? "razem (jak dudnienia)" : "naprzemiennie (jak przełączanie trybów)") : "bez związku (jak niezależne)"
    vn = r.t1n.p < 0.05 ? (r.t1n.r > 0 ? "razem" : "naprzemiennie") : "bez związku"
    @sprintf("%-22s T1 moc obu pików: r=%+.2f (p=%.3f) -> %s; po podzieleniu przez jasność: r=%+.2f (p=%.3f) -> %s | T3 nachylenia w składowych: %s (największe |z|=%.1f) | T4 kształt amplitudy: r=%.2f, max różnica %.0f%% | T5 S/N przy M=%.4f: %.1f, przy D=%.4f: %.1f",
             name, r.t1.r, r.t1.p, v, r.t1n.r, r.t1n.p, vn,
             join([@sprintf("skł.%d %.2f±%.2f vs %.2f±%.2f", t.c, t.s1..., t.s2...) for t in r.T3], "; ") * " °/bin", r.zslope, r.rA, 100r.dA, r.M, r.snrM, r.D, r.snrD)
end

"Tylko dane syntetyczne (dziala lokalnie, bez danych pulsarow)."
function beat_test_synth(; W=64, seeds=1:3)
    L = String[]
    for kind in (:beat, :indep, :indep_scint, :modes), s in seeds
        X, b1, b2 = beat_synth(kind; seed=s)
        push!(L, _bt_line("synt. $(kind) #$s", beat_test_core(X, b1, b2, 0.0508, 0.0938; W=W)))
        println(L[end]); flush(stdout)
    end
    L
end

# piki do testu: J1555 (nasze z V1: 0.0508 i 0.0938); J1901: F3 i F1 (F1 - F3 = F2, wiec F2 = "2M"... kontrola
# ukladu dudnien z TPA XX w innym ukladzie - interpretacja ostrozna)
BT_PSRS = Dict("J1555-3134" => (0.0508, 0.0938), "J1901-0906" => (0.1406, 0.3281))

function beat_test_run(; W=64, outdir="/home/psr/data/OUTPUT/beat_test/")
    isdefined(Main, :_sl_data) || include(joinpath(@__DIR__, "sliding_lrfs.jl"))
    mkpath(outdir)
    L = ["# Test systemu dudnień (okna $W impulsów, krok 1)", "",
         "Oczekiwania: DUDNIENIA - T1 r > 0 istotne, T3 te same nachylenia (|z| < 2), T4 ten sam kształt (r ≈ 1);",
         "NIEZALEŻNE - T1 r ≈ 0 (nieistotne); PRZEŁĄCZANIE TRYBÓW - T1 r < 0 istotne. Różne nachylenia (|z| > 3) wykluczają dudnienia.",
         "Rozstrzyga wersja T1 po podzieleniu przez jasność (wspólna scyntylacja / nulle dają pozorne 'razem' w zwykłej T1).",
         "J1901-0906 (TPA XX): para F3 i F1 - F3 to dudnienie F1 z F2, więc jego moc powinna iść razem z mocą F1 (r > 0).",
         "", "## Dane syntetyczne (znana odpowiedź)"]
    append!(L, beat_test_synth(; W=W))
    push!(L, "", "## Pulsary")
    for (psr, (f1, f2)) in BT_PSRS
        try
            D, p = Base.invokelatest(_sl_data, psr)
            r = beat_test_core(D, Int(p["bin_st"]), Int(p["bin_end"]), f1, f2; W=W)
            push!(L, _bt_line(psr, r)); println(L[end])
            # czulosc testu przy TEJ jasnosci pikow, liczbie impulsow i czestotliwosciach: syntetyczne dudnienia
            # i niezalezne dryfy z sila pikow (moc/szum) dopasowana do pulsara
            N = size(D, 1)
            cand = [(sn, beat_test_core(beat_synth(:beat; N=N, f1=f1, f2=f2, snr=sn, seed=99)...; f1=f1, f2=f2, W=W, nperm=100).snrP)
                    for sn in (0.03, 0.05, 0.08, 0.12, 0.2, 0.3, 0.5, 0.8, 1.2)]
            sn = cand[argmin([abs(log(max(c[2], 1e-6) / max(r.snrP, 1e-6))) for c in cand])][1]
            for kind in (:beat, :indep)
                rs = [beat_test_core(beat_synth(kind; N=N, f1=f1, f2=f2, snr=sn, seed=s)...; f1=f1, f2=f2, W=W, nperm=300).t1n for s in 1:5]
                push!(L, @sprintf("    czułość: syntetyczne %-6s przy sile pików jak w %s (moc/szum %.1f, snr=%.2f): T1 po podzieleniu r = %s; p = %s",
                                  kind, psr, r.snrP, sn, join([@sprintf("%+.2f", x.r) for x in rs], ", "), join([@sprintf("%.3f", x.p) for x in rs], ", ")))
                println(L[end])
            end
        catch e
            e isa InterruptException && rethrow()
            push!(L, "$psr BŁĄD: " * first(sprint(showerror, e), 200))
        end
    end
    write(outdir * "WYNIK.txt", join(L, "\n") * "\n")
    println("\nZapisane: $(outdir)WYNIK.txt")
end
