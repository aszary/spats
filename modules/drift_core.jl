# Rdzen analizy dryfu (bez PyPlot - testowalny osobno): wlaczany do modulu Plot (plot.jl) i do testow.
# Metody wg TPA VIII (Song+2023, MNRAS 520, 4562) tam, gdzie sa opisane:
#  - dlugosc FFT: potega 2 od 512 w dol (min. 64), polowiona, dopoki rosnie liczba komorek LRFS istotnych > 3σ
#    (u nich: "liczba binow LRFS z istotna moca", szczegoly w materiale online D5 - tu nasze przyblizenie);
#  - szerokosc cechy σ_1/P3: srednia i wariancja wazona moca po odjeciu tla, bez dopasowania (ich r. 2-3);
#  - 2DFS skladowych, srodek ciezkosci 1/P2 w prostokacie wokol cechy (po odjeciu sredniego widma impulsow
#    przetasowanych), blad = rozrzut miedzy kilkoma prostokatami + bootstrap blokow; DRYF gdy znak P2 jest
#    wyznaczony w granicach bledu, inaczej TYLKO P3 (ich 2.3.5 i dodatek B);
#  - P_asym (ich r. 4) z bledem z tasowania impulsow.
# Nasze: tlo widma pelnej rozdzielczosci (biezaca mediana z maska cech, okno rosnace z czestotliwoscia - szerokie
# i wedrujace cechy nie sa juz wchlaniane przez tlo), nachylenie fazy ze SCIEZKI SREDNIEJ (spojna suma faz z okien,
# wazona amplituda; rozwijanie tylko w skladowej, od jej najsilniejszego binu) zamiast sredniej z nachylen okien.
# Konwencja: + = dryf ku pozniejszym dlugosciom (faza rosnie z dlugoscia, 1/P2 > 0).
using FFTW, Statistics

dc_live(X) = [any(!iszero, @view X[i, :]) for i in 1:size(X, 1)]

"Stos z odjeta srednia w kazdym binie; wyciete (zerowe) impulsy zostaja zerami."
function dc_demean(X)
    good = dc_live(X); Y = Float64.(X)
    any(good) || return Y, good
    Y .-= mean(Y[good, :], dims=1)
    Y[.!good, :] .= 0.0
    return Y, good
end

dc_box(x, k) = k <= 1 ? copy(x) : (c = vcat(0.0, cumsum(x)); h = k ÷ 2;
                                   [(c[min(length(x), i + h) + 1] - c[max(1, i - h)]) / (min(length(x), i + h) - max(1, i - h) + 1) for i in eachindex(x)])
dc_mad(u) = 1.4826 * median(abs.(u .- median(u)))

"Widmo pelnej rozdzielczosci: R(f) = moc w oknie impulsu / sredni szum (moc poza impulsem), f = k/N, k = 1..N/2."
function dc_spectrum(X, on, off)
    N = size(X, 1)
    Fo = rfft(X[:, on] .- mean(X[:, on], dims=1), 1); Ff = rfft(X[:, off] .- mean(X[:, off], dims=1), 1)
    nf = N ÷ 2
    Pon = vec(sum(abs2.(Fo[2:nf+1, :]), dims=2)); Poff = vec(sum(abs2.(Ff[2:nf+1, :]), dims=2)) .* (length(on) / length(off))
    s = mean(Poff)
    return Pon ./ s, Poff ./ s, collect(1:nf) ./ N
end

"""
Tlo widma: biezaca mediana z oknem rosnacym z czestotliwoscia (±50% indeksu, min. 5 binow, max. ±0.05 cyklu/P -
czerwony szum przy zerze zmienia sie szybko, dalej tlo jest plaskie), z maska cech (iteracyjnie: biny > 3σ po
wygladzeniu, poszerzone o 4) - szerokie i wedrujace cechy (J2053, J0304) nie podnosza juz tla.
"""
function dc_background(R; iters=4, zmask=3.0, frac=0.5)
    nf = length(R); maxw = max(10, nf ÷ 10)
    hw = [clamp(round(Int, frac * i), 5, maxw) for i in 1:nf]
    mask = falses(nf); B = zeros(nf)
    for it in 1:iters
        for i in 1:nf
            r = max(1, i - hw[i]):min(nf, i + hw[i])
            v = [R[q] for q in r if !mask[q]]
            B[i] = length(v) >= 3 ? median(v) : median(@view R[r])
        end
        it == iters && break
        d = dc_box(R .- B, 4)
        u = d[.!mask]; length(u) < 10 && break
        σ = dc_mad(u)
        m2 = falses(nf)
        for i in findall(d .> zmask * σ); m2[max(1, i - 4):min(nf, i + 4)] .= true; end
        m2 == mask && break
        mask = m2
    end
    return B, mask
end

"Szerokosc cechy jak TPA VIII (r. 2-3): moc po odjeciu tla w obszarze cechy (tlo podniesione, gdy ujemne)."
function dc_width(D, freq, i, k, N)
    sm = dc_box(D, k)
    cap = max(5, 3k)
    # obszar cechy: do zejscia do tla ALBO do dolka (gdy wygladzone widmo znow rosnie - sasiednia cecha;
    # J1555: dwie cechy 0.057 i 0.098 zlewaly sie w jedna o P3 = 12.5)
    # (sasiednia cecha = wzrost o > 30% wysokosci szczytu ponad dotychczasowe minimum; drobne wahania szumu nie tna)
    function edge(dir)
        e = i; mn = sm[i]
        while true
            n = e + dir
            (1 <= n <= length(D) && abs(n - i) <= cap && sm[n] > 0) || break
            sm[n] > mn + 0.3 * (sm[i] - mn) && sm[n] > sm[e] && break
            e = n; mn = min(mn, sm[n])
        end
        e
    end
    lo = edge(-1); hi = edge(1)
    P = D[lo:hi]; mn = minimum(P); mn < 0 && (P = P .- mn)
    s = sum(P)
    s > 0 || return (mu=freq[i], sigma=1 / N, lo=lo, hi=hi)
    μ = sum(P .* freq[lo:hi]) / s
    σ = sqrt(sum(P .* (freq[lo:hi] .- μ) .^ 2) / s)
    return (mu=μ, sigma=max(σ, 1 / N), lo=lo, hi=hi)
end

"""
Wieloskalowe wyszukiwanie cech (wygladzanie k = 1, 2, 4, ... binow; istotnosc = sygnal / szum MAD przy tej skali,
z binow poza maska cech). Zwraca liste (i, f = szczyt, mu, sigma = σ_1/P3, k, snr), najwyzej `maxfeat`.
"""
function dc_features(R, freq; thr=6.0, fmin_bins=3, maxk=nothing, onlymax=false, maxfeat=4)
    nf = length(R); N = 2nf
    B, mask = dc_background(R)
    D = R .- B
    maxk = something(maxk, max(1, round(Int, 0.05 * N)))
    ks = [k for k in (1, 2, 4, 8, 16, 32, 64, 128) if k <= maxk]
    best = fill(-Inf, nf); bk = ones(Int, nf)
    for k in ks
        sm = dc_box(D, k)
        u = sm[.!mask]; length(u) < 10 && (u = sm)
        σ = dc_mad(u)
        σ > 0 || continue
        z = sm ./ σ
        upd = z .> best
        best[upd] .= z[upd]; bk[upd] .= k
    end
    onlymax && return maximum(best[fmin_bins:end])
    loc = [i for i in fmin_bins:nf-1 if best[i] >= best[i-1] && best[i] >= best[i+1] && best[i] >= thr]
    sort!(loc, by=i -> -best[i])
    feats = Int[]
    for i in loc
        sep = max(3, 2bk[i])
        smi = dc_box(D, bk[i])
        ok = all(abs(i - q) >= sep && minimum(smi[min(i, q):max(i, q)]) < 0.7 * smi[i] for q in feats)
        ok && push!(feats, i)
        length(feats) >= maxfeat && break
    end
    out = []; regs = UnitRange{Int}[]
    for i in feats                                   # od najsilniejszej
        w = dc_width(D, freq, i, bk[i], N)
        # ta sama (np. wedrujaca) cecha rozbita na kilka szczytow: srodek w obszarze silniejszej -> pomijamy
        im = clamp(round(Int, w.mu * N), 1, nf)
        any(r -> i in r || im in r, regs) && continue
        push!(regs, w.lo:w.hi)
        push!(out, (i=i, f=freq[i], mu=w.mu, sigma=w.sigma, k=bk[i], snr=best[i]))
    end
    return out
end

"Srednie LRFS (moc) z niezachodzacych blokow po nfft impulsow: (on, off), wiersze f = 1..nfft/2 (bez DC)."
function dc_lrfs_power(X, on, off, nfft)
    nbk = size(X, 1) ÷ nfft
    h = nfft ÷ 2
    Pon = zeros(h, length(on)); Poff = zeros(h, length(off))
    for b in 1:nbk
        r = (b - 1) * nfft + 1:b * nfft
        Fo = rfft(X[r, on] .- mean(X[r, on], dims=1), 1); Ff = rfft(X[r, off] .- mean(X[r, off], dims=1), 1)
        Pon .+= abs2.(Fo[2:h+1, :]); Poff .+= abs2.(Ff[2:h+1, :])
    end
    return Pon ./ max(nbk, 1), Poff ./ max(nbk, 1)
end

"""
Dlugosc FFT (TPA VIII, dodatek B: najwyzej 512 albo najwieksza potega 2 miesci sie w obserwacji, najmniej 64):
od najdluzszej, polowiona dopoki rosnie liczba komorek LRFS istotnych (> srednia + nsig·σ szumu poza impulsem).
Zwraca (nfft, historia [(nfft, liczba)]).
"""
function dc_nfft(X, on, off; nsig=3.0, nmax=512, nmin=64)
    N = size(X, 1)
    n = min(nmax, prevpow(2, max(N, 2)))
    n <= nmin && return n, [(n, -1)]
    cnt(n) = (P = dc_lrfs_power(X, on, off, n); m = mean(P[2]); s = std(P[2]); count(>(m + nsig * s), P[1]))
    c = cnt(n); hist = [(n, c)]
    while n ÷ 2 >= nmin
        c2 = cnt(n ÷ 2); push!(hist, (n ÷ 2, c2))
        c2 > c || break
        n ÷= 2; c = c2
    end
    return n, hist
end

"""
Sciezka srednia fazy przy czestotliwosci f: okna nfft impulsow co `step`; w kazdym oknie jedna czestotliwosc
(najsilniejsza w f ± hw), faza obrocona tak, by w binie `jref` (szczyt profilu) byla 0 (okna, w ktorych jref jest
nieistotny, pomijane); sumowanie zespolone po oknach (waga = amplituda). Biny istotne: |Z| > nsig·σ (σ z tej samej
sumy dla binow poza impulsem). W kazdej skladowej rozwijanie od jej najsilniejszego binu, prosta wazona (1/blad²),
blad powiekszony o sqrt(χ²/(n-2)), gdy rozrzut wiekszy od szumu. Zwraca nachylenia [°/°] i bledy na skladowa.
"""
function dc_track(X, on, off, f, nfft; hw=0.0, step=max(1, nfft ÷ 4), comps=nothing, jref=1, nsig=3.0, dlon=1.0)
    N = size(X, 1)
    nc = isnothing(comps) ? 1 : length(comps)
    empty = (phase=fill(NaN, length(on)), amp=zeros(length(on)), sig=falses(length(on)), slopes=fill(NaN, nc),
             errs=fill(NaN, nc), nwin=0, npts=zeros(Int, nc))
    N >= nfft || return empty
    h = nfft ÷ 2; fr = collect(1:h) ./ nfft
    band = [k for k in 1:h if abs(fr[k] - f) <= max(hw, 1 / nfft) + 0.5 / nfft]
    isempty(band) && (band = [clamp(round(Int, f * nfft), 1, h)])
    Z = zeros(ComplexF64, length(on)); Zo = zeros(ComplexF64, length(off))
    nw = 0
    for s in 1:step:N-nfft+1
        r = s:s+nfft-1
        Fo = rfft(X[r, on] .- mean(X[r, on], dims=1), 1); Ff = rfft(X[r, off] .- mean(X[r, off], dims=1), 1)
        k = band[argmax([sum(abs2, @view Fo[k + 1, :]) for k in band])]
        c = conj.(Fo[k + 1, :]); co = conj.(Ff[k + 1, :])          # conj: jadro e^(+2πi f n)
        rmsw = std(vcat(real.(co), imag.(co)))
        abs(c[jref]) > nsig * rmsw || continue
        rot = conj(c[jref]) / abs(c[jref])
        Z .+= c .* rot; Zo .+= co .* rot; nw += 1
    end
    nw >= 1 || return empty
    σ = std(vcat(real.(Zo), imag.(Zo)))
    amp = abs.(Z); sig = amp .> nsig * σ
    ph = rad2deg.(angle.(Z))
    comps = something(comps, [1:length(on)])
    slopes = fill(NaN, length(comps)); errs = fill(NaN, length(comps)); npts = zeros(Int, length(comps))
    U = fill(NaN, length(on))
    for (c, r) in enumerate(comps)
        js = [j for j in r if sig[j]]
        length(js) >= 3 || continue
        j0 = js[argmax(amp[js])]
        U[j0] = ph[j0]
        for dir in (1, -1)
            prev = ph[j0]; j = j0 + dir
            while first(r) <= j <= last(r)
                if sig[j]
                    U[j] = ph[j] - 360 * round((ph[j] - prev) / 360); prev = U[j]
                end
                j += dir
            end
        end
        x = (js .- jref) .* dlon; y = U[js]
        e = rad2deg.(σ ./ amp[js])                                  # blad fazy w binie [°]
        w = 1 ./ e .^ 2
        xm = sum(w .* x) / sum(w); ym = sum(w .* y) / sum(w)
        Sxx = sum(w .* (x .- xm) .^ 2)
        Sxx > 0 || continue
        b = sum(w .* (x .- xm) .* (y .- ym)) / Sxx
        χ2 = sum(w .* (y .- ym .- b .* (x .- xm)) .^ 2)
        slopes[c] = b
        errs[c] = sqrt(1 / Sxx) * sqrt(max(1.0, χ2 / max(length(js) - 2, 1)))
        npts[c] = length(js)
    end
    return (phase=U, amp=amp ./ σ, sig=sig, slopes=slopes, errs=errs, nwin=nw, npts=npts)
end

"Ciagly kawalek `m` binow poza impulsem, najdalej od okna impulsu (do 2DFS szumu)."
function dc_offblock(nb, bin_st, bin_end, m)
    c = (bin_st + bin_end) / 2
    best = 0; bd = -1.0
    for s in 1:nb
        r = [mod1(s + q, nb) for q in 0:m-1]
        any(j -> bin_st - 20 <= j <= bin_end + 20, r) && continue
        d = minimum(min(abs(j - c), nb - abs(j - c)) for j in r)
        d > bd && (bd = d; best = s)
    end
    best == 0 && return Int[]
    return [mod1(best + q, nb) for q in 0:m-1]
end

"""
2DFS skladowej (Edwards & Stappers 2002; TPA VIII 2.3.3): bloki po nfft impulsow, FFT wzdluz impulsow i dlugosci
(skladowa dopelniona zerami do M). Zwraca moce blokow P[f3, q, blok], os 1/P3 (0..0.5) i 1/P2 [cykli/P] odbita tak,
by + = dryf ku pozniejszym dlugosciom.
"""
function dc_2dfs(X, bins, nfft, nb)
    nbk = size(X, 1) ÷ nfft
    m0 = length(bins); M = nextpow(2, max(2m0, 16))
    P = zeros(nfft ÷ 2 + 1, M, max(nbk, 0))
    A = zeros(nfft, M)
    for b in 1:nbk
        r = (b - 1) * nfft + 1:b * nfft
        A .= 0.0; A[:, 1:m0] .= X[r, bins] .- mean(X[r, bins], dims=1)
        P[:, :, b] .= abs2.(fft(rfft(A, 1), 2))
    end
    q = -fftfreq(M) .* nb
    f3 = collect(0:nfft÷2) ./ nfft
    return P, f3, q
end

"P_asym (TPA VIII r. 4) ze znakiem: + = przewaga dryfu ku pozniejszym dlugosciom."
function dc_pasym(S, q)
    pos = findall(>(0), q); num = 0.0
    for j in pos
        jm = findfirst(x -> isapprox(x, -q[j]; atol=1e-9), q)
        isnothing(jm) && continue
        num += sum(@view S[:, j]) - sum(@view S[:, jm])
    end
    den = sum(S)
    return den > 0 ? num / den : NaN
end

function _dc_centroid(S, f3, q, f, h, qmax)
    rows = findall(x -> abs(x - f) <= h, f3)
    cols = findall(x -> abs(x) <= qmax, q)
    (isempty(rows) || isempty(cols)) && return NaN
    w = vec(sum(max.(S[rows, cols], 0.0), dims=1))
    s = sum(w)
    return s > 0 ? sum(w .* q[cols]) / s : NaN
end

"""
Klasyfikacja cechy (f, σ_f) w skladowej `bins` jak TPA VIII: srodek ciezkosci 1/P2 w prostokacie f ± h w 2DFS
po odjeciu sredniego 2DFS impulsow przetasowanych (`perms` - permutacje impulsow; usuwa pionowy pas zmiennosci
stochastycznej). Blad = najwieksza roznica miedzy prostokatami (h × 0.5/1/1.5, pelny / polowa zakresu 1/P2) +
bootstrap blokow. DRYF, gdy wszystkie srodki maja ten sam znak i |1/P2| > blad. P_asym z bledem z tasowania.
"""
function dc_classify(X, bins, f, σf, nfft, nb, perms; nboot=100, rng_idx=nothing)
    P, f3, q = dc_2dfs(X, bins, nfft, nb)
    nbk = size(P, 3)
    nbk >= 1 || return nothing
    Sh = zeros(size(P, 1), size(P, 2)); Sh2 = zeros(size(Sh)); pas_sh = Float64[]
    for p in perms
        Ps, _, _ = dc_2dfs(X[p, :], bins, nfft, nb)
        Ss = dropdims(mean(Ps, dims=3), dims=3)
        Sh .+= Ss; Sh2 .+= Ss .^ 2
        push!(pas_sh, dc_pasym(Ss, q))
    end
    np = max(length(perms), 1)
    Sh ./= np
    σcell = sqrt.(max.(Sh2 ./ np .- Sh .^ 2, 0.0))       # rms komorki w widmach przetasowanych (jak TPA VIII 2.3.4)
    S = dropdims(mean(P, dims=3), dims=3)
    pas = dc_pasym(S, q)
    pas_err = length(pas_sh) >= 2 ? std(pas_sh) : NaN
    Sd = S .- Sh
    h0 = max(2σf, 1.5 / nfft); qm = maximum(abs.(q))
    c0 = _dc_centroid(Sd, f3, q, f, h0, qm)
    cs = [_dc_centroid(Sd, f3, q, f, hh, qq) for hh in (max(0.5h0, 1 / nfft), h0, 1.5h0) for qq in (qm, qm / 2)]
    cs = filter(isfinite, cs)
    syst = isempty(cs) || !isfinite(c0) ? NaN : maximum(abs.(cs .- c0))
    bs = Float64[]
    if nbk >= 2
        idx = isnothing(rng_idx) ? [rand(1:nbk, nbk) for _ in 1:nboot] : rng_idx(nbk, nboot)
        for id in idx
            Sb = dropdims(mean(P[:, :, id], dims=3), dims=3) .- Sh
            c = _dc_centroid(Sb, f3, q, f, h0, qm)
            isfinite(c) && push!(bs, c)
        end
    end
    berr = length(bs) >= 5 ? std(bs) : 0.0
    # szum srodka ciezkosci z rms komorek widm przetasowanych (zmiennosc stochastyczna + szum; przy 2 blokach
    # bootstrap blokow prawie nic nie mowi - szumowa cecha dostawala klase dryf)
    nerr = 0.0
    if isfinite(c0)
        rows = findall(x -> abs(x - f) <= h0, f3)
        w = vec(sum(max.(Sd[rows, :], 0.0), dims=1)); sw = sum(w)
        σcol2 = vec(sum(σcell[rows, :] .^ 2, dims=1))
        sw > 0 && (nerr = sqrt(sum((q .- c0) .^ 2 .* σcol2)) / sw)
    end
    err = (isfinite(syst) ? syst : 0.0) + sqrt(berr^2 + nerr^2)
    drift = isfinite(c0) && !isempty(cs) && all(sign.(cs) .== sign(c0)) && abs(c0) > err
    return (inv_p2=c0, err=err, P2=isfinite(c0) && c0 != 0 ? 360 / c0 : NaN, drift=drift,
            sign=drift ? (c0 > 0 ? "+" : "-") : "0", pasym=pas, pasym_err=pas_err, nblocks=nbk)
end
