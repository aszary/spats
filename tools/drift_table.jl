# Tabela cech dryfu z flagami + kalibracja na prawdziwym szumie (pomysly, ktore sie oplacaja):
#  1. widmo pelnej rozdzielczosci (jedno FFT po wszystkich impulsach, suma mocy po oknie impulsu), tlo czerwone
#     (nulle, scyntylacja) usuniete biezaca mediana, szum odpornie z samego widma (MAD; nie z chi^2 - po sumowaniu
#     skorelowanych binow ten rozklad nie obowiazuje);
#  2. przeszukiwanie wieloskalowe: wygladzanie oknem k = 1, 2, 4, ... binow, istotnosc = sygnal / szum przy tej
#     samej skali; skala z najwieksza istotnoscia daje FWHM; cechy = osobne szczyty (dolek miedzy nimi);
#  3. dlugosc okna dobrana do piku N = clamp(4/FWHM, 32, Nimp): nachylenie fazy w skladowych i kierunek dryfu;
#  4. stabilnosc: f w polowkach (tercjach dla >= 3000 impulsow) obserwacji;
#  5. flagi: harmoniczna, blisko aliasu (0.5) / zera, malo impulsow, niestabilna, niezgodna z literatura;
#  6. literatura: input/drift_pulsars_P3.txt i p3only_pulsars_P3.txt (cenzus TPA, Song+2023);
#  7. kalibracja: impulsy w losowej kolejnosci (prawdziwy szum i rozklad jasnosci, bez okresowosci) -> falszywe
#     wykrycia; + wstrzykniety dryf o znanym f i sile -> kompletnosc i dokladnosc f.
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/drift_table.jl"); drift_table_run()
# Wynik: /home/psr/data/OUTPUT/drift_table/{cechy.csv, WYNIK.txt}
isdefined(Main, :_sl_data) || include(joinpath(@__DIR__, "sliding_lrfs.jl"))
using Random, Statistics, Printf, FFTW
DT_DIR = F6_ROOT * "drift_table/"
DT_THR = 6.0            # prog istotnosci cechy (skalibrowany: falszywe wykrycia na przetasowanych danych)

"Widmo pelnej rozdzielczosci: R(f) = moc w oknie impulsu / sredni szum (moc poza impulsem), f = k/N, k = 1..N/2."
function _dt_spectrum(X, on, off)
    N = size(X, 1)
    Fo = rfft(X[:, on] .- mean(X[:, on], dims=1), 1); Ff = rfft(X[:, off] .- mean(X[:, off], dims=1), 1)
    nf = N ÷ 2
    Pon = vec(sum(abs2.(Fo[2:nf+1, :]), dims=2)); Poff = vec(sum(abs2.(Ff[2:nf+1, :]), dims=2)) .* (length(on) / length(off))
    s = mean(Poff)
    return Pon ./ s, Poff ./ s, collect(1:nf) ./ N
end

_dt_box(x, k) = k <= 1 ? copy(x) : (c = vcat(0.0, cumsum(x)); h = k ÷ 2;
                                    [(c[min(length(x), i + h) + 1] - c[max(1, i - h)]) / (min(length(x), i + h) - max(1, i - h) + 1) for i in eachindex(x)])
_dt_runmed(x, w) = [median(@view x[max(1, i - w):min(length(x), i + w)]) for i in eachindex(x)]

"Najwieksza istotnosc w widmie (do kalibracji progu na przetasowanych impulsach)."
_dt_maxz(R, Roff, freq; fmin_bins=3) = _dt_features(R, Roff, freq; thr=-Inf, fmin_bins=fmin_bins, onlymax=true)

"""
Wieloskalowe wyszukiwanie cech. Zwraca liste (f, k_best, snr, fwhm) - osobne szczyty: lokalne maksima istotnosci
>= thr, oddzielone dolkiem (minimum wygladzonego widma miedzy nimi < 0.7 slabszego) i >= 3 biny skali od siebie.
"""
function _dt_features(R, Roff, freq; thr=DT_THR, fmin_bins=3, maxk=nothing, onlymax=false)
    nf = length(R); N = 2nf
    maxk = something(maxk, max(1, round(Int, 0.05 * N)))
    B = _dt_runmed(R, max(5, nf ÷ 20))                        # tlo czerwone + bialy szum
    ks = [k for k in (1, 2, 4, 8, 16, 32, 64, 128) if k <= maxk]
    best = fill(-Inf, nf); bk = ones(Int, nf); sm_best = zeros(nf)
    for k in ks
        sm = _dt_box(R .- B, k)
        # szum przy tej skali z samego widma w oknie impulsu, odpornie (MAD; pojedyncze piki nie wplywaja):
        # szum spoza impulsu nie zawiera zmiennosci samej emisji od impulsu do impulsu i zanizal szum
        # (falszywe cechy na przetasowanych impulsach)
        σ = 1.4826 * median(abs.(sm .- median(sm)))
        σ > 0 || continue
        z = sm ./ σ
        upd = z .> best
        best[upd] .= z[upd]; bk[upd] .= k; sm_best[upd] .= sm[upd]
    end
    onlymax && return maximum(best[fmin_bins:end])
    loc = [i for i in fmin_bins:nf-1 if best[i] >= best[i-1] && best[i] >= best[i+1] && best[i] >= thr]
    sort!(loc, by=i -> -best[i])
    feats = Int[]
    for i in loc
        sep = max(3, 2bk[i])
        ok = all(abs(i - q) >= sep && minimum(_dt_box(R .- B, bk[i])[min(i, q):max(i, q)]) < 0.7 * _dt_box(R .- B, bk[i])[i]
                 for q in feats)
        ok && push!(feats, i)
        length(feats) >= 4 && break
    end
    out = []
    for i in feats
        k = bk[i]; sm = _dt_box(R .- B, k)
        h = sm[i] / 2
        lo = i; while lo > 1 && sm[lo-1] > h; lo -= 1; end
        hi = i; while hi < nf && sm[hi+1] > h; hi += 1; end
        wmeas = (hi - lo + 1) / N
        fwhm = sqrt(max(wmeas^2 - ((k - 1) / N)^2, (1 / N)^2))
        push!(out, (i=i, f=freq[i], k=k, snr=best[i], fwhm=fwhm))
    end
    return out
end

"""
Prog istotnosci dla pulsara: najwieksza istotnosc w widmach jego impulsow w losowej kolejnosci (prawdziwy szum i
rozklad jasnosci, bez okresowosci), z zapasem 20%, nie mniej niz DT_THR. Moc w widmie ma dlugi (wykladniczy) ogon,
wiec staly prog dawal falszywe cechy.
"""
function _dt_threshold(D, on, off; nshuf=10, seed=11)
    PL = SpaTs.Plot
    rng = MersenneTwister(seed)
    N = size(D, 1)
    z = Float64[]
    for _ in 1:nshuf
        Y = D[shuffle(rng, 1:N), :]
        X, _ = PL._demeaned(Y)
        R, Ro, f = _dt_spectrum(X, on, off)
        push!(z, _dt_maxz(R, Ro, f))
    end
    return max(DT_THR, 1.2 * maximum(z))
end

"Lista P3 z literatury (cenzus TPA): nazwa => (P3, blad)."
function _dt_literature()
    lit = Dict{String,Tuple{Float64,Float64}}()
    for f in ("drift_pulsars_P3.txt", "p3only_pulsars_P3.txt")
        pf = joinpath(@__DIR__, "..", "input", f); isfile(pf) || continue
        for l in eachline(pf)
            m = match(r"^(\S+)\s+([\d.]+)\((\d+)\)", l); isnothing(m) && continue
            v = m.captures[2]; e = m.captures[3]
            dec = occursin('.', v) ? length(split(v, '.')[2]) : 0
            lit[m.captures[1]] = (parse(Float64, v), parse(Float64, e) * 10.0^(-dec))
        end
    end
    lit
end

"Wszystko dla jednego pulsara: cechy, nachylenia, stabilnosc, flagi."
function _dt_one(psr, D, p; lit=Dict())
    PL = SpaTs.Plot
    N, nb = size(D)
    bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
    X, _ = PL._demeaned(D)
    on = bst:bend; off = PL._far_offpulse(nb, bst, bend)
    R, Roff, freq = _dt_spectrum(X, on, off)
    thr = _dt_threshold(D, on, off)
    feats = _dt_features(R, Roff, freq; thr=thr)
    rows = []
    for (q, ft) in enumerate(feats)
        flags = String[]
        N < 300 && push!(flags, "mało impulsów")
        ft.f > 0.45 && push!(flags, "blisko aliasu 0.5")
        ft.f < 4 / N + 0.01 && push!(flags, "blisko zera")
        for o in feats
            o === ft && continue
            abs(ft.f - 2o.f) <= max(ft.fwhm + o.fwhm, 2 / N) && push!(flags, @sprintf("harmoniczna %.4f", o.f))
        end
        # okno dobrane do piku: nachylenie fazy w skladowych (nasza konwencja: + = ku pozniejszym dlugosciom)
        Wn = clamp(round(Int, 4 / ft.fwhm), 32, N)
        cs = String[]; sgn = String[]
        if Wn < N
            t = PL.sliding_lrfs_tracks(D, bst, bend; window=Wn, step=max(1, Wn ÷ 4), p3=1 / ft.f)
            neff = max(1, (N - Wn + 1) ÷ Wn)
            for c in axes(t.cslope, 2)
                v = filter(isfinite, t.cslope[:, c])
                length(v) < 3 && (push!(cs, "-"); continue)
                m, e = mean(v), std(v) / sqrt(neff)
                push!(cs, @sprintf("%.1f±%.1f", m, e))
                push!(sgn, abs(m) < 2e ? "0" : (m > 0 ? "+" : "-"))
            end
        end
        # stabilnosc: polowki (tercje)
        nparts = N >= 3000 ? 3 : 2
        fs = Float64[]
        for k in 1:nparts
            r = (k - 1) * (N ÷ nparts) + 1:k * (N ÷ nparts)
            Rk, _, fk = _dt_spectrum(X[r, :], on, off)
            Bk = _dt_runmed(Rk, max(5, length(Rk) ÷ 20))
            w = [i for i in eachindex(fk) if abs(fk[i] - ft.f) <= 3ft.fwhm + 2 / length(r)]
            isempty(w) && continue
            sm = _dt_box(Rk .- Bk, max(1, ft.k ÷ nparts))
            push!(fs, fk[w[argmax(sm[w])]])
        end
        length(fs) >= 2 && (maximum(fs) - minimum(fs) > max(ft.fwhm, 2nparts / N)) && push!(flags, "niestabilna (części: " * join([@sprintf("%.4f", x) for x in fs], ", ") * ")")
        push!(rows, (psr=psr, q=q, f=ft.f, P3=1 / ft.f, fwhm=ft.fwhm, sf=ft.fwhm / 2.355, snr=ft.snr, N=N, W=Wn,
                     slopes=join(cs, " | "), sign=join(sgn, ""), flags=flags))
    end
    # literatura: czy jakas cecha zgadza sie z P3 z cenzusu
    litflag = ""
    if haskey(lit, psr)
        P, e = lit[psr]
        if isempty(rows)
            litflag = @sprintf("literatura P3=%.3g(%.2g): brak cech", P, e)
        else
            d = [abs(r.P3 - P) / sqrt(e^2 + (r.sf / r.f^2)^2) for r in rows]
            j = argmin(d)
            litflag = d[j] <= 3 || abs(rows[j].P3 / P - 1) <= 0.05 ? @sprintf("literatura P3=%.3g(%.2g): zgodna cecha %d", P, e, rows[j].q) :
                                  @sprintf("literatura P3=%.3g(%.2g): NIEZGODNA (najbliższa P3=%.3g, %.1fσ)", P, e, rows[j].P3, d[j])
        end
    end
    return rows, litflag, (R=R, freq=freq, thr=thr)
end

"Kalibracja na prawdziwym szumie pulsara: przetasowane impulsy (falszywe wykrycia) + wstrzykniety dryf."
function _dt_calibrate(D, p; nshuf=5, amps=(0.03, 0.06, 0.12, 0.25), nf0=3, seed=1, thr=DT_THR)
    PL = SpaTs.Plot
    rng = MersenneTwister(seed)
    N, nb = size(D)
    bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
    on = bst:bend; off = PL._far_offpulse(nb, bst, bend)
    good = PL._live_rows(D)
    prof = vec(mean(Float64.(D[good, on]), dims=1))
    nfalse = 0
    det = Dict(a => Int[] for a in amps); ferr = Dict(a => Float64[] for a in amps); snrs = Dict(a => Float64[] for a in amps)
    for s in 1:nshuf
        Y = Float64.(D[shuffle(rng, 1:N), :])
        X, _ = PL._demeaned(Y)
        R, Roff, freq = _dt_spectrum(X, on, off)
        nfalse += length(_dt_features(R, Roff, freq; thr=thr))
        for a in amps, _ in 1:nf0
            f0 = 0.05 + 0.4rand(rng); P2 = 30 + 40rand(rng)
            Z = copy(Y)
            for n in 1:N, (jj, j) in enumerate(on)
                Z[n, j] += a * prof[jj] * 2cos(2π * (f0 * (n - 1) - (jj - 1) / P2))
            end
            Xz, _ = PL._demeaned(Z)
            Rz, Rzo, fz = _dt_spectrum(Xz, on, off)
            fs = _dt_features(Rz, Rzo, fz; thr=thr)
            hit = findfirst(x -> abs(x.f - f0) <= max(3 / N, 2x.fwhm), fs)
            push!(det[a], isnothing(hit) ? 0 : 1)
            isnothing(hit) || (push!(ferr[a], fs[hit].f - f0); push!(snrs[a], fs[hit].snr))
        end
    end
    return (nfalse=nfalse, nshuf=nshuf, det=det, ferr=ferr, snrs=snrs, amps=amps)
end

function drift_table_run(psrs=SL_PSRS; outdir=DT_DIR, calibrate=true)
    mkpath(outdir)
    lit = _dt_literature()
    L = ["# Tabela cech dryfu (próg istotności z przetasowanych impulsów każdego pulsara, min. $(DT_THR); nasza konwencja: + = dryf ku późniejszym długościom)", ""]
    csv = ["pulsar,cecha,f,P3,FWHM,sigma_f,SN,N,okno,nachylenia_skladowych,kierunek,flagi"]
    review = String[]
    cal = String[]
    for psr in psrs
        try
            D, p = Base.invokelatest(_sl_data, psr)
            rows, litflag, info = Base.invokelatest(_dt_one, psr, D, p; lit=lit)
            push!(L, @sprintf("%s (N=%d, próg %.1f)", psr, size(D, 1), info.thr) * (isempty(litflag) ? "" : "  -  " * litflag))
            isempty(rows) && push!(L, "    brak istotnych cech")
            for r in rows
                push!(L, @sprintf("    F%d: f=%.4f (P3=%.2f) FWHM=%.4f S/N=%.1f okno=%d nachylenia [%s] kierunek %s%s",
                                  r.q, r.f, r.P3, r.fwhm, r.snr, r.W, r.slopes, isempty(r.sign) ? "?" : r.sign,
                                  isempty(r.flags) ? "" : "  FLAGI: " * join(r.flags, "; ")))
                push!(csv, join([psr, r.q, @sprintf("%.5f", r.f), @sprintf("%.3f", r.P3), @sprintf("%.5f", r.fwhm),
                                 @sprintf("%.5f", r.sf), @sprintf("%.1f", r.snr), r.N, r.W, "\"" * r.slopes * "\"", r.sign,
                                 "\"" * join(r.flags, "; ") * "\""], ","))
            end
            (any(!isempty(r.flags) for r in rows) || occursin("NIEZGODNA", litflag) || occursin("brak cech", litflag)) && push!(review, psr)
            if calibrate
                c = Base.invokelatest(_dt_calibrate, D, p; thr=info.thr)
                push!(cal, @sprintf("%-12s fałszywe cechy na przetasowanych danych: %d w %d próbach | kompletność (siła wstrzyknięcia): ", psr, c.nfalse, c.nshuf) *
                           join([@sprintf("%.2f: %d%% (S/N %s, błąd f %s)", a, round(Int, 100 * mean(c.det[a])),
                                          isempty(c.snrs[a]) ? "-" : @sprintf("%.0f", median(c.snrs[a])),
                                          isempty(c.ferr[a]) ? "-" : @sprintf("%.4f", median(abs.(c.ferr[a])))) for a in c.amps], ", "))
                println(cal[end])
            end
        catch e
            e isa InterruptException && rethrow()
            push!(L, "$psr BŁĄD: " * first(sprint(showerror, e), 200))
        end
        println(L[end]); flush(stdout)
    end
    push!(L, "", "## Do ręcznego przejrzenia (flagi albo niezgodność z literaturą): " * (isempty(review) ? "brak" : join(review, ", ")))
    calibrate && append!(L, ["", "## Kalibracja na prawdziwym szumie (impulsy w losowej kolejności; siła = amplituda modulacji / średni profil)"], cal)
    write(outdir * "WYNIK.txt", join(L, "\n") * "\n")
    write(outdir * "cechy.csv", join(csv, "\n") * "\n")
    println("\nZapisane: $(outdir)WYNIK.txt, cechy.csv")
end
