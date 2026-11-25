# Tabela cech dryfu z klasa (dryf / tylko P3) jak TPA VIII (Song+2023) + kalibracja na prawdziwym szumie.
# Rdzen obliczen: modules/drift_core.jl (funkcje dc_*, w module Plot). Kroki:
#  1. widmo pelnej rozdzielczosci, tlo = biezaca mediana z maska cech i oknem rosnacym z czestotliwoscia (szerokie
#     i wedrujace cechy nie sa juz wchlaniane przez tlo - J2053, J0304 wczesniej "brak cech"), szum MAD;
#  2. przeszukiwanie wieloskalowe; prog z impulsow przetasowanych kazdego pulsara (×1.2, min. DT_THR);
#  3. szerokosc σ_1/P3 i P3 ze srodka ciezkosci mocy cechy (TPA VIII r. 2-3), bez dopasowania;
#  4. dlugosc FFT jak TPA VIII: od 512 (min. 64) polowiona, dopoki rosnie liczba istotnych (> 3σ) komorek LRFS;
#  5. klasa w kazdej skladowej z 2DFS: srodek ciezkosci 1/P2 w prostokacie wokol cechy (po odjeciu 2DFS impulsow
#     przetasowanych), blad z kilku prostokatow + bootstrap blokow; DRYF gdy znak 1/P2 wyznaczony, inaczej TYLKO P3;
#     P_asym (TPA VIII r. 4) z bledem z tasowania;
#  6. nachylenie fazy ze sciezki sredniej w kazdej skladowej (+ polowki obserwacji); w naszych jednostkach
#     nachylenie [°/°] = 1/P2 [cykli/P], wiec oba sposoby mozna wprost porownac (niezgodny znak = flaga);
#  7. stabilnosc f w polowkach (tercjach), flagi: harmoniczna, alias 0.5, wolna (P3 > 100 albo blisko zera: nulle,
#     zmiany jasnosci), malo impulsow;
#  8. literatura: input/drift_pulsars_P3.txt (klasa dryf) i p3only_pulsars_P3.txt (tylko P3), cenzus TPA VIII;
#  9. kalibracja: przetasowane impulsy (falszywe cechy) + wstrzykniety dryf o znanym f, P2 (znak losowy), stalym
#     albo wedrujacym P3 -> kompletnosc, poprawna klasa i znak.
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/drift_table.jl"); drift_table_run()
# Wynik: /home/psr/data/OUTPUT/drift_table/{cechy.csv, WYNIK.txt}
isdefined(Main, :_sl_data) || include(joinpath(@__DIR__, "sliding_lrfs.jl"))
using Random, Statistics, Printf, FFTW
DT_DIR = F6_ROOT * "drift_table/"
DT_THR = 6.0            # najnizszy prog istotnosci cechy
DT_NPERM = 10           # tasowania do 2DFS (odjecie tla stochastycznego, blad P_asym)

"Prog istotnosci pulsara: najwieksza istotnosc w widmach jego impulsow w losowej kolejnosci ×1.2, min. DT_THR."
function _dt_threshold(X, on, off, rng; nshuf=10)
    PL = SpaTs.Plot
    N = size(X, 1)
    z = Float64[]
    for _ in 1:nshuf
        R, _, f = PL.dc_spectrum(X[shuffle(rng, 1:N), :], on, off)
        push!(z, PL.dc_features(R, f; onlymax=true))
    end
    return max(DT_THR, 1.2 * maximum(z))
end

"Lista P3 z literatury (cenzus TPA VIII): nazwa => (P3, blad, klasa)."
function _dt_literature()
    lit = Dict{String,Tuple{Float64,Float64,String}}()
    for (f, cl) in (("drift_pulsars_P3.txt", "dryf"), ("p3only_pulsars_P3.txt", "tylko P3"))
        pf = joinpath(@__DIR__, "..", "input", f); isfile(pf) || continue
        for l in eachline(pf)
            m = match(r"^(\S+)\s+([\d.]+)\((\d+)\)", l); isnothing(m) && continue
            v = m.captures[2]; e = m.captures[3]
            dec = occursin('.', v) ? length(split(v, '.')[2]) : 0
            lit[m.captures[1]] = (parse(Float64, v), parse(Float64, e) * 10.0^(-dec), cl)
        end
    end
    lit
end

_dt_rngidx(rng) = (n, k) -> [rand(rng, 1:n, n) for _ in 1:k]
_dt_pm(v, e) = isfinite(v) ? @sprintf("%.1f±%.1f", v, e) : "-"

"Klasa i nachylenia jednej cechy (f = mu, σ) we wszystkich skladowych."
function _dt_feature(X, on, off, nb, mu, σf, nfft, comps, jref, perms, rng; halves=true)
    PL = SpaTs.Plot
    tr = PL.dc_track(X, on, off, mu, nfft; hw=2σf, comps=comps, jref=jref, dlon=360 / nb)
    hs = []
    if halves
        N = size(X, 1)
        for r in (1:N÷2, N÷2+1:N)
            push!(hs, PL.dc_track(X[r, :], on, off, mu, min(nfft, prevpow(2, length(r))); hw=2σf, comps=comps, jref=jref, dlon=360 / nb))
        end
    end
    cl = [PL.dc_classify(X, on[r], mu, σf, nfft, nb, perms; rng_idx=_dt_rngidx(rng)) for r in comps]
    return tr, hs, cl
end

"Wszystko dla jednego pulsara: cechy, klasa, nachylenia, stabilnosc, flagi."
function _dt_one(psr, D, p; lit=Dict(), seed=11)
    PL = SpaTs.Plot
    rng = MersenneTwister(seed)
    N, nb = size(D)
    bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
    X, good = PL.dc_demean(D)
    on = bst:bend; off = PL._far_offpulse(nb, bst, bend)
    R, _, freq = PL.dc_spectrum(X, on, off)
    thr = _dt_threshold(X, on, off, rng)
    feats = PL.dc_features(R, freq; thr=thr)
    nfft, nhist = PL.dc_nfft(X, on, off)
    prof = vec(mean(Float64.(D[good, on]), dims=1))
    jref = argmax(prof)
    comps = PL._profile_components(prof)
    perms = [shuffle(rng, 1:N) for _ in 1:DT_NPERM]
    rows = []
    for (q, ft) in enumerate(feats)
        flags = String[]
        N < 300 && push!(flags, "mało impulsów")
        ft.mu > 0.45 && push!(flags, "blisko aliasu 0.5")
        slow = ft.mu < 4 / N + 0.01 || 1 / ft.mu > 100
        slow && push!(flags, "wolna (nulle / zmiany jasności?)")
        for o in feats
            o === ft && continue
            abs(ft.mu - 2o.mu) <= max(2(ft.sigma + o.sigma), 2 / N) && push!(flags, @sprintf("harmoniczna %.4f", o.mu))
        end
        tr, hs, cl = _dt_feature(X, on, off, nb, ft.mu, ft.sigma, nfft, comps, jref, perms, rng)
        # klasa: dryf, gdy w ktorejkolwiek skladowej znak 1/P2 jest wyznaczony
        isdrift = any(c -> !isnothing(c) && c.drift, cl)
        klasa = isdrift ? "dryf" : "tylko P3"
        dirs = join([isnothing(c) ? "?" : c.sign for c in cl], "")
        for (c, cc) in enumerate(cl)
            isnothing(cc) && continue
            s, e = tr.slopes[c], tr.errs[c]
            if cc.drift && isfinite(s) && abs(s) > 2e && sign(s) != sign(cc.inv_p2)
                push!(flags, "skł.$c: znak nachylenia fazy ≠ znak 1/P2")
            end
        end
        # stabilnosc f: polowki (tercje)
        nparts = N >= 3000 ? 3 : 2
        fs = Float64[]
        for k in 1:nparts
            r = (k - 1) * (N ÷ nparts) + 1:k * (N ÷ nparts)
            Rk, _, fk = PL.dc_spectrum(X[r, :], on, off)
            Bk, _ = PL.dc_background(Rk)
            w = [i for i in eachindex(fk) if abs(fk[i] - ft.mu) <= 3ft.sigma + 2 / length(r)]
            isempty(w) && continue
            sm = PL.dc_box(Rk .- Bk, max(1, ft.k ÷ nparts))
            push!(fs, fk[w[argmax(sm[w])]])
        end
        length(fs) >= 2 && (maximum(fs) - minimum(fs) > max(2.355ft.sigma, 2nparts / N)) &&
            push!(flags, "niestabilna (części: " * join([@sprintf("%.4f", x) for x in fs], ", ") * ")")
        slopes = join([_dt_pm(tr.slopes[c], tr.errs[c]) for c in eachindex(tr.slopes)], " | ")
        hslopes = join([join([_dt_pm(h.slopes[c], h.errs[c]) for h in hs], " / ") for c in eachindex(tr.slopes)], " | ")
        p2s = join([isnothing(c) ? "-" : _dt_pm(c.inv_p2, c.err) for c in cl], " | ")
        pas = join([isnothing(c) ? "-" : @sprintf("%.2f±%.2f", c.pasym, c.pasym_err) for c in cl], " | ")
        push!(rows, (psr=psr, q=q, f=ft.mu, P3=1 / ft.mu, sf=ft.sigma, snr=ft.snr, N=N, nfft=nfft, klasa=klasa, dirs=dirs,
                     slopes=slopes, hslopes=hslopes, p2=p2s, pasym=pas, flags=flags, slow=slow))
    end
    # literatura: najblizsza cecha (P3 w granicach bledu albo 5%) i jej klasa
    litflag = ""
    if haskey(lit, psr)
        P, e, lcl = lit[psr]
        if isempty(rows)
            litflag = @sprintf("literatura P3=%.3g(%.2g) %s: brak cech", P, e, lcl)
        else
            d = [abs(r.P3 - P) / sqrt(e^2 + (r.sf / r.f^2)^2) for r in rows]
            j = argmin(d)
            if d[j] <= 3 || abs(rows[j].P3 / P - 1) <= 0.05
                litflag = @sprintf("literatura P3=%.3g(%.2g) %s: zgodna cecha %d, klasa %s", P, e, lcl, rows[j].q,
                                   rows[j].klasa == lcl ? "zgodna" : "INNA ($(rows[j].klasa))")
            else
                litflag = @sprintf("literatura P3=%.3g(%.2g) %s: NIEZGODNA (najbliższa P3=%.3g, %.1fσ)", P, e, lcl, rows[j].P3, d[j])
            end
        end
    end
    return rows, litflag, (thr=thr, nfft=nfft, nhist=nhist, comps=comps)
end

"""
Kalibracja na prawdziwym szumie pulsara: przetasowane impulsy (falszywe cechy) + wstrzykniety dryf (f0, P2 o
losowym znaku; staly albo wedrujacy P3: f(n) = f0 + 0.015 sin(2πn/300)) -> wykrycie, poprawna klasa i znak.
"""
function _dt_calibrate(D, p; nshuf=5, amps=(0.03, 0.06, 0.12, 0.25), nf0=3, seed=1, thr=DT_THR, nfft=256)
    PL = SpaTs.Plot
    rng = MersenneTwister(seed)
    N, nb = size(D)
    bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
    on = bst:bend; off = PL._far_offpulse(nb, bst, bend)
    good = PL._live_rows(D)
    prof = vec(mean(Float64.(D[good, on]), dims=1))
    comps = PL._profile_components(prof)
    nfalse = 0
    res = Dict((a, w) => (det=Int[], ok=Int[], ferr=Float64[]) for a in amps, w in (false, true))
    for s in 1:nshuf
        Y = Float64.(D[shuffle(rng, 1:N), :])
        X, _ = PL.dc_demean(Y)
        R, _, freq = PL.dc_spectrum(X, on, off)
        nfalse += length(PL.dc_features(R, freq; thr=thr))
        for a in amps, wander in (false, true), _ in 1:nf0
            f0 = 0.05 + 0.4rand(rng); P2 = (30 + 40rand(rng)) * rand(rng, (-1, 1))
            Z = copy(Y); φ = 0.0
            for n in 1:N
                φ += 2π * (f0 + (wander ? 0.015sin(2π * n / 300) : 0.0))
                good[n] || continue
                for (jj, j) in enumerate(on)
                    Z[n, j] += a * prof[jj] * 2cos(φ - 2π * (jj - 1) / P2)
                end
            end
            Xz, _ = PL.dc_demean(Z)
            Rz, _, fz = PL.dc_spectrum(Xz, on, off)
            fs = PL.dc_features(Rz, fz; thr=thr)
            hit = findfirst(x -> abs(x.mu - f0) <= max(3 / N, 3x.sigma, wander ? 0.02 : 0.0), fs)
            r = res[(a, wander)]
            push!(r.det, isnothing(hit) ? 0 : 1)
            isnothing(hit) && continue
            push!(r.ferr, fs[hit].mu - f0)
            perms = [shuffle(rng, 1:N) for _ in 1:3]
            cl = [PL.dc_classify(Xz, on[c], fs[hit].mu, fs[hit].sigma, nfft, nb, perms; rng_idx=_dt_rngidx(rng)) for c in comps]
            want = P2 > 0 ? "+" : "-"
            push!(r.ok, any(c -> !isnothing(c) && c.drift && c.sign == want, cl) ? 1 : 0)
        end
    end
    return (nfalse=nfalse, nshuf=nshuf, res=res, amps=amps)
end

function drift_table_run(psrs=SL_PSRS; outdir=DT_DIR, calibrate=true)
    mkpath(outdir)
    lit = _dt_literature()
    L = ["# Tabela cech dryfu (TPA VIII: klasa dryf / tylko P3 z 2DFS; P3 i σ_1/P3 ze środka ciężkości mocy; nasza konwencja: + = dryf ku późniejszym długościom;",
         "# nachylenie fazy [°/°] i 1/P2 [cykli/P] to ta sama wielkość - dwa niezależne pomiary)", ""]
    csv = ["pulsar,cecha,f,P3,sigma_f,SN,N,nfft,klasa,kierunki,nachylenia,nachylenia_polowki,invP2,Pasym,flagi"]
    review = String[]
    cal = String[]
    for psr in psrs
        try
            D, p = Base.invokelatest(_sl_data, psr)
            rows, litflag, info = Base.invokelatest(_dt_one, psr, D, p; lit=lit)
            push!(L, @sprintf("%s (N=%d, próg %.1f, nfft %d, składowe %d)", psr, size(D, 1), info.thr, info.nfft, length(info.comps)) *
                     (isempty(litflag) ? "" : "  -  " * litflag))
            isempty(rows) && push!(L, "    brak istotnych cech")
            for r in rows
                push!(L, @sprintf("    F%d: P3=%.2f (f=%.4f, σ_f=%.4f) S/N=%.1f  klasa: %s [%s]  1/P2: %s  nachylenie: %s (połówki %s)  P_asym: %s%s",
                                  r.q, r.P3, r.f, r.sf, r.snr, r.klasa, r.dirs, r.p2, r.slopes, r.hslopes, r.pasym,
                                  isempty(r.flags) ? "" : "  FLAGI: " * join(r.flags, "; ")))
                push!(csv, join([psr, r.q, @sprintf("%.5f", r.f), @sprintf("%.3f", r.P3), @sprintf("%.5f", r.sf), @sprintf("%.1f", r.snr),
                                 r.N, r.nfft, r.klasa, r.dirs, "\"" * r.slopes * "\"", "\"" * r.hslopes * "\"", "\"" * r.p2 * "\"",
                                 "\"" * r.pasym * "\"", "\"" * join(r.flags, "; ") * "\""], ","))
            end
            (any(!isempty(r.flags) for r in rows) || occursin(r"NIEZGODNA|INNA|brak cech", litflag)) && push!(review, psr)
            if calibrate
                c = Base.invokelatest(_dt_calibrate, D, p; thr=info.thr, nfft=info.nfft)
                fmt(a, w) = (r = c.res[(a, w)]; @sprintf("%.2f: %d%%/%s", a, round(Int, 100 * mean(r.det)),
                                                          isempty(r.ok) ? "-" : string(round(Int, 100 * mean(r.ok))) * "%"))
                push!(cal, @sprintf("%-12s fałszywe: %d w %d | wykryte/klasa i znak poprawne - stały P3: ", psr, c.nfalse, c.nshuf) *
                           join([fmt(a, false) for a in c.amps], ", ") * " | wędrujący P3: " * join([fmt(a, true) for a in c.amps], ", "))
                println(cal[end])
            end
        catch e
            e isa InterruptException && rethrow()
            push!(L, "$psr BŁĄD: " * first(sprint(showerror, e), 300))
        end
        println(L[end]); flush(stdout)
    end
    push!(L, "", "## Do ręcznego przejrzenia (flagi albo niezgodność z literaturą): " * (isempty(review) ? "brak" : join(review, ", ")))
    calibrate && append!(L, ["", "## Kalibracja na prawdziwym szumie (impulsy w losowej kolejności; siła = amplituda modulacji / średni profil)"], cal)
    write(outdir * "WYNIK.txt", join(L, "\n") * "\n")
    write(outdir * "cechy.csv", join(csv, "\n") * "\n")
    println("\nZapisane: $(outdir)WYNIK.txt, cechy.csv")
end
