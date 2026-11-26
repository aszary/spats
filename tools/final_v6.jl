# OSTATECZNE PRZEJSCIE (v6): wszystkie testy na wszystkich pulsarach + wybor najlepszego wariantu
# wycinania nulli. pazur2 (v5) NIE jest zmieniany; wszystko w /home/psr/data/OUTPUT/final_v6/.
#
# Uruchom w Julii na kopierniku (moze trwac kilka godzin; wznawialne - zrobione etapy sa pomijane):
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/final_v6.jl"); final_v6()
#
# Etapy:
#   1. energie impulsow wszystkich pulsarow: 3 linie bazowe (srednia / Fourier 1 / Fourier 2 harmoniki,
#      w kazdym pulsarze wybierana ta z najmniejszym rozrzutem E_off), faza z filtra dopasowanego (CHIME),
#      dolna granica NF (CHIME/Esamdin)
#   2. model nulli (TPA XIX) dla wszystkich: sigma empiryczne, bootstrap blokowy, P_null kazdego impulsu
#   3. testy syntetyczne ze znana prawda: benchmark wycinania, 8 wariantow regul nulli, czulosc progow
#   4. skan wszystkich pulsarow: 8 wariantow wycinania nulli na tych samych danych (nic nie zapisuje)
#   5. ranking wariantow (testy syntetyczne + prawdziwe dane + metoda reczna) -> zalecany wariant
#   6. wersja koncowa zalecanego wariantu -> final_v6/pazur2_final/ (+ weryfikacja, porownanie z metoda
#      reczna, model nulli z koncowymi etykietami)
#   7. niepewnosc NF: wrazliwosc na dlugosc bloku (10 / 50) i nested sampling (dynesty) na wszystkich
#      pulsarach z wiarygodnym NF - najdluzszy etap, na koncu
# Potem (po decyzji): final_v6_adopt()  ->  pazur2 -> pazur2_przed_v6, final_v6/pazur2_final -> pazur2
# Inny wariant zamiast zalecanego:  final_v6_build("automat bez marg. + P>=0.95")
using Printf, Statistics, LinearAlgebra, FFTW, JSON

F6_ROOT = "/home/psr/data/OUTPUT/"
F6_DIR = F6_ROOT * "final_v6/"
F6_DATA = "/home/psr/data/new/"
# _f6_load potrzebuje _debased_spcf (modules/zap_verify.jl) - wczytywane od razu, zeby dzialalo w kazdej sesji
isdefined(Main, :_debased_spcf) || Base.include(Main, joinpath(@__DIR__, "..", "modules", "zap_verify.jl"))
F6_GP = Ref(false)          # v7: linia bazowa GP jako dodatkowa opcja (domyslnie jak v6)
# nazwa => (margines przy nullach, tryb, prog P_null). Tryb :auto = nulle automatu, :model = JEDNA REGULA
# (tylko model energii, jak w publikacji), :union = automat + model. Bez wiarygodnego modelu -> automat.
F6_VARIANTS = [
    ("automat (v5)", true, :auto, 0.0),
    ("automat bez marg.", false, :auto, 0.0),
    ("model P>=0.5 (jedna reguła)", true, :model, 0.5),
    ("model P>=0.8", true, :model, 0.8),
    ("model P>=0.95", true, :model, 0.95),
    ("automat + P>=0.95", true, :union, 0.95),
    ("automat bez marg. + P>=0.95", false, :union, 0.95),
    ("automat bez marg. + P>=0.8", false, :union, 0.8),
]
F6_SCAN_HEADER = "pulsar,variant,N,n_cut,n_null,snr_none,snr,snr_manual,null_E_frac,null_P_lt005,null_phase_in,phase_random,cut_NF,NF_model,status"

_f6_mark(s) = F6_DIR * ".etap_$s"
_f6_done(s) = isfile(_f6_mark(s))
_f6_set(s) = (mkpath(F6_DIR); write(_f6_mark(s), Libc.strftime("%Y-%m-%d %H:%M", time())))

function _f6_psrs()
    st = Dict{String,String}()
    for l in Iterators.drop(eachline(F6_ROOT * "pazur2/zap_summary.csv"), 1)
        r = split(l, ","); length(r) >= 2 && (st[r[1]] = r[2])
    end
    sort([k for (k, v) in st if v in ("ok", "reczny")])
end

F6_CACHE = F6_ROOT * "cache_debased/"     # stosy impulsow po odjeciu tla (pamiec podreczna, dokladna kopia)

"""
Klucz pamieci podrecznej: obserwacja, okno impulsu, kanaly wyciete przed sumowaniem, nazwy, rozmiary i czasy
modyfikacji plikow zrodlowych - kazda zmiana danych albo ustawien daje nowy klucz (stara kopia nieuzywana).
"""
function _f6_key(psr, p, obsdir)
    fs = sort(readdir(obsdir))
    src = join(["$f:$(filesize(joinpath(obsdir, f))):$(mtime(joinpath(obsdir, f)))" for f in fs], ";")
    "$psr|$(basename(obsdir))|$(Int(p["bin_st"]))-$(Int(p["bin_end"]))|$(something(get(p, "zap_chans16", nothing), ""))|$src"
end

function _f6_load(psr, p; cache=true)
    obs = sort(filter(d -> isdir(joinpath(F6_DATA, psr, d)), readdir(joinpath(F6_DATA, psr))))
    obsdir = joinpath(F6_DATA, psr, obs[1])
    # pamiec podreczna: raz przerobione dane (PSRCHIVE: suma, kanaly, odjecie tla, konwersja do tekstu) zapisane
    # binarnie - kolejne testy wczytuja je w ulamku sekundy; te same liczby co bez pamieci (zapis bit w bit)
    key = cache ? _f6_key(psr, p, obsdir) : ""
    cf = F6_CACHE * psr * ".bin"; kf = F6_CACHE * psr * ".key"
    if cache && isfile(cf) && isfile(kf) && read(kf, String) == key
        try
            return open(cf) do io
                T = read(io, UInt8) == 0x04 ? Float32 : Float64
                N = read(io, Int64); nb = read(io, Int64)
                M = Matrix{T}(undef, N, nb); read!(io, M); M
            end
        catch e
            println("pamięć podręczna $psr nieczytelna - liczę od nowa: ", first(sprint(showerror, e), 100))
        end
    end
    D = _f6_load_raw(psr, p, obsdir)
    if cache && D isa AbstractMatrix && eltype(D) in (Float32, Float64)
        try
            mkpath(F6_CACHE)
            tmpf = cf * ".tmp$(getpid())"
            open(tmpf, "w") do io
                write(io, eltype(D) == Float32 ? 0x04 : 0x08); write(io, Int64(size(D, 1)), Int64(size(D, 2))); write(io, D)
            end
            mv(tmpf, cf; force=true); write(kf, key)
        catch e
            println("nie zapisano pamięci podręcznej $psr: ", first(sprint(showerror, e), 100))
        end
    end
    return D
end

function _f6_load_raw(psr, p, obsdir)
    tmp = mktempdir()
    try
        return Main._debased_spcf(obsdir * "/", tmp, Int(p["bin_st"]), Int(p["bin_end"]);
                                  chans16=something(get(p, "zap_chans16", nothing), ""))
    finally
        rm(tmp; recursive=true, force=true)
    end
end

# ---------------------------------------------------------------- etap 1: energie
"Energie jednego pulsara (etykiety z `labels_root`/<psr>/params.json). Zwraca \"\" albo powód pominięcia."
function _f6_energy_one(psr; labels_root=F6_ROOT * "pazur2/", out=F6_DIR * "model/")
    _VD = SpaTs.Data
    p = SpaTs.Tools.read_params(joinpath(labels_root, psr, "params.json"))
    D = _f6_load(psr, p)
    N, nb = size(D)
    bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
    w = Main.ne_windows(D, bst, bend)
    isnothing(w) && return "okno on-pulse za szerokie na dwa rozłączne okna off-pulse"
    on, W, B = w
    lab = zeros(Int, N)
    lab[_VD.zaps_to_rows(_VD.zaps_rfi_of(p), N)] .= 1
    lab[_VD.zaps_to_rows(_VD.zaps_null_of(p), N)] .= 2
    lab[_VD.zaps_to_rows(_VD.zaps_other_of(p), N)] .= 3
    good = [lab[i] != 1 && any(!iszero, @view D[i, :]) for i in 1:N]
    count(good) < 20 && return "za mało impulsów"
    e = _f6_energy_core(D, on, W, B, lab, good)
    mkpath(joinpath(out, "energies"))
    open(joinpath(out, "energies", "$psr.csv"), "w") do io
        println(io, "# pulsar=$psr\n# N=$N\n# n_on=$(length(on))\n# n_B=$(length(B))\n# sigma_bin=$(e.σb)\n# bin_st=$bst\n# bin_end=$bend")
        println(io, "# W=$(first(W))-$(last(W))\n# nbin=$nb\n# baseline=$(e.name)\n# std_Eoff=$(join(e.sds, " "))\n# NF_low=$(e.nflow)")
        println(io, "i,E_on,E_off,label,dphi")
        for i in 1:N
            any(!iszero, @view D[i, :]) || continue
            println(io, "$(i - 1),$(e.Eon[i]),$(e.Eoff[i]),$(lab[i]),$(e.dphi[i])")
        end
    end
    return ""
end

"Energie, wybór linii bazowej, faza z filtra dopasowanego i dolna granica NF (bez plików - testowalne)."
function _f6_energy_core(D, on, W, B, lab, good)
    N, nb = size(D)
    φ = 2π .* (0:nb-1) ./ nb
    Xs = [("srednia", ones(nb, 1)), ("fourier1", hcat(ones(nb), cos.(φ), sin.(φ))),
          ("fourier2", hcat(ones(nb), cos.(φ), sin.(φ), cos.(2φ), sin.(2φ)))]
    # linia bazowa impulsu = P * (dane w binach B), P liniowe (nb x |B|): Fourier = rzut na harmoniki,
    # GP (F6_GP[]) = srednia warunkowa procesu gaussowskiego (jadro SE na okregu fazy + stala) ze stalymi
    # parametrami - ten sam operator dla wszystkich impulsow pulsara, wiec szybko. Wybor jedna regula:
    # najmniejszy rozrzut E_off w oknie W, ktore nie bierze udzialu w dopasowaniu.
    Ps = [(name, X * pinv(X[B, :])) for (name, X) in Xs]
    if F6_GP[]
        dist = [min(abs(i - j), nb - abs(i - j)) for i in 1:nb, j in 1:nb]
        for ell in (nb / 16, nb / 8, nb / 4), lam in (0.3, 3.0)
            K = exp.(-dist .^ 2 ./ (2ell^2)) .+ 10.0
            push!(Ps, (@sprintf("gp_l%d_s%.1f", round(Int, ell), lam), K[:, B] / (K[B, B] + lam * I)))
        end
    end
    best, sds = nothing, String[]
    DB = D[:, B]
    for (name, P) in Ps
        pon, pw = vec(sum(P[on, :], dims=1)), vec(sum(P[W, :], dims=1))
        Eon = vec(sum(D[:, on], dims=2)) .- DB * pon
        Eoff = vec(sum(D[:, W], dims=2)) .- DB * pw
        sd = std(Eoff[good])
        push!(sds, @sprintf("%s=%.4g", name, sd))
        (isnothing(best) || sd < best[2]) && (best = (name, sd, Eon, Eoff, P))
    end
    name, sdE, Eon, Eoff, P = best
    R = D .- DB * P'                                       # dane po odjeciu linii bazowej
    σb = std(vec(R[good, W]))
    # faza z filtra dopasowanego (CHIME, Ng+2020): przesuniecie maksimum korelacji z szablonem
    tmpl = vec(mean(R[good .& (lab .== 0), :], dims=1))
    Tf = conj.(rfft(tmpl))
    CC = irfft(rfft(R, 2) .* reshape(Tf, 1, :), nb, 2)
    dphi = [(k = argmax(@view CC[i, :]) - 1; k > nb ÷ 2 ? k - nb : k) for i in 1:N]
    # dolna granica NF (CHIME/Esamdin): najwiecej najslabszych impulsow, ktorych suma nie ma emisji (< 3 sigma)
    Eg = sort(Eon[good]); n = length(Eg)
    sn = cumsum(Eg) ./ (sqrt.(1:n) .* sdE)
    k = findlast(<(3), sn)
    nflow = isnothing(k) ? 0.0 : k / n
    return (name=name, sds=sds, Eon=Eon, Eoff=Eoff, σb=σb, dphi=dphi, nflow=nflow, sdE=sdE)
end

function _f6_energies_part(k, n; labels_root=F6_ROOT * "pazur2/", out=F6_DIR * "model/", force=false)
    psrs = _f6_psrs()[k:n:end]
    for (j, psr) in enumerate(psrs)
        f = joinpath(out, "energies", "$psr.csv")
        (!force && isfile(f)) && continue
        t0 = time()
        msg = try
            _f6_energy_one(psr; labels_root=labels_root, out=out)
        catch e
            e isa InterruptException && rethrow()
            "BŁĄD: " * first(sprint(showerror, e), 150)
        end
        isempty(msg) || open(io -> println(io, "$psr,$(replace(msg, ',' => ';'))"), joinpath(out, "skipped_$k.txt"), "a")
        @printf("[%d/%d] %-12s %s (%.0f s)\n", j, length(psrs), psr, isempty(msg) ? "ok" : msg, time() - t0)
        flush(stdout)
    end
end

function _f6_parallel(fn, n, tag; args="", extra="", logs=F6_DIR * "logs")
    pre = "include(\"spats.jl\"); Base.include(Main, \"modules/zap_verify.jl\"); ENV[\"NULLING_OUT\"]=\"$(F6_DIR)model/\"; " *
          "Base.include(Main, \"tools/nulling_energies.jl\"); include(\"tools/final_v6.jl\"); " * extra
    codes = [pre * "Base.invokelatest($fn, $k, $n$args)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, logs, tag; every=120)
end

# ---------------------------------------------------------------- etap 4: skan wariantow
"P_null z modelu (wiersze od 1) i status modelu dla pulsara."
function _f6_model(psr)
    st = Dict{String,String}()
    f = F6_DIR * "model/nulls_model.csv"
    isfile(f) && for l in Iterators.drop(eachline(f), 1)
        r = split(l, ","); st[r[1]] = r[2]
    end
    pf = F6_DIR * "model/nulls_model/$psr.csv"
    isfile(pf) || return nothing, get(st, psr, "brak")
    rows = Int[]; P = Float64[]
    for l in Iterators.drop(eachline(pf), 1)
        r = split(l, ","); push!(rows, parse(Int, r[1]) + 1); push!(P, parse(Float64, r[3]))
    end
    return (rows=rows, P=P), get(st, psr, "brak")
end

"Lista pewnych nulli z modelu (P_null >= prog) dla wszystkich pulsarów z wiarygodnym modelem."
function _f6_extra(thr)
    out = Dict{String,Vector{Int}}()
    for psr in _f6_psrs()
        m, st = _f6_model(psr)
        (isnothing(m) || st != "ok") && continue
        out[psr] = m.rows[m.P .>= thr]
    end
    out
end

function _f6_scan_part(k, n)
    _VD = SpaTs.Data
    out = F6_DIR * "scan/"; mkpath(out)
    f = out * "scan_$k.csv"
    done = isfile(f) ? Set(split(l, ",")[1] for l in Iterators.drop(eachline(f), 1)) : Set{String}()
    isfile(f) || open(io -> println(io, F6_SCAN_HEADER), f, "w")
    nfm = Dict{String,Float64}()
    sf = F6_DIR * "model/summary.csv"
    if isfile(sf)
        L = readlines(sf); H = split(L[1], ","); ix = findfirst(==("NF_model"), H)
        for l in L[2:end]; r = split(l, ","); nfm[r[1]] = something(tryparse(Float64, r[ix]), NaN); end
    end
    psrs = filter(!in(done), _f6_psrs()[k:n:end])
    for (j, psr) in enumerate(psrs)
        t0 = time()
        try
            p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
            c = get(p, "zap_cmp", Dict())
            D = _f6_load(psr, p)
            N = size(D, 1)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            on = bst:bend
            off = _VD.offpulse_bins(D, bst, bend)
            base = [mean(@view D[i, off]) for i in 1:N]
            E = [sum(@view D[i, on]) - length(on) * base[i] for i in 1:N]
            m, st = _f6_model(psr)
            Pr = fill(NaN, N)
            isnothing(m) || (Pr[m.rows] .= m.P)
            dph = fill(NaN, N); nbin = size(D, 2)
            ef = F6_DIR * "model/energies/$psr.csv"
            if isfile(ef)
                for l in eachline(ef)
                    (startswith(l, "#") || startswith(l, "i,")) && continue
                    r = split(l, ","); length(r) >= 5 && (dph[parse(Int, r[1]) + 1] = parse(Float64, r[5]))
                end
            end
            hw = length(on) / 2
            snr_none = _VD.profile_snr(D, bst, bend)
            for (name, marg, mode, thr) in F6_VARIANTS
                extra = (mode == :auto || isnothing(m) || st != "ok") ? nothing : m.rows[m.P .>= thr]
                a = _VD.auto_zap(D, bst, bend; null_margin=marg, null_mode=mode, null_extra=extra)
                nul = findall(x -> x in (3, 6), a.reason)
                kept = setdiff(1:N, a.bad)
                rfi = _VD.split_zaps(a.reason).rfi
                nonrfi = N - length(rfi)
                ef_ = length(nul) >= 5 && !isempty(kept) ? mean(E[nul]) / mean(E[kept]) : NaN
                pe = length(nul) >= 5 && any(isfinite, Pr[nul]) ? mean(filter(isfinite, Pr[nul]) .< 0.05) : NaN
                ph = length(nul) >= 5 && any(isfinite, dph[nul]) ? mean(abs.(filter(isfinite, dph[nul])) .<= hw) : NaN
                row = [psr, name, N, length(a.bad), length(nul), snr_none, _VD.profile_snr(D, bst, bend; drop=a.bad),
                       something(get(c, "snr_maciej", nothing), NaN), ef_, pe, ph, 2hw / nbin,
                       nonrfi > 0 ? length(nul) / nonrfi : NaN, get(nfm, psr, NaN), st]
                open(io -> println(io, join([x isa AbstractFloat ? @sprintf("%.5g", x) : string(x) for x in row], ",")), f, "a")
            end
            @printf("[%d/%d] %-12s %d wariantów (%.0f s)\n", j, length(psrs), psr, length(F6_VARIANTS), time() - t0)
        catch e
            e isa InterruptException && rethrow()
            println("[$j/$(length(psrs))] $psr BŁĄD: ", first(sprint(showerror, e), 200))
        end
        flush(stdout)
    end
end

# ---------------------------------------------------------------- etap 5: ranking
function f6_ranking(; tie=1.0)
    H = split(F6_SCAN_HEADER, ",")
    R = Dict{String,Vector{Dict{String,String}}}()
    for f in filter(f -> occursin(r"^scan_\d+\.csv$", basename(f)), readdir(F6_DIR * "scan/"; join=true))
        for l in Iterators.drop(eachline(f), 1)
            v = split(l, ","); length(v) == length(H) || continue
            push!(get!(R, v[2], Dict{String,String}[]), Dict(zip(H, v)))
        end
    end
    syn = isfile(F6_DIR * "tests/null_rule/null_rule.json") ? JSON.parsefile(F6_DIR * "tests/null_rule/null_rule.json") : Dict()
    num(r, k) = something(tryparse(Float64, r[k]), NaN)
    fin(x) = filter(isfinite, x)
    L = String[]; pr(s) = (println(s); push!(L, s))
    pr("# Ranking wariantów wycinania nulli (v6)\n")
    pr("Wszystkie warianty: te same dane, te same wycięcia RFI (automat v5); różni je tylko to, które impulsy są wycinane jako nulle.\n")
    pr("| wariant | syntetyczne: kompletność | emisja wycięta | max bias profilu | lepszy : gorszy : remis vs ręcznie | średni zysk S/N | wycięte | energia nulli (mediana) | nulle z P<0.05 | faza nulli w oknie (losowo) | |NF wycięte − NF model| |")
    pr("|---|---|---|---|---|---|---|---|---|---|---|")
    score = Dict{String,Float64}()
    for (name, _, mode, _) in F6_VARIANTS
        rr = get(R, name, Dict{String,String}[])
        isempty(rr) && continue
        s = get(syn, name, Dict()); sa = get(s, "_srednio", Dict())
        k_ = get(sa, "kompletnosc", NaN); e_ = get(sa, "emisja_wycieta", NaN); b_ = get(sa, "bias_max", NaN)
        nb = nw = ne = 0; g = Float64[]
        for r in rr
            s0, s1, sm = num(r, "snr_none"), num(r, "snr"), num(r, "snr_manual")
            isfinite(sm) || continue
            ga, gm = 100 * (s1 - s0) / abs(s0), 100 * (sm - s0) / abs(s0)
            ga - gm >= tie ? (nb += 1) : gm - ga >= tie ? (nw += 1) : (ne += 1)
            s0 >= 3 && push!(g, ga)
        end
        cutp = 100 * mean(num(r, "n_cut") / num(r, "N") for r in rr)
        efr = fin([num(r, "null_E_frac") for r in rr])
        pe = fin([num(r, "null_P_lt005") for r in rr])
        ph = fin([num(r, "null_phase_in") for r in rr]); phr = fin([num(r, "phase_random") for r in rr])
        dnf = fin([abs(num(r, "cut_NF") - num(r, "NF_model")) for r in rr if r["status"] == "ok"])
        md(x) = isempty(x) ? NaN : median(x)
        pr(@sprintf("| %s | %.1f%% | %.2f%% | %+.2f%% | %d : %d : %d | %+.1f%% | %.1f%% | %.1f%% | %.1f%% | %.2f (%.2f) | %.1f pp |",
                    name, 100 * k_, 100 * e_, 100 * b_, nb, nw, ne, isempty(g) ? NaN : mean(g), cutp, 100 * md(efr),
                    100 * md(pe), md(ph), md(phr), 100 * md(dnf)))
        # kryterium: syntetyczne (prawda znana) - profil bez znieksztalcen (bias <= 1%, emisja <= 1%);
        # wsrod nich najwieksza kompletnosc nulli; jedna regula wygrywa przy roznicy < 1 pp
        # + prawdziwe dane: wyciete nulle bez emisji (energia usrednionych nulli, mediana < 2% dobrego impulsu)
        ok = isfinite(b_) && b_ <= 0.01 && e_ <= 0.01 && !(md(efr) >= 0.02)
        score[name] = ok ? k_ + (mode == :model ? 0.01 : 0.0) : -1.0     # jedna regula: premia 1 pp
    end
    best = isempty(score) ? "automat (v5)" : first(sort(collect(score), by=x -> -x[2]))[1]
    pr("\n**Kryterium wyboru:** tylko warianty, które na danych syntetycznych ze znaną prawdą nie zniekształcają profilu " *
       "(bias ≤ 1%, wycięta emisja ≤ 1%); spośród nich największa kompletność nulli. Jedna reguła (tylko model, jak w publikacji) wygrywa " *
       "przy różnicy kompletności poniżej 1 pp. Na prawdziwych danych energia uśrednionych nulli musi być < 2%.")
    pr("\n**Zalecany wariant: $best**")
    write(F6_DIR * "RAPORT_ranking.md", join(L, "\n") * "\n")
    return best
end

# ---------------------------------------------------------------- etap 6: wersja koncowa
"Wersja końcowa wybranego wariantu w final_v6/pazur2_final/ (pazur2 bez zmian) + weryfikacja i porównania."
function final_v6_build(variant=nothing; n=32)
    isnothing(variant) && (variant = f6_ranking())
    i = findfirst(v -> v[1] == variant, F6_VARIANTS)
    isnothing(i) && error("nieznany wariant: $variant; dostępne: $(join(first.(F6_VARIANTS), " | "))")
    _, marg, mode, thr = F6_VARIANTS[i]
    out = F6_DIR * "pazur2_final/"
    if isdir(out) && isfile(out * ".wariant") && strip(read(out * ".wariant", String)) != variant
        error("$out zawiera inny wariant ($(strip(read(out * ".wariant", String)))) - usuń go albo przenieś przed budową innego")
    end
    mkpath(out); write(out * ".wariant", variant)
    println("\n===== 6. WERSJA KOŃCOWA: $variant → $out =====")
    psrs = _f6_psrs()
    cfg = "SpaTs.Data.NULL_MARGIN[] = $marg; SpaTs.Data.NULL_MODE[] = :$mode; " *
          (mode == :auto ? "" : "merge!(SpaTs.Data.NULL_EXTRA, _f6_extra($thr)); ")
    pre = "include(\"spats.jl\"); Base.include(Main, \"modules/zap_verify.jl\"); include(\"tools/final_v6.jl\"); " * cfg
    cf = "SpaTs.chan_first_list([\"$(F6_ROOT)zap_verify_v3/\", \"$(F6_ROOT)zap_verify/\"])"
    codes = [pre * "SpaTs.autozap_all(out=\"$out\", only=$(repr(psrs[k:n:end])), skip_done=true, keep_manual=false, " *
             "reuse_window=true, window_root=\"$(F6_ROOT)pazur2/\", chan_first=$cf)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, F6_DIR * "logs", "build"; every=120)
    vd = F6_DIR * "zap_verify_final/"
    codes = [pre * "Base.invokelatest(Main.zap_verify; pazur2=\"$out\", out=\"$vd\", part=\"$k/$n\", chans=false, save_rfi=true)"
             for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, F6_DIR * "logs", "verify"; every=120)
    Base.invokelatest(Main.zap_verify_summary; out=vd)
    txt = mktemp() do path, io                     # wynik na ekran i do pliku (do WYNIK.txt)
        redirect_stdout(io) do
            SpaTs.compare_zap_methods(mine=out, outdir=F6_DIR * "zap_comparison_final/")
        end
        flush(io); read(path, String)
    end
    print(txt); mkpath(F6_DIR * "zap_comparison_final/")
    write(F6_DIR * "zap_comparison_final/summary.txt", txt)
    # model nulli z koncowymi etykietami (te same impulsy RFI -> bez ponownego dopasowania)
    _f6_parallel("_f6_energies_part", n, "energie_final"; args=", labels_root=\"$out\", force=true")
    try
        run(addenv(`python3 tools/nulling_model.py`, "NULLING_OUT" => F6_DIR * "model/", "NULLING_RELABEL" => "1"))
    catch e
        println("model nulli: ", e)
    end
end

"Przyjęcie wersji końcowej: pazur2 → pazur2_przed_v6, final_v6/pazur2_final → pazur2."
function final_v6_adopt()
    src, dst, bk = F6_DIR * "pazur2_final", F6_ROOT * "pazur2", F6_ROOT * "pazur2_przed_v6"
    isdir(src) || error("brak $src - najpierw final_v6() albo final_v6_build()")
    k = 1
    while isdir(bk)                      # nic nie nadpisujemy: kolejna kopia dostaje numer
        k += 1; bk = F6_ROOT * "pazur2_przed_v6_$k"
    end
    v = strip(read(src * "/.wariant", String))
    mv(dst, bk); mv(src, dst)
    println("Przyjęte: $v. Stara wersja: $bk. Domyślnie w kodzie: Data.NULL_MODE[] = :union, NULL_THR[] = 0.95 (v5: NULL_MODE[] = :auto).")
end

# ---------------------------------------------------------------- calosc
function final_v6(; n=32, dynesty=true, build=true)
    t0 = time()
    mkpath(F6_DIR)
    Base.include(Main, "modules/zap_verify.jl")
    py(args...; env=Dict{String,String}()) = try
        run(addenv(`python3 $args`, merge(Dict("NULLING_OUT" => F6_DIR * "model/"), env)...)); true
    catch e
        println("python: ", e); false
    end
    if !_f6_done("1_energie")
        println("\n===== 1. ENERGIE wszystkich pulsarów ($n procesów) =====")
        _f6_parallel("_f6_energies_part", n, "energie")
        _f6_set("1_energie")
    end
    if !_f6_done("2_model")
        println("\n===== 2. MODEL NULLI (wszystkie pulsary; sigma empiryczne, bootstrap blokowy) =====")
        ctl = String[]
        for f in readdir(F6_DIR * "model/energies"; join=true)
            any(l -> !startswith(l, "#") && !startswith(l, "i,") && split(l, ",")[4] == "2", eachline(f)) ||
                push!(ctl, basename(f)[1:end-4])
        end
        write(F6_DIR * "model/control.txt", join(ctl, "\n") * "\n")
        py("tools/nulling_model.py")
        _f6_set("2_model")
    end
    if !_f6_done("3_testy")
        println("\n===== 3. TESTY SYNTETYCZNE =====")
        td = F6_DIR * "tests/"; mkpath(td)
        for (t, cmd) in [("benchmark", `julia --project=. tests/zap_benchmark.jl`),
                         ("czulosc_progow", `julia --project=. tests/threshold_sensitivity.jl`)]
            try
                write(td * "$t.txt", read(cmd, String)); println("  $t: zapisane $(td)$t.txt")
            catch e
                println("  $t: ", e)
            end
        end
        try
            run(`julia --project=. tests/null_rule_gen.jl $(td)null_rule`)
            run(addenv(`python3 tests/null_rule_eval.py $(td)null_rule`, "NR_OUT" => td * "null_rule"))
        catch e
            println("  reguły nulli: ", e)
        end
        _f6_set("3_testy")
    end
    if !_f6_done("4_skan")
        println("\n===== 4. SKAN: $(length(F6_VARIANTS)) wariantów wycinania nulli na wszystkich pulsarach =====")
        _f6_parallel("_f6_scan_part", n, "skan")
        _f6_set("4_skan")
    end
    println("\n===== 5. RANKING =====")
    best = f6_ranking()
    if build && !_f6_done("6_build")
        final_v6_build(best; n=n)
        _f6_set("6_build")
    end
    if !_f6_done("7_niepewnosc")
        println("\n===== 7. NIEPEWNOŚĆ NF: długość bloku 10 / 50 =====")
        cdir = F6_DIR * "model/checks/"
        for b in ("10", "50")
            py("tools/nulling_checks.py"; env=Dict("NC_BLOCK" => b))
            isfile(cdir * "RAPORT_checks.md") && cp(cdir * "RAPORT_checks.md", cdir * "RAPORT_checks_blok$b.md"; force=true)
        end
        _f6_set("7_niepewnosc")
    end
    if dynesty && !_f6_done("8_dynesty")
        println("\n===== 8. NESTED SAMPLING (dynesty): wszystkie pulsary z wiarygodnym NF, 64 procesy =====")
        py("tools/nulling_dynesty.py", "--all"; env=Dict("NULLING_NPROC" => "64")) ? _f6_set("8_dynesty") :
            println("ETAP 8 NIEUDANY - nie oznaczam jako zrobiony; po naprawie uruchom final_v6() ponownie")
    end
    if isfile("input/tpa19_nf.csv")
        println("\n===== 9. PORÓWNANIE Z TPA XIX (Keith+2026) =====")
        py("tools/compare_tpa19.py")
    end
    # wszystko do wklejenia w jednym pliku
    rd(f, n=10^6) = isfile(f) ? join(first(readlines(f), n), "\n") : "(brak $f)"
    md = F6_DIR * "model/"
    wyn = join(["=== RANKING", rd(F6_DIR * "RAPORT_ranking.md"),
                "=== MODEL", join(filter(l -> startswith(l, "- "), split(rd(md * "RAPORT_nulle.md"), "\n")), "\n"),
                "=== BLOK 10", rd(md * "checks/RAPORT_checks_blok10.md", 30), "=== BLOK 50", rd(md * "checks/RAPORT_checks_blok50.md", 30),
                "=== DYNESTY", join(last(split(rd(md * "dynesty/RAPORT_dynesty.md"), "\n"), 4), "\n"),
                "=== PORÓWNANIE Z METODĄ RĘCZNĄ", rd(F6_DIR * "zap_comparison_final/summary.txt", 20),
                "=== TPA XIX", rd(md * "RAPORT_tpa19.md", 70)], "\n")
    write(F6_DIR * "WYNIK.txt", wyn * "\n")
    println(@sprintf("\nKONIEC final_v6 (%.1f h). pazur2 bez zmian.", (time() - t0) / 3600))
    println("Do wklejenia (całość): cat $(F6_DIR)WYNIK.txt  oraz blok AUTOMAT vs METODA RĘCZNA z etapu 6, jeśli był.")
end

"""Ponowne przejście po poprawce f_i3 (2026-11-18: f_i3 liczone z rozkładu emisji bez nulli, jak w TPA XIX;
wcześniej razem z nullami, co zaniżało f_i3 pulsarom z dużym NF). Energie i testy syntetyczne zostają
(nie zależą od f_i3); od nowa: model nulli (bez ponownego dopasowania - tylko f_i3), skan wariantów, ranking,
wersja końcowa (final_v6/pazur2_final), bootstrap blokowy, dynesty. Stary skan -> final_v6/scan_przed_fi3/.
pazur2 bez zmian; potem final_v6_adopt()."""
function final_v6_redo(; kw...)
    sc = F6_DIR * "scan/"
    if isdir(sc)
        dst = F6_DIR * "scan_przed_fi3/"
        isdir(dst) || mv(sc, dst)
        isdir(sc) && rm(sc; recursive=true)
    end
    for s in ("2_model", "4_skan", "6_build", "7_niepewnosc", "8_dynesty")
        rm(_f6_mark(s); force=true)
    end
    final_v6(; kw...)
end
