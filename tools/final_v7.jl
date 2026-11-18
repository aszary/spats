# PRZEJSCIE v7: ostatnie testy automatu. pazur2 (v6) NIE jest zmieniany; wszystko w /home/psr/data/OUTPUT/final_v7/.
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/final_v7.jl"); final_v7()
# Wznawialne (zrobione etapy sa pomijane). Etapy:
#   1. kroki automatu, dane syntetyczne: kazdy krok wylaczony osobno (tests/zap_ablation.jl)
#   2. kroki automatu, wszystkie pulsary: to samo na prawdziwych danych (+ metoda reczna)
#   3. wersja uproszczona: wszystkie zbedne kroki wylaczone naraz - sprawdzenie syntetyczne i na wszystkich
#   4. linia bazowa GP (jak w TPA XIX) jako dodatkowa opcja: energie, model nulli, porownanie z v6 i TPA XIX
#   5. nulle w 16 kanalach: pewne nulle musza miec zerowa energie w KAZDYM kanale
#   6. decyzja (jedna regula: krok zostaje, jesli bez niego cos sie pogarsza) -> WYNIK_v7.txt
# Potem (po decyzji): final_v7_build() -> final_v7/pazur2_final, final_v7_adopt().
include(joinpath(@__DIR__, "final_v6.jl"))
include(joinpath(@__DIR__, "v7_steps.jl"))
F7_DIR = F6_ROOT * "final_v7/"
F7_SCAN_HEADER = "pulsar,variant,N,n_cut,n_null,snr_none,snr,snr_manual,snr_rfi,kept_rfi3,null_E_frac"
_f7_mark(s) = F7_DIR * ".etap_$s"
_f7_done(s) = isfile(_f7_mark(s))
_f7_set(s) = (mkpath(F7_DIR); write(_f7_mark(s), Libc.strftime("%Y-%m-%d %H:%M", time())))
const F7_PRE = "include(\"spats.jl\"); Base.include(Main, \"modules/zap_verify.jl\"); include(\"tools/final_v7.jl\"); "

# ---------------------------------------------------------------- etap 2/3: skan krokow
function _f7_scan_part(k, n; variants=V7_STEPS, tag="scan")
    _VD = SpaTs.Data
    out = F7_DIR * "$tag/"; mkpath(out)
    f = out * "scan_$k.csv"
    done = isfile(f) ? Set(split(l, ",")[1] for l in Iterators.drop(eachline(f), 1)) : Set{String}()
    isfile(f) || open(io -> println(io, F7_SCAN_HEADER), f, "w")
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
            E = [sum(@view D[i, on]) - length(on) * mean(@view D[i, off]) for i in 1:N]
            z = _VD.rfi_rms_z(D, bst, bend)
            extra = _VD.null_extra_for(psr)
            snr_none = _VD.profile_snr(D, bst, bend)
            rows = String[]
            for (name, st) in variants
                a = _VD.auto_zap(D, bst, bend; null_extra=extra, _v7_kw(st)...)
                nul = findall(x -> x in (3, 6), a.reason)
                kept = setdiff(1:N, a.bad)
                rfi = _VD.split_zaps(a.reason).rfi
                ef_ = length(nul) >= 5 && !isempty(kept) ? mean(E[nul]) / mean(E[kept]) : NaN
                row = [psr, name, N, length(a.bad), length(nul), snr_none, _VD.profile_snr(D, bst, bend; drop=a.bad),
                       something(get(c, "snr_maciej", nothing), NaN), _VD.profile_snr(D, bst, bend; drop=rfi),
                       count(i -> z[i] > 3, kept), ef_]
                push!(rows, join([x isa AbstractFloat ? @sprintf("%.6g", x) : string(x) for x in row], ","))
            end
            open(io -> foreach(r -> println(io, r), rows), f, "a")
            @printf("[%d/%d] %-12s %d wariantów (%.0f s)\n", j, length(psrs), psr, length(variants), time() - t0)
        catch e
            e isa InterruptException && rethrow()
            println("[$j/$(length(psrs))] $psr BŁĄD: ", first(sprint(showerror, e), 200))
        end
        flush(stdout)
    end
end

function _f7_read_scan(tag)
    H = split(F7_SCAN_HEADER, ",")
    R = Dict{String,Dict{String,Dict{String,Float64}}}()        # wariant => pulsar => kolumny
    d = F7_DIR * "$tag/"
    isdir(d) || return R
    for f in filter(f -> occursin(r"^scan_\d+\.csv$", basename(f)), readdir(d; join=true))
        for l in Iterators.drop(eachline(f), 1)
            v = split(l, ","); length(v) == length(H) || continue
            get!(R, v[2], Dict{String,Dict{String,Float64}}())[v[1]] =
                Dict(H[i] => something(tryparse(Float64, v[i]), NaN) for i in 3:length(H))
        end
    end
    R
end

"Ocena wariantow wzgledem v6 na prawdziwych danych. Zwraca Dict wariant => (zbedny, opis)."
function _f7_judge(R, names, base; pr=println, tie=1.0)
    B = get(R, base, Dict())
    out = Dict{String,Any}()
    pr("| wariant | gorzej niż v6 o >1% S/N | lepiej o >1% | najgorzej | vs ręcznie lepszy : gorszy : remis | średni zysk S/N | wycięte | RFI zostawione (z>3) | S/N danych tylko po RFI (najgorzej) | ocena |")
    pr("|---|---|---|---|---|---|---|---|---|---|")
    for name in names
        V = get(R, name, Dict()); isempty(V) && continue
        nw = nb = 0; worst = 0.0; worst_p = ""; b = w = t = 0; g = Float64[]; cut = Float64[]; kr = kr0 = 0.0; wr = 0.0
        for (psr, r) in V
            r0 = get(B, psr, nothing); isnothing(r0) && continue
            dS = 100 * (r["snr"] / r0["snr"] - 1)
            dS < -tie && (nw += 1); dS > tie && (nb += 1)
            dS < worst && (worst = dS; worst_p = psr)
            s0, sm = r["snr_none"], r["snr_manual"]
            if isfinite(sm)
                ga, gm = 100 * (r["snr"] - s0) / abs(s0), 100 * (sm - s0) / abs(s0)
                ga - gm >= tie ? (b += 1) : gm - ga >= tie ? (w += 1) : (t += 1)
                s0 >= 3 && push!(g, ga)
            end
            push!(cut, r["n_cut"] / r["N"]); kr += r["kept_rfi3"]; kr0 += r0["kept_rfi3"]
            wr = min(wr, 100 * (r["snr_rfi"] / r0["snr_rfi"] - 1))
        end
        dkr = kr0 > 0 ? 100 * (kr / kr0 - 1) : 0.0
        ok = name != base && nw <= 5 && worst >= -3 && dkr <= 2 && wr >= -1
        out[name] = (zbedny=ok, nw=nw, worst=worst)
        pr(@sprintf("| %s | %d | %d | %+.1f%% %s | %d : %d : %d | %+.1f%% | %.1f%% | %+.1f%% | %+.1f%% | %s |", name, nw, nb, worst,
                    worst_p, b, w, t, isempty(g) ? NaN : mean(g), 100 * mean(cut), dkr, wr,
                    name == base ? "-" : ok ? "bez zmian - można usunąć" : "potrzebny"))
    end
    out
end

# ---------------------------------------------------------------- etap 4: linia bazowa GP
"Porownanie linii bazowej v6 (srednia/Fourier) i v7 (+GP) z naglowkow plikow energii."
function _f7_gp_summary(; pr=println)
    d6, d7 = F6_DIR * "model/energies/", F7_DIR * "model/energies/"
    rat = Float64[]; ngp = 0; n = 0; top = Tuple{Float64,String,String}[]
    for f in readdir(d7)
        endswith(f, ".csv") || continue
        h7 = Dict(m[1] => m[2] for m in (match(r"^# (\w+)=(.*)$", l) for l in Iterators.take(eachline(d7 * f), 12)) if !isnothing(m))
        sds = Dict(split(x, "=")[1] => parse(Float64, split(x, "=")[2]) for x in split(get(h7, "std_Eoff", "")))
        isempty(sds) && continue
        s6 = minimum(v for (k, v) in sds if !startswith(k, "gp"))
        s7 = minimum(values(sds))
        n += 1; push!(rat, s7 / s6)
        startswith(get(h7, "baseline", ""), "gp") && (ngp += 1; push!(top, (s7 / s6, f[1:end-4], h7["baseline"])))
    end
    pr("- pulsarów: $n; GP wybrane (mniejszy rozrzut E_off niż średnia/Fourier) w $ngp")
    isempty(rat) || pr(@sprintf("- rozrzut E_off v7 / v6: mediana %.3f, w pulsarach z GP mediana %.3f, najmniej %.3f",
                               median(rat), isempty(top) ? NaN : median(first.(top)), minimum(rat)))
    sort!(top)
    isempty(top) || pr("- największy zysk: " * join([@sprintf("%s %.2f (%s)", p, r, b) for (r, p, b) in top[1:min(15, end)]], ", "))
    (n=n, ngp=ngp, med=isempty(top) ? 1.0 : median(first.(top)))
end

"Zgodnosc z TPA XIX z pliku tpa19_compare.csv: (n, mediana |roznicy| pp, odsetek w 2 sigma)."
function _f7_tpa_stats(dir)
    f = dir * "tpa19_compare.csv"; isfile(f) || return nothing
    L = readlines(f); H = split(L[1], ","); ix(k) = findfirst(==(k), H)
    d = Float64[]; z = Bool[]
    for l in L[2:end]
        v = split(l, ","); g(k) = something(tryparse(Float64, v[ix(k)]), NaN)
        (g("fi_t") > 0.3 && g("fi_o") >= 0.3 && v[ix("degen")] != "True" && isfinite(g("nf_t")) && isfinite(g("nf_o"))) || continue
        e = hypot(isfinite(g("e_t")) ? g("e_t") : 0.0, isfinite(g("e_o")) ? g("e_o") : 0.0)
        push!(d, g("nf_o") - g("nf_t")); push!(z, abs(g("nf_o") - g("nf_t")) <= 2 * max(e, 1e-3))
    end
    isempty(d) ? nothing : (n=length(d), mad=100 * median(abs.(d)), med=100 * median(d), in2=mean(z))
end

# ---------------------------------------------------------------- etap 5: nulle w 16 kanalach
const F7_CH_HEADER = "pulsar,set,n_null,n_chan,frac_sum,max_frac,chan_max_frac,max_z,chan_max_z,n_z3"
function _f7_chan_part(k, n)
    _VD = SpaTs.Data
    out = F7_DIR * "chan16/"; mkpath(out)
    f = out * "chan_$k.csv"
    done = isfile(f) ? Set(split(l, ",")[1] for l in Iterators.drop(eachline(f), 1)) : Set{String}()
    isfile(f) || open(io -> println(io, F7_CH_HEADER), f, "w")
    for psr in filter(!in(done), _f6_psrs()[k:n:end])
        t0 = time()
        tmp = mktempdir()
        try
            p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
            obs = sort(filter(d -> isdir(joinpath(F6_DATA, psr, d)), readdir(joinpath(F6_DATA, psr))))
            cu = Main._cube16(joinpath(F6_DATA, psr, obs[1]) * "/", tmp)
            isnothing(cu) && (println("$psr: brak .spCf16"); continue)
            d = cu[1]; N, C, nb = size(d)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            on = bst:bend
            S = dropdims(sum(d, dims=2), dims=2)
            off = _VD.offpulse_bins(S, bst, bend)
            E = [sum(@view d[i, c, on]) - length(on) * mean(@view d[i, c, off]) for i in 1:N, c in 1:C]
            m, st = _f6_model(psr)
            cut = falses(N); cut[_VD.zaps_to_rows(_VD.zaps_of(p), N)] .= true
            zc = Set(parse(Int, x) + 1 for x in split(something(get(p, "zap_chans16", nothing), "")) if !isnothing(tryparse(Int, x)))
            sets = Pair{String,Vector{Int}}[]
            if !isnothing(m) && st == "ok" && maximum(m.rows) <= N
                P = fill(NaN, N); P[m.rows] .= m.P
                push!(sets, "model P>=0.95" => findall(x -> x >= 0.95, P))
                emis = findall(i -> P[i] < 0.05 && !cut[i], 1:N)
            else
                emis = findall(.!cut)
            end
            push!(sets, "nulle automatu" => _VD.zaps_to_rows(_VD.zaps_null_of(p), N))
            rows = String[]
            for (nm, nl) in sets
                length(nl) >= 10 || continue
                fr = Float64[]; zz = Float64[]; cs = Int[]
                for c in 1:C
                    c in zc && continue
                    me = mean(E[emis, c]); sn = std(E[nl, c]) / sqrt(length(nl))
                    (me > 0 && sn > 0) || continue
                    push!(fr, mean(E[nl, c]) / me); push!(zz, mean(E[nl, c]) / sn); push!(cs, c - 1)
                end
                isempty(fr) && continue
                fs = sum(E[nl, cc + 1] for cc in cs) |> x -> mean(x) / mean(sum(E[emis, cc + 1] for cc in cs))
                im, iz = argmax(fr), argmax(zz)
                push!(rows, join([psr, nm, length(nl), length(cs), @sprintf("%.4f", fs), @sprintf("%.4f", fr[im]), cs[im],
                                  @sprintf("%.2f", zz[iz]), cs[iz], count(>(3), zz)], ","))
            end
            open(io -> foreach(r -> println(io, r), rows), f, "a")
            @printf("%-12s %d zestawów (%.0f s)\n", psr, length(rows), time() - t0)
        catch e
            e isa InterruptException && rethrow()
            println("$psr BŁĄD: ", first(sprint(showerror, e), 200))
        finally
            rm(tmp; recursive=true, force=true)
        end
        flush(stdout)
    end
end

function _f7_chan_summary(; pr=println)
    d = F7_DIR * "chan16/"
    rows = Dict{String,Vector{Vector{SubString{String}}}}()
    isdir(d) && for f in filter(f -> startswith(basename(f), "chan_"), readdir(d; join=true))
        for l in Iterators.drop(eachline(f), 1); v = split(l, ","); push!(get!(rows, v[2], []), v); end
    end
    for (nm, R) in sort(collect(rows))
        mz = [parse(Float64, r[8]) for r in R]; mf = [parse(Float64, r[6]) for r in R]; fs = [parse(Float64, r[5]) for r in R]
        nc = [parse(Int, r[4]) for r in R]
        # bez emisji max z z ~16 kanalow ~ maksimum 16 N(0,1): > 3.5 w ok. 0.4% pulsarow
        bad = [r for r in R if parse(Float64, r[8]) > 3.5 && parse(Float64, r[6]) > 0.1]
        pr(@sprintf("- %s: %d pulsarów; energia nulli we wszystkich kanałach razem: mediana %+.1f%% emisji; najjaśniejszy kanał: mediana %.1f%%; kanał z emisją w nullach (z > 3.5 i > 10%%): %d pulsarów",
                    nm, length(R), 100 * median(fs), 100 * median(mf), length(bad)))
        isempty(bad) || pr("  - " * join([@sprintf("%s (kanał %s: %.0f%%, z %.1f)", r[1], r[7], 100 * parse(Float64, r[6]), parse(Float64, r[8]))
                                          for r in sort(bad, by=r -> -parse(Float64, r[8]))[1:min(20, end)]], ", "))
    end
end

# ---------------------------------------------------------------- calosc
function final_v7(; n=32)
    t0 = time(); mkpath(F7_DIR)
    Base.include(Main, "modules/zap_verify.jl")
    par(fn, tag; args="", extra="") = Base.invokelatest(Main._run_parallel,
        [F7_PRE * extra * "Base.invokelatest($fn, $k, $n$args)" for k in 1:n], F7_DIR * "logs", tag; every=120)
    py(args...; env=Dict{String,String}()) = try
        run(addenv(`python3 $args`, env...)); true
    catch e
        println("python: ", e); false
    end
    L = String[]; pr(s) = (println(s); push!(L, s))
    if !_f7_done("1_syn")
        println("\n===== 1. KROKI AUTOMATU: dane syntetyczne =====")
        write(F7_DIR * "ablation_syn.txt", read(`julia --project=. tests/zap_ablation.jl $(F7_DIR)ablation_syn.json`, String))
        _f7_set("1_syn")
    end
    if !_f7_done("2_scan")
        println("\n===== 2. KROKI AUTOMATU: wszystkie pulsary ($(length(V7_STEPS)) wariantów) =====")
        par("_f7_scan_part", "scan")
        _f7_set("2_scan")
    end
    syn = JSON.parsefile(F7_DIR * "ablation_syn.json")
    pr("# v7: czy każdy krok automatu jest potrzebny\n")
    pr("Jedna reguła: krok zostaje, jeśli bez niego coś się pogarsza (na danych syntetycznych ze znaną prawdą albo na prawdziwych).\n")
    pr("## Dane syntetyczne (najgorszy scenariusz względem pełnego automatu)\n")
    pr(join(filter(l -> occursin("ocena", l) || occursin(r"(zbędny|potrzebny|^v6)", l),
                   split(read(F7_DIR * "ablation_syn.txt", String), "\n")), "\n"))
    pr("\n## Wszystkie pulsary (zmiana względem v6)\n")
    R = _f7_read_scan("scan")
    names = first.(V7_STEPS)
    J = _f7_judge(R, names, names[1]; pr=pr)
    drop = [(nm, st) for (nm, st) in V7_STEPS[2:end-1]
            if get(J, nm, (zbedny=false,)).zbedny && get(get(syn, nm, Dict()), "zbedny", false)]
    pr("\n**Zbędne (bez zmian na obu rodzajach danych):** " * (isempty(drop) ? "brak - wszystkie kroki są potrzebne" : join(first.(drop), ", ")))
    simple = NamedTuple()
    if !isempty(drop)
        simple = merge((last(x) for x in drop)...)
        write(F7_DIR * "uproszczony.txt", repr(simple))
        if !_f7_done("3_simple")
            println("\n===== 3. WERSJA UPROSZCZONA: $(repr(simple)) =====")
            vs = "[(\"uproszczony\", $(repr(simple)))]"
            write(F7_DIR * "ablation_simple.txt", read(addenv(`julia --project=. tests/zap_ablation.jl $(F7_DIR)ablation_simple.json`,
                                                               "V7_EXTRA" => repr(simple)), String))
            par("_f7_scan_part", "scan_simple"; args=", variants=$vs, tag=\"scan_simple\"")
            _f7_set("3_simple")
        end
        Rs = _f7_read_scan("scan_simple"); Rs[names[1]] = R[names[1]]
        pr("\n### Wszystkie zbędne kroki wyłączone naraz\n")
        js = _f7_judge(Rs, ["uproszczony"], names[1]; pr=pr)
        sj = get(JSON.parsefile(F7_DIR * "ablation_simple.json"), "uproszczony", Dict())
        pr(@sprintf("- syntetyczne: śmieci %+.1f pp, nulle %+.1f pp, fałszywe %+.2f pp, S/N vs ideał %+.2f pp",
                    get(sj, "smieci_pp", NaN), get(sj, "nulle_pp", NaN), get(sj, "falszywe_pp", NaN), get(sj, "vs_ideal_pp", NaN)))
        ok = get(js, "uproszczony", (zbedny=false,)).zbedny && get(sj, "zbedny", false)
        pr(ok ? "\n**Zalecane: automat bez tych kroków (prostszy, wyniki bez zmian).**" :
                "\n**Razem pogarszają wynik - zostawiamy automat v6.**")
        ok || (simple = NamedTuple())
    end
    write(F7_DIR * "zalecane_kroki.txt", repr(simple))
    if !_f7_done("4_gp")
        println("\n===== 4. LINIA BAZOWA GP: energie + model nulli =====")
        md = F7_DIR * "model/"; mkpath(md * "energies")
        cp(F6_DIR * "model/control.txt", md * "control.txt"; force=true)
        par("_f6_energies_part", "energie_gp"; args=", out=\"$md\"",
            extra="ENV[\"NULLING_OUT\"]=\"$md\"; Base.include(Main, \"tools/nulling_energies.jl\"); F6_GP[] = true; ")
        py("tools/nulling_model.py"; env=Dict("NULLING_OUT" => md))
        isfile("input/tpa19_nf.csv") && py("tools/compare_tpa19.py"; env=Dict("NULLING_OUT" => md))
        _f7_set("4_gp")
    end
    pr("\n## Linia bazowa GP (dodatkowa opcja, wybór jedną regułą: najmniejszy rozrzut E_off)\n")
    g = _f7_gp_summary(; pr=pr)
    t6, t7 = _f7_tpa_stats(F6_DIR * "model/"), _f7_tpa_stats(F7_DIR * "model/")
    for (lab, t) in (("v6", t6), ("v7 z GP", t7))
        isnothing(t) || pr(@sprintf("- zgodność z TPA XIX, %s: %d pulsarów, |różnica| mediana %.2f pp, mediana %+.2f pp, w 2σ %.1f%%",
                                    lab, t.n, t.mad, t.med, 100 * t.in2))
    end
    gp_ok = !isnothing(t6) && !isnothing(t7) && g.ngp > 0 && g.med <= 0.95 && t7.in2 >= t6.in2 - 0.005 && t7.mad <= t6.mad + 0.05
    pr(gp_ok ? "\n**Zalecane: GP jako opcja linii bazowej (mniejszy rozrzut, zgodność z TPA XIX nie gorsza).**" :
               "\n**GP nie poprawia wyniku - zostaje linia bazowa v6.**")
    write(F7_DIR * "zalecane_gp.txt", string(gp_ok))
    if !_f7_done("5_chan")
        println("\n===== 5. NULLE W 16 KANAŁACH =====")
        par("_f7_chan_part", "chan16")
        _f7_set("5_chan")
    end
    pr("\n## Nulle w 16 kanałach (prawdziwy null nie ma emisji w żadnym kanale)\n")
    _f7_chan_summary(; pr=pr)
    pr(@sprintf("\nKONIEC final_v7 (%.1f h). pazur2 bez zmian.", (time() - t0) / 3600))
    write(F7_DIR * "WYNIK_v7.txt", join(L, "\n") * "\n")
    println("Do wklejenia: cat $(F7_DIR)WYNIK_v7.txt")
end

"Wersja koncowa v7 (zalecane kroki i linia bazowa) w final_v7/pazur2_final; pazur2 bez zmian."
function final_v7_build(; steps=eval(Meta.parse(read(F7_DIR * "zalecane_kroki.txt", String))),
                        gp=parse(Bool, strip(read(F7_DIR * "zalecane_gp.txt", String))), n=32)
    out = F7_DIR * "pazur2_final/"; mkpath(out)
    md = gp ? F7_DIR * "model/" : F6_DIR * "model/"
    write(out * ".wariant", "kroki=$(repr(steps)) gp=$gp")
    println("\n===== WERSJA KOŃCOWA v7: kroki wyłączone $(repr(steps)), GP $gp, model nulli $md → $out =====")
    cfg = "SpaTs.Data.ZAP_STEPS[] = $(repr(steps)); SpaTs.Data.NULL_MODEL_DIR[] = \"$md\"; "
    cf = "SpaTs.chan_first_list([\"$(F6_ROOT)zap_verify_v3/\", \"$(F6_ROOT)zap_verify/\"])"
    codes = [F7_PRE * cfg * "SpaTs.autozap_all(out=\"$out\", only=$(repr(_f6_psrs()[k:n:end])), skip_done=true, keep_manual=false, " *
             "reuse_window=true, window_root=\"$(F6_ROOT)pazur2/\", chan_first=$cf)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, F7_DIR * "logs", "build"; every=120)
    vd = F7_DIR * "zap_verify_final/"
    codes = [F7_PRE * cfg * "Base.invokelatest(Main.zap_verify; pazur2=\"$out\", out=\"$vd\", part=\"$k/$n\", chans=false, save_rfi=true)"
             for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, F7_DIR * "logs", "verify"; every=120)
    Base.invokelatest(Main.zap_verify_summary; out=vd)
    SpaTs.compare_zap_methods(mine=out, outdir=F7_DIR * "zap_comparison_final/")
    println("Potem: final_v7_adopt()")
end

"Przyjęcie v7: pazur2 → pazur2_przed_v7 (albo _2, _3...), final_v7/pazur2_final → pazur2."
function final_v7_adopt()
    src, dst, bk = F7_DIR * "pazur2_final", F6_ROOT * "pazur2", F6_ROOT * "pazur2_przed_v7"
    isdir(src) || error("brak $src - najpierw final_v7_build()")
    k = 1
    while isdir(bk); k += 1; bk = F6_ROOT * "pazur2_przed_v7_$k"; end
    mv(dst, bk); mv(src, dst)
    println("Przyjęte v7: $(strip(read(dst * "/.wariant", String))). Stara wersja: $bk.")
end
