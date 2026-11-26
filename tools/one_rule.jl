# Test JEDNEJ metody wycinania nulli dla WSZYSTKICH pulsarow (bez mieszania: model w czesci pulsarow, automat
# w reszcie). Kazdy wariant to jedna cala regula, ta sama dla kazdego pulsara; RFI we wszystkich tak samo (jak teraz).
#   A  automat                    - nulle automatu (bez modelu energii)
#   B  model energii P >= 0.5     - regula z publikacji (TPA XIX) dla wszystkich, takze tam, gdzie model jest
#                                   niepewny (bez dopasowanego modelu: bez nulli)
#   C  impuls bez emisji          - energia impulsu < 3σ szumu (σ z energii poza impulsem)
#   D  seria bez emisji           - kazde 5 kolejnych impulsow, ktorych suma energii < 3σ (szum sumy 5)
# Porownanie z obecnym v8 (model w czesci pulsarow + automat w reszcie) i z metoda reczna; pazur2 bez zmian.
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/one_rule.jl"); one_rule()
# ~20-30 min na 32 procesach, wznawialne. Wynik: final_v7/one_rule/WYNIK.txt
isdefined(Main, :_s8_part) || include(joinpath(@__DIR__, "simple_v8.jl"))
S8_DIR = F7_DIR * "one_rule/"
# nazwa => regula nulli (symbol); ustawienia automatu bez zmian
OR_RULES = [
    ("v8 (obecny): model P>=0.5, bez modelu automat", :v8),
    ("A: automat", :auto),
    ("B: model energii P>=0.5 dla wszystkich", :model_all),
    ("C: impuls bez emisji (E < 3σ)", :pulse3),
    ("D: seria 5 impulsów bez emisji (suma < 3σ)", :run5),
]
S8_VARIANTS = [(r[1], NamedTuple(), :x, 0.0) for r in OR_RULES]     # tylko nazwy - do podsumowania

"Model energii pulsara niezaleznie od oceny wiarygodnosci (plik z P_null, jesli dopasowanie istnieje)."
function _or_model_any(psr)
    pf = NR_MODEL * "nulls_model/$psr.csv"
    isfile(pf) || return nothing
    rows = Int[]; P = Float64[]
    for l in Iterators.drop(eachline(pf), 1)
        r = split(l, ","); push!(rows, parse(Int, r[1]) + 1); push!(P, parse(Float64, r[3]))
    end
    return (rows=rows, P=P)
end

"Energie impulsow (on) i szum energii: σ z energii w ciaglym kawalku poza impulsem tej samej szerokosci (MAD)."
function _or_energy(D, bst, bend)
    PL = SpaTs.Plot
    N, nb = size(D)
    on = bst:bend
    offb = PL.dc_offblock(nb, bst, bend, length(on))
    rest = setdiff(SpaTs.Data.offpulse_bins(D, bst, bend), offb)
    isempty(rest) && (rest = offb)
    base = [mean(@view D[i, rest]) for i in 1:N]
    E = [sum(@view D[i, on]) - length(on) * base[i] for i in 1:N]
    Eo = [sum(@view D[i, offb]) - length(offb) * base[i] for i in 1:N]
    live = [any(!iszero, @view D[i, :]) for i in 1:N]
    u = Eo[live]
    σ = 1.4826 * median(abs.(u .- median(u)))
    return E, σ, live
end

function _or_nulls(rule, E, σ, live, m_ok, m_any)
    N = length(E)
    if rule == :pulse3
        return [i for i in 1:N if live[i] && E[i] < 3σ]
    elseif rule == :run5
        out = falses(N)
        for s in 1:N-4
            r = s:s+4
            all(live[r]) || continue
            sum(E[r]) < 3σ * sqrt(5) && (out[r] .= true)
        end
        return findall(out)
    elseif rule == :model_all
        return isnothing(m_any) ? Int[] : m_any.rows[m_any.P .>= 0.5]
    elseif rule == :model_ok
        return isnothing(m_ok) ? Int[] : m_ok.rows[m_ok.P .>= 0.5]     # opcja 3: tylko pewny model, reszta bez nulli
    elseif rule == :v8
        return isnothing(m_ok) ? nothing : m_ok.rows[m_ok.P .>= 0.5]
    end
    return nothing                                       # :auto
end

function _or_part(k, n)
    _VD = SpaTs.Data
    mkpath(S8_DIR)
    f = S8_DIR * "part_$k.csv"
    done = isfile(f) ? Set(split(l, ",")[1] for l in Iterators.drop(eachline(f), 1)) : Set{String}()
    isfile(f) || open(io -> println(io, S8_HEADER), f, "w")
    psrs = filter(!in(done), _f6_psrs()[k:n:end])
    for (j, psr) in enumerate(psrs)
        t0 = time()
        try
            p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
            D = _f6_load(psr, p)
            N = size(D, 1)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            on = bst:bend
            off = _VD.offpulse_bins(D, bst, bend)
            Ebase = [sum(@view D[i, on]) - length(on) * mean(@view D[i, off]) for i in 1:N]
            E, σ, live = _or_energy(D, bst, bend)
            zl = _VD.rfi_rms_z(D, bst, bend)
            zg = _VD.rfi_rms_z(D, bst, bend; local_win=0)
            m = _nr_model(psr); ma = _or_model_any(psr)
            P = fill(NaN, N); isnothing(m) || (P[m.rows] .= m.P)       # ocena wycietej emisji: tylko wiarygodny model
            snr_none = _VD.profile_snr(D, bst, bend)
            snr_man = something(get(get(p, "zap_cmp", Dict()), "snr_maciej", nothing), NaN)
            for rr in OR_RULES
                name, rule = rr[1], rr[2]; kw = length(rr) >= 3 ? rr[3] : NamedTuple()
                extra = _or_nulls(rule, E, σ, live, m, ma)
                a = isnothing(extra) ? _VD.auto_zap(D, bst, bend; null_mode=:auto, kw...) :
                                       _VD.auto_zap(D, bst, bend; null_mode=:model, null_extra=extra, kw...)
                sp = _VD.split_zaps(a.reason)
                rfi = sp.rfi; nul = sp.null
                cut = Set(a.bad); rs = Set(rfi)
                U = [i for i in 1:N if isfinite(P[i]) && !(i in rs)]
                C = [i for i in U if i in cut]; K = [i for i in U if !(i in cut)]
                ea, na = isempty(U) ? (NaN, NaN) : (sum(1 .- P[U]), sum(P[U]))
                Ee = isempty(U) || !(ea > 0) ? NaN : sum((1 .- P[U]) .* Ebase[U]) / ea
                row = [psr, name, N, length(a.bad), length(rfi), length(nul), snr_none,
                       _VD.profile_snr(D, bst, bend; drop=a.bad), snr_man, _VD.profile_snr(D, bst, bend; drop=rfi),
                       count(i -> !(i in cut) && zl[i] > 3, 1:N), count(i -> !(i in cut) && zg[i] > 3, 1:N),
                       isempty(C) ? 0.0 : sum(1 .- P[C]), ea, isempty(C) ? 0.0 : sum(P[C]), na,
                       isempty(K) ? NaN : mean(Ebase[K]), Ee]
                open(io -> println(io, join([x isa AbstractFloat ? @sprintf("%.6g", x) : string(x) for x in row], ",")), f, "a")
            end
            @printf("[%d/%d] %-12s ok (%.0f s)\n", j, length(psrs), psr, time() - t0)
        catch e
            e isa InterruptException && rethrow()
            println("[$j/$(length(psrs))] $psr BŁĄD: ", first(sprint(showerror, e), 200))
        end
        flush(stdout)
    end
end

function one_rule(; n=32)
    isdefined(Main, :_run_parallel) || Base.include(Main, "modules/zap_verify.jl")
    mkpath(S8_DIR)
    codes = [F7_PRE * "include(\"tools/one_rule.jl\"); Base.invokelatest(_or_part, $k, $n)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, S8_DIR * "logs", "or"; every=120)
    one_rule_summary()
end

one_rule_summary() = simple_v8_summary(; title="Jedna metoda wycinania nulli dla wszystkich pulsarów (porównanie z v8)", basename_="v8")

# ---- Margines przy RFI na prawdziwych danych: regula B z marginesem (1 impuls z kazdej strony odcinka RFI) i bez.
# Ocena: S/N profilu (bez marginesu w RFI zostaja sasiednie, czesciowo zaklocone impulsy -> S/N spada) i liczba
# niewycietych impulsow z szumem poza impulsem > 3σ ("RFI zostawione").
#   include("tools/one_rule.jl"); one_rule_margin()         # ~20-30 min; wynik: final_v7/one_rule_margin/WYNIK.txt
function _or_margin_setup!()
    global S8_DIR = F7_DIR * "one_rule_margin/"
    global OR_RULES = [("B z marginesem przy RFI", :model_all, NamedTuple()), ("B bez marginesu przy RFI", :model_all, (pad=0,))]
    global S8_VARIANTS = [(r[1], NamedTuple(), :x, 0.0) for r in OR_RULES]
end

function one_rule_margin(; n=32)
    isdefined(Main, :_run_parallel) || Base.include(Main, "modules/zap_verify.jl")
    _or_margin_setup!(); mkpath(S8_DIR)
    codes = [F7_PRE * "include(\"tools/one_rule.jl\"); _or_margin_setup!(); Base.invokelatest(_or_part, $k, $n)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, S8_DIR * "logs", "orm"; every=120)
    one_rule_margin_summary()
end

one_rule_margin_summary() = (_or_margin_setup!(); simple_v8_summary(; title="Reguła B: margines przy RFI czy bez (prawdziwe dane)", basename_="B z marginesem"))

# ---- Wersja v9 = regula B dla wszystkich pulsarow (po tescie i decyzji). pazur2 bez zmian do one_rule_adopt().
#   include("tools/one_rule.jl"); one_rule_build()                # margines przy RFI jak dotad (margin=false: bez)
#   one_rule_adopt()                                              # pazur2 -> pazur2_przed_v9, nowa wersja -> pazur2
OR_OUT = F7_DIR * "one_rule/"
function one_rule_build(; rule=:model_ok, margin=true, n=32)
    # rule = :model_ok (opcja 3, decyzja 2026-11-26): nulle z modelu tylko w pulsarach z wiarygodnym modelem, reszta
    # bez nulli; :model_all (B): model kazdego pulsara z dopasowaniem
    isdefined(Main, :_run_parallel) || Base.include(Main, "modules/zap_verify.jl")
    out = OR_OUT * "pazur2_final/"; mkpath(out)
    kw = margin ? NamedTuple() : (pad=0,)
    anym = rule == :model_all
    write(out * ".wariant", "v9: nulle = model energii P>=0.5 " * (anym ? "dla wszystkich pulsarow z modelem" : "tylko gdy model wiarygodny") *
                            " (reszta bez nulli); RFI jak v8" * (margin ? ", margines przy RFI" : ", bez marginesu przy RFI"))
    println("\n===== WERSJA v9 (" * (anym ? "reguła B" : "opcja 3: model tylko pewny") * ") → $out =====")
    cfg = "SpaTs.Data.ZAP_STEPS[] = $(repr(kw)); SpaTs.Data.NULL_MODE[] = :model; SpaTs.Data.NULL_THR[] = 0.5; " *
          "SpaTs.Data.NULL_MODEL_ANY[] = $anym; SpaTs.Data.NULL_MODEL_FALLBACK[] = :none; " *
          "SpaTs.Data.NULL_MODEL_DIR[] = \"$(F7_DIR)model/\"; "
    cf = "SpaTs.chan_first_list([\"$(F6_ROOT)zap_verify_v3/\", \"$(F6_ROOT)zap_verify/\"])"
    codes = [F7_PRE * cfg * "SpaTs.autozap_all(out=\"$out\", only=$(repr(_f6_psrs()[k:n:end])), skip_done=true, keep_manual=false, " *
             "reuse_window=true, window_root=\"$(F6_ROOT)pazur2/\", chan_first=$cf)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, OR_OUT * "logs", "build"; every=120)
    vd = OR_OUT * "zap_verify_final/"
    codes = [F7_PRE * cfg * "Base.invokelatest(Main.zap_verify; pazur2=\"$out\", out=\"$vd\", part=\"$k/$n\", chans=false, save_rfi=true)"
             for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, OR_OUT * "logs", "verify"; every=120)
    Base.invokelatest(Main.zap_verify_summary; out=vd)
    SpaTs.compare_zap_methods(mine=out, outdir=OR_OUT * "zap_comparison_final/")
    println("Potem, po decyzji: one_rule_adopt()")
end

"Przyjecie v9: pazur2 -> pazur2_przed_v9 (albo _2...), one_rule/pazur2_final -> pazur2."
function one_rule_adopt()
    src, dst, bk = OR_OUT * "pazur2_final", F6_ROOT * "pazur2", F6_ROOT * "pazur2_przed_v9"
    isdir(src) || error("brak $src - najpierw one_rule_build()")
    k = 1
    while isdir(bk); k += 1; bk = F6_ROOT * "pazur2_przed_v9_$k"; end
    mv(dst, bk); mv(src, dst)
    println("Przyjęte: $(strip(read(dst * "/.wariant", String))). Stara wersja: $bk.")
    println("Domyślne w kodzie są już v9 (Data.NULL_MODEL_FALLBACK[] = :none, commit ee57001).")
end

# ---- Szczegoly: B wobec v8 dla kazdego pulsara (S/N, nulle, status modelu energii).
#   include("tools/one_rule.jl"); one_rule_details()      # wynik: final_v7/one_rule/SZCZEGOLY.txt
function _or_rows(dir)
    R = Dict{String,Dict{String,Dict{String,Float64}}}()
    for f in filter(f -> occursin(r"^part_\d+\.csv$", basename(f)), readdir(dir; join=true))
        H = split(first(eachline(f)), ",")
        for l in Iterators.drop(eachline(f), 1)
            v = split(l, ",")
            length(v) > length(H) && (v = vcat(v[1:1], [join(v[2:2+length(v)-length(H)], ",")], v[3+length(v)-length(H):end]))
            length(v) == length(H) || continue
            get!(R, v[2], Dict{String,Dict{String,Float64}}())[v[1]] =
                Dict(h => something(tryparse(Float64, x), NaN) for (h, x) in zip(H[3:end], v[3:end]))
        end
    end
    R
end

function one_rule_details(; a="v8 (obecny): model P>=0.5, bez modelu automat", b="B: model energii P>=0.5 dla wszystkich", tie=1.0)
    R = _or_rows(F7_DIR * "one_rule/")
    A = R[a]; B = R[b]
    st = Dict{String,String}()
    for l in Iterators.drop(eachline(NR_MODEL * "nulls_model.csv"), 1)
        r = split(l, ","); length(r) >= 2 && (st[r[1]] = r[2])
    end
    status(q) = (s = get(st, q, "brak w liście"); s == "ok" ? "wiarygodny" : isfile(NR_MODEL * "nulls_model/$q.csv") ? "niepewny ($s)" : "brak modelu ($s)")
    L = String[]; pr(s) = (println(s); push!(L, s))
    ps = sort([q for q in keys(B) if haskey(A, q)])
    pr("# Szczegóły: $b wobec $a ($(length(ps)) pulsarów)\n")
    # liczby wg statusu modelu
    pr("## Status modelu energii\n")
    pr("| status | pulsary | B gorsze o >$(tie)% | B lepsze o >$(tie)% | nulle v8 | nulle B |")
    pr("|---|---|---|---|---|---|")
    groups = Dict{String,Vector{String}}()
    for q in ps; push!(get!(groups, first(split(status(q), " ")), String[]), q); end
    for (g, qs) in sort(collect(groups), by=x -> x[1])
        d = [100 * (B[q]["snr"] / A[q]["snr"] - 1) for q in qs if A[q]["snr"] > 0]
        pr(@sprintf("| %s | %d | %d | %d | %d | %d |", g, length(qs), count(<(-tie), d), count(>(tie), d),
                    round(Int, sum(A[q]["n_null"] for q in qs)), round(Int, sum(B[q]["n_null"] for q in qs))))
    end
    for (title, sel) in (("B gorsze od v8", x -> x < -tie), ("B lepsze od v8", x -> x > tie))
        pr("\n## $title (S/N o >$(tie)%)\n")
        pr("| pulsar | model | S/N v8 | S/N B | zmiana | S/N ręcznie | nulle v8 | nulle B | wycięte v8 | wycięte B | N |")
        pr("|---|---|---|---|---|---|---|---|---|---|---|")
        rows = [(q, 100 * (B[q]["snr"] / A[q]["snr"] - 1)) for q in ps if A[q]["snr"] > 0]
        for (q, d) in sort(filter(x -> sel(x[2]), rows), by=x -> sel(-1) ? x[2] : -x[2])
            pr(@sprintf("| %s | %s | %.1f | %.1f | %+.1f%% | %.1f | %d | %d | %d | %d | %d |", q, status(q), A[q]["snr"], B[q]["snr"], d,
                        B[q]["snr_manual"], A[q]["n_null"], B[q]["n_null"], A[q]["n_cut"], B[q]["n_cut"], B[q]["N"]))
        end
    end
    write(F7_DIR * "one_rule/SZCZEGOLY.txt", join(L, "\n") * "\n")
    println("\nZapisane: $(F7_DIR)one_rule/SZCZEGOLY.txt")
end

# ---- Trzy opcje do wyboru na tych samych danych: v8, B (model wszedzie), 3 (model tylko pewny, reszta bez nulli)
#   include("tools/one_rule.jl"); one_rule_three()          # ~20-30 min; wynik: final_v7/one_rule_three/WYNIK.txt
function _or_three_setup!()
    global S8_DIR = F7_DIR * "one_rule_three/"
    global OR_RULES = [("v8 (obecny) model pewny + automat", :v8), ("B model wszedzie", :model_all),
                       ("3 model tylko pewny reszta bez nulli", :model_ok)]
    global S8_VARIANTS = [(r[1], NamedTuple(), :x, 0.0) for r in OR_RULES]
end

function one_rule_three(; n=32)
    isdefined(Main, :_run_parallel) || Base.include(Main, "modules/zap_verify.jl")
    _or_three_setup!(); mkpath(S8_DIR)
    codes = [F7_PRE * "include(\"tools/one_rule.jl\"); _or_three_setup!(); Base.invokelatest(_or_part, $k, $n)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, S8_DIR * "logs", "or3"; every=120)
    one_rule_three_summary()
end

one_rule_three_summary() = (_or_three_setup!(); simple_v8_summary(; title="Trzy opcje wycinania nulli: v8, B, 3", basename_="v8"))

"Po przyjeciu v9 (one_rule_adopt): wykresy koncowe sliding LRFS i tabela dryfu na nowych danych."
function one_rule_replot()
    Base.include(Main, joinpath(@__DIR__, "sliding_lrfs.jl")); Base.include(Main, joinpath(@__DIR__, "drift_table.jl"))
    Base.invokelatest(Main.sliding_lrfs_final)
    Base.invokelatest(Main.drift_table_run)
end
