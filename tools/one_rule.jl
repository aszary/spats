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
isdefined(Main, :_tz_tpa) || include(joinpath(@__DIR__, "tpa_zap_compare.jl"))     # lista wyciec TPA (Xiaoxi) do oceny
OR_HEADER = S8_HEADER * ",rfi_clean,left_bad,tpa_n,tpa_hit,rfi_clean_odc"
isdefined(Main, :_OR_TPA) || (global _OR_TPA = Ref{Any}(nothing))
_or_tpa() = (isnothing(_OR_TPA[]) && (_OR_TPA[] = try _tz_tpa() catch; Dict{String,Any}() end); _OR_TPA[])
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
    isfile(f) || open(io -> println(io, OR_HEADER), f, "w")
    psrs = filter(!in(done), _f6_psrs()[k:n:end])
    T = _or_tpa()
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
            F = _VD.snr_fast(D, bst, bend)                # szybkie S/N: profil i biny poza impulsem raz (te same liczby)
            snr_none = _VD.profile_snr(F)
            snr_man = something(get(get(p, "zap_cmp", Dict()), "snr_maciej", nothing), NaN)
            # jakosc wyciec: ksztalt tla, czyste / wyraznie zle tlo, lista TPA (Xiaoxi)
            ks = _VD.rfi_struct_ratio(D, bst, bend)
            clean(i) = zl[i] < 2 && zg[i] < 2 && ks[i] < 1.5
            badbg(i) = live[i] && (ks[i] > 3 || zl[i] > 10 || zg[i] > 10)
            tz = haskey(T, psr) ? Set(i + 1 for i in T[psr].zaps if 0 <= i < N) : Set{Int}()
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
                       _VD.profile_snr(F; drop=a.bad), snr_man, _VD.profile_snr(F; drop=rfi),
                       count(i -> !(i in cut) && zl[i] > 3, 1:N), count(i -> !(i in cut) && zg[i] > 3, 1:N),
                       isempty(C) ? 0.0 : sum(1 .- P[C]), ea, isempty(C) ? 0.0 : sum(P[C]), na,
                       isempty(K) ? NaN : mean(Ebase[K]), Ee,
                       count(clean, rfi), count(i -> badbg(i) && !(i in cut), 1:N), length(tz), length(intersect(tz, cut)),
                       count(i -> clean(i) && a.reason[i] in (4, 5), 1:N)]
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
                            " (reszta bez nulli); RFI jak v8" * (margin ? ", margines przy RFI" : ", bez marginesu przy RFI") *
                            "; silne RFI (z > $(SpaTs.Data.RFI_HARD[])) zawsze wycinane")
    println("\n===== WERSJA v9 (" * (anym ? "reguła B" : "opcja 3: model tylko pewny") * ") → $out =====")
    cfg = "SpaTs.Data.ZAP_STEPS[] = $(repr(kw)); SpaTs.Data.NULL_MODE[] = :model; SpaTs.Data.NULL_THR[] = 0.5; " *
          "SpaTs.Data.NULL_MODEL_ANY[] = $anym; SpaTs.Data.NULL_MODEL_FALLBACK[] = :none; " *
          "SpaTs.Data.NULL_MODEL_DIR[] = \"$(F7_DIR)model/\"; "
    cf = "SpaTs.chan_first_list([\"$(F6_ROOT)zap_verify_v3/\", \"$(F6_ROOT)zap_verify/\"])"
    codes = [F7_PRE * cfg * "SpaTs.autozap_all(out=\"$out\", only=$(repr(_f6_psrs()[k:n:end])), skip_done=true, keep_manual=false, " *
             "reuse_window=true, window_root=\"$(F6_ROOT)pazur2/\", chan_first=$cf)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, OR_OUT * "logs", "build"; every=120)
    vd = OR_OUT * "zap_verify_final/"
    # nowa budowa: weryfikacja od zera (jest wznawialna - inaczej pominelaby pulsary z poprzedniej budowy)
    isfile(out * ".verify_started") || rm(vd; recursive=true, force=true); touch(out * ".verify_started")
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

# ---- Silne RFI zawsze wycinane (z > prog, bez bezpiecznika S/N), na v9; test na wszystkich pulsarach.
#   include("tools/one_rule.jl"); one_rule_hard()            # ~20-30 min; wynik: final_v7/one_rule_hard7/WYNIK.txt
# (one_rule_hard/ = pierwszy test bez wykluczenia jasnych impulsow: z>6 psul 162 pulsary)
function _or_hard_setup!()
    global S8_DIR = F7_DIR * "one_rule_hard7/"
    global OR_RULES = [("v9", :model_ok, (rfi_struct=Inf,)),
                       ("v9 + RFI: struktura tła >2x i z>6", :model_ok, (rfi_struct=2.0, rfi_struct_z=6.0)),
                       ("v9 + RFI: struktura tła >2x i z>4", :model_ok, (rfi_struct=2.0, rfi_struct_z=4.0)),
                       ("v9 + RFI jak TPA VIII: rms > mediana + 4σ, iteracyjnie", :model_ok, (rfi_tpa_k=4.0,))]
    global S8_VARIANTS = [(r[1], NamedTuple(), :x, 0.0) for r in OR_RULES]
end

function one_rule_hard(; n=32)
    isdefined(Main, :_run_parallel) || Base.include(Main, "modules/zap_verify.jl")
    _or_hard_setup!(); mkpath(S8_DIR)
    codes = [F7_PRE * "include(\"tools/one_rule.jl\"); _or_hard_setup!(); Base.invokelatest(_or_part, $k, $n)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, S8_DIR * "logs", "orh"; every=120)
    one_rule_hard_summary()
end

one_rule_hard_summary() = (_or_hard_setup!(); simple_v8_summary(; title="Silne RFI zawsze wycinane (bez bezpiecznika S/N) - porównanie z v9", basename_="v9"))

# ---- Podglad: ktore impulsy dochodza przy silnym RFI (rfi_hard) wzgledem v9 - stos impulsow z zaznaczeniami,
# z(RFI) i jasnosc tych impulsow, profil sredni dodanych impulsow (czy maja emisje / paski).
#   include("tools/one_rule.jl"); hard_rfi_view(["J1012-5857", "J0820-1350", "J1110-5637", "J0932-3217"])
#   -> final_v7/one_rule_hard5/podglad.pdf
using PyPlot, PyCall
function hard_rfi_view(psrs; struct_=3.0, z=6.0, excess=0.0, zkw=NamedTuple(), out=F7_DIR * "one_rule_hard5/podglad.pdf")
    _VD = SpaTs.Data
    PyPlot = Main.PyPlot
    mkpath(dirname(out)); PdfPages = Main.PyCall.pyimport("matplotlib.backends.backend_pdf").PdfPages
    pdf = PdfPages(out)
    for psr in psrs
        try
            p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
            D = Float64.(_f6_load(psr, p)); N, nb = size(D)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            m = _nr_model(psr); extra = isnothing(m) ? Int[] : m.rows[m.P .>= 0.5]
            a0 = _VD.auto_zap(D, bst, bend; null_mode=:model, null_extra=extra, rfi_struct=Inf)
            a1 = _VD.auto_zap(D, bst, bend; null_mode=:model, null_extra=extra, rfi_struct=struct_, zkw...)
            add = sort(setdiff(a1.bad, a0.bad)); zz = _VD.rfi_rms_z(D, bst, bend)
            sn = _VD.pulse_signal_noise(D, bst, bend)
            sg = sn.s[sn.good]; md = median(sg); sd = 1.4826 * median(abs.(sg .- md))
            s0 = _VD.profile_snr(D, bst, bend; drop=a0.bad); s1 = _VD.profile_snr(D, bst, bend; drop=a1.bad)
            fig = PyPlot.figure(figsize=(11, 7.5))
            fig.suptitle(@sprintf("%s: przy strukturze >%.1fx dochodzi %d impulsów; S/N v9 %.1f -> %.1f (%+.1f%%)", psr, struct_, length(add), s0, s1, 100(s1 / s0 - 1)), fontsize=9)
            live = [any(!iszero, @view D[i, :]) for i in 1:N]; v = vec(D[live, :])
            ax = fig.add_axes([0.06, 0.32, 0.50, 0.6])
            ax.imshow(D, origin="lower", aspect="auto", cmap="Greys_r", vmin=quantile(v, 0.01), vmax=quantile(v, 0.995),
                      interpolation="nearest", extent=[-0.5, nb - 0.5, 0.5, N + 0.5])
            # znaczniki: niebieskie - wyciete juz w v9 (kazdy powod), czerwone - dochodza z nowej reguly
            isempty(a0.bad) || ax.scatter(fill(-0.045nb, length(a0.bad)), a0.bad, marker="_", s=10, c="royalblue", clip_on=false)
            isempty(add) || ax.scatter(fill(-0.02nb, length(add)), add, marker=">", s=12, c="red", clip_on=false)
            ax.axvline(bst, c="orange", lw=0.5); ax.axvline(bend, c="orange", lw=0.5)
            ax.set_ylabel("impuls (od 1); niebieskie: wycięte w v9, czerwone: dochodzą", fontsize=7); ax.tick_params(labelsize=6)
            # struktura tla poza oknem w kazdym impulsie (to mierzy regula) - widac cale zle odcinki
            ksv = _VD.rfi_struct_ratio(D, bst, bend)
            ax2 = fig.add_axes([0.62, 0.62, 0.34, 0.3])
            ax2.scatter(1:N, ksv, s=2, c="grey", lw=0)
            isempty(a0.bad) || ax2.scatter(a0.bad, ksv[a0.bad], s=3, c="royalblue", lw=0)
            isempty(add) || ax2.scatter(add, ksv[add], s=8, c="red", lw=0)
            ax2.axhline(struct_, c="k", ls="--", lw=0.6); ax2.set_yscale("log")
            ax2.set_xlabel("impuls", fontsize=7); ax2.set_ylabel("struktura poza oknem / typowa", fontsize=7); ax2.tick_params(labelsize=6)
            ax3 = fig.add_axes([0.62, 0.32, 0.34, 0.22])
            good = [i for i in 1:N if live[i] && !(i in a1.bad)]
            pg = vec(mean(D[good, :], dims=1)); ax3.plot(pg ./ max(maximum(pg), eps()), c="grey", lw=0.8, label="zostawione")
            if !isempty(add)
                pa = vec(mean(D[add, :], dims=1)); ax3.plot(pa ./ max(maximum(pg), eps()), c="red", lw=0.8, label="dochodzące (średnia)")
            end
            ax3.legend(fontsize=6); ax3.tick_params(labelsize=6); ax3.set_xlabel("bin", fontsize=7)
            ax4 = fig.add_axes([0.06, 0.06, 0.50, 0.18])
            for (k, i) in enumerate(add[1:min(end, 6)])
                ax4.plot(D[i, :] .+ 3k * std(v), lw=0.5)
            end
            ax4.set_title("kilka dochodzących impulsów (przesunięte w pionie)", fontsize=7); ax4.tick_params(labelsize=6)
            pdf.savefig(fig); PyPlot.close(fig)
            @printf("%-12s dochodzi %d imp., S/N %+.1f%%, jasność (σ) mediana %.1f, z mediana %.1f\n", psr, length(add), 100(s1 / s0 - 1),
                    isempty(add) ? NaN : median((sn.s[add] .- md) ./ sd), isempty(add) ? NaN : median(zz[add]))
            # rozbicie: wprost z reguly silnego RFI (te same warunki co w auto_zap) a posrednio (inny wybor reszty)
            zg = _VD.rfi_rms_z(D, bst, bend; local_win=0); q = _VD.rfi_rms_ratio(D, bst, bend)
            ks = _VD.rfi_struct_ratio(D, bst, bend)
            direct = [i for i in add if ks[i] > struct_]
            other = setdiff(add, direct)
            rs = Dict{Int,Int}(); for i in add; rs[a1.reason[i]] = get(rs, a1.reason[i], 0) + 1; end
            @printf("    wprost z reguły: %d (szum / typowy: mediana %.2f, min %.2f) | pośrednio: %d (szum / typowy: mediana %.2f, z lok. %.1f, z całość %.1f) | powody w a1: %s\n",
                    length(direct), isempty(direct) ? NaN : median(q[direct]), isempty(direct) ? NaN : minimum(q[direct]),
                    length(other), isempty(other) ? NaN : median(q[other]), isempty(other) ? NaN : median(zz[other]),
                    isempty(other) ? NaN : median(zg[other]), join(["$k:$v" for (k, v) in sort(collect(rs))], " "))
            @printf("    szum / typowy we wszystkich impulsach: mediana 1.00, 90%%: %.2f, 99%%: %.2f\n",
                    quantile(filter(isfinite, q), 0.9), quantile(filter(isfinite, q), 0.99))
            fm(v) = (w = filter(isfinite, v); isempty(w) ? "-" : @sprintf("%.2f", median(w)))
            @printf("    struktura poza oknem / typowa: wprost %s, pośrednio %s; wszystkie impulsy 90%%: %.2f, 99%%: %.2f\n",
                    fm(ks[direct]), fm(ks[other]), quantile(filter(isfinite, ks), 0.9), quantile(filter(isfinite, ks), 0.99))
            # okno: czy jest tam sygnal pulsara (korelacja z profilem) i ile (energia wzgledem typowej)
            wc = _VD.window_signal_check(D, bst, bend)
            @printf("    okno - korelacja z profilem: wprost %s, pośrednio %s, wszystkie %s (10%% najniższych < %.2f); energia / typowa: wprost %s, pośrednio %s\n",
                    fm(wc.r[direct]), fm(wc.r[other]), fm(wc.r), quantile(filter(isfinite, wc.r), 0.1), fm(wc.e[direct]), fm(wc.e[other]))
        catch e
            e isa InterruptException && rethrow()
            println("$psr BŁĄD: ", first(sprint(showerror, e), 200))
        end
    end
    pdf.close(); println("Podgląd: $out")
end

# ---- Tabela impulsow: struktura poza oknem, z, szum, korelacja okna z profilem i czy wyciete (v9 / nowa regula).
#   include("tools/one_rule.jl"); rfi_pulse_table("J0932-3217", 1:25)
function rfi_pulse_table(psr, rng; struct_=3.0)
    _VD = SpaTs.Data
    p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
    D = Float64.(_f6_load(psr, p)); N = size(D, 1)
    bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
    m = _nr_model(psr); extra = isnothing(m) ? Int[] : m.rows[m.P .>= 0.5]
    a0 = _VD.auto_zap(D, bst, bend; null_mode=:model, null_extra=extra, rfi_struct=Inf)
    a1 = _VD.auto_zap(D, bst, bend; null_mode=:model, null_extra=extra, rfi_struct=struct_)
    ks = _VD.rfi_struct_ratio(D, bst, bend); zl = _VD.rfi_rms_z(D, bst, bend); q = _VD.rfi_rms_ratio(D, bst, bend)
    wc = _VD.window_signal_check(D, bst, bend)
    powod = Dict(0 => "-", 1 => "pusty", 2 => "RFI/śmieć", 3 => "null", 4 => "zaszumiony", 5 => "margines", 6 => "wygaszony")
    println("$psr: impulsy od 0 (jak w pazi); próg struktury $(struct_)x")
    println(" impuls | struktura/typowa | z (RFI) | szum/typowy | korelacja okna | v9 | nowa reguła")
    for i in rng
        1 <= i + 1 <= N || continue
        j = i + 1
        @printf(" %6d | %16.2f | %7.1f | %11.2f | %14.2f | %-10s | %s\n", i, ks[j], zl[j], q[j], wc.r[j],
                get(powod, Int(a0.reason[j]), string(a0.reason[j])), get(powod, Int(a1.reason[j]), string(a1.reason[j])))
    end
end

# ---- Szybki test kilku progow na kilku pulsarach (minuty): ile dochodzi, S/N, i co ZLEGO zostaje niewyciete.
# "Zle zostawione" = niewyciety impuls z wyrazna struktura poza oknem (> 2x typowej) i silnym szumem (z > 6 wzgledem
# sasiadow albo calosci) - takie impulsy powinny wypasc. "Wyciete z sygnalem" = dochodzace impulsy z typowa korelacja
# okna z profilem i bez struktury > 2x (podejrzenie, ze to dobre impulsy).
#   include("tools/one_rule.jl"); quick_rfi_test()        # -> final_v7/one_rule_quick/{WYNIK.txt, podglad_<prog>.pdf}
function quick_rfi_test(psrs=["J1012-5857", "J0820-1350", "J1110-5637", "J0932-3217", "J1915+0752"];
                        variants=[("2x", (rfi_struct=2.0,)), ("2.5x", (rfi_struct=2.5,)), ("1.5x i z>6", (rfi_struct=1.5, rfi_struct_z=6.0)),
                                  ("2x i z>6", (rfi_struct=2.0, rfi_struct_z=6.0)), ("2x i z>4", (rfi_struct=2.0, rfi_struct_z=4.0)),
                                  ("TPA 4σ", (rfi_tpa_k=4.0,))],
                        out=F7_DIR * "one_rule_quick/", plots=true)
    _VD = SpaTs.Data
    mkpath(out)
    L = String[]; pr(s) = (println(s); push!(L, s))
    pr("# Szybki test reguły RFI (struktura tła poza oknem > próg × typowa): $(join(psrs, ", "))\n")
    pr("| pulsar | próg | dochodzi | odpada (v9 wycinał) | S/N wobec v9 | złe zostawione (struktura>2x i z>6) | wycięte, a wyglądają dobrze |")
    pr("|---|---|---|---|---|---|---|")
    for psr in psrs
        try
            p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
            D = Float64.(_f6_load(psr, p)); N = size(D, 1)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            m = _nr_model(psr); extra = isnothing(m) ? Int[] : m.rows[m.P .>= 0.5]
            ks = _VD.rfi_struct_ratio(D, bst, bend)
            zl = _VD.rfi_rms_z(D, bst, bend); zg = _VD.rfi_rms_z(D, bst, bend; local_win=0)
            wc = _VD.window_signal_check(D, bst, bend)
            rtyp = median(filter(isfinite, wc.r))
            live = [any(!iszero, @view D[i, :]) for i in 1:N]
            bad(i) = live[i] && ks[i] > 2 && (zl[i] > 6 || zg[i] > 6)
            F = _VD.snr_fast(D, bst, bend)
            a0 = _VD.auto_zap(D, bst, bend; null_mode=:model, null_extra=extra, rfi_struct=Inf)
            s0 = _VD.profile_snr(F; drop=a0.bad)
            left0 = [i - 1 for i in 1:N if bad(i) && !(i in a0.bad)]
            pr(@sprintf("| %s | v9 | - | - | - | %d %s | - |", psr, length(left0), isempty(left0) ? "" : string(left0[1:min(end, 12)])))
            for (lab, kw) in variants
                a1 = _VD.auto_zap(D, bst, bend; null_mode=:model, null_extra=extra, kw...)
                cut = Set(a1.bad); add = sort(setdiff(a1.bad, a0.bad))
                s1 = _VD.profile_snr(F; drop=a1.bad)
                left = [i - 1 for i in 1:N if bad(i) && !(i in cut)]
                good = [i - 1 for i in add if ks[i] <= 2 && isfinite(wc.r[i]) && wc.r[i] >= rtyp - 0.1]
                pr(@sprintf("| %s | %s | %d | %d | %+.1f%% | %d %s | %d %s |", psr, lab, length(add), length(setdiff(a0.bad, a1.bad)), 100(s1 / s0 - 1),
                            length(left), isempty(left) ? "" : string(left[1:min(end, 12)]),
                            length(good), isempty(good) ? "" : string(good[1:min(end, 12)])))
            end
        catch e
            e isa InterruptException && rethrow()
            pr("| $psr | BŁĄD: $(first(sprint(showerror, e), 150)) |")
        end
    end
    pr("\nNumery impulsów od 0 (jak w pazi). Struktura = tło poza oknem wygładzone, względem typowego impulsu.")
    write(out * "WYNIK.txt", join(L, "\n") * "\n")
    if plots
        for (k, (lab, kw)) in enumerate(variants)
            hard_rfi_view(psrs; struct_=get(kw, :rfi_struct, Inf), zkw=Base.structdiff(kw, (rfi_struct=0,)), out=out * "podglad_$(k)_" * replace(lab, r"[^0-9a-z.]" => "") * ".pdf")
        end
    end
    println("\nZapisane: $(out)WYNIK.txt" * (plots ? " i podglad_<próg>.pdf" : ""))
end

# ---- CALE RFI jedna prosta regula (jak nulle: jedna zasada). Wylaczone wszystkie kroki RFI automatu (pojedyncze smieci,
# skoki szumu, zaszumione odcinki / zly poczatek, regula z bezpiecznikiem S/N) - zostaje tylko wybrana regula.
# Nulle jak v9 (model, tylko wiarygodny). Kanaly wycinane przed sumowaniem bez zmian (to przygotowanie danych).
OR_SIMPLE = (use_junk=false, use_jump=false, use_noisy=false, rfi_safe=false)
OR_SIMPLE_VARIANTS = [
    ("v9 (obecne: kilka krokow RFI)", NamedTuple()),
    ("tylko TPA VIII: rms > mediana + 4σ, iteracyjnie", merge(OR_SIMPLE, (rfi_tpa_k=4.0,))),
    ("tylko TPA VIII z 6σ", merge(OR_SIMPLE, (rfi_tpa_k=6.0,))),
    ("tylko z > 6 (sasiedzi albo calosc)", merge(OR_SIMPLE, (rfi_struct=0.0, rfi_struct_z=6.0))),
    ("tylko struktura tla > 2x i z > 6", merge(OR_SIMPLE, (rfi_struct=2.0, rfi_struct_z=6.0))),
]
function _or_simple_setup!()
    global S8_DIR = F7_DIR * "one_rule_simple/"
    global OR_RULES = [(n, :model_ok, kw) for (n, kw) in OR_SIMPLE_VARIANTS]
    global S8_VARIANTS = [(r[1], NamedTuple(), :x, 0.0) for r in OR_RULES]
end
"Test na wszystkich pulsarach: kazda prosta regula RFI zamiast wszystkich krokow RFI automatu. -> final_v7/one_rule_simple/WYNIK.txt"
function one_rule_simple(; n=32)
    isdefined(Main, :_run_parallel) || Base.include(Main, "modules/zap_verify.jl")
    _or_simple_setup!(); mkpath(S8_DIR)
    codes = [F7_PRE * "include(\"tools/one_rule.jl\"); _or_simple_setup!(); Base.invokelatest(_or_part, $k, $n)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, S8_DIR * "logs", "ors"; every=120)
    one_rule_simple_summary()
end
one_rule_simple_summary() = (_or_simple_setup!(); simple_v8_summary(; title="Całe RFI jedną prostą regułą (nulle jak v9) - porównanie z v9", basename_="v9"))
"Szybko na 5 pulsarach (tabela + podglady): proste reguly RFI."
quick_simple() = quick_rfi_test(; variants=OR_SIMPLE_VARIANTS[2:end], out=F7_DIR * "one_rule_quick_simple/")

# ---- Prosta zasada RFI z DWOCH czesci: (1) zle odcinki (krok automatu "zaszumione odcinki / zly poczatek" -
# odpowiednik obcinania zlych fragmentow u TPA VIII) + (2) jedna regula dla pojedynczych impulsow. Wylaczone: pojedyncze
# smieci, skoki szumu, regula z bezpiecznikiem S/N. Sama regula pojedynczych impulsow (one_rule_simple) gubila dlugie
# zle odcinki (J1844 -101%, sredni zysk S/N +5% zamiast +19%).
OR_SIMPLE2 = (use_junk=false, use_jump=false, use_noisy=true, rfi_safe=false)
OR_SIMPLE2_VARIANTS = [
    ("v9 (obecne: kilka krokow RFI)", NamedTuple()),
    ("odcinki + TPA VIII 4σ", merge(OR_SIMPLE2, (rfi_tpa_k=4.0,))),
    ("odcinki + TPA VIII 6σ", merge(OR_SIMPLE2, (rfi_tpa_k=6.0,))),
    ("odcinki + z > 6", merge(OR_SIMPLE2, (rfi_struct=0.0, rfi_struct_z=6.0))),
    ("odcinki + struktura > 2x i z > 6", merge(OR_SIMPLE2, (rfi_struct=2.0, rfi_struct_z=6.0))),
]
function _or_simple2_setup!()
    global S8_DIR = F7_DIR * "one_rule_simple2/"
    global OR_RULES = [(n, :model_ok, kw) for (n, kw) in OR_SIMPLE2_VARIANTS]
    global S8_VARIANTS = [(r[1], NamedTuple(), :x, 0.0) for r in OR_RULES]
end
"Test na wszystkich pulsarach: zle odcinki + jedna regula dla pojedynczych impulsow. -> final_v7/one_rule_simple2/WYNIK.txt"
function one_rule_simple2(; n=32)
    isdefined(Main, :_run_parallel) || Base.include(Main, "modules/zap_verify.jl")
    _or_simple2_setup!(); mkpath(S8_DIR)
    codes = [F7_PRE * "include(\"tools/one_rule.jl\"); _or_simple2_setup!(); Base.invokelatest(_or_part, $k, $n)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, S8_DIR * "logs", "ors2"; every=120)
    one_rule_simple2_summary()
end
one_rule_simple2_summary() = (_or_simple2_setup!(); simple_v8_summary(; title="Prosta zasada RFI: złe odcinki + jedna reguła dla pojedynczych impulsów - porównanie z v9", basename_="v9"))
quick_simple2() = quick_rfi_test(["J1012-5857", "J0820-1350", "J1110-5637", "J0932-3217", "J1915+0752", "J1844-0433", "J1842+0257"];
                                 variants=OR_SIMPLE2_VARIANTS[2:end], out=F7_DIR * "one_rule_quick_simple2/")

# ---- Proste zasady RFI BEZ marginesow, oceniane po tym, CO wycinaja (nie po S/N):
#   wyciete z czystym tlem (ksztalt < 1.5x typowego i szum w normie, z < 2) - prawdopodobnie niepotrzebnie;
#   zostawione z wyraznie zlym tlem (ksztalt > 3x albo z > 10) - niewyciete RFI;
#   zgodnosc z lista wyciec TPA (Xiaoxi).
#   include("tools/one_rule.jl"); one_rule_eval()           # -> final_v7/one_rule_eval/{WYNIK.txt, OCENA.txt}
OR_A = (use_junk=false, use_jump=false, use_noisy=true, rfi_safe=false, pad=0)     # A = zle odcinki, bez marginesu
OR_AB = (use_junk=true, use_jump=false, use_noisy=true, rfi_safe=false, pad=0)     # A + B (pojedyncze smieci)
OR_EVAL_VARIANTS = [
    ("v9 (obecne, z marginesem)", NamedTuple()),
    ("A + TPA 4σ", merge(OR_A, (rfi_tpa_k=4.0,))),
    ("A + TPA 6σ", merge(OR_A, (rfi_tpa_k=6.0,))),
    ("A + TPA 4σ i kształt >2x", merge(OR_A, (rfi_tpa_k=4.0, rfi_tpa_struct=2.0))),
    ("A + kształt >2x i z>6", merge(OR_A, (rfi_struct=2.0, rfi_struct_z=6.0))),
    ("A + B + TPA 4σ", merge(OR_AB, (rfi_tpa_k=4.0,))),
    ("A + B + kształt >2x i z>6", merge(OR_AB, (rfi_struct=2.0, rfi_struct_z=6.0))),
]
function _or_eval_setup!()
    global S8_DIR = F7_DIR * "one_rule_eval/"
    global OR_RULES = [(n, :model_ok, kw) for (n, kw) in OR_EVAL_VARIANTS]
    global S8_VARIANTS = [(r[1], NamedTuple(), :x, 0.0) for r in OR_RULES]
end
function one_rule_eval(; n=32)
    isdefined(Main, :_run_parallel) || Base.include(Main, "modules/zap_verify.jl")
    _or_eval_setup!(); mkpath(S8_DIR)
    codes = [F7_PRE * "include(\"tools/one_rule.jl\"); _or_eval_setup!(); Base.invokelatest(_or_part, $k, $n)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, S8_DIR * "logs", "ore"; every=120)
    one_rule_eval_summary()
end
"Ocena jakosci wyciec (co wycinaja) + standardowe podsumowanie."
function one_rule_eval_summary(; tie=1.0)
    _or_eval_setup!()
    simple_v8_summary(; title="Proste zasady RFI bez marginesów - porównanie z v9", basename_="v9")
    R = _or_rows(S8_DIR)
    L = ["# Co realnie wycinają proste zasady RFI (wszystkie pulsary)\n",
         "Wycięte z czystym tłem = impulsy wycięte jako RFI, a tło bez kształtu (< 1,5× typowego) i szum w normie (z < 2): prawdopodobnie niepotrzebnie.",
         "Zostawione złe = niewycięte impulsy z wyraźnie złym tłem (kształt > 3× typowego albo z > 10): niewycięte RFI.",
         "TPA = ile impulsów z listy Xiaoxi (stats.csv) wycinamy też my.\n",
         "| zasada | wycięte RFI | z czystym tłem | % | zostawione złe | TPA odtworzone | lepiej : gorzej : remis vs ręcznie | średni zysk S/N |",
         "|---|---|---|---|---|---|---|---|"]
    for (name, _) in OR_EVAL_VARIANTS
        rr = get(R, name, nothing); isnothing(rr) && continue
        sumk(k) = sum(get(r, k, 0.0) for r in values(rr))
        nr, nc, lb, tn, th = sumk("n_rfi"), sumk("rfi_clean"), sumk("left_bad"), sumk("tpa_n"), sumk("tpa_hit")
        dm = [100 * (r["snr"] / r["snr_manual"] - 1) for r in values(rr) if isfinite(get(r, "snr_manual", NaN)) && r["snr_manual"] > 0]
        g = [100 * (r["snr"] / r["snr_none"] - 1) for r in values(rr) if r["snr_none"] >= 3 && isfinite(r["snr"])]
        push!(L, @sprintf("| %s | %d | %d | %.0f%% | %d | %.0f%% | %d : %d : %d | %+.1f%% |", name, nr, nc, 100nc / max(nr, 1), lb,
                          100th / max(tn, 1), count(>(tie), dm), count(<(-tie), dm), count(x -> abs(x) <= tie, dm), mean(g)))
    end
    write(S8_DIR * "OCENA.txt", join(L, "\n") * "\n"); println(join(L, "\n")); println("\nZapisane: $(S8_DIR)OCENA.txt")
end
quick_eval() = quick_rfi_test(["J1012-5857", "J0820-1350", "J1110-5637", "J0932-3217", "J1915+0752", "J1844-0433", "J1842+0257"];
                              variants=OR_EVAL_VARIANTS[2:end], out=F7_DIR * "one_rule_quick_eval/")

# ---- To samo z poprawionym krokiem A (mediana szumu w oknie zamiast sredniej: pojedynczy glosny impuls nie wycina
# calego okna czystych sasiadow). Kolumna "z czystym tłem w odcinkach" pokazuje, ile czystych wycina sam krok A.
#   include("tools/one_rule.jl"); one_rule_eval2()          # -> final_v7/one_rule_eval2/OCENA.txt
OR_Am = merge(OR_A, (noisy_med=true,)); OR_ABm = merge(OR_AB, (noisy_med=true,))
OR_EVAL2_VARIANTS = [
    ("v9 (obecne, z marginesem)", NamedTuple()),
    ("A(med) + TPA 4σ", merge(OR_Am, (rfi_tpa_k=4.0,))),
    ("A(med) + B + TPA 4σ", merge(OR_ABm, (rfi_tpa_k=4.0,))),
    ("A(med) + B + TPA 4σ i kształt >2x", merge(OR_ABm, (rfi_tpa_k=4.0, rfi_tpa_struct=2.0))),
    ("A(med) + B + kształt >2x i z>6", merge(OR_ABm, (rfi_struct=2.0, rfi_struct_z=6.0))),
    ("A(śr.) + B + TPA 4σ (jak poprzednio)", merge(OR_AB, (rfi_tpa_k=4.0,))),
]
function _or_eval2_setup!()
    global S8_DIR = F7_DIR * "one_rule_eval2/"
    global OR_RULES = [(n, :model_ok, kw) for (n, kw) in OR_EVAL2_VARIANTS]
    global S8_VARIANTS = [(r[1], NamedTuple(), :x, 0.0) for r in OR_RULES]
    global OR_EVAL_VARIANTS_CUR = OR_EVAL2_VARIANTS
end
function one_rule_eval2(; n=32)
    isdefined(Main, :_run_parallel) || Base.include(Main, "modules/zap_verify.jl")
    _or_eval2_setup!(); mkpath(S8_DIR)
    codes = [F7_PRE * "include(\"tools/one_rule.jl\"); _or_eval2_setup!(); Base.invokelatest(_or_part, $k, $n)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, S8_DIR * "logs", "ore2"; every=120)
    one_rule_eval2_summary()
end
function one_rule_eval2_summary(; tie=1.0)
    _or_eval2_setup!()
    simple_v8_summary(; title="Proste zasady RFI (krok A z mediana) - porównanie z v9", basename_="v9")
    R = _or_rows(S8_DIR)
    L = ["# Co realnie wycinają proste zasady RFI, krok A z medianą (wszystkie pulsary)\n",
         "| zasada | wycięte RFI | z czystym tłem | w tym w odcinkach (A) | % czystych | zostawione złe | TPA odtworzone | lepiej : gorzej : remis vs ręcznie | średni zysk S/N |",
         "|---|---|---|---|---|---|---|---|---|"]
    for (name, _) in OR_EVAL2_VARIANTS
        rr = get(R, name, nothing); isnothing(rr) && continue
        sumk(k) = sum(get(r, k, 0.0) for r in values(rr))
        nr, nc, no, lb, tn, th = sumk("n_rfi"), sumk("rfi_clean"), sumk("rfi_clean_odc"), sumk("left_bad"), sumk("tpa_n"), sumk("tpa_hit")
        dm = [100 * (r["snr"] / r["snr_manual"] - 1) for r in values(rr) if isfinite(get(r, "snr_manual", NaN)) && r["snr_manual"] > 0]
        g = [100 * (r["snr"] / r["snr_none"] - 1) for r in values(rr) if r["snr_none"] >= 3 && isfinite(r["snr"])]
        push!(L, @sprintf("| %s | %d | %d | %d | %.0f%% | %d | %.0f%% | %d : %d : %d | %+.1f%% |", name, nr, nc, no, 100nc / max(nr, 1), lb,
                          100th / max(tn, 1), count(>(tie), dm), count(<(-tie), dm), count(x -> abs(x) <= tie, dm), mean(g)))
    end
    write(S8_DIR * "OCENA.txt", join(L, "\n") * "\n"); println(join(L, "\n")); println("\nZapisane: $(S8_DIR)OCENA.txt")
end

# ---- Wykres porownawczy zasad RFI: jeden stos impulsow, obok kolumna znacznikow dla kazdej zasady
# (pomaranczowe - zly odcinek, czerwone - pojedynczy impuls RFI, szare - null, jasnoczerwone - margines).
#   include("tools/one_rule.jl"); compare_rules_view()      # -> final_v7/porownanie_zasad.pdf
OR_VIEW_VARIANTS = [
    ("v9", NamedTuple()),
    ("A+B+TPA4σ", merge(OR_ABm, (rfi_tpa_k=4.0,))),
    ("A+B+TPA4σ i kształt", merge(OR_ABm, (rfi_tpa_k=4.0, rfi_tpa_struct=2.0))),
    ("A+B+kształt i z>6", merge(OR_ABm, (rfi_struct=2.0, rfi_struct_z=6.0))),
]
function compare_rules_view(psrs=["J0932-3217", "J1012-5857", "J0820-1350", "J1110-5637", "J1915+0752", "J1844-0433",
                                  "J1842+0257", "J1057-7914", "J1404+1159", "J0211-8159"];
                            variants=OR_VIEW_VARIANTS, out=F7_DIR * "porownanie_zasad.pdf")
    _VD = SpaTs.Data
    PyPlot = Main.PyPlot
    PdfPages = Main.PyCall.pyimport("matplotlib.backends.backend_pdf").PdfPages
    mkpath(dirname(out)); pdf = PdfPages(out)
    cols = Dict(4 => "darkorange", 2 => "red", 1 => "red", 3 => "grey", 6 => "grey", 5 => "salmon")
    for psr in psrs
        try
            p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
            D = Float64.(_f6_load(psr, p)); N, nb = size(D)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            m = _nr_model(psr); extra = isnothing(m) ? Int[] : m.rows[m.P .>= 0.5]
            F = _VD.snr_fast(D, bst, bend); s0 = _VD.profile_snr(F)
            res = [(l, _VD.auto_zap(D, bst, bend; null_mode=:model, null_extra=extra, rfi_struct=Inf, kw...)) for (l, kw) in variants]
            nv = length(res)
            fig = PyPlot.figure(figsize=(12, 8))
            fig.suptitle("$psr: co wycina każda zasada (pomarańczowe: zły odcinek, czerwone: pojedynczy impuls RFI, szare: null, jasne: margines)", fontsize=8)
            wcol = 0.045
            live = [any(!iszero, @view D[i, :]) for i in 1:N]; v = vec(D[live, :])
            for (k, (l, a)) in enumerate(res)
                ax = fig.add_axes([0.04 + (k - 1) * wcol, 0.22, wcol * 0.8, 0.7])
                for r in (4, 2, 1, 3, 6, 5)
                    ys = findall(==(r), a.reason)
                    isempty(ys) || ax.scatter(fill(0, length(ys)), ys, marker="_", s=40, c=cols[r], lw=1.2)
                end
                ax.set_xlim(-1, 1); ax.set_ylim(0.5, N + 0.5); ax.set_xticks([])
                ax.set_title(l, fontsize=5, rotation=60, ha="left")
                k == 1 ? ax.tick_params(labelsize=6) : ax.set_yticks([])
            end
            axs = fig.add_axes([0.04 + nv * wcol + 0.01, 0.22, 0.62 - nv * wcol, 0.7])
            axs.imshow(D, origin="lower", aspect="auto", cmap="Greys_r", vmin=quantile(v, 0.01), vmax=quantile(v, 0.995),
                       interpolation="nearest", extent=[-0.5, nb - 0.5, 0.5, N + 0.5])
            axs.axvline(bst, c="orange", lw=0.5); axs.axvline(bend, c="orange", lw=0.5)
            axs.tick_params(labelleft=false, labelsize=6); axs.set_xlabel("bin (pomarańczowe linie: okno impulsu)", fontsize=7)
            # struktura tla i szum (pomocniczo)
            ks = _VD.rfi_struct_ratio(D, bst, bend); zg = _VD.rfi_rms_z(D, bst, bend; local_win=0)
            axk = fig.add_axes([0.69, 0.22, 0.12, 0.7], sharey=axs)
            axk.plot(ks, 1:N, c="grey", lw=0.5); axk.axvline(2, c="k", ls="--", lw=0.5); axk.set_xscale("log")
            axk.set_xlabel("kształt tła / typowy", fontsize=6); axk.tick_params(labelsize=5, labelleft=false)
            axz = fig.add_axes([0.83, 0.22, 0.12, 0.7], sharey=axs)
            axz.plot(clamp.(zg, -5, 1000), 1:N, c="grey", lw=0.5); axz.axvline(4, c="k", ls="--", lw=0.5); axz.set_xscale("symlog", linthresh=10)
            axz.set_xlabel("szum tła z (cała obs.)", fontsize=6); axz.tick_params(labelsize=5, labelleft=false)
            # tabela liczb
            axt = fig.add_axes([0.04, 0.02, 0.9, 0.15]); axt.axis("off")
            lines = [@sprintf("%-22s  odcinki %4d   pojedyncze RFI %4d   nulle %4d   margines %3d   razem %4d   S/N %.1f -> %.1f (%+.1f%%)",
                              l, count(==(4), a.reason), count(r -> r in (1, 2), a.reason), count(r -> r in (3, 6), a.reason),
                              count(==(5), a.reason), length(a.bad), s0, _VD.profile_snr(F; drop=a.bad), 100(_VD.profile_snr(F; drop=a.bad) / s0 - 1))
                     for (l, a) in res]
            axt.text(0, 1, join(lines, "\n"), fontsize=7, family="monospace", va="top", transform=axt.transAxes)
            pdf.savefig(fig); PyPlot.close(fig)
            println(psr, " ok")
        catch e
            e isa InterruptException && rethrow()
            println("$psr BŁĄD: ", first(sprint(showerror, e), 200))
        end
    end
    pdf.close(); println("Wykres: $out")
end

# ---- Lepszy wykres porownawczy: osobny stos impulsow dla kazdej zasady, wyciete impulsy przykryte pasem przez cala
# szerokosc (turkus - pojedynczy RFI, zielony - zly odcinek, fiolet - null, zolty - margines). Dane w mapie "inferno".
#   include("tools/one_rule.jl"); compare_rules_view2()     # -> final_v7/porownanie_zasad2.pdf
function compare_rules_view2(psrs=["J0932-3217", "J1012-5857", "J0820-1350", "J1110-5637", "J1915+0752", "J1844-0433",
                                   "J1842+0257", "J1057-7914", "J1404+1159", "J0211-8159"];
                             variants=OR_VIEW_VARIANTS, out=F7_DIR * "porownanie_zasad2.pdf", zoom=nothing)
    _VD = SpaTs.Data
    PyPlot = Main.PyPlot
    PdfPages = Main.PyCall.pyimport("matplotlib.backends.backend_pdf").PdfPages
    patches = Main.PyCall.pyimport("matplotlib.patches")
    mkpath(dirname(out)); pdf = PdfPages(out)
    rgba = Dict(2 => (0.0, 0.9, 0.9, 0.55), 1 => (0.0, 0.9, 0.9, 0.55), 4 => (0.2, 0.9, 0.2, 0.55),
                3 => (0.7, 0.3, 1.0, 0.55), 6 => (0.7, 0.3, 1.0, 0.55), 5 => (1.0, 0.95, 0.2, 0.55))
    for psr in psrs
        try
            p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
            D = Float64.(_f6_load(psr, p)); N, nb = size(D)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            m = _nr_model(psr); extra = isnothing(m) ? Int[] : m.rows[m.P .>= 0.5]
            F = _VD.snr_fast(D, bst, bend); s0 = _VD.profile_snr(F)
            res = [(l, _VD.auto_zap(D, bst, bend; null_mode=:model, null_extra=extra, rfi_struct=Inf, kw...)) for (l, kw) in variants]
            live = [any(!iszero, @view D[i, :]) for i in 1:N]; v = vec(D[live, :])
            lo, hi = quantile(v, 0.02), quantile(v, 0.995)
            rows = isnothing(zoom) ? (1:N) : zoom
            nv = length(res)
            fig = PyPlot.figure(figsize=(4.2 * nv, 9))
            fig.suptitle("$psr - co wycina każda zasada", fontsize=11)
            for (k, (l, a)) in enumerate(res)
                ax = fig.add_subplot(1, nv, k)
                ax.imshow(D, origin="lower", aspect="auto", cmap="inferno", vmin=lo, vmax=hi, interpolation="nearest",
                          extent=[-0.5, nb - 0.5, 0.5, N + 0.5])
                ov = zeros(N, 1, 4)
                for i in 1:N
                    r = Int(a.reason[i]); haskey(rgba, r) && (ov[i, 1, :] .= rgba[r])
                end
                ax.imshow(ov, origin="lower", aspect="auto", interpolation="nearest", extent=[-0.5, nb - 0.5, 0.5, N + 0.5])
                ax.axvline(bst, c="white", lw=0.6, ls="--"); ax.axvline(bend, c="white", lw=0.6, ls="--")
                ax.set_ylim(first(rows) - 0.5, last(rows) + 0.5)
                s1 = _VD.profile_snr(F; drop=a.bad)
                ax.set_title(@sprintf("%s\nRFI %d, odcinki %d, nulle %d%s\nS/N %.0f -> %.0f (%+.1f%%)", l,
                                      count(r -> r in (1, 2), a.reason), count(==(4), a.reason), count(r -> r in (3, 6), a.reason),
                                      count(==(5), a.reason) > 0 ? ", margines $(count(==(5), a.reason))" : "",
                                      s0, s1, 100(s1 / s0 - 1)), fontsize=9)
                ax.set_xlabel("bin", fontsize=8); k == 1 && ax.set_ylabel("impuls (od 1)", fontsize=8); ax.tick_params(labelsize=7)
            end
            hs = [patches.Patch(color=c[1:3], alpha=0.7, label=t) for (c, t) in ((rgba[2], "pojedynczy impuls RFI"), (rgba[4], "zły odcinek"),
                                                                                (rgba[3], "null"), (rgba[5], "margines"))]
            fig.legend(handles=hs, loc="lower center", ncol=4, fontsize=9, frameon=false)
            fig.subplots_adjust(left=0.05, right=0.99, top=0.88, bottom=0.08, wspace=0.12)
            pdf.savefig(fig); PyPlot.close(fig)
            println(psr, " ok")
        catch e
            e isa InterruptException && rethrow()
            println("$psr BŁĄD: ", first(sprint(showerror, e), 200))
        end
    end
    pdf.close(); println("Wykres: $out")
end

# ---- Wykres "przed / po" dla kazdej zasady (styl jak pulsar_autozap): jedna zasada na strone, lewy panel oryginal,
# prawy dane po wycieciu (wyciete impulsy biale), miedzy nimi pasek z powodem. Strony pulsar po pulsarze: v9 i
# kolejne zasady dla tego samego pulsara - przewijajac widac, co ZOSTALO w danych.
#   include("tools/one_rule.jl"); compare_rules_pages()      # -> final_v7/porownanie_zasad3.pdf
function compare_rules_pages(psrs=["J0932-3217", "J1012-5857", "J0820-1350", "J1110-5637", "J1915+0752", "J1844-0433",
                                   "J1842+0257", "J1057-7914", "J1404+1159", "J0211-8159"];
                             variants=OR_VIEW_VARIANTS, out=F7_DIR * "porownanie_zasad3.pdf", zoom=nothing)
    _VD = SpaTs.Data
    PyPlot = Main.PyPlot
    PdfPages = Main.PyCall.pyimport("matplotlib.backends.backend_pdf").PdfPages
    patches = Main.PyCall.pyimport("matplotlib.patches")
    mkpath(dirname(out)); pdf = PdfPages(out)
    col = Dict(1 => "red", 2 => "red", 4 => "darkorange", 3 => "grey", 6 => "grey", 5 => "gold")
    for psr in psrs
        try
            p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
            D = Float64.(_f6_load(psr, p)); N, nb = size(D)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            m = _nr_model(psr); extra = isnothing(m) ? Int[] : m.rows[m.P .>= 0.5]
            F = _VD.snr_fast(D, bst, bend); s0 = _VD.profile_snr(F)
            live = [any(!iszero, @view D[i, :]) for i in 1:N]; v = vec(D[live, :])
            lo, hi = quantile(v, 0.01), quantile(v, 0.995)
            rows = isnothing(zoom) ? (1:N) : zoom
            w = max(1, round(Int, 0.15 * (bend - bst + 1)))
            xl = (max(1, bst - 4w), min(nb, bend + 4w))          # przyblizenie na okno impulsu z kawalkiem tla
            for (l, kw) in variants
                a = _VD.auto_zap(D, bst, bend; null_mode=:model, null_extra=extra, rfi_struct=Inf, kw...)
                cut = falses(N); cut[a.bad] .= true
                A = copy(D); A[cut, :] .= NaN
                s1 = _VD.profile_snr(F; drop=a.bad)
                fig = PyPlot.figure(figsize=(11, 8))
                fig.suptitle(@sprintf("%s - zasada: %s\nwycięte %d z %d  (RFI %d, złe odcinki %d, nulle %d%s)   S/N profilu %.1f -> %.1f",
                                      psr, l, length(a.bad), N, count(r -> r in (1, 2), a.reason), count(==(4), a.reason),
                                      count(r -> r in (3, 6), a.reason), count(==(5), a.reason) > 0 ? ", margines $(count(==(5), a.reason))" : "",
                                      s0, s1), fontsize=10)
                cm = PyPlot.get_cmap("afmhot").copy(); cm.set_bad("white")
                ax1 = fig.add_axes([0.07, 0.30, 0.40, 0.58]); ax2 = fig.add_axes([0.53, 0.30, 0.40, 0.58], sharey=ax1)
                axm = fig.add_axes([0.475, 0.30, 0.02, 0.58], sharey=ax1)
                for (ax, M, t) in ((ax1, D, "PRZED (oryginał)"), (ax2, A, "PO (wycięte impulsy białe)"))
                    ax.imshow(M, origin="lower", aspect="auto", cmap=cm, vmin=lo, vmax=hi, interpolation="nearest",
                              extent=[0.5, nb + 0.5, 0.5, N + 0.5])
                    ax.set_xlim(xl...); ax.set_ylim(first(rows) - 0.5, last(rows) + 0.5)
                    ax.set_title(t, fontsize=10); ax.set_xlabel("bin", fontsize=8); ax.tick_params(labelsize=7)
                end
                ax1.set_ylabel("impuls (od 1)", fontsize=8); ax2.tick_params(labelleft=false)
                for r in (4, 1, 2, 3, 6, 5)
                    ys = findall(==(r), a.reason)
                    isempty(ys) || axm.barh(ys, 1, height=1.0, color=col[r], lw=0)
                end
                axm.set_xlim(0, 1); axm.set_xticks([]); axm.tick_params(labelleft=false)
                # profile przed / po
                axp = fig.add_axes([0.07, 0.07, 0.86, 0.16])
                pr0 = vec(mean(D[live, :], dims=1)); keep = live .& .!cut
                pr1 = any(keep) ? vec(mean(D[keep, :], dims=1)) : zeros(nb)
                pm = max(maximum(pr0), eps())
                axp.plot(1:nb, pr0 ./ pm, c="grey", lw=0.8, label="przed"); axp.plot(1:nb, pr1 ./ pm, c="black", lw=0.9, label="po")
                axp.set_xlim(xl...); axp.legend(fontsize=7); axp.tick_params(labelsize=7); axp.set_xlabel("bin", fontsize=8)
                hs = [patches.Patch(color=c, label=t) for (c, t) in (("red", "RFI"), ("darkorange", "zły odcinek"), ("grey", "null"), ("gold", "margines"))]
                fig.legend(handles=hs, loc="upper right", fontsize=8, frameon=false)
                pdf.savefig(fig); PyPlot.close(fig)
            end
            println(psr, " ok")
        catch e
            e isa InterruptException && rethrow()
            println("$psr BŁĄD: ", first(sprint(showerror, e), 200))
        end
    end
    pdf.close(); println("Wykres: $out  (dla każdego pulsara $(length(variants)) strony: ", join(first.(variants), ", "), ")")
end
