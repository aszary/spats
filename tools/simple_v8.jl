# Test "prostszego" wycinania na WSZYSTKICH pulsarach (pazur2 bez zmian):
#   1. RFI z calej obserwacji (odniesienie = mediana szumu wszystkich impulsow) zamiast 31 sasiadow
#   2. bez marginesu (tylko impulsy realnie wykryte jako RFI albo null)
#   3. nulle tylko z modelu energii, P_null >= 0.8 (jedna regula); pulsar bez wiarygodnego modelu: bez nulli
#   4. wszystko naraz ("prosty")
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/simple_v8.jl"); simple_v8()
# ~20-30 min na 32 procesach, wznawialne. Wynik: final_v7/simple_v8/WYNIK.txt
isdefined(Main, :F7_DIR) || include(joinpath(@__DIR__, "final_v7.jl"))
isdefined(Main, :_nr_model) || include(joinpath(@__DIR__, "null_rules_real.jl"))
S8_DIR = F7_DIR * "simple_v8/"
# nazwa => (ustawienia auto_zap, tryb nulli, prog P_null)
S8_VARIANTS = [
    ("v7 (obecny)", NamedTuple(), :union, 0.95),
    ("RFI z całej obserwacji", (rfi_win=0,), :union, 0.95),
    ("bez marginesu", (pad=0,), :union, 0.95),
    ("nulle tylko z modelu P>=0.8", NamedTuple(), :model, 0.8),
    ("prosty (wszystko naraz)", (rfi_win=0, pad=0), :model, 0.8),
    ("nulle tylko z modelu P>=0.5 (publikacja)", NamedTuple(), :model, 0.5),
    ("publikacja P>=0.5 + RFI z całości + bez marginesu", (rfi_win=0, pad=0), :model, 0.5),
]
S8_HEADER = "pulsar,variant,N,n_cut,n_rfi,n_null,snr_none,snr,snr_manual,snr_rfi,kept_z3_local,kept_z3_global," *
            "emis_cut,emis_all,null_cut,null_all,E_kept,E_emis"

function _s8_part(k, n)
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
            E = [sum(@view D[i, on]) - length(on) * mean(@view D[i, off]) for i in 1:N]
            zl = _VD.rfi_rms_z(D, bst, bend)
            zg = _VD.rfi_rms_z(D, bst, bend; local_win=0)
            m = _nr_model(psr)                           # nothing = brak wiarygodnego modelu
            P = fill(NaN, N); isnothing(m) || (P[m.rows] .= m.P)
            snr_none = _VD.profile_snr(D, bst, bend)
            snr_man = something(get(get(p, "zap_cmp", Dict()), "snr_maciej", nothing), NaN)
            for (name, kw, mode, thr) in S8_VARIANTS
                # tryb :model bez wiarygodnego modelu: pusta lista -> zadnych nulli (jedna regula, bez automatu)
                extra = isnothing(m) ? (mode == :model ? Int[] : nothing) : m.rows[m.P .>= thr]
                a = _VD.auto_zap(D, bst, bend; null_mode=mode, null_extra=extra, kw...)
                sp = _VD.split_zaps(a.reason)
                rfi = sp.rfi; nul = sp.null
                cut = Set(a.bad); rs = Set(rfi)
                U = [i for i in 1:N if isfinite(P[i]) && !(i in rs)]
                C = [i for i in U if i in cut]; K = [i for i in U if !(i in cut)]
                ea, na = isempty(U) ? (NaN, NaN) : (sum(1 .- P[U]), sum(P[U]))
                Ee = isempty(U) || !(ea > 0) ? NaN : sum((1 .- P[U]) .* E[U]) / ea
                row = [psr, name, N, length(a.bad), length(rfi), length(nul), snr_none,
                       _VD.profile_snr(D, bst, bend; drop=a.bad), snr_man, _VD.profile_snr(D, bst, bend; drop=rfi),
                       count(i -> !(i in cut) && zl[i] > 3, 1:N), count(i -> !(i in cut) && zg[i] > 3, 1:N),
                       isempty(C) ? 0.0 : sum(1 .- P[C]), ea, isempty(C) ? 0.0 : sum(P[C]), na,
                       isempty(K) ? NaN : mean(E[K]), Ee]
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

function simple_v8_summary(; tie=1.0)
    R = Dict{String,Dict{String,Dict{String,Float64}}}()          # wariant => pulsar => wiersz
    for f in filter(f -> occursin(r"^part_\d+\.csv$", basename(f)), readdir(S8_DIR; join=true))
        H = split(first(eachline(f)), ",")
        for l in Iterators.drop(eachline(f), 1)
            v = split(l, ","); length(v) == length(H) || continue
            get!(R, v[2], Dict{String,Dict{String,Float64}}())[v[1]] =
                Dict(h => something(tryparse(Float64, x), NaN) for (h, x) in zip(H[3:end], v[3:end]))
        end
    end
    base = get(R, S8_VARIANTS[1][1], Dict{String,Dict{String,Float64}}())
    L = String[]; pr(s) = (println(s); push!(L, s))
    pr("# Prostsze wycinanie: RFI z całej obserwacji, bez marginesu, nulle tylko z modelu (P>=0.8)\n")
    pr("Pulsary: $(length(base)). Zmiana względem v7 na tych samych danych; „RFI zostawione” = impulsy niewycięte z szumem poza impulsem > 3σ (lokalnie / względem całej obserwacji).\n")
    pr("| wariant | gorzej niż v7 o >1% S/N | lepiej o >1% | najgorzej | vs ręcznie lepszy : gorszy : remis | średni zysk S/N | wycięte | RFI | nulle | RFI zostawione: lokalnie / całość | dane tylko po RFI: tracą >1% S/N |")
    pr("|---|---|---|---|---|---|---|---|---|---|---|")
    num(r, k) = get(r, k, NaN)
    for (name, _, _, _) in S8_VARIANTS
        rr = get(R, name, nothing); isnothing(rr) && continue
        ps = [q for q in keys(rr) if haskey(base, q)]
        d = [(q, 100 * (rr[q]["snr"] / base[q]["snr"] - 1)) for q in ps if base[q]["snr"] > 0 && isfinite(rr[q]["snr"])]
        worse = count(x -> x[2] < -1, d); better = count(x -> x[2] > 1, d)
        wq = isempty(d) ? ("", NaN) : d[argmin(last.(d))]
        dm = [100 * (r["snr"] / r["snr_manual"] - 1) for r in values(rr) if isfinite(num(r, "snr_manual")) && r["snr_manual"] > 0]
        g = [100 * (r["snr"] / r["snr_none"] - 1) for r in values(rr) if r["snr_none"] >= 3 && isfinite(r["snr"])]
        Nt = sum(r["N"] for r in values(rr))
        lost = count(r -> r["snr_none"] >= 3 && r["snr_rfi"] < 0.99 * r["snr_none"], values(rr))
        pr(@sprintf("| %s | %d | %d | %+.1f%% %s | %d : %d : %d | %+.1f%% | %.1f%% | %.2f%% | %.1f%% | %d / %d | %d |", name, worse, better,
                    wq[2], wq[1], count(>(tie), dm), count(<(-tie), dm), count(x -> abs(x) <= tie, dm), mean(g),
                    100 * sum(r["n_cut"] for r in values(rr)) / Nt, 100 * sum(r["n_rfi"] for r in values(rr)) / Nt,
                    100 * sum(r["n_null"] for r in values(rr)) / Nt,
                    sum(r["kept_z3_local"] for r in values(rr)), sum(r["kept_z3_global"] for r in values(rr)), lost))
    end
    pr("\n## Nulle (tylko pulsary z wiarygodnym modelem; szacunek z P_null, dla reguł z automatem za surowy)\n")
    pr("| wariant | wycięta emisja | nulle znalezione | profil jaśniejszy >1% | najgorzej (+) | profil ciemniejszy >1% | najgorzej (−) |")
    pr("|---|---|---|---|---|---|---|")
    for (name, _, _, _) in S8_VARIANTS
        rr = [r for r in values(get(R, name, Dict{String,Dict{String,Float64}}())) if isfinite(num(r, "emis_all"))]
        isempty(rr) && continue
        ec = sum(r["emis_cut"] for r in rr) / sum(r["emis_all"] for r in rr)
        nc = sum(r["null_cut"] for r in rr) / max(sum(r["null_all"] for r in rr), eps())
        b = [100 * (r["E_kept"] / r["E_emis"] - 1) for r in rr if r["null_all"] >= 5 && r["E_emis"] > 0 && isfinite(r["E_kept"])]
        pr(@sprintf("| %s | %.2f%% | %.1f%% | %d z %d | %+.1f%% | %d z %d | %+.1f%% |", name, 100ec, 100nc,
                    count(>(1), b), length(b), isempty(b) ? NaN : maximum(b), count(<(-1), b), length(b),
                    isempty(b) ? NaN : minimum(b)))
    end
    write(S8_DIR * "WYNIK.txt", join(L, "\n") * "\n")
    println("\nZapisane: $(S8_DIR)WYNIK.txt")
end

"""
Wersja końcowa wybranego wariantu w final_v7/simple_v8/pazur2_final (pazur2 bez zmian) + weryfikacja
i porównanie z metodą ręczną. Domyślnie: publikacja P>=0.5 + RFI z całości + bez marginesu. ~30-40 min.
"""
function simple_v8_build(variant="publikacja P>=0.5 + RFI z całości + bez marginesu"; n=32)
    isdefined(Main, :_run_parallel) || Base.include(Main, "modules/zap_verify.jl")
    i = findfirst(v -> v[1] == variant, S8_VARIANTS)
    isnothing(i) && error("nieznany wariant: $variant; dostępne: $(join(first.(S8_VARIANTS), " | "))")
    _, kw, mode, thr = S8_VARIANTS[i]
    out = S8_DIR * "pazur2_final/"; mkpath(out)
    write(out * ".wariant", "$variant: ustawienia=$(repr(kw)) nulle=$mode P>=$thr")
    println("\n===== WERSJA v8: $variant → $out =====")
    cfg = "SpaTs.Data.ZAP_STEPS[] = $(repr(kw)); SpaTs.Data.NULL_MODE[] = :$mode; SpaTs.Data.NULL_THR[] = $thr; " *
          "SpaTs.Data.NULL_MODEL_DIR[] = \"$(F7_DIR)model/\"; "
    cf = "SpaTs.chan_first_list([\"$(F6_ROOT)zap_verify_v3/\", \"$(F6_ROOT)zap_verify/\"])"
    codes = [F7_PRE * cfg * "SpaTs.autozap_all(out=\"$out\", only=$(repr(_f6_psrs()[k:n:end])), skip_done=true, keep_manual=false, " *
             "reuse_window=true, window_root=\"$(F6_ROOT)pazur2/\", chan_first=$cf)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, S8_DIR * "logs", "build"; every=120)
    vd = S8_DIR * "zap_verify_final/"
    codes = [F7_PRE * cfg * "Base.invokelatest(Main.zap_verify; pazur2=\"$out\", out=\"$vd\", part=\"$k/$n\", chans=false, save_rfi=true)"
             for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, S8_DIR * "logs", "verify"; every=120)
    Base.invokelatest(Main.zap_verify_summary; out=vd)
    SpaTs.compare_zap_methods(mine=out, outdir=S8_DIR * "zap_comparison_final/")
    println("Potem, po decyzji: simple_v8_adopt()")
end

"Przyjęcie v8: pazur2 → pazur2_przed_v8 (albo _2...), simple_v8/pazur2_final → pazur2."
function simple_v8_adopt()
    src, dst, bk = S8_DIR * "pazur2_final", F6_ROOT * "pazur2", F6_ROOT * "pazur2_przed_v8"
    isdir(src) || error("brak $src - najpierw simple_v8_build()")
    k = 1
    while isdir(bk); k += 1; bk = F6_ROOT * "pazur2_przed_v8_$k"; end
    mv(dst, bk); mv(src, dst)
    println("Przyjęte v8: $(strip(read(dst * "/.wariant", String))). Stara wersja: $bk.")
    println("Domyślne w kodzie trzeba jeszcze zmienić (napisz Claude'owi).")
end

function simple_v8(; n=32)
    isdefined(Main, :_run_parallel) || Base.include(Main, "modules/zap_verify.jl")
    mkpath(S8_DIR)
    codes = [F7_PRE * "include(\"tools/simple_v8.jl\"); Base.invokelatest(_s8_part, $k, $n)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, S8_DIR * "logs", "s8"; every=120)
    simple_v8_summary()
end
