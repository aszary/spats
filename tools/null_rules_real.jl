# Ile prawdziwej EMISJI i ile NULLI wycina kazda regula nulli NA PRAWDZIWYCH DANYCH (bez danych syntetycznych).
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/null_rules_real.jl"); null_rules_real()
# Tylko pulsary z wiarygodnym modelem nulli (final_v7/model, status ok). Model daje dla kazdego impulsu
# P = prawdopodobienstwo nulla, czyli 1-P = prawdopodobienstwo emisji. Dla kazdej reguly (te same 8 co w v6):
#   wycieta emisja   = suma (1-P) po wycietych impulsach (bez RFI) / suma (1-P) po wszystkich (bez RFI)
#   kompletnosc      = suma P po wycietych / suma P po wszystkich
#   zmiana profilu   = srednia energia impulsow, ktore zostaly / srednia energia samej emisji - 1
#                      (energia emisji z modelu: suma (1-P)E / suma (1-P)); > 0: wycieta slaba emisja,
#                      < 0: zostaly nulle
# Wznawialne; ~10-15 min na 32 procesach. Wynik: final_v7/null_rules_real/WYNIK.txt
isdefined(Main, :F7_DIR) || include(joinpath(@__DIR__, "final_v7.jl"))
NR_DIR = F7_DIR * "null_rules_real/"
NR_MODEL = F7_DIR * "model/"
NR_HEADER = "pulsar,variant,N,n_cut,emis_cut,emis_all,null_cut,null_all,E_kept,E_emis,snr_none,snr"

"P_null (wiersze od 1) i status modelu."
function _nr_model(psr)
    st = "brak"
    for l in Iterators.drop(eachline(NR_MODEL * "nulls_model.csv"), 1)
        r = split(l, ","); r[1] == psr && (st = r[2]; break)
    end
    pf = NR_MODEL * "nulls_model/$psr.csv"
    (st != "ok" || !isfile(pf)) && return nothing
    rows = Int[]; P = Float64[]
    for l in Iterators.drop(eachline(pf), 1)
        r = split(l, ","); push!(rows, parse(Int, r[1]) + 1); push!(P, parse(Float64, r[3]))
    end
    return (rows=rows, P=P)
end

function _nr_part(k, n)
    _VD = SpaTs.Data
    mkpath(NR_DIR)
    f = NR_DIR * "part_$k.csv"
    done = isfile(f) ? Set(split(l, ",")[1] for l in Iterators.drop(eachline(f), 1)) : Set{String}()
    isfile(f) || open(io -> println(io, NR_HEADER), f, "w")
    psrs = filter(!in(done), _f6_psrs()[k:n:end])
    for (j, psr) in enumerate(psrs)
        t0 = time()
        try
            m = _nr_model(psr)
            if isnothing(m)
                open(io -> println(io, "$psr,-,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN"), f, "a")
                continue
            end
            p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
            D = _f6_load(psr, p)
            N = size(D, 1)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            on = bst:bend
            off = _VD.offpulse_bins(D, bst, bend)
            E = [sum(@view D[i, on]) - length(on) * mean(@view D[i, off]) for i in 1:N]
            P = fill(NaN, N); P[m.rows] .= m.P
            snr_none = _VD.profile_snr(D, bst, bend)
            for (name, marg, mode, thr) in F6_VARIANTS
                extra = mode == :auto ? nothing : m.rows[m.P .>= thr]
                a = _VD.auto_zap(D, bst, bend; null_margin=marg, null_mode=mode, null_extra=extra)
                rfi = Set(_VD.split_zaps(a.reason).rfi)
                U = [i for i in 1:N if isfinite(P[i]) && !(i in rfi)]       # impulsy z modelem, bez RFI
                cut = Set(a.bad)
                C = [i for i in U if i in cut]
                K = [i for i in U if !(i in cut)]
                ea, na = sum(1 .- P[U]), sum(P[U])
                Ee = ea > 0 ? sum((1 .- P[U]) .* E[U]) / ea : NaN
                row = [psr, name, N, length(C), sum(1 .- P[C]; init=0.0), ea, sum(P[C]; init=0.0), na,
                       isempty(K) ? NaN : mean(E[K]), Ee, snr_none, _VD.profile_snr(D, bst, bend; drop=a.bad)]
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

function null_rules_real_summary()
    H = split(NR_HEADER, ",")
    R = Dict{String,Vector{Dict{String,Float64}}}()
    for f in filter(f -> occursin(r"^part_\d+\.csv$", basename(f)), readdir(NR_DIR; join=true))
        for l in Iterators.drop(eachline(f), 1)
            v = split(l, ","); (length(v) == length(H) && v[2] != "-") || continue
            push!(get!(R, v[2], Dict{String,Float64}[]),
                  Dict(h => something(tryparse(Float64, x), NaN) for (h, x) in zip(H[3:end], v[3:end])))
        end
    end
    L = String[]; pr(s) = (println(s); push!(L, s))
    npsr = isempty(R) ? 0 : length(first(values(R)))
    pr("# Reguły nulli na PRAWDZIWYCH danych: ile emisji i ile nulli wycinają (z P_null modelu)\n")
    pr("Pulsary z wiarygodnym modelem: $npsr. Wycięta emisja = Σ(1−P) wyciętych / Σ(1−P) wszystkich (bez RFI); " *
       "kompletność = ΣP wyciętych / ΣP wszystkich; zmiana profilu = energia impulsów, które zostały / energia samej emisji − 1.\n")
    pr("| reguła | wycięta emisja (razem) | kompletność nulli (razem) | zmiana profilu: mediana | profil jaśniejszy >1% (wycięta emisja) | najgorzej (+) | profil ciemniejszy >1% (zostały nulle) | najgorzej (−) | średni zysk S/N |")
    pr("|---|---|---|---|---|---|---|---|---|")
    for (name, _, _, _) in F6_VARIANTS
        rr = get(R, name, Dict{String,Float64}[]); isempty(rr) && continue
        ec = sum(r["emis_cut"] for r in rr) / sum(r["emis_all"] for r in rr)
        nc = sum(r["null_cut"] for r in rr) / max(sum(r["null_all"] for r in rr), eps())
        # zmiana profilu tylko tam, gdzie sa nulle (>= 5 oczekiwanych) i emisja jest wyrazna
        b = [100 * (r["E_kept"] / r["E_emis"] - 1) for r in rr if r["null_all"] >= 5 && r["E_emis"] > 0 && isfinite(r["E_kept"])]
        g = [100 * (r["snr"] / r["snr_none"] - 1) for r in rr if r["snr_none"] >= 3 && isfinite(r["snr"])]
        bp = isempty(b) ? NaN : maximum(b); bn = isempty(b) ? NaN : minimum(b)
        pr(@sprintf("| %s | %.2f%% | %.1f%% | %+.2f%% | %d z %d | %+.1f%% | %d z %d | %+.1f%% | %+.1f%% |", name, 100ec, 100nc,
                    isempty(b) ? NaN : median(b), count(>(1), b), length(b), bp, count(<(-1), b), length(b), bn,
                    isempty(g) ? NaN : mean(g)))
    end
    write(NR_DIR * "WYNIK.txt", join(L, "\n") * "\n")
    println("\nZapisane: $(NR_DIR)WYNIK.txt")
end

function null_rules_real(; n=32)
    isdefined(Main, :_run_parallel) || Base.include(Main, "modules/zap_verify.jl")
    mkpath(NR_DIR)
    codes = [F7_PRE * "include(\"tools/null_rules_real.jl\"); Base.invokelatest(_nr_part, $k, $n)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, NR_DIR * "logs", "nr"; every=120)
    null_rules_real_summary()
end
