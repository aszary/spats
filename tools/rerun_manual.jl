# 13 pulsarow zrobionych wczesniej recznie (manual_zap_all) -> automat v4 z oknem on-pulse
# zapisanym przy recznym przerobieniu. Wtedy wszystkie 533 pulsary sa z automatu.
# Uruchom w Julii na kopierniku:
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/rerun_manual.jl"); rerun_manual()
# Kroki: kopia zapasowa params.json i zap_summary.csv -> automat v4 (reuse_window=true; stare reczne
# wyciecia zostaja w params jako zaps2_manual_old) -> weryfikacja tych pulsarow (+ pulsar_rfi.spCF,
# kanaly bez limitu impulsow) -> energie i model nulli dla nich od nowa -> porownanie + tabela
# "recznie wczesniej vs automat teraz". pazur nie jest ruszany.
using Printf

function _rm_dedupe_summary(f)
    L = readlines(f)
    last_ = Dict{String,Int}()
    for (i, l) in enumerate(L[2:end]); last_[split(l, ",")[1]] = i + 1; end
    keep = [1; sort(collect(values(last_)))]
    length(keep) < length(L) && write(f, join(L[keep], "\n") * "\n")
end

function rerun_manual(; root="/home/psr/data/OUTPUT/")
    t0 = time()
    pz = root * "pazur2/"
    summ = pz * "zap_summary.csv"
    st = Dict{String,String}()
    for l in Iterators.drop(eachline(summ), 1)
        r = split(l, ","); length(r) >= 2 && (st[r[1]] = r[2])
    end
    psrs = sort([k for (k, v) in st if v == "reczny"])
    println("PULSARY ZROBIONE RĘCZNIE ($(length(psrs))): ", join(psrs, " "))
    isempty(psrs) && return

    # 1. kopie zapasowe
    stamp = Dates_stamp()
    cp(summ, pz * "zap_summary_przed_rerun_$stamp.csv"; force=true)
    old = Dict{String,Any}()
    for psr in psrs
        pf = joinpath(pz, psr, "params.json")
        isfile(pf) || continue
        cp(pf, joinpath(pz, psr, "params_reczne.json"); force=true)
        old[psr] = SpaTs.Tools.read_params(pf)
    end
    println("Kopie: $(pz)zap_summary_przed_rerun_$stamp.csv, <psr>/params_reczne.json")

    # 2. automat v4 z zapisanym oknem
    cf = SpaTs.chan_first_list([root * "zap_verify_v3/", root * "zap_verify/"])
    SpaTs.autozap_all(only=psrs, skip_done=false, keep_manual=false, reuse_window=true, chan_first=cf)
    _rm_dedupe_summary(summ)

    # 3. weryfikacja (spojnosc, nulle, kanaly) + pulsar_rfi.spCF
    Base.include(Main, "modules/zap_verify.jl")
    vd = root * "zap_verify_v4_manual/"
    rm(vd; force=true, recursive=true)
    Base.invokelatest(Main.zap_verify; out=vd, only=psrs, save_rfi=true, max_chan_pulses=10^6)

    # 4. energie i model nulli od nowa dla tych pulsarow (etykiety sie zmienily)
    nout = root * "nulling_model_v4/"
    for psr in psrs
        rm(joinpath(nout, "energies", "$psr.csv"); force=true)
        rm(joinpath(nout, "fits", "$psr.v3.json"); force=true)
    end
    ENV["NULLING_OUT"] = nout
    Base.include(Main, "tools/nulling_energies.jl")
    Base.invokelatest(Main.nulling_energies_parallel)
    try
        run(addenv(`python3 tools/nulling_model.py`, "NULLING_OUT" => nout))
    catch e
        println("model nulli: ", e)
    end

    # 5. porownanie z metoda reczna (pazur) i tabela: recznie wczesniej vs automat teraz
    SpaTs.compare_zap_methods()
    st2 = Dict{String,String}()
    for l in Iterators.drop(eachline(summ), 1)
        r = split(l, ","); length(r) >= 2 && (st2[r[1]] = r[2])
    end
    r1(x) = x isa Number ? @sprintf("%.1f", x) : "-"
    println("\n===== 13 PULSARÓW: ręcznie wcześniej vs automat teraz =====")
    println(rpad("pulsar", 13), rpad("status", 9), lpad("imp.", 6), lpad("ręcz.", 7), lpad("auto", 6),
            lpad("S/N bez", 9), lpad("S/N ręcz.", 10), lpad("S/N auto", 10), lpad("zmiana", 8), lpad("RFI", 5), lpad("nulle", 6))
    nb = nw = ne = 0
    for psr in psrs
        pf = joinpath(pz, psr, "params.json")
        isfile(pf) || (println(rpad(psr, 13), "brak params.json"); continue)
        p = SpaTs.Tools.read_params(pf)
        c, co = get(p, "zap_cmp", Dict()), get(get(old, psr, Dict()), "zap_cmp", Dict())
        sa, sm = get(c, "snr_mine", nothing), get(co, "snr_mine", nothing)
        d = (sa isa Number && sm isa Number && sm != 0) ? 100 * (sa - sm) / abs(sm) : NaN
        isfinite(d) && (d >= 1 ? (nb += 1) : d <= -1 ? (nw += 1) : (ne += 1))
        nrfi = get(c, "n_rfi", nothing); nnul = get(c, "n_null", nothing)
        println(rpad(psr, 13), rpad(get(st2, psr, "?"), 9), lpad(string(get(c, "n_pulses", "-")), 6),
                lpad(string(get(co, "n_mine", "-")), 7), lpad(string(get(c, "n_mine", "-")), 6),
                lpad(r1(get(c, "snr_none", nothing)), 9), lpad(r1(sm), 10), lpad(r1(sa), 10),
                lpad(isfinite(d) ? @sprintf("%+.1f%%", d) : "-", 8), lpad(string(something(nrfi, "-")), 5),
                lpad(string(something(nnul, "-")), 6))
    end
    println("automat lepszy niż wcześniejsze ręczne: $nb, gorszy: $nw, remis (<1%): $ne")
    println("Jeśli automat gdzieś wyraźnie gorszy: stare wycięcia są w <psr>/params_reczne.json (zaps2) i w params.json (zaps2_manual_old).")
    @printf("\nKONIEC rerun_manual (%.0f min). Wklej powyższą tabelę, blok AUTOMAT vs METODA RĘCZNA i: cat %ssummary.txt\n",
            (time() - t0) / 60, vd)
end

Dates_stamp() = replace(string(round(Int, time())), r"\D" => "")
