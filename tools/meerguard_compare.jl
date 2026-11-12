# Test MeerGuard, kroki 6-8: porownanie flag MeerGuard (tools/meerguard_test.py)
# z automatem. S/N liczone TA SAMA funkcja co w automacie (Data.profile_snr, stale
# okno on-pulse) na tych samych danych (_debased_spcf z zap_verify.jl, z kanalami
# zap_chans16 jak w automacie). Nic nie zmienia w pazur/pazur2.
#
#     cd ~/software/spats && SPATS_NO_MAIN=1 julia tools/meerguard_compare.jl
#
# Wyniki: /home/psr/data/OUTPUT/meerguard_test/summary.csv, <psr>.pdf, RAPORT.md
using Statistics, Printf, DelimitedFiles
ENV["SPATS_NO_MAIN"] = "1"                 # tylko modul SpaTs, bez main()
include(joinpath(@__DIR__, "..", "spats.jl"))
Base.include(Main, joinpath(@__DIR__, "..", "modules", "zap_verify.jl"))
import JSON

const MGROOT = "/home/psr/data/OUTPUT/meerguard_test"
const MG_PAZUR2 = "/home/psr/data/OUTPUT/pazur2"
const MG_DATAROOT = "/home/psr/data/new"

"Suma czestotliwosci przez PSRCHIVE z wagami (jak .spCF): 1 polaryzacja, pam -F, debase -> (impuls, bin)."
function _mg_fsum(ar, tmp, bst, bend, tag)
    a = joinpath(tmp, "mg_$tag.spCF")
    cp(ar, a; force=true)
    try run(pipeline(`pam -p -m $a`, stdout=devnull, stderr=devnull)) catch; end
    run(pipeline(`pam -F -m $a`, stdout=devnull, stderr=devnull))
    run(pipeline(`pmod -onpulse "$bst $bend" -device "/NULL" -debase $a`; stdout=devnull, stderr=devnull))
    txt = joinpath(tmp, "mg_$tag.debase.txt")
    _VD.convert_psrfit_ascii(joinpath(tmp, "mg_$tag.debase.gg"), txt)
    return _VD.load_ascii(txt)
end

prf(a, ref) = (tp = length(intersect(a, ref));
               (jac = (u = length(union(a, ref)); u == 0 ? NaN : tp / u),
                prec = isempty(a) ? NaN : tp / length(a), rec = isempty(ref) ? NaN : tp / length(ref)))

function mg_compare()
    mf = joinpath(MGROOT, "meta.json")
    isfile(mf) || (println("Brak $mf — najpierw python3 tools/meerguard_test.py (bez błędu STOP)."); return)
    meta = JSON.parsefile(mf)
    rows = Dict{String,Any}[]
    notes = String[]
    plt = _VP.PyPlot
    PdfPages = _VP.PyCall.pyimport("matplotlib.backends.backend_pdf").PdfPages
    for (psr, st) in sort(collect(meta["pulsars"]), by=first)
        if haskey(st, "skip") || !haskey(st, "flags")
            push!(notes, "$psr: pominięty — $(get(st, "skip", get(st, "flags_error", "brak flag")))")
            continue
        end
        p = _VT.read_params(joinpath(MG_PAZUR2, psr, "params.json"))
        bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
        obs = sort(filter(d -> isdir(joinpath(MG_DATAROOT, psr, d)), readdir(joinpath(MG_DATAROOT, psr))))
        indir = joinpath(MG_DATAROOT, psr, obs[1]) * "/"
        tmp = mktempdir()
        try
            println("== $psr")
            D = _debased_spcf(indir, tmp, bst, bend; chans16=something(get(p, "zap_chans16", nothing), ""))
            N = size(D, 1)
            auto = _VD.zaps_to_rows(_VD.zaps_of(p), N)
            rfi = _VD.zaps_to_rows(_VD.zaps_rfi_of(p), N)
            null = _VD.zaps_to_rows(_VD.zaps_null_of(p), N)
            snr(drop) = _VD.profile_snr(D, bst, bend; drop=drop)
            s_none, s_auto = snr(Int[]), snr(auto)
            ar_in = joinpath(MGROOT, psr, psr * ".ar")
            s_pairs_none = _VD.profile_snr(_mg_fsum(ar_in, tmp, bst, bend, "in"), bst, bend)
            chans16 = [parse(Int, c) for c in split(something(get(p, "zap_chans16", nothing), ""))]
            live = readdlm(joinpath(MGROOT, psr, "live_in.txt"), Int)
            fl = Dict{String,Any}()
            for v in ("A", "B", "C")
                f = st["flags"][v]
                haskey(f, "error") && (push!(notes, "$psr $v: $(f["error"])"); continue)
                full = Int.(f["subints_full"]) .+ 1
                cfrac = Float64.(f["chan_frac"])
                nz = readdlm(joinpath(MGROOT, psr, v, "newzero.txt"), Int)
                # kanaly wyciete przez MeerGuard w CALOSCI nie mowia nic o impulsach:
                # bez nich "impuls z flaga" (anyf) i "impuls zly" (half: >= polowa
                # zywych kanalow oflagowana) - inaczej kazdy impuls ma flage (J0630)
                chfull = Int.(f["chans_full"]) .+ 1
                nzx = copy(nz); nzx[:, chfull] .= 0
                nlive = vec(sum(live, dims=2))
                anyf = findall(vec(any(nzx .> 0, dims=2)))
                half = findall(i -> nlive[i] > 0 && sum(@view nzx[i, :]) >= 0.5 * nlive[i], 1:size(nzx, 1))
                s_pairs = _VD.profile_snr(_mg_fsum(joinpath(MGROOT, psr, v, psr * ".mg"), tmp, bst, bend, v), bst, bend)
                mf, ma, mh = prf(full, rfi), prf(anyf, rfi), prf(half, rfi)
                top = sortperm(cfrac, rev=true)[1:3] .- 1
                only_auto = setdiff(rfi, half)
                only_mg = sort(setdiff(anyf, auto), by=i -> -sum(nzx[i, :]))
                r = Dict{String,Any}(
                    "pulsar" => psr, "wariant" => v, "N" => N, "n_rfi_auto" => length(rfi), "n_null_auto" => length(null),
                    "frac_pairs_pct" => 100 * f["frac_pairs"], "n_subints_full" => length(full),
                    "n_subints_any" => length(anyf), "chans_full" => join(Int.(f["chans_full"]), " "),
                    "jac_full" => mf.jac, "prec_full" => mf.prec, "rec_full" => mf.rec,
                    "jac_any" => ma.jac, "prec_any" => ma.prec, "rec_any" => ma.rec,
                    "n_subints_half" => length(half), "jac_half" => mh.jac, "prec_half" => mh.prec, "rec_half" => mh.rec,
                    "null_in_full" => count(in(Set(full)), null), "null_in_any" => count(in(Set(anyf)), null),
                    "null_in_half" => count(in(Set(half)), null), "snr_mg_half" => snr(half),
                    "top3_chans" => join(top, " "), "zap_chans16" => join(chans16, " "),
                    "chans16_in_top3" => count(in(Set(top)), chans16),
                    "snr_none" => s_none, "snr_auto" => s_auto, "snr_mg_full" => snr(full), "snr_mg_any" => snr(anyf),
                    "snr_pairs_none" => s_pairs_none, "snr_mg_pairs" => s_pairs,
                    "only_auto" => length(only_auto), "only_mg" => length(only_mg),
                    "only_auto_list" => _VD.pulses_to_zaps(only_auto), "only_mg_list" => _VD.pulses_to_zaps(only_mg[1:min(50, end)]),
                    "chan_frac" => join([@sprintf("%.3f", x) for x in cfrac], " "),
                    "input_unchanged" => get(st, "input_unchanged", false))
                push!(rows, r)
                fl[v] = (nz=nz, cfrac=cfrac, only_auto=only_auto, only_mg=only_mg)
                @printf("   %s: pary %.2f%%, kanały całe: %s | impulsy pełne %d / ≥połowa kanałów %d / z flagą %d | nulle (z %d): %d / %d / %d | Jaccard z RFI automatu (≥połowa): %.2f | S/N bez %.1f, automat %.1f, MG ≥połowa %.1f, MG pary %.1f (bez %.1f)\n",
                        v, r["frac_pairs_pct"], isempty(chfull) ? "-" : join(chfull .- 1, " "), length(full), length(half),
                        length(anyf), length(null), r["null_in_full"], r["null_in_half"], r["null_in_any"], mh.jac,
                        s_none, s_auto, r["snr_mg_half"], s_pairs, s_pairs_none)
            end
            haskey(fl, "A") || continue
            # ---- PDF
            pdf = PdfPages(joinpath(MGROOT, "$psr.pdf"))
            nz = fl["A"].nz
            fig, ax = plt.subplots(1, 1, figsize=(9, 11))
            ax.imshow(nz, aspect="auto", origin="lower", cmap="Greys", interpolation="nearest",
                      extent=(-0.5, size(nz, 2) - 0.5, -0.5, size(nz, 1) - 0.5))
            ax.scatter(fill(-1.2, length(rfi)), rfi .- 1, s=4, c="C3", marker="s", label="zaps2_rfi (automat)")
            ax.scatter(fill(-2.0, length(null)), null .- 1, s=4, c="C0", marker="s", label="zaps2_null (automat)")
            ax.set_xlim(-2.6, size(nz, 2) - 0.5)
            ax.set_xlabel("kanał"); ax.set_ylabel("impuls (subint)")
            ax.set_title("$psr — nowe flagi MeerGuard, wariant A (czarne); po lewej wycięcia automatu", fontsize=10)
            ax.legend(loc="upper right", fontsize=8)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

            fig, ax = plt.subplots(1, 1, figsize=(11, 5))
            x = 0:length(fl["A"].cfrac)-1
            for (k, v) in enumerate(("A", "B", "C"))
                haskey(fl, v) || continue
                ax.bar(x .+ (k - 2) * 0.27, 100 .* fl[v].cfrac, width=0.27, label="wariant $v")
            end
            for c in chans16
                ax.axvline(c, color="C3", ls="--", lw=0.8)
            end
            ax.set_xticks(x); ax.set_xlabel("kanał (czerwona linia = zap_chans16 automatu)")
            ax.set_ylabel("% oflagowanych impulsów w kanale"); ax.legend()
            ax.set_title("$psr — odsetek flag MeerGuard na kanał", fontsize=10)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

            w = bend - bst + 1
            lo, hi = max(1, bst - w), min(size(D, 2), bend + w)
            fig, axs = plt.subplots(1, 2, figsize=(12, 7))
            for (k, (lab, list)) in enumerate((("tylko MeerGuard (A, najwięcej kanałów)", fl["A"].only_mg[1:min(10, end)]),
                                                ("tylko automat (A)", fl["A"].only_auto[1:min(10, end)])))
                a = axs[k]
                if isempty(list)
                    a.text(0.5, 0.5, "brak", ha="center"); a.set_title(lab, fontsize=9); continue
                end
                M = D[list, lo:hi]
                v = vec(M)
                a.imshow(M, aspect="auto", origin="lower", cmap="viridis", interpolation="nearest",
                         vmin=quantile(v, 0.01), vmax=quantile(v, 0.995), extent=(lo, hi, -0.5, length(list) - 0.5))
                a.set_yticks(0:length(list)-1); a.set_yticklabels(string.(list .- 1), fontsize=7)
                a.axvline(bst, color="w", ls="--", lw=0.6); a.axvline(bend, color="w", ls="--", lw=0.6)
                a.set_title(lab, fontsize=9); a.set_xlabel("bin"); a.set_ylabel("impuls")
            end
            fig.suptitle("$psr — impulsy rozbieżne (dane jak w automacie, skala 1–99,5%)", fontsize=10)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
            pdf.close()
        catch e
            e isa InterruptException && rethrow()
            push!(notes, "$psr: błąd porównania — $(first(sprint(showerror, e), 200))")
            println("   BŁĄD: ", first(sprint(showerror, e), 300))
        finally
            rm(tmp; recursive=true, force=true)
        end
    end

    cols = ["pulsar", "wariant", "N", "n_rfi_auto", "n_null_auto", "frac_pairs_pct", "n_subints_full", "n_subints_any",
            "chans_full", "jac_full", "prec_full", "rec_full", "jac_any", "prec_any", "rec_any",
            "n_subints_half", "jac_half", "prec_half", "rec_half", "null_in_full", "null_in_half", "null_in_any", "snr_mg_half", "top3_chans", "zap_chans16", "chans16_in_top3", "snr_none", "snr_auto", "snr_mg_full",
            "snr_mg_any", "snr_pairs_none", "snr_mg_pairs", "only_auto", "only_mg", "input_unchanged",
            "chan_frac", "only_auto_list", "only_mg_list"]
    f2(x) = x isa AbstractFloat ? (isfinite(x) ? @sprintf("%.3f", x) : "") : string(x)
    open(joinpath(MGROOT, "summary.csv"), "w") do io
        println(io, join(cols, ","))
        for r in rows
            println(io, join([occursin(" ", f2(get(r, c, ""))) ? "\"$(f2(r[c]))\"" : f2(get(r, c, "")) for c in cols], ","))
        end
    end
    write_report(meta, rows, notes)
    println("\nGotowe: $(joinpath(MGROOT, "summary.csv")), $(joinpath(MGROOT, "RAPORT.md")), $(MGROOT)/<psr>.pdf")
end

function write_report(meta, rows, notes)
    g(r, k) = get(r, k, NaN)
    pct(a, b) = 100 * (a / b - 1)
    io = IOBuffer()
    p(s="") = println(io, s)
    p("# Test MeerGuard na 7 pulsarach — porównanie z automatem")
    p()
    p("Część automatyczna (liczby z `summary.csv`). Wnioski i rekomendacja: sekcja na końcu, do uzupełnienia po obejrzeniu PDF-ów.")
    p()
    p("## Środowisko")
    p("- MeerGuard: commit `$(get(meta, "meerguard_commit", "?"))`, skrypt `$(get(meta, "clean_archive", "?"))`")
    p("- Python $(get(meta, "python", "?")), numpy $(get(meta, "numpy", "?")), scipy $(get(meta, "scipy", "?")), psrchive: `$(get(meta, "psrchive", "?"))`")
    miss = get(meta, "help_missing_options", [])
    p("- Opcje z opisu w `clean_archive.py -h`: " * (isempty(miss) ? "wszystkie są" : "BRAK $(join(miss, ", "))"))
    unch = [r["input_unchanged"] for r in rows if r["wariant"] == "A"]
    p("- Pliki wejściowe po MeerGuard bez zmian (md5): " * (all(==(true), unch) ? "tak, wszystkie" : "NIE dla części — sprawdzić"))
    p("- Warianty: A `-C mg16.cfg -c 7 -s 7 -bc 0.95 -bs 0.95`; B domyślne; C `-C mg16.cfg -ag`")
    isempty(notes) || (p(); p("Pominięte / błędy:"); foreach(n -> p("- $n"), notes))
    p()
    p("## Wyniki")
    p()
    p("Impuls „zły według MeerGuard” = oflagowana co najmniej połowa jego kanałów (bez kanałów wyciętych przez MeerGuard w całości). „Z flagą” = choć jedna para.")
    p()
    p("| pulsar | war. | pary % | kanały całe | impulsy pełne / ≥połowa / z flagą | Jaccard z RFI automatu (≥połowa) | nulle trafione: pełne / ≥połowa / z flagą (z n_null) | top-3 kanały | zap_chans16 | S/N: bez / automat / MG ≥połowa / MG pary (bez wycięć w tej sumie) |")
    p("|---|---|---|---|---|---|---|---|---|---|")
    for r in rows
        p(@sprintf("| %s | %s | %.2f | %s | %d / %d / %d | %.2f | %d / %d / %d (z %d) | %s | %s | %.1f / %.1f / %.1f / %.1f (%.1f) |",
                   r["pulsar"], r["wariant"], r["frac_pairs_pct"], isempty(r["chans_full"]) ? "-" : r["chans_full"],
                   r["n_subints_full"], r["n_subints_half"], r["n_subints_any"], g(r, "jac_half"),
                   r["null_in_full"], r["null_in_half"], r["null_in_any"], r["n_null_auto"],
                   r["top3_chans"], isempty(r["zap_chans16"]) ? "-" : r["zap_chans16"],
                   r["snr_none"], r["snr_auto"], r["snr_mg_half"], r["snr_mg_pairs"], r["snr_pairs_none"]))
    end
    p()
    p("## Odpowiedzi (z liczb)")
    A = filter(r -> r["wariant"] == "A", rows)
    B = Dict(r["pulsar"] => r for r in rows if r["wariant"] == "B")
    p("1. **Czy MeerGuard rusza nulle?** Wariant A: " * join([@sprintf("%s: %d/%d pełnych, %d ≥połowa kanałów, %d z jakąkolwiek flagą",
            r["pulsar"], r["null_in_full"], r["n_null_auto"], r["null_in_half"], r["null_in_any"]) for r in A if r["n_null_auto"] > 0], "; ") *
      ". Wariant B (domyślny): " * join([@sprintf("%s: %d ≥połowa, %d z flagą", r["pulsar"], r["null_in_half"], r["null_in_any"])
            for r in rows if r["wariant"] == "B" && r["n_null_auto"] > 0], "; ") * ".")
    p("2. **Te same impulsy RFI?** Jaccard impulsów ≥połowa kanałów z `zaps2_rfi` (A): " *
      join([@sprintf("%s %.2f (precyzja %.2f, czułość %.2f)", r["pulsar"], g(r, "jac_half"), g(r, "prec_half"), g(r, "rec_half")) for r in A], ", ") *
      ". Kanały: zap_chans16 w top-3 kanałów MeerGuard: " *
      join([string(r["pulsar"], " ", r["chans16_in_top3"], "/", length(split(r["zap_chans16"]))) for r in A if !isempty(r["zap_chans16"])], ", ") * ".")
    p("3. **Domyślne ustawienia (B) wobec A, odsetek oflagowanych par:** " *
      join([@sprintf("%s %.2f%% vs %.2f%%", r["pulsar"], B[r["pulsar"]]["frac_pairs_pct"], r["frac_pairs_pct"]) for r in A if haskey(B, r["pulsar"])], ", ") * ".")
    p("4. **Coś, czego automat nie wyciął?** Impulsy tylko MeerGuard (A, z jakąkolwiek flagą): " *
      join([string(r["pulsar"], " ", r["only_mg"]) for r in A], ", ") * ". Czy to RFI — strona 3 PDF-ów.")
    j = filter(r -> r["pulsar"] == "J1057-5226", A)
    p("5. **Interpuls J1057−5226:** " * (isempty(j) ? "brak wyniku" :
      @sprintf("pary %.2f%%, subinty pełne %d, S/N bez %.1f / automat %.1f / MG pełne %.1f.", j[1]["frac_pairs_pct"],
               j[1]["n_subints_full"], j[1]["snr_none"], j[1]["snr_auto"], j[1]["snr_mg_full"])))
    p()
    p("## Rozbieżności z opisem działania (z kodu, commit wyżej)")
    p("- `nthreads` nie jest opcją `clean_archive.py`; jest tylko w `default.cfg`, a dopasowanie w surgical i tak liczy w 1 wątku.")
    p("- Z `-C` i bez `-c`/`-s`/`-bc`/`-bs` progi biorą się z pliku konfiguracji; w `default.cfg` to 3σ i 0,8, nie 7σ i 0,95. Dlatego w A progi są podane jawnie.")
    p("- Z `-ag` progi 5 / 0,8 obowiązują także przy `-C`; `-ag` zmienia też łączenie statystyk ze średniej na maksimum.")
    p("- Statystyk są 4 (komentarz w kodzie mówi o 5); wyzerowanie przy wyniku ≥ 1.")
    p("- `bandwagon`: ostra nierówność (>); mianownik bez kanałów/subintów już w pełni wyzerowanych.")
    p("- Szablon: jeśli ma mniej kanałów niż dane (i więcej niż 1), MeerGuard go f-scrunchuje; szablonu z wieloma subintami nie scrunchuje sam (błąd).")
    p("- Nazwa wyjścia z `-e` w tym samym katalogu co wejście z tym samym rozszerzeniem nadpisałaby wejście (kod tego nie sprawdza); tu wyjście idzie do osobnych katalogów.")
    p()
    p("## Wnioski i rekomendacja")
    p("_Do uzupełnienia po obejrzeniu `<psr>.pdf`._")
    write(joinpath(MGROOT, "RAPORT.md"), String(take!(io)))
end

mg_compare()
