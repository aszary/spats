# Automat v5 = v4 + nowy bezpiecznik RFI (Data.auto_zap(rfi_greedy=true)): gdy caly zestaw
# kandydatow RFI (> 3 sigma) obniza S/N o > 1%, stary bezpiecznik odrzucal WSZYSTKICH (podnosil prog);
# nowy dokłada kandydatow od tych, ktore najmniej szkodza S/N, i bierze tylko te, ktore S/N nie obnizaja.
# Gdzie stary nie podnosil progu - wynik identyczny. Test syntetyczny: falszywe 1.0% / 1.0%,
# S/N 100.6% / 100.8% ideału; J2136-like: RFI wyciete 7% -> 15%.
#
# DWIE WERSJE: pazur2 (v4) NIE jest zmieniany. Uruchom w Julii na kopierniku:
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/rfi_v5.jl"); rfi_v5()
#     1. skan wszystkich pulsarow (16 procesow, nic nie zapisuje w pazur2): v4 vs v5 -> lista zmienionych
#     2. dla zmienionych: pelny automat v5 -> /home/psr/data/OUTPUT/pazur2_v5/ (te same okna on-pulse i kanaly)
#     3. weryfikacja v5 (zap_verify_v5/) + porownanie z metoda reczna (zap_comparison_v5/) + tabela v4 vs v5
#     4. diagnostyka J2111+2106 (zap_checks_j2111/)
#   Jesli v5 lepsze - DOPIERO WTEDY:  rfi_v5_adopt()   (pazur2/<psr> -> pazur2_v4_backup/<psr>, v5 -> pazur2,
#   energie i model nulli dla tych pulsarow od nowa, porownanie)
using Printf, Statistics

V5_ROOT = "/home/psr/data/OUTPUT/"
V5_HEADER = "pulsar,N,n_old,n_new,snr_old,snr_new,snr_rfi_old,snr_rfi_new,n_rfi_old,n_rfi_new,snr_none,zaps_old,zaps_new"

function _v5_psrs(root=V5_ROOT)
    st = Dict{String,String}()
    for l in Iterators.drop(eachline(root * "pazur2/zap_summary.csv"), 1)
        r = split(l, ","); length(r) >= 2 && (st[r[1]] = r[2])
    end
    sort([k for (k, v) in st if v in ("ok", "reczny")])
end

"Część k/n skanu (w procesie potomnym). Liczy automat v4 i v5 na tych samych danych; nic nie zapisuje w pazur2."
function _v5_scan_part(k, n; root=V5_ROOT, dataroot="/home/psr/data/new/")
    out = root * "rfi_v5_scan/"; mkpath(out)
    f = out * "scan_$k.csv"
    done = isfile(f) ? Set(split(l, ",")[1] for l in Iterators.drop(eachline(f), 1)) : Set{String}()
    isfile(f) || open(io -> println(io, V5_HEADER), f, "w")
    psrs = filter(!in(done), _v5_psrs(root)[k:n:end])
    println("[0/$(length(psrs))] skan v4/v5, część $k/$n"); flush(stdout)
    D_ = SpaTs.Data
    for (j, psr) in enumerate(psrs)
        t0 = time()
        try
            p = SpaTs.Tools.read_params(joinpath(root, "pazur2", psr, "params.json"))
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            obs = sort(filter(d -> isdir(joinpath(dataroot, psr, d)), readdir(joinpath(dataroot, psr))))
            tmp = mktempdir()
            D = try
                Main._debased_spcf(joinpath(dataroot, psr, obs[1]) * "/", tmp, bst, bend;
                                   chans16=something(get(p, "zap_chans16", nothing), ""))
            finally
                rm(tmp; recursive=true, force=true)
            end
            a = D_.auto_zap(D, bst, bend; rfi_greedy=false)
            b = D_.auto_zap(D, bst, bend; rfi_greedy=true)
            ra, rb = D_.split_zaps(a.reason).rfi, D_.split_zaps(b.reason).rfi
            s(x) = D_.profile_snr(D, bst, bend; drop=x)
            row = [psr, size(D, 1), length(a.bad), length(b.bad), s(a.bad), s(b.bad), s(ra), s(rb), length(ra), length(rb),
                   s(Int[]), "\"$(D_.pulses_to_zaps(a.bad))\"", "\"$(D_.pulses_to_zaps(b.bad))\""]
            open(io -> println(io, join([x isa AbstractFloat ? @sprintf("%.3f", x) : string(x) for x in row], ",")), f, "a")
            @printf("[%d/%d] %-12s v4 %d cięć S/N %.1f | v5 %d cięć S/N %.1f  (%.0f s)\n", j, length(psrs), psr,
                    length(a.bad), s(a.bad), length(b.bad), s(b.bad), time() - t0)
        catch e
            e isa InterruptException && rethrow()
            println("[$j/$(length(psrs))] $psr BŁĄD: ", first(sprint(showerror, e), 200))
        end
        flush(stdout)
    end
end

function _v5_scan_read(root=V5_ROOT)
    R = Dict{String,Dict{String,String}}()
    H = split(V5_HEADER, ",")
    for f in filter(f -> occursin(r"^scan_\d+\.csv$", basename(f)), readdir(root * "rfi_v5_scan/"; join=true))
        for l in Iterators.drop(eachline(f), 1)
            # dwie ostatnie kolumny w cudzyslowach (zakresy ze spacjami, bez przecinkow)
            v = split(l, ","; limit=length(H))
            R[v[1]] = Dict(zip(H, v))
        end
    end
    R
end

function rfi_v5_scan_summary(; root=V5_ROOT)
    R = _v5_scan_read(root)
    num(r, k) = parse(Float64, r[k])
    ch = sort([r for r in values(R) if r["zaps_old"] != r["zaps_new"]], by=r -> -(num(r, "snr_new") / num(r, "snr_old")))
    println("\n===== SKAN v4 vs v5 (nowy bezpiecznik RFI): $(length(R)) pulsarów =====")
    println("zmienione wycięcia: $(length(ch)) pulsarów;  identyczne: $(length(R) - length(ch))")
    isempty(ch) && return String[]
    g = [100 * (num(r, "snr_new") / num(r, "snr_old") - 1) for r in ch]
    gr = [100 * (num(r, "snr_rfi_new") / num(r, "snr_rfi_old") - 1) for r in ch]
    @printf("zmiana S/N v5/v4 w zmienionych: mediana %+.2f%%, maks. %+.2f%%, min. %+.2f%%;  >1%%: %d,  <-1%%: %d\n",
            median(g), maximum(g), minimum(g), count(>(1), g), count(<(-1), g))
    @printf("dane tylko po RFI (do nullingu): mediana %+.2f%%, min. %+.2f%%\n", median(gr), minimum(gr))
    println(rpad("pulsar", 13), lpad("imp.", 6), lpad("cięć v4", 9), lpad("v5", 6), lpad("S/N v4", 10), lpad("v5", 10),
            lpad("zmiana", 9), lpad("RFI v4", 8), lpad("v5", 6), lpad("tylko-RFI", 11))
    for (r, x, y) in zip(ch, g, gr)
        @printf("%-13s%6s%9s%6s%10.1f%10.1f%+8.2f%%%8s%6s%+10.2f%%\n", r["pulsar"], r["N"], r["n_old"], r["n_new"],
                num(r, "snr_old"), num(r, "snr_new"), x, r["n_rfi_old"], r["n_rfi_new"], y)
    end
    return [r["pulsar"] for r in ch]
end

function rfi_v5(; root=V5_ROOT, n=16)
    t0 = time()
    Base.include(Main, "modules/zap_verify.jl")
    println("===== 1. SKAN: automat v4 i v5 na wszystkich pulsarach ($n procesów; pazur2 bez zmian) =====")
    codes = ["include(\"spats.jl\"); Base.include(Main, \"modules/zap_verify.jl\"); include(\"tools/rfi_v5.jl\"); " *
             "Base.invokelatest(_v5_scan_part, $k, $n)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, root * "rfi_v5_scan/logs", "scan"; every=120)
    changed = rfi_v5_scan_summary(root=root)

    if !isempty(changed)
        println("\n===== 2. AUTOMAT v5 dla $(length(changed)) zmienionych → $(root)pazur2_v5/ =====")
        SpaTs.Data.RFI_GREEDY[] = true
        try
            cf = SpaTs.chan_first_list([root * "zap_verify_v3/", root * "zap_verify/"])
            SpaTs.autozap_all(out=root * "pazur2_v5/", only=changed, skip_done=false, keep_manual=false,
                              reuse_window=true, window_root=root * "pazur2/", chan_first=cf)
            println("\n===== 3. WERYFIKACJA v5 → $(root)zap_verify_v5/ =====")
            Base.invokelatest(Main.zap_verify; pazur2=root * "pazur2_v5/", out=root * "zap_verify_v5/", only=changed,
                              save_rfi=true, max_chan_pulses=10^6)
        finally
            SpaTs.Data.RFI_GREEDY[] = false
        end
        SpaTs.compare_zap_methods(mine=root * "pazur2_v5/", outdir=root * "zap_comparison_v5/")
        println("\n===== v4 (pazur2) vs v5 (pazur2_v5), ta sama miara S/N =====")
        println(rpad("pulsar", 13), lpad("S/N bez", 9), lpad("v4", 9), lpad("v5", 9), lpad("zmiana", 9),
                lpad("ręcznie", 9), lpad("nulle v4", 9), lpad("v5", 5))
        for psr in changed
            c4 = get(SpaTs.Tools.read_params(joinpath(root, "pazur2", psr, "params.json")), "zap_cmp", Dict())
            f5 = joinpath(root, "pazur2_v5", psr, "params.json")
            isfile(f5) || (println(rpad(psr, 13), " brak wyniku v5"); continue)
            c5 = get(SpaTs.Tools.read_params(f5), "zap_cmp", Dict())
            r1(x) = x isa Number ? @sprintf("%.1f", x) : "-"
            a, b = get(c4, "snr_mine", NaN), get(c5, "snr_mine", NaN)
            println(rpad(psr, 13), lpad(r1(get(c4, "snr_none", nothing)), 9), lpad(r1(a), 9), lpad(r1(b), 9),
                    lpad(@sprintf("%+.2f%%", 100 * (b / a - 1)), 9), lpad(r1(get(c4, "snr_maciej", nothing)), 9),
                    lpad(string(get(c4, "n_null", "-")), 9), lpad(string(get(c5, "n_null", "-")), 5))
        end
    end

    println("\n===== 4. DIAGNOSTYKA J2111+2106 =====")
    Base.include(Main, "tools/zap_checks.jl")
    Base.invokelatest(Main.zap_checks; verify_dir=root * "zap_verify_v4_manual/", out=root * "zap_checks_j2111/",
                      psrs=["J2111+2106"])
    @printf("\nKONIEC rfi_v5 (%.0f min). pazur2 (v4) bez zmian; v5 w %spazur2_v5/. Wklej: tabelę skanu, tabelę v4 vs v5, blok AUTOMAT vs METODA RĘCZNA (v5) i część o J2111+2106.\n",
            (time() - t0) / 60, root)
end

"""
Szybka wersja kroków 2-4 (skan musi być już zrobiony): automat v5 i weryfikacja równolegle w `n`
procesach; weryfikacja bez części o kanałach (nie zależą od bezpiecznika RFI, sprawdzone w v4).
Wznawialna: pulsary już przeliczone do pazur2_v5 / zap_verify_v5 są pomijane.
"""
function rfi_v5_fast(; root=V5_ROOT, n=32)
    t0 = time()
    Base.include(Main, "modules/zap_verify.jl")
    changed = rfi_v5_scan_summary(root=root)
    isempty(changed) && return
    n = min(n, length(changed))
    chg = repr(changed)
    pre = "include(\"spats.jl\"); SpaTs.Data.RFI_GREEDY[] = true; "
    println("\n===== 2. AUTOMAT v5: $(length(changed)) pulsarów, $n procesów → $(root)pazur2_v5/ =====")
    codes = [pre * "ch = $chg; cf = SpaTs.chan_first_list([\"$(root)zap_verify_v3/\", \"$(root)zap_verify/\"]); " *
             "SpaTs.autozap_all(out=\"$(root)pazur2_v5/\", only=ch[$k:$n:end], skip_done=true, keep_manual=false, " *
             "reuse_window=true, window_root=\"$(root)pazur2/\", chan_first=cf)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, root * "pazur2_v5/logs", "auto"; every=60)
    println("\n===== 3. WERYFIKACJA v5 (bez kanałów): $n procesów → $(root)zap_verify_v5/ =====")
    codes = [pre * "Base.include(Main, \"modules/zap_verify.jl\"); " *
             "Base.invokelatest(Main.zap_verify; pazur2=\"$(root)pazur2_v5/\", out=\"$(root)zap_verify_v5/\", " *
             "only=$chg, part=\"$k/$n\", chans=false, save_rfi=true)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, root * "zap_verify_v5/logs", "verify"; every=60)
    Base.invokelatest(Main.zap_verify_summary; out=root * "zap_verify_v5/")
    _v5_tail(root, changed)
    @printf("\nKONIEC rfi_v5_fast (%.0f min). pazur2 (v4) bez zmian; v5 w %spazur2_v5/. Wklej: summary weryfikacji v5, tabelę v4 vs v5, blok AUTOMAT vs METODA RĘCZNA (v5) i część o J2111+2106.\n",
            (time() - t0) / 60, root)
end

"Porównanie z metodą ręczną dla v5, tabela v4 vs v5 i diagnostyka J2111+2106."
function _v5_tail(root, changed)
    SpaTs.compare_zap_methods(mine=root * "pazur2_v5/", outdir=root * "zap_comparison_v5/")
    println("\n===== v4 (pazur2) vs v5 (pazur2_v5), ta sama miara S/N =====")
    r1(x) = x isa Number ? @sprintf("%.1f", x) : "-"
    d = Float64[]; worse = String[]
    for psr in changed
        c4 = get(SpaTs.Tools.read_params(joinpath(root, "pazur2", psr, "params.json")), "zap_cmp", Dict())
        f5 = joinpath(root, "pazur2_v5", psr, "params.json")
        isfile(f5) || (println(rpad(psr, 13), " brak wyniku v5"); continue)
        c5 = get(SpaTs.Tools.read_params(f5), "zap_cmp", Dict())
        a, b = get(c4, "snr_mine", NaN), get(c5, "snr_mine", NaN)
        x = 100 * (b - a) / abs(a); push!(d, x)
        x < -0.5 && push!(worse, @sprintf("%s (%+.2f%%)", psr, x))
    end
    isempty(d) && return
    @printf("S/N v5 względem v4 w %d pulsarach: mediana %+.2f%%, średnia %+.2f%%, maks. %+.2f%%, min. %+.2f%%; lepiej >1%%: %d, gorzej <-0.5%%: %d\n",
            length(d), median(d), mean(d), maximum(d), minimum(d), count(>(1), d), length(worse))
    isempty(worse) || println("gorzej: ", join(worse, ", "))
    println("\n===== 4. DIAGNOSTYKA J2111+2106 =====")
    Base.include(Main, "tools/zap_checks.jl")
    Base.invokelatest(Main.zap_checks; verify_dir=root * "zap_verify_v4_manual/", out=root * "zap_checks_j2111/",
                      psrs=["J2111+2106"])
end

"Przyjęcie v5: tylko dla pulsarów z pazur2_v5. Kopia v4 w pazur2_v4_backup/. Potem energie i model nulli od nowa."
function rfi_v5_adopt(; root=V5_ROOT)
    v5, pz, bk = root * "pazur2_v5/", root * "pazur2/", root * "pazur2_v4_backup/"
    psrs = sort(filter(d -> isfile(joinpath(v5, d, "params.json")), readdir(v5)))
    isempty(psrs) && return println("Brak wyników v5 w $v5")
    mkpath(bk)
    cp(pz * "zap_summary.csv", bk * "zap_summary.csv"; force=true)
    rows5 = Dict(split(l, ",")[1] => l for l in Iterators.drop(eachline(v5 * "zap_summary.csv"), 1))
    for psr in psrs
        isdir(bk * psr) || cp(pz * psr, bk * psr)            # kopia v4 (raz)
        rm(pz * psr; recursive=true, force=true)
        cp(v5 * psr, pz * psr)
    end
    L = readlines(pz * "zap_summary.csv")
    L = [i == 1 || !haskey(rows5, split(l, ",")[1]) ? l : rows5[split(l, ",")[1]] for (i, l) in enumerate(L)]
    write(pz * "zap_summary.csv", join(L, "\n") * "\n")
    println("v5 przyjęte w $(length(psrs)) pulsarach; kopia v4: $bk")
    nout = root * "nulling_model_v4/"
    for psr in psrs
        rm(joinpath(nout, "energies", "$psr.csv"); force=true)
        rm(joinpath(nout, "fits", "$psr.v3.json"); force=true)
    end
    ENV["NULLING_OUT"] = nout
    Base.include(Main, "tools/nulling_energies.jl")
    Base.invokelatest(Main.nulling_energies_parallel)
    try run(addenv(`python3 tools/nulling_model.py`, "NULLING_OUT" => nout)) catch e; println("model nulli: ", e) end
    SpaTs.compare_zap_methods()
end
