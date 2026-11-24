# LRFS z poprawionym panelem fazy + sliding LRFS (okno 256 impulsów, krok 1) dla wybranych pulsarów.
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/sliding_lrfs.jl"); sliding_lrfs_run()
# Dane: surowe impulsy (kanały jak w pazur2), wycięte TYLKO RFI (zaps2_rfi), nulle zostają (jak w publikacjach).
# Dla każdego pulsara w /home/psr/data/OUTPUT/sliding_lrfs/<psr>/:
# trzy wersje:
#   V1 nasza:       <psr>_v1_nasza_lrfs.pdf (LRFS, poprawiony panel fazy), _v1_nasza_psrsalsa_lrfs.pdf (plik pspec),
#                   _v1_nasza_sliding.pdf (sliding LRFS krok 1: moc/szum, jasność impulsów, faza przy P3 w oknach)
#   V2 Basu&Mitra:  _v2_basu_mitra_sliding.pdf (rys. 2 lewy: okno 256, krok 50, średni szczyt = 1),
#                   _v2_basu_mitra_tracks.pdf (rys. 2 prawy: amplituda i faza przy fp dla okien, 0 w szczycie profilu)
#   V3 TPA XX:      _v3_tpa_lrfs.pdf (rys. 2: FFT 200, DC = 0, max = 1, skala log),
#                   _v3_tpa_phase.pdf (rys. 4: fazy w blokach 64, wyrównane, histogram + średnia, reszta dudnień)
# Na koniec wszystko w jednym pliku: sliding_lrfs/wszystkie.pdf (pdfunite).
isdefined(Main, :F7_DIR) || include(joinpath(@__DIR__, "final_v7.jl"))
# wczytane tutaj (na poziomie pliku), nie wewnatrz funkcji: inaczej "method too new" (world age)
isdefined(Main, :_debased_spcf) || Base.include(Main, joinpath(@__DIR__, "..", "modules", "zap_verify.jl"))
SL_DIR = F6_ROOT * "sliding_lrfs/"
# cztery pulsary z Basu & Mitra (2018) do porównania + znane dryfujące / z nullami / zmianami P3
SL_PSRS = ["J0820-1350", "J1034-3224", "J1555-3134", "J1720-2933",
           "J0034-0721", "J1750-3503", "J1727-2739", "J2007+0910", "J1539-6322", "J2053-7200",
           "J1944+1755", "J0304+1932",
           "J1901-0906", "J1514-4834"]                       # TPA XX (Hsu+2026) i H25: wspolistniejace okresowosci
# znane cechy z publikacji (P3 w okresach) - fazy w blokach 64 impulsow jak w TPA XX; pozostale: szczyty z widma
SL_P3S = Dict("J1901-0906" => [3.05, 5.31, 7.33])          # F1, F2, F3 (TPA XX, Tabela 1, C1); F1 = F2 + F3

"Impulsy z wyciętym tylko RFI (wiersze zerowe) i parametry z pazur2."
function _sl_data(psr)
    p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
    D = Float64.(_f6_load(psr, p))
    for (a, b) in SpaTs.Data.zaps_to_ranges(SpaTs.Data.zaps_rfi_of(p))
        D[max(1, a + 1):min(size(D, 1), b + 1), :] .= 0.0
    end
    return D, p
end

function sliding_lrfs_run(psrs=SL_PSRS; window=256, step=1, basu_step=50, outdir=SL_DIR)
    mkpath(outdir)
    PL = SpaTs.Plot
    pdfs = String[]; summary = String[]
    for psr in psrs
        t0 = time()
        try
            D, p = Base.invokelatest(_sl_data, psr)
            od = joinpath(outdir, psr); mkpath(od)
            N = size(D, 1)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            push!(summary, "$psr  N=$N  P3 (params) = $(get(p, "p3", "brak"))")
            # ---------- V1: nasza wersja - LRFS z poprawionym panelem fazy + sliding LRFS (krok 1)
            Y, _ = PL._demeaned(D)
            F = rfft(Y, 1)
            arr = cat(real.(F), imag.(F); dims=3)
            Base.invokelatest(PL.lrfs, arr, od, p; name_mod="$(psr)_v1_nasza", pulses=D)
            push!(pdfs, joinpath(od, "$(psr)_v1_nasza_lrfs.pdf"))
            sf = joinpath(F6_ROOT, "pazur2", psr, "pulsar.debase.lrfs")
            if isfile(sf)
                d3 = SpaTs.Data.load_ascii_all(sf)
                Base.invokelatest(PL.lrfs, d3, od, p; name_mod="$(psr)_v1_nasza_psrsalsa", pulses=D)
                push!(pdfs, joinpath(od, "$(psr)_v1_nasza_psrsalsa_lrfs.pdf"))
            end
            N > window || (push!(summary, "    za mało impulsów na okno $window"); continue)
            Base.invokelatest(PL.sliding_lrfs, D, od, p; window=window, step=step, name_mod=psr, style=:nasza)
            push!(pdfs, joinpath(od, "$(psr)_v1_nasza_sliding.pdf"))
            # ---------- V2: jak Basu & Mitra (2018) - okno 256, krok jak u nich (50), rys. 2 lewy i prawy
            res = Base.invokelatest(PL.sliding_lrfs, D, od, p; window=window, step=basu_step, name_mod=psr, style=:basu)
            push!(pdfs, joinpath(od, "$(psr)_v2_basu_mitra_sliding.pdf"), joinpath(od, "$(psr)_v2_basu_mitra_tracks.pdf"))
            push!(summary, "    V2 Basu & Mitra: " * res.info)
            # ---------- V3: jak TPA XX (Hsu+2026) - LRFS (FFT 200, rys. 2) + fazy w blokach 64 (rys. 4)
            Base.invokelatest(PL.lrfs_tpa_plot, D, bst, bend, joinpath(od, "$(psr)_v3_tpa_lrfs.pdf");
                              nfft=200, title_="$psr: LRFS jak TPA XX (FFT 200, DC = 0, max = 1, skala log)")
            push!(pdfs, joinpath(od, "$(psr)_v3_tpa_lrfs.pdf"))
            p3s = get(SL_P3S, psr, [1 / t.fp for t in res.tracks])
            bl = Base.invokelatest(PL.lrfs_phase_blocks, D, bst, bend, p3s)
            if bl.nblocks >= 3
                dl = 360 * (bend - bst + 1) / size(D, 2)
                lon = collect(range(-dl / 2, dl / 2, length=bend - bst + 1)); lon .-= lon[bl.jref]
                tl = "$psr: jak TPA XX rys. 4 - bloki $(bl.block) impulsów, użyte $(bl.nblocks) (S/N ≥ 7 w ≥ 60/64)\n" *
                     "faza: konwencja TPA XX (rosnąca = dryf ku późniejszym długościom)" *
                     (length(p3s) == 3 ? (bl.beat_ok ? ", dół: reszta dudnień φ1−(φ2+φ3)" : ", biny k1 ≠ k2+k3 - bez reszty dudnień") : "")
                Base.invokelatest(PL.lrfs_phase_blocks_plot, bl, lon, joinpath(od, "$(psr)_v3_tpa_phase.pdf"); title_=tl)
                push!(pdfs, joinpath(od, "$(psr)_v3_tpa_phase.pdf"))
                bs = isnothing(bl.beat) ? "" : @sprintf(", reszta dudnień: mediana |φ1−(φ2+φ3)| = %.0f° (losowo 90°)", median(abs.(filter(isfinite, bl.beat))))
                push!(summary, "    V3 TPA XX: bloki 64: $(bl.nblocks) użytych, cechy P3 = " *
                               join([@sprintf("%.2f", 1 / f.f) for f in bl.feats], ", ") * bs)
            else
                push!(summary, "    V3 TPA XX: za mało bloków 64 impulsów z S/N ≥ 7 w ≥ 60/64 ($(bl.nblocks))")
            end
            @printf("%-12s ok (%.0f s)\n", psr, time() - t0)
        catch e
            e isa InterruptException && rethrow()
            println("$psr BŁĄD: ", first(sprint(showerror, e), 300))
            push!(summary, "$psr  BŁĄD: " * first(sprint(showerror, e), 200))
        end
        flush(stdout)
    end
    write(joinpath(outdir, "podsumowanie.txt"), join(summary, "\n") * "\n")
    println("\n", join(summary, "\n"))
    all_ = joinpath(outdir, "wszystkie.pdf")
    ex = filter(isfile, pdfs)
    if !isempty(ex) && Sys.which("pdfunite") !== nothing
        run(`pdfunite $ex $all_`); println("\nWszystko w jednym pliku: $all_")
    end
end
