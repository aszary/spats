# LRFS z poprawionym panelem fazy + sliding LRFS (okno 256 impulsów, krok 1) dla wybranych pulsarów.
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/sliding_lrfs.jl"); sliding_lrfs_run()
# Dane: surowe impulsy (kanały jak w pazur2), wycięte TYLKO RFI (zaps2_rfi), nulle zostają (jak w publikacjach).
# Dla każdego pulsara w /home/psr/data/OUTPUT/sliding_lrfs/<psr>/:
#   <psr>_lrfs.pdf            LRFS z całej obserwacji, panel fazy: istotne biny, rozwinięta faza, błędy
#   <psr>_psrsalsa_lrfs.pdf   to samo dla pliku pspec z pazur2 (jeśli jest) - faza liczona z impulsów
#   <psr>_sliding_lrfs.pdf    sliding LRFS: widmo vs numer impulsu + faza przy P3 w każdym oknie
#   <psr>_lrfs_tracks.pdf     amplituda i faza przy fp w każdym oknie (jak Basu & Mitra 2018, Fig. 2)
# Na koniec wszystko w jednym pliku: sliding_lrfs/wszystkie.pdf (pdfunite).
isdefined(Main, :F7_DIR) || include(joinpath(@__DIR__, "final_v7.jl"))
SL_DIR = F6_ROOT * "sliding_lrfs/"
# cztery pulsary z Basu & Mitra (2018) do porównania + znane dryfujące / z nullami / zmianami P3
SL_PSRS = ["J0820-1350", "J1034-3224", "J1555-3134", "J1720-2933",
           "J0034-0721", "J1750-3503", "J1727-2739", "J2007+0910", "J1539-6322", "J2053-7200",
           "J1944+1755", "J0304+1932"]

"Impulsy z wyciętym tylko RFI (wiersze zerowe) i parametry z pazur2."
function _sl_data(psr)
    p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
    D = Float64.(_f6_load(psr, p))
    for (a, b) in SpaTs.Data.zaps_to_ranges(SpaTs.Data.zaps_rfi_of(p))
        D[max(1, a + 1):min(size(D, 1), b + 1), :] .= 0.0
    end
    return D, p
end

function sliding_lrfs_run(psrs=SL_PSRS; window=256, step=1, outdir=SL_DIR)
    isdefined(Main, :_debased_spcf) || Base.include(Main, "modules/zap_verify.jl")
    PL = SpaTs.Plot
    pdfs = String[]; summary = String[]
    for psr in psrs
        t0 = time()
        try
            D, p = _sl_data(psr)
            od = joinpath(outdir, psr); mkpath(od)
            N = size(D, 1)
            # LRFS z całej obserwacji (układ jak plik pspec: wiersze 0..0.5 cykla/P, pol 1/2 = Re/Im)
            Y, _ = PL._demeaned(D)
            F = rfft(Y, 1)
            arr = cat(real.(F), imag.(F); dims=3)
            Base.invokelatest(PL.lrfs, arr, od, p; name_mod=psr, pulses=D)
            push!(pdfs, joinpath(od, "$(psr)_lrfs.pdf"))
            sf = joinpath(F6_ROOT, "pazur2", psr, "pulsar.debase.lrfs")
            if isfile(sf)
                d3 = SpaTs.Data.load_ascii_all(sf)
                Base.invokelatest(PL.lrfs, d3, od, p; name_mod="$(psr)_psrsalsa", pulses=D)
                push!(pdfs, joinpath(od, "$(psr)_psrsalsa_lrfs.pdf"))
            end
            if N > window
                res = Base.invokelatest(PL.sliding_lrfs, D, od, p; window=window, step=step, name_mod=psr)
                push!(pdfs, joinpath(od, "$(psr)_sliding_lrfs.pdf"), joinpath(od, "$(psr)_lrfs_tracks.pdf"))
                push!(summary, "$psr  N=$N  " * res.info)
            else
                push!(summary, "$psr  N=$N  za mało impulsów na okno $window")
            end
            @printf("%-12s ok (%.0f s)\n", psr, time() - t0)
        catch e
            e isa InterruptException && rethrow()
            println("$psr BŁĄD: ", first(sprint(showerror, e), 300))
            push!(summary, "$psr  BŁĄD")
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
