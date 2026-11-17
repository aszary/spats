# Diagnostyka otwartych spraw po weryfikacji automatu (tylko odczyt, nic nie zmienia):
#  1. pulsary, w ktorych metoda reczna jest lepsza (J1843+2024, J1632-4621): ktore impulsy
#     wycieto tylko recznie, czym sie roznia (sygnal, szum poza impulsem wzgledem sasiadow)
#     i czy dociecie ich przez automat podnioslo by S/N (kazdy osobno i razem);
#  2. pulsary, ktore traca > 2% S/N na danych tylko po RFI (z verify*.csv);
#  3. pulsary bez sprawdzenia kanalow (z verify*.csv, z powodem).
# Wynik: /home/psr/data/OUTPUT/zap_checks/checks.txt (+ wydruk).
#
#     cd ~/software/spats && SPATS_NO_MAIN=1 julia --project=. tools/zap_checks.jl
using Statistics, Printf
if !isdefined(Main, :SpaTs)
    ENV["SPATS_NO_MAIN"] = "1"
    include(joinpath(@__DIR__, "..", "spats.jl"))
end
isdefined(Main, :_debased_spcf) || Base.include(Main, joinpath(@__DIR__, "..", "modules", "zap_verify.jl"))

"rms poza impulsem kazdego impulsu jako odporne z wzgledem mediany `win` sasiadow (jak regula RFI automatu)."
function _rms_local_z(D, off; win=31)
    N = size(D, 1)
    r = [std(@view D[i, off]) for i in 1:N]
    hw = win ÷ 2
    res = [r[t] - median(@view r[max(1, t - hw):min(N, t + hw)]) for t in 1:N]
    sd = 1.4826 * median(abs.(res .- median(res)))
    return sd > 0 ? res ./ sd : zeros(N)
end

function _manual_wins(io, psr; pazur2="/home/psr/data/OUTPUT/pazur2/", pazur="/home/psr/data/OUTPUT/pazur/",
                      dataroot="/home/psr/data/new/")
    p = _VT.read_params(joinpath(pazur2, psr, "params.json"))
    bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
    obs = sort(filter(d -> isdir(joinpath(dataroot, psr, d)), readdir(joinpath(dataroot, psr))))
    tmp = mktempdir()
    try
        D = _debased_spcf(joinpath(dataroot, psr, obs[1]) * "/", tmp, bst, bend;
                          chans16=something(get(p, "zap_chans16", nothing), ""))
        N = size(D, 1)
        auto = _VD.zaps_to_rows(_VD.zaps_of(p), N)
        rp = _VT.read_params(joinpath(pazur, psr, "params.json"))
        ref = _VD.zapped_rows(joinpath(pazur, psr, "pulsar.spCF"), bst, bend)
        (isnothing(ref) || ref.N != N) && (println(io, "$psr: brak archiwum metody ręcznej albo inna liczba impulsów"); return)
        manual = sort(unique(vcat(ref.rows, _VD.zaps_to_rows(get(rp, "zaps", nothing), N))))
        snr(rows) = _VD.profile_snr(D, bst, bend; drop=rows)
        s0, sa, sm = snr(Int[]), snr(auto), snr(manual)
        mo = setdiff(manual, auto)
        ao = setdiff(auto, manual)
        sn = _VD.pulse_signal_noise(D, bst, bend)
        zr = _rms_local_z(D, _VD.offpulse_bins(D, bst, bend))
        gain = [100 * (snr(vcat(auto, i)) / sa - 1) for i in mo]
        pos = mo[gain .> 0]
        println(io, "\n== $psr: N=$N, wycięte automat $(length(auto)), ręcznie $(length(manual)), wspólne $(length(intersect(auto, manual)))")
        println(io, @sprintf("S/N: bez %.2f | automat %.2f | ręcznie %.2f | automat + wszystkie tylko-ręczne %.2f | automat + tylko-ręczne z zyskiem (%d) %.2f | ręcznie + tylko-automat %.2f",
                             s0, sa, sm, snr(vcat(auto, mo)), length(pos), snr(vcat(auto, pos)), snr(vcat(manual, ao))))
        println(io, "impulsy wycięte TYLKO ręcznie (numeracja od 0): sygnał/szum impulsu, rms poza impulsem (z wzgl. 31 sąsiadów; automat tnie > 3), zmiana S/N po dołożeniu do automatu:")
        for (i, g) in sort(collect(zip(mo, gain)), by=x -> -x[2])[1:min(40, end)]
            println(io, @sprintf("   %5d   S/N impulsu %6.2f   rms z %6.2f   %+6.3f%%", i - 1, sn.s[i] / sqrt(max(sn.v[i], eps())), zr[i], g))
        end
        if !isempty(mo)
            println(io, @sprintf("   mediana rms z tylko-ręcznych %.2f (wszystkich %.2f); z rms z > 3: %d z %d; z zyskiem S/N: %d",
                                 median(zr[mo]), median(zr), count(>(3), zr[mo]), length(mo), length(pos)))
        end
        println(io, "wycięte tylko przez automat: $(length(ao)) — S/N ręcznie + one: $(@sprintf("%.2f", snr(vcat(manual, ao)))) (wobec ręcznie $(@sprintf("%.2f", sm)))")
    finally
        rm(tmp; recursive=true, force=true)
    end
end

function zap_checks(; verify_dir="/home/psr/data/OUTPUT/zap_verify_v4/", out="/home/psr/data/OUTPUT/zap_checks/",
                    psrs=["J1843+2024", "J1632-4621"])
    mkpath(out)
    io = IOBuffer()
    println(io, "===== DIAGNOSTYKA OTWARTYCH SPRAW =====")
    println(io, "\n--- 1. Metoda ręczna lepsza niż automat ---")
    for psr in psrs
        try
            _manual_wins(io, psr)
        catch e
            e isa InterruptException && rethrow()
            println(io, "$psr: błąd — $(first(sprint(showerror, e), 200))")
        end
    end
    H = split(VERIFY_HEADER, ",")
    R = Dict{SubString{String},SubString{String}}[]
    for f in _verify_csvs(verify_dir), l in Iterators.drop(eachline(f), 1)
        push!(R, Dict(zip(H, split(l, ","; limit=length(H)))))
    end
    num(r, k) = (x = tryparse(Float64, get(r, k, "")); isnothing(x) ? NaN : x)
    println(io, "\n--- 2. Dane tylko po RFI tracą > 2% S/N (z $verify_dir, $(length(R)) pulsarów) ---")
    for r in R
        d = 100 * (num(r, "new_snr_rfi") / num(r, "snr_none") - 1)
        isfinite(d) && d < -2 && println(io, @sprintf("   %-12s %+.1f%%   wycięte RFI %s, nulle %s, S/N bez %.1f, tylko po RFI %.1f",
                                                     r["pulsar"], d, r["n_rfi"], r["n_null"], num(r, "snr_none"), num(r, "new_snr_rfi")))
    end
    println(io, "\n--- 3. Pulsary bez sprawdzenia kanałów ---")
    for r in R
        isfinite(num(r, "C")) || println(io, "   $(r["pulsar"]): $(isempty(r["note"]) ? "(brak notatki)" : r["note"])")
    end
    s = String(take!(io))
    write(joinpath(out, "checks.txt"), s)
    println(s)
    println("Zapisane: $(joinpath(out, "checks.txt"))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    zap_checks()
end
