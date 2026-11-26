# Porownanie naszych wyciec z lista wycietych impulsow TPA VIII (Song+2023, plik stats.csv w katalogu repo:
# kolumny NrZappenPulses, ZappedPulsesList). TPA wycinalo TYLKO pojedyncze impulsy z RFI (AutoClean), nulli nie -
# wiec porownujemy glownie z naszym RFI (zaps2_rfi); osobno: czy ich impuls wycielismy jako null.
# Numeracja: oba od 0 (sprawdzane takze przesuniecie o 1 - kolumna "przes.+1").
# Dla impulsow TPA, ktorych nie wycielismy: nasza miara RFI (szum poza impulsem, z wzgledem ~31 sasiadow).
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/tpa_zap_compare.jl"); tpa_zap_compare()
# Wynik: /home/psr/data/OUTPUT/tpa_zap/WYNIK.txt (+ tabela.csv)
isdefined(Main, :_f6_load) || include(joinpath(@__DIR__, "final_v6.jl"))
using Printf, Statistics
TZ_DIR = F6_ROOT * "tpa_zap/"
TZ_STATS = joinpath(@__DIR__, "..", "stats.csv")
# pulsary do szczegolowego wypisania (wszystkie wspolne z wycieciami TPA i tak trafiaja do tabeli)
TZ_DETAIL = ["J1931+1439", "J2136-1606", "J1843+2024"]

"CSV z polami w cudzyslowach obejmujacymi nowe linie (stats.csv TPA)."
function _tz_parsecsv(s)
    rows = Vector{Vector{String}}(); row = String[]; cur = IOBuffer(); q = false
    for c in s
        if c == '"'; q = !q
        elseif c == ',' && !q; push!(row, String(take!(cur)))
        elseif c == '\n' && !q; push!(row, String(take!(cur))); push!(rows, row); row = String[]
        elseif c == '\r'
        else; write(cur, c); end
    end
    isempty(row) || push!(rows, row)
    rows
end

"TPA: nazwa => (N, wyciete impulsy (od 0), AutoClean, obserwacja); przy kilku obserwacjach - ta z cenzusu."
function _tz_tpa(file=TZ_STATS)
    R = _tz_parsecsv(read(file, String)); H = R[1]
    c(n) = findfirst(==(n), H)
    iJ, iC, iO, iN, iA, iZ, iL = c("JName"), c("Census"), c("Obsname"), c("NrPulsesTot"), c("AutoClean"), c("NrZappenPulses"), c("ZappedPulsesList")
    T = Dict{String,Any}()
    for v in R[3:end]
        length(v) == length(H) || continue
        haskey(T, v[iJ]) && v[iC] != "YES" && continue
        z = [x for x in (tryparse(Int, w) for w in split(v[iL])) if !isnothing(x)]
        T[v[iJ]] = (N=something(tryparse(Int, v[iN]), -1), zaps=z, auto=v[iA], obs=v[iO], nz=something(tryparse(Int, v[iZ]), 0))
    end
    T
end

_tz_rows0(zaps, N) = isnothing(zaps) || zaps == "" ? Int[] : SpaTs.Data.zaps_to_rows(zaps, N) .- 1     # od 0

function tpa_zap_compare(; out=TZ_DIR, detail=TZ_DETAIL)
    _VD = SpaTs.Data
    mkpath(out)
    T = _tz_tpa()
    ours = Set(_f6_psrs())
    psrs = sort([q for (q, t) in T if q in ours && t.nz > 0])
    L = String[]; pr(s) = (println(s); push!(L, s))
    pr("# Wycięte impulsy: TPA VIII (stats.csv) a nasze (pazur2)\n")
    pr("TPA wycina tylko pojedyncze impulsy z RFI; pulsary wspólne z wycięciami TPA: $(length(psrs)) " *
       "(TPA ma wycięcia w $(count(t -> t.nz > 0, values(T))) z $(length(T)); naszych pulsarów $(length(ours)))\n")
    pr("| pulsar | N TPA / nasze | TPA wycięło | my RFI | wspólne z naszym RFI | wycięte u nas jako null | u nas niewycięte | przes.+1 | nasze RFI spoza TPA | z (RFI) niewyciętych: mediana / max |")
    pr("|---|---|---|---|---|---|---|---|---|---|")
    csv = ["pulsar,N_tpa,N,tpa,ours_rfi,common_rfi,as_null,not_cut,shift1,ours_only,z_med,z_max"]
    tot = zeros(Int, 6)
    for psr in psrs
        try
            t = T[psr]
            p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
            D = _f6_load(psr, p); N = size(D, 1)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            rfi = Set(_tz_rows0(_VD.zaps_rfi_of(p), N)); nul = Set(_tz_rows0(_VD.zaps_null_of(p), N))
            allc = Set(_tz_rows0(_VD.zaps_of(p), N))
            tz = Set(t.zaps)
            com = length(intersect(tz, rfi)); asn = length(intersect(setdiff(tz, rfi), nul))
            notc = sort(collect(setdiff(tz, allc)))
            sh1 = length(intersect(Set(t.zaps .+ 1), rfi))
            only = length(setdiff(rfi, tz))
            z = _VD.rfi_rms_z(D, bst, bend)
            zn = [z[i + 1] for i in notc if 0 <= i < N]
            pr(@sprintf("| %s | %d / %d | %d | %d | %d | %d | %d | %d | %d | %s |", psr, t.N, N, length(tz), length(rfi), com, asn,
                        length(notc), sh1, only, isempty(zn) ? "-" : @sprintf("%.1f / %.1f", median(zn), maximum(zn))))
            push!(csv, join([psr, t.N, N, length(tz), length(rfi), com, asn, length(notc), sh1, only,
                             isempty(zn) ? "" : @sprintf("%.2f", median(zn)), isempty(zn) ? "" : @sprintf("%.2f", maximum(zn))], ","))
            tot .+= [length(tz), length(rfi), com, asn, length(notc), only]
            if psr in detail
                push!(L, "    $psr: TPA $(sort(collect(tz)))")
                push!(L, "    $psr: nasze RFI w pobliżu impulsów TPA (±3): " *
                         string(sort([i for i in rfi if any(abs(i - j) <= 3 for j in tz)])))
                push!(L, "    $psr: niewycięte przez nas, z(RFI): " * join([@sprintf("%d (z=%.1f)", i, z[i + 1]) for i in notc if 0 <= i < N], ", "))
            end
        catch e
            e isa InterruptException && rethrow()
            pr("| $psr | BŁĄD: $(first(sprint(showerror, e), 120)) |")
        end
    end
    pr("\nRazem: TPA wycięło $(tot[1]) impulsów; z nich nasze RFI $(tot[3]) ($(round(Int, 100tot[3] / max(tot[1], 1)))%), " *
       "u nas jako null $(tot[4]), u nas niewycięte $(tot[5]). Nasze RFI w tych pulsarach: $(tot[2]), z tego spoza listy TPA $(tot[6]).")
    pr("z(RFI) = odchylenie szumu poza impulsem od mediany sąsiadów w σ; nasza reguła RFI tnie od z > 3.")
    write(out * "WYNIK.txt", join(L, "\n") * "\n"); write(out * "tabela.csv", join(csv, "\n") * "\n")
    println("\nZapisane: $(out)WYNIK.txt, tabela.csv")
end
