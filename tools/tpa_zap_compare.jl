# Porownanie naszych wyciec z lista wycietych impulsow TPA VIII (Song+2023, plik stats.csv w katalogu repo:
# kolumny NrZappenPulses, ZappedPulsesList). TPA wycinalo TYLKO pojedyncze impulsy z RFI (AutoClean), nulli nie -
# wiec porownujemy glownie z naszym RFI (zaps2_rfi); osobno: czy ich impuls wycielismy jako null.
# Numeracja: oba od 0 (sprawdzane takze przesuniecie o 1 - kolumna "przes.+1").
# Dla impulsow TPA, ktorych nie wycielismy: nasza miara RFI (szum poza impulsem, z wzgledem ~31 sasiadow).
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/tpa_zap_compare.jl"); tpa_zap_compare()
# Wynik: /home/psr/data/OUTPUT/tpa_zap/WYNIK.txt (+ tabela.csv)
isdefined(Main, :_f6_load) || include(joinpath(@__DIR__, "final_v6.jl"))
isdefined(Main, :_nr_model) || include(joinpath(@__DIR__, "null_rules_real.jl"))      # model nulli (tpa_zap_plots z zkw)
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

# ---- Wykresy porownawcze: strona zbiorcza + strona na pulsar (stos impulsow z zaznaczeniami, z(RFI), profile).
#   include("tools/tpa_zap_compare.jl"); tpa_zap_plots()          # -> tpa_zap/porownanie.pdf
using PyPlot, PyCall
function tpa_zap_plots(; out=TZ_DIR, only=nothing, zkw=nothing, name="porownanie.pdf")
    # zkw = ustawienia automatu (np. (rfi_struct=2.0, rfi_struct_z=6.0)): nasze wyciecia liczone od nowa ta regula
    # (v9 + zkw) zamiast z pazur2 - porownanie proponowanej reguly z TPA przed jej przyjeciem
    _VD = SpaTs.Data
    mkpath(out)
    T = _tz_tpa(); ours = Set(_f6_psrs())
    psrs = sort([q for (q, t) in T if q in ours && t.nz > 0])
    isnothing(only) || (psrs = filter(in(Set(only)), psrs))
    PdfPages = pyimport("matplotlib.backends.backend_pdf").PdfPages
    pf = out * name; pdf = PdfPages(pf)
    S = []                                   # do strony zbiorczej
    pages = []
    for psr in psrs
        try
            t = T[psr]
            p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
            D = Float64.(_f6_load(psr, p)); N, nb = size(D)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            if isnothing(zkw)
                rfi = Set(_tz_rows0(_VD.zaps_rfi_of(p), N)); nul = Set(_tz_rows0(_VD.zaps_null_of(p), N))
                allc0 = Set(_tz_rows0(_VD.zaps_of(p), N))
            else
                m = _nr_model(psr); extra = isnothing(m) ? Int[] : m.rows[m.P .>= 0.5]
                a = _VD.auto_zap(D, bst, bend; null_mode=:model, null_extra=extra, zkw...)
                sp = _VD.split_zaps(a.reason)
                rfi = Set(sp.rfi .- 1); nul = Set(sp.null .- 1); allc0 = Set(a.bad .- 1)
            end
            allc = allc0; tz = Set(i for i in t.zaps if 0 <= i < N)
            z = _VD.rfi_rms_z(D, bst, bend)
            com = intersect(tz, rfi); tonly = setdiff(tz, allc); oonly = setdiff(rfi, tz)
            snr0 = _VD.profile_snr(D, bst, bend)
            snr_o = _VD.profile_snr(D, bst, bend; drop=sort(collect(allc)) .+ 1)
            snr_t = _VD.profile_snr(D, bst, bend; drop=sort(collect(tz)) .+ 1)
            snr_ot = _VD.profile_snr(D, bst, bend; drop=sort(collect(union(allc, tz))) .+ 1)
            push!(S, (psr=psr, com=length(com), tonly=length(tonly), oonly=length(oonly), asnull=length(intersect(setdiff(tz, rfi), nul)),
                      zt=[z[i + 1] for i in tonly], zc=[z[i + 1] for i in com], snr0=snr0, snr_o=snr_o, snr_t=snr_t, snr_ot=snr_ot))
            push!(pages, (psr, D, bst, bend, z, tz, rfi, nul, allc, S[end]))
        catch e
            e isa InterruptException && rethrow()
            println("$psr BŁĄD: ", first(sprint(showerror, e), 150))
        end
    end
    # strona zbiorcza
    fig = figure(figsize=(11, 7.5))
    ax = fig.add_axes([0.08, 0.08, 0.42, 0.85])
    y = collect(1:length(S))
    ax.barh(y, [s.com for s in S], color="purple", label="wycięte przez oba")
    ax.barh(y, [s.tonly for s in S], left=[s.com for s in S], color="red", label="tylko TPA (u nas zostały)")
    ax.barh(y, [s.oonly for s in S], left=[s.com + s.tonly for s in S], color="royalblue", label="tylko nasze RFI")
    ax.set_yticks(y); ax.set_yticklabels([s.psr for s in S], fontsize=6); ax.invert_yaxis()
    ax.set_xscale("symlog", linthresh=10); ax.set_xlabel("liczba impulsów"); ax.legend(fontsize=7, loc="lower right")
    ax.set_title(@sprintf("Wycięte impulsy: TPA %d, wspólne %d, tylko TPA %d, tylko nasze RFI %d", sum(s.com + s.tonly for s in S; init=0) + sum(s.asnull for s in S; init=0),
                          sum(s.com for s in S; init=0), sum(s.tonly for s in S; init=0), sum(s.oonly for s in S; init=0)), fontsize=8)
    ax2 = fig.add_axes([0.58, 0.58, 0.38, 0.35])
    zt = filter(isfinite, vcat([s.zt for s in S]...)); zc = filter(isfinite, vcat([s.zc for s in S]...))
    bins = vcat(-5:1:10, 15, 20, 50, 100, 300)
    ax2.hist(clamp.(zc, -5, 300), bins=bins, histtype="step", color="purple", label="wycięte przez oba")
    ax2.hist(clamp.(zt, -5, 300), bins=bins, histtype="step", color="red", label="tylko TPA")
    ax2.axvline(3, c="k", ls="--", lw=0.8); ax2.set_xscale("symlog", linthresh=10)
    ax2.set_xlabel("z (RFI): szum poza impulsem względem sąsiadów [σ]", fontsize=7); ax2.set_ylabel("impulsy", fontsize=7)
    ax2.legend(fontsize=7); ax2.tick_params(labelsize=6)
    ax2.set_title(@sprintf("impulsy tylko TPA: z>3 w %d z %d", count(>(3), zt), length(zt)), fontsize=8)
    ax3 = fig.add_axes([0.58, 0.08, 0.38, 0.38])
    r(a, b) = 100 * (a / b - 1)
    ax3.scatter([r(s.snr_o, s.snr0) for s in S], [r(s.snr_t, s.snr0) for s in S], s=12, c="k")
    for s in S; ax3.annotate(s.psr, (r(s.snr_o, s.snr0), r(s.snr_t, s.snr0)), fontsize=4); end
    lim = [min(-1, minimum([r(s.snr_t, s.snr0) for s in S]; init=0.0)), max(1, maximum([r(s.snr_o, s.snr0) for s in S]; init=0.0))]
    ax3.plot(lim, lim, c="grey", ls="--", lw=0.6)
    ax3.set_xlabel("zmiana S/N: nasze wycięcia [%]", fontsize=7); ax3.set_ylabel("zmiana S/N: wycięcia TPA [%]", fontsize=7)
    ax3.tick_params(labelsize=6)
    pdf.savefig(fig); close(fig)
    # strony pulsarow
    for (psr, D, bst, bend, z, tz, rfi, nul, allc, s) in pages
        N, nb = size(D)
        fig = figure(figsize=(11, 7.5))
        fig.suptitle(@sprintf("%s: TPA wycięło %d (wspólne z naszym RFI %d, u nas null %d, u nas zostały %d) | nasze RFI %d (spoza TPA %d)",
                              psr, length(tz), s.com, s.asnull, s.tonly, length(rfi), s.oonly), fontsize=8)
        live = [any(!iszero, @view D[i, :]) for i in 1:N]
        v = vec(D[live, :]); lo, hi = quantile(v, 0.01), quantile(v, 0.995)
        # znaczniki: kolumny TPA | nasze RFI | nasze nulle
        axm = fig.add_axes([0.05, 0.30, 0.04, 0.62])
        for (x, set, c) in ((0, tz, "red"), (1, rfi, "royalblue"), (2, nul, "grey"))
            ys = sort(collect(set)); isempty(ys) || axm.scatter(fill(x, length(ys)), ys, marker="_", s=30, c=c, lw=1.2)
        end
        axm.set_xlim(-0.5, 2.5); axm.set_ylim(-0.5, N - 0.5); axm.set_xticks([0, 1, 2])
        axm.set_xticklabels(["TPA", "RFI", "null"], fontsize=6, rotation=90); axm.set_ylabel("impuls (od 0)", fontsize=7); axm.tick_params(labelsize=6)
        axs = fig.add_axes([0.10, 0.30, 0.50, 0.62], sharey=axm)
        axs.imshow(D, origin="lower", aspect="auto", cmap="Greys_r", vmin=lo, vmax=hi, interpolation="nearest", extent=[-0.5, nb - 0.5, -0.5, N - 0.5])
        axs.axvline(bst, c="orange", lw=0.5); axs.axvline(bend, c="orange", lw=0.5)
        axs.tick_params(labelleft=false, labelsize=6); axs.set_xlabel("bin (pomarańczowe: okno impulsu)", fontsize=7)
        axz = fig.add_axes([0.62, 0.30, 0.18, 0.62], sharey=axm)
        axz.plot(clamp.(z, -5, 1000), 0:N-1, c="grey", lw=0.5)
        axz.axvline(3, c="k", ls="--", lw=0.7)
        for (set, c, m) in ((setdiff(tz, allc), "red", "o"), (intersect(tz, rfi), "purple", "o"), (setdiff(rfi, tz), "royalblue", "x"))
            ys = sort(collect(set)); isempty(ys) || axz.scatter(clamp.([z[i + 1] for i in ys], -5, 1000), ys, s=8, c=c, marker=m, lw=0.8)
        end
        axz.set_xscale("symlog", linthresh=10); axz.tick_params(labelleft=false, labelsize=6)
        axz.set_xlabel("z (RFI) [σ], próg 3", fontsize=7)
        axl = fig.add_axes([0.82, 0.30, 0.16, 0.62]); axl.axis("off")
        for (k, (c, lab)) in enumerate((("red", "o tylko TPA (u nas zostały)"), ("purple", "o wycięte przez oba"), ("royalblue", "x tylko nasze RFI"), ("grey", "— nasze nulle")))
            axl.text(0, 1 - 0.07k, lab, color=c, fontsize=7, transform=axl.transAxes)
        end
        axl.text(0, 0.55, @sprintf("S/N profilu:\n bez wycięć %.1f\n nasze %.1f\n TPA %.1f\n nasze + TPA %.1f", s.snr0, s.snr_o, s.snr_t, s.snr_ot),
                 fontsize=7, transform=axl.transAxes, va="top")
        # profile
        axp = fig.add_axes([0.10, 0.06, 0.50, 0.18])
        prof(drop) = vec(mean(D[[i for i in 1:N if live[i] && !(i - 1 in drop)], :], dims=1))
        p0 = prof(Set{Int}()); pm = max(maximum(p0), eps())
        axp.plot(p0 ./ pm, c="grey", lw=0.8, label="bez wycięć")
        axp.plot(prof(allc) ./ pm, c="royalblue", lw=0.8, label="nasze")
        axp.plot(prof(tz) ./ pm, c="red", lw=0.8, ls="--", label="TPA")
        axp.set_xlim(-0.5, nb - 0.5); axp.legend(fontsize=6); axp.tick_params(labelsize=6); axp.set_xlabel("bin", fontsize=7)
        pdf.savefig(fig); close(fig)
    end
    pdf.close()
    println("Wykresy: $pf  ($(length(pages)) pulsarów + strona zbiorcza)")
    return pf
end

"""
Porownanie z TPA (Xiaoxi, stats.csv) na wybranych pulsarach: obecne v9 (pazur2) i proponowana regula RFI
(struktura tla poza oknem > 2x i z > 6). Dwa pliki: tpa_zap/porownanie_v9.pdf i porownanie_regula.pdf.
    include("tools/tpa_zap_compare.jl"); tpa_rule_compare()
"""
function tpa_rule_compare(psrs=["J0932-3217", "J1915+0752", "J1057-7914", "J1926+1648", "J1808-0813", "J1404+1159"];
                          zkw=(rfi_struct=2.0, rfi_struct_z=6.0))
    tpa_zap_plots(; only=psrs, name="porownanie_v9.pdf")
    tpa_zap_plots(; only=psrs, zkw=zkw, name="porownanie_regula.pdf")
end

"""
Ktora zasada najblizej odtwarza wyciecia TPA (Xiaoxi, stats.csv): na wspolnych pulsarach z wycieciami TPA, dla kazdej
zasady - ile impulsow TPA wycinamy (odtworzone), ile naszych wyciec RFI jest poza lista TPA. Uwaga: TPA wycinalo tylko
w 79 pulsarach (tam, gdzie nie bylo dlugiego czystego odcinka), wiec "poza lista TPA" nie musi znaczyc blednie.
    include("tools/tpa_zap_compare.jl"); tpa_rule_scores()
"""
function tpa_rule_scores(; variants=[("v9", NamedTuple()), ("struktura >2x i z>6", (rfi_struct=2.0, rfi_struct_z=6.0)),
                                     ("struktura >2x i z>4", (rfi_struct=2.0, rfi_struct_z=4.0)), ("TPA VIII 4σ", (rfi_tpa_k=4.0,))],
                         out=TZ_DIR)
    _VD = SpaTs.Data
    mkpath(out)
    T = _tz_tpa(); ours = Set(_f6_psrs())
    psrs = sort([q for (q, t) in T if q in ours && t.nz > 0])
    tot = Dict(l => zeros(Int, 4) for (l, _) in variants)      # TPA, odtworzone, nasze RFI, nasze poza TPA
    per = Dict{String,Vector{String}}()
    for psr in psrs
        try
            t = T[psr]
            p = SpaTs.Tools.read_params(joinpath(F6_ROOT, "pazur2", psr, "params.json"))
            D = Float64.(_f6_load(psr, p)); N = size(D, 1)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            m = _nr_model(psr); extra = isnothing(m) ? Int[] : m.rows[m.P .>= 0.5]
            tz = Set(i + 1 for i in t.zaps if 0 <= i < N)               # od 1
            row = String[]
            for (l, kw) in variants
                a = _VD.auto_zap(D, bst, bend; null_mode=:model, null_extra=extra, rfi_struct=Inf, kw...)
                rfi = Set(_VD.split_zaps(a.reason).rfi); cut = Set(a.bad)
                v = [length(tz), length(intersect(tz, cut)), length(rfi), length(setdiff(rfi, tz))]
                tot[l] .+= v
                push!(row, @sprintf("%d/%d", v[2], v[1]))
            end
            per[psr] = row
            println(rpad(psr, 12), " ", join(row, "  "))
        catch e
            e isa InterruptException && rethrow()
            println("$psr BŁĄD: ", first(sprint(showerror, e), 150))
        end
    end
    L = ["# Która zasada najbliżej odtwarza wycięcia TPA (Xiaoxi) - $(length(psrs)) wspólnych pulsarów z wycięciami TPA\n",
         "| zasada | impulsy TPA | wycinamy też my | odtworzone | nasze RFI | nasze RFI poza listą TPA |", "|---|---|---|---|---|---|"]
    for (l, _) in variants
        v = tot[l]
        push!(L, @sprintf("| %s | %d | %d | %.0f%% | %d | %d |", l, v[1], v[2], 100v[2] / max(v[1], 1), v[3], v[4]))
    end
    push!(L, "", "## Pulsary (wycięte z listy TPA / lista TPA): " * join([l for (l, _) in variants], " | "))
    for q in sort(collect(keys(per))); push!(L, "- $q: " * join(per[q], " | ")); end
    write(out * "zasady_vs_TPA.txt", join(L, "\n") * "\n")
    println(join(L[1:3+length(variants)], "\n")); println("\nZapisane: $(out)zasady_vs_TPA.txt")
end
