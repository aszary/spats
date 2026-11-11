# Weryfikacja automatu wycinania w JEDNYM przebiegu (tylko odczyt surowych danych;
# w pazur2 dopisywany jest jedynie pulsar_rfi.spCF = archiwum z wycietym tylko RFI,
# nulle zostawione - do badan nullingu). Dla kazdego pulsara z pazur2:
#   A. spojnosc: S/N bez wyciec vs zap_cmp, czy automat odtwarza zapisane zaps2;
#   B. nulle: profil sredni z impulsow zaps2_null -> czy jest w nim sygnal
#      (S/N stosu, ulamek energii dobrego impulsu); frakcja nulli: nasza vs
#      estymator z rozkladu energii (N(E_on<0)/N(E_off<0)) vs literatura;
#   C. nowy bezpiecznik RFI (rfi_ref=:both) vs stary (:cut) na tych samych danych;
#   D. kanaly (.spCf16): chan_rfi_core -> kanaly do wyciecia; automat na danych
#      zsumowanych ze wszystkich kanalow vs bez zlych kanalow (czy najpierw kanaly).
# Wyniki: /home/psr/data/OUTPUT/zap_verify/verify.csv (dopisywane, wznawialne),
# nulls.pdf (profile nulli), summary.txt (do wklejenia). Uruchamiane z main().
using Statistics, Printf
isdefined(Main, :SpaTs) || error("Najpierw include(\"spats.jl\").")
include(joinpath(@__DIR__, "chan_rfi.jl"))   # zawsze: po ponownym include("spats.jl") SpaTs jest nowym modulem
_VD, _VT, _VP = SpaTs.Data, SpaTs.Tools, SpaTs.Plot

# frakcje nullingu z literatury (przeglad Claude online, 2026-11-11); rozne
# teleskopy/metody - porownanie orientacyjne
LIT_NF = Dict(
    "J0034-0721" => (44.0, "44±1% Gajjar 2017"), "J1049-5833" => (47.0, "47±3% Wang+2007"),
    "J1255-6131" => (75.0, "~70-80% (10-100% miedzy obs.) Burke-Spolaor+2012"),
    "J1537-4912" => (20.0, "dlugie nulle do ~20% TPA XVIII"), "J1727-2739" => (52.0, "52±3% Wang+2007"),
    "J1850-0006" => (48.0, "48±3% Yan+2026 FAST"), "J1851-0053" => (45.0, "47±4% / 42±4% Yan+2026"),
    "J1926-0652" => (75.0, "~75% (b. niepewne) Zhang+2019"), "J1931+1439" => (42.0, "~42% Yan+2026"),
    "J1944+1755" => (60.0, ">=60% Lorimer+2002; FAST ~48%"), "J1946+1805" => (55.0, "55±5% Yang+2014"),
    "J2037+1942" => (44.0, "44±4% Herfindal&Rankin 2009 (albo 24%)"),
    "J2253+1516" => (49.0, "<=49% Redman&Rankin 2009 (gorna granica)"))

VERIFY_HEADER = "pulsar,status,N,n_cut,n_rfi,n_null,snr_none,snr_none_cmp,repro," *
    "null_snr,null_frac,nf_ours,nf_energy,nf_lit,n_null_new,null_snr_new,null_frac_new,nf_new," *
    "old_ncut,old_snr,old_snr_rfi,new_ncut,new_snr,new_snr_rfi," *
    "C,chans_bad,chans_loo,chans_cut,n_sub_bad,pair_pct," *
    "ch_all_none,ch_all_ncut,ch_all_snr,ch_good_none,ch_good_ncut,ch_good_snr,note"

# --- szybki odczyt `pdv -t` (16 kanalow = ~16 mln linii; split() byl za wolny) ---
@inline function _pint(buf, j)
    while buf[j] == UInt8(' ') || buf[j] == UInt8('\t'); j += 1; end
    v = 0
    while UInt8('0') <= buf[j] <= UInt8('9'); v = 10v + (buf[j] - UInt8('0')); j += 1; end
    return v, j
end
_load_cube_fast(txt::AbstractString) = _load_cube_fast(read(txt))
function _load_cube_fast(buf::Vector{UInt8})
    n = length(buf); pos = 1
    data = nothing
    GC.@preserve buf begin
        p0 = pointer(buf)
        while pos <= n
            e = something(findnext(==(UInt8('\n')), buf, pos), n + 1)
            c = buf[pos]
            if !isnothing(data) && UInt8('0') <= c <= UInt8('9')
                i, j = _pint(buf, pos); ch, j = _pint(buf, j); b, j = _pint(buf, j)
                v = ccall(:strtod, Float64, (Ptr{UInt8}, Ptr{Ptr{UInt8}}), p0 + j - 1, C_NULL)
                data[i+1, ch+1, b+1] = v
            elseif isnothing(data) && c == UInt8('F')
                r = split(String(buf[pos:e-1]))
                if r[1] == "File:"
                    data = zeros(Float32, parse(Int, r[6]), parse(Int, r[8]), parse(Int, r[12]))
                end
            end
            pos = e + 1
        end
    end
    return data
end

# S/N profilu zlozonego z `rows` (off-pulse z calych danych)
function _stack_snr(D, rows, on, off)
    isempty(rows) && return NaN
    P = vec(sum(@view(D[rows, :]), dims=1))
    μ, σ = mean(P[off]), std(P[off])
    return σ > 0 ? sum(P[on] .- μ) / (σ * sqrt(length(on))) : NaN
end

"Frakcja nulli z rozkładu energii: N(E_on<0)/N(E_off<0), okno off tej samej szerokości."
function _nf_energy(D, rows, on, off)
    n_on = length(on)
    segs = UnitRange{Int}[]
    for j in off
        (!isempty(segs) && last(segs[end]) == j - 1) ? (segs[end] = first(segs[end]):j) : push!(segs, j:j)
    end
    sg = segs[argmax(length.(segs))]
    length(sg) >= n_on || return NaN
    w = first(sg):first(sg)+n_on-1
    Eon = [sum(@view D[i, on]) for i in rows]
    Eoff = [sum(@view D[i, w]) for i in rows]
    k = count(<(0), Eoff)
    return k >= 10 ? clamp(count(<(0), Eon) / k, 0.0, 1.0) : NaN
end

_jacc(a, b) = (u = length(union(a, b)); u == 0 ? 1.0 : length(intersect(a, b)) / u)

function _debased_spcf(indir, tmp, bst, bend; chans16="")
    ar = joinpath(tmp, "pulsar.spCF")
    redirect_stdout(devnull) do
        if isempty(chans16)
            _VD.add_psrfiles(indir, ar)
        else                                  # jak autozap_all z kanalami wycietymi przed sumowaniem
            _VD.add_psrfiles(indir, ar; sixteen=true)
            try run(pipeline(`pam -p -m $ar`, stdout=devnull, stderr=devnull)) catch; end
            run(pipeline(`paz -m -z $chans16 $ar`, stdout=devnull, stderr=devnull))
            run(pipeline(`pam -F -m $ar`, stdout=devnull, stderr=devnull))
        end
    end
    cp(ar, joinpath(tmp, "raw.spCF"))
    run(pipeline(`pmod -onpulse "$bst $bend" -device "/NULL" -debase $ar`; stdout=devnull, stderr=devnull))
    txt = joinpath(tmp, "pulsar.debase.txt")
    _VD.convert_psrfit_ascii(joinpath(tmp, "pulsar.debase.gg"), txt)
    return _VD.load_ascii(txt)
end

function _cube16(indir, tmp)
    any(f -> endswith(f, ".spCf16"), readdir(indir)) || return nothing
    ar = joinpath(tmp, "p16.ar")
    redirect_stdout(devnull) do
        _VD.add_psrfiles(indir, ar; sixteen=true)
    end
    try run(pipeline(`pam -p -m $ar`, stdout=devnull, stderr=devnull)) catch; end
    d = _load_cube_fast(read(pipeline(`pdv -t $ar`, stderr=devnull)))   # prosto do pamieci, bez pliku
    return d, ar
end

# suma czestotliwosci PRZEZ PSRCHIVE (z wagami kanalow, jak .spCF): opcjonalnie
# `paz -z` zlych kanalow, `pam -F`, debase z oknem -> macierz (impuls, bin)
function _fsum_debased(ar16, tmp, chans, bst, bend, tag)
    a = joinpath(tmp, "f_$tag.spCF")
    cp(ar16, a; force=true)
    isempty(chans) || run(pipeline(`paz -m -z $(join(chans, " ")) $a`, stdout=devnull, stderr=devnull))
    run(pipeline(`pam -F -m $a`, stdout=devnull, stderr=devnull))
    run(pipeline(`pmod -onpulse "$bst $bend" -device "/NULL" -debase $a`; stdout=devnull, stderr=devnull))
    txt = joinpath(tmp, "f_$tag.debase.txt")
    _VD.convert_psrfit_ascii(joinpath(tmp, "f_$tag.debase.gg"), txt)
    return _VD.load_ascii(txt)
end

"Środkowe częstotliwości kanałów (MHz) z `psredit -c freq,bw`; `nothing`, gdy się nie uda."
function _chan_freqs(ar, C)
    try
        o = read(pipeline(`psredit -c freq,bw $ar`, stderr=devnull), String)
        f = parse(Float64, match(r"freq=\s*([-+0-9.eE]+)", o).captures[1])
        bw = parse(Float64, match(r"bw=\s*([-+0-9.eE]+)", o).captures[1])
        return [f + (c - (C + 1) / 2) * bw / C for c in 1:C]
    catch
        return nothing
    end
end

"""
PDF z uzasadnieniem decyzji o kanałach dla jednego pulsara (strona 1: cechy
każdego kanału z progami reguł + S/N kanału + zmiana S/N sumy bez kanału;
strona 2: stosy impulsów oznaczonych kanałów obok najlepszego kanału).
"""
function _chan_pdf(path, psr, d, a, bst, bend, freqs, snr_c, loo, snr_all, cut)
    plt = _VP.PyPlot
    PdfPages = _VP.PyCall.pyimport("matplotlib.backends.backend_pdf").PdfPages
    N, C, nb = size(d)
    g, off = a.g, a.off
    chl = [isnothing(freqs) ? "$(c-1)" : @sprintf("%d\n%.0f", c - 1, freqs[c]) for c in 1:C]
    col = [c in cut ? "C3" : c in a.chan_all ? "C1" : a.dead[c] ? "0.8" : "C0" for c in 1:C]
    x = 0:C-1
    pdf = PdfPages(path)

    fig, ax = plt.subplots(2, 3, figsize=(17, 10.5))
    # 1. czestotliwosc x faza: profil sredni kazdego kanalu w jednostkach szumu
    w = bend - bst + 1
    lo, hi = max(1, bst - 2w), min(nb, bend + 2w)
    M = zeros(C, hi - lo + 1)
    for c in 1:C
        a.dead[c] && continue
        P = vec(sum(Float64.(@view d[g, c, :]), dims=1))
        σ = std(P[off]); σ > 0 && (M[c, :] = (P[lo:hi] .- median(P[off])) ./ σ)
    end
    ax[1, 1].imshow(M, aspect="auto", origin="lower", cmap="viridis", extent=(lo, hi, -0.5, C - 0.5),
                    vmin=-3, vmax=max(5, quantile(vec(M), 0.995)), interpolation="nearest")
    ax[1, 1].axvline(bst, color="w", lw=0.6, ls="--"); ax[1, 1].axvline(bend, color="w", lw=0.6, ls="--")
    ax[1, 1].set_yticks(x); ax[1, 1].set_yticklabels(chl, fontsize=6)
    ax[1, 1].set_title("profil każdego kanału [σ szumu]; kanał / MHz", fontsize=9); ax[1, 1].set_xlabel("bin")
    # 2. S/N kanalu
    ax[1, 2].bar(x, snr_c, color=col)
    ax[1, 2].set_title("S/N profilu kanału (sygnał pulsara w kanale)", fontsize=9)
    # 3. zmiana S/N sumy bez kanalu (PSRCHIVE, z wagami) - KRYTERIUM wyciecia
    lv = [isnan(v) ? 0.0 : v for v in loo]
    ax[1, 3].bar(x, lv, color=col); ax[1, 3].axhline(1.0, color="C3", ls="--", lw=1)
    ax[1, 3].axhline(0, color="k", lw=0.5)
    ax[1, 3].set_title("zmiana S/N sumy BEZ kanału [%] (tylko oznaczone)\n" *
                       "wycinamy, gdy > +1% (czerwona linia)", fontsize=9)
    # 4-6. reguly detektora z progami
    ax[2, 1].bar(x, 100 .* a.fch, color=col); ax[2, 1].axhline(100 * a.fchan, color="C3", ls="--", lw=1)
    ax[2, 1].set_title("reguła 'czas': % impulsów z zakłóceniem w kanale (próg $(round(Int, 100a.fchan))%)", fontsize=9)
    ax[2, 2].bar(x, [isnan(v) ? 0.0 : v for v in a.rmsc], color=col)
    isfinite(a.rms_thr) && ax[2, 2].axhline(a.rms_thr, color="C3", ls="--", lw=1)
    ax[2, 2].set_title("reguła 'rms': szum off-pulse kanału (próg Q3+1.5 IQR)\n" *
                       "surowe dane: zależy też od wzmocnienia kanału", fontsize=9)
    ax[2, 3].bar(x, [isnan(v) ? 0.0 : v for v in a.acfc], color=col)
    ax[2, 3].axhline(a.acf_thr, color="C3", ls="--", lw=1)
    ax[2, 3].set_title("reguła 'acf': |ACF1| off-pulse (paski, próg mediana+4σ)", fontsize=9)
    for k in (ax[1, 2], ax[1, 3], ax[2, 1], ax[2, 2], ax[2, 3])
        k.set_xticks(x); k.set_xticklabels(chl, fontsize=6); k.set_xlabel("kanał / MHz", fontsize=8)
    end
    oz = join(["$(c-1)($(a.why[c]))" for c in a.chan_all], " ")
    fig.suptitle("$psr — kanały.  niebieski = OK, pomarańczowy = oznaczony przez detektor, czerwony = wycięcie " *
                 "uzasadnione (S/N sumy rośnie > 1%)\noznaczone: $(isempty(oz) ? "brak" : oz)   →   " *
                 "do wycięcia: $(isempty(cut) ? "BRAK (wycięcie obniżyłoby S/N)" : join(cut .- 1, " "))   " *
                 @sprintf("(S/N sumy wszystkich kanałów: %.1f)", snr_all), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93)); pdf.savefig(fig); plt.close(fig)

    # strona 2: stosy impulsow
    ref = argmax([a.dead[c] || c in a.chan_all ? -Inf : snr_c[c] for c in 1:C])
    show = vcat(a.chan_all[1:min(5, end)], ref)
    fig, ax = plt.subplots(1, length(show), figsize=(3.4 * length(show), 9), squeeze=false)
    for (k, c) in enumerate(show)
        A = Float64.(d[:, c, lo:hi]) .- a.F[:, c, 4]
        A[.!g, :] .= NaN
        v = filter(isfinite, vec(A))
        vmin, vmax = isempty(v) ? (0.0, 1.0) : (quantile(v, 0.01), quantile(v, 0.995))
        ax[1, k].imshow(A, aspect="auto", origin="lower", cmap="viridis", vmin=vmin, vmax=vmax,
                        extent=(lo, hi, 0, N), interpolation="nearest")
        lab = c == ref ? "najlepszy kanał (wzór)" : "oznaczony: $(a.why[c])"
        fr = isnothing(freqs) ? "" : @sprintf(" (%.0f MHz)", freqs[c])
        ax[1, k].set_title(@sprintf("kanał %d%s\n%s\nS/N kanału %.1f%s", c - 1, fr, lab, snr_c[c],
                           isnan(loo[c]) ? "" : @sprintf("\nS/N sumy bez niego %+.1f%%", loo[c])), fontsize=8)
        ax[1, k].set_xlabel("bin")
    end
    ax[1, 1].set_ylabel("impuls")
    fig.suptitle("$psr — impulsy w oznaczonych kanałach obok najlepszego kanału (skala 1-99.5% każdego panelu)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); pdf.savefig(fig); plt.close(fig)
    pdf.close()
end

function _verify_one(psr, p, st, indir, outpsr; chans=true, max_chan_pulses=4000, save_rfi=true,
                     chan_pdf_dir=nothing)
    bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
    on = bst:bend
    tmp = mktempdir()
    row = Dict{String,Any}("pulsar" => psr, "status" => st)
    note = String[]
    fig = nothing
    try
        D = _debased_spcf(indir, tmp, bst, bend; chans16=something(get(p, "zap_chans16", nothing), ""))
        N = size(D, 1)
        row["N"] = N
        rows_all = _VD.zaps_to_rows(_VD.zaps_of(p), N)
        rows_rfi = _VD.zaps_to_rows(_VD.zaps_rfi_of(p), N)
        has_split = haskey(p, "zaps2_rfi") || haskey(p, "zaps2_null")
        rows_null = has_split ? _VD.zaps_to_rows(_VD.zaps_null_of(p), N) : Int[]
        has_split || push!(note, "brak podzialu rfi/null")
        row["n_cut"], row["n_rfi"], row["n_null"] = length(rows_all), length(rows_rfi), length(rows_null)
        c = get(p, "zap_cmp", Dict())
        get(c, "n_pulses", N) == N || push!(note, "N inne niz w zap_cmp ($(c["n_pulses"]))")
        row["snr_none"] = _VD.profile_snr(D, bst, bend)
        row["snr_none_cmp"] = something(get(c, "snr_none", nothing), NaN)

        # B. nulle
        off = _VD.offpulse_bins(D, bst, bend)
        good = setdiff(1:N, rows_all)
        base = [mean(@view D[i, off]) for i in 1:N]
        E = [sum(@view D[i, on]) - length(on) * base[i] for i in 1:N]
        row["null_snr"] = _stack_snr(D, rows_null, on, off)
        row["null_frac"] = (isempty(rows_null) || isempty(good)) ? NaN : mean(E[rows_null]) / mean(E[good])
        nonrfi = setdiff(1:N, rows_rfi)
        row["nf_ours"] = has_split && !isempty(nonrfi) ? 100 * length(rows_null) / length(nonrfi) : NaN
        nfe = _nf_energy(D, nonrfi, on, off)
        row["nf_energy"] = 100 * nfe
        row["nf_lit"] = haskey(LIT_NF, psr) ? LIT_NF[psr][1] : NaN

        # C. bezpiecznik RFI: stary vs nowy na tych samych danych
        for (tag, ref) in (("old", :cut), ("new", :both))
            a = _VD.auto_zap(D, bst, bend; rfi_ref=ref)
            row["$(tag)_ncut"] = length(a.bad)
            row["$(tag)_snr"] = _VD.profile_snr(D, bst, bend; drop=a.bad)
            row["$(tag)_snr_rfi"] = _VD.profile_snr(D, bst, bend; drop=_VD.split_zaps(a.reason).rfi)
            tag == "old" && (row["repro"] = _jacc(Set(a.bad), Set(rows_all)))
            if tag == "new"                   # nulle wg nowego automatu (z poprawka 'nulli' z sygnalem)
                pn = _VD.split_zaps(a.reason)
                gn = setdiff(1:N, a.bad)
                row["n_null_new"] = length(pn.null)
                row["null_snr_new"] = _stack_snr(D, pn.null, on, off)
                row["null_frac_new"] = (isempty(pn.null) || isempty(gn)) ? NaN : mean(E[pn.null]) / mean(E[gn])
                nr = setdiff(1:N, pn.rfi)
                row["nf_new"] = isempty(nr) ? NaN : 100 * length(pn.null) / length(nr)
                if !isempty(pn.null) && !isempty(gn)
                    fig = (psr, vec(mean(D[gn, :], dims=1)), vec(mean(D[pn.null, :], dims=1)), bst, bend,
                           row["null_snr_new"], row["null_frac_new"], length(pn.null))
                end
            end
        end

        # archiwum tylko z RFI (nulle zostawione) dla badan nullingu
        if save_rfi && has_split && !isempty(rows_null)
            out = joinpath(outpsr, "pulsar_rfi.spCF")
            cp(joinpath(tmp, "raw.spCF"), out; force=true)
            redirect_stdout(devnull) do
                _VD.apply_zaps(out, _VD.zaps_rfi_of(p); chans=_VD.zap_chans_of(p))
            end
        end

        # D. kanaly (blad tutaj nie kasuje A-C)
        if chans && N <= max_chan_pulses
          try
            dc = _cube16(indir, tmp)
            if isnothing(dc)
                push!(note, "brak .spCf16")
            else
                d, ar16 = dc
                Nc, C, nb = size(d)
                row["C"] = C
                if Nc != N || nb != size(D, 2)
                    push!(note, "16ch: $Nc imp./$nb binow vs $N/$(size(D, 2))")
                else
                    a = chan_rfi_core(d, bst, bend)
                    row["chans_bad"] = join(["$(c-1):$(a.why[c])" for c in a.chan_all], " ")
                    row["n_sub_bad"] = length(a.sub_all)
                    live = count(.!a.dead)
                    row["pair_pct"] = 100 * length(a.pairs) / max(count(a.g) * live, 1)
                    # S/N kazdego kanalu (niezalezne od wzmocnienia)
                    snr_c = fill(NaN, C)
                    for c in findall(.!a.dead)
                        P = vec(sum(Float64.(@view d[a.g, c, :]), dims=1))
                        σ = std(P[a.off])
                        σ > 0 && (snr_c[c] = sum(P[on] .- mean(P[a.off])) / (σ * sqrt(length(on))))
                    end
                    # KRYTERIUM: S/N sumy (PSRCHIVE, z wagami) bez kanalu; wycinamy, gdy > +1%
                    S_all = _fsum_debased(ar16, tmp, Int[], bst, bend, "all")
                    s_all = _VD.profile_snr(S_all, bst, bend)
                    row["ch_all_none"] = s_all
                    loo = fill(NaN, C)
                    for c in a.chan_all
                        loo[c] = 100 * (_VD.profile_snr(_fsum_debased(ar16, tmp, [c - 1], bst, bend, "loo"),
                                                        bst, bend) / s_all - 1)
                    end
                    cut = [c for c in a.chan_all if loo[c] > 1.0]
                    row["chans_loo"] = join([@sprintf("%d:%+.1f", c - 1, loo[c]) for c in a.chan_all], " ")
                    row["chans_cut"] = join(cut .- 1, " ")
                    if !isnothing(chan_pdf_dir)
                        try
                            _chan_pdf(joinpath(chan_pdf_dir, "$psr.pdf"), psr, d, a, bst, bend,
                                      _chan_freqs(ar16, C), snr_c, loo, s_all, cut)
                        catch e
                            push!(note, "pdf kanalow: " * first(replace(sprint(showerror, e), ',' => ';', '\n' => ' '), 80))
                        end
                    end
                    d = nothing                       # kostka niepotrzebna (pamiec)
                    if !isempty(cut)                  # automat: wszystkie kanaly vs bez uzasadnionych
                        z = _VD.auto_zap(S_all, bst, bend)
                        row["ch_all_ncut"], row["ch_all_snr"] = length(z.bad), _VD.profile_snr(S_all, bst, bend; drop=z.bad)
                        S = _fsum_debased(ar16, tmp, cut .- 1, bst, bend, "good")
                        z = _VD.auto_zap(S, bst, bend)
                        row["ch_good_none"] = _VD.profile_snr(S, bst, bend)
                        row["ch_good_ncut"], row["ch_good_snr"] = length(z.bad), _VD.profile_snr(S, bst, bend; drop=z.bad)
                    end
                end
            end
          catch e
            e isa InterruptException && rethrow()
            push!(note, "kanaly: " * first(replace(sprint(showerror, e), '"' => '\'', '\n' => ' ', ',' => ';'), 120))
          end
        elseif chans
            push!(note, "kanaly pominiete (N>$max_chan_pulses)")
        end
    catch e
        e isa InterruptException && rethrow()
        row["status"] = "blad"
        push!(note, first(replace(sprint(showerror, e), '"' => '\'', '\n' => ' ', ',' => ';'), 150))
    finally
        rm(tmp; recursive=true, force=true)
    end
    row["note"] = join(note, "; ")
    return row, fig
end

function _fmtv(x)
    (isnothing(x) || (x isa AbstractFloat && !isfinite(x))) && return ""
    x isa AbstractFloat && return string(round(x, digits=3))
    return string(x)
end

function _save_nullprof(dir, fig)
    psr, Pg, Pn, bst, bend, s, f, nn = fig
    open(joinpath(dir, "$psr.txt"), "w") do io
        println(io, join([psr, bst, bend, s, f, nn], " "))
        foreach(j -> println(io, Pg[j], " ", Pn[j]), eachindex(Pg))
    end
end

function _null_pages(dir, path)
    files = isdir(dir) ? sort(filter(f -> endswith(f, ".txt"), readdir(dir))) : String[]
    isempty(files) && return
    figs = map(files) do f
        L = readlines(joinpath(dir, f))
        h = split(L[1])
        v = [parse.(Float64, split(l)) for l in L[2:end]]
        (h[1], first.(v), last.(v), parse(Int, h[2]), parse(Int, h[3]), parse(Float64, h[4]),
         parse(Float64, h[5]), parse(Int, h[6]))
    end
    sort!(figs, by=x -> -x[6])                       # najsilniejszy sygnal w nullach najpierw
    plt = _VP.PyPlot
    PdfPages = _VP.PyCall.pyimport("matplotlib.backends.backend_pdf").PdfPages
    pdf = PdfPages(path)
    for k in 1:12:length(figs)
        fig, axs = plt.subplots(4, 3, figsize=(15, 16))
        for (t, ax) in enumerate(vec(permutedims(axs)))
            idx = k + t - 1
            if idx > length(figs); ax.axis("off"); continue; end
            psr, Pg, Pn, bst, bend, s, f, nn = figs[idx]
            w = bend - bst + 1
            lo, hi = max(1, bst - w), min(length(Pg), bend + w)
            x = lo:hi
            ax.plot(x, Pg[x], color="0.6", lw=1, label="dobre (średnia)")
            ax.plot(x, Pn[x], color="C0", lw=1, label="nulle (średnia, $nn imp.)")
            ax.axvspan(bst, bend, color="C1", alpha=0.08)
            ax.set_title(@sprintf("%s  S/N nulli=%.1f  energia=%.1f%%", psr, s, 100f), fontsize=9)
            ax.legend(fontsize=7)
        end
        fig.suptitle("Profile średnie (posortowane: najwięcej sygnału w nullach najpierw): impulsy dobre vs wycięte jako nulle (okno on-pulse na pomarańczowo)")
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)
    end
    pdf.close()
end

"""
    zap_verify(; pazur2, out, only, chans=true)

Pełna weryfikacja (patrz nagłówek pliku). Wznawialna: pulsary już w
`out/verify.csv` są pomijane. Kolejność: najpierw pulsary z nullami i
podejrzane, potem reszta. Na końcu `zap_verify_summary`.
"""
function zap_verify(; pazur2="/home/psr/data/OUTPUT/pazur2/", dataroot="/home/psr/data/new/",
                    out="/home/psr/data/OUTPUT/zap_verify/", only=nothing, chans=true, save_rfi=true,
                    part=nothing)
    mkpath(out)
    kp, np_ = isnothing(part) ? (1, 1) : parse.(Int, split(part, "/"))
    csv = joinpath(out, np_ == 1 ? "verify.csv" : "verify_$kp.csv")
    done = Set{String}()
    for f in _verify_csvs(out), l in Iterators.drop(eachline(f), 1)
        push!(done, split(l, ",")[1])
    end
    isfile(csv) || open(f -> println(f, VERIFY_HEADER), csv, "w")
    summ = joinpath(pazur2, "zap_summary.csv")
    st = Dict{String,String}(); prio = Dict{String,Float64}()
    for l in Iterators.drop(eachline(summ), 1)
        r = split(l, ",")
        length(r) < 13 && continue
        st[r[1]] = r[2]
        v(i) = (x = tryparse(Float64, r[i]); isnothing(x) ? 0.0 : x)
        prio[r[1]] = v(13) + v(12) + 10 * v(5)          # nulle, RFI, % wycietych
    end
    psrs = sort([k for (k, s) in st if s in ("ok", "reczny")], by=k -> -prio[k])
    isnothing(only) || (psrs = filter(in(Set(only)), psrs))
    psrs = psrs[kp:np_:end]                              # czesc kp z np_ (rownolegle procesy)
    todo = filter(!in(done), psrs)
    println("WERYFIKACJA (część $kp/$np_): $(length(todo)) pulsarów do zrobienia ($(length(done)) już zrobionych w $out)")
    flush(stdout)                                      # log w pliku: bez flush widac go dopiero po ~kilkudziesieciu KB
    profdir = joinpath(out, "nullprof"); mkpath(profdir)
    chandir = joinpath(out, "kanaly"); mkpath(chandir)
    keys_ = split(VERIFY_HEADER, ",")
    t0 = time()
    for (k, psr) in enumerate(todo)
        pf = joinpath(pazur2, psr, "params.json")
        isfile(pf) || continue
        p = _VT.read_params(pf)
        psr_path = joinpath(dataroot, psr)
        obs = sort(filter(d -> isdir(joinpath(psr_path, d)), readdir(psr_path)))
        isempty(obs) && continue
        row, fig = _verify_one(psr, p, st[psr], joinpath(psr_path, obs[1]) * "/", joinpath(pazur2, psr);
                               chans=chans, save_rfi=save_rfi, chan_pdf_dir=chandir)
        isnothing(fig) || _save_nullprof(profdir, fig)
        open(f -> println(f, join([_fmtv(get(row, kk, nothing)) for kk in keys_], ",")), csv, "a")
        el = (time() - t0) / k
        g(x) = (v = get(row, x, nothing); v isa Number ? v : NaN)
        @printf("[%d/%d] %-12s %s  nulle: energia %s%%  kanały oznaczone: %s  zmiana S/N bez: %s  wyciąć: %s  %s  (%.0f s/pulsar, zostało ~%.1f h)\n",
                k, length(todo), psr, row["status"], _fmtv(100 * g("null_frac")), _fmtv(get(row, "chans_bad", "-")),
                _fmtv(get(row, "chans_loo", "-")), _fmtv(get(row, "chans_cut", "-")), row["note"], el,
                el * (length(todo) - k) / 3600)
        flush(stdout)
        k % 50 == 0 && (try _null_pages(profdir, joinpath(out, "nulls.pdf")) catch e; println("PDF: $e"); end)
    end
    try _null_pages(profdir, joinpath(out, "nulls.pdf")) catch e; println("PDF nulli nie powstał: $e"); end
    return zap_verify_summary(out=out)
end

"""
Podsumowanie z `verify.csv` (działa w każdej chwili, także w trakcie
przebiegu). Wypisuje i zapisuje `summary.txt` — to wklej.
"""
_verify_csvs(out) = sort(filter(f -> occursin(r"^verify(_\d+)?\.csv$", basename(f)),
                                readdir(out; join=true)))

function zap_verify_summary(; out="/home/psr/data/OUTPUT/zap_verify/")
    H = split(VERIFY_HEADER, ",")
    R = Dict{SubString{String},SubString{String}}[]
    for f in _verify_csvs(out), l in Iterators.drop(eachline(f), 1)
        push!(R, Dict(zip(H, split(l, ","; limit=length(H)))))
    end
    num(r, k) = (x = tryparse(Float64, get(r, k, "")); isnothing(x) ? NaN : x)
    ok = filter(r -> r["status"] != "blad", R)
    io = IOBuffer()
    pr(s="") = (println(s); println(io, s))
    fin(v) = filter(isfinite, v)
    pr("===== WERYFIKACJA AUTOMATU: $(length(R)) pulsarów (błędy: $(length(R) - length(ok))) =====")
    for r in R
        r["status"] == "blad" && pr("  BŁĄD $(r["pulsar"]): $(r["note"])")
    end

    pr("\n--- A. Spójność ---")
    rat = fin([num(r, "snr_none") / num(r, "snr_none_cmp") for r in ok])
    rep = fin([num(r, "repro") for r in ok if r["status"] == "ok"])
    pr(@sprintf("S/N bez wycięć vs zapisane: mediana %.3f, |różnica|>2%%: %d z %d", median(rat),
                count(x -> abs(x - 1) > 0.02, rat), length(rat)))
    pr(@sprintf("automat odtwarza zapisane zaps2 (Jaccard): mediana %.3f, <0.9: %d z %d", median(rep),
                count(<(0.9), rep), length(rep)))

    pr("\n--- B. Nulle: czy w profilu średnim wyciętych nulli jest sygnał ---")
    NR = filter(r -> isfinite(num(r, "null_snr")), ok)
    s = [num(r, "null_snr") for r in NR]
    en = [num(r, "null_frac") for r in NR]
    pr(@sprintf("pulsary z nullami: %d;  S/N profilu nulli <3 (nic): %d,  >=3 (coś jest): %d", length(NR),
                count(<(3), s), count(>=(3), s)))
    pr(@sprintf("energia nulli vs dobry impuls: <2%% (czyste nulle): %d,  2-10%% (domieszka emisji): %d,  >10%% (to nie nulle): %d",
                count(<(0.02), en), count(x -> 0.02 <= x < 0.1, en), count(>=(0.1), en)))
    pr("(S/N rośnie z jasnością pulsara - o czystości nulli mówi energia)  Najwięcej energii w nullach:")
    pr(@sprintf("  %-12s %7s %9s %8s %8s   | po poprawce: %6s %8s %8s", "pulsar", "n_null", "S/N nulli",
                "energia", "NF nasz", "n_null", "energia", "NF"))
    for r in sort(NR, by=r -> -num(r, "null_frac"))[1:min(25, end)]
        pr(@sprintf("  %-12s %7s %9.1f %7.1f%% %7.1f%%   |              %6s %7.1f%% %7.1f%%", r["pulsar"],
                    r["n_null"], num(r, "null_snr"), 100num(r, "null_frac"), num(r, "nf_ours"),
                    get(r, "n_null_new", ""), 100num(r, "null_frac_new"), num(r, "nf_new")))
    end
    en2 = fin([num(r, "null_frac_new") for r in ok])
    pr(@sprintf("PO POPRAWCE ('nulle' z sygnałem -> zaszumione): energia <2%%: %d,  2-10%%: %d,  >10%%: %d",
                count(<(0.02), en2), count(x -> 0.02 <= x < 0.1, en2), count(>=(0.1), en2)))
    fr = fin([num(r, "null_frac") for r in NR])
    isempty(fr) || pr(@sprintf("energia nulli jako %% dobrego impulsu: mediana %.1f%%, >10%%: %d", 100median(fr),
                               count(>(0.1), fr)))

    pr("\n--- B2. Frakcja nulli: nasza vs z rozkładu energii vs literatura ---")
    pr(@sprintf("  %-12s %8s %9s %8s  %s", "pulsar", "NF nasz", "NF energ.", "NF lit.", "źródło"))
    for r in sort(filter(r -> haskey(LIT_NF, r["pulsar"]), ok), by=r -> r["pulsar"])
        pr(@sprintf("  %-12s %7.1f%% %8.1f%% %7.1f%%  %s  (po poprawce: %.1f%%)", r["pulsar"], num(r, "nf_ours"),
                    num(r, "nf_energy"), num(r, "nf_lit"), LIT_NF[r["pulsar"]][2], num(r, "nf_new")))
    end
    d1 = fin([num(r, "nf_ours") - num(r, "nf_lit") for r in ok if haskey(LIT_NF, r["pulsar"])])
    d2 = fin([num(r, "nf_ours") - num(r, "nf_energy") for r in NR if num(r, "nf_ours") > 5])
    isempty(d1) || pr(@sprintf("nasz - literatura: mediana %+.1f pp (z %d)", median(d1), length(d1)))
    isempty(d2) || pr(@sprintf("nasz - energia (pulsary z NF>5%%): mediana %+.1f pp (z %d)", median(d2), length(d2)))

    pr("\n--- C. Bezpiecznik RFI: stary (:cut) vs nowy (:both), te same dane ---")
    drop(r, t) = 100 * (num(r, "$(t)_snr_rfi") / num(r, "snr_none") - 1)
    for t in ("old", "new")
        v = fin([drop(r, t) for r in ok])
        pr(@sprintf("%s: dane tylko po RFI tracą >1%% S/N: %d pulsarów; >2%%: %d", t == "old" ? "stary" : "nowy ",
                    count(<(-1), v), count(<(-2), v)))
    end
    ch = filter(r -> num(r, "old_ncut") != num(r, "new_ncut"), ok)
    pr("zmienione przez nowy bezpiecznik: $(length(ch))")
    for r in ch[1:min(25, end)]
        pr(@sprintf("  %-12s cięć %4d -> %4d   S/N %.1f -> %.1f   tylko-RFI vs bez: %+.1f%% -> %+.1f%%", r["pulsar"],
                    num(r, "old_ncut"), num(r, "new_ncut"), num(r, "old_snr"), num(r, "new_snr"),
                    drop(r, "old"), drop(r, "new")))
    end

    pr("\n--- D. Kanały (.spCf16): czy najpierw wycinać kanały ---")
    pr("Kryterium: kanał wycinamy tylko, gdy S/N sumy (PSRCHIVE, z wagami jak .spCF) BEZ niego rośnie > 1%.")
    CR = filter(r -> isfinite(num(r, "C")), ok)
    BC = filter(r -> !isempty(strip(get(r, "chans_bad", ""))), CR)
    JC = filter(r -> !isempty(strip(get(r, "chans_cut", ""))), CR)
    san = fin([num(r, "ch_all_none") / num(r, "snr_none") for r in CR])
    isempty(san) || pr(@sprintf("kontrola: suma przez PSRCHIVE vs .spCF (S/N bez wycięć): mediana %.3f", median(san)))
    pr("sprawdzone: $(length(CR));  z kanałami oznaczonymi przez detektor: $(length(BC));  " *
       "z UZASADNIONYM wycięciem: $(length(JC))")
    rs = Dict{String,Vector{Float64}}()
    for r in BC, (t, l) in zip(split(r["chans_bad"]), split(get(r, "chans_loo", "")))
        for why in split(split(t, ":")[2], "+")
            push!(get!(rs, why, Float64[]), parse(Float64, split(l, ":")[2]))
        end
    end
    for (why, v) in sort(collect(rs))
        pr(@sprintf("  reguła %-5s oznaczyła %3d kanałów; S/N sumy bez nich: mediana %+.1f%%, rośnie >1%%: %d",
                    why, length(v), median(v), count(>(1), v)))
    end
    cnt = Dict{String,Int}()
    for r in JC, c in split(r["chans_cut"]); cnt[c] = get(cnt, c, 0) + 1; end
    isempty(cnt) || pr("najczęściej uzasadnione kanały: " *
                       join(["$c ($(n)x)" for (c, n) in sort(collect(cnt), by=x -> -x[2])[1:min(8, end)]], ", "))
    pp = fin([num(r, "pair_pct") for r in CR])
    isempty(pp) || pr(@sprintf("pojedyncze pary (impuls×kanał): mediana %.2f%%, >1%%: %d", median(pp), count(>(1), pp)))
    gain(r) = 100 * (num(r, "ch_good_snr") / num(r, "ch_all_snr") - 1)
    pr(@sprintf("  %-12s %-14s %-24s %13s %24s", "pulsar", "wyciąć", "zmiany S/N bez kanału", "cięć wsz->bez",
                "S/N po automacie wsz->bez"))
    for r in sort(JC, by=r -> -gain(r))[1:min(30, end)]
        pr(@sprintf("  %-12s %-14s %-24s %5d -> %5d   %7.1f -> %7.1f (%+.1f%%)", r["pulsar"], r["chans_cut"],
                    r["chans_loo"], num(r, "ch_all_ncut"), num(r, "ch_good_ncut"), num(r, "ch_all_snr"),
                    num(r, "ch_good_snr"), gain(r)))
    end
    pr("Wykresy z uzasadnieniem dla każdego pulsara: $(joinpath(out, "kanaly"))/<pulsar>.pdf")
    write(joinpath(out, "summary.txt"), String(take!(io)))
    println("\nZapisane: $(joinpath(out, "summary.txt"))  (wklej),  profile nulli: $(joinpath(out, "nulls.pdf"))")
end


"""
    _run_parallel(codes, logdir, tag; every=300, progress=k -> "")

Uruchamia `codes` (kod Julii, jeden proces na element) w tle z `SPATS_NO_MAIN=1`
w katalogu repozytorium, logi w `logdir/<tag>_<k>.log` (zapisywane na bieżąco),
co `every` s wypisuje `progress(k)` i ostatnią istotną linię logu każdej części.
Ctrl+C zatrzymuje procesy w tle. Zwraca, czy wszystkie zakończyły się sukcesem.
"""
function _run_parallel(codes, logdir, tag; every=300, progress=k -> "")
    mkpath(logdir)
    n = length(codes)
    repo = dirname(@__DIR__)
    ios, procs = IO[], Base.Process[]
    for (k, code) in enumerate(codes)
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) -e $code`
        cmd = setenv(cmd, merge(Dict(ENV), Dict("SPATS_NO_MAIN" => "1")); dir=repo)
        io = open(joinpath(logdir, "$(tag)_$k.log"), "w")
        push!(ios, io)
        push!(procs, run(pipeline(cmd, stdout=io, stderr=io), wait=false))
    end
    function lastline(k)
        log = joinpath(logdir, "$(tag)_$k.log")
        L = isfile(log) ? filter(l -> startswith(l, "[") || startswith(l, "WERYFIKACJA") ||
                                      occursin("BŁĄD", l) || occursin("ERROR", l), readlines(log)) : String[]
        return progress(k) * (isempty(L) ? "(ładuje pakiety / liczy pierwszy pulsar)" : last(L))
    end
    try
        t = 0
        while any(process_running, procs)
            sleep(10); t += 10
            if t % every == 0
                println("--- ", Libc.strftime("%H:%M", time()), " ---")
                for k in 1:n
                    println("  część $k/$n: ", first(lastline(k), 170))
                end
            end
        end
    catch e
        e isa InterruptException && (println("Przerywam procesy w tle..."); foreach(kill, procs))
        rethrow()
    finally
        foreach(close, ios)
    end
    ok = true
    for (k, p) in enumerate(procs)
        success(p) && continue
        ok = false
        println("UWAGA: część $k/$n zakończona błędem - zobacz $(joinpath(logdir, "$(tag)_$k.log"))")
    end
    return ok
end

_nproc(n) = something(n, clamp(min(Sys.CPU_THREADS ÷ 4, floor(Int, Sys.free_memory() / 2^30 / 3)), 1, 16))

"""
    zap_verify_parallel(; n=nothing, out=...)

Cała weryfikacja z JEDNEGO miejsca: `n` procesów Julii w tle (każdy
`zap_verify(part="k/n")`, log `out/part_k.log`), postęp co `every` s, na końcu
`nulls.pdf` i `summary.txt` ze wszystkich części. `n` domyślnie z liczby rdzeni
i wolnej pamięci (~3 GB na proces, najwyżej 16). Wznawialne.
"""
function zap_verify_parallel(; n=nothing, out="/home/psr/data/OUTPUT/zap_verify/", every=300, _code=nothing)
    mkpath(out)
    n = _nproc(n)
    println("WERYFIKACJA RÓWNOLEGLE: $n procesów (rdzenie: $(Sys.CPU_THREADS), wolna pamięć: " *
            "$(round(Sys.free_memory() / 2^30, digits=1)) GB)  wyniki: $out")
    codes = [something(_code, "include(\"spats.jl\"); Base.include(Main, \"modules/zap_verify.jl\"); " *
                              "Base.invokelatest(Main.zap_verify; part=\"$k/$n\", out=\"$out\")") for k in 1:n]
    nd(k) = (csv = joinpath(out, "verify_$k.csv"); "gotowe $(isfile(csv) ? countlines(csv) - 1 : 0) | ")
    _run_parallel(codes, out, "part"; every=every, progress=nd)
    isnothing(_code) || return nothing
    try _null_pages(joinpath(out, "nullprof"), joinpath(out, "nulls.pdf")) catch e; println("PDF nulli: $e"); end
    return zap_verify_summary(out=out)
end

"""
    autozap_parallel(; n=nothing, out=pazur2, verify_dir=...)

Ponowny PEŁNY przebieg automatu (`autozap_all`) w `n` procesach w tle, z
poprawkami z weryfikacji: bezpiecznik RFI `:both`, 'nulle' z sygnałem ->
zaszumione, kanały wycięte przed sumowaniem tam, gdzie weryfikacja w
`verify_dir` wykazała zysk >= 1% (`SpaTs.chan_first_list`). Pulsary zrobione
ręcznie (`zap_manual`) zostają. Stara tabela -> `zap_summary_<data>.csv`.
Stan w `out/.autozap_v3`: brak = start od nowa, "running" = wznowienie
(zrobione pomijane), "done" = już zrobione, nic się nie dzieje.
"""
function autozap_parallel(; n=nothing, out="/home/psr/data/OUTPUT/pazur2/", dataroot="/home/psr/data/new/",
                          verify_dir="/home/psr/data/OUTPUT/zap_verify/", every=300)
    marker = joinpath(out, ".autozap_v3")
    state = isfile(marker) ? strip(read(marker, String)) : ""
    state == "done" && (println("AUTOMAT v3: już zrobiony ($marker) - pomijam."); return true)
    summ = joinpath(out, "zap_summary.csv")
    if state != "running"
        stamp = Libc.strftime("%Y%m%d-%H%M", time())
        for f in ("zap_summary.csv", "do_recznie.csv")
            pth = joinpath(out, f)
            isfile(pth) && mv(pth, joinpath(out, replace(f, ".csv" => "_$stamp.csv")))
        end
        open(f -> println(f, SpaTs.AUTOZAP_HEADER), summ, "w")
        open(f -> println(f, "pulsar,powod,obs"), joinpath(out, "do_recznie.csv"), "w")
        write(marker, "running")
    end
    psrs = sort(filter(d -> isdir(joinpath(dataroot, d)), readdir(dataroot)))
    cf = SpaTs.chan_first_list(verify_dir)
    n = _nproc(n)
    println("AUTOMAT v3 RÓWNOLEGLE: $n procesów, $(length(psrs)) pulsarów" *
            (state == "running" ? " (wznowienie)" : "") * ";  kanały przed sumowaniem w $(length(cf)) pulsarach")
    codes = ["include(\"spats.jl\"); Base.invokelatest(SpaTs.autozap_all; only=$(repr(psrs[k:n:end])), " *
             "chan_first=$(repr(cf)))" for k in 1:n]
    tot() = (isfile(summ) ? countlines(summ) - 1 : 0)
    ok = _run_parallel(codes, joinpath(out, "logs_v3"), "automat"; every=every,
                       progress=k -> "razem gotowe $(tot())/$(length(psrs)) | ")
    st = Dict{String,Int}()
    for l in Iterators.drop(eachline(summ), 1)
        r = split(l, ","); length(r) >= 2 && (st[r[2]] = get(st, r[2], 0) + 1)
    end
    println("AUTOMAT v3 gotowy: ", join(["$k=$v" for (k, v) in sort(collect(st))], "  "))
    ok && tot() >= length(psrs) && write(marker, "done")
    return ok
end
