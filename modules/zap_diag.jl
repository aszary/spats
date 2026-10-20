# Diagnostyka automatu zapow na serwerze (dane nie wychodza z kopernika).
# W tym samym REPL, w ktorym byl include("spats.jl"):
#     include("modules/zap_diag.jl"); zap_diag(["J1750-3503", "J0459-0210"])
# Wypisuje krotkie tabelki (do wklejenia) i zapisuje <psr>/pulsar_zapdiag.png.
# Numery impulsow od 0 (jak w pazi).
using Statistics
isdefined(Main, :SpaTs) || error("Najpierw include(\"spats.jl\") (moze byc przerwane Ctrl+C po starcie testu).")
const _D, _T, _P = SpaTs.Data, SpaTs.Tools, SpaTs.Plot

function zap_diag(psrs; vpmout="/home/psr/output/", ranges=nothing)
    for psr in psrs
        outdir = vpmout * psr
        p = _T.read_params(joinpath(outdir, "params.json"))
        d = _D.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        N, nb = size(d)
        bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
        on = bst:bend
        off = vcat(1:bst-1, bend+1:nb)
        sn = _D.pulse_signal_noise(d, bst, bend)
        a = _D.auto_zap_core(sn)
        g = sn.good
        nz = _D.noise_jump_z(sn.vwhite, g)
        med(x) = median(x[g])
        # per pulse
        sig = [std(@view d[i, off]) for i in 1:N]
        base = [mean(@view d[i, off]) for i in 1:N]
        basez = base ./ (med(sig) / sqrt(length(off)))            # level of the whole off-pulse row
        ac1 = [(x = d[i, off] .- base[i]; sum(x[1:end-1] .* x[2:end]) / max(sum(abs2, x), eps())) for i in 1:N]
        spike = [maximum(abs.(@view(d[i, off]) .- base[i])) / max(sig[i], eps()) for i in 1:N]
        h = length(off) ÷ 2
        slope = [(mean(@view d[i, off[1:h]]) - mean(@view d[i, off[h+1:end]])) / (med(sig) / sqrt(h / 2)) for i in 1:N]
        vr, vv = sn.vwhite ./ med(sn.vwhite), sn.v ./ med(sn.v)
        snr1 = sn.s ./ sqrt.(max.(sn.vwhite, eps()))
        r(x, k=2) = round(x, digits=k)

        println("\n================ $psr  N=$N nbin=$nb on=$bst-$bend (" *
                "$(length(on)) binów)  okna off: $(sn.nwin) użyte / $(sn.nleak) z emisją ================")
        println("automat: $(length(a.bad)) wyciętych  \"$(_D.pulses_to_zaps(a.bad))\"")
        println("typowo: σ_off=$(r(med(sig), 4))  ac1=$(r(med(ac1), 3)) (biały szum ≈ 0±$(r(1 / sqrt(length(off)), 3)))  " *
                "S/N 1 impulsu=$(r(med(snr1)))")
        rg = isnothing(ranges) ? [((k - 1) * N ÷ 4 + 1):(k * N ÷ 4) for k in 1:4] : ranges
        println("\nodcinek (od 0)   szum_biały  szum_okien   ac1    |poziom|z  S/N1  frac(szum>1.5)  frac(|poziom|z>5)  wycięte")
        for q in rg
            k = filter(i -> g[i], q)
            println(rpad("$(first(q)-1)-$(last(q)-1)", 16),
                    lpad(r(median(vr[k])), 10), lpad(r(median(vv[k])), 12), lpad(r(median(ac1[k]), 3), 8),
                    lpad(r(median(abs.(basez[k])), 1), 10), lpad(r(median(snr1[k]), 1), 7),
                    lpad(r(mean(vr[k] .> 1.5)), 14), lpad(r(mean(abs.(basez[k]) .> 5)), 18),
                    lpad(count(in(Set(a.bad)), q), 10))
        end
        for (name, x) in (("szum biały (×typowy)", vr), ("szum okien (×typowy)", vv), ("|poziom tła| z", abs.(basez)),
                          ("ac1 (paski/falowanie)", ac1), ("spike max/σ", spike), ("nachylenie z", abs.(slope)),
                          ("skok szumu nz", nz))
            idx = sort(findall(g), by=i -> -x[i])[1:min(15, count(g))]
            println("\nTOP $name: " * join(["$(i-1):$(r(x[i], 1))" for i in idx], "  "))
        end

        # figure: metrics vs pulse number, cut pulses shaded
        plt = _P.PyPlot
        fig, ax = plt.subplots(6, 1, figsize=(11, 12), sharex=true)
        x0 = 0:N-1
        for (k, (lab, y)) in enumerate((("szum biały ×typ", vr), ("szum okien ×typ", vv), ("poziom tła z", basez),
                                        ("ac1", ac1), ("spike/σ", spike), ("S/N 1 imp.", snr1)))
            ax[k].plot(x0, y, lw=0.6, color="k")
            ax[k].set_ylabel(lab, fontsize=8)
            for b in a.bad; ax[k].axvspan(b - 1.5, b - 0.5, color="tab:red", alpha=0.15, lw=0); end
        end
        ax[6].set_xlabel("impuls (od 0)   czerwone tło = wycięte przez automat")
        fig.suptitle(psr)
        fig.tight_layout()
        png = joinpath(outdir, "pulsar_zapdiag.png")
        fig.savefig(png, dpi=90)
        plt.close(fig)
        println("\nwykres: $png")
    end
end
