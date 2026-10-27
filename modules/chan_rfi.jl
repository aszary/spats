# Detektor RFI w kanalach (impuls x kanal) - TYLKO diagnostyka, archiwa nietkniete.
# Rekomendacja nr 1 z przegladu narzedzi RFI: statystyki wylacznie w off-pulse
# (nulle i jasne impulsy nie moga byc "anomalia"), odporne z = (x - mediana) /
# (IQR/1.349) liczone w KAZDYM kanale osobno wzdluz impulsow (logika clfd/TPA).
# W tym samym REPL, w ktorym byl include("spats.jl"):
#     include("modules/chan_rfi.jl"); chan_rfi(["J1110-5637"])
# Wypisuje tabelke kanalow (do wklejenia), zapisuje pazur2/<psr>/chan_rfi.csv
# (oflagowane pary) i chan_rfi.png. Numery impulsow i kanalow od 0 (jak w pazi/paz).
using Statistics
isdefined(Main, :SpaTs) || error("Najpierw include(\"spats.jl\") (moze byc przerwane Ctrl+C po starcie).")
const _CD, _CT, _CP = SpaTs.Data, SpaTs.Tools, SpaTs.Plot

"Kostka (impuls, kanal, bin) z `pdv -t` archiwum o jednej polaryzacji."
function _load_cube(txt)
    nsub = nch = nb = 0
    data = nothing
    for line in eachline(txt)
        r = split(line)
        isempty(r) && continue
        if !isdigit(r[1][1])
            if r[1] == "File:" && isnothing(data)
                nsub, nch, nb = parse(Int, r[6]), parse(Int, r[8]), parse(Int, r[12])
                data = zeros(Float32, nsub, nch, nb)
            end
            continue
        end
        data[parse(Int, r[1]) + 1, parse(Int, r[2]) + 1, parse(Int, r[3]) + 1] = parse(Float32, r[4])
    end
    return data
end

robz_t(x, g) = (m = median(x[g]); s = max((quantile(x[g], 0.75) - quantile(x[g], 0.25)) / 1.349, eps()); (x .- m) ./ s)
# cechy rosnace przy zakloceniach (std, ptp, |ACF1|): poziom i skala z CICHEJ
# czesci rozkladu (kwantyle 10/25%) - kanal zaszumiony w 35% impulsow
# przesuwa mediane i IQR, a dolnego ogona nie
function robz_lo(x, g)
    q10, q25 = quantile(x[g], 0.1), quantile(x[g], 0.25)
    s = max((q25 - q10) / 0.6071, eps())
    return (x .- (q25 + 0.6745s)) ./ s
end

"Rdzeń `chan_rfi` na kostce `d` (impuls, kanał, bin); indeksy od 1."
function chan_rfi_core(d, bst, bend; q=4.0, nmin=2, fchan=0.3, fsub=0.5, margin=0.05)
    N, C, nb = size(d)
    # off-pulse = poza oknem on-pulse z marginesem (5% obrotu)
    mg = round(Int, margin * nb)
    off = filter(j -> j < bst - mg || j > bend + mg, 1:nb)

    F = zeros(N, C, 4); fz = zeros(N, C)
    for i in 1:N, c in 1:C
        x = Float64.(@view d[i, c, off])
        m = mean(x); y = x .- m; v = sum(abs2, y)
        F[i, c, 1] = std(y)
        F[i, c, 2] = maximum(y) - minimum(y)
        F[i, c, 3] = v > 0 ? abs(sum(y[1:end-1] .* y[2:end]) / v) : 0.0
        F[i, c, 4] = m
        fz[i, c] = count(iszero, x) / length(x)
    end
    empty = [all(iszero, @view d[i, :, :]) for i in 1:N]       # juz wyciete / puste impulsy
    dead = [all(iszero, @view d[.!empty, c, :]) for c in 1:C]  # martwe kanaly
    g = .!empty
    Z = zeros(N, C, 4)
    for c in 1:C
        dead[c] && continue
        gc = g .& (fz[:, c] .<= 0.5)
        count(gc) < 10 && continue
        for k in 1:3; Z[:, c, k] = robz_lo(F[:, c, k], gc); end
        Z[:, c, 4] = robz_t(F[:, c, 4], gc)
    end
    # kanal z paskami przez cala obserwacje: wzdluz czasu nic nie odstaje, ale
    # |ACF1| (niezalezne od wzmocnienia kanalu) jest wyzsze niz w innych kanalach
    acfc = [dead[c] ? NaN : median(F[g, c, 3]) for c in 1:C]
    lv = findall(.!dead)
    ma = median(acfc[lv]); sa = max(1.4826 * median(abs.(acfc[lv] .- ma)), 1 / sqrt(length(off)))
    acf_bad = [!dead[c] && (acfc[c] - ma) / sa > q for c in 1:C]
    flag = falses(N, C)
    for i in findall(g), c in findall(.!dead)
        flag[i, c] = count(abs.(Z[i, c, :]) .> q) >= nmin || fz[i, c] > 0.5
    end
    live = count(.!dead)
    fch = [count(flag[g, c]) / max(count(g), 1) for c in 1:C]
    chan_all = findall(((fch .> fchan) .| acf_bad) .& .!dead)
    fsb = [count(flag[i, .!dead]) / max(live, 1) for i in 1:N]
    sub_all = findall(g .& (fsb .> fsub))
    pairs = [(i, c) for i in findall(g) for c in 1:C if flag[i, c] && !(c in chan_all) && !(i in sub_all)]

    return (F=F, Z=Z, fz=fz, flag=flag, g=g, dead=dead, fch=fch, chan_all=chan_all, sub_all=sub_all,
            pairs=pairs, off=off)
end

"""
    chan_rfi(psrs; q=4.0, nmin=2, fchan=0.3, fsub=0.5, ...)

Para (impuls, kanal) oflagowana, gdy >= `nmin` z cech off-pulse (std, ptp,
|ACF1|, srednia) ma |z| > `q`, albo ponad polowa probek to zera.
Klasy: kanal_caly (> `fchan` impulsow w kanale oflagowanych - jak w TPA: te
same kanaly we wszystkich impulsach), impuls_caly (> `fsub` kanalow), para.
"""
function chan_rfi(psrs; dataroot="/home/psr/data/new/", params_root="/home/psr/data/OUTPUT/pazur2/",
                  q=4.0, nmin=2, fchan=0.3, fsub=0.5, margin=0.05, plot=true)
    res = Dict{String,Any}()
    for psr in psrs
        pdir = joinpath(params_root, psr)
        pf = joinpath(pdir, "params.json")
        isfile(pf) || (println("$psr: brak $pf - pomijam"); continue)
        p = _CT.read_params(pf)
        psr_path = joinpath(dataroot, psr)
        obs = sort(filter(d -> isdir(joinpath(psr_path, d)), readdir(psr_path)))
        isempty(obs) && (println("$psr: brak danych"); continue)
        tmp = mktempdir()
        ar = joinpath(tmp, "p16.ar")
        _CD.add_psrfiles(joinpath(psr_path, obs[1]) * "/", ar; sixteen=true)
        try run(pipeline(`pam -p -m $ar`, stderr=devnull)) catch; end   # jedna polaryzacja
        txt = joinpath(tmp, "p16.txt")
        run(pipeline(`pdv -t $ar`, stdout=txt, stderr=devnull))
        d = _load_cube(txt)
        rm(tmp, recursive=true, force=true)
        isnothing(d) && (println("$psr: nie wczytano danych z pdv"); continue)
        N, C, nb = size(d)

        nb != 1024 && println("UWAGA: $nb binów, okno on-pulse z params zakłada 1024")
        a = chan_rfi_core(d, Int(p["bin_st"]), Int(p["bin_end"]); q=q, nmin=nmin, fchan=fchan, fsub=fsub, margin=margin)
        F, Z, fz, flag, g, dead, fch, chan_all, sub_all, pairs, off = a.F, a.Z, a.fz, a.flag, a.g, a.dead,
            a.fch, a.chan_all, a.sub_all, a.pairs, a.off
        empty, live = .!g, count(.!dead)

        # porownanie z wycietymi przez automat calymi impulsami (zaps2)
        mine = Set(_CD.zaps_to_rows(_CD.zaps_of(p), N))
        np_in = count(t -> t[1] in mine, pairs)
        r(x, k=2) = round(x, digits=k)
        println("\n================ $psr  N=$N kanały=$C biny=$nb  off-pulse: $(length(off)) binów  " *
                "q=$q nmin=$nmin ================")
        println("puste impulsy: $(count(empty))   martwe kanały: $(join(findall(dead) .- 1, " "))")
        println("kanał  oflag%   med_std    med|ACF1|  klasa")
        for c in 1:C
            dead[c] && continue
            gc = g .& (fz[:, c] .<= 0.5)
            ms = count(gc) > 0 ? median(F[gc, c, 1]) : NaN
            ma = count(gc) > 0 ? median(F[gc, c, 3]) : NaN
            println(lpad(c - 1, 4), lpad(r(100 * fch[c], 1), 9), lpad(r(ms, 4), 11), lpad(r(ma, 3), 11),
                    c in chan_all ? "   KANAŁ DO WYCIĘCIA" : "")
        end
        println("kanały do wycięcia w całej obserwacji (paz -z): \"$(join(chan_all .- 1, " "))\"")
        println("całe impulsy (> $(Int(100fsub))% kanałów): $(length(sub_all))  " *
                "(z tego już w zaps2 automatu: $(count(in(mine), sub_all)))  \"$(_CD.pulses_to_zaps(sub_all))\"")
        println("pojedyncze pary (impuls, kanał): $(length(pairs))  " *
                "($(r(100 * length(pairs) / max(count(g) * live, 1), 3))% par; w impulsach już wyciętych: $np_in)")

        open(joinpath(pdir, "chan_rfi.csv"), "w") do f
            println(f, "isub,ichan,klasa,z_std,z_ptp,z_acf1,z_mean,frac_zero")
            rows = vcat([(i, c, "kanal_caly") for c in chan_all for i in findall(g)],
                        [(i, c, "impuls_caly") for i in sub_all for c in findall(.!dead)],
                        [(i, c, "para") for (i, c) in pairs])
            for (i, c, k) in rows
                println(f, join([i - 1, c - 1, k, r(Z[i, c, 1]), r(Z[i, c, 2]), r(Z[i, c, 3]), r(Z[i, c, 4]),
                                 r(fz[i, c], 3)], ","))
            end
        end
        if plot
            plt = _CP.PyPlot
            fig, ax = plt.subplots(1, 3, figsize=(15, 9), sharey=true)
            zs = clamp.(maximum(abs.(Z), dims=3)[:, :, 1], 0, 2q)
            ax[1].imshow(zs, aspect="auto", origin="lower", cmap="magma", interpolation="nearest")
            ax[1].set_title("max |z| cech off-pulse (obcięte do $(2q))")
            ax[2].imshow(flag, aspect="auto", origin="lower", cmap="Greys", interpolation="nearest")
            ax[2].set_title("oflagowane pary; kanały do wycięcia: $(join(chan_all .- 1, " "))")
            ax[3].imshow(dropdims(sum(d, dims=2), dims=2), aspect="auto", origin="lower", cmap="viridis",
                         interpolation="nearest")
            ax[3].set_title("dane zsumowane w częstotliwości")
            ax[1].set_ylabel("impuls"); ax[1].set_xlabel("kanał"); ax[2].set_xlabel("kanał"); ax[3].set_xlabel("bin")
            fig.suptitle("$psr — detektor RFI w kanałach (tylko off-pulse)")
            fig.tight_layout()
            png = joinpath(pdir, "chan_rfi.png")
            fig.savefig(png, dpi=80); plt.close(fig)
            println("wykres: $png")
        end
        res[psr] = (chan_all=chan_all .- 1, sub_all=sub_all .- 1, n_pairs=length(pairs), frac_chan=fch)
    end
    return res
end
