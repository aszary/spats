# Animowany sliding LRFS: okno `window` impulsow przesuwane po obserwacji; dla kazdego polozenia pelny LRFS
# (mapa |LRFS|, widmo, faza przy fp, profil okna) + podglad calej obserwacji z zaznaczonym oknem.
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/sliding_lrfs_anim.jl")
#   sliding_lrfs_anim("J1555-3134")              # PDF (strona = polozenie okna, w xpdf strzalki = suwak) + GIF
#   sliding_lrfs_player("J1555-3134")            # interaktywna animacja: Odtwórz/Pauza, suwak, strzałki (VNC)
# Pliki: /home/psr/data/OUTPUT/sliding_lrfs/<psr>/<psr>_lrfs_anim.pdf i .gif
isdefined(Main, :_sl_data) || include(joinpath(@__DIR__, "sliding_lrfs.jl"))
using PyPlot, PyCall

"Dane klatek: dla kazdego okna |LRFS| (czestotliwosc x dlugosc), widmo, faza przy fp, profil okna."
function _anim_prepare(D, p; window=256, step=8, nsig=3.0)
    PL = SpaTs.Plot
    N, nb = size(D)
    bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
    X, good = PL._demeaned(D)
    on = bst:bend; off = PL._far_offpulse(nb, bst, bend)
    nf = window ÷ 2
    freq = collect(1:nf) ./ window
    starts = collect(1:step:N-window+1)
    A = Array{Float64}(undef, length(starts), nf, length(on))       # |LRFS| w oknie
    Ph = fill(NaN, length(starts), length(on))                        # faza przy fp (nasza konwencja)
    spec = zeros(length(starts), nf); noise = zeros(length(starts))
    prof = zeros(length(starts), length(on))
    for (b, s) in enumerate(starts)
        r = s:s+window-1
        B = X[r, on]; B .-= mean(B, dims=1)
        Bo = X[r, off]; Bo .-= mean(Bo, dims=1)
        F = conj.(rfft(B, 1))[2:nf+1, :]; Fo = conj.(rfft(Bo, 1))[2:nf+1, :]
        A[b, :, :] = abs.(F)
        spec[b, :] = vec(sum(abs2.(F), dims=2))
        noise[b] = mean(sum(abs2.(Fo), dims=2)) * length(on) / length(off)
        gr = [i for i in r if good[i]]
        prof[b, :] = isempty(gr) ? zeros(length(on)) : vec(mean(Float64.(D[gr, on]), dims=1))
    end
    # fp: szczyt sredniego widma przy P3 z params (+-10%), jak w pozostalych wykresach
    ms = vec(mean(spec ./ noise, dims=1))
    kp = 1 + argmax(ms[2:end])
    p3 = get(p, "p3", nothing)
    if p3 isa Number && p3 > 0
        w = [k for k in 1:nf if abs(freq[k] - 1 / p3) <= max(0.1 / p3, 1 / window)]
        isempty(w) || (kp = w[argmax(ms[w])])
    end
    # faza przy fp w kazdym oknie (wzgledem bezwzglednego numeru impulsu: staly dryf = ten sam wzor)
    car = [cis(2π * freq[kp] * (n - 1)) for n in 1:N]
    C = cumsum(X[:, on] .* car, dims=1); Co = cumsum(X[:, off] .* car, dims=1)
    for (b, s) in enumerate(starts)
        e = s + window - 1
        L = s == 1 ? C[e, :] : C[e, :] .- C[s-1, :]
        Lo = s == 1 ? Co[e, :] : Co[e, :] .- Co[s-1, :]
        σ = std(vcat(real.(Lo), imag.(Lo)))
        ok = abs.(L) .> nsig * σ
        ph = rad2deg.(angle.(L))
        jr = argmax(vec(mean(Float64.(D[good, on]), dims=1)))
        ph .-= ph[jr]; ph = mod.(ph .+ 180, 360) .- 180
        Ph[b, ok] .= ph[ok]
    end
    lon = collect(range(-180 * length(on) / nb, 180 * length(on) / nb, length=length(on)))
    jr = argmax(vec(mean(Float64.(D[good, on]), dims=1))); lon .-= lon[jr]
    return (A=A, Ph=Ph, spec=spec ./ noise, ms=ms, prof=prof, freq=freq, starts=starts, window=window,
            lon=lon, fp=freq[kp], N=N, profall=vec(mean(Float64.(D[good, on]), dims=1)))
end

"Buduje osie i elementy wykresu raz; klatki zmieniaja tylko dane (_anim_update!) - plynna animacja."
function _anim_setup(fig, d, psr; bottom=0.08)
    vmax = quantile(vec(d.A), 0.998)
    ttl = fig.suptitle("", fontsize=9)
    gs = fig.add_gridspec(5, 6, hspace=0.0, wspace=0.0, left=0.08, right=0.92, bottom=bottom, top=0.9)
    # faza przy fp w oknie (nasza konwencja, 0 w szczycie profilu)
    axp = fig.add_subplot(get(gs, (0, pyslice(1, 4))))
    sc = axp.scatter(d.lon, d.Ph[1, :], s=6, c="black")
    axp.set_ylim(-190, 190); axp.set_yticks([-180, 0, 180]); axp.set_xlim(d.lon[1], d.lon[end])
    axp.set_ylabel("faza (°)", fontsize=7); axp.tick_params(labelbottom=false, labelsize=7)
    # widmo okna (czarne) i srednie z calej obserwacji (szare)
    axs = fig.add_subplot(get(gs, (pyslice(1, 4), 0)))
    ls, = axs.plot(d.spec[1, :], d.freq, c="black", lw=0.7)
    axs.plot(d.ms, d.freq, c="grey", lw=0.6, ls="--")
    axs.set_ylim(0, 0.5); axs.set_xlim(1.05 * maximum(d.spec), 0); axs.set_ylabel("częstotliwość (1/P)", fontsize=7)
    axs.tick_params(labelsize=7); axs.set_xticks([])
    # mapa |LRFS| w oknie
    axm = fig.add_subplot(get(gs, (pyslice(1, 4), pyslice(1, 4))))
    im = axm.imshow(d.A[1, :, :], origin="lower", aspect="auto", cmap="viridis", interpolation="nearest",
                    extent=[d.lon[1], d.lon[end], d.freq[1], d.freq[end]], vmin=0, vmax=vmax)
    axm.axhline(d.fp, c="red", lw=0.5, ls=":")
    axm.tick_params(labelleft=false, labelbottom=false)
    # profil okna (czarny) i calej obserwacji (szary)
    axb = fig.add_subplot(get(gs, (4, pyslice(1, 4))))
    pm = max(maximum(d.profall), eps())
    lp, = axb.plot(d.lon, d.prof[1, :] ./ pm, c="black", lw=0.8)
    axb.plot(d.lon, d.profall ./ pm, c="grey", lw=0.6, ls="--")
    axb.set_xlim(d.lon[1], d.lon[end]); axb.set_ylim(-0.1, 1.05 * maximum(d.prof) / pm)
    axb.set_xlabel("długość (°), 0 = szczyt profilu", fontsize=7); axb.tick_params(labelsize=7); axb.set_yticks([])
    # podglad calej obserwacji z zaznaczonym oknem
    axo = fig.add_subplot(get(gs, (pyslice(0, 5), pyslice(4, 6))))
    axo.imshow(d.spec, origin="lower", aspect="auto", cmap="viridis", interpolation="nearest",
               extent=[0, 0.5, d.starts[1] - 1, d.starts[end] - 1], vmin=0, vmax=quantile(vec(d.spec), 0.995))
    hl = axo.axhline(d.starts[1] - 1, c="red", lw=1.2)
    axo.set_xlabel("częstotliwość (1/P)", fontsize=7); axo.yaxis.set_label_position("right"); axo.yaxis.tick_right()
    axo.set_ylabel("początek okna (impuls)", fontsize=7); axo.tick_params(labelsize=7); axo.set_title("cała obserwacja", fontsize=7)
    return (ttl=ttl, sc=sc, ls=ls, im=im, lp=lp, hl=hl, psr=psr)
end

function _anim_update!(art, d, b)
    s = d.starts[b]
    art.ttl.set_text(@sprintf("%s: impulsy %d–%d (okno %d), fp = %.4f (P3 = %.2f P)", art.psr, s - 1, s + d.window - 2,
                              d.window, d.fp, 1 / d.fp))
    art.sc.set_offsets(hcat(d.lon, d.Ph[b, :]))
    art.ls.set_xdata(d.spec[b, :])
    art.im.set_data(d.A[b, :, :])
    art.lp.set_ydata(d.prof[b, :] ./ max(maximum(d.profall), eps()))
    art.hl.set_ydata([s - 1, s - 1])
end

"PDF (strona = polozenie okna; w xpdf strzalki / PgDn = przesuwanie) i GIF dla pulsara."
function sliding_lrfs_anim(psr="J1555-3134"; window=256, step=8, outdir=SL_DIR, gif=true, fps=8)
    D, p = Base.invokelatest(_sl_data, psr)
    d = _anim_prepare(D, p; window=window, step=step)
    od = joinpath(outdir, psr); mkpath(od)
    PdfPages = pyimport("matplotlib.backends.backend_pdf").PdfPages
    pf = joinpath(od, "$(psr)_lrfs_anim.pdf")
    pdf = PdfPages(pf)
    fig = figure(figsize=(9, 6))
    art = _anim_setup(fig, d, psr)
    pngs = String[]; tmp = mktempdir()
    for b in eachindex(d.starts)
        _anim_update!(art, d, b)
        pdf.savefig(fig)
        gif && (f = joinpath(tmp, @sprintf("f%04d.png", b)); fig.savefig(f, dpi=70); push!(pngs, f))
    end
    pdf.close(); close(fig)
    println("PDF: $pf  ($(length(d.starts)) stron; w xpdf strzałki / PgDn = przesuwanie okna)")
    if gif
        try
            Image = pyimport("PIL.Image")
            fr = [Image.open(f) for f in pngs]
            gf = joinpath(od, "$(psr)_lrfs_anim.gif")
            fr[1].save(gf, save_all=true, append_images=fr[2:end], duration=round(Int, 1000 / fps), loop=0)
            println("GIF: $gf")
        catch e
            println("GIF pominięty (brak Pillow?): ", first(sprint(showerror, e), 120))
        end
    end
    rm(tmp; recursive=true, force=true)
    return pf
end

"""
    sliding_lrfs_player(psr; window=256, step=2, interval=120)

Interaktywna animacja w jednym oknie (wymaga ekranu, np. VNC): przycisk Odtwórz / Pauza (albo spacja) -
okno samo przesuwa sie od poczatku do konca obserwacji; w pauzie suwak przesuwa okno recznie, strzalki <- ->
o jeden krok; Odtwórz po pauzie kontynuuje od miejsca suwaka, na koncu zaczyna od poczatku.
`interval` = czas jednej klatki w ms, `step` = przesuniecie okna miedzy klatkami (impulsy).
"""
function sliding_lrfs_player(psr="J1555-3134"; window=256, step=2, interval=120)
    D, p = Base.invokelatest(_sl_data, psr)
    d = _anim_prepare(D, p; window=window, step=step)
    widgets = pyimport("matplotlib.widgets")
    pygui(true)
    fig = figure(figsize=(10, 7))
    art = _anim_setup(fig, d, psr; bottom=0.13)
    nb = length(d.starts)
    axsl = fig.add_axes([0.22, 0.025, 0.55, 0.03])
    axbt = fig.add_axes([0.06, 0.02, 0.12, 0.045])
    sl = widgets.Slider(axsl, "okno", 1, nb, valinit=1, valstep=1)
    bt = widgets.Button(axbt, "Odtwórz")
    state = Dict(:playing => false)
    timer = fig.canvas.new_timer(interval=interval)
    show_frame(v) = (_anim_update!(art, d, clamp(round(Int, v), 1, nb)); fig.canvas.draw_idle())
    sl.on_changed(show_frame)
    function tick()
        state[:playing] || return
        v = round(Int, sl.val)
        if v >= nb
            state[:playing] = false; timer.stop(); bt.label.set_text("Odtwórz"); fig.canvas.draw_idle()
        else
            sl.set_val(v + 1)                                  # przesuwa suwak i rysuje klatke
        end
    end
    timer.add_callback(tick)
    function toggle(_...)
        if state[:playing]
            state[:playing] = false; timer.stop(); bt.label.set_text("Odtwórz")
        else
            round(Int, sl.val) >= nb && sl.set_val(1)          # od poczatku, jesli doszlo do konca
            state[:playing] = true; timer.start(); bt.label.set_text("Pauza")
        end
        fig.canvas.draw_idle()
    end
    bt.on_clicked(toggle)
    function onkey(ev)
        k = ev.key
        if k == " "
            toggle()
        elseif k in ("right", "left") && !state[:playing]
            sl.set_val(clamp(round(Int, sl.val) + (k == "right" ? 1 : -1), 1, nb))
        end
    end
    fig.canvas.mpl_connect("key_press_event", onkey)
    # referencje, zeby Python nie usunal suwaka / przycisku / timera
    global _SL_PLAYER = (fig=fig, sl=sl, bt=bt, timer=timer, tick=tick, toggle=toggle, onkey=onkey)
    show()
    println("Odtwórz / Pauza: przycisk albo spacja; w pauzie suwak albo strzałki ← →. Zamknij okno, żeby skończyć.")
    return fig
end
