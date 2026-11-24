# Animowany sliding LRFS: okno `window` impulsow przesuwane po obserwacji; dla kazdego polozenia pelny LRFS
# (mapa |LRFS|, widmo, fazy cech, profil okna) + podglad calej obserwacji z zaznaczonym oknem.
# Cechy (szczyty widma, do 3) wykrywane OSOBNO W KAZDYM OKNIE i sledzone w czasie (ten sam kolor, dopoki szczyt
# jest blisko poprzedniego) - linie w widmie i fazy zmieniaja sie razem z danymi.
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/sliding_lrfs_anim.jl")
#   sliding_lrfs_player("J1555-3134")            # interaktywna animacja: Odtwórz/Pauza, suwak, strzałki (VNC)
#   sliding_lrfs_anim("J1555-3134")              # PDF (strona = polozenie okna, w xpdf strzalki = suwak) + GIF
# Pliki: /home/psr/data/OUTPUT/sliding_lrfs/<psr>/<psr>_lrfs_anim.pdf i .gif
isdefined(Main, :_sl_data) || include(joinpath(@__DIR__, "sliding_lrfs.jl"))
using PyPlot, PyCall
# wycinek Pythona (slice) do indeksowania GridSpec, np. gs[0, 1:4] -> get(gs, (0, _ps(1, 4)))
_ps(a, b) = pybuiltin("slice")(a, b)

const _ANIM_COLS = ["royalblue", "green", "crimson"]
const _ANIM_MARK = ["o", "^", "s"]

"""
Szczyty widma jednego okna: lokalne maksima (bez 2 najnizszych czestotliwosci - powolne zmiany, nulle), moc/szum
>= `thr` i >= `rel` najwyzszego, co najmniej `sep` binow od silniejszych i oddzielone od nich wyraznym dolkiem
(minimum widma miedzy nimi < `dip` x wysokosc slabszego) - inaczej to zbocze tego samego szczytu (J1555: "cecha"
0.109 obok 0.094); najwyzej `nmax`, od najsilniejszego.
"""
function _anim_peaks(sp; thr=5.0, rel=0.3, sep=6, dip=0.7, nmax=3)
    nf = length(sp)
    loc = [k for k in 3:nf-1 if sp[k] >= sp[k-1] && sp[k] >= sp[k+1]]
    isempty(loc) && return Int[]
    top = maximum(sp[loc])
    sort!(loc, by=k -> -sp[k])
    ks = Int[]
    for k in loc
        length(ks) >= nmax && break
        sp[k] >= thr && sp[k] >= rel * top || continue
        all(abs(k - q) >= sep && minimum(sp[min(k, q):max(k, q)]) < dip * sp[k] for q in ks) && push!(ks, k)
    end
    return ks
end

"""
Sledzenie cech w kolejnych oknach: szczyt dostaje etykiete (1-3) cechy, ktorej czestotliwosc w poprzednim
oknie byla najblizej (do `tol` binow); nowy szczyt - pierwsza wolna etykieta. Zwraca K[etykieta, okno] (bin albo 0).
"""
function _anim_track(peaks; nlab=3, tol=4)
    K = zeros(Int, nlab, length(peaks))
    prev = zeros(Int, nlab)
    for (b, ks) in enumerate(peaks)
        used = falses(nlab)
        for k in ks                                         # od najsilniejszego
            cand = [l for l in 1:nlab if !used[l] && prev[l] > 0 && abs(prev[l] - k) <= tol]
            l = isempty(cand) ? findfirst(l -> !used[l] && prev[l] == 0, 1:nlab) : cand[argmin([abs(prev[q] - k) for q in cand])]
            isnothing(l) && (l = findfirst(!, used))
            isnothing(l) && continue
            K[l, b] = k; used[l] = true
        end
        prev = [K[l, b] > 0 ? K[l, b] : (b > 1 && prev[l] > 0 && any(K[l, max(1, b - 5):b] .> 0) ? prev[l] : 0) for l in 1:nlab]
    end
    return K
end

"Dane klatek: |LRFS| okna, widmo, cechy okna (sledzone) z fazami, profil okna."
function _anim_prepare(D, p; window=256, step=8, nsig=3.0, thr=5.0)
    PL = SpaTs.Plot
    N, nb = size(D)
    bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
    X, good = PL._demeaned(D)
    on = bst:bend; off = PL._far_offpulse(nb, bst, bend)
    nf = window ÷ 2
    freq = collect(1:nf) ./ window
    starts = collect(1:step:N-window+1)
    nw = length(starts)
    A = Array{Float64}(undef, nw, nf, length(on))
    Fc = Array{ComplexF64}(undef, nw, nf, length(on))          # zespolone LRFS (faza wzgledem poczatku obserwacji)
    σo = zeros(nw, nf)
    spec = zeros(nw, nf); noise = zeros(nw)
    prof = zeros(nw, length(on))
    for (b, s) in enumerate(starts)
        r = s:s+window-1
        B = X[r, on]; B .-= mean(B, dims=1)
        Bo = X[r, off]; Bo .-= mean(Bo, dims=1)
        F = conj.(rfft(B, 1))[2:nf+1, :]; Fo = conj.(rfft(Bo, 1))[2:nf+1, :]    # conj: kernel e^(+2πi f n)
        # faza wzgledem bezwzglednego numeru impulsu: staly dryf daje ten sam wzor w kazdym oknie
        F .*= cis.(2π .* freq .* (s - 1))
        A[b, :, :] = abs.(F); Fc[b, :, :] = F
        σo[b, :] = [std(vcat(real.(Fo[k, :]), imag.(Fo[k, :]))) for k in 1:nf]
        spec[b, :] = vec(sum(abs2.(F), dims=2))
        noise[b] = mean(sum(abs2.(Fo), dims=2)) * length(on) / length(off)
        gr = [i for i in r if good[i]]
        prof[b, :] = isempty(gr) ? zeros(length(on)) : vec(mean(Float64.(D[gr, on]), dims=1))
    end
    snr = spec ./ noise
    K = _anim_track([_anim_peaks(snr[b, :]; thr=thr) for b in 1:nw])
    nlab = size(K, 1)
    Ph = fill(NaN, nlab, nw, length(on)); Fq = fill(NaN, nlab, nw)
    for l in 1:nlab, b in 1:nw
        k = K[l, b]; k == 0 && continue
        Fq[l, b] = freq[k]
        c = Fc[b, k, :]
        ok = abs.(c) .> nsig * σo[b, k]
        Ph[l, b, ok] .= rad2deg.(angle.(c[ok]))
    end
    profall = vec(mean(Float64.(D[good, on]), dims=1))
    lon = collect(range(-180 * length(on) / nb, 180 * length(on) / nb, length=length(on)))
    lon .-= lon[argmax(profall)]
    return (A=A, Ph=Ph, Fq=Fq, spec=snr, ms=vec(mean(snr, dims=1)), prof=prof, freq=freq, starts=starts,
            window=window, lon=lon, N=N, profall=profall)
end

"Uklad: podglad calej obserwacji z pasem okna i sladami cech | widmo okna z liniami cech | fazy + mapa LRFS + profil."
function _anim_setup(fig, d, psr; bottom=0.08, animated=false)
    # animated=true (odtwarzacz): ruchome elementy rysowane tylko przy blittingu - plynna animacja
    art = Dict{Symbol,Any}(:psr => psr)
    nlab = size(d.Fq, 1)
    gs = fig.add_gridspec(6, 12, hspace=0.0, wspace=0.0, left=0.07, right=0.98, bottom=bottom, top=0.94)
    # podglad calej obserwacji (sliding LRFS) + slady cech + pas biezacego okna
    axo = fig.add_subplot(get(gs, (_ps(0, 6), _ps(0, 5))))
    y0 = d.starts[1] - 1; y1 = d.starts[end] - 1 + d.window
    axo.imshow(d.spec, origin="lower", aspect="auto", cmap="viridis", interpolation="nearest",
               extent=[0, 0.5, y0, y1], vmin=0, vmax=quantile(vec(d.spec), 0.995))
    cen = d.starts .- 1 .+ d.window / 2
    for l in 1:nlab
        axo.scatter(d.Fq[l, :], cen, s=1.5, c=_ANIM_COLS[l], lw=0, alpha=0.8)
    end
    axo.set_xlim(0, 0.5); axo.set_ylim(y0, y1); axo.set_title(@sprintf("SLRFS (okno %d)", d.window), fontsize=8)
    axo.set_xlabel("Frequency 1/P₃ (c/P)", fontsize=7); axo.set_ylabel("Pulse number", fontsize=7); axo.tick_params(labelsize=6)
    art[:band] = pyimport("matplotlib.patches").Rectangle((0, d.starts[1] - 1), 0.5, d.window, fc="white", ec="white",
                                                          alpha=0.3, lw=0.8, animated=animated)
    axo.add_patch(art[:band])
    # legenda cech
    axl = fig.add_subplot(get(gs, (0, _ps(5, 7)))); axl.axis("off")
    for l in 1:nlab; axl.scatter([], [], marker=_ANIM_MARK[l], c=_ANIM_COLS[l], label="F$l"); end
    axl.legend(loc="center", frameon=false, fontsize=7)
    # widmo okna (moc, max = 1) z liniami cech tego okna
    axs = fig.add_subplot(get(gs, (_ps(1, 5), _ps(5, 7))))
    ls, = axs.plot(d.spec[1, :] ./ maximum(d.spec[1, :]), d.freq, c="grey", lw=0.9, animated=animated)
    art[:fl] = [axs.axhline(0.0, c=_ANIM_COLS[l], ls="--", lw=0.8, animated=animated) for l in 1:nlab]
    axs.set_ylim(0, 0.5); axs.set_xlim(0, 1.1); axs.set_xlabel("Power", fontsize=7)
    axs.set_ylabel("Frequency (c/P)", fontsize=7); axs.tick_params(labelsize=6)
    art[:ls] = ls
    # fazy cech w oknie
    axp = fig.add_subplot(get(gs, (0, _ps(7, 12))))
    art[:sc] = [axp.scatter(d.lon, d.Ph[l, 1, :], s=7, marker=_ANIM_MARK[l], c=_ANIM_COLS[l], lw=0, animated=animated)
                for l in 1:nlab]
    axp.set_ylim(-185, 185); axp.set_yticks([-180, 0, 180]); axp.set_xlim(d.lon[1], d.lon[end])
    axp.set_ylabel("Phase (deg)", fontsize=7); axp.tick_params(labelbottom=false, labelsize=6)
    # mapa |LRFS| okna
    axm = fig.add_subplot(get(gs, (_ps(1, 5), _ps(7, 12))))
    art[:im] = axm.imshow(d.A[1, :, :], origin="lower", aspect="auto", cmap="viridis", interpolation="nearest",
                          extent=[d.lon[1], d.lon[end], 0, 0.5], vmin=0, vmax=quantile(vec(d.A), 0.998), animated=animated)
    axm.tick_params(labelleft=false, labelbottom=false)
    # profil okna (czarny) i calej obserwacji (szary)
    axb = fig.add_subplot(get(gs, (5, _ps(7, 12))))
    pm = max(maximum(d.profall), eps())
    lp, = axb.plot(d.lon, d.prof[1, :] ./ pm, c="black", lw=1.0, animated=animated)
    axb.plot(d.lon, d.profall ./ pm, c="grey", lw=0.7, ls="--")
    axb.set_xlim(d.lon[1], d.lon[end]); axb.set_ylim(-0.05, 1.1 * maximum(d.prof) / pm)
    axb.set_xlabel("Longitude (°)", fontsize=7); axb.set_ylabel("Intensity", fontsize=7); axb.tick_params(labelsize=6)
    art[:lp] = lp
    art[:status] = fig.text(0.07, 0.01, "", fontsize=7, animated=animated)
    art[:dyn] = vcat([art[:band], ls], art[:fl], art[:sc], [art[:im], lp, art[:status]])    # ruchome elementy
    return art
end

function _anim_update!(art, d, b)
    s = d.starts[b]
    art[:band].set_y(s - 1)
    for (l, sc) in enumerate(art[:sc]); sc.set_offsets(hcat(d.lon, d.Ph[l, b, :])); end
    for (l, fl) in enumerate(art[:fl])
        f = d.Fq[l, b]
        fl.set_visible(isfinite(f)); isfinite(f) && fl.set_ydata([f, f])
    end
    art[:ls].set_xdata(d.spec[b, :] ./ max(maximum(d.spec[b, :]), eps()))
    art[:im].set_data(d.A[b, :, :])
    art[:lp].set_ydata(d.prof[b, :] ./ max(maximum(d.profall), eps()))
    fs = [@sprintf("F%d=%.3f cpp (P3=%.1f)", l, d.Fq[l, b], 1 / d.Fq[l, b]) for l in axes(d.Fq, 1) if isfinite(d.Fq[l, b])]
    art[:status].set_text(@sprintf("%s   okno: impulsy %d–%d   |   ", art[:psr], s - 1, s + d.window - 2) *
                          (isempty(fs) ? "brak wyraźnych cech" : join(fs, "   ")))
end

"PDF (strona = polozenie okna; w xpdf strzalki / PgDn = przesuwanie) i GIF dla pulsara."
function sliding_lrfs_anim(psr="J1555-3134"; window=256, step=8, outdir=SL_DIR, gif=true, fps=8)
    D, p = Base.invokelatest(_sl_data, psr)
    d = _anim_prepare(D, p; window=window, step=step)
    od = joinpath(outdir, psr); mkpath(od)
    PdfPages = pyimport("matplotlib.backends.backend_pdf").PdfPages
    pf = joinpath(od, "$(psr)_lrfs_anim.pdf")
    pdf = PdfPages(pf)
    fig = figure(figsize=(11, 6))
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
    sliding_lrfs_player(psr; window=256, speed=2, interval=30, size=(9.5, 5.2), dpi=80, thr=5.0)

Interaktywna animacja w jednym oknie (wymaga ekranu, np. VNC): Odtwórz / Pauza (przycisk albo spacja) - okno
samo przesuwa sie od poczatku do konca obserwacji; w pauzie suwak przesuwa okno recznie, strzalki <- -> o jeden
krok; Odtwórz po pauzie kontynuuje od miejsca suwaka, na koncu zaczyna od poczatku. Predkosc: klawisze + / -
(albo strzalki gora / dol) - ile impulsow okno przesuwa sie na klatke (`speed`, 1-32). `interval` = przerwa
miedzy klatkami (ms). Animacje napedza petla w Julii (timer okna graficznego przez PyCall/VNC bywa martwy);
funkcja trzyma konsole, dopoki okno jest otwarte - zamknij okno, zeby wrocic.
"""
function sliding_lrfs_player(psr="J1555-3134"; window=256, speed=2, interval=30, size=(9.5, 5.2), dpi=80, thr=5.0)
    D, p = Base.invokelatest(_sl_data, psr)
    d = _anim_prepare(D, p; window=window, step=1, thr=thr)
    widgets = pyimport("matplotlib.widgets"); plt = pyimport("matplotlib.pyplot")
    pygui(true)
    fig = figure(figsize=size, dpi=dpi)
    art = _anim_setup(fig, d, psr; bottom=0.17, animated=true)
    nb = length(d.starts)
    axsl = fig.add_axes([0.18, 0.065, 0.78, 0.03])
    axbt = fig.add_axes([0.07, 0.055, 0.09, 0.05])
    sl = widgets.Slider(axsl, "", 1, nb, valinit=1, valstep=1)
    sl.drawon = false                                          # suwak nie przerysowuje calej figury
    bt = widgets.Button(axbt, "Odtwórz")
    state = Dict{Symbol,Any}(:playing => false, :speed => clamp(speed, 1, 32))
    # blitting: tlo (bez ruchomych elementow) zapamietane po pelnym rysowaniu; klatka = tlo + ruchome elementy
    bg = Ref{Any}(nothing)
    function blit_frame()
        c = fig.canvas
        isnothing(bg[]) && return c.draw_idle()
        c.restore_region(bg[])
        for a in art[:dyn]; fig.draw_artist(a); end
        fig.draw_artist(axsl); fig.draw_artist(axbt)
        c.blit(fig.bbox)
    end
    function ondraw(_)
        bg[] = fig.canvas.copy_from_bbox(fig.bbox)
        for a in art[:dyn]; fig.draw_artist(a); end
        fig.canvas.blit(fig.bbox)
    end
    fig.canvas.mpl_connect("draw_event", ondraw)
    _anim_update!(art, d, 1)
    show_frame(v) = (_anim_update!(art, d, clamp(round(Int, v), 1, nb)); blit_frame())
    sl.on_changed(show_frame)
    function toggle(_...)
        if state[:playing]
            state[:playing] = false; bt.label.set_text("Odtwórz")
        else
            round(Int, sl.val) >= nb && sl.set_val(1)          # od poczatku, jesli doszlo do konca
            state[:playing] = true; bt.label.set_text("Pauza")
        end
        blit_frame()
    end
    bt.on_clicked(toggle)
    function onkey(ev)
        k = ev.key
        if k == " "
            toggle()
        elseif k in ("+", "=", "up")
            state[:speed] = min(32, 2 * state[:speed]); println("prędkość: $(state[:speed]) imp./klatkę")
        elseif k in ("-", "down")
            state[:speed] = max(1, state[:speed] ÷ 2); println("prędkość: $(state[:speed]) imp./klatkę")
        elseif k in ("right", "left") && !state[:playing]
            sl.set_val(clamp(round(Int, sl.val) + (k == "right" ? 1 : -1), 1, nb))
        end
    end
    fig.canvas.mpl_connect("key_press_event", onkey)
    global _SL_PLAYER = (fig=fig, sl=sl, bt=bt, toggle=toggle, onkey=onkey, ondraw=ondraw, show_frame=show_frame)
    show(); fig.canvas.draw()
    println("Odtwórz / Pauza: przycisk albo spacja; w pauzie suwak albo strzałki ← →; prędkość: + / - (↑ ↓).")
    println("Zamknij okno, żeby wrócić do konsoli.")
    # petla animacji: przesuwa okno i obsluguje zdarzenia okna (klikniecia, klawisze) miedzy klatkami
    while plt.fignum_exists(fig.number)
        if state[:playing]
            v = round(Int, sl.val)
            if v >= nb
                state[:playing] = false; bt.label.set_text("Odtwórz"); blit_frame()
            else
                sl.set_val(min(nb, v + state[:speed]))         # przesuwa suwak i rysuje klatke (blitting)
            end
        end
        fig.canvas.start_event_loop(interval / 1000)          # zdarzenia okna bez pelnego przerysowania
    end
    return nothing
end
