# Animowany sliding LRFS: okno `window` impulsow przesuwane po obserwacji; dla kazdego polozenia pelny LRFS
# (mapa |LRFS|, widmo, faza przy fp, profil okna) + podglad calej obserwacji z zaznaczonym oknem.
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/sliding_lrfs_anim.jl")
#   sliding_lrfs_anim("J1555-3134")              # PDF (strona = polozenie okna, w xpdf strzalki = suwak) + GIF
#   sliding_lrfs_player("J1555-3134")            # interaktywna animacja: Odtwórz/Pauza, suwak, strzałki (VNC)
# Pliki: /home/psr/data/OUTPUT/sliding_lrfs/<psr>/<psr>_lrfs_anim.pdf i .gif
isdefined(Main, :_sl_data) || include(joinpath(@__DIR__, "sliding_lrfs.jl"))
using PyPlot, PyCall
# wycinek Pythona (slice) do indeksowania GridSpec, np. gs[0, 1:4] -> get(gs, (0, _ps(1, 4)))
_ps(a, b) = pybuiltin("slice")(a, b)

"Cechy widma (do 3): z SL_P3S (np. J1901: F1-F3 z TPA XX), inaczej szczyty sredniego widma - pierwszy przy P3
z params (+-10%), potem kolejne najwyzsze lokalne maksima (>= 25% najwyzszego, >= 3 biny od siebie, nie najnizsze
czestotliwosci)."
function _anim_features(ms, freq, window, p, psr; nmax=3)
    haskey(SL_P3S, psr) && return [freq[argmin(abs.(freq .- 1 / q))] for q in SL_P3S[psr]][1:min(nmax, end)]
    nf = length(ms)
    loc = [k for k in 3:nf-1 if ms[k] >= ms[k-1] && ms[k] >= ms[k+1]]
    isempty(loc) && return [freq[argmax(ms)]]
    sort!(loc, by=k -> -ms[k])
    p3 = get(p, "p3", nothing)
    first_ = loc[1]
    if p3 isa Number && p3 > 0
        w = [k for k in loc if abs(freq[k] - 1 / p3) <= max(0.1 / p3, 1 / window)]
        isempty(w) || (first_ = w[1])
    end
    ks = [first_]
    for k in loc
        length(ks) >= nmax && break
        all(abs(k - q) >= 3 for q in ks) && ms[k] >= 0.25 * ms[loc[1]] && push!(ks, k)
    end
    return freq[ks]
end

"Dane klatek: dla kazdego okna |LRFS| (czestotliwosc x dlugosc), widmo, fazy cech, profil okna."
function _anim_prepare(D, p; window=256, step=8, nsig=3.0, psr="")
    PL = SpaTs.Plot
    N, nb = size(D)
    bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
    X, good = PL._demeaned(D)
    on = bst:bend; off = PL._far_offpulse(nb, bst, bend)
    nf = window ÷ 2
    freq = collect(1:nf) ./ window
    starts = collect(1:step:N-window+1)
    A = Array{Float64}(undef, length(starts), nf, length(on))
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
    snr = spec ./ noise
    ms = vec(mean(snr, dims=1))
    feats = _anim_features(ms, freq, window, p, psr)
    # fazy cech w kazdym oknie (wzgledem bezwzglednego numeru impulsu: staly dryf = ten sam wzor), zawiniete do +-180
    Ph = fill(NaN, length(feats), length(starts), length(on))
    for (q, f) in enumerate(feats)
        car = [cis(2π * f * (n - 1)) for n in 1:N]
        C = cumsum(X[:, on] .* car, dims=1); Co = cumsum(X[:, off] .* car, dims=1)
        for (b, s) in enumerate(starts)
            e = s + window - 1
            L = s == 1 ? C[e, :] : C[e, :] .- C[s-1, :]
            Lo = s == 1 ? Co[e, :] : Co[e, :] .- Co[s-1, :]
            σ = std(vcat(real.(Lo), imag.(Lo)))
            ok = abs.(L) .> nsig * σ
            Ph[q, b, ok] .= rad2deg.(angle.(L[ok]))
        end
    end
    profall = vec(mean(Float64.(D[good, on]), dims=1))
    lon = collect(range(-180 * length(on) / nb, 180 * length(on) / nb, length=length(on)))
    lon .-= lon[argmax(profall)]
    return (A=A, Ph=Ph, spec=snr, ms=ms, prof=prof, freq=freq, starts=starts, window=window, lon=lon,
            feats=feats, fp=feats[1], N=N, profall=profall)
end

const _ANIM_COLS = ["royalblue", "green", "crimson"]
const _ANIM_MARK = ["o", "^", "s"]

"Uklad jak w narzedziu GLMakie: podglad calej obserwacji z pasem okna | widmo okna z cechami | faza + mapa LRFS + profil."
function _anim_setup(fig, d, psr; bottom=0.08)
    art = Dict{Symbol,Any}(:psr => psr)
    gs = fig.add_gridspec(6, 12, hspace=0.0, wspace=0.0, left=0.06, right=0.98, bottom=bottom, top=0.95)
    # podglad calej obserwacji (sliding LRFS) z pasem biezacego okna
    axo = fig.add_subplot(get(gs, (_ps(0, 6), _ps(0, 5))))
    axo.imshow(d.spec, origin="lower", aspect="auto", cmap="viridis", interpolation="nearest",
               extent=[0, 0.5, d.starts[1] - 1, d.starts[end] - 1 + d.window], vmin=0, vmax=quantile(vec(d.spec), 0.995))
    axo.set_xlim(0, 0.5); axo.set_title(@sprintf("SLRFS Overview (win=%d)", d.window), fontsize=8)
    axo.set_xlabel("Frequency 1/P₃ (c/P)", fontsize=7); axo.set_ylabel("Pulse number", fontsize=7); axo.tick_params(labelsize=7)
    art[:axo] = axo
    art[:band] = axo.axhspan(d.starts[1] - 1, d.starts[1] - 1 + d.window, color="white", alpha=0.35, lw=0.8, ec="white")
    # legenda cech nad widmem
    axl = fig.add_subplot(get(gs, (0, _ps(5, 7)))); axl.axis("off")
    for (q, f) in enumerate(d.feats)
        axl.scatter([], [], marker=_ANIM_MARK[q], c=_ANIM_COLS[q], label="F$q")
    end
    axl.legend(loc="center", frameon=false, fontsize=7)
    # widmo okna (moc, max = 1) z liniami cech
    axs = fig.add_subplot(get(gs, (_ps(1, 5), _ps(5, 7))))
    ls, = axs.plot(d.spec[1, :] ./ maximum(d.spec[1, :]), d.freq, c="grey", lw=0.9)
    for (q, f) in enumerate(d.feats); axs.axhline(f, c=_ANIM_COLS[q], ls="--", lw=0.7); end
    axs.set_ylim(0, 0.5); axs.set_xlim(0, 1.1); axs.set_xlabel("Power", fontsize=7)
    axs.set_ylabel("Frequency (c/P)", fontsize=7); axs.tick_params(labelsize=7)
    art[:ls] = ls
    # fazy cech w oknie
    axp = fig.add_subplot(get(gs, (0, _ps(7, 12))))
    art[:sc] = [axp.scatter(d.lon, d.Ph[q, 1, :], s=8, marker=_ANIM_MARK[q], c=_ANIM_COLS[q], lw=0) for q in eachindex(d.feats)]
    axp.set_ylim(-185, 185); axp.set_yticks([-180, -90, 0, 90, 180]); axp.set_xlim(d.lon[1], d.lon[end])
    axp.set_ylabel("Phase (deg)", fontsize=7); axp.tick_params(labelbottom=false, labelsize=6)
    # mapa |LRFS| okna
    axm = fig.add_subplot(get(gs, (_ps(1, 5), _ps(7, 12))))
    art[:im] = axm.imshow(d.A[1, :, :], origin="lower", aspect="auto", cmap="viridis", interpolation="nearest",
                          extent=[d.lon[1], d.lon[end], 0, 0.5], vmin=0, vmax=quantile(vec(d.A), 0.998))
    axm.tick_params(labelleft=false, labelbottom=false)
    # profil okna (czarny) i calej obserwacji (szary)
    axb = fig.add_subplot(get(gs, (5, _ps(7, 12))))
    pm = max(maximum(d.profall), eps())
    lp, = axb.plot(d.lon, d.prof[1, :] ./ pm, c="black", lw=1.0)
    axb.plot(d.lon, d.profall ./ pm, c="grey", lw=0.7, ls="--")
    axb.set_xlim(d.lon[1], d.lon[end]); axb.set_ylim(-0.05, 1.1 * maximum(d.prof) / pm)
    axb.set_xlabel("Longitude (°)", fontsize=7); axb.set_ylabel("Intensity", fontsize=7); axb.tick_params(labelsize=6)
    art[:lp] = lp
    art[:status] = fig.text(0.06, 0.01, "", fontsize=7)
    return art
end

function _anim_update!(art, d, b)
    s = d.starts[b]
    art[:band].remove()
    art[:band] = art[:axo].axhspan(s - 1, s - 1 + d.window, color="white", alpha=0.35, lw=0.8, ec="white")
    for (q, sc) in enumerate(art[:sc]); sc.set_offsets(hcat(d.lon, d.Ph[q, b, :])); end
    art[:ls].set_xdata(d.spec[b, :] ./ max(maximum(d.spec[b, :]), eps()))
    art[:im].set_data(d.A[b, :, :])
    art[:lp].set_ydata(d.prof[b, :] ./ max(maximum(d.profall), eps()))
    art[:status].set_text(@sprintf("%s   Window: pulses %d–%d   |   ", art[:psr], s - 1, s + d.window - 2) *
                          join([@sprintf("F%d=%.3f cpp (P3=%.1f)", q, f, 1 / f) for (q, f) in enumerate(d.feats)], "   "))
end

"PDF (strona = polozenie okna; w xpdf strzalki / PgDn = przesuwanie) i GIF dla pulsara."
function sliding_lrfs_anim(psr="J1555-3134"; window=256, step=8, outdir=SL_DIR, gif=true, fps=8)
    D, p = Base.invokelatest(_sl_data, psr)
    d = _anim_prepare(D, p; window=window, step=step, psr=psr)
    od = joinpath(outdir, psr); mkpath(od)
    PdfPages = pyimport("matplotlib.backends.backend_pdf").PdfPages
    pf = joinpath(od, "$(psr)_lrfs_anim.pdf")
    pdf = PdfPages(pf)
    fig = figure(figsize=(12, 6.5))
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
    d = _anim_prepare(D, p; window=window, step=step, psr=psr)
    widgets = pyimport("matplotlib.widgets")
    pygui(true)
    fig = figure(figsize=(13, 7))
    art = _anim_setup(fig, d, psr; bottom=0.17)
    nb = length(d.starts)
    axsl = fig.add_axes([0.16, 0.06, 0.8, 0.03])
    axbt = fig.add_axes([0.06, 0.05, 0.07, 0.045])
    sl = widgets.Slider(axsl, "", 1, nb, valinit=1, valstep=1)
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
