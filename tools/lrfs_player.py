#!/usr/bin/env python3
"""Odtwarzacz animowanego sliding LRFS (dane z sliding_lrfs_export w tools/sliding_lrfs_anim.jl).
    python3 tools/lrfs_player.py /home/psr/data/OUTPUT/sliding_lrfs/J1555-3134/player [--speed 2] [--interval 30]
Odtwórz / Pauza: przycisk albo spacja; w pauzie suwak albo strzałki <- ->; prędkość: + / - (albo strzałki góra / dół).
Okno samo przesuwa się od początku do końca obserwacji; po pauzie Odtwórz kontynuuje od miejsca suwaka."""
import argparse
import json
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button, Slider

ap = argparse.ArgumentParser()
ap.add_argument("dir")
ap.add_argument("--speed", type=int, default=2, help="impulsy na klatkę")
ap.add_argument("--interval", type=int, default=30, help="ms między klatkami")
ap.add_argument("--size", type=float, nargs=2, default=(9.5, 5.2))
ap.add_argument("--dpi", type=int, default=80)
a = ap.parse_args()

m = json.load(open(os.path.join(a.dir, "meta.json")))
nw, nf, non, nlab = m["nw"], m["nf"], m["non"], m["nlab"]


def load(name, shape):
    return np.fromfile(os.path.join(a.dir, name + ".f32"), dtype="<f4").reshape(shape, order="F")   # Julia: kolumnowo


A = load("A", (nw, nf, non))
spec = load("spec", (nw, nf))
Ph = load("Ph", (nlab, nw, non))
Fq = load("Fq", (nlab, nw))
prof = load("prof", (nw, non))
profall = load("profall", (non,))
lon = np.array(m["lon"]); freq = np.array(m["freq"]); starts = np.array(m["starts"]); W = m["window"]; psr = m["psr"]
cols = ["royalblue", "green", "crimson"]; marks = ["o", "^", "s"]

fig = plt.figure(figsize=tuple(a.size), dpi=a.dpi)
gs = GridSpec(6, 12, figure=fig, hspace=0.0, wspace=0.0, left=0.07, right=0.98, bottom=0.17, top=0.94)
# podgląd całej obserwacji + ślady cech + pas bieżącego okna
axo = fig.add_subplot(gs[0:6, 0:5])
y0, y1 = starts[0] - 1, starts[-1] - 1 + W
axo.imshow(spec, origin="lower", aspect="auto", cmap="viridis", interpolation="nearest",
           extent=[0, 0.5, y0, y1], vmin=0, vmax=np.nanquantile(spec, 0.995))
cen = starts - 1 + W / 2
for l in range(nlab):
    axo.scatter(Fq[l], cen, s=1.5, c=cols[l], lw=0, alpha=0.8)
axo.set_xlim(0, 0.5); axo.set_ylim(y0, y1); axo.set_title("SLRFS (okno %d)" % W, fontsize=8)
axo.set_xlabel("Frequency 1/P$_3$ (c/P)", fontsize=7); axo.set_ylabel("Pulse number", fontsize=7); axo.tick_params(labelsize=6)
band = Rectangle((0, starts[0] - 1), 0.5, W, fc="white", ec="white", alpha=0.3, lw=0.8)
axo.add_patch(band)
# legenda cech
axl = fig.add_subplot(gs[0, 5:7]); axl.axis("off")
for l in range(nlab):
    axl.scatter([], [], marker=marks[l], c=cols[l], label="F%d" % (l + 1))
axl.legend(loc="center", frameon=False, fontsize=7)
# widmo okna z liniami cech
axs = fig.add_subplot(gs[1:5, 5:7])
ls, = axs.plot(spec[0] / max(spec[0].max(), 1e-12), freq, c="grey", lw=0.9)
fl = [axs.axhline(0.0, c=cols[l], ls="--", lw=0.8) for l in range(nlab)]
axs.set_ylim(0, 0.5); axs.set_xlim(0, 1.1); axs.set_xlabel("Power", fontsize=7)
axs.set_ylabel("Frequency (c/P)", fontsize=7); axs.tick_params(labelsize=6)
# fazy cech
axp = fig.add_subplot(gs[0, 7:12])
sc = [axp.scatter(lon, Ph[l, 0], s=7, marker=marks[l], c=cols[l], lw=0) for l in range(nlab)]
axp.set_ylim(-185, 185); axp.set_yticks([-180, 0, 180]); axp.set_xlim(lon[0], lon[-1])
axp.set_ylabel("Phase (deg)", fontsize=7); axp.tick_params(labelbottom=False, labelsize=6)
# mapa |LRFS|
axm = fig.add_subplot(gs[1:5, 7:12])
im = axm.imshow(A[0], origin="lower", aspect="auto", cmap="viridis", interpolation="nearest",
                extent=[lon[0], lon[-1], 0, 0.5], vmin=0, vmax=np.nanquantile(A, 0.998))
axm.tick_params(labelleft=False, labelbottom=False)
# profil okna i całej obserwacji
axb = fig.add_subplot(gs[5, 7:12])
pm = max(profall.max(), 1e-12)
lp, = axb.plot(lon, prof[0] / pm, c="black", lw=1.0)
axb.plot(lon, profall / pm, c="grey", lw=0.7, ls="--")
axb.set_xlim(lon[0], lon[-1]); axb.set_ylim(-0.05, 1.1 * prof.max() / pm)
axb.set_xlabel("Longitude (deg)", fontsize=7); axb.set_ylabel("Intensity", fontsize=7); axb.tick_params(labelsize=6)
status = fig.text(0.07, 0.01, "", fontsize=7)
# sterowanie
axsl = fig.add_axes([0.18, 0.065, 0.78, 0.03])
axbt = fig.add_axes([0.07, 0.055, 0.09, 0.05])
sl = Slider(axsl, "", 1, nw, valinit=1, valstep=1)
bt = Button(axbt, "Odtwórz")
st = {"play": False, "i": 0, "speed": max(1, a.speed)}


def draw(i):
    s = starts[i]
    band.set_y(s - 1)
    for l in range(nlab):
        sc[l].set_offsets(np.c_[lon, Ph[l, i]])
        f = Fq[l, i]
        fl[l].set_visible(bool(np.isfinite(f)))
        if np.isfinite(f):
            fl[l].set_ydata([f, f])
    ls.set_xdata(spec[i] / max(spec[i].max(), 1e-12))
    im.set_data(A[i])
    lp.set_ydata(prof[i] / pm)
    fs = ["F%d=%.3f cpp (P3=%.1f)" % (l + 1, Fq[l, i], 1 / Fq[l, i]) for l in range(nlab) if np.isfinite(Fq[l, i])]
    status.set_text("%s   okno: impulsy %d-%d   |   %s   |   prędkość %d imp./klatkę" %
                    (psr, s - 1, s + W - 2, "   ".join(fs) if fs else "brak wyraźnych cech", st["speed"]))


def on_slider(v):
    st["i"] = int(round(v)) - 1
    if not st["play"]:
        draw(st["i"]); fig.canvas.draw_idle()


sl.on_changed(on_slider)


def toggle(_=None):
    if st["play"]:
        st["play"] = False; bt.label.set_text("Odtwórz")
    else:
        if st["i"] >= nw - 1:
            st["i"] = 0
        st["play"] = True; bt.label.set_text("Pauza")
    fig.canvas.draw_idle()


bt.on_clicked(toggle)


def on_key(ev):
    k = ev.key
    if k == " ":
        toggle()
    elif k in ("+", "=", "up"):
        st["speed"] = min(32, 2 * st["speed"])
    elif k in ("-", "down"):
        st["speed"] = max(1, st["speed"] // 2)
    elif k in ("right", "left") and not st["play"]:
        sl.set_val(min(nw, max(1, st["i"] + 1 + (1 if k == "right" else -1))))
    draw(st["i"]); fig.canvas.draw_idle()


fig.canvas.mpl_connect("key_press_event", on_key)


def step(_):
    if st["play"]:
        st["i"] = min(nw - 1, st["i"] + st["speed"])
        if st["i"] >= nw - 1:
            st["play"] = False; bt.label.set_text("Odtwórz")
        sl.eventson = False; sl.set_val(st["i"] + 1); sl.eventson = True
        draw(st["i"])
    return []


draw(0)
try:
    ani = FuncAnimation(fig, step, interval=a.interval, blit=False, cache_frame_data=False)
except TypeError:                                   # starsze matplotlib bez cache_frame_data
    ani = FuncAnimation(fig, step, interval=a.interval, blit=False)
print("Odtwórz / Pauza: przycisk albo spacja; w pauzie suwak albo strzałki <- ->; prędkość: + / -. Zamknij okno, żeby skończyć.")
plt.show()
