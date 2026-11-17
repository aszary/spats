#!/usr/bin/env python3
"""Frakcja nullingu z modelu mieszaniny energii (Keith+2026, TPA XIX, wzor 2-3)
i porownanie z nullami wycietymi przez automat (zaps2_null).

Wejscie: energie z tools/nulling_energies.jl (/home/psr/data/OUTPUT/nulling_model/energies/).
Model:   f(E_on) = NF * f_off(E) + (1-NF) * (f_int (*) f_off)(E),  f_off = N(mu0, s0),
         lnL = sum ln f(E_on) + sum ln f_off(E_off)   (szum wiazany z E_off).
         f_int: G, LN, G+G, G+LN, LN+LN, LN+PL; wybor modelu: BIC.
Roznica wobec Keith+2026: najwieksza wiarygodnosc (Nelder-Mead, kilka startow) +
bootstrap po impulsach zamiast nested sampling (dynesty); szum gaussowski, nie empiryczny.
Impulsy wyciete przez automat jako RFI sa pomijane; nulle automatu zostaja w danych.
Etykiety: 0 zostaje, 1 RFI, 2 null automatu, 3 inne wyciete (wygaszona emisja, marginesy; v4).

    python3 tools/nulling_model.py              # wszystkie pliki energii, rownolegle
    python3 tools/nulling_model.py J0034-0721   # wybrane
    python3 tools/nulling_model.py --test       # test na danych syntetycznych ze znanym NF
    NULLING_OUT=/home/psr/data/OUTPUT/nulling_model_v4 python3 tools/nulling_model.py   # inny katalog
"""
import glob
import json
import math
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy import optimize, special

OUT = os.environ.get("NULLING_OUT", "/home/psr/data/OUTPUT/nulling_model").rstrip("/")
NBOOT = 30
KMED = 0.5          # mediana kazdego skladnika emisji >= KMED x szum energii impulsu (s0)
MODEL_VERSION = 2   # zmiana modelu -> stare wyniki z fits/ nie sa uzywane
MODELS = [("G",), ("LN",), ("G", "G"), ("G", "LN"), ("LN", "LN"), ("LN", "PL")]
NPAR = {"G": 2, "LN": 2, "PL": 2}
LIT_NF = {"J0034-0721": 44.0, "J1049-5833": 47.0, "J1255-6131": 75.0, "J1537-4912": 20.0, "J1727-2739": 52.0,
          "J1850-0006": 48.0, "J1851-0053": 45.0, "J1926-0652": 75.0, "J1931+1439": 42.0, "J1944+1755": 60.0,
          "J1946+1805": 55.0, "J2037+1942": 44.0, "J2253+1516": 49.0}
SQ2PI = math.sqrt(2 * math.pi)


# ------------------------------------------------------------------ model
def comp_pdf(kind, q, x):
    """gestosc skladnika intrinsic na siatce x (E w jednostkach <E_on>)."""
    if kind == "G":
        mu, sig = q[0], math.exp(q[1])
        return np.exp(-0.5 * ((x - mu) / sig) ** 2) / (sig * SQ2PI)
    out = np.zeros_like(x)
    if kind == "LN":
        m, s = q[0], math.exp(q[1])
        pos = x > 0
        out[pos] = np.exp(-0.5 * ((np.log(x[pos]) - m) / s) ** 2) / (x[pos] * s * SQ2PI)
        return out
    if kind == "PL":
        xmin, a = math.exp(q[0]), 1.0 + math.exp(q[1])
        pos = x >= xmin
        out[pos] = (a - 1) / xmin * (x[pos] / xmin) ** (-a)
        return out
    raise ValueError(kind)


def unpack(model, th):
    mu0, s0, nf = th[0], math.exp(th[1]), special.expit(th[2])
    k = 3
    comps = []
    for kind in model:
        comps.append((kind, th[k:k + NPAR[kind]]))
        k += NPAR[kind]
    nc = len(model)
    w = np.ones(1) if nc == 1 else special.softmax(np.concatenate([[0.0], th[k:k + nc - 1]]))
    return mu0, s0, nf, comps, w


def prior_ok(comps, s0=0.0):
    """Keith+2026: mu_N > 2.5 sigma_N. Do tego mediana kazdego skladnika emisji >= KMED*s0:
    skladnik z emisja ~0 po splocie z szumem wyglada jak null i przejmuje jego wage
    (J0034-0721, J1931+1439, J1901-0312: NF = 0 przy NF z literatury ~40%)."""
    for kind, q in comps:
        med = q[0] if kind == "G" else math.exp(q[0])
        if med < KMED * s0:
            return False
        if kind == "G" and q[0] < 2.5 * math.exp(q[1]):        # Keith+2026: mu_N > 2.5 sigma_N
            return False
        if kind == "LN" and not (-8 < q[0] < 4 and -4 < q[1] < 2):
            return False
        if kind == "PL" and not (-8 < q[0] < 3 and math.exp(q[1]) < 6):
            return False
    return True


def f_on(model, th, E):
    """gestosc E_on oraz skladnik null (do P_null)."""
    mu0, s0, nf, comps, w = unpack(model, th)
    lo, hi = min(E.min(), 0.0) - 6 * s0, E.max() + 6 * s0
    dx = max(s0 / 6, (hi - lo) / 20000)
    x = np.arange(lo, hi + dx, dx)
    fint = sum(wi * comp_pdf(kind, q, x) for wi, (kind, q) in zip(w, comps))
    norm = fint.sum() * dx
    if not np.isfinite(norm) or norm <= 0:
        return None, None
    fint = fint / norm
    hk = int(math.ceil(6 * s0 / dx))
    t = np.arange(-hk, hk + 1) * dx
    ker = np.exp(-0.5 * (t / s0) ** 2) / (s0 * SQ2PI)
    conv = np.convolve(fint, ker, mode="same") * dx
    fem = np.interp(E - mu0, x, conv, left=0.0, right=0.0)
    fnull = np.exp(-0.5 * ((E - mu0) / s0) ** 2) / (s0 * SQ2PI)
    return nf * fnull + (1 - nf) * fem, nf * fnull


def nll(th, model, Eon, Eoff):
    try:
        _, s0u, _, comps, _ = unpack(model, th)
        if not prior_ok(comps, s0u):
            return 1e12
        f, _ = f_on(model, th, Eon)
        if f is None:
            return 1e12
        mu0, s0 = th[0], math.exp(th[1])
        loff = -0.5 * ((Eoff - mu0) / s0) ** 2 - math.log(s0 * SQ2PI)
        v = -(np.sum(np.log(np.maximum(f, 1e-300))) + np.sum(loff))
        return v if np.isfinite(v) else 1e12
    except (FloatingPointError, OverflowError, ValueError):
        return 1e12


def starts(model, Eon, Eoff, nf0s=(0.05, 0.3, 0.6)):
    mu0, s0 = float(np.mean(Eoff)), float(np.std(Eoff))
    em = Eon[Eon > 3 * s0 + mu0]
    if len(em) < 10:
        em = Eon[Eon > np.median(Eon)]
    lo, hi = np.percentile(em, [30, 80]) if len(em) > 3 else (0.7, 1.5)
    sd = max(float(np.std(em)), s0)
    guess = {"G": [float(np.mean(em)), math.log(max(min(sd, np.mean(em) / 3), 1e-3))],
             "LN": [math.log(max(np.median(em), 1e-3)), math.log(0.5)],
             "PL": [math.log(max(hi, 1e-3)), math.log(1.5)]}
    out = []
    for nf0 in nf0s:
        th = [mu0, math.log(max(s0, 1e-6)), float(special.logit(nf0))]
        for j, kind in enumerate(model):
            g = list(guess[kind])
            if len(model) > 1 and kind in ("G", "LN"):             # rozsun skladniki
                c = lo if j == 0 else hi
                g[0] = c if kind == "G" else math.log(max(c, 1e-3))
                if kind == "G":
                    g[1] = math.log(max(min(sd / 2, c / 3), 1e-3))
            th += g
        th += [0.0] * (len(model) - 1)
        out.append(np.array(th))
    return out


def fit(model, Eon, Eoff, x0s=None):
    best = None
    for x0 in (x0s or starts(model, Eon, Eoff)):
        r = optimize.minimize(nll, x0, args=(model, Eon, Eoff), method="Nelder-Mead",
                              options={"maxiter": 6000, "maxfev": 9000, "xatol": 1e-5, "fatol": 1e-4, "adaptive": True})
        r = optimize.minimize(nll, r.x, args=(model, Eon, Eoff), method="Nelder-Mead",
                              options={"maxiter": 3000, "xatol": 1e-6, "fatol": 1e-5, "adaptive": True})
        if best is None or r.fun < best.fun:
            best = r
    return best


# ---------------------------------------------------------------- dane i analiza
def load(path):
    meta, rows = {}, []
    for line in open(path):
        if line.startswith("#"):
            k, v = line[1:].strip().split("=", 1)
            meta[k.strip()] = v.strip()
        elif not line.startswith("i,"):
            rows.append([float(x) for x in line.split(",")])
    a = np.array(rows)
    return meta, a[:, 0].astype(int), a[:, 1], a[:, 2], a[:, 3].astype(int)


def analyse(path, nboot=NBOOT, seed=1):
    meta, idx, Eon, Eoff, lab = load(path)
    psr = meta["pulsar"]
    keep = lab != 1                                   # bez RFI wycietego przez automat
    Eon_k, Eoff_k, lab_k, idx_k = Eon[keep], Eoff[keep], lab[keep], idx[keep]
    norm = float(np.mean(Eon_k))                      # <E_on> razem z nullami
    if not norm > 0:
        return {"pulsar": psr, "error": "<E_on> <= 0"}
    e_on, e_off = Eon_k / norm, Eoff_k / norm
    n = len(e_on)
    res = {"pulsar": psr, "N": int(len(Eon)), "n_fit": n, "n_rfi_auto": int((lab == 1).sum()),
           "n_null_auto": int((lab == 2).sum()), "control": False}
    # NF z progu 3 sigma (Bhattacharyya+2010): E_on < 3 sqrt(n_on) sigma_bin
    sE = math.sqrt(float(meta["n_on"])) * float(meta["sigma_bin"]) / norm
    nn3 = int((e_on < 3 * sE).sum())
    res.update(NF_3sigma=nn3 / n, NF_3sigma_err=math.sqrt(nn3) / n, NF_auto=res["n_null_auto"] / n,
               s0_over_sE=float(np.std(e_off)) / sE)
    # modele
    fits = {}
    for model in MODELS:
        r = fit(model, e_on, e_off)
        k = len(r.x)
        fits["+".join(model)] = {"nll": float(r.fun), "k": k, "bic": float(2 * r.fun + k * math.log(2 * n)),
                                 "x": r.x.tolist(), "model": model}
    bestname = min(fits, key=lambda m: fits[m]["bic"])
    b = fits[bestname]
    for m in fits:
        fits[m]["dbic"] = fits[m]["bic"] - b["bic"]
    th = np.array(b["x"])
    mu0, s0, nf, comps, w = unpack(b["model"], th)
    f, fn = f_on(b["model"], th, e_on)
    pnull = fn / np.maximum(f, 1e-300)
    res.update(model="N+" + bestname, NF_model=nf, mu0=mu0, s0=s0,
               dbic={("N+" + m): round(fits[m]["dbic"], 2) for m in fits},
               nf_by_model={("N+" + m): round(float(special.expit(fits[m]["x"][2])), 4) for m in fits})
    # bootstrap po impulsach (pary on/off razem), najlepszy model; dlugie obserwacje
    # (np. J1057-5226, ~22 tys. impulsow) - mniej losowan, inaczej liczy sie godzinami
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(nboot if n <= 5000 else min(nboot, 8)):
        j = rng.integers(0, n, n)
        rb = fit(b["model"], e_on[j], e_off[j], x0s=[th])
        bs.append(float(special.expit(rb.x[2])))
    if bs:
        res.update(NF_model_lo=float(np.percentile(bs, 16)), NF_model_hi=float(np.percentile(bs, 84)))
    # porownanie impuls po impulsie z automatem
    an = lab_k == 2
    mn = pnull >= 0.5
    res.update(n_null_model=int(mn.sum()), ambiguous_frac=float(((pnull > 0.05) & (pnull < 0.95)).mean()),
               auto_null_with_p_lt_005=int((an & (pnull < 0.05)).sum()),
               auto_null_with_p_ge_05=int((an & mn).sum()),
               model_null_cut_by_auto=int((mn & an).sum()), model_null_not_cut=int((mn & ~an).sum()),
               mean_E_auto_nulls=float(e_on[an].mean()) if an.any() else float("nan"),
               mean_E_model_nulls=float(e_on[mn].mean()) if mn.any() else float("nan"))
    res["NF_lit"] = LIT_NF.get(psr)
    # dane do wykresow
    os.makedirs(os.path.join(OUT, "fits"), exist_ok=True)
    np.savez(os.path.join(OUT, "fits", psr + ".npz"), idx=idx_k, e_on=e_on, e_off=e_off, lab=lab_k, pnull=pnull,
             th=th, model=np.array(b["model"]), sE=sE)
    return res


def plot(psr, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    z = np.load(os.path.join(OUT, "fits", psr + ".npz"))
    model = tuple(str(m) for m in z["model"])
    th = z["th"]
    e_on, e_off, lab, pn, idx = z["e_on"], z["e_off"], z["lab"], z["pnull"], z["idx"]
    mu0, s0, nf, comps, w = unpack(model, th)
    fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    a = ax[0, 0]
    lo, hi = min(e_on.min(), e_off.min()), e_on.max()
    bins = np.linspace(lo, hi, 80)
    a.hist(e_off, bins=bins, density=True, color="0.7", label="E_off")
    a.hist(e_on, bins=bins, density=True, histtype="step", color="k", label="E_on")
    xx = np.linspace(lo, hi, 800)
    f, fn = f_on(model, th, xx)
    a.plot(xx, f, "C0", label="model N+%s" % "+".join(model))
    a.plot(xx, fn, "C1", ls="--", label="składnik null (NF=%.3f)" % nf)
    a.axvline(3 * float(z["sE"]), color="C3", ls=":", label="3σ_E")
    a.set_xlabel("E / <E_on>"); a.set_yscale("log"); a.set_ylim(1e-3, None); a.legend(fontsize=8)
    a.set_title("%s: NF model %.3f, 3σ %.3f, automat %.3f" % (psr, nf, res["NF_3sigma"], res["NF_auto"]), fontsize=10)
    a = ax[0, 1]
    a.hist(np.log10(np.clip(pn, 1e-6, 1)), bins=60, color="C0")
    a.set_xlabel("log10 P_null"); a.set_title("P_null; niejednoznaczne (0.05–0.95): %.0f%%" % (100 * res["ambiguous_frac"]),
                                             fontsize=10)
    a = ax[1, 0]
    m = min(500, len(e_on))
    a.plot(idx[:m], e_on[:m], ".", ms=3, color="0.4")
    a.plot(idx[:m][pn[:m] >= 0.5], e_on[:m][pn[:m] >= 0.5], "o", mfc="none", color="C1", ms=6, label="null wg modelu")
    a.plot(idx[:m][lab[:m] == 2], e_on[:m][lab[:m] == 2], "x", color="C3", ms=4, label="null wg automatu")
    a.axhline(0, color="k", lw=0.5); a.set_xlabel("impuls"); a.set_ylabel("E / <E_on>"); a.legend(fontsize=8)
    a.set_title("sekwencja energii (pierwsze %d impulsów bez RFI)" % m, fontsize=10)
    a = ax[1, 1]
    an = lab == 2
    a.hist([pn[an], pn[~an]], bins=np.linspace(0, 1, 21), stacked=True, color=["C3", "0.7"],
           label=["wycięte przez automat jako null", "pozostałe"])
    a.set_xlabel("P_null"); a.set_yscale("log"); a.legend(fontsize=8)
    a.set_title("czy nulle automatu to nulle wg modelu", fontsize=10)
    fig.tight_layout()
    os.makedirs(os.path.join(OUT, "plots"), exist_ok=True)
    fig.savefig(os.path.join(OUT, "plots", psr + ".pdf"))
    plt.close(fig)


def work(path):
    psr = os.path.basename(path)[:-4]
    cache = os.path.join(OUT, "fits", "%s.v%d.json" % (psr, MODEL_VERSION))   # gotowy pulsar: wynik z pliku
    if os.path.isfile(cache):
        return json.load(open(cache, encoding="utf-8"))
    try:
        nlines = sum(1 for _ in open(path))
        if nlines > 5000:
            print("  %s: długa obserwacja (%d impulsów) - to potrwa dłużej" % (psr, nlines - 9), flush=True)
        res = analyse(path)
    except Exception as e:
        res = {"pulsar": psr, "error": repr(e)}
    if "error" not in res:
        try:
            plot(psr, res)
        except Exception as e:                              # np. brak matplotlib: wynik zostaje
            res["plot_error"] = repr(e)
        with open(cache, "w", encoding="utf-8") as fh:
            json.dump(res, fh)
    print("  %s: %s" % (psr, res.get("error") or "NF model %.3f [%s] 3σ %.3f automat %.3f  model %s" % (
        res["NF_model"], "%.3f-%.3f" % (res.get("NF_model_lo", float("nan")), res.get("NF_model_hi", float("nan"))),
        res["NF_3sigma"], res["NF_auto"], res["model"])), flush=True)
    return res


def report(results, control):
    ok = [r for r in results if "error" not in r]
    for r in ok:
        r["control"] = r["pulsar"] in control
    cols = ["pulsar", "control", "N", "n_fit", "n_rfi_auto", "n_null_auto", "n_null_model", "NF_auto", "NF_3sigma",
            "NF_3sigma_err", "NF_model", "NF_model_lo", "NF_model_hi", "NF_lit", "model", "ambiguous_frac",
            "auto_null_with_p_lt_005", "auto_null_with_p_ge_05", "model_null_not_cut", "mean_E_auto_nulls",
            "mean_E_model_nulls", "s0_over_sE"]
    with open(os.path.join(OUT, "summary.csv"), "w", encoding="utf-8") as fh:
        fh.write(",".join(cols) + ",dbic\n")
        for r in sorted(ok, key=lambda r: r["pulsar"]):
            fh.write(",".join("" if r.get(c) is None else (("%.4f" % r[c]) if isinstance(r.get(c), float) else str(r[c]))
                              for c in cols) + ',"%s"\n' % json.dumps(r["dbic"]).replace('"', "'"))
    nul = [r for r in ok if not r["control"]]
    ctl = [r for r in ok if r["control"]]
    L = []
    p = L.append
    p("# Frakcja nullingu z modelu mieszaniny energii a nulle automatu")
    p("")
    p("Model jak w TPA XIX (Keith+2026, wzory 2-3), dopasowany metodą największej wiarygodności (nie nested sampling), "
      "wybór rozkładu przez BIC, niepewność z bootstrapu po impulsach (%d losowań). Impulsy wycięte przez automat jako RFI "
      "pominięte. Energie: okno E_off rozłączne z oknem linii bazowej." % NBOOT)
    p("")
    p("- pulsary z nullami: %d, kontrolne bez nulli: %d, błędy: %d" % (len(nul), len(ctl), len(results) - len(ok)))
    for r in results:
        if "error" in r:
            p("  - %s: %s" % (r["pulsar"], r["error"]))
    if ctl:
        p("- **kontrola** (pulsary, w których automat nie wyciął nulli): NF z modelu mediana %.3f, maks. %.3f; "
          "NF 3σ mediana %.3f" % (np.median([r["NF_model"] for r in ctl]), max(r["NF_model"] for r in ctl),
                                  np.median([r["NF_3sigma"] for r in ctl])))
    if nul:
        d = np.array([r["NF_auto"] - r["NF_model"] for r in nul])
        d3 = np.array([r["NF_3sigma"] - r["NF_model"] for r in nul])
        p("- **automat − model:** mediana %+.3f, automat niżej o >0.05 w %d, wyżej o >0.05 w %d pulsarach" %
          (np.median(d), int((d < -0.05).sum()), int((d > 0.05).sum())))
        p("- **3σ − model:** mediana %+.3f (oczekiwane > 0: próg liczy słabą emisję jako nulle)" % np.median(d3))
        tot_auto = sum(r["n_null_auto"] for r in nul)
        p("- **nulle automatu wg modelu:** P_null ≥ 0.5: %d z %d (%.0f%%); P_null < 0.05 (emisja wycięta jako null): %d (%.1f%%)" %
          (sum(r["auto_null_with_p_ge_05"] for r in nul), tot_auto,
           100 * sum(r["auto_null_with_p_ge_05"] for r in nul) / max(tot_auto, 1),
           sum(r["auto_null_with_p_lt_005"] for r in nul), 100 * sum(r["auto_null_with_p_lt_005"] for r in nul) / max(tot_auto, 1)))
        p("- **nulle wg modelu niewycięte przez automat:** %d (automat zostawia krótkie nulle, gdy ich wycięcie nie podnosi S/N)" %
          sum(r["model_null_not_cut"] for r in nul))
        p("- impulsy niejednoznaczne (0.05 < P_null < 0.95): mediana %.0f%% na pulsar" %
          (100 * np.median([r["ambiguous_frac"] for r in nul])))
        p("")
        p("## Pulsary z literaturą")
        p("")
        p("| pulsar | NF automat | NF 3σ | NF model [16–84%] | literatura | model |")
        p("|---|---|---|---|---|---|")
        for r in sorted([r for r in nul if r["NF_lit"] is not None], key=lambda r: r["pulsar"]):
            p("| %s | %.1f%% | %.1f%% | %.1f%% [%.1f–%.1f] | %.0f%% | %s |" % (
                r["pulsar"], 100 * r["NF_auto"], 100 * r["NF_3sigma"], 100 * r["NF_model"],
                100 * r.get("NF_model_lo", float("nan")), 100 * r.get("NF_model_hi", float("nan")), r["NF_lit"], r["model"]))
        lit = [r for r in nul if r["NF_lit"] is not None]
        if lit:
            p("")
            p("Mediana różnicy z literaturą: model %+.1f pp, automat %+.1f pp, 3σ %+.1f pp." % (
                np.median([100 * r["NF_model"] - r["NF_lit"] for r in lit]),
                np.median([100 * r["NF_auto"] - r["NF_lit"] for r in lit]),
                np.median([100 * r["NF_3sigma"] - r["NF_lit"] for r in lit])))
        p("")
        p("## Najwięcej emisji wyciętej jako null (P_null < 0.05)")
        p("")
        p("| pulsar | nulle automatu | z P_null < 0.05 | średnia E nulli automatu | NF automat / model |")
        p("|---|---|---|---|---|")
        for r in sorted(nul, key=lambda r: -r["auto_null_with_p_lt_005"])[:20]:
            if r["auto_null_with_p_lt_005"] == 0:
                break
            p("| %s | %d | %d | %.2f | %.1f%% / %.1f%% |" % (r["pulsar"], r["n_null_auto"], r["auto_null_with_p_lt_005"],
                                                       r["mean_E_auto_nulls"], 100 * r["NF_auto"], 100 * r["NF_model"]))
    with open(os.path.join(OUT, "RAPORT_nulle.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


# ------------------------------------------------------------------ test syntetyczny
def synth(nf_true, n=1500, snr=30.0, seed=0):
    """impulsy: null (E=0) z prawd. nf, emisja lognormalna; szum o sigma = 1/snr (w <E_emisji>)."""
    rng = np.random.default_rng(seed)
    isnull = rng.random(n) < nf_true
    em = rng.lognormal(mean=-0.1, sigma=0.5, size=n)
    em[isnull] = 0.0
    s = 1.0 / snr
    eon = em + rng.normal(0, s, n)
    eoff = rng.normal(0, s, n)
    return eon, eoff, isnull


def selftest():
    print("TEST: dane syntetyczne, znane NF")
    for nf_true, snr in ((0.0, 30), (0.2, 30), (0.5, 30), (0.5, 5), (0.3, 3), (0.0, 3)):
        eon, eoff, isnull = synth(nf_true, snr=snr)
        norm = eon.mean()
        e_on, e_off = eon / norm, eoff / norm
        fits = {"+".join(m): fit(m, e_on, e_off) for m in MODELS}
        n = len(e_on)
        bic = {m: 2 * r.fun + len(r.x) * math.log(2 * n) for m, r in fits.items()}
        b = min(bic, key=bic.get)
        model = tuple(b.split("+"))
        th = fits[b].x
        nf = float(special.expit(th[2]))
        f, fn = f_on(model, th, e_on)
        pn = fn / f
        sE = np.std(e_off)
        nf3 = float((e_on < 3 * sE).mean())
        acc = float(((pn >= 0.5) == isnull).mean())
        print("  NF prawdziwe %.2f, S/N impulsu %4.0f: model %.3f (N+%s), 3σ %.3f, trafność P_null %.1f%%, niejednoznaczne %.0f%%" %
              (nf_true, snr, nf, b, nf3, 100 * acc, 100 * float(((pn > 0.05) & (pn < 0.95)).mean())))


def main():
    if "--test" in sys.argv:
        selftest()
        return
    files = sorted(glob.glob(os.path.join(OUT, "energies", "*.csv")))
    want = [a for a in sys.argv[1:] if not a.startswith("-")]
    if want:
        files = [f for f in files if os.path.basename(f)[:-4] in want]
    if not files:
        print("Brak plików energii w %s - najpierw: SPATS_NO_MAIN=1 julia --project=. tools/nulling_energies.jl" % OUT)
        return
    cf = os.path.join(OUT, "control.txt")
    control = set(open(cf).read().split()) if os.path.isfile(cf) else set()
    print("MODEL NULLI: %d pulsarów, %d procesów" % (len(files), min(32, os.cpu_count() or 1)))
    with Pool(min(32, os.cpu_count() or 1)) as pool:
        results = pool.map(work, files)
    report(results, control)
    print("\nWyniki: %s/summary.csv, RAPORT_nulle.md, plots/<psr>.pdf" % OUT)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
