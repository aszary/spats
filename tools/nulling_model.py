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
MODEL_VERSION = 4   # zmiana modelu -> stare wyniki z fits/ nie sa uzywane
# v4: sigma energii empiryczne (std E_off) w f_i3 i NF_3sigma (s0/sE ~2 na danych TPA - rozrzut z samego
#     szumu odbiornika zanizal sigma), bootstrap blokowy (nulle seriami), test fazy (filtr dopasowany,
#     CHIME Ng+2020), dolna granica NF (CHIME/Esamdin), tryb NULLING_RELABEL=1 (nowe etykiety automatu)
BLOCK = int(os.environ.get("NULLING_BLOCK", "20"))   # dlugosc bloku bootstrapu (impulsy); 1 = zwykly
FI3_MIN = 0.3       # TPA XIX: NF wiarygodne tylko, gdy > 3 sigma jest >= 30% impulsow (f_i3)
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
    dphi = a[:, 4] if a.shape[1] > 4 else np.full(len(a), np.nan)   # przesuniecie fazy z filtra dopasowanego
    return meta, a[:, 0].astype(int), a[:, 1], a[:, 2], a[:, 3].astype(int), dphi


def auto_compare(res, e_on, lab_k, pnull):
    """porownanie impuls po impulsie: nulle automatu (etykieta 2) a P_null z modelu."""
    an = lab_k == 2
    mn = pnull >= 0.5
    res.update(n_null_auto=int(an.sum()), NF_auto=float(an.sum()) / len(an),
               n_null_model=int(mn.sum()), ambiguous_frac=float(((pnull > 0.05) & (pnull < 0.95)).mean()),
               auto_null_with_p_lt_005=int((an & (pnull < 0.05)).sum()),
               auto_null_with_p_ge_05=int((an & mn).sum()),
               model_null_cut_by_auto=int((mn & an).sum()), model_null_not_cut=int((mn & ~an).sum()),
               mean_E_auto_nulls=float(e_on[an].mean()) if an.any() else float("nan"),
               mean_E_model_nulls=float(e_on[mn].mean()) if mn.any() else float("nan"))


def relabel(path, res):
    """NULLING_RELABEL=1: nowe etykiety automatu (te same impulsy RFI) bez ponownego dopasowania.
    Zwraca None, gdy zbior impulsow sie zmienil (inne RFI) - wtedy trzeba dopasowac od nowa."""
    meta, idx, Eon, Eoff, lab, dphi = load(path)
    keep = lab != 1
    f = os.path.join(OUT, "fits", res["pulsar"] + ".npz")
    if not os.path.isfile(f):
        return None
    z = dict(np.load(f))
    if len(z["idx"]) != keep.sum() or not np.array_equal(z["idx"], idx[keep]):
        return None
    z["lab"] = lab[keep]
    np.savez(f, **z)
    auto_compare(res, z["e_on"], z["lab"], z["pnull"])
    return res


def analyse(path, nboot=NBOOT, seed=1):
    meta, idx, Eon, Eoff, lab, dphi = load(path)
    psr = meta["pulsar"]
    keep = lab != 1                                   # bez RFI wycietego przez automat
    Eon_k, Eoff_k, lab_k, idx_k, dphi_k = Eon[keep], Eoff[keep], lab[keep], idx[keep], dphi[keep]
    norm = float(np.mean(Eon_k))                      # <E_on> razem z nullami
    if not norm > 0:
        return {"pulsar": psr, "error": "<E_on> <= 0"}
    e_on, e_off = Eon_k / norm, Eoff_k / norm
    n = len(e_on)
    res = {"pulsar": psr, "N": int(len(Eon)), "n_fit": n, "n_rfi_auto": int((lab == 1).sum()),
           "n_null_auto": int((lab == 2).sum()), "control": False}
    # NF z progu 3 sigma (Bhattacharyya+2010). sigma energii: EMPIRYCZNE std(E_off) (v4); rozrzut z samego
    # szumu odbiornika sqrt(n_on)*sigma_bin zanizal je ~2x (linia bazowa / skorelowany szum)
    sE_rad = math.sqrt(float(meta["n_on"])) * float(meta["sigma_bin"]) / norm
    sE = float(np.std(e_off))
    nn3 = int((e_on < 3 * sE).sum())
    res.update(NF_3sigma=nn3 / n, NF_3sigma_err=math.sqrt(nn3) / n, NF_auto=res["n_null_auto"] / n,
               s0_over_sE=float(np.std(e_off)) / sE_rad, baseline=meta.get("baseline", "srednia"),
               NF_low=float(meta["NF_low"]) if "NF_low" in meta else float("nan"))
    # f_i3: odsetek impulsow > 3 sigma szumu. Ponizej FI3_MIN emisji nie da sie odroznic od
    # szumu i "null" znaczy tylko "niewykryty" (J1755-2025, J1652-1400: NF 67-81% przy
    # ogonie kilku jasnych impulsow) - NF raportowane, ale oznaczone jako niewiarygodne
    res["f_i3"] = float((e_on > 3 * sE).mean())
    res["reliable"] = res["f_i3"] >= FI3_MIN
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
    nblk = int(math.ceil(n / BLOCK))
    for _ in range(nboot if n <= 5000 else min(nboot, 8)):
        if BLOCK <= 1:
            j = rng.integers(0, n, n)
        else:                                   # bootstrap blokowy: nulle wystepuja seriami
            st = rng.integers(0, max(n - BLOCK, 1), nblk)
            j = np.concatenate([np.arange(a, min(a + BLOCK, n)) for a in st])[:n]
        rb = fit(b["model"], e_on[j], e_off[j], x0s=[th])
        bs.append(float(special.expit(rb.x[2])))
    if bs:
        res.update(NF_model_lo=float(np.percentile(bs, 16)), NF_model_hi=float(np.percentile(bs, 84)))
    auto_compare(res, e_on, lab_k, pnull)
    res["boot_block"] = BLOCK
    # test fazy (CHIME, Ng+2020): filtr dopasowany - w nullu faza dopasowania jest losowa, w emisji ~0.
    # Odsetek impulsow z |dphi| <= polowa okna on-pulse; losowo ~ n_on/nbin
    if np.isfinite(dphi_k).any() and "nbin" in meta:
        hw, nbin = float(meta["n_on"]) / 2, float(meta["nbin"])
        inw = np.abs(dphi_k) <= hw
        g = {"P95": pnull >= 0.95, "P50": pnull >= 0.5, "auto": lab_k == 2, "emisja": pnull < 0.05}
        res["phase_random"] = 2 * hw / nbin
        for k, m in g.items():
            res["phase_in_" + k] = float(inw[m].mean()) if m.sum() >= 5 else float("nan")
    res["NF_lit"] = LIT_NF.get(psr)
    res["n_on"], res["nbin"] = int(meta["n_on"]), int(meta.get("nbin", 0))
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
        res = json.load(open(cache, encoding="utf-8"))
        if os.environ.get("NULLING_RELABEL") != "1":
            return res
        r2 = relabel(path, res)
        if r2 is not None:
            with open(cache, "w", encoding="utf-8") as fh:
                json.dump(r2, fh)
            return r2
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


def degenerate(r):
    """Model nie odroznia nulli od slabej emisji: NF_model ~0, choc automat wycial >=5 nulli
    o sredniej energii zgodnej z zerem (< 3 bledy sredniej szumu; J1947+0915: 86 nulli, E -0.03,
    NF model 0, NF 3 sigma 0.57). Szeroki rozklad emisji obejmuje wtedy zero i 'zjada' skladnik null."""
    n = r.get("n_null_auto", 0)
    # model musi byc wyraznie ponizej automatu: przy NF ~2% model i automat zgodne (J1820-0427: 1.9% / 1.9%)
    if r.get("control") or n < 5 or r["NF_model"] >= 0.02 or r["NF_model"] >= 0.5 * r["NF_auto"]:
        return False
    m = r.get("mean_E_auto_nulls", float("nan"))
    return bool(np.isfinite(m) and m < 3 * r["s0"] / math.sqrt(n))


def report(results, control):
    ok = [r for r in results if "error" not in r]
    for r in ok:
        r["control"] = r["pulsar"] in control
        r["degenerate"] = degenerate(r)          # liczone z wynikow (cache), bez ponownego dopasowania
        if r["degenerate"]:
            r["reliable"] = False
    cols = ["pulsar", "control", "N", "n_fit", "n_rfi_auto", "n_null_auto", "n_null_model", "NF_auto", "NF_3sigma",
            "NF_3sigma_err", "NF_model", "NF_model_lo", "NF_model_hi", "NF_lit", "model", "f_i3", "reliable", "degenerate", "ambiguous_frac",
            "auto_null_with_p_lt_005", "auto_null_with_p_ge_05", "model_null_not_cut", "mean_E_auto_nulls",
            "mean_E_model_nulls", "s0_over_sE", "baseline", "NF_low", "phase_random", "phase_in_P95", "phase_in_P50",
            "phase_in_auto", "phase_in_emisja"]
    with open(os.path.join(OUT, "summary.csv"), "w", encoding="utf-8") as fh:
        fh.write(",".join(cols) + ",dbic\n")
        for r in sorted(ok, key=lambda r: r["pulsar"]):
            fh.write(",".join("" if r.get(c) is None else (("%.4f" % r[c]) if isinstance(r.get(c), float) else str(r[c]))
                              for c in cols) + ',"%s"\n' % json.dumps(r["dbic"]).replace('"', "'"))
    allnul = [r for r in ok if not r["control"]]
    unrel = [r for r in allnul if not r.get("reliable", True) and not r["degenerate"]]
    degen = [r for r in allnul if r["degenerate"]]
    nul = [r for r in allnul if r.get("reliable", True)]
    ctl = [r for r in ok if r["control"]]
    L = []
    p = L.append
    p("# Frakcja nullingu z modelu mieszaniny energii a nulle automatu")
    p("")
    p("Model jak w TPA XIX (Keith+2026, wzory 2-3), dopasowany metodą największej wiarygodności (nie nested sampling), "
      "wybór rozkładu przez BIC, niepewność z bootstrapu po impulsach (%d losowań). Impulsy wycięte przez automat jako RFI "
      "pominięte. Energie: okno E_off rozłączne z oknem linii bazowej." % NBOOT)
    p("")
    p("- pulsary z nullami: %d, w tym %d z wiarygodnym NF (f_i3 ≥ %.1f: co najmniej %.0f%% impulsów > 3σ, kryterium TPA XIX); "
      "kontrolne bez nulli: %d; błędy: %d" % (len(allnul), len(nul), FI3_MIN, 100 * FI3_MIN, len(ctl), len(results) - len(ok)))
    p("- **dalsze liczby tylko dla pulsarów z wiarygodnym NF.** Niewiarygodne (za słabe, emisji nie da się odróżnić od szumu): %s" %
      ", ".join("%s (f_i3 %.2f)" % (r["pulsar"], r["f_i3"]) for r in sorted(unrel, key=lambda r: r["f_i3"])))
    p("- **model nie odróżnia nulli od słabej emisji** (NF model < 0.02, choć nulle automatu mają energię zgodną z zerem; "
      "NF jest tu co najmniej NF automatu): %d pulsarów: %s" % (len(degen), ", ".join(
          "%s (automat %.1f%%, 3σ %.1f%%)" % (r["pulsar"], 100 * r["NF_auto"], 100 * r["NF_3sigma"])
          for r in sorted(degen, key=lambda r: -r["NF_auto"]))))
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
        p("- **sigma energii:** empiryczne std(E_off); rozrzut E_off / rozrzut z samego szumu odbiornika: mediana %.2f "
          "(1 = stabilna linia bazowa)" % np.median([r["s0_over_sE"] for r in ok]))
        bl = {}
        for r in ok:
            bl[r.get("baseline", "srednia")] = bl.get(r.get("baseline", "srednia"), 0) + 1
        p("- **linia bazowa** (wybrana w każdym pulsarze: najmniejszy rozrzut E_off): " +
          ", ".join("%s %d" % kv for kv in sorted(bl.items())))
        p("- **niepewność NF:** bootstrap blokowy, bloki po %d impulsów (nulle występują seriami)" % BLOCK)
        lo = [(r["NF_low"], r["NF_model"]) for r in nul if np.isfinite(r.get("NF_low", float("nan")))]
        if lo:
            d = np.array([100 * (a - b) for a, b in lo])
            p("- **NF metodą CHIME (Ng+2020):** najwięcej najsłabszych impulsów bez wykrywalnej emisji w sumie (< 3σ); "
              "NF CHIME − NF model: mediana %+.1f pp; różnica > 5 pp: %d. Uwaga: na danych syntetycznych ta metoda "
              "ZAWYŻA NF (0.35 przy prawdziwym 0.30), bo najsłabsze impulsy z emisją też nie dają wykrywalnej sumy" %
              (np.median(d), int((np.abs(d) > 5).sum())))
        ph = [r for r in nul if np.isfinite(r.get("phase_in_P95", float("nan")))]
        if ph:
            p("- **test fazy (filtr dopasowany, CHIME):** odsetek impulsów z fazą dopasowania w oknie on-pulse; "
              "losowo mediana %.2f. Nulle P≥0.95: %.2f, P≥0.5: %.2f, nulle automatu: %.2f, emisja: %.2f "
              "(prawdziwe nulle ≈ losowo, emisja ≈ 1)" % (
                  np.median([r["phase_random"] for r in ph]), np.nanmedian([r["phase_in_P95"] for r in ph]),
                  np.nanmedian([r.get("phase_in_P50", np.nan) for r in ph]),
                  np.nanmedian([r.get("phase_in_auto", np.nan) for r in ph]),
                  np.nanmedian([r.get("phase_in_emisja", np.nan) for r in ph])))
        p("")
        p("## Pulsary z literaturą")
        p("")
        p("| pulsar | NF automat | NF 3σ | NF model [16–84%] | literatura | model | f_i3 |")
        p("|---|---|---|---|---|---|---|")
        for r in sorted([r for r in allnul if r["NF_lit"] is not None], key=lambda r: r["pulsar"]):
            p("| %s | %.1f%% | %.1f%% | %.1f%% [%.1f–%.1f]%s | %.0f%% | %s | %.2f |" % (
                r["pulsar"], 100 * r["NF_auto"], 100 * r["NF_3sigma"], 100 * r["NF_model"],
                100 * r.get("NF_model_lo", float("nan")), 100 * r.get("NF_model_hi", float("nan")),
                "" if r.get("reliable", True) else " (model nie rozróżnia)" if r["degenerate"] else " (niewiarygodne)", r["NF_lit"], r["model"], r["f_i3"]))
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
    export_nulls(ok)


def to_ranges(ix):
    """[3,4,5,9] -> '3-5 9' (jak zaps2, numeracja impulsow od 0)."""
    ix = sorted(int(i) for i in ix)
    out, k = [], 0
    while k < len(ix):
        j = k
        while j + 1 < len(ix) and ix[j + 1] == ix[j] + 1:
            j += 1
        out.append(str(ix[k]) if j == k else "%d-%d" % (ix[k], ix[j]))
        k = j + 1
    return " ".join(out)


def export_nulls(ok):
    """Lista nulli wg modelu (P_null >= 0.5) dla kazdego pulsara - do badan nullingu (dlugosci
    nulli, granice, R1-R3). nulls_model/<psr>.csv: impuls po impulsie (E_on, P_null, null, niejednoznaczny,
    etykieta automatu); nulls_model.csv: zbiorczo (NF, liczba i srednia dlugosc odcinkow nulli,
    nulle w zakresach jak zaps2). Impulsy wyciete jako RFI nie sa w listach (brak energii)."""
    d = os.path.join(OUT, "nulls_model")
    os.makedirs(d, exist_ok=True)
    rows = []
    for r in sorted(ok, key=lambda r: r["pulsar"]):
        f = os.path.join(OUT, "fits", r["pulsar"] + ".npz")
        if not os.path.isfile(f):
            continue
        z = np.load(f)
        idx, pn, e_on, lab = z["idx"], z["pnull"], z["e_on"], z["lab"]
        isn = pn >= 0.5
        with open(os.path.join(d, r["pulsar"] + ".csv"), "w", encoding="utf-8") as fh:
            fh.write("i,E_on,P_null,null,ambiguous,label_auto\n")
            for i, e, p, l in zip(idx, e_on, pn, lab):
                fh.write("%d,%.5f,%.4f,%d,%d,%d\n" % (i, e, p, p >= 0.5, 0.05 < p < 0.95, l))
        # odcinki nulli (kolejne impulsy; impuls wyciety jako RFI przerywa odcinek)
        runs, cur, prev = [], 0, None
        for i, s in zip(idx, isn):
            if s and prev is not None and i == prev + 1 and cur > 0:
                cur += 1
            elif s:
                if cur:
                    runs.append(cur)
                cur = 1
            else:
                if cur:
                    runs.append(cur)
                cur = 0
            prev = i
        if cur:
            runs.append(cur)
        # status: czy P_null z modelu nadaje sie do wycinania (kontrolne tez, np. J2111+2106 - automat bez nulli)
        status = ("model_nie_rozroznia" if r.get("degenerate") else
                  "za_slaby" if r.get("f_i3", 1.0) < FI3_MIN else "ok")
        rows.append([r["pulsar"], status, len(idx), int(isn.sum()), "%.4f" % r["NF_model"],
                     len(runs), "%.2f" % (np.mean(runs) if runs else 0.0), max(runs) if runs else 0,
                     int(((pn > 0.05) & (pn < 0.95)).sum()), '"%s"' % to_ranges(idx[isn])])
    with open(os.path.join(OUT, "nulls_model.csv"), "w", encoding="utf-8") as fh:
        fh.write("pulsar,status,n_pulses,n_null,NF_model,n_null_runs,mean_null_len,max_null_len,n_ambiguous,nulls\n")
        for row in rows:
            fh.write(",".join(str(x) for x in row) + "\n")
    nok = [x for x in rows if x[1] == "ok"]
    print("\nLISTY NULLI WG MODELU: %d pulsarów (%d z wiarygodnym NF) -> %s/nulls_model.csv, nulls_model/<psr>.csv"
          % (len(rows), len(nok), OUT))


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
