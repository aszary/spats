#!/usr/bin/env python3
"""Szybkie sprawdzenia modelu nulli (z gotowych wynikow nulling_model.py, bez nowych energii):

1. Linia bazowa: czy potrzebna korekta GP (jak w TPA XIX). s0/sE = rozrzut energii poza impulsem
   (E_off) wobec rozrzutu z samego szumu radiometru (sqrt(n_on) * sigma_bin). Gdy linia bazowa
   "plywa" w obrebie impulsu, szum w binach jest skorelowany i s0/sE > 1. ~1 = GP nic nie da.
2. Bootstrap blokowy: nulle wystepuja seriami, a zwykly bootstrap zaklada niezalezne impulsy, wiec
   moze zanizac niepewnosc NF. Losujemy bloki po BLOCK kolejnych impulsow i porownujemy szerokosc
   przedzialu 16-84% NF z bootstrapem zwyklym (z nulling_model.py).

    NULLING_OUT=/home/psr/data/OUTPUT/nulling_model_v4 python3 tools/nulling_checks.py
Wyniki: $NULLING_OUT/checks/RAPORT_checks.md, block_bootstrap.csv
"""
import csv
import importlib.util
import json
import math
import os
from multiprocessing import Pool

import numpy as np
from scipy import special

_spec = importlib.util.spec_from_file_location("nm", os.path.join(os.path.dirname(os.path.abspath(__file__)), "nulling_model.py"))
nm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nm)
OUT = nm.OUT
BLOCK = int(os.environ.get("NC_BLOCK", "20"))   # porownanie z przedzialem z summary.csv (v4: juz blokowy 20)
NB = 30


def block_boot(psr):
    try:
        z = np.load(os.path.join(OUT, "fits", psr + ".npz"))
        r = json.load(open(os.path.join(OUT, "fits", "%s.v%d.json" % (psr, nm.MODEL_VERSION)), encoding="utf-8"))
        e_on, e_off, th, model = z["e_on"], z["e_off"], z["th"], tuple(str(m) for m in z["model"])
        n = len(e_on)
        nb = NB if n <= 5000 else 8
        rng = np.random.default_rng(7)
        nblk = int(math.ceil(n / BLOCK))
        bs = []
        for _ in range(nb):
            starts = rng.integers(0, max(n - BLOCK, 1), nblk)
            j = np.concatenate([np.arange(s, min(s + BLOCK, n)) for s in starts])[:n]
            rb = nm.fit(model, e_on[j], e_off[j], x0s=[th])
            bs.append(float(special.expit(rb.x[2])))
        lo, hi = np.percentile(bs, 16), np.percentile(bs, 84)
        return {"pulsar": psr, "NF": r["NF_model"], "iid_lo": r.get("NF_model_lo"), "iid_hi": r.get("NF_model_hi"),
                "blk_lo": float(lo), "blk_hi": float(hi)}
    except Exception as e:
        return {"pulsar": psr, "error": repr(e)}


def main():
    rows = list(csv.DictReader(open(os.path.join(OUT, "summary.csv"), encoding="utf-8")))
    os.makedirs(os.path.join(OUT, "checks"), exist_ok=True)
    L = []
    p = L.append
    f = lambda r, k: float(r[k]) if r.get(k) not in (None, "") else float("nan")
    # 1. linia bazowa
    ratio = [(r["pulsar"], f(r, "s0_over_sE")) for r in rows if np.isfinite(f(r, "s0_over_sE"))]
    v = np.array([x for _, x in ratio])
    p("# Szybkie sprawdzenia modelu nulli")
    p("")
    p("## 1. Linia bazowa: czy potrzebna korekta GP")
    p("")
    p("s0/sE = rozrzut energii poza impulsem / rozrzut z samego szumu. ~1: linia bazowa stabilna, GP nic nie da.")
    p("")
    p("- pulsarów: %d;  mediana s0/sE %.3f;  16–84%%: %.3f–%.3f" % (len(v), np.median(v), *np.percentile(v, [16, 84])))
    p("- s0/sE > 1.1: %d,  > 1.2: %d,  > 1.5: %d" % ((v > 1.1).sum(), (v > 1.2).sum(), (v > 1.5).sum()))
    big = sorted([(x, ps) for ps, x in ratio if x > 1.2], reverse=True)
    if big:
        nf = {r["pulsar"]: r for r in rows}
        p("- największe (tu GP mogłoby pomóc): " + ", ".join("%s %.2f (NF %.1f%%)" % (ps, x, 100 * f(nf[ps], "NF_model"))
                                                          for x, ps in big[:20]))
    # 2. bootstrap blokowy - pulsary z wiarygodnym NF
    rel = [r["pulsar"] for r in rows if r.get("reliable") == "True" and r.get("control") == "False"
           and r.get("degenerate", "False") != "True"]
    print("BOOTSTRAP BLOKOWY (bloki po %d impulsów): %d pulsarów, %d procesów" % (BLOCK, len(rel), min(32, os.cpu_count() or 1)),
          flush=True)
    with Pool(min(32, os.cpu_count() or 1)) as pool:
        res = pool.map(block_boot, rel)
    ok = [r for r in res if "error" not in r and r["iid_lo"] is not None]
    w_iid = np.array([r["iid_hi"] - r["iid_lo"] for r in ok])
    w_blk = np.array([r["blk_hi"] - r["blk_lo"] for r in ok])
    q = w_blk / np.maximum(w_iid, 1e-4)
    with open(os.path.join(OUT, "checks", "block_bootstrap.csv"), "w", encoding="utf-8") as fh:
        fh.write("pulsar,NF,iid_lo,iid_hi,blk_lo,blk_hi,width_ratio\n")
        for r, x in zip(ok, q):
            fh.write("%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f\n" % (r["pulsar"], r["NF"], r["iid_lo"], r["iid_hi"], r["blk_lo"], r["blk_hi"], x))
    p("")
    p("## 2. Bootstrap blokowy (bloki po %d kolejnych impulsów) a zwykły" % BLOCK)
    p("")
    p("- pulsarów: %d (błędy: %d)" % (len(ok), len(res) - len(ok)))
    p("- szerokość przedziału 16–84%% NF: zwykły mediana %.1f pp, blokowy mediana %.1f pp" % (100 * np.median(w_iid), 100 * np.median(w_blk)))
    p("- blokowy / zwykły: mediana %.2f;  > 1.5×: %d,  > 2×: %d pulsarów" % (np.median(q), (q > 1.5).sum(), (q > 2).sum()))
    top = sorted(zip(q, ok), key=lambda t: -t[0])[:12]
    p("- najwięcej: " + ", ".join("%s %.1f× (NF %.1f%% [%.1f–%.1f] → [%.1f–%.1f])" % (
        r["pulsar"], x, 100 * r["NF"], 100 * r["iid_lo"], 100 * r["iid_hi"], 100 * r["blk_lo"], 100 * r["blk_hi"]) for x, r in top))
    with open(os.path.join(OUT, "checks", "RAPORT_checks.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))
    print("\nWyniki: %s/checks/RAPORT_checks.md, block_bootstrap.csv" % OUT)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
