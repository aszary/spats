#!/usr/bin/env python3
"""Porownanie naszej frakcji nullingu z wynikami TPA XIX (Keith+2026, te same dane MeerKAT TPA).
Wyniki TPA XIX: Zenodo 10.5281/zenodo.18982781 (tpa_pulse_energy_dist.zip, 9.6 GB) - z archiwum wyciagniete
tylko najlepsze modele (najwyzsze ln Z) do input/tpa19_nf.csv (NF, NF_e, Null5/Null95 = 5/95%, int3s = f_i3).
    NULLING_OUT=/home/psr/data/OUTPUT/final_v6/model/ python3 tools/compare_tpa19.py
Wynik: $NULLING_OUT/RAPORT_tpa19.md, tpa19_compare.csv"""
import csv
import math
import os
import sys

import numpy as np

OUT = os.environ.get("NULLING_OUT", "/home/psr/data/OUTPUT/final_v7/model/")
HERE = os.path.dirname(os.path.abspath(__file__))
TPA = os.path.join(HERE, "..", "input", "tpa19_nf.csv")
FI3 = 0.3
AB = {"LN": "LN", "Nm": "G", "PL": "PL", "&delta;": "delta"}


def fnum(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def main():
    tpa = {r["psr"]: r for r in csv.DictReader(open(TPA, encoding="utf-8")) if not r["model"].startswith("BLAD")}
    ours = {r["pulsar"]: r for r in csv.DictReader(open(os.path.join(OUT, "summary.csv"), encoding="utf-8"))}
    dyf = os.path.join(OUT, "dynesty", "summary.csv")
    dyn = {r["pulsar"]: r for r in csv.DictReader(open(dyf, encoding="utf-8"))} if os.path.isfile(dyf) else {}
    common = sorted(set(tpa) & set(ours))
    L = []
    p = L.append
    p("# Frakcja nullingu: nasz model a TPA XIX (Keith+2026)\n")
    p("Te same obserwacje MeerKAT TPA; TPA XIX: dynesty, najlepszy model wg ln Z, f_i3 z rozkładu emisji. "
      "My: największa wiarygodność + BIC, bootstrap blokowy, RFI wycięte automatem.\n")
    p("- pulsarów w TPA XIX: %d, u nas: %d, wspólnych: %d" % (len(tpa), len(ours), len(common)))
    rows = []
    for psr in common:
        t, o = tpa[psr], ours[psr]
        nf_t, e_t, fi_t = fnum(t["NF"]), fnum(t["NF_e"]), fnum(t["int3s"])
        nf_o, lo, hi, fi_o = fnum(o["NF_model"]), fnum(o.get("NF_model_lo")), fnum(o.get("NF_model_hi")), fnum(o["f_i3"])
        e_o = (hi - lo) / 2 if math.isfinite(hi - lo) else float("nan")
        d = dyn.get(psr, {})
        rows.append(dict(psr=psr, nf_t=nf_t, e_t=e_t, n5=fnum(t["Null5"]), n95=fnum(t["Null95"]), fi_t=fi_t,
                         model_t="+".join(AB.get(x, x) for x in t["model"].split("+")),
                         nf_o=nf_o, e_o=e_o, fi_o=fi_o, f3s_o=fnum(o.get("f_3sigma")), model_o=o["model"],
                         auto=fnum(o["NF_auto"]), s3=fnum(o["NF_3sigma"]), degen=o.get("degenerate") == "True",
                         rel_o=o.get("reliable") == "True", nf_d=fnum(d.get("NF_med")), lit=fnum(o.get("NF_lit"))))
    R = rows
    # f_i3 i wiarygodnosc
    fi = np.array([[r["fi_t"], r["fi_o"]] for r in R if math.isfinite(r["fi_t"]) and math.isfinite(r["fi_o"])])
    p("\n## f_i3 (czy pulsar jest dość jasny na wiarygodne NF)\n")
    if len(fi):
        rel_t, rel_o = fi[:, 0] > FI3, fi[:, 1] >= FI3
        p("- f_i3 nasze − TPA XIX: mediana %+.3f, 16–84%%: %+.3f…%+.3f (z %d)" % (
            np.median(fi[:, 1] - fi[:, 0]), *np.percentile(fi[:, 1] - fi[:, 0], [16, 84]), len(fi)))
        p("- wiarygodne (f_i3 > 0.3): TPA XIX %d, my %d; zgodna klasyfikacja: %d z %d (%.0f%%); tylko TPA: %d, tylko my: %d" % (
            rel_t.sum(), rel_o.sum(), (rel_t == rel_o).sum(), len(fi), 100 * (rel_t == rel_o).mean(),
            (rel_t & ~rel_o).sum(), (~rel_t & rel_o).sum()))
    # NF dla wiarygodnych w obu
    G = [r for r in R if r["fi_t"] > FI3 and r["fi_o"] >= FI3 and not r["degen"]
         and math.isfinite(r["nf_t"]) and math.isfinite(r["nf_o"])]
    p("\n## Frakcja nullingu (pulsary wiarygodne w obu: %d)\n" % len(G))
    if G:
        d = np.array([r["nf_o"] - r["nf_t"] for r in G])
        err = np.array([math.hypot(r["e_t"] if math.isfinite(r["e_t"]) else 0, r["e_o"] if math.isfinite(r["e_o"]) else 0) for r in G])
        p("- NF nasze − TPA XIX: mediana %+.1f pp, 16–84%%: %+.1f…%+.1f pp, |różnica| mediana %.1f pp" % (
            100 * np.median(d), *(100 * np.percentile(d, [16, 84])), 100 * np.median(np.abs(d))))
        p("- zgodne w granicach 2σ (łączna niepewność): %d z %d (%.0f%%); |różnica| < 2 pp: %d, < 5 pp: %d" % (
            (np.abs(d) <= 2 * np.maximum(err, 1e-3)).sum(), len(G), 100 * (np.abs(d) <= 2 * np.maximum(err, 1e-3)).mean(),
            (np.abs(d) < 0.02).sum(), (np.abs(d) < 0.05).sum()))
        det_t = np.array([r["nf_t"] > 3 * r["e_t"] and r["nf_t"] > 0.01 for r in G])
        det_o = np.array([r["nf_o"] > 3 * max(r["e_o"], 1e-4) and r["nf_o"] > 0.01 for r in G])
        p("- nulling wykryty (NF > 3σ i > 1%%): TPA XIX %d, my %d, zgodnie %d z %d (%.0f%%)" % (
            det_t.sum(), det_o.sum(), (det_t == det_o).sum(), len(G), 100 * (det_t == det_o).mean()))
        for lab, sel in (("NF TPA < 5%", [r["nf_t"] < 0.05 for r in G]), ("5–30%", [0.05 <= r["nf_t"] < 0.3 for r in G]),
                         ("≥ 30%", [r["nf_t"] >= 0.3 for r in G])):
            sel = np.array(sel)
            if sel.any():
                p("  - %s: %d pulsarów, różnica mediana %+.1f pp" % (lab, sel.sum(), 100 * np.median(d[sel])))
        same = sum(r["model_t"] == r["model_o"][2:] for r in G)
        p("- ten sam najlepszy kształt rozkładu emisji: %d z %d (TPA ma też modele PL, δ i trójskładnikowe)" % (same, len(G)))
        D = [r for r in G if math.isfinite(r["nf_d"])]
        if D:
            dd = np.array([r["nf_d"] - r["nf_t"] for r in D])
            p("- NF z naszego dynesty − TPA XIX: mediana %+.1f pp, |różnica| mediana %.1f pp (z %d)" % (
                100 * np.median(dd), 100 * np.median(np.abs(dd)), len(D)))
        p("\n### Największe różnice\n")
        p("| pulsar | NF TPA XIX [5–95%] | NF nasz [16–84%] | dynesty | automat | 3σ | f_i3 TPA / nasz | model TPA / nasz |")
        p("|---|---|---|---|---|---|---|---|")
        for r in sorted(G, key=lambda r: -abs(r["nf_o"] - r["nf_t"]))[:25]:
            p("| %s | %.1f%% [%.1f–%.1f] | %.1f%% ±%.1f | %s | %.1f%% | %.1f%% | %.2f / %.2f | %s / %s |" % (
                r["psr"], 100 * r["nf_t"], 100 * r["n5"], 100 * r["n95"], 100 * r["nf_o"], 100 * r["e_o"],
                "%.1f%%" % (100 * r["nf_d"]) if math.isfinite(r["nf_d"]) else "-", 100 * r["auto"], 100 * r["s3"],
                r["fi_t"], r["fi_o"], r["model_t"], r["model_o"][2:]))
    p("\n## Pulsary z literaturą\n")
    p("| pulsar | literatura | TPA XIX | nasz model | automat v6 | f_i3 TPA / nasz |")
    p("|---|---|---|---|---|---|")
    for r in R:
        if math.isfinite(r["lit"]) or r["psr"] == "J1850-0006":
            p("| %s | %s | %.1f%% ±%.1f | %.1f%%%s | %.1f%% | %.2f / %.2f |" % (
                r["psr"], "%.0f%%" % r["lit"] if math.isfinite(r["lit"]) else "-", 100 * r["nf_t"], 100 * r["e_t"],
                100 * r["nf_o"], "" if r["rel_o"] else " (niewiar.)", 100 * r["auto"], r["fi_t"], r["fi_o"]))
    miss = sorted(set(ours) - set(tpa))
    if miss:
        p("\nU nas, a nie w TPA XIX (albo błąd pobrania): %d: %s" % (len(miss), " ".join(miss[:40])))
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(OUT, "RAPORT_tpa19.md"), "w", encoding="utf-8").write(txt + "\n")
    keys = list(R[0].keys()) if R else []
    with open(os.path.join(OUT, "tpa19_compare.csv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        [w.writerow(r) for r in R]
    print("\nWyniki: %sRAPORT_tpa19.md, tpa19_compare.csv" % OUT)


if __name__ == "__main__":
    sys.exit(main())
