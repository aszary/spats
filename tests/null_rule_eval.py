#!/usr/bin/env python3
"""Test wyboru reguly wycinania nulli (krok 2 z 2, po tests/null_rule_gen.jl): model energii na
danych syntetycznych ze znana prawda i porownanie wariantow wycinania nulli.
    python3 tests/null_rule_eval.py <katalog>
Miary (srednio po zestawach): kompletnosc = wyciete prawdziwe nulle / wszystkie nulle; czystosc =
nulle / wyciete; emisja wycieta = wyciete impulsy z emisja / wszystkie z emisja; slabe wyciete = to samo
dla emisji < 0.3 sredniej; bias profilu = zmiana sredniej emisji pozostalych impulsow (0 = profil
bez znieksztalcenia). Wynik: <katalog>/null_rule.json, null_rule.txt"""
import glob
import importlib.util
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

np.seterr(all="ignore")
out = os.environ.get("NR_OUT") or sys.argv[1]
os.environ["NULLING_OUT"] = out
os.environ["NR_OUT"] = out
spec = importlib.util.spec_from_file_location("nm", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools", "nulling_model.py"))
nm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(nm)
nm.OUT = out
VARIANTS = ["automat (v5)", "automat bez marg.", "model P>=0.5 (jedna reguła)", "model P>=0.8", "model P>=0.95",
            "automat + P>=0.95", "automat bez marg. + P>=0.95", "automat bez marg. + P>=0.8", "automat + P>=0.8"]


def one(f):
    psr = os.path.basename(f)[:-4]
    r = nm.analyse(f)
    z = np.load(os.path.join(out, "fits", psr + ".npz"))
    idx, pn = z["idx"], z["pnull"]
    tr = np.loadtxt(os.path.join(out, psr + ".truth"), delimiter=",")[idx]
    isn, amp, am, an = tr[:, 0] > 0.5, tr[:, 1], tr[:, 2] > 0.5, tr[:, 3] > 0.5
    weak = (~isn) & (amp < 0.3)
    cuts = dict(zip(VARIANTS, [am, an, pn >= 0.5, pn >= 0.8, pn >= 0.95, am | (pn >= 0.95), an | (pn >= 0.95),
                               an | (pn >= 0.8), am | (pn >= 0.8)]))
    res = {"psr": psr, "nf_true": float(isn.mean()), "nf_model": r["NF_model"], "f_i3": r["f_i3"]}
    for name, cut in cuts.items():
        res[name] = dict(kompletnosc=float((cut & isn).sum() / max(isn.sum(), 1)),
                         czystosc=float((cut & isn).sum() / max(cut.sum(), 1)) if cut.sum() else float("nan"),
                         emisja_wycieta=float((cut & ~isn).sum() / max((~isn).sum(), 1)),
                         slabe_wyciete=float((cut & weak).sum() / max(weak.sum(), 1)),
                         bias=float(amp[~cut & ~isn].mean() / amp[~isn].mean() - 1),
                         # to samo oszacowane z P_null modelu (jak tools/null_rules_real.jl na prawdziwych danych)
                         kompletnosc_szac=float(pn[cut].sum() / max(pn.sum(), 1e-12)),
                         emisja_szac=float((1 - pn[cut]).sum() / max((1 - pn).sum(), 1e-12)))
    return res


def main():
    files = sorted(glob.glob(os.path.join(out, "energies", "*.csv")))
    with Pool(min(16, os.cpu_count() or 1)) as p:
        R = p.map(one, files)
    keys = ["kompletnosc", "czystosc", "emisja_wycieta", "slabe_wyciete", "bias", "kompletnosc_szac", "emisja_szac"]
    L, summ = [], {}
    scen = sorted(set(r["psr"].rsplit("_", 1)[0] for r in R))
    for s in scen:
        rr = [r for r in R if r["psr"].startswith(s + "_")]
        L.append("\n== %s: NF prawdziwe %.3f, NF model %.3f, f_i3 %.2f" % (
            s[1:], np.mean([r["nf_true"] for r in rr]), np.mean([r["nf_model"] for r in rr]), np.mean([r["f_i3"] for r in rr])))
        L.append("   %-30s %11s %9s %14s %13s %12s" % ("wariant", "kompletność", "czystość", "emisja wycięta", "słabe wycięte", "bias profilu"))
        for v in VARIANTS:
            m = {k: float(np.nanmean([r[v][k] for r in rr])) for k in keys}
            summ.setdefault(v, {})[s[1:]] = m
            L.append("   %-30s %10.1f%% %8.1f%% %13.2f%% %12.1f%% %+11.2f%%" % (v, 100 * m["kompletnosc"], 100 * m["czystosc"],
                     100 * m["emisja_wycieta"], 100 * m["slabe_wyciete"], 100 * m["bias"]))
    L.append("\n== ŚREDNIO po scenariuszach (bez 'modulowany' w kompletności: brak nulli)")
    for v in VARIANTS:
        sc = summ[v]
        k = np.mean([sc[s]["kompletnosc"] for s in sc if s != "modulowany"])
        e = np.mean([sc[s]["emisja_wycieta"] for s in sc]); b = max(sc[s]["bias"] for s in sc)
        summ[v]["_srednio"] = {"kompletnosc": float(k), "emisja_wycieta": float(e), "bias_max": float(b)}
        ks = np.mean([sc[s]["kompletnosc_szac"] for s in sc if s not in ("modulowany", "_srednio")])
        es = np.mean([sc[s]["emisja_szac"] for s in sc if s != "_srednio"])
        summ[v]["_srednio"].update({"kompletnosc_szac": float(ks), "emisja_szac": float(es)})
        L.append("   %-30s kompletność %5.1f%% (z modelu %5.1f%%)  emisja wycięta %.2f%% (z modelu %.2f%%)  największy bias profilu %+.2f%%" %
                 (v, 100 * k, 100 * ks, 100 * e, 100 * es, 100 * b))
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(out, "null_rule.txt"), "w", encoding="utf-8").write(txt + "\n")
    clean = lambda o: ({k: clean(v) for k, v in o.items()} if isinstance(o, dict) else
                       (None if isinstance(o, float) and not np.isfinite(o) else o))     # NaN -> null (poprawny JSON)
    json.dump(clean(summ), open(os.path.join(out, "null_rule.json"), "w", encoding="utf-8"), indent=1)


if __name__ == "__main__":
    main()
