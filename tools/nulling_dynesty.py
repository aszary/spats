#!/usr/bin/env python3
"""Sprawdzenie modelu nulli nested samplingiem (dynesty), jak w TPA XIX (Keith+2026)
i w metodyce z publikacji w recenzji. Nasz model (tools/nulling_model.py) dopasowuje
metoda najwiekszej wiarygodnosci, wybiera rozklad przez BIC i liczy niepewnosc bootstrapem.
Tu ten sam model (ta sama wiarygodnosc, te same ograniczenia) liczony dynesty:
  - NF: mediana a posteriori i 16-84%  vs  nasze NF i przedzial z bootstrapu,
  - wybor rozkladu: ln Z (Keith: odrzucamy modele z ln Z - ln Z_best < -2.5)  vs  BIC.
Wejscie: fits/<psr>.npz z nulling_model.py (te same energie, RFI pominiete).

    NULLING_OUT=/home/psr/data/OUTPUT/nulling_model_v4 python3 tools/nulling_dynesty.py          # 10 pulsarow z literatura
    NULLING_OUT=... python3 tools/nulling_dynesty.py J0034-0721 J1727-2739                        # wybrane
    python3 tools/nulling_dynesty.py --test                                                       # dane syntetyczne
Wyniki: $NULLING_OUT/dynesty/summary.csv, RAPORT_dynesty.md
"""
import importlib.util
import json
import math
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy import special

try:
    import dynesty
    from dynesty import utils as dyfunc
except ImportError:
    sys.exit("Brak dynesty. Zainstaluj bez roota:  pip install --user dynesty   i uruchom ponownie.")

_spec = importlib.util.spec_from_file_location("nm", os.path.join(os.path.dirname(os.path.abspath(__file__)), "nulling_model.py"))
nm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nm)

OUT = nm.OUT
DEFAULT = ["J0034-0721", "J1049-5833", "J1537-4912", "J1727-2739", "J1851-0053", "J1926-0652",
           "J1944+1755", "J1946+1805", "J2037+1942", "J2253+1516"]
NLIVE = 250
DLOGZ = 0.5


def make_prior(model, e_off):
    """Szerokie priory jednostajne w parametrach modelu (jednostki <E_on>), jak [uzup.] w metodyce:
    NF ~ U(0,1), szum wokol E_off, polozenia/szerokosci skladnikow szeroko, wagi przez softmax."""
    m, s = float(np.mean(e_off)), float(np.std(e_off))
    lo, hi = [m - 0.5 * s, math.log(0.5 * s), 0.0], [m + 0.5 * s, math.log(2.0 * s), 1.0]   # nf w [0,1]
    box = {"G": ([0.0, math.log(0.01)], [5.0, math.log(5.0)]),
           "LN": ([-4.0, -3.0], [2.0, 1.5]),
           "PL": ([-5.0, math.log(0.01)], [2.0, math.log(5.0)])}
    for kind in model:
        lo += box[kind][0]; hi += box[kind][1]
    lo += [-5.0] * (len(model) - 1); hi += [5.0] * (len(model) - 1)
    lo, hi = np.array(lo), np.array(hi)

    def transform(u):
        x = lo + u * (hi - lo)
        x[2] = special.logit(min(max(x[2], 1e-9), 1 - 1e-9))     # NF ~ U(0,1) -> logit, jak w nm.unpack
        return x
    return transform, len(lo)


def run_one(args):
    psr, model, e_on, e_off, seed = args
    transform, nd = make_prior(model, e_off)

    def loglike(th):
        v = nm.nll(th, model, e_on, e_off)
        return -v if v < 1e11 else -1e100

    rng = np.random.default_rng(seed)
    s = dynesty.NestedSampler(loglike, transform, nd, nlive=NLIVE, bound="multi", sample="rwalk", rstate=rng)
    s.run_nested(dlogz=DLOGZ, print_progress=False)
    r = s.results
    w = np.exp(r.logwt - r.logz[-1])
    post = dyfunc.resample_equal(r.samples, w / w.sum(), rstate=rng)
    nf = special.expit(post[:, 2])
    return {"pulsar": psr, "model": "N+" + "+".join(model), "lnZ": float(r.logz[-1]), "lnZ_err": float(r.logzerr[-1]),
            "NF_med": float(np.median(nf)), "NF_lo": float(np.percentile(nf, 16)), "NF_hi": float(np.percentile(nf, 84)),
            "neval": int(sum(r.ncall))}


def analyse(psrs):
    jobs, ml = [], {}
    for psr in psrs:
        f = os.path.join(OUT, "fits", psr + ".npz")
        c = os.path.join(OUT, "fits", "%s.v%d.json" % (psr, nm.MODEL_VERSION))
        if not (os.path.isfile(f) and os.path.isfile(c)):
            print("  %s: brak wyniku nulling_model.py (fits/%s.npz) - pomijam" % (psr, psr))
            continue
        z = np.load(f)
        ml[psr] = json.load(open(c, encoding="utf-8"))
        for k, model in enumerate(nm.MODELS):
            jobs.append((psr, model, z["e_on"], z["e_off"], 1000 + k))
    print("DYNESTY: %d pulsarów x %d modeli = %d przebiegów, %d procesów" %
          (len(ml), len(nm.MODELS), len(jobs), min(32, os.cpu_count() or 1)), flush=True)
    with Pool(min(32, os.cpu_count() or 1)) as pool:
        res = pool.map(run_one, jobs)
    return ml, res


def report(ml, res, outdir):
    os.makedirs(outdir, exist_ok=True)
    by = {}
    for r in res:
        by.setdefault(r["pulsar"], []).append(r)
    L = []
    p = L.append
    p("# Model nulli: nested sampling (dynesty) a nasze dopasowanie (największa wiarygodność + BIC + bootstrap)")
    p("")
    p("Ta sama wiarygodność i ograniczenia co w nulling_model.py; priory jednostajne; nlive=%d, dlogz=%.1f." % (NLIVE, DLOGZ))
    p("")
    p("| pulsar | NF nasze [bootstrap 16–84%] | NF dynesty [16–84%] | różnica | model BIC | model ln Z | modele w ΔlnZ > −2.5 | literatura |")
    p("|---|---|---|---|---|---|---|---|")
    rows, diffs, same = [], [], 0
    for psr in sorted(by):
        R = sorted(by[psr], key=lambda r: -r["lnZ"])
        best = R[0]
        ok_models = [r["model"] for r in R if r["lnZ"] - best["lnZ"] > -2.5]
        m = ml[psr]
        # NF z dynesty dla modelu wybranego przez BIC (porownanie tego samego modelu) i dla najlepszego wg ln Z
        same_model = next((r for r in R if r["model"] == m["model"]), best)
        d = 100 * (same_model["NF_med"] - m["NF_model"])
        diffs.append(d)
        same += best["model"] == m["model"]
        lit = m.get("NF_lit")
        p("| %s | %.1f%% [%.1f–%.1f] | %.1f%% [%.1f–%.1f] | %+.1f pp | %s | %s (NF %.1f%%) | %s | %s |" % (
            psr, 100 * m["NF_model"], 100 * m.get("NF_model_lo", float("nan")), 100 * m.get("NF_model_hi", float("nan")),
            100 * same_model["NF_med"], 100 * same_model["NF_lo"], 100 * same_model["NF_hi"], d,
            m["model"], best["model"], 100 * best["NF_med"], ", ".join(ok_models), "%.0f%%" % lit if lit else "-"))
        for r in R:
            rows.append([psr, r["model"], "%.2f" % r["lnZ"], "%.2f" % (r["lnZ"] - best["lnZ"]), "%.4f" % r["NF_med"],
                         "%.4f" % r["NF_lo"], "%.4f" % r["NF_hi"], r["neval"], m["model"], "%.4f" % m["NF_model"]])
    if diffs:
        p("")
        p("- NF dynesty − NF nasze (ten sam model): mediana %+.1f pp, maks. |różnica| %.1f pp (z %d pulsarów)" %
          (float(np.median(diffs)), float(np.max(np.abs(diffs))), len(diffs)))
        p("- ten sam najlepszy model wg BIC i wg ln Z: %d z %d" % (same, len(diffs)))
    with open(os.path.join(outdir, "summary.csv"), "w", encoding="utf-8") as fh:
        fh.write("pulsar,model,lnZ,dlnZ,NF_med,NF_lo,NF_hi,neval,model_BIC,NF_ML\n")
        for r in rows:
            fh.write(",".join(str(x) for x in r) + "\n")
    with open(os.path.join(outdir, "RAPORT_dynesty.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))
    print("\nWyniki: %s/summary.csv, RAPORT_dynesty.md" % outdir)


def selftest():
    """Dane syntetyczne ze znanym NF: dynesty i nasze dopasowanie na tych samych danych."""
    print("TEST dynesty: dane syntetyczne")
    for nf_true, snr in ((0.3, 30), (0.5, 5)):
        eon, eoff, isnull = nm.synth(nf_true, n=800, snr=snr)
        norm = eon.mean()
        e_on, e_off = eon / norm, eoff / norm
        model = ("LN",)
        r = nm.fit(model, e_on, e_off)
        d = run_one(("synt", model, e_on, e_off, 1))
        print("  NF prawdziwe %.2f (w próbce %.3f), S/N %d: nasze %.3f, dynesty %.3f [%.3f–%.3f], ln Z %.1f" %
              (nf_true, isnull.mean(), snr, float(special.expit(r.x[2])), d["NF_med"], d["NF_lo"], d["NF_hi"], d["lnZ"]))


def main():
    if "--test" in sys.argv:
        selftest()
        return
    want = [a for a in sys.argv[1:] if not a.startswith("-")] or DEFAULT
    ml, res = analyse(want)
    if res:
        report(ml, res, os.path.join(OUT, "dynesty"))


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
