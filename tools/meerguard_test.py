#!/usr/bin/env python3
"""Test MeerGuard (kroki 1-5): srodowisko, dane 16-kanalowe, szablon 2D,
warianty A/B/C, wagi -> flagi. Porownanie z automatem: tools/meerguard_compare.jl.

Nie zmienia niczego w data/new, pazur ani pazur2: dane sa laczone (psradd) do
/home/psr/data/OUTPUT/meerguard_test/<psr>/<psr>.ar, MeerGuard pisze do <psr>/<A|B|C>/.
Ponowne uruchomienie pomija gotowe kroki (istniejace pliki).

    python3 tools/meerguard_test.py            # wszystko
    python3 tools/meerguard_test.py J0630-2834 # wybrane pulsary
"""
import concurrent.futures as cf
import glob
import hashlib
import json
import os
import re
import subprocess
import sys
import time

ROOT = "/home/psr/data/OUTPUT/meerguard_test"
DATAROOT = "/home/psr/data/new"
PAZUR2 = "/home/psr/data/OUTPUT/pazur2"
MGDIR = os.path.expanduser("~/src/MeerGuard")
MGURL = "https://github.com/danielreardon/MeerGuard"
PSRS = ["J0630-2834",   # bardzo jasny; automat wycina kanal 15
        "J0034-0721",   # duzo nulli
        "J1057-5226",   # interpuls
        "J1843+2024",   # slaby, silne RFI; metoda reczna wycina wiecej
        "J1543-0620",   # automat wycina kanal 3
        "J1727-2739",   # nulle + kanaly 6, 15
        "J1850-0006"]   # slaby nuller
MG16_CFG = """surgical_default_params = 'template=None,chan_breakpoints=None,chan_numpieces=1,chan_order=1,chanthresh=7,subint_breakpoints=None,subint_numpieces=1,subint_order=1,subintthresh=7,plot=None'
rcvrstd_default_params = 'badchans=None,badfreqs=None,badsubints=None,response=None,trimbw=0,trimfrac=0,trimnum=0'
bandwagon_default_params = 'badchantol=0.95,badsubtol=0.95'
"""


def log(msg, f=None):
    print(msg, flush=True)
    if f is not None:
        with open(f, "a") as fh:
            fh.write(msg + "\n")


def stop(msg):
    log("STOP: " + msg, os.path.join(ROOT, "status.txt"))
    sys.exit(1)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd, cwd=None, logfile=None):
    t0 = time.time()
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    dt = time.time() - t0
    if logfile:
        with open(logfile, "w") as fh:
            fh.write("$ " + " ".join(cmd) + "\n")
            fh.write(p.stdout)
            fh.write("\n[exit %d, %.1f s]\n" % (p.returncode, dt))
    return p.returncode, p.stdout, dt


def zaps_ranges(z):
    """'0 5-9 q 49'' -> [(0,0),(5,9)] (indeksy subintow od 0, jak paz -w)"""
    out = []
    for tok in (z or "").split():
        m = re.match(r"^(\d+)(?:-(\d+))?", tok)
        if m:
            lo = int(m.group(1))
            out.append((lo, int(m.group(2)) if m.group(2) else lo))
    return out


# ---------------------------------------------------------------- 1. srodowisko
def find_psrchive_path():
    """Katalog z psrchive.py + _psrchive.so dla TEJ wersji Pythona (PSRCHIVE instaluje
    do .../site-packages, a Debian/Ubuntu szuka w dist-packages)."""
    ver = "python%d.%d" % sys.version_info[:2]
    for base in ("/usr/local/lib", "/usr/lib", os.path.expanduser("~/.local/lib")):
        d = os.path.join(base, ver, "site-packages")
        if os.path.isfile(os.path.join(d, "psrchive.py")) and glob.glob(os.path.join(d, "_psrchive*.so")):
            return d
    return None


def step_env(meta):
    os.makedirs(ROOT, exist_ok=True)
    try:
        import psrchive
    except Exception as e:
        d = find_psrchive_path()
        if d is None:
            stop("brak bindingow PSRCHIVE dla Pythona %d.%d (import psrchive: %s; nie ma ich tez w "
                 "/usr/local/lib/python%d.%d/site-packages). MeerGuard ich wymaga." % (sys.version_info[:2] + (e,) + sys.version_info[:2]))
        sys.path.insert(0, d)
        os.environ["PYTHONPATH"] = d + os.pathsep + os.environ.get("PYTHONPATH", "")   # dla procesow MeerGuard
        meta["pythonpath_added"] = d
        log("psrchive nie byl na sciezce Pythona - dodaje %s (sys.path i PYTHONPATH dla MeerGuard)" % d)
        try:
            import psrchive
        except Exception as e2:
            stop("psrchive z %s nie laduje sie: %s" % (d, e2))
    meta["psrchive"] = psrchive.__file__
    meta["python"] = sys.version.split()[0]
    for mod in ("numpy", "scipy"):
        try:
            meta[mod] = __import__(mod).__version__
        except Exception as e:
            stop("brak pakietu %s (%s): pip install --user %s" % (mod, e, mod))
    if not os.path.isdir(MGDIR):
        os.makedirs(os.path.dirname(MGDIR), exist_ok=True)
        rc, out, _ = run(["git", "clone", MGURL, MGDIR])
        if rc:
            stop("git clone %s nie powiodl sie:\n%s" % (MGURL, out))
    meta["meerguard_commit"] = run(["git", "rev-parse", "HEAD"], cwd=MGDIR)[1].strip()
    rc, _, _ = run([sys.executable, "-c", "import coast_guard"])
    if rc:
        rc, out, _ = run([sys.executable, "-m", "pip", "install", "--user", "."], cwd=MGDIR,
                         logfile=os.path.join(ROOT, "pip_install.txt"))
        if rc:
            stop("pip install --user MeerGuard nie powiodl sie, zob. %s/pip_install.txt" % ROOT)
    script = os.path.join(MGDIR, "clean_archive.py")
    if not os.path.isfile(script):
        stop("brak %s" % script)
    meta["clean_archive"] = script
    rc, out, _ = run([sys.executable, script, "-h"], logfile=os.path.join(ROOT, "clean_archive_help.txt"))
    if rc:
        stop("clean_archive.py -h nie dziala, zob. %s/clean_archive_help.txt" % ROOT)
    missing = [o for o in ("-T", "-C", "-c", "-s", "-bc", "-bs", "-ag", "-e", "-O") if o + " " not in out and o + "," not in out]
    meta["help_missing_options"] = missing
    if missing:
        log("UWAGA: w clean_archive.py -h brak opcji: %s" % missing)
    with open(os.path.join(ROOT, "mg16.cfg"), "w") as fh:
        fh.write(MG16_CFG)
    log("srodowisko OK: psrchive %s, MeerGuard %s" % (meta["psrchive"], meta["meerguard_commit"][:10]))


# ---------------------------------------------------------- 2.-3. dane i szablon
def step_data(psr, meta):
    """psradd plikow .spCf16 (pierwsza obserwacja, jak w automacie) + szablon 2D."""
    import psrchive
    st = meta["pulsars"].setdefault(psr, {})
    pf = os.path.join(PAZUR2, psr, "params.json")
    if not os.path.isfile(pf):
        st["skip"] = "brak %s" % pf
        return False
    p = json.load(open(pf))
    obs = sorted(d for d in os.listdir(os.path.join(DATAROOT, psr)) if os.path.isdir(os.path.join(DATAROOT, psr, d)))
    if not obs:
        st["skip"] = "brak katalogu obserwacji"
        return False
    indir = os.path.join(DATAROOT, psr, obs[0])
    files = glob.glob(os.path.join(indir, "*.spCf16"))
    if not files:
        st["skip"] = "brak plikow .spCf16 w %s" % indir
        return False

    def start(f):
        m = re.search(r"_(\d+)-(\d+)\.spCf16$", f)
        return int(m.group(1)) if m else 1 << 60
    files.sort(key=start)
    d = os.path.join(ROOT, psr)
    os.makedirs(d, exist_ok=True)
    ar_path = os.path.join(d, psr + ".ar")
    if not os.path.isfile(ar_path):
        rc, out, _ = run(["psradd", "-o", ar_path] + files, logfile=os.path.join(d, "psradd.txt"))
        if rc:
            st["skip"] = "psradd nie powiodl sie (zob. %s/psradd.txt)" % d
            return False
    ar = psrchive.Archive_load(ar_path)
    nsub, nchan, npol, nbin = ar.get_nsubint(), ar.get_nchan(), ar.get_npol(), ar.get_nbin()
    n_spcf = (p.get("zap_cmp") or {}).get("n_pulses")
    st.update(obs=obs[0], nfiles=len(files), nsub=nsub, nchan=nchan, npol=npol, nbin=nbin, n_spcf=n_spcf,
              bin_st=p.get("bin_st"), bin_end=p.get("bin_end"),
              zaps2=p.get("zaps2"), zaps2_rfi=p.get("zaps2_rfi"), zaps2_null=p.get("zaps2_null"),
              zap_chans16=p.get("zap_chans16"), md5_in=md5(ar_path))
    if n_spcf is not None and nsub != n_spcf:
        st["skip"] = "plik 16-kanalowy ma %d subintow, .spCF %d - bez mapowania po czasie nie porownuje" % (nsub, n_spcf)
        return False
    if nchan != 16:
        log("UWAGA %s: %d kanalow (oczekiwane 16)" % (psr, nchan))

    tpl = os.path.join(d, "template.ar")
    if not os.path.isfile(tpl):
        t = psrchive.Archive_load(ar_path)
        rows = [i for lo, hi in zaps_ranges(p.get("zaps2")) for i in range(lo, hi + 1) if 0 <= i < nsub]
        chans = [int(c) for c in (p.get("zap_chans16") or "").split() if c.isdigit() and int(c) < nchan]
        for i in rows:
            integ = t.get_Integration(i)
            for c in range(nchan):
                integ.set_weight(c, 0.0)
        if chans:
            for i in range(nsub):
                integ = t.get_Integration(i)
                for c in chans:
                    integ.set_weight(c, 0.0)
        t.tscrunch()
        t.pscrunch()
        t.unload(tpl)
        st.update(template_zeroed_subints=len(rows), template_zeroed_chans=chans)
    tt = psrchive.Archive_load(tpl)
    st["template_shape"] = [tt.get_nsubint(), tt.get_nchan(), tt.get_npol(), tt.get_nbin()]
    if tt.get_nsubint() != 1:
        st["skip"] = "szablon ma %d subintow (musi byc 1)" % tt.get_nsubint()
        return False
    return True


# -------------------------------------------------------------- 4. uruchomienia
def variants(psr):
    d = os.path.join(ROOT, psr)
    ar, tpl, cfg = os.path.join(d, psr + ".ar"), os.path.join(d, "template.ar"), os.path.join(ROOT, "mg16.cfg")
    return {
        "A": [ar, "-T", tpl, "-C", cfg, "-c", "7", "-s", "7", "-bc", "0.95", "-bs", "0.95", "-e", "mg", "-O", os.path.join(d, "A")],
        "B": [ar, "-T", tpl, "-e", "mg", "-O", os.path.join(d, "B")],
        "C": [ar, "-T", tpl, "-C", cfg, "-ag", "-e", "mg", "-O", os.path.join(d, "C")],
    }


def run_variant(psr, v, args, script):
    d = os.path.join(ROOT, psr, v)
    os.makedirs(d, exist_ok=True)
    out = os.path.join(d, psr + ".mg")
    if os.path.isfile(out):
        return psr, v, 0, None, "juz jest"
    cmd = [sys.executable, script] + args
    rc, _, dt = run(cmd, cwd=d, logfile=os.path.join(d, "log.txt"))
    return psr, v, rc, dt, ("brak pliku wyjsciowego" if rc == 0 and not os.path.isfile(out) else "")


# ------------------------------------------------------------------ 5. wyniki
def step_flags(psr, meta):
    import numpy as np
    import psrchive
    st = meta["pulsars"][psr]
    d = os.path.join(ROOT, psr)
    w0 = psrchive.Archive_load(os.path.join(d, psr + ".ar")).get_weights() > 0
    np.savetxt(os.path.join(d, "live_in.txt"), w0.astype(int), fmt="%d")
    res = {}
    for v in ("A", "B", "C"):
        out = os.path.join(d, v, psr + ".mg")
        if not os.path.isfile(out):
            res[v] = {"error": "brak %s" % out}
            continue
        w1 = psrchive.Archive_load(out).get_weights() > 0
        new = w0 & ~w1                                   # nowe zera (bez zer z wejscia)
        np.savetxt(os.path.join(d, v, "newzero.txt"), new.astype(int), fmt="%d")
        live_sub = w0.any(axis=1)
        live_ch = w0.any(axis=0)
        sub_full = np.flatnonzero(live_sub & ~w1.any(axis=1))
        ch_full = np.flatnonzero(live_ch & ~w1.any(axis=0))
        sub_any = np.flatnonzero(new.any(axis=1))
        nlive = int(w0.sum())
        chan_frac = [float(new[:, c].sum()) / max(int(w0[:, c].sum()), 1) for c in range(new.shape[1])]
        res[v] = {"pairs": int(new.sum()), "frac_pairs": float(new.sum()) / max(nlive, 1),
                  "subints_full": sub_full.tolist(), "chans_full": ch_full.tolist(),
                  "subints_any": sub_any.tolist(), "chan_frac": chan_frac}
    st["flags"] = res
    st["md5_after"] = md5(os.path.join(d, psr + ".ar"))
    st["input_unchanged"] = st["md5_after"] == st.get("md5_in")


def main():
    psrs = sys.argv[1:] or PSRS
    os.makedirs(ROOT, exist_ok=True)
    mf = os.path.join(ROOT, "meta.json")
    meta = json.load(open(mf)) if os.path.isfile(mf) else {}
    meta.setdefault("pulsars", {})
    step_env(meta)
    ok = []
    for psr in psrs:
        try:
            good = step_data(psr, meta)
        except Exception as e:
            meta["pulsars"].setdefault(psr, {})["skip"] = "blad danych/szablonu: %r" % e
            good = False
        log("%s: %s" % (psr, "dane i szablon OK" if good else "POMINIETY - " + meta["pulsars"][psr]["skip"]))
        if good:
            ok.append(psr)
    json.dump(meta, open(mf, "w"), indent=1)

    jobs = [(psr, v, a) for psr in ok for v, a in variants(psr).items()]
    log("uruchamiam MeerGuard: %d przebiegow (%d pulsarow x 3 warianty)" % (len(jobs), len(ok)))
    with cf.ThreadPoolExecutor(max_workers=min(7, max(1, len(jobs)))) as ex:
        futs = [ex.submit(run_variant, psr, v, a, meta["clean_archive"]) for psr, v, a in jobs]
        for f in cf.as_completed(futs):
            psr, v, rc, dt, note = f.result()
            meta["pulsars"][psr].setdefault("runs", {})[v] = {"exit": rc, "seconds": dt, "note": note,
                                                             "cmd": " ".join(["python3", "clean_archive.py"] + variants(psr)[v])}
            log("  %s %s: exit %s%s %s" % (psr, v, rc, "" if dt is None else ", %.0f s" % dt, note))
    for psr in ok:
        try:
            step_flags(psr, meta)
            fl = meta["pulsars"][psr]["flags"]
            log("%s: pary A/B/C = %s   plik wejsciowy bez zmian: %s" % (
                psr, " / ".join("%.2f%%" % (100 * fl[v]["frac_pairs"]) if "frac_pairs" in fl[v] else "-" for v in "ABC"),
                meta["pulsars"][psr]["input_unchanged"]))
        except Exception as e:
            meta["pulsars"][psr]["flags_error"] = repr(e)
            log("%s: blad odczytu wag: %r" % (psr, e))
    json.dump(meta, open(mf, "w"), indent=1)
    log("Gotowe: %s. Dalej: SPATS_NO_MAIN=1 julia tools/meerguard_compare.jl" % mf)


if __name__ == "__main__":
    main()
