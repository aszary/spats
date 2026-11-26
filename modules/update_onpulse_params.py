#!/usr/bin/env python3
"""
update_onpulse_params.py

Updates on-pulse regions in `params.json` files for pulsars based on `input/onpulsefixing.csv`:
- MP mode: sets `bin_st` to lowest bin value, `bin_end` to highest bin value across components (C1-C4).
- IP mode: sets `IP_start` to lowest bin value, `IP_end` to highest bin value across components (C1-C4).
- Creates a backup `params_pre_reonpuls.json` in the same directory before overwriting `params.json`.
"""

import os
import sys
import csv
import json
import shutil
from typing import Dict, Any, Optional

def parse_float_safe(val: str) -> Optional[float]:
    if not val:
        return None
    s = val.strip().replace(',', '.')
    if not s:
        return None
    try:
        v = float(s)
        return v if v != 0.0 else None
    except ValueError:
        return None

def parse_onpulse_csv(csv_path: str = "input/onpulsefixing.csv") -> Dict[str, Dict[str, int]]:
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    results = {}
    # Component column index triplets in onpulsefixing.csv (0-indexed):
    # C1: reg=16, st=17, en=18
    # C2: reg=19, st=20, en=21
    # C3: reg=22, st=23, en=24
    # C4: reg=25, st=26, en=27
    comp_cols = [
        (16, 17, 18),
        (19, 20, 21),
        (22, 23, 24),
        (25, 26, 27)
    ]

    with open(csv_path, mode='r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader, None)

        for row in reader:
            if not row or not row[0].strip():
                continue
            jname = row[0].strip()

            mp_starts, mp_ends = [], []
            ip_starts, ip_ends = [], []

            # FFT length used (column index 2: FFT_length_used)
            fft_len = parse_float_safe(row[2]) if len(row) > 2 else None

            for reg_idx, st_idx, en_idx in comp_cols:
                if reg_idx >= len(row):
                    continue
                reg = row[reg_idx].strip().upper()
                st = parse_float_safe(row[st_idx]) if st_idx < len(row) else None
                en = parse_float_safe(row[en_idx]) if en_idx < len(row) else None

                if reg == "MP":
                    if st is not None: mp_starts.append(st)
                    if en is not None: mp_ends.append(en)
                elif reg == "IP":
                    if st is not None: ip_starts.append(st)
                    if en is not None: ip_ends.append(en)

            entry = {}
            if mp_starts and mp_ends:
                entry["bin_st"] = int(round(min(mp_starts)))
                entry["bin_end"] = int(round(max(mp_ends)))

            if ip_starts and ip_ends:
                entry["IP_start"] = int(round(min(ip_starts)))
                entry["IP_end"] = int(round(max(ip_ends)))

            if fft_len is not None:
                entry["fft"] = int(round(fft_len))
                entry["nfft"] = int(round(fft_len))

            if entry:
                results[jname] = entry

    return results

def update_single_params(params_path: str, onpulse_info: Dict[str, int], dry_run: bool = False, only_fft: bool = True) -> bool:
    if not os.path.isfile(params_path):
        return False

    with open(params_path, 'r', encoding='utf-8') as f:
        params = json.load(f)

    # Only update MP/IP if only_fft is False
    if not only_fft:
        # Apply MP updates
        if "bin_st" in onpulse_info and "bin_end" in onpulse_info:
            params["bin_st"] = onpulse_info["bin_st"]
            params["bin_end"] = onpulse_info["bin_end"]

        # Apply IP updates
        if "IP_start" in onpulse_info and "IP_end" in onpulse_info:
            params["IP_start"] = onpulse_info["IP_start"]
            params["IP_end"] = onpulse_info["IP_end"]

    # Apply fft / nfft updates
    if "fft" in onpulse_info:
        params["fft"] = onpulse_info["fft"]
        params["nfft"] = onpulse_info["nfft"]

    backup_path = os.path.join(os.path.dirname(params_path), "params_pre_reonpuls.json")

    if not dry_run:
        shutil.copyfile(params_path, backup_path)
        with open(params_path, 'w', encoding='utf-8') as f:
            json.dump(params, f, indent=4)

    return True

def update_params_files(search_dirs=None, csv_path="input/onpulsefixing.csv", dry_run=False, only_fft=True):
    if search_dirs is None:
        search_dirs = ["/home/psr/output", "output"]

    print(f"Loading data from {csv_path} (only_fft={only_fft}) ...")
    onpulse_data = parse_onpulse_csv(csv_path)
    print(f"Loaded {len(onpulse_data)} pulsar definitions.")

    updated_count = 0
    not_found = []

    for jname, info in onpulse_data.items():
        found = False
        folder_variants = [jname, f"{jname}_16", f"{jname}_32"]

        for sdir in search_dirs:
            if not os.path.isdir(sdir):
                continue
            for fvar in folder_variants:
                candidate_dir = os.path.join(sdir, fvar)
                pfile = os.path.join(candidate_dir, "params.json")
                if os.path.isfile(pfile):
                    update_single_params(pfile, info, dry_run=dry_run, only_fft=only_fft)
                    updated_count += 1
                    found = True
                    if only_fft:
                        print(f"  [OK] Updated {pfile} -> fft={info.get('fft')}")
                    else:
                        ip_str = f", IP_start={info['IP_start']}, IP_end={info['IP_end']}" if "IP_start" in info else ""
                        print(f"  [OK] Updated {pfile} -> bin_st={info.get('bin_st')}, bin_end={info.get('bin_end')}, fft={info.get('fft')}{ip_str}")

        if not found:
            not_found.append(jname)

    print("\nSummary:")
    print(f"  Updated params.json files: {updated_count}")
    print(f"  Pulsars not located in search dirs: {len(not_found)}")
    return updated_count, not_found

if __name__ == "__main__":
    dirs = sys.argv[1:] if len(sys.argv) > 1 else ["/home/psr/output", "output"]
    update_params_files(search_dirs=dirs, only_fft=True)
