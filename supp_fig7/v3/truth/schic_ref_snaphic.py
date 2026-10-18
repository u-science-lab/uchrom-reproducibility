"""The independent reference of R7: the authors' SnapHiC (github.com/HuMingLab/SnapHiC @ 782dc2b, GPL-3.0; environment
venv-snaphic, the newest numpy / pandas / scipy for Python 3.12, no U-Chrom) on the 50 cells of
data_v3/nagano2017_dip_chr19 with its defaults, chr19, 8 processes (threaded mode); compared with U-Chrom's calls.

    python schic_ref_snaphic.py --data3 /data_v3 --shared /shared --work DIR --out schic_ref_snaphic.json
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

SNAP = Path("/opt/venv/src/SnapHiC")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data3", type=Path, required=True)
    ap.add_argument("--shared", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--nproc", type=int, default=8)
    a = ap.parse_args()
    d = a.data3 / "nagano2017_dip_chr19"
    outdir = Path(f"/tmp/snaphic_out_{os.environ.get('SLURM_JOB_ID', os.getpid())}")
    if outdir.exists():                                  # SnapHiC resumes from (and trips over) an old output folder
        import shutil
        shutil.rmtree(outdir)
    cmd = [sys.executable, str(SNAP / "snap.py"), "-i", str(d / "contacts"), "-s", "contacts.txt", "-o", str(outdir),
           "-c", "2", "4", "-p", "3", "5", "-l", str(SNAP / "ext/mm10.chrom.sizes"), "-g", "mm10", "--chrom", "chr19",
           "--filter-file", str(a.shared / "snaphic/mm10_filter_regions.txt"), "--threaded", "-n", str(a.nproc),
           "--no-hic", "--no-cool", "--prefix", "ref"]
    t0 = time.perf_counter()
    p = subprocess.run(cmd, cwd=str(SNAP), capture_output=True, text=True)
    secs = time.perf_counter() - t0
    head = subprocess.run(["git", "-C", str(SNAP), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    res = {"commit": head, "exit": p.returncode, "runtime_s": secs, "cmd": " ".join(cmd)}
    summits = outdir / "postprocessed" / "ref.postprocessed.summits.bedpe"
    if p.returncode != 0 or not summits.exists():
        res["stderr"] = p.stderr[-4000:]
        res["stdout"] = p.stdout[-2000:]
        log = outdir / "snapHiC.log"
        if log.exists():
            res["log_tail"] = log.read_text()[-3000:]
    else:
        t = pd.read_csv(summits, sep="\t") if summits.stat().st_size else pd.DataFrame()
        res["n_loops"] = int(len(t))
        t.to_csv(a.work / "schic_R7_snaphic_summits.tsv", sep="\t", index=False)
        ours = a.work / "schic_R7_uchrom_loops.tsv"
        if ours.exists() and len(t):
            u = pd.read_csv(ours, sep="\t")
            c = list(t.columns)
            theirs = set(zip(t[c[1]].astype(int), t[c[4]].astype(int)))
            mine = set(zip(u["start1"].astype(int), u["start2"].astype(int)))
            near = sum(any(abs(x - x2) <= 20_000 and abs(y - y2) <= 20_000 for x2, y2 in mine) for x, y in theirs)
            res["vs_uchrom"] = {"n_uchrom": len(u), "identical": len(theirs & mine),
                                "theirs_within_20kb_of_ours": int(near)}
    a.out.write_text(json.dumps(res, indent=1))
    print(json.dumps(res)[:2000])


if __name__ == "__main__":
    main()
