"""The independent reference of R1 / R2: the authors' NucDynamics (github.com/tjs23/nuc_dynamics @ 00eaf0b, its own
environment venv-nucdyn: numpy + Cython, no U-Chrom, no scipy) on the GEO contacts of a Stevens et al. 2017 cell,
1-Mb final particles, compared with the published model 1 as the question asks.

    python schic_ref_nucdyn.py --cell-dir DIR --label R1 --out schic_ref_nucdyn.json [--seeds 0,1]

The GEO table (``chr_A pos_A chr_B pos_B``) is written as a gzipped NCC file (one active, unambiguous contact per
line, '+' strand, the position as fragment start and end; gzipped because the code opens plain files with mode 'rU',
which Python >= 3.11 rejects, and gzip files with 'rt'), then ``nuc_dynamics.py <ncc> -m 1 -s 8 4 2 1 -ran <seed>``
(the README's command with the question's particle sizes; defaults otherwise: 500 temperature steps x 100 dynamics
steps per stage).  Its N3D output (chr1 particles) is averaged into 1-Mb bins and the pairwise distances are
compared with those of the published PDB model 1 (Spearman, ties averaged, numpy only).
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

ND = Path("/opt/venv/src/nuc_dynamics")


def write_ncc(pairs: Path, ncc: Path) -> int:
    n = 0
    with gzip.open(pairs, "rt") as fh, gzip.open(ncc, "wt") as out:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            ca, pa, cb, pb = line.split()[:4]
            out.write(f"{ca}\t{pa}\t{pa}\t{pa}\t{pa}\t+\t{cb}\t{pb}\t{pb}\t{pb}\t{pb}\t+\t1.1\t{n}\t0\n")
            n += 1
    return n


def read_n3d(path: Path, chrom: str = "chr1"):
    pos, xyz, cur = [], [], None
    for line in open(path):
        f = line.split()
        if len(f) == 3:
            cur = f[0]
        elif cur == chrom and len(f) >= 4:
            pos.append(int(f[0]))
            xyz.append([float(v) for v in f[1:4]])
    return np.array(pos, np.int64), np.array(xyz)


def read_pdb(path: Path, model: int = 1, chain: str = "a"):
    pos, xyz, m = [], [], 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("MODEL"):
                m = int(line.split()[1])
            elif m == model and line.startswith(("HETATM", "ATOM")) and line[21] == chain:
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                pos.append(int(line.split()[-1]))
            elif m > model:
                break
    return np.array(pos, np.int64), np.array(xyz)


def binned(pos, xyz, res=1_000_000):
    b = pos // res
    keys = np.unique(b)
    mean = np.array([xyz[b == k].mean(axis=0) for k in keys])
    d = np.sqrt(((mean[:, None] - mean[None]) ** 2).sum(-1))
    i, j = np.triu_indices(len(keys), 1)
    return {(int(keys[p]), int(keys[q])): float(d[p, q]) for p, q in zip(i, j)}


def rank(x):
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x))
    xs = x[order]
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and xs[j + 1] == xs[i]:
            j += 1
        r[order[i:j + 1]] = (i + j) / 2.0 + 1
        i = j + 1
    return r


def spearman(a: dict, b: dict) -> tuple:
    keys = sorted(set(a) & set(b))
    x = np.array([a[k] for k in keys])
    y = np.array([b[k] for k in keys])
    return float(np.corrcoef(rank(x), rank(y))[0, 1]), len(keys)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cell-dir", type=Path, required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--seeds", default="0,1")
    ap.add_argument("--extra", default="", help="more nuc_dynamics arguments, e.g. '-temps 100' (a shorter anneal)")
    a = ap.parse_args()
    pairs = next(a.cell_dir.glob("*contact_pairs.txt.gz"))
    pdb = next(a.cell_dir.glob("*.pdb.gz"))
    ref = binned(*read_pdb(pdb, 1))
    work = Path(tempfile.mkdtemp())
    ncc = work / "cell.ncc.gz"
    n = write_ncc(pairs, ncc)
    head = subprocess.run(["git", "-C", str(ND), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    rows = []
    for seed in map(int, a.seeds.split(",")):
        o = work / f"model_s{seed}"
        cmd = [sys.executable, str(ND / "nuc_dynamics.py"), str(ncc), "-m", "1", "-s", "8", "4", "2", "1",
               "-o", str(o), "-ran", str(seed), "-f", "n3d"] + a.extra.split()
        t0 = time.perf_counter()
        p = subprocess.run(cmd, cwd=str(ND), capture_output=True, text=True,
                           env={**os.environ, "PYTHONPATH": str(ND)})
        secs = time.perf_counter() - t0
        outs = sorted(work.glob(f"model_s{seed}*.n3d"))
        row = {"seed": seed, "seconds": secs, "exit": p.returncode, "outputs": [x.name for x in outs]}
        if p.returncode != 0 or not outs:
            row["stderr"] = p.stderr[-3000:]
            row["stdout"] = p.stdout[-2000:]
        else:
            ours = binned(*read_n3d(outs[-1]))
            row["rho"], row["n_pairs"] = spearman(ours, ref)
            row["log_tail"] = p.stdout[-800:]
        rows.append(row)
        print(json.dumps(row)[:1500], flush=True)
    res = {"label": a.label, "extra_args": a.extra, "nuc_dynamics_commit": head, "n_contacts": n, "runs": rows,
           "rho": float(np.mean([r["rho"] for r in rows if "rho" in r])) if any("rho" in r for r in rows) else None,
           "cpus": os.cpu_count(), "slurm_cpus": os.environ.get("SLURM_CPUS_PER_TASK")}
    out = json.loads(a.out.read_text()) if a.out.exists() else {}
    out[a.label] = res
    a.out.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
