"""The data folders of the schic family of version 3 (data_v3/), built from registered sources only.

    python schic_build_data.py --shared /shared --data /data --out /data_v3

``--shared`` is the shared data directory of Sherlock (``$SCRATCH/uchrom-nd/data``: the files of ``fetch``),
``--data`` the shared data root (``agentbench/data``).  Writes:

* ``stevens2017_cell2/``: the contact pairs and the published structure of Stevens et al. 2017 Cell 2 (fetch
  ``stevens2017``, GSM2219498), copied unchanged (question R2).
* ``mm10_cpg/mm10_cpg_1Mb.tsv``: the CpG ratio of every 1-Mb bin of the scHiCAR per-cell maps
  (``schicar_mop/contacts_cells_1Mb.scool``, shared data root), as scHiCluster computes it (``bedtools nuc -pattern CG -C``:
  CG occurrences, any case, within the bin over the bin's length minus its N bases; 0 where the bin has no base),
  from the UCSC mm10 FASTA (fetch ``mm10_genome``) (question R5).
* ``nagano2017_hap_2i/``: the 395 haploid mES cells of Nagano et al. 2017 grown in 2i (condition ``2i_all``) that
  passed the authors' QC (``passed_qc == 1`` in GEO's feature table): their fend-pair tables as the authors
  distribute them (``adj/<cell>.adj``: ``fend1 fend2 count``, from fetch ``nagano2017_hap``), the mm9 GATC fragment
  ends (``GATC.fends``: ``fend chr coord``) and ``cells.tsv`` (cell, condition, contacts in the adj table);
  the authors' cell-cycle features are left out (question R6).
* ``schicar_cell_types/schicar_cell_types.tsv`` (``--only types``): the cell id (= the cell name in
  ``schicar_mop/contacts_cells_1Mb.scool``) and ``cell_type`` of the 5,313 scHiCAR cells, from the cell table of the
  shared store (question R4; a plain table so that the question does not depend on reading the store).
* ``nagano2017_dip_chr19/`` (``--only snaphic``): 50 of the 742 diploid serum mES cells SnapHiC used (Nagano et al.
  2017; loader ``nagano2017_dip_serum``: the authors' fend pairs lifted to mm10 and binned at 10 kb), drawn with
  ``numpy.random.default_rng(0).permutation`` of the sorted ids of ``cells.snaphic_top742``: per cell
  ``contacts/<cell>.contacts.txt`` (SnapHiC's input: ``cell chr19 pos chr19 pos``, one line per non-empty 10-kb
  pixel of chr19, positions at the bin centres) and ``cells.txt`` (question R7).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def stevens_cell2(shared: Path, out: Path) -> dict:
    d = out / "stevens2017_cell2"
    d.mkdir(parents=True, exist_ok=True)
    files = {}
    for name in ("GSM2219498_Cell_2_contact_pairs.txt.gz", "GSM2219498_Cell_2_genome_structure_model.pdb.gz"):
        shutil.copyfile(shared / "stevens2017_mesc" / name, d / name)
        files[name] = md5(d / name)
    return files


def _fasta(path: Path):
    """(chrom, sequence bytes) of a gzipped FASTA, one chromosome at a time."""
    name, parts = None, []
    with gzip.open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(parts)
                name, parts = line[1:].split()[0].decode(), []
            else:
                parts.append(line.rstrip())
    if name is not None:
        yield name, b"".join(parts)


def mm10_cpg(shared: Path, data: Path, out: Path) -> dict:
    import cooler

    scool = data / "schicar_mop" / "contacts_cells_1Mb.scool"
    first = cooler.fileops.list_scool_cells(str(scool))[0]
    bins = cooler.Cooler(f"{scool}::{first}").bins()[["chrom", "start", "end"]][:]
    bins["chrom"] = bins["chrom"].astype(str)
    want = set(bins["chrom"])
    ratio = np.full(len(bins), np.nan)
    seen = set()
    for chrom, seq in _fasta(shared / "genomes" / "mm10" / "mm10.fa.gz"):
        if chrom not in want:
            continue
        seen.add(chrom)
        up = seq.upper()
        idx = np.flatnonzero(bins["chrom"].to_numpy() == chrom)
        for i in idx:
            s, e = int(bins.at[i, "start"]), int(bins.at[i, "end"])
            w = up[s:e]
            n = len(w) - w.count(b"N")
            ratio[i] = w.count(b"CG") / n if n > 0 else np.nan
    missing = sorted(want - seen)
    if missing:
        raise SystemExit(f"chromosomes of the maps not in the FASTA: {missing}")
    bins["cpg_ratio"] = np.nan_to_num(ratio, nan=0.0)          # scHiCluster: fillna(0)
    d = out / "mm10_cpg"
    d.mkdir(parents=True, exist_ok=True)
    path = d / "mm10_cpg_1Mb.tsv"
    bins.to_csv(path, sep="\t", index=False, float_format="%.8g")
    return {"bins": int(len(bins)), "zero_bins": int((bins["cpg_ratio"] == 0).sum()), "md5": md5(path),
            "chroms": list(dict.fromkeys(bins["chrom"]))}


def nagano_2i(shared: Path, out: Path) -> dict:
    src = shared / "nagano2017_hap"
    feat = pd.read_csv(src / "GSE94489_haploids_features_table.txt.gz", sep="\t")
    sel = feat[(feat["cond"] == "2i_all") & (feat["passed_qc"] == 1)].copy()
    sel["cell"] = sel["cell_nm"].str.replace("_", ".", regex=False)
    d = out / "nagano2017_hap_2i"
    (d / "adj").mkdir(parents=True, exist_ok=True)
    rows, missing = [], []
    for cell, cond, total in zip(sel["cell"], sel["cond"], sel["total_contacts"]):
        p = src / "adj_2i" / cell / "adj"
        if not p.exists():
            missing.append(cell)
            continue
        shutil.copyfile(p, d / "adj" / f"{cell}.adj")
        n = sum(1 for _ in open(p)) - 1
        rows.append({"cell": cell, "condition": cond, "n_fend_pairs": n, "total_contacts_geo": int(total)})
    if missing:
        raise SystemExit(f"{len(missing)} selected cells without an adj file, e.g. {missing[:5]}")
    shutil.copyfile(src / "GATC.fends", d / "GATC.fends")
    cells = pd.DataFrame(rows)
    cells[["cell", "condition", "n_fend_pairs"]].to_csv(d / "cells.tsv", sep="\t", index=False)
    diff = (cells["n_fend_pairs"] - cells["total_contacts_geo"]).abs()
    return {"cells": int(len(cells)), "fends_md5": md5(d / "GATC.fends"),
            "n_fend_pairs_median": float(cells["n_fend_pairs"].median()),
            "n_fend_pairs_range": [int(cells["n_fend_pairs"].min()), int(cells["n_fend_pairs"].max())],
            "adj_rows_minus_geo_total_contacts_max_abs": int(diff.max()),
            "cells_md5": md5(d / "cells.tsv")}


def nagano_dip_chr19(shared: Path, out: Path, n: int = 50, seed: int = 0) -> dict:
    import cooler
    from chromdata import ChromData

    src = shared / "nagano2017_dip_serum"
    cd = ChromData.read(str(src / "nagano2017_dip_serum.chromdata.zarr"))
    top = sorted(map(str, cd.cells.index[cd.cells["snaphic_top742"].to_numpy(bool)]))
    perm = np.random.default_rng(seed).permutation(len(top))
    cells = [top[k] for k in perm[:n]]
    scool = src / "nagano2017_dip_serum_10kb.scool"
    d = out / "nagano2017_dip_chr19"
    (d / "contacts").mkdir(parents=True, exist_ok=True)
    lines = []
    for c in cells:
        px = cooler.Cooler(f"{scool}::/cells/{c}").matrix(balance=False, as_pixels=True, join=True).fetch("chr19")
        t = pd.DataFrame({"cell": c, "c1": "chr19", "p1": px["start1"].to_numpy() + 5000, "c2": "chr19",
                          "p2": px["start2"].to_numpy() + 5000})
        t.to_csv(d / "contacts" / f"{c}.contacts.txt", sep="\t", header=False, index=False)
        lines.append(len(t))
    (d / "cells.txt").write_text("\n".join(cells) + "\n")
    return {"cells": len(cells), "pixels_median": float(np.median(lines)), "pixels_range": [min(lines), max(lines)],
            "seed": seed, "first": cells[:3]}


def schicar_types(data: Path, out: Path) -> dict:
    from chromdata import ChromData

    cd = ChromData.read(str(data / "schicar_mop" / "schicar_mop.chromdata.zarr"))
    t = pd.DataFrame({"cell": cd.cells.index.astype(str), "cell_type": cd.cells["cell_type"].astype(str).to_numpy()})
    d = out / "schicar_cell_types"
    d.mkdir(parents=True, exist_ok=True)
    t.to_csv(d / "schicar_cell_types.tsv", sep="\t", index=False)
    return {"cells": int(len(t)), "types": int(t["cell_type"].nunique()), "md5": md5(d / "schicar_cell_types.tsv")}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--shared", type=Path, required=True)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--only", default="stevens,cpg,nagano")
    a = ap.parse_args()
    only = set(a.only.split(","))
    log = {}
    if "stevens" in only:
        log["stevens2017_cell2"] = stevens_cell2(a.shared, a.out)
        print("stevens2017_cell2", log["stevens2017_cell2"], flush=True)
    if "nagano" in only:
        log["nagano2017_hap_2i"] = nagano_2i(a.shared, a.out)
        print("nagano2017_hap_2i", log["nagano2017_hap_2i"], flush=True)
    if "types" in only:
        log["schicar_cell_types"] = schicar_types(a.data, a.out)
        print("schicar_cell_types", log["schicar_cell_types"], flush=True)
    if "snaphic" in only:
        log["nagano2017_dip_chr19"] = nagano_dip_chr19(a.shared, a.out)
        print("nagano2017_dip_chr19", log["nagano2017_dip_chr19"], flush=True)
    if "cpg" in only:
        log["mm10_cpg"] = mm10_cpg(a.shared, a.data, a.out)
        print("mm10_cpg", log["mm10_cpg"], flush=True)
    logf = a.out / "schic_build_log.json"
    old = json.loads(logf.read_text()) if logf.exists() else {}
    logf.write_text(json.dumps({**old, **log}, indent=1))


if __name__ == "__main__":
    main()
