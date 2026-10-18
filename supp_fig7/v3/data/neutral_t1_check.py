"""T1 from the FOF-CT export instead of the store: the store keeps float32 coordinates, and the core table writes their
shortest decimals, which read back as float64 differ from the float32 values in the last bits.  This builds T1's input
(chr15, 1-Mb bins, the question's recipe) from the core table with pandas in float64 -- the read of an agent without
U-Chrom -- next to the reference's input (benchmarks/arcfish_parity/inputs.py from the store), and writes a manifest
for run_upstream.py / run_ours.py; the job then compares the domain calls with the reference.

    python neutral_t1_check.py prepare --core CORE.csv --store STORE --manifest IN/manifest.json --out DIR
    python neutral_t1_check.py compare --runs RUN_DIR --reference results/arcfish.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "truth"))
NAME = "takei2025_chr15_1mb"


def prepare(a) -> None:
    from chromdata import ChromData

    from neutral_build import read_fofct

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cd = ChromData.read(a.store)
    xyz = np.asarray(cd.coords)
    t, _ = read_fofct(Path(a.core))
    f64 = t[["X", "Y", "Z"]].to_numpy(np.float64)
    d = np.abs(f64 - xyz.astype(np.float64))
    rep = {"store_dtype": str(xyz.dtype), "max_abs_diff_um": float(np.nanmax(d)),
           "max_rel_diff": float(np.nanmax(d / np.maximum(np.abs(xyz.astype(np.float64)), 1e-12))),
           "equal_as_float32": bool(np.array_equal(f64.astype(np.float32), xyz.astype(np.float32), equal_nan=True))}

    # the question's recipe on the core table (pandas, float64)
    on = t["Chrom"].astype(str) == "chr15"
    df = pd.DataFrame({"trace_id": t.loc[on, "Trace_ID"].astype(str).to_numpy(),
                       "b": (t.loc[on, "Chrom_Start"].to_numpy(np.int64) + t.loc[on, "Chrom_End"].to_numpy(np.int64))
                       // 2 // 1_000_000})
    df[["x", "y", "z"]] = f64[on.to_numpy()]
    g = df.groupby(["trace_id", "b"], sort=True).agg(x=("x", "mean"), y=("y", "mean"), z=("z", "mean")).reset_index()
    spots = pd.DataFrame({"chrom": "chr15", "start": (g["b"] * 1_000_000).astype(np.int64),
                          "end": ((g["b"] + 1) * 1_000_000).astype(np.int64), "trace_id": g["trace_id"]})
    b = ChromData(coords=g[["x", "y", "z"]].to_numpy(float), spots=spots)
    b.uns["xyz_unit"] = cd.uns.get("xyz_unit", "micron")
    path = out / f"{NAME}.csv"
    b.to_fofct(path)

    # the reference's input from the store, for the size of the difference
    man = json.load(open(a.manifest))
    ref_t, _ = read_fofct(Path(man[NAME]["files"]))
    key = ["Trace_ID", "Chrom_Start"]
    new_t, _ = read_fofct(path)
    m = ref_t.assign(Trace_ID=ref_t["Trace_ID"].astype(str)).merge(
        new_t.assign(Trace_ID=new_t["Trace_ID"].astype(str)), on=key, suffixes=("_ref", "_csv"))
    rep["binned_rows"] = {"reference": int(len(ref_t)), "from_csv": int(len(g)), "matched": int(len(m))}
    rep["binned_max_abs_diff_um"] = float(np.nanmax(np.abs(m[["X_ref", "Y_ref", "Z_ref"]].to_numpy()
                                                           - m[["X_csv", "Y_csv", "Z_csv"]].to_numpy())))
    man = {NAME: {**man[NAME], "files": str(path)}}
    (out / "manifest.json").write_text(json.dumps(man, indent=1))
    (out / "prepare.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1), flush=True)


def compare(a) -> None:
    from arcfish import t1

    ours, theirs = t1([Path(a.runs)])
    ref = json.load(open(a.reference))["T1"]["values"]
    res = {"from_csv": {"uchrom": ours, "arcfish_a21d889": theirs}, "reference": ref,
           "uchrom_equals_reference": ours == ref, "authors_equal_reference": theirs == ref}
    Path(a.runs, "t1_from_csv.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1), flush=True)
    sys.exit(0 if ours == ref and theirs == ref else 1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=["prepare", "compare"])
    ap.add_argument("--core")
    ap.add_argument("--store")
    ap.add_argument("--manifest")
    ap.add_argument("--out")
    ap.add_argument("--runs")
    ap.add_argument("--reference")
    a = ap.parse_args()
    (prepare if a.what == "prepare" else compare)(a)


if __name__ == "__main__":
    main()
