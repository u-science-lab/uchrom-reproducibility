"""Multi-way contacts in the store (``ChromData.clusters``, format 2.4) on real SPRITE data.

Quinodoz et al. 2018 (Cell 174:744), GEO GSE114242: the mouse ES cell clusters of the authors' pipeline, one
cluster per line (``barcode<TAB>chrom:pos<TAB>chrom:pos ...``, one token per read).  A read position ``pos``
(0-based: the pipeline writes pysam's ``reference_start``; positions 0 occur) becomes the fragment ``[pos, pos + 1)``.

Measured:

* **round trip** — every cluster of the file, as a multiset of (chrom, pos), against the store read back: an
  order-independent 64-bit checksum over clusters (each cluster hashed from its reads), cluster / read counts and
  the size distribution;
* **size** — the store's parts (members, offsets, locus index) against the gzipped text;
* **time / memory** — parsing, ``ClusterTable.from_frame``, writing, opening;
* **region queries** — random windows of 40 kb, 1 Mb and 10 Mb: in memory, backed from the local store, and backed
  over HTTP (a local server with range requests: requests and bytes per window), the backed results checked
  against the in-memory ones.

    python extended_data_fig1/g_clusters_sprite.py --archive $UCHROM_DATA/quinodoz2018_sprite/GSE114242_mouse_combined_clusters.tar.gz \\
        --member mouse.mapq-ge10.clusters.gz --work $L_SCRATCH/clusters
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import resource
import shutil
import sys
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "shared"))
from _record import save  # noqa: E402

M64 = np.uint64(0xFFFFFFFFFFFFFFFF)


def _mix(x: np.ndarray) -> np.ndarray:
    """splitmix64 finaliser (uint64 in, uint64 out)."""
    with np.errstate(over="ignore"):
        z = x.astype(np.uint64) + np.uint64(0x9E3779B97F4A7C15)
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return z ^ (z >> np.uint64(31))


def cluster_checksum(cid: np.ndarray, pos: np.ndarray, sizes: np.ndarray) -> tuple:
    """(sum, xor) over clusters of hash(cluster); a cluster's hash = mix(size + sum of mix(chrom, pos)) — independent
    of the order of the clusters and of the reads in a cluster.  ``cid``: the reads' chromosome ids in one fixed
    name -> id table, so both sides hash the same values; every cluster has at least one read."""
    cid = cid.astype(np.uint64)
    with np.errstate(over="ignore"):
        h = _mix((cid << np.uint64(40)) | pos.astype(np.uint64))
        starts = np.r_[0, np.cumsum(sizes)[:-1]].astype(np.int64)
        per = np.add.reduceat(h, starts) if len(h) else np.zeros(0, np.uint64)
        per = _mix(per + sizes.astype(np.uint64))
        return int(per.sum(dtype=np.uint64)), int(np.bitwise_xor.reduce(per)) if len(per) else 0


def parse(archive: Path, member: str, limit: int | None, chunk: int = 1 << 20):
    """The member's clusters: (cluster sizes, chrom names (categorical codes + names), positions); a checksum is
    taken from the text as it is parsed."""
    sizes, codes, poss = [], [], []
    names: dict = {}
    n_lines = 0
    with tarfile.open(archive, "r|gz") as tar:
        for m in tar:
            if Path(m.name).name != member:
                continue
            with gzip.open(tar.extractfile(m), "rt") as fh:
                buf = []
                for line in fh:
                    buf.append(line)
                    if len(buf) == chunk or (limit and n_lines + len(buf) >= limit):
                        n_lines += _chunk(buf, sizes, codes, poss, names)
                        buf = []
                        if limit and n_lines >= limit:
                            break
                if buf:
                    n_lines += _chunk(buf, sizes, codes, poss, names)
            break
        else:
            raise SystemExit(f"{member} not in {archive}")
    chroms = [c for c, _ in sorted(names.items(), key=lambda kv: kv[1])]
    return (np.concatenate(sizes), np.concatenate(codes), np.concatenate(poss), chroms)


def _chunk(lines, sizes, codes, poss, names) -> int:
    toks, n = [], np.empty(len(lines), dtype=np.int64)
    for i, line in enumerate(lines):
        reads = [x for x in line.rstrip("\n").split("\t")[1:] if x]    # a trailing tab is no read
        n[i] = len(reads)
        toks.extend(reads)
    s = pd.Series(toks, dtype=object).str.rpartition(":")
    ch = s[0].to_numpy()
    for c in pd.unique(ch):
        names.setdefault(c, len(names))
    codes.append(np.asarray(pd.Series(ch).map(names), dtype=np.int16))
    poss.append(s[2].astype(np.int64).to_numpy())
    sizes.append(n)
    return len(lines)


def rss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1 << 30) if sys.platform == "darwin" else r / (1 << 20)


def dir_bytes(p: Path) -> int:
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file())


def windows(lengths: dict, n: int, width: int, seed: int):
    """Random windows on the main chromosomes (inside the reads' extent)."""
    rng = np.random.default_rng(seed)
    big = [c for c, n_ in lengths.items() if c not in ("chrM", "chrY") and "_" not in c and n_ > width + 4_000_000]
    out = []
    for _ in range(n):
        c = big[rng.integers(len(big))]
        s = int(rng.integers(3_000_000, lengths[c] - width))
        out.append((c, s, s + width))
    return out


def timed(f, wins):
    ts, rows = [], 0
    for w in wins:
        t0 = time.perf_counter()
        r = f(*w)
        ts.append(time.perf_counter() - t0)
        rows += len(r)
    return {"median_ms": 1e3 * float(np.median(ts)), "p90_ms": 1e3 * float(np.quantile(ts, 0.9)),
            "mean_fragments": rows / len(wins)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--member", default="mouse.mapq-ge10.clusters.gz")
    ap.add_argument("--limit", type=int, default=None, help="first N clusters only")
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--n-windows", type=int, default=50)
    ap.add_argument("--out", type=Path, default=HERE / "results" / "clusters_sprite.json")
    args = ap.parse_args()

    from chromdata import ChromData, ClusterTable
    from chromdata.testing import RangeRequestHandler, serve_directory

    res: dict = {"dataset": "quinodoz2018_sprite", "member": args.member}
    t0 = time.perf_counter()
    sizes, codes, pos, chroms = parse(args.archive, args.member, args.limit)
    res["parse_s"] = time.perf_counter() - t0
    ids = {c: i for i, c in enumerate(chroms)}
    n_empty = int((sizes == 0).sum())                # barcodes without reads: not clusters of the table
    sizes = sizes[sizes > 0]
    src_sum = cluster_checksum(codes, pos, sizes)
    lengths = {chroms[c]: int(v) for c, v in pd.Series(pos).groupby(codes).max().items()}
    res["file"] = {"clusters": int(len(sizes)), "lines_without_reads": n_empty, "reads": int(sizes.sum()),
                   "chromosomes": len(chroms),
                   "size_hist": {str(k): int(v) for k, v in zip(*np.unique(np.minimum(sizes, 1001), return_counts=True))
                                 if k in (0, 1, 2, 3, 4, 5, 10, 100, 1000, 1001)},
                   "checksum": list(src_sum)}
    res["file"]["size_quantiles"] = {str(q): float(np.quantile(sizes, q)) for q in (0.5, 0.9, 0.99, 0.999)}
    print(json.dumps(res["file"])[:400], flush=True)

    t0 = time.perf_counter()
    df = pd.DataFrame({"cluster": np.repeat(np.arange(len(sizes), dtype=np.int64), sizes),
                       "chrom": pd.Categorical.from_codes(codes, categories=chroms),
                       "start": pos, "end": pos + 1})
    del codes
    t = ClusterTable.from_frame(df, keep_names=False, attrs={
        "assay": "SPRITE", "source": "GEO GSE114242 " + args.member,
        "fragments": "read positions as written (0-based read starts) -> [pos, pos + 1)"})
    del df
    res["from_frame_s"] = time.perf_counter() - t0
    res["peak_rss_gb_after_build"] = rss_gb()
    print(t, res["from_frame_s"], flush=True)

    args.work.mkdir(parents=True, exist_ok=True)
    store = args.work / "sprite.chromdata.zarr"
    shutil.rmtree(store, ignore_errors=True)
    empty = pd.DataFrame({"chrom": pd.Series([], dtype=str), "start": pd.Series([], dtype=np.int64),
                          "end": pd.Series([], dtype=np.int64), "trace_id": pd.Series([], dtype=str)})
    cd = ChromData(np.zeros((0, 3)), empty, uns={"genome_assembly": "mm9"})
    cd.clusters["sprite"] = t
    t0 = time.perf_counter()
    cd.write(store)
    res["write_s"] = time.perf_counter() - t0
    parts = {p.name: dir_bytes(p) if p.is_dir() else p.stat().st_size
             for p in (store / "clusters" / "sprite").iterdir()}
    res["store_bytes"] = {"total": dir_bytes(store), **parts}
    res["archive_member_bytes"] = None
    with tarfile.open(args.archive, "r:gz") as tar:
        for m in tar.getmembers():
            if Path(m.name).name == args.member:
                res["archive_member_bytes"] = int(m.size)
    res["bytes_per_read"] = res["store_bytes"]["total"] / max(1, res["file"]["reads"])

    # round trip: the store read back (backed, chunk by chunk) against the file
    t0 = time.perf_counter()
    b = ChromData.read(store, backed=True).clusters["sprite"]
    res["open_backed_s"] = time.perf_counter() - t0
    off = b.member_offsets
    mem = b.members_array.astype(np.int64)
    code_map = np.array([ids[c] for c in b.chroms], dtype=np.uint64)
    back_sum = cluster_checksum(code_map[mem[0]], mem[1], np.diff(off))     # start == pos
    res["round_trip"] = {"checksum_store": list(back_sum), "equal": bool(back_sum == src_sum),
                         "clusters_equal": int(b.n_clusters) == res["file"]["clusters"],
                         "reads_equal": int(b.n_members) == res["file"]["reads"],
                         "one_bp_fragments": bool((mem[2] - mem[1] == 1).all())}
    del mem
    print(res["round_trip"], flush=True)

    # region queries
    q = {}
    for width in (40_000, 1_000_000, 10_000_000):
        wins = windows(lengths, args.n_windows, width, seed=width)
        mem_t = timed(t.region, wins)
        back_t = timed(b.region, wins)
        same = all(b.region(*w).equals(t.region(*w)) for w in wins[:10])
        q[str(width)] = {"in_memory": mem_t, "backed_local": back_t, "backed_equals_memory": same}
    with serve_directory(args.work) as base:
        h = ChromData.read(f"{base}/{store.name}", backed=True).clusters["sprite"]
        for width in (40_000, 1_000_000, 10_000_000):
            wins = windows(lengths, args.n_windows, width, seed=width + 1)
            h.region(*wins[0])            # the chromosome's fence once
            RangeRequestHandler.log.clear()
            ht = timed(h.region, wins[1:])
            log = list(RangeRequestHandler.log)
            ht["requests_per_window"] = len(log) / (len(wins) - 1)
            ht["kb_per_window"] = sum(n for _, _, n in log) / 1024 / (len(wins) - 1)
            ht["only_locus_index"] = all("/by_locus/" in p for p, _, _ in log)
            ht["equals_memory"] = all(h.region(*w).equals(t.region(*w)) for w in wins[1:6])
            q[str(width)]["backed_http_local_server"] = ht
    res["region"] = q
    res["peak_rss_gb"] = rss_gb()
    save(args.out, res, script=__file__, datasets=["quinodoz2018_sprite"], params=vars(args))
    print(json.dumps(res, indent=1, default=str)[:4000])


if __name__ == "__main__":
    main()
