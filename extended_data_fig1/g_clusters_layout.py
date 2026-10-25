"""The layout of the clusters' locus index (``ChromData.clusters``, format 2.4) on real SPRITE data: chunk size and
codec against the store's size and what a window reads.

The same table as ``g_clusters_sprite.py`` (Quinodoz et al. 2018, GSE114242, the MAPQ >= 10 clusters) is written
once per variant with ``chromdata.clusters.add_clusters(locus_chunk=, locus_compressor=)``; for each: the locus
index's size, and 40 kb / 1 Mb windows (50 each) backed from local disk and over HTTP (a local server with range
requests: requests and bytes per window), checked against the in-memory table, whose own query time is measured too.

    python extended_data_fig1/g_clusters_layout.py --archive $UCHROM_DATA/quinodoz2018_sprite/GSE114242_mouse_combined_clusters.tar.gz \\
        --work $L_SCRATCH/layout
"""

from __future__ import annotations

import argparse
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "shared"))
sys.path.insert(0, str(HERE))
from _record import save  # noqa: E402
from g_clusters_sprite import dir_bytes, parse, rss_gb, timed, windows  # noqa: E402


def variants():
    import zarr

    zstd = zarr.codecs.ZstdCodec(level=9)
    blosc = zarr.codecs.BloscCodec(cname="zstd", clevel=5, shuffle="bitshuffle")
    return {"zstd9_16k": (1 << 14, zstd), "zstd9_8k": (1 << 13, zstd), "zstd9_4k": (1 << 12, zstd),
            "bitshuffle_16k": (1 << 14, blosc), "bitshuffle_8k": (1 << 13, blosc), "bitshuffle_4k": (1 << 12, blosc)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--member", default="mouse.mapq-ge10.clusters.gz")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--n-windows", type=int, default=50)
    ap.add_argument("--only", nargs="*", help="variant names")
    ap.add_argument("--out", type=Path, default=HERE / "results" / "clusters_layout.json")
    args = ap.parse_args()

    from chromdata import ChromData, ClusterTable
    from chromdata.clusters import add_clusters
    from chromdata.testing import RangeRequestHandler, serve_directory

    sizes, codes, pos, chroms = parse(args.archive, args.member, args.limit)
    sizes = sizes[sizes > 0]
    lengths = {chroms[c]: int(v) for c, v in pd.Series(pos).groupby(codes).max().items()}
    df = pd.DataFrame({"cluster": np.repeat(np.arange(len(sizes), dtype=np.int64), sizes),
                       "chrom": pd.Categorical.from_codes(codes, categories=chroms), "start": pos, "end": pos + 1})
    del codes, pos
    t = ClusterTable.from_frame(df, keep_names=False)
    del df
    res = {"dataset": "quinodoz2018_sprite", "member": args.member, "clusters": t.n_clusters, "reads": t.n_members,
           "windows": {}, "variants": {}}
    wins = {w: windows(lengths, args.n_windows, w, seed=w + 7) for w in (40_000, 1_000_000)}
    for w, ws in wins.items():                       # the in-memory table (locus built on first use)
        for c in {x[0] for x in ws}:
            t.region(c, 0, 1)
        res["windows"][str(w)] = {"in_memory": timed(t.region, ws)}
    print(t, res["windows"], flush=True)

    empty = pd.DataFrame({"chrom": pd.Series([], dtype=str), "start": pd.Series([], dtype=np.int64),
                          "end": pd.Series([], dtype=np.int64), "trace_id": pd.Series([], dtype=str)})
    for name, (chunk, codec) in variants().items():
        if args.only and name not in args.only:
            continue
        store = args.work / f"{name}.chromdata.zarr"
        shutil.rmtree(store, ignore_errors=True)
        ChromData(np.zeros((0, 3)), empty).write(store)
        t0 = time.perf_counter()
        add_clusters(store, "sprite", t, locus_chunk=chunk, locus_compressor=codec)
        v = {"locus_chunk": chunk, "codec": type(codec).__name__, "write_s": time.perf_counter() - t0,
             "locus_bytes": dir_bytes(store / "clusters" / "sprite" / "by_locus"),
             "members_bytes": dir_bytes(store / "clusters" / "sprite" / "members"), "windows": {}}
        b = ChromData.read(store, backed=True).clusters["sprite"]
        for w, ws in wins.items():
            v["windows"][str(w)] = {"backed_local": timed(b.region, ws),
                                    "backed_equals_memory": all(b.region(*x).equals(t.region(*x)) for x in ws[:8])}
        with serve_directory(args.work) as base:
            h = ChromData.read(f"{base}/{store.name}", backed=True).clusters["sprite"]
            for w, ws in wins.items():
                for c in {x[0] for x in ws}:
                    h.region(c, 0, 1)                # each chromosome's fence once
                RangeRequestHandler.log.clear()
                ht = timed(h.region, ws)
                log = list(RangeRequestHandler.log)
                ht["requests_per_window"] = len(log) / len(ws)
                ht["kb_per_window"] = sum(n for _, _, n in log) / 1024 / len(ws)
                ht["equals_memory"] = all(h.region(*x).equals(t.region(*x)) for x in ws[:5])
                v["windows"][str(w)]["backed_http_local_server"] = ht
        res["variants"][name] = v
        print(name, {k: v[k] for k in ("locus_bytes", "write_s")},
              {w: (round(x["backed_http_local_server"]["kb_per_window"]), round(x["backed_http_local_server"]["requests_per_window"], 1),
                   round(x["backed_http_local_server"]["median_ms"], 1)) for w, x in v["windows"].items()}, flush=True)
        shutil.rmtree(store, ignore_errors=True)
    res["peak_rss_gb"] = rss_gb()
    save(args.out, res, script=__file__, datasets=["quinodoz2018_sprite"], params=vars(args))


if __name__ == "__main__":
    main()
