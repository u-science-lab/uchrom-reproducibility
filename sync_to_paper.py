"""Copy the rendered figures and generated tables into the paper checkout.

The paper repository keeps only what its LaTeX reads; this repository makes those files, one folder per item.
``FILES`` maps each of them to its place in the paper checkout (``MANUSCRIPT``: the manuscript's folder there).

    python sync_to_paper.py [--paper ../u-chrom-paper] [--dry-run]
"""

from __future__ import annotations

import argparse
import filecmp
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paper_paths import PAPER_CHECKOUT  # noqa: E402

MANUSCRIPT = "manuscripts/nmeth-bc"          # the manuscript's folder in the paper checkout
FIG = f"{MANUSCRIPT}/figures"
FILES = {
    "figures/fig1/fig1.pdf": f"{FIG}/fig1.pdf",
    "figures/fig2/fig2.pdf": f"{FIG}/fig2.pdf",
    "figures/extended_data_fig1/ed_store.pdf": f"{FIG}/ed_store.pdf",
    "figures/extended_data_fig2/ed_melanoma.pdf": f"{FIG}/ed_melanoma.pdf",
    "figures/extended_data_fig3/ed_brain.pdf": f"{FIG}/ed_brain.pdf",
    "figures/supp_fig1/supp_access.pdf": f"{FIG}/supp_access.pdf",
    "figures/supp_fig2/supp_linked.pdf": f"{FIG}/supp_linked.pdf",
    **{f"figures/{item}/supp_browser_{name}_framed.{ext}": f"{FIG}/captures/supp_browser_{name}_framed.{ext}"
       for item, name in (("supp_fig3", "spatial"), ("supp_fig4", "scrna"), ("supp_fig5", "tracing"))
       for ext in ("png", "tex")},
    "figures/supp_fig6/supp_atlas.pdf": f"{FIG}/supp_atlas.pdf",
    "figures/supp_fig7/supp_agents.pdf": f"{FIG}/supp_agents.pdf",
    "tables/supp_table2/table_atlas.tex": f"{MANUSCRIPT}/supp/table_atlas.tex",
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paper", type=Path, default=PAPER_CHECKOUT, help="the paper checkout ($UCHROM_PAPER)")
    ap.add_argument("--dry-run", action="store_true", help="list what would change")
    a = ap.parse_args()
    paper = a.paper.expanduser().resolve()
    if not (paper / MANUSCRIPT).is_dir():
        raise SystemExit(f"{paper}: no {MANUSCRIPT} (set UCHROM_PAPER or pass --paper)")
    changed = []
    for src_rel, dst_rel in FILES.items():
        src, dst = HERE / src_rel, paper / dst_rel
        if not src.exists():
            print(f"not made yet: {src_rel}", file=sys.stderr)
        elif not dst.exists() or not filecmp.cmp(src, dst, shallow=False):
            changed.append(dst_rel)
            if not a.dry_run:
                shutil.copy2(src, dst)
    for rel in changed:
        print(("would copy " if a.dry_run else "copied ") + rel)
    print(f"{len(changed)} of {len(FILES)} file(s) {'to copy' if a.dry_run else 'copied'} into {paper}")


if __name__ == "__main__":
    main()
