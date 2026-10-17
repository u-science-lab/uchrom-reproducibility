"""Copy the rendered figures and generated tables into the paper checkout.

The paper repository (``u-science-lab/u-chrom-paper``) keeps only what its LaTeX reads: the figure files under
``manuscripts/<ms>/figures/`` and the generated tables (``manuscripts/nmeth-bc/supp/table_atlas.tex``,
``manuscripts/full/supp/tableS1_*.tex``).  This repository makes them, at the same paths.  The files to copy are the
ones the paper checkout tracks there, so a new figure is added to the paper once (``git add`` there) and kept in
sync from then on.

    python sync_to_paper.py [--paper ../u-chrom-paper] [--dry-run]
"""

from __future__ import annotations

import argparse
import filecmp
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paper_paths import PAPER_CHECKOUT  # noqa: E402

PATTERNS = [":(glob)manuscripts/*/figures/**", "manuscripts/nmeth-bc/supp/table_atlas.tex", "manuscripts/full/supp/tableS1_*.tex"]


def tracked(paper: Path) -> list[str]:
    out = subprocess.run(["git", "-C", str(paper), "ls-files", "--", *PATTERNS], capture_output=True, text=True,
                         check=True).stdout
    return [line for line in out.splitlines() if line]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paper", type=Path, default=PAPER_CHECKOUT, help="the u-chrom-paper checkout ($UCHROM_PAPER)")
    ap.add_argument("--dry-run", action="store_true", help="list what would change")
    a = ap.parse_args()
    paper = a.paper.expanduser().resolve()
    if not (paper / "manuscripts").is_dir():
        raise SystemExit(f"{paper}: not a u-chrom-paper checkout (set UCHROM_PAPER or pass --paper)")
    changed, missing = [], []
    for rel in tracked(paper):
        src, dst = HERE / rel, paper / rel
        if not src.exists():
            missing.append(rel)
        elif not dst.exists() or not filecmp.cmp(src, dst, shallow=False):
            changed.append(rel)
            if not a.dry_run:
                shutil.copy2(src, dst)
    for rel in changed:
        print(("would copy " if a.dry_run else "copied ") + rel)
    for rel in missing:
        print(f"not made here: {rel}", file=sys.stderr)
    print(f"{len(changed)} file(s) {'to copy' if a.dry_run else 'copied'} into {paper}")


if __name__ == "__main__":
    main()
