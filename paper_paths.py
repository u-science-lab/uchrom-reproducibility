"""Where this repository finds the U-Chrom code, the data and the paper.

This is the reproducibility repository of the U-Chrom papers (``u-science-lab/uchrom-reproducibility``): the
benchmarks that make the papers' numbers (``benchmarks/``), their results (``manuscripts/<ms>/figures/<fig>/results``)
and the figure scripts (``manuscripts/<ms>/figures``), in the layout of the paper repository, whose LaTeX reads the
rendered figures (``sync_to_paper.py`` copies them there).

* ``UCHROM_REPO``: the U-Chrom checkout (default: a sibling folder ``U-Chrom`` next to this repository): its
  packages when the installed ones should not be used, ``skills/``, ``docs/``, test fixtures, ``apps/atlas``;
* ``UCHROM_DATA``: the data directory of ``uchrom.datasets`` (default ``~/.cache/uchrom``);
* ``UCHROM_PAPER``: the u-chrom-paper checkout that ``sync_to_paper.py`` writes into (default: a sibling folder
  ``u-chrom-paper``).

Scripts import it after putting the repository root on ``sys.path``::

    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))   # manuscripts/full/figures/<fig>/x.py -> root
    from paper_paths import UCHROM_REPO, DATA
"""

from __future__ import annotations

import os
from pathlib import Path

REPRO = Path(__file__).resolve().parent            # this repository
PAPER = REPRO                                       # former name: the figure folders' root (same layout as the paper)
UCHROM_REPO = Path(os.environ.get("UCHROM_REPO", PAPER.parent / "U-Chrom")).expanduser().resolve()
DATA = (Path(os.environ["UCHROM_DATA"]) if os.environ.get("UCHROM_DATA") else Path.home() / ".cache" / "uchrom").expanduser().resolve()
EXAMPLE_DATA = DATA                                  # former name
FIGURES = REPRO / "manuscripts" / "full" / "figures"   # the full manuscript's figure folders
PAPER_CHECKOUT = Path(os.environ.get("UCHROM_PAPER", REPRO.parent / "u-chrom-paper")).expanduser().resolve()
