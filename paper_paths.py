"""Where this repository finds its own parts, the U-Chrom code, the data and the paper.

* ``FIGURES / <item>``: one folder per figure (``fig1``, ``extended_data_fig1``, ``supp_fig1``, ...), with its script,
  inputs and rendered figure; ``TABLES / <item>`` likewise for generated tables;
* ``results(<topic>)``: ``benchmarks/<topic>/results``, the numbers a benchmark writes and the figures read
  (topics: store, browser, agents, atlas, cross_assay, spatial_variable, toolkit_parity/...);
* ``SHARED``: drawing code and illustrations several figures use;
* ``UCHROM_REPO``: the U-Chrom checkout (default: a sibling folder ``U-Chrom``): its packages when the installed ones
  should not be used, ``skills/``, ``docs/``, test fixtures, ``apps/atlas``;
* ``DATA``: the data directory of ``uchrom.datasets`` (``$UCHROM_DATA``, default ``~/.cache/uchrom``);
* ``PAPER_CHECKOUT``: the paper repository that ``sync_to_paper.py`` writes into (``$UCHROM_PAPER``, default a sibling
  folder ``u-chrom-paper``).

Scripts import it after putting the repository root on ``sys.path``::

    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # figures/<item>/x.py -> the root
    from paper_paths import SHARED, results
"""

from __future__ import annotations

import os
from pathlib import Path

REPRO = Path(__file__).resolve().parent            # this repository
PAPER = REPRO                                       # former name: the figure folders' root (same layout as the paper)
UCHROM_REPO = Path(os.environ.get("UCHROM_REPO", PAPER.parent / "U-Chrom")).expanduser().resolve()
DATA = (Path(os.environ["UCHROM_DATA"]) if os.environ.get("UCHROM_DATA") else Path.home() / ".cache" / "uchrom").expanduser().resolve()
EXAMPLE_DATA = DATA                                  # former name
FIGURES = REPRO / "figures"
TABLES = REPRO / "tables"
SHARED = REPRO / "shared"
BENCHMARKS = REPRO / "benchmarks"


def results(topic: str) -> Path:
    """``benchmarks/<topic>/results``: what a benchmark writes and the figures read."""
    return BENCHMARKS / topic / "results"


PAPER_CHECKOUT = Path(os.environ.get("UCHROM_PAPER", REPRO.parent / "u-chrom-paper")).expanduser().resolve()
