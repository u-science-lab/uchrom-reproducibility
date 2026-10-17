"""Where things are: this repository's folders, the U-Chrom code and the data.

One folder per figure or table whose numbers come from data; a benchmark that several items use sits in the folder
of the first one.  ``folder(topic)`` finds a benchmark by its topic, ``results(topic)`` its results:

    cross_assay       fig2/                  store            extended_data_fig1/
    spatial_variable  extended_data_fig2/    browser          supp_fig6/browser/
    atlas             supp_fig6/atlas/       agents           supp_fig7/
    toolkit_parity/<check>  supp_table3/<check>/

* ``REPRO``: this repository, or ``$UCHROM_REPRO`` (another checkout of it, e.g. a copy on a cluster: results go
  there);
* ``UCHROM_REPO``: the U-Chrom checkout (``$UCHROM_REPO``, default a sibling folder ``U-Chrom``): its packages when
  the installed ones should not be used, ``skills/``, ``docs/``, test fixtures, ``apps/atlas``;
* ``DATA``: the data directory of ``uchrom.datasets`` (``$UCHROM_DATA``, default ``~/.cache/uchrom``).

Scripts put this folder on ``sys.path``::

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))   # <folder>/x.py
    from paper_paths import results
"""

from __future__ import annotations

import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPRO = Path(os.environ.get("UCHROM_REPRO") or HERE.parent).expanduser().resolve()
UCHROM_REPO = Path(os.environ.get("UCHROM_REPO", HERE.parent.parent / "U-Chrom")).expanduser().resolve()
DATA = (Path(os.environ["UCHROM_DATA"]) if os.environ.get("UCHROM_DATA") else Path.home() / ".cache" / "uchrom").expanduser().resolve()

FOLDERS = {
    "cross_assay": "fig2",
    "store": "extended_data_fig1",
    "spatial_variable": "extended_data_fig2",
    "browser": "supp_fig6/browser",
    "atlas": "supp_fig6/atlas",
    "agents": "supp_fig7",
    "toolkit_parity": "supp_table3",
}


def folder(topic: str) -> Path:
    """The folder of a benchmark: ``store``, ``toolkit_parity/nucdyn``, ..."""
    head, _, rest = topic.partition("/")
    if head not in FOLDERS:
        raise KeyError(f"unknown topic {topic!r} (one of {', '.join(FOLDERS)})")
    return REPRO / FOLDERS[head] / rest if rest else REPRO / FOLDERS[head]


def results(topic: str) -> Path:
    """``<folder>/results``: what a benchmark writes and the figures read (not created here; writers create it)."""
    return folder(topic) / "results"
