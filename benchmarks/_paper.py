"""Where the benchmarks write their results: this repository.

A figure's numbers are ``manuscripts/<ms>/figures/<fig>/results/*.json`` here (the layout of the paper repository,
where they lived until 2026-10-17), the only source the figure scripts read.  ``$UCHROM_REPRO`` points the benchmarks
at another checkout of this repository (a copy synced to a cluster, say).
"""

from __future__ import annotations

import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def paper_root() -> Path:
    """The root that holds ``manuscripts/``: this repository, or ``$UCHROM_REPRO``."""
    env = os.environ.get("UCHROM_REPRO")
    return Path(env).expanduser().resolve() if env else REPO


def results_dir(fig: str) -> Path:
    """``manuscripts/full/figures/<fig>/results`` (not created here; writers create it)."""
    return paper_root() / "manuscripts" / "full" / "figures" / fig / "results"


def paper_file(rel: str) -> Path:
    """A file under the root, e.g. ``manuscripts/full/figures/fig2/tableS_capabilities.md``."""
    return paper_root() / rel
