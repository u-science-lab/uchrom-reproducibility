"""Where the benchmarks write their results: ``benchmarks/<topic>/results`` of this repository.

The figure scripts read only those files.  ``$UCHROM_REPRO`` points the benchmarks at another checkout of this
repository (a copy synced to a cluster, say).  The names of the earlier layout still work: ``fig2`` is the store
benchmark, ``fig6`` the browser benchmark, ``fig3`` the NucDynamics engine checks.
"""

from __future__ import annotations

import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TOPICS = {"fig2": "store", "fig6": "browser", "fig3": "toolkit_parity/nucdyn"}


def paper_root() -> Path:
    """The root of this repository, or ``$UCHROM_REPRO``."""
    env = os.environ.get("UCHROM_REPRO")
    return Path(env).expanduser().resolve() if env else REPO


def results_dir(topic: str) -> Path:
    """``benchmarks/<topic>/results`` (not created here; writers create it)."""
    return paper_root() / "benchmarks" / TOPICS.get(topic, topic) / "results"


def paper_file(rel: str) -> Path:
    """A file under the root."""
    return paper_root() / rel
