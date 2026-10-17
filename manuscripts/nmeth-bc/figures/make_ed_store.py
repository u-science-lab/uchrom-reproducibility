"""Brief Communication, Extended Data Fig. 1 -- ChromData storage: capabilities and access at scale.

The full manuscript's Fig. 2 (``../../full/figures/fig2/make_fig2.py``, data in its ``results/``) with panel letters
only, as in the journal's published figures: the panel titles, subtitles and the notes inside d and e are in the
legend (``../../full/captions/fig2.tex``, shared, with this paper's title ``\\storefigtitle`` from ``main.tex``).

    python manuscripts/nmeth-bc/figures/make_ed_store.py      # -> ed_store.pdf / ed_store.png
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "full" / "figures" / "fig2"))
import make_fig2  # noqa: E402

if __name__ == "__main__":
    print(make_fig2.make_figure(letters_only=True, out=HERE / "ed_store"))
