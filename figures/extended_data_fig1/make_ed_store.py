"""Extended Data Fig. 1 -- ChromData storage: capabilities and access at scale.

The store figure (``shared/store_panels/make_fig2.py``, data in ``benchmarks/store/results``) with panel letters
only, as in the journal's published figures: the panel titles, subtitles and the notes inside d and e are in the
legend (``../../full/captions/fig2.tex``, shared, with this paper's title ``\\storefigtitle`` from ``main.tex``).

    python figures/extended_data_fig1/make_ed_store.py      # -> ed_store.pdf / ed_store.png
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "shared" / "store_panels"))   # the store figure's panels
import make_fig2  # noqa: E402

if __name__ == "__main__":
    print(make_fig2.make_figure(letters_only=True, out=HERE / "ed_store"))
