"""Brief Communication, Supplementary Fig. 3 -- researchers and coding agents.

Panels a, d and e of the full manuscript's Fig. 6, relettered a, b, c, drawn by that figure's own panel functions
(``../../full/figures/fig6/make_fig6.py``) from the same data: b reads ``agent_run/`` and the Takei 2025 field of view
(``_cache/takei_fov0.npz``, ``extract_data.py``), c the agent benchmark (``results/agent_bench.json``).  Fig. 6b is
Fig. 1b of the main text and Fig. 6c Supplementary Fig. 2b, so they are left out.

Names that changed since Fig. 6 was drawn are relabelled after drawing: the browser is now the package
``uchrom_browser`` (``uchrom_browser.connect()``, was ``uchrom.browser.connect()``), the analysis modules are the six
of Fig. 1c (``spatial`` in place of ``pl``), and contact maps can be embedded in the store as well as linked.  The
recorded session of b keeps the code it ran (``ub.connect()``, ``ws.selected_cells()``; the same calls today).

    python manuscripts/nmeth-bc/figures/make_supp_agents.py          # -> supp_agents.pdf / supp_agents.png
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
FIG6 = HERE.parents[1] / "full" / "figures" / "fig6"
H = 121.0                                   # layout height (mm) of the three panels; the width stays 170 mm

_spec = importlib.util.spec_from_file_location("make_fig6", FIG6 / "make_fig6.py")
m6 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m6)                # draws nothing; makes the full-size canvas, replaced below

# a shorter canvas in the same layout units (mm), as make_fig6 sets it up
plt.close(m6.fig)
m6.H = H
m6.PW, m6.PH = m6.W, H * m6.SCALE + 2 * m6.MARGIN
m6.fig = plt.figure(figsize=(m6.PW * m6.MM, m6.PH * m6.MM))
m6.A = m6.fig.add_axes([0, 0, 1, 1])
m6.A.set_xlim(-m6.MARGIN / m6.SCALE, m6.PW / m6.SCALE - m6.MARGIN / m6.SCALE)
m6.A.set_ylim(m6.PH / m6.SCALE - m6.MARGIN / m6.SCALE, -m6.MARGIN / m6.SCALE)
m6.A.axis("off")

# the panels keep their drawing code; their letters change and their titles go (the legend says what each shows)
LETTER = {"a": "a", "d": "b", "e": "c"}
LETTER_DY = 5.0                         # the content moves up into the titles' room; the letter stays at its place
m6.panel = lambda x, y, letter, title=None, sub=None: m6.text(x, y + LETTER_DY, LETTER[letter], 8, weight="bold",
                                                              va="top")

RELABEL = {
    "WebMCP · uchrom.browser.connect()": "WebMCP · uchrom_browser.connect()",
    "uchrom.browser.connect() → views as evidence": "uchrom_browser.connect() → views as evidence",
    "dashed: linked, not copied": "dashed: linked (or embedded)",
    "pl": "spatial",
}


def relabel(fig) -> None:
    """Replace the outdated labels (``RELABEL``) wherever make_fig6 drew them; fail if one is missing."""
    seen = set()
    for ax in fig.axes:
        for t in ax.texts:
            s = t.get_text()
            if s in RELABEL:
                t.set_text(RELABEL[s])
                seen.add(s)
                if s == "pl":                   # a longer name in the same module box
                    t.set_fontsize(t.get_fontsize() * 0.82)
    missing = set(RELABEL) - seen
    if missing:
        raise SystemExit(f"labels not found in make_fig6's drawing (changed?): {sorted(missing)}")


def main() -> None:
    m6.panel_a(0, 2 - LETTER_DY)
    m6.panel_d({}, x0=0, y0=59 - 2 * LETTER_DY)
    r = m6.load_results()
    if "agent_bench" not in r:
        raise SystemExit("../../full/figures/fig6/results/agent_bench.json missing")
    m6.panel_e(r["agent_bench"], x0=108, y0=59 - 2 * LETTER_DY)
    relabel(m6.fig)
    m6.fig.savefig(HERE / "supp_agents.pdf")
    m6.fig.savefig(HERE / "supp_agents.png", dpi=300)


if __name__ == "__main__":
    main()
