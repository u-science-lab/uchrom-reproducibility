"""Cut the generated illustration sheets (``gen_image.py``, prompts in ``PROMPTS.md``) into one PNG per object.

Objects are found automatically, separated by white gutters (rows of the sheet first, then columns within a
row), in reading order; their number must match the names given below.  Each object is cropped
with a margin; the generator's off-white paper (the near-white pixels connected to the border) is set to pure
white, fills inside the objects are kept.  The figure panels are white, so no transparency is needed.

    python shared/art/cut.py
"""
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

HERE = Path(__file__).resolve().parent
SHEETS = {  # sheet -> names in reading order (row by row, left to right)
    "a_kinds_25sunburst": ["a_cells", "a_tissue", "a_nucleus"],
    "f_icons_25": ["f_notebook", "f_laptop", "f_web", "f_agent"],
    "d_icons_25": ["d_subset", "d_aggregate", "d_convert", "d_call", "d_features", "d_link", "c_store"],
    "hero_25_1": ["hero"],
    "fig2a_icons": ["a2_atlas", "a2_contacts", "a2_features", "a2_destripe", "a2_explain"],
    "ed3a_icons": ["x3_recon", "x3_trace_v1", "x3_separation", "x3_classes", "x3_cells_v1"],
    "ed3a_icons2": ["x3_trace", "x3_cells"],      # an open trace; the heat map whole, the check mark beside it
}


def runs(ink: np.ndarray, gap: int):
    """(start, stop) of the runs of True in a 1-D mask, joining runs separated by fewer than gap False."""
    idx = np.where(ink)[0]
    if not len(idx):
        return []
    cut = np.where(np.diff(idx) > gap)[0]
    return list(zip(np.r_[idx[0], idx[cut + 1]], np.r_[idx[cut], idx[-1]] + 1))


def objects(a: np.ndarray, gap: int):
    """Objects separated by white gutters: rows of the sheet first, then columns within each row."""
    ink = a.min(2) < 225
    out = []
    for r0, r1 in runs(ink.any(1), gap):
        for c0, c1 in runs(ink[r0:r1].any(0), gap):
            sub = ink[r0:r1, c0:c1]
            rr = np.where(sub.any(1))[0]
            out.append((r0 + rr[0], c0, r0 + rr[-1] + 1, c1))
    return [b for b in out if (b[2] - b[0]) * (b[3] - b[1]) > 0.003 * a.shape[0] * a.shape[1]]


def clean(im: Image.Image) -> Image.Image:
    a = np.asarray(im).copy()
    lab, _ = ndimage.label(a.min(2) >= 236)
    edge = np.unique(np.r_[lab[0], lab[-1], lab[:, 0], lab[:, -1]])
    a[np.isin(lab, edge[edge > 0])] = 255
    return Image.fromarray(a)


def main() -> None:
    for sheet, names in SHEETS.items():
        im = Image.open(HERE / f"{sheet}.png").convert("RGB")
        a = np.asarray(im)
        for gap in (12, 20, 30, 45, 8):
            boxes = objects(a, gap)
            if len(boxes) == len(names):
                break
        else:
            raise SystemExit(f"{sheet}: found {len(boxes)} objects, expected {len(names)}")
        m = 8
        for name, (r0, c0, r1, c1) in zip(names, boxes):
            out = clean(im.crop((max(c0 - m, 0), max(r0 - m, 0), min(c1 + m, im.width), min(r1 + m, im.height))))
            out.save(HERE / f"{name}.png")
            print(sheet, name, out.size, f"(gap {gap})")


if __name__ == "__main__":
    main()
