# Generated illustrations

Generated 2026-10-06 with `openai/gpt-image-2.5-sunburst` through OpenRouter's images endpoint
(`../gen_image.py`, about 0.04 USD per sheet).  Every sheet has its request next to it (`<sheet>.json`: model,
prompt, size, quality).  `cut.py` cuts the sheets into one PNG per object.  The illustrations are not data.

    ( set -a; . ~/.env; set +a; python ../gen_image.py art/<sheet> "<prompt>" --background opaque )

| Sheet | Objects | Figure |
|---|---|---|
| `a_kinds_25sunburst.png` | `a_cells`, `a_tissue`, `a_nucleus` | Fig. 1a |
| `d_icons_25.png` | `c_store` (and `d_subset` ... `d_link`, not used since version 3) | Fig. 1d |
| `f_icons_25.png` | `f_notebook`, `f_laptop`, `f_web`, `f_agent` | Fig. 1d |
| `hero_25_1.png` (of `hero_25_0/1`, n = 2) | `hero`: ChromData as layers on one cells x genomic-bins base | Fig. 1a |
| `fig2a_icons.png` (2026-10-09) | `a2_atlas`, `a2_contacts`, `a2_features`, `a2_destripe`, `a2_explain` | Fig. 2a (the five steps of the workflow) |
| `fig2a_sv_0.png` (of `fig2a_sv_0/1`, n = 2; 2026-10-14) | `a2_spatial` (cropped, near-white set to white): a spot grid with a smooth pattern and a null distribution with the value in its tail | Extended Data Fig. 2a, step 5 (`spatially_variable`; replaces `a2_explain`) |

The first version of the figure used `openai/gpt-5.4-image-2` (commit 75ddd00); gpt-image-2.5 draws contact maps
as proper 45° triangles and the brain section with a more plausible hippocampus.
