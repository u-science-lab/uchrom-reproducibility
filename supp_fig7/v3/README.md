# Agent benchmark, version 3: building the questions

Plan and question table: `../questions_v3_plan.md` (approved 2026-10-13, extended Part A F1-F4 and Part B to 20).
The agent of version 3 is Codex CLI (gpt-6-astra / gpt-6-sol / gpt-6-luna) in the same four conditions (none, bare,
docs, skills) as version 2, 1,800-s limit per session.  Each family of questions is built in its own files, merged
into `questions_v3.yaml` / `answers_v3.json` at the end:

```
v3/questions/<family>.yaml      the questions (schema below)
v3/truth/<family>.py            reference answers with U-Chrom + the independent check (one function per question)
v3/sherlock/<family>.sbatch     the job that runs it on Sherlock
v3/results/<family>.json        {id: {"values": ..., "verify": {<implementation>: ...}, "truth": fn}} as answers_v2.json
v3/data/<family>.md             what each new data folder holds and how it is built (from registered sources)
v3/reports/<family>.md          per question: reference, independent value, agreement, runtime, decision
```

## Rules for a question

1. **Name the method, never U-Chrom.**  A method question names the published method, its paper, every parameter
   that changes the answer and how the input is prepared, so that the authors' implementation (or any correct one)
   gives the reference.  Words like "U-Chrom", "uchrom", "ChromData API" never appear in a Part B prompt (a Part A
   prompt may say "ChromData store", the format's name).
2. **Two implementations.**  The reference is computed with U-Chrom (the version the agents get) **and** with an
   independent implementation: the authors' code (GitHub / PyPI, pinned to a commit or version) or, for a short
   definition, a separate numpy implementation written from the paper.  Grade types:
   - `numbers` with `tol` (relative): the two must agree within `tol` (state both values in `verify`);
   - `numbers_min` with `<name>_floor` in `values`: quality against an external truth (published structure, published
     loops, cell types); the floor is set below what the named method reaches (U-Chrom's value minus a margin
     justified in the report), and the independent implementation must also reach it;
   - `text`: a name (cell type, chromosome); `rubric`: judged items (use sparingly).
   A question whose implementations disagree beyond the tolerance is fixed (tolerance covering both, with reason) or
   dropped; say which in the report.
3. **Feasible without U-Chrom.**  The authors' code must be installable by an agent (pip / git clone, Python 3.12,
   Debian bookworm with gcc, no root, internet) and finish within ~10 min on 8 CPUs for the reference input; record
   its runtime.  If it cannot, the question is conditional and the report says so.
4. **Real data only, registered.**  Inputs come from sources registered in
   `packages/uchrom/uchrom/datasets/_sources.py` / `apps/atlas/recipes/README.md` or the public atlas
   (`https://uchrom-atlas-r2.u-science.org/<id>.chromdata.zarr`).  A new source is registered in the same work (entry
   + README section).  Never fabricate or simulate an input.  A derived input (a slice, a table) is built by a script
   in `v3/` and described in `v3/data/<family>.md`.
5. **Unambiguous.**  Version 2 lost its discrimination to underspecified wording (mean or median, binning, which
   measure, how to fit).  State units, bin edges, half-open ranges, which statistic, which subset, and the exact
   quantity to report.

## Question schema (as version 2's `questions.yaml`)

```yaml
- id: R3
  category: schic            # imaging | schic | bulk | crossmodal | spatial | format
  split: v3
  kind: method               # method | format | anchor
  data: schicar              # data key(s) of the runner ('a+b'); new keys: give folder + file in v3/data/<family>.md
  prompt: '...'
  truth: rwr_pair_corr       # function in v3/truth/<family>.py
  grade:
    numbers: [r]
    tol: 0.02
  method_ref: random walk with restart (scHiCluster; Zhou et al. 2019, PNAS); schicluster 1.3.x
  check: parity              # parity | floor
```

## Sherlock

Follow `/etc/agents/AGENTS.md` (read it first).  No pip / compile on login nodes (jobs only), `--time` on every job,
polls of the queue at least 60 s apart, downloads inside jobs, transfers through `sherlock-dtn`.  Partition
`xiaojie,normal`; logs to `$SCRATCH/uchrom-nd/agentbench/logs/`.  Image `$GROUP_HOME/$USER/uchrom-nd/agentbench/py312full.sif`
(python:3.12-bookworm); the agents' U-Chrom environment `.../agentbench/venv-uchrom` (commit fb43ac4; bind at `/opt/venv`,
as `sherlock/truth2.sbatch` does); a cooltools environment `.../agentbench/venv-verify`; other existing envs in
`$GROUP_HOME/$USER/uchrom-nd/` (e.g. `venv-snaphic-ref`: the authors' SnapHiC; `minimds`).  New environments for
authors' code go to `$GROUP_HOME/$USER/uchrom-nd/agentbench/venv-<tool>` (built inside a job).  Data root of version 2:
`$SCRATCH/uchrom-nd/agentbench/data/`; version 3 adds folders under `$SCRATCH/uchrom-nd/agentbench/data_v3/`.  Shared
data on Sherlock: `$SCRATCH/uchrom-nd/data/` (read only).  Sync code with
`rsync -a <files> sherlock-dtn:$SCRATCH/uchrom-nd/agentbench/bench/v3/...`.
