# Supplementary Fig. 7: the agent benchmark

Also used by Supplementary Note 6.

* The benchmark (version 1): `questions.yaml`, `ground_truth.py` -> `answers.json` (the answer key), `run_bench.py`
  (the agents' runs into `runs/`: transcripts stay local, `grades.json` are tracked), `grade.py`, `summarize.py`
  -> `results/agent_bench.json`.
* `session/`: the handoff session (one agent session, its notebook and its record).
* Version 3 (`v3/`, 2026-10-13/14): 27 questions in two parts -- Part A reads the data (anchors A1, A6, H1, S3 and the
  store-reading questions F1-F4), Part B asks for 19 published methods by name, each with a reference computed by
  U-Chrom **and** by the authors' code (`v3/reports/*.md`).  Codex CLI (gpt-6-astra / sol / luna) in four conditions:
  `none` (no U-Chrom), `bare` (U-Chrom installed, no documentation), `docs`, `skills`; 3 runs each, 972 sessions,
  1,800 s per session, run on Sherlock (`v3/sherlock/codex_run.sbatch`, snapshot of U-Chrom fb43ac4).  Graded by
  `v3/grade_v2.py` with `answers_v3.json` / `questions_v3.yaml` (numbers tested in code, names and rubric items by a
  `claude -p` judge, claude-opus-5-5): `v3/results/main_v3_grades.json`; `v3/summarize_v3.py` ->
  `results/agent_bench_v3.json` (question as the unit, Wilcoxon signed-rank over questions).  The plan is
  `v3/questions_v3_plan.md`; how each question and its data were built, `v3/README.md`, `v3/data/*.md`.
  Session transcripts, pilots and the inputs built for v3 are archived on Oak,
  `$OAK/uchrom/results/agent_bench/2026-10-14_v3/` (`v3/sherlock/archive_v3.sbatch`).
* Version 2 (2026-10-13; 43 questions, Claude models, same four conditions) showed no difference between conditions
  with the question as the unit; most of its discriminating questions turned on what their wording left open, which
  version 3 was built to remove (`v3/questions_v3_plan.md`).  Its summary `results/agent_bench_v2.json` stays local
  (not tracked); its sessions are archived on Oak (`agent_bench/2026-10-13_fb43ac4`).
