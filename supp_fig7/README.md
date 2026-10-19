# Supplementary Fig. 7: interfaces for agents and the agent benchmark

Also used by Supplementary Note 6.

* `session/` (panel b): the handoff session, one agent session started from a selection made in the browser
  (`session/run_d.py`), with its transcript, answer, screenshots and record.
* `v3/` (panel c): the agent benchmark.  27 questions in two parts: Part A reads the data (anchors A1, A6, H1, S3 and
  the store-reading questions F1-F4), Part B asks for 19 published methods by name, each with a reference computed by
  U-Chrom **and** by the authors' code (`v3/reports/*.md`).  Codex CLI 0.156 (gpt-6-astra / sol / luna) in four
  conditions: `none` (no U-Chrom), `bare` (U-Chrom installed, no documentation), `docs`, `skills`; 3 runs each, 972
  sessions, 1,800 s per session, run on Sherlock (`v3/sherlock/codex_run.sbatch`; the agents get U-Chrom at b6d8873,
  `v3/build_snapshot.py`).  Graded by `v3/grade.py` with `v3/answers_v3.json` / `v3/questions_v3.yaml` (numbers tested
  in code, names and rubric items by a `claude -p` judge, claude-opus-5-5, which does not see the condition):
  `v3/results/main_v3_grades.json`; `v3/summarize_v3.py` -> `results/agent_bench_v3.json` (the question as the unit,
  Wilcoxon signed-rank over questions).  How each question and its data were built: `v3/README.md`,
  `v3/data/*.md`; the anchors' references: `v3/truth/anchors.py` -> `v3/answers_anchors.json`.
  Session transcripts, pilots and the inputs built for the benchmark are archived on Oak,
  `$OAK/uchrom/results/agent_bench/2026-10-14_v3/` (`v3/sherlock/archive_v3.sbatch`).
