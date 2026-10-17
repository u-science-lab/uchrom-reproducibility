# Supplementary Fig. 7: the agent benchmark

Also used by Supplementary Note 6.

* The benchmark (version 1): `questions.yaml`, `ground_truth.py` -> `answers.json` (the answer key), `run_bench.py`
  (the agents' runs into `runs/`: transcripts stay local, `grades.json` are tracked), `grade.py`, `summarize.py`
  -> `results/agent_bench.json`.
* `session/`: the handoff session (one agent session, its notebook and its record).
