"""What the agents of the benchmark get from U-Chrom, frozen at one commit: the documentation, the skills and the
wheels (``sherlock/setup_v3.sbatch`` builds the agents' environment from them).

    python supp_fig7/v3/build_snapshot.py --commit b6d8873 --out SNAPSHOT

Writes ``SNAPSHOT/{COMMIT, docs/, skills/, wheels/}``:

* ``docs/`` -- the repository README and every Markdown page of ``docs/source`` that the docs build publishes (the
  chapters by data type, the common chapters, the reference; ``discovery/`` is excluded from the build), and every tutorial notebook as Markdown with its code and text but without its
  outputs (the docs show them executed; the outputs are figures and tables of the atlas data);
* ``skills/`` -- the five skills of ``skills/`` (SKILL.md and scripts), as written before the questions of split
  "new";
* ``wheels/`` -- ``chromdata`` and ``u-chrom`` built from a clean ``git archive`` of the commit (not editable: the
  agent finds the package in site-packages, as a user who installed it would).

Nothing else of the repository goes to the agents: no tests, no benchmarks (they hold the questions' answers),
no notebooks with outputs.  Documentation that states a question's reference value is left out (``EXCLUDE``,
``REDACT``: the release notes report this benchmark's own scHiCluster numbers, R4; a validation table of the bulk
callers gives BH2's recall and precision on the same map), listed in ``SNAPSHOT/REDACTIONS.json`` (outside ``docs/``:
the agents do not see it).  Each redaction must match, so that a change of the page is noticed.
"""

from __future__ import annotations

import argparse
import io
import json
import re
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]

#: docs pages left out (path under docs/source -> reason)
EXCLUDE = {"reference/changelog.md": "release notes: report this benchmark's scHiCluster reference (R4 purity)"}
#: passages removed from docs pages: path -> [(regex, replacement, reason)]
REDACT = {
    "bulk/structures.md": [(
        r"  Against the published Rao 2014 calls on IMR90\nchr21 \(92 HiCCUPS loops, 86 Arrowhead domains\):\n\n"
        r"\| Caller \|.*?\n\n",
        "\n\n",
        "validation table: recall / precision of HiCCUPS and Mustache against Rao 2014 on IMR90 chr21 at 10 kb "
        "(BH2's reference values; BH1's loop count follows from it)")],
}


def notebook_markdown(nb: dict) -> str:
    """A notebook as Markdown: text cells as they are, code cells fenced, outputs left out."""
    out = []
    for c in nb.get("cells", []):
        src = "".join(c.get("source", []))
        if not src.strip():
            continue
        out.append(src if c.get("cell_type") == "markdown" else f"```python\n{src}\n```")
    return "\n\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--commit", required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    sha = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", a.commit], capture_output=True, text=True,
                         check=True).stdout.strip()
    if a.out.exists():
        shutil.rmtree(a.out)
    (a.out / "docs").mkdir(parents=True)
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "src"
        src.mkdir()
        tar = subprocess.run(["git", "-C", str(REPO), "archive", sha], capture_output=True, check=True).stdout
        tarfile.open(fileobj=io.BytesIO(tar)).extractall(src, filter="data")
        shutil.copy(src / "README.md", a.out / "docs" / "README.md")
        n_md = n_nb = 0
        redactions = []
        for p in sorted((src / "docs" / "source").rglob("*")):
            rel = p.relative_to(src / "docs" / "source")
            if any(part.startswith(("_", ".")) for part in rel.parts):      # _static, _ext, _templates, _build
                continue
            if rel.parts[0] == "discovery":                    # excluded from the docs build (conf.py)
                continue
            dst = a.out / "docs" / rel
            if str(rel) in EXCLUDE:
                redactions.append({"page": str(rel), "removed": "page", "reason": EXCLUDE[str(rel)]})
                continue
            if p.suffix == ".md":
                dst.parent.mkdir(parents=True, exist_ok=True)
                text = p.read_text()
                for pattern, repl, reason in REDACT.get(str(rel), []):
                    text, n = re.subn(pattern, repl, text, count=1, flags=re.S)
                    if n != 1:
                        raise SystemExit(f"redaction of {rel} did not match: {reason}")
                    redactions.append({"page": str(rel), "removed": "passage", "reason": reason})
                dst.write_text(text)
                n_md += 1
            elif p.suffix == ".ipynb":
                dst = dst.with_suffix(".md")
                dst.parent.mkdir(parents=True, exist_ok=True)
                dst.write_text(notebook_markdown(json.loads(p.read_text())))
                n_nb += 1
        for s in sorted((src / "skills").iterdir()):
            if (s / "SKILL.md").exists():
                shutil.copytree(s, a.out / "skills" / s.name, ignore=shutil.ignore_patterns("__pycache__"))
        for pkg in ("chromdata", "uchrom"):
            subprocess.run(["uv", "build", "--wheel", "-o", str(a.out / "wheels"), str(src / "packages" / pkg)],
                           check=True, capture_output=True)
    (a.out / "COMMIT").write_text(sha + "\n")
    (a.out / "REDACTIONS.json").write_text(json.dumps(redactions, indent=1) + "\n")
    print(json.dumps({"commit": sha, "markdown_pages": n_md, "tutorials": n_nb,
                      "skills": sorted(p.name for p in (a.out / "skills").iterdir()), "redactions": redactions,
                      "wheels": sorted(p.name for p in (a.out / "wheels").glob("*.whl"))}, indent=1))


if __name__ == "__main__":
    main()
