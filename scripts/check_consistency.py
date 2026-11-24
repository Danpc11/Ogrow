#!/usr/bin/env python
"""Check that the version and the stamped numbers of the paper agree.

Every model number the paper uses is in PAPER_NUMBERS.md, written by
`make paper` with the model version and the commit that produced it. This
checks that the version is the same in pyproject.toml, the package,
PAPER.md, CLOSED_MODEL.md and PAPER_NUMBERS.md; that the stamp names a
commit rather than "unknown" and carries no uncommitted changes; and that
the stamped commit is in the history of the current one. It runs in CI on
every push and needs only the standard library and git.
"""
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
fail = []


def read(name):
    return (ROOT / name).read_text(encoding="utf-8")


def find(pattern, name):
    m = re.search(pattern, read(name))
    return m.group(1) if m else None


version = re.search(r'^version = "([^"]+)"', read("pyproject.toml"), re.M).group(1)
found = {
    "src/organgrow/__init__.py": find(r'__version__ = "([^"]+)"', "src/organgrow/__init__.py"),
    "PAPER.md": find(r"Reference model: version ([0-9.]+)", "PAPER.md"),
    "CLOSED_MODEL.md": find(r"# The closed model \(version ([0-9.]+)\)", "CLOSED_MODEL.md"),
    "PAPER_NUMBERS.md": find(r"Model version ([0-9.]+)", "PAPER_NUMBERS.md"),
}
# No figure may name a colour of its own: one palette lives in
# src/organgrow/figstyle.py, and a literal here is how Figure 4 came to draw
# the liver green while every other figure drew it blue.
_lit = re.compile(r"#[0-9a-fA-F]{6}")
for _f in sorted((ROOT / "figures").glob("fig_*.py")):
    _hits = _lit.findall(_f.read_text(encoding="utf-8"))
    if _hits:
        fail.append(f"{_f.relative_to(ROOT)} names {len(_hits)} colour(s) of its "
                    f"own ({', '.join(sorted(set(_hits))[:4])}); use figstyle.C")

print(f"pyproject.toml version {version}")
for name, v in found.items():
    ok = v == version
    print(f"  {'ok  ' if ok else 'FAIL'} {name}: {v}")
    if not ok:
        if name == "PAPER_NUMBERS.md":
            fail.append(f"{name} holds numbers from version {v}, while the code is "
                        f"{version}. Run `make paper` and commit the regenerated "
                        f"PAPER_NUMBERS.md and source_data together with the version "
                        f"change, so that no commit carries numbers from two versions.")
        else:
            fail.append(f"{name} says {v}, pyproject.toml says {version}")

stamp = re.search(r"Model version [0-9.]+, commit ([^,]+),", read("PAPER_NUMBERS.md"))
commit = stamp.group(1).strip() if stamp else "missing"
if commit in ("missing", "unknown") or "uncommitted" in commit:
    print(f"  FAIL stamp commit: {commit}")
    fail.append(f"PAPER_NUMBERS.md stamp is '{commit}'; run `make paper` on a clean commit")
else:
    sha = commit.split()[0]
    known = subprocess.run(["git", "merge-base", "--is-ancestor", sha, "HEAD"],
                           cwd=ROOT).returncode == 0
    print(f"  {'ok  ' if known else 'FAIL'} stamp commit {sha} is in the history of HEAD")
    if not known:
        fail.append(f"stamped commit {sha} is not an ancestor of HEAD")

if fail:
    print("\n" + "\n".join(fail))
    sys.exit(1)
print("\nversion and stamp agree")
