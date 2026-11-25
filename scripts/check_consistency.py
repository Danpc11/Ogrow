#!/usr/bin/env python
"""Check that the version and the stamped numbers of the paper agree.

Every model number the paper uses is in PAPER_NUMBERS.md, written by
`make paper` with the model version and the commit that produced it. This
checks that the version is the same in pyproject.toml, the package,
PAPER.md, CLOSED_MODEL.md and PAPER_NUMBERS.md; that the stamp names a
commit rather than "unknown" and carries no uncommitted changes; that
the stamped commit is in the history of the current one; and that nothing
under src, config or scripts has changed since it, which ancestry alone does
not give. It runs in CI on every push and needs only the standard library
and git.
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
    # The manuscript names the version of the model its numbers came from, and
    # CITATION.cff names the version a release is tagged as. The first had
    # drifted seven versions behind the code before it was checked here.
    "MANUSCRIPT.md": find(r"Reference model version ([0-9.]+)", "MANUSCRIPT.md"),
    "CITATION.cff": find(r"(?m)^version: ([0-9.]+)", "CITATION.cff"),
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

# The regime parameters are written out by hand in every script that grows a
# construct, and `config/params.yaml` says that nothing is set elsewhere. They
# are not read from one place, so the only thing that keeps them equal is that
# someone notices. One copy had already drifted: the kidney's loss rate read
# 0.0 in scripts/27_pool_fate.py and 0.05 everywhere else. This is the cheap
# half of the repair -- the copies stay, and a disagreement between them is now
# an error rather than something to notice.
REGIMES = {
    "liver":    dict(P0="2600.0", rho="0.10", delta="0.10", q_set="0.6", renewal="1.0"),
    "kidney":   dict(P0="300.0",  rho="0.00", delta="0.05", q_set="0.6", renewal="0.0"),
    "pancreas": dict(P0="60.0",   rho="0.00", delta="0.00", q_set="0.6", renewal="0.0"),
}


def num(s):
    """Compare by value, so that 0.0 and 0.00 are the same number."""
    try:
        return float(s)
    except ValueError:
        return s


_defn = re.compile(r'"(liver|kidney|pancreas)"\s*:\s*dict\(([^)]*)\)')
_bad = []
for _p in sorted((ROOT / "scripts").glob("*.py")):
    if _p.name == pathlib.Path(__file__).name:
        continue                      # the canonical table above is not a copy
    for _m in _defn.finditer(_p.read_text()):
        _reg = _m.group(1)
        _got = dict(kv.split("=", 1) for kv in
                    (x.strip() for x in " ".join(_m.group(2).split()).split(","))
                    if "=" in kv)
        for _k, _v in _got.items():
            _want = REGIMES[_reg].get(_k.strip())
            if _want is not None and num(_v.strip()) != num(_want):
                _bad.append(f"{_p.name}: {_reg} {_k.strip()}={_v.strip()}, "
                            f"the other scripts say {_want}")
print(f"  {'ok  ' if not _bad else 'FAIL'} the regime tables agree across the scripts")
fail += _bad

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
    else:
        # Ancestry alone does not mean the numbers are current: a later commit
        # can change the model and leave the stamp standing. It has happened
        # here: 5fe3834, whose message is that PAPER_NUMBERS was stamped with
        # the v4.49 commit, also changed closed.py by 57 lines. What has to
        # hold is that nothing the runs depend on has moved since the stamp.
        same = subprocess.run(["git", "diff", "--quiet", sha, "HEAD",
                               "--", "src", "config", "scripts"],
                              cwd=ROOT).returncode == 0
        print(f"  {'ok  ' if same else 'FAIL'} src, config and scripts unchanged since {sha}")
        if not same:
            moved = subprocess.run(["git", "diff", "--name-only", sha, "HEAD",
                                    "--", "src", "config", "scripts"],
                                   cwd=ROOT, capture_output=True, text=True).stdout.split()
            fail.append(f"the model has changed since the stamped commit {sha} "
                        f"({', '.join(moved[:6])}"
                        f"{f' and {len(moved) - 6} more' if len(moved) > 6 else ''}); "
                        f"run `make paper` and commit the regenerated numbers")

if fail:
    print("\n" + "\n".join(fail))
    sys.exit(1)
print("\nversion and stamp agree")
