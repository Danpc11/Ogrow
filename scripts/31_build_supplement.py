#!/usr/bin/env python
"""Build the supplement of the second manuscript, and its Word file.

Every table of the supplement is a section of `PAPER_NUMBERS.md`, copied as it
stands, so the supplement cannot hold a number the stamped file does not. The
figure legend and the note come from `manuscript/MANUSCRIPT_A.md`.

    manuscript/SUPPLEMENT_A.md    always written
    manuscript/supplement_A.docx  if pandoc is installed; set in landscape
                                  with small type in the tables if python-docx
                                  is installed too, because several tables
                                  have a dozen columns and break their words
                                  on an upright page

Run by `make paper` after `scripts/20_paper_numbers.py`. Usage:
    python scripts/31_build_supplement.py
"""
import re
import shutil
import subprocess

from organgrow.config import ROOT

PN = ROOT / "PAPER_NUMBERS.md"
MS = ROOT / "manuscript" / "MANUSCRIPT_A.md"
OUT = ROOT / "manuscript" / "SUPPLEMENT_A.md"
DOCX = ROOT / "manuscript" / "supplement_A.docx"

#: (number, title, heading in PAPER_NUMBERS.md, heading at which the copy stops or None)
TABLES = (
    ("S1", "The factorial of renewal, seeded pool, and turnover, with the elasticity of "
           "each source, related to Figure 2",
     "What separates the regimes: renewal against pool and turnover", None),
    ("S2", "Nearness to capacity against size: the perfusion threshold and the domain "
           "changed inside one source, related to Figure 2",
     "Nearness to capacity against size", None),
    ("S3", "Elasticities of the final unit count under alternative settings of the "
           "nucleation law, related to Figure 2",
     "Elasticities of the final count (mean ± s.e., paired replicates)", None),
    ("S4", "Robustness of the three sources to the auxiliary rules of the model, related "
           "to Figures 1 and 2",
     "Robustness to the auxiliary rules", "Interventions on the closed state"),
    ("S5", "The reference model beside published progenitor ablations and paired human "
           "kidneys, related to Figure 3",
     "Against published perturbations and paired human kidneys", None),
    ("S6", "The extension in which progenitors drive branching, in every setting "
           "examined, related to Figure 4 and Table 2",
     "An extension: progenitors that make their tip branch", None),
    ("S7", "The dividing kidney-like source with and without the coupling, related to "
           "Figure 4 and Table 2",
     "The kidney-like dividing source with and without the coupling", None),
    ("S8", "Predictions for ablations that have not been performed, related to Figure 5",
     "Predictions for ablations that have not been done", None),
)

#: What the four shares of the pool are, because they are not a partition.
SHARES = (
    "**How to read the shares of the pool.** The columns *supply*, *room* and *active* "
    "divide the progenitors a construct still holds at arrest among tips whose perfusion "
    "has fallen to the threshold, tips above it with no free site, and tips that can "
    "still form a unit. They are mutually exclusive, and each is given as a share of the "
    "pool that was seeded. *Spent* is the part of the seeded pool no longer held, one "
    "minus the pool held over the pool seeded, and is set to zero when the pool has "
    "grown. In a source that is spent the four add to one. In a source that divides, the "
    "pool held can exceed the pool seeded, so the three held shares can add to more than "
    "one while *spent* is zero: 0.21 + 1.03 + 0.18 in one row is a pool that has grown "
    "to about 1.2 times what was seeded, not an error of accounting. The shares are "
    "averaged over seeds after *spent* is set to zero where the pool has grown, so a row "
    "can show a spent share above zero beside held shares that add to more than one. "
    "The bound named in "
    "the last column is the largest of the four.")


def section(pn, title, stop=None):
    """One section of PAPER_NUMBERS.md, without its own heading."""
    m = re.search(r"^## " + re.escape(title) + r".*$", pn, flags=re.M)
    if not m:
        return None
    j = pn.find("\n## " + stop, m.end()) if stop else pn.find("\n## ", m.end())
    body = pn[m.end():j if j > 0 else len(pn)].strip()
    return re.sub(r"^## (.+)$", r"**\1**", body, flags=re.M)


def build():
    pn, ms = PN.read_text(encoding="utf-8"), MS.read_text(encoding="utf-8")
    title = ms.split("\n", 1)[0].lstrip("# ").strip()
    stamp = pn[:pn.index("\n## ")].strip().split("\n", 1)[1].strip()
    supp = ms[ms.index("## SUPPLEMENTAL INFORMATION"):ms.index("## STAR METHODS")]
    supp = re.sub(r"<!--.*?-->", "", supp)
    fig = supp[supp.index("**Figure S1."):supp.index("**Note S1.")].strip()
    note = supp[supp.index("**Note S1."):].strip().split("\n", 1)[1].strip()
    # the supplement has its own short reference list; equations as plain text
    refs = ms[ms.index("## REFERENCES"):ms.index("## Notes for the authors")]
    cited = sorted({int(x) for m in re.finditer(r"\[(\d+(?:, \d+)*)\]", fig + note)
                    for x in m.group(1).split(", ")})
    renum = {n: f"S{i + 1}" for i, n in enumerate(cited)}
    def swap(t):
        return re.sub(r"\[(\d+(?:, \d+)*)\]",
                      lambda m: "[" + ", ".join(renum[int(x)] for x in m.group(1).split(", ")) + "]", t)
    fig, note = swap(fig), swap(note)
    note = re.sub(r"\$\\rho = ([+\-0-9.]+)\$, \$P = ([0-9.]+)\$", r"ρ = \1, P = \2", note)
    reflist = []
    for n in cited:
        m = re.search(rf"^\[{n}\] (.+)$", refs, flags=re.M)
        if m:
            reflist.append(f"[{renum[n]}] {m.group(1)}")
    # the author line is the first paragraph after the title; affiliation marks off
    authors = ms.split("\n\n", 2)[1]
    authors = re.sub(r"\^[^^]*\^|\\\*", "", authors).replace(" ,", ",")
    out = ["# Supplemental information\n", f"**{title}**\n", " ".join(authors.split()) + "\n",
           f"This file contains Figure S1 (provided as a separate file), Notes S1 and S2, and "
           f"Tables S1 to S{len(TABLES)}. Every table is a section of `PAPER_NUMBERS.md` as "
           f"stamped, and is rebuilt by `make paper`.\n", stamp + "\n", "---\n",
           "## Figure S1\n", fig + "\n", "## Note S1. Variation between two organs grown alike\n",
           note + "\n", "## Note S2. How to read the tables\n",
           "All values are from 24 independent seeds of the stochastic model unless stated. "
           "Counts and shares are means with their standard errors; contrasts between a "
           "perturbed construct and its reference use the same seeds and are ratios of means "
           "with 95% bootstrap intervals. The three cell sources are named liver-like, "
           "kidney-like, and pancreas-like for the tissues whose cell source they resemble; no "
           "parameter was optimized against the measurements with which the model is "
           "compared.\n", SHARES + "\n"]
    # numbered in sequence over the tables that exist, so that a section not yet
    # in PAPER_NUMBERS.md leaves no gap and the count announced above is true
    missing, present = [], []
    for _, ttl, head, stop in TABLES:
        body = section(pn, head, stop)
        if body is None:
            missing.append(head)
        else:
            present.append((ttl, body))
    for i, (ttl, body) in enumerate(present, 1):
        out += [f"## Table S{i}. {ttl}\n", body + "\n"]
    out[3] = out[3].replace(f"Tables S1 to S{len(TABLES)}", f"Tables S1 to S{len(present)}")
    if reflist:
        out += ["## Supplemental references\n", "\n\n".join(reflist) + "\n"]
    OUT.write_text("\n".join(out), encoding="utf-8")
    return missing


def to_docx():
    if shutil.which("pandoc") is None:
        return "pandoc is not installed, so supplement_A.docx was not built"
    subprocess.run(["pandoc", str(OUT), "-f", "markdown+superscript+tex_math_dollars",
                    "-t", "docx", "-o", str(DOCX), "--standalone"], check=True)
    try:
        from docx import Document
        from docx.enum.section import WD_ORIENT
        from docx.shared import Inches, Pt
    except ImportError:
        return ("wrote supplement_A.docx upright: `pip install python-docx` and run again "
                "for the landscape layout the wide tables need")
    d = Document(str(DOCX))
    for sec in d.sections:
        sec.orientation = WD_ORIENT.LANDSCAPE
        sec.page_width, sec.page_height = Inches(11), Inches(8.5)
        for m in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
            setattr(sec, m, Inches(0.7))
    for t in d.tables:
        t.autofit = True
        for row in t.rows:
            for c in row.cells:
                for p in c.paragraphs:
                    p.paragraph_format.space_after = Pt(1)
                    p.paragraph_format.space_before = Pt(1)
                    for r in p.runs:
                        r.font.size = Pt(7.5)
    d.save(str(DOCX))
    return "wrote manuscript/supplement_A.docx, landscape"


if __name__ == "__main__":
    if not (PN.exists() and MS.exists()):
        raise SystemExit("needs PAPER_NUMBERS.md and manuscript/MANUSCRIPT_A.md")
    gone = build()
    print(f"wrote {OUT.relative_to(ROOT)}"
          + (f"; not in PAPER_NUMBERS.md yet: {', '.join(gone)}" if gone else ""))
    print(to_docx())
