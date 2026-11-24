#!/usr/bin/env python
"""Build the submission document from MANUSCRIPT.md.

MANUSCRIPT.md is the source of truth, because that is the file
`scripts/29_check_manuscript.py` holds to the regenerated tables. The Word
file is derived from it, the way the figures are derived from the results, so
it is never edited by hand: edit the markdown and build again.

The repository header is stripped on the way out. It says which model version
and which commit the numbers come from, which belongs in the repository and
not in a file sent to a journal.

Needs pandoc. If pandoc is missing the script says so and stops rather than
writing a half-converted file.

Usage:  python scripts/30_build_docx.py [--keep-notes]

    --keep-notes  keep the "Notes for the authors" section, which is what you
                  want when circulating the draft to coauthors. Without it the
                  section is dropped, which is what you want for submission.
"""
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from organgrow.config import ROOT

SRC = ROOT / "MANUSCRIPT.md"
OUT = ROOT / "manuscript"


def strip_repo_header(text):
    """Drop the two paragraphs that are about the repository, not the paper."""
    for pat in (r"\nTarget journal: \*Biofabrication\*.*?\n\n",
                r"\nReference model version .*?is written by hand\.\n\n"):
        text = re.sub(pat, "\n", text, flags=re.S)
    return text


def strip_notes(text):
    """Drop the section that exists to be removed before submission."""
    i = text.find("## Notes for the authors")
    if i < 0:
        return text
    j = text.find("\n## ", i + 1)
    return text[:i] + (text[j + 1:] if j > 0 else "")


def main():
    keep = "--keep-notes" in sys.argv[1:]
    if not SRC.exists():
        raise SystemExit(f"no manuscript at {SRC}")
    if not shutil.which("pandoc"):
        raise SystemExit(
            "pandoc is not installed, so the Word file cannot be built.\n"
            "  macOS:  brew install pandoc\n"
            "  Debian: sudo apt install pandoc\n"
            "MANUSCRIPT.md is the manuscript either way; this only converts it.")

    text = strip_repo_header(SRC.read_text(encoding="utf-8"))
    if not keep:
        text = strip_notes(text)
    OUT.mkdir(parents=True, exist_ok=True)
    name = "manuscript_with_notes" if keep else "manuscript"
    dest = OUT / f"{name}.docx"

    with tempfile.NamedTemporaryFile("w", suffix=".md", delete=False,
                                     encoding="utf-8") as f:
        f.write(text)
        tmp = f.name
    try:
        subprocess.run(["pandoc", tmp, "-f", "markdown+superscript+tex_math_dollars",
                        "-t", "docx", "-o", str(dest), "--standalone"], check=True)
    finally:
        Path(tmp).unlink(missing_ok=True)

    print(f"wrote {dest.relative_to(ROOT)}"
          + ("  (with the notes section)" if keep else "  (notes removed)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
