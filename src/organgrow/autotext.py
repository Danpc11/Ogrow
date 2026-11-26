"""Passages of the manuscript that a script writes from its results.

A number typed into prose goes stale the first time the results move, and
`scripts/29_check_manuscript.py` only finds out afterwards. The passages that
quote the most numbers are therefore written by the script that computes
them, between a pair of markers:

    <!-- auto:NAME -->  ...  <!-- /auto:NAME -->

Writing the same text twice changes nothing, so a run that finds the
manuscript already in step leaves the file, and the stamp, untouched.
"""
import re


def fill(text, name, body):
    """Replace what stands between the two markers of one named block."""
    pat = re.compile(rf"(<!-- auto:{re.escape(name)} -->)(.*?)(<!-- /auto:{re.escape(name)} -->)",
                     re.S)
    if not pat.search(text):
        return text
    return pat.sub(lambda m: m.group(1) + body + m.group(3), text)


def write(doc, blocks):
    """Bring the named passages of `doc` in step. True if the file changed."""
    if not doc.exists():
        return False
    old = doc.read_text(encoding="utf-8")
    new = old
    for name, body in blocks.items():
        new = fill(new, name, body)
    if new == old:
        return False
    doc.write_text(new, encoding="utf-8")
    return True
