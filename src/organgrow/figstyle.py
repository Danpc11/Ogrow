"""One visual style for every figure, and panel letters that really align.

Panel letters are set in upper case, as the journal requires.

Panel letters are placed in FIGURE coordinates, taken from each axes'
bounding box, so a row of panels shares one baseline and a column shares one
left margin whatever the individual axes widths are.

Type sizes are set for the final printed size. A figure declared at 180 mm
is placed at 180 mm, so 7 pt here is 7 pt on the page; a figure declared at
90 mm would be, so the same rcParams apply and nothing is rescaled.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

MM = 1 / 25.4
WIDTH_FULL = 180 * MM          # Nature double column
WIDTH_HALF = 88 * MM

# One palette, used by every figure so a colour means the same thing throughout.
C = dict(
    ink="#1a1a1a", grey="#8a8a8a", pale="#dcdcdc", box="#f2f4f5",
    liver="#c0392b", kidney="#2874a6", pancreas="#b7950b", lung="#117a65",
    sim="#7d5bb8", warn="#b7950b", good="#117a65", bad="#c0392b",
)

# Sizes in points at final print size. Nothing below 6 pt.
SZ = dict(panel=8, label=7, tick=6.5, legend=6.5, annot=6.5, small=6)


def apply():
    """Set the shared rcParams. Call once at the top of a figure script."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": SZ["label"],
        "axes.labelsize": SZ["label"],
        "axes.titlesize": SZ["label"],
        "xtick.labelsize": SZ["tick"],
        "ytick.labelsize": SZ["tick"],
        "legend.fontsize": SZ["legend"],
        "legend.frameon": False,
        "legend.handlelength": 1.3,
        "legend.borderpad": 0.15,
        "legend.labelspacing": 0.25,
        "legend.handletextpad": 0.5,
        "legend.borderaxespad": 0.3,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.2, "ytick.major.size": 2.2,
        "xtick.major.pad": 2, "ytick.major.pad": 2,
        "axes.labelpad": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.0,
        "patch.linewidth": 0.7,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.01,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def panel(fig, ax, letter, dx=0.030, dy=0.012):
    """Place a panel letter in figure coordinates, above and left of `ax`.

    Because the position comes from the axes bounding box, every panel in a
    row lands on the same baseline and every panel in a column on the same
    left margin, regardless of tick labels or axis titles.
    """
    fig.canvas.draw()                       # bounding boxes need a draw
    bb = ax.get_position()
    fig.text(bb.x0 - dx, bb.y1 + dy, letter.upper(), fontsize=SZ["panel"],
             fontweight="bold", va="bottom", ha="left")


def panel_row(fig, axes, letters, dx=0.030, dy=0.012):
    """Letter a row of panels on one shared baseline."""
    fig.canvas.draw()
    top = max(a.get_position().y1 for a in axes)
    for a, L in zip(axes, letters):
        fig.text(a.get_position().x0 - dx, top + dy, L.upper(),
                 fontsize=SZ["panel"], fontweight="bold", va="bottom",
                 ha="left")


def panel_grid(fig, rows, letters, dx=0.030, dy=0.012):
    """Letter a whole figure: `rows` is a list of lists of axes.

    Each row gets one baseline; the first panel of every row gets the same
    left margin, so the letters form a clean column down the page.
    """
    fig.canvas.draw()
    left = min(r[0].get_position().x0 for r in rows)
    k = 0
    for row in rows:
        top = max(a.get_position().y1 for a in row)
        for i, a in enumerate(row):
            x = left if i == 0 else a.get_position().x0
            fig.text(x - dx, top + dy, letters[k].upper(),
                     fontsize=SZ["panel"], fontweight="bold", va="bottom",
                     ha="left")
            k += 1


def tidy(ax, hide_left=False):
    """Small consistency pass: no stray spines, no clipped tick labels."""
    if hide_left:
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    ax.set_axisbelow(True)
