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

MM = 1 / 25.4
WIDTH_FULL = 180 * MM          # Nature double column

# One palette, used by every figure so a colour means the same thing throughout.
#
# The six hues are Okabe-Ito, which is the usual choice for a figure that has
# to survive colour-blind readers and greyscale printing. The three regimes
# take blue, vermillion and green, which are the three furthest apart in that
# set. The three resources of the schematic panels take purple, orange and
# ink; they never share a panel with a regime colour.
C = dict(
    ink="#1a1a1a", grey="#6e6e6e", pale="#d9d9d9", box="#f2f4f5",
    liver="#0072b2", kidney="#d55e00", pancreas="#009e73",
    supply="#cc79a7", cells="#e69f00", room="#1a1a1a",
    lung="#009e73", sim="#cc79a7",            # kept: older figures use these
    warn="#e69f00", good="#009e73", bad="#d55e00",
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


def panel_grid(fig, rows, letters, dx=0.008, dy=0.006):
    """Letter a whole figure: `rows` is a list of lists of axes.

    Each row gets one baseline and the first panel of every row the same left
    margin, so the letters form a clean column down the page.

    The letter is placed against each axes' tight bounding box, which includes
    its tick labels and its y label, rather than against the axes frame. Put
    at a fixed offset from the frame instead it lands on whichever tick label
    happens to be widest or tallest, which is what `audit` kept reporting.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.canvas.get_width_height()

    def box(a):
        bb = a.get_tightbbox(r)
        return bb.x0 / W, bb.y1 / H

    left = min(box(row[0])[0] for row in rows)
    k = 0
    for row in rows:
        top = max(box(a)[1] for a in row)
        for i, a in enumerate(row):
            x = left if i == 0 else box(a)[0]
            fig.text(max(x - dx, 0.002), min(top + dy, 0.998),
                     letters[k].upper(), fontsize=SZ["panel"],
                     fontweight="bold", va="bottom", ha="left")
            k += 1


def audit(fig, name="figure", pad=1.0, quiet=False):
    """Report text and legends that overlap, which no figure should ever have.

    Checking this by eye fails: a label that clears one panel at one font size
    runs into its neighbour at another. So every figure script calls this
    before it saves, and a hit is printed with both strings and the number of
    overlapping pixels. Tick labels count, since a crowded axis is the usual
    way two panels collide.

    Only text that is actually drawn is compared. Collecting every Text in the
    figure instead picks up tick labels matplotlib never places, which sit at
    the origin and overlap each other wholesale.

    Returns the list of overlapping pairs, so a caller can fail on it.
    """
    import itertools

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    drawn = list(fig.texts)
    for lg in getattr(fig, "legends", []):
        drawn += lg.get_texts()
    for ax in fig.axes:
        if not ax.get_visible():
            continue
        drawn += [ax.title, ax.xaxis.label, ax.yaxis.label]
        drawn += list(ax.texts)
        # only the ticks inside the view: matplotlib keeps label objects for
        # ticks it never draws, and they sit wherever they last were
        if ax.axison:
            for axis, lim in ((ax.xaxis, ax.get_xlim()),
                              (ax.yaxis, ax.get_ylim())):
                lo, hi = sorted(lim)
                drawn += [lab for loc, lab in zip(axis.get_ticklocs(),
                                                  axis.get_ticklabels())
                          if lo - 1e-9 <= loc <= hi + 1e-9]
        lg = ax.get_legend()
        if lg is not None:
            drawn += lg.get_texts()

    W, H = fig.canvas.get_width_height()
    items = []
    for t in drawn:
        if t is None or not t.get_visible() or not str(t.get_text()).strip():
            continue
        try:
            bb = t.get_window_extent(renderer=r)
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        if bb.x1 < 0 or bb.y1 < 0 or bb.x0 > W or bb.y0 > H:
            continue
        items.append((t, bb))

    hits = []
    for (t1, b1), (t2, b2) in itertools.combinations(items, 2):
        dx = min(b1.x1, b2.x1) - max(b1.x0, b2.x0) - pad
        dy = min(b1.y1, b2.y1) - max(b1.y0, b2.y0) - pad
        if dx > 0 and dy > 0:
            hits.append((t1.get_text(), t2.get_text(), dx * dy))
    if hits and not quiet:
        print(f"  {name}: {len(hits)} overlapping label pair(s)")
        for a, b, area in sorted(hits, key=lambda h: -h[2])[:12]:
            print(f"    {area:7.0f} px  {a[:38]!r}  x  {b[:38]!r}")
    return hits


def tidy(ax, hide_left=False):
    """Small consistency pass: no stray spines, no clipped tick labels."""
    if hide_left:
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    ax.set_axisbelow(True)
