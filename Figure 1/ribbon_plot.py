import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import numpy as np

# ── ADJUSTABLE PARAMETERS ────────────────────────────────────────────────────
BAR_W_LEFT  = 0.25   # width of the LEFT (WGS) bar
BAR_W_RIGHT = 0.25   # width of the RIGHT (FISH) bar
FONT_SIZE   = 25     # base font size (labels, percentages, titles scale from this)
# ─────────────────────────────────────────────────────────────────────────────

# ── DATA ────────────────────────────────────────────────────────────────────
WGS = ['No Data', 'No amp', 'Amplified', 'Focal Amp']  # reversed

FISH = ['No amp (6)', 'Subclonal gain (9)', 'CN gain (34)', 'ecDNA (11)', 'HSR (5)']  # reversed


WGS_COLORS = {
    'Focal Amp': '#7b2d8b',
    'Amplified': '#c0392b',
    'No amp': '#888780',
    'No Data': '#36454F'
}

FISH_COLORS = {
    'HSR (5)': '#BEA2FE',
    'ecDNA (11)': '#b5585b',
    'CN gain (34)': '#2ca02c',
    'Subclonal gain (9)': '#e58d71',
    'No amp (6)': '#7f8c8d'
}


# Flows: data[wgs][fish] = fraction (%)
data = {
    'Focal Amp': {
        'HSR (5)': 7.7,
        'ecDNA (11)': 13.8,
        'CN gain (34)': 0.0,
        'Subclonal gain (9)': 0.0,
        'No amp (6)': 0.0
    },

    'Amplified': {
        'HSR (5)': 0.0,
        'ecDNA (11)': 1.54,
        'CN gain (34)': 16.1,
        'Subclonal gain (9)': 1.54,
        'No amp (6)': 0.0
    },

    'No amp': {
        'HSR (5)': 0.0,
        'ecDNA (11)': 0.0,
        'CN gain (34)': 25.4,
        'Subclonal gain (9)': 10.7,
        'No amp (6)': 6.2
    },

    'No Data': {
        'HSR (5)': 0.0,
        'ecDNA (11)': 1.54,
        'CN gain (34)': 10.8,
        'Subclonal gain (9)': 1.54,
        'No amp (6)': 3.1
    },
}

# ── LAYOUT CONSTANTS ────────────────────────────────────────────────────────
GAP     = 0.025     # gap between segments in each column
ALPHA   = 0.38      # ribbon fill opacity

# ── HELPERS ─────────────────────────────────────────────────────────────────
def build_segments(keys, totals, colormap, gap=GAP):
    """Return list of {key, y0, y1, color} stacked bottom-up."""
    total_height = sum(totals[k] for k in keys if totals[k] > 0)
    segs = []
    cursor = 0.0
    for k in keys:
        v = totals[k]
        if v <= 0:
            continue
        h = v / total_height
        segs.append({'key': k, 'y0': cursor, 'y1': cursor + h, 'color': colormap[k]})
        cursor += h
    # re-normalise to fill 0–1 minus inter-segment gaps
    n = len(segs)
    total_gaps = (n - 1) * GAP
    scale = 1.0 - total_gaps
    cursor = 0.0
    for s in segs:
        h = (s['y1'] - s['y0'])
        s['y0'] = cursor
        s['y1'] = cursor + h * scale
        cursor = s['y1'] + GAP
    return segs

def bezier_ribbon(ax, lx, ly0, ly1, rx, ry0, ry1, color, alpha):
    """Draw a cubic-Bezier ribbon between two vertical spans."""
    cx1, cx2 = lx + (rx - lx) * 0.45, lx + (rx - lx) * 0.55
    verts = [
        (lx, ly0),
        (cx1, ly0), (cx2, ry0), (rx, ry0),
        (rx, ry1),
        (cx2, ry1), (cx1, ly1), (lx, ly1),
        (lx, ly0),
    ]
    codes = [Path.MOVETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CLOSEPOLY]
    patch = PathPatch(Path(verts, codes), facecolor=color,
                      edgecolor='none', alpha=alpha, zorder=1)
    ax.add_patch(patch)

# ── COMPUTE TOTALS ───────────────────────────────────────────────────────────
wgs_totals  = {w: sum(data[w][f] for f in FISH) for w in WGS}
fish_totals = {f: sum(data[w][f] for w in WGS)  for f in FISH}

l_segs = build_segments(WGS,  wgs_totals,  WGS_COLORS)
r_segs = build_segments(FISH, fish_totals, FISH_COLORS)

l_map = {s['key']: s for s in l_segs}
r_map = {s['key']: s for s in r_segs}

# Track fill progress within each bar
l_off = {w: l_map[w]['y0'] for w in WGS  if w in l_map}
r_off = {f: r_map[f]['y0'] for f in FISH if f in r_map}

total_flow = sum(data[w][f] for w in WGS for f in FISH)

# ── FIGURE ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))
ax.set_xlim(0, 1)
ax.set_ylim(-0.05, 1.08)
ax.axis('off')

LX = 0.18   # left bar x-start
RX = 0.76   # right bar x-start

# Draw ribbons
for w in WGS:
    if w not in l_map:
        continue
    l_seg = l_map[w]
    for f in FISH:
        v = data[w][f]
        if v <= 0 or f not in r_map:
            continue
        r_seg = r_map[f]

        lh = (v / wgs_totals[w])  * (l_seg['y1'] - l_seg['y0'])
        rh = (v / fish_totals[f]) * (r_seg['y1'] - r_seg['y0'])

        ly0 = l_off[w]
        ly1 = ly0 + lh
        ry0 = r_off[f]
        ry1 = ry0 + rh

        bezier_ribbon(ax, LX + BAR_W_LEFT, ly0, ly1, RX, ry0, ry1, WGS_COLORS[w], ALPHA)

        l_off[w] += lh
        r_off[f] += rh

# Draw WGS bars
for s in l_segs:
    ax.add_patch(mpatches.FancyBboxPatch(
        (LX, s['y0']), BAR_W_LEFT, s['y1'] - s['y0'],
        boxstyle='round,pad=0.002', facecolor=s['color'], edgecolor='white', lw=0.5, zorder=3))
    mid = (s['y0'] + s['y1']) / 2
    pct = f"{wgs_totals[s['key']]:.1f}%"
    ax.text(LX - 0.012, mid, s['key'],
            ha='right', va='center', fontsize=FONT_SIZE, color='#222')
    ax.text(LX + BAR_W_LEFT / 2, mid, pct,
            ha='center', va='center', fontsize=FONT_SIZE * 0.85, color='white')

# Draw FISH bars
for s in r_segs:
    ax.add_patch(mpatches.FancyBboxPatch(
        (RX, s['y0']), BAR_W_RIGHT, s['y1'] - s['y0'],
        boxstyle='round,pad=0.002', facecolor=s['color'], edgecolor='white', lw=0.5, zorder=3))
    mid = (s['y0'] + s['y1']) / 2
    pct = f"{fish_totals[s['key']]:.1f}%"
    ax.text(RX + BAR_W_RIGHT + 0.012, mid, s['key'],
            ha='left', va='center', fontsize=FONT_SIZE, color='#222')
    ax.text(RX + BAR_W_RIGHT / 2, mid, pct,
            ha='center', va='center', fontsize=FONT_SIZE * 0.85, color='white')

# Column titles
ax.text(LX + BAR_W_LEFT  / 2, 1.02, 'WGS',
        ha='center', va='bottom', fontsize=FONT_SIZE * 1.17, fontweight='500', color='#111')
ax.text(RX + BAR_W_RIGHT / 2, 1.02, 'FISH',
        ha='center', va='bottom', fontsize=FONT_SIZE * 1.17, fontweight='500', color='#111')

fig.patch.set_facecolor('white')
plt.tight_layout(pad=0.5)

out_png = 'ribbon_plot.png'
out_pdf = 'ribbon_plot.pdf'
fig.savefig(out_png, dpi=180, bbox_inches='tight', facecolor='white')
fig.savefig(out_pdf,           bbox_inches='tight', facecolor='white')
print("Saved:", out_png)
print("Saved:", out_pdf)
