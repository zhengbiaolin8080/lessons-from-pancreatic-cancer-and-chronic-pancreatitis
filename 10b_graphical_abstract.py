"""Graphical abstract: gut microbial metabolites -> pancreatic TME AHR/STAT3 axis.
Layout: top row (gut -> metabolites -> TME), bottom row (evidence + outcome)."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Ellipse
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")

plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.bbox": "tight",
                     "savefig.pad_inches": 0.04})

fig, ax = plt.subplots(figsize=(7.09, 4.0))
ax.set_xlim(0, 100); ax.set_ylim(0, 56); ax.axis("off")

# palette
C_GUT  = "#EAF1EA"
C_MET  = "#FDF2E3"
C_TME  = "#EAE6F2"
C_AXIS = {"AHR": "#B03A2E", "STAT3": "#2E6DA4", "NFkB": "#4F9A6D"}


def box(x, y, w, h, fc, title, tfs=7.6, ec="#555555"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.45",
                                fc=fc, ec=ec, lw=0.9))
    if title:
        ax.text(x + w / 2, y + h - 2.2, title, ha="center", va="top",
                fontsize=tfs, fontweight="bold")


def arrow(x0, y0, x1, y1, lw=1.5, color="#333333", style="-|>", con="arc3,rad=0.0"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                                 mutation_scale=13, lw=lw, color=color,
                                 connectionstyle=con, zorder=5,
                                 shrinkA=2, shrinkB=2))


# =====================================================================
# Stage 1: Gut microbiota
# =====================================================================
box(1, 28, 21, 26, C_GUT, "Gut microbiota")
rng = np.random.default_rng(3)
for i in range(16):
    cx = 3.2 + rng.random() * 17
    cy = 33.8 + rng.random() * 14.7
    r  = 0.9 + rng.random() * 1.3
    col = ["#7FAF7A", "#A5C8A0", "#5F9160"][i % 3]
    if i % 3 == 0:
        ax.add_patch(Ellipse((cx, cy), r * 2.4, r * 1.1, angle=rng.random() * 180,
                             fc=col, ec="white", lw=0.4, zorder=2))
    else:
        ax.add_patch(Circle((cx, cy), r * 0.8, fc=col, ec="white", lw=0.4, zorder=2))
ax.text(11.5, 31.8, "dysbiosis in chronic\npancreatitis & PDAC",
        ha="center", va="top", fontsize=5.5, style="italic")

# =====================================================================
# Stage 2: Microbial metabolites
# =====================================================================
box(26, 28, 20, 26, C_MET, "Microbial\nmetabolites", tfs=7.0)
mets = [("Indoles / Trp", "#7E57C2"), ("Bile acids", "#C08A3E"),
        ("SCFAs", "#4F9A6D"), ("Polyamines", "#C05A7E")]
for i, (m, c) in enumerate(mets):
    yy = 44.5 - i * 4.6
    ax.add_patch(Circle((29.6, yy), 1.0, fc=c, ec="white", lw=0.5, zorder=3))
    ax.text(31.2, yy, m, ha="left", va="center", fontsize=6.2)
ax.text(36, 30.0, "69 focus-class metabolites", ha="center", va="top",
        fontsize=5.7, style="italic")

# arrows stage1 -> stage2 -> stage3
arrow(22.5, 41, 25.5, 41)
arrow(46.5, 41, 50.5, 41)
ax.text(48.5, 44.0, "systemic\ncirculation", ha="center", va="center",
        fontsize=4.8, style="italic", color="#666666")

# =====================================================================
# Stage 3: Pancreatic TME
# =====================================================================
box(51.5, 10, 47.5, 44, C_TME, "Pancreatic tumour microenvironment")

# tumour / PSC cell
ax.add_patch(Ellipse((74, 32), 30, 24, fc="#F5E1DA", ec="#B03A2E", lw=1.1, zorder=2))
ax.add_patch(Circle((74, 33.5), 7.0, fc="#F9F0EA", ec="#B03A2E", lw=1.0, zorder=3))
ax.text(74, 33.5, "nucleus", ha="center", va="center", fontsize=6.0, zorder=4)
ax.text(74, 23.6, "CCND1 / MYC / VEGFA / HIF1A", ha="center", fontsize=5.6,
        fontweight="bold", color="#6B2E22", zorder=4)
ax.text(74, 18.4, "PDAC / PSC cell", ha="center", fontsize=6.0, style="italic",
        color="#888888")

# AHR path
ax.add_patch(Circle((62, 40), 1.5, fc=C_AXIS["AHR"], ec="white", lw=0.6, zorder=4))
ax.text(62, 43.2, "AHR", ha="center", fontsize=6.8, fontweight="bold",
        color=C_AXIS["AHR"])
ax.text(59.8, 40, "indoles,\nbile acids", fontsize=5.0, ha="right", va="center",
        color=C_AXIS["AHR"])
arrow(63.6, 39.2, 68.6, 35.6, lw=1.3, color=C_AXIS["AHR"])

# STAT3 path
ax.add_patch(Circle((84, 40), 1.5, fc=C_AXIS["STAT3"], ec="white", lw=0.6, zorder=4))
ax.text(84, 43, "STAT3", ha="center", fontsize=6.8, fontweight="bold",
        color=C_AXIS["STAT3"])
ax.text(85.8, 40, "P-STAT3\ndimerisation", fontsize=4.8, ha="left", va="center",
        color=C_AXIS["STAT3"])
arrow(82.6, 39.2, 78.6, 35.6, lw=1.3, color=C_AXIS["STAT3"])

# TME immune / stromal components
for cx, cy, c, lab in ((57, 17.5, "#4F9A6D", "Treg"),
                       (65.5, 17.5, "#8E8E8E", "M2 mac."),
                       (89, 17.5, "#DF8A34", "CAF"),
                       (94, 17.5, "#5B7FA6", "CD8 T\nexhausted")):
    ax.add_patch(Circle((cx, cy), 2.1, fc=c, ec="white", lw=0.6, alpha=0.9, zorder=3))
    ax.text(cx, cy - 3.3, lab, ha="center", fontsize=5.2, zorder=4)
ax.text(75.5, 11.2, "immunosuppressive, desmoplastic niche", ha="center",
        fontsize=5.4, style="italic", color="#666666")

# =====================================================================
# Outcome banner (bottom-right)
# =====================================================================
ax.add_patch(FancyBboxPatch((51.5, 2.5), 47.5, 5.0, boxstyle="round,pad=0.4",
                            fc="#8E44AD", ec="#6C3483", lw=0.9))
ax.text(75.2, 5.0, "tumour progression - therapy resistance - poor survival",
        ha="center", va="center", fontsize=5.6, color="white", fontweight="bold")
arrow(74, 19.4, 74, 8.2, lw=1.4)

# =====================================================================
# Evidence strip (bottom-left)
# =====================================================================
box(1, 2.5, 48.5, 22, "#F4F4F4", "Multi-layer evidence (this study)", tfs=7.0)
ev = ["23 core targets at the metabolite-disease interface",
      "hubs TP53 / JUN / HIF1A (STRING PPI >= 0.9)",
      "AHR / STAT3 / NF-kB axis engagement",
      "TCGA-PAAD: module score -> OS HR 1.91, PFI HR 1.94",
      "core module correlates with Treg / CAF / M2",
      "Vina docking: indole-3-acrylic acid -8.5 (AHR)"]
for i, t in enumerate(ev):
    ax.text(2.8, 18.6 - i * 2.8, "\u25B6 " + t, fontsize=5.3, va="center",
            color="#333333")
# connecting arrows
arrow(23.5, 28.5, 23.5, 25.2, lw=1.2)
arrow(49.8, 13, 52.5, 18.5, lw=1.1, con="arc3,rad=-0.25")

fig.savefig(os.path.join(FIG, "Fig0_graphical_abstract.png"), dpi=300)
fig.savefig(os.path.join(FIG, "Fig0_graphical_abstract.tif"), dpi=300,
            pil_kwargs={"compression": "tiff_lzw"})
print("saved Fig0_graphical_abstract")
