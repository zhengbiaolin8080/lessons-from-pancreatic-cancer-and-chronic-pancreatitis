"""Step 10: main figures Fig1-Fig6, Fig8 and graphical abstract.
Publication style, 300 dpi PNG + LZW TIFF."""
import os
import re
import textwrap
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe
from matplotlib_venn import venn2
from lifelines import KaplanMeierFitter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.5,
    "axes.titlesize": 8.5, "axes.labelsize": 7.5,
    "xtick.labelsize": 6.8, "ytick.labelsize": 6.8,
    "axes.linewidth": 0.8, "axes.edgecolor": "#333333",
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "legend.fontsize": 6.8, "legend.frameon": False,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
})

C_PDAC, C_CP, C_SH = "#8E44AD", "#DF8A34", "#2E6DA4"
C_MIC, C_MET, C_TGT = "#7FAF7A", "#E8A33D", "#5B7FA6"
C_AXIS = {"AHR": "#B03A2E", "STAT3": "#2E6DA4", "NFkB": "#4F9A6D"}
CLS_COL = {"indole_trp": "#7E57C2", "bile_acid": "#C08A3E", "scfa": "#4F9A6D",
           "polyamine": "#C05A7E", "other": "#999999"}


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=300)
    fig.savefig(os.path.join(FIG, name + ".tif"), dpi=300,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("saved", name)


def ttl(ax, letter, text, pad=0.014):
    bb = ax.get_subplotspec().get_position(ax.figure)
    ax.figure.text(bb.x0, bb.y1 + pad, r"$\mathbf{%s}$   %s" % (letter, text),
                   fontsize=7.8, va="bottom", ha="left")


core = pd.read_csv(os.path.join(RES, "core_targets.csv"))["core_target"].tolist()
union = pd.read_csv(os.path.join(RES, "disease_genes_union.csv"))
shared = pd.read_csv(os.path.join(RES, "pdac_cp_shared_genes.csv"))
ev = pd.read_csv(os.path.join(RES, "core_target_evidence.csv"))
hub = pd.read_csv(os.path.join(RES, "hub_genes.csv"))["hub_gene"].tolist()
votes = pd.read_csv(os.path.join(RES, "hub_votes.csv"))
axm = pd.read_csv(os.path.join(RES, "axis_membership.csv"))
axe = pd.read_csv(os.path.join(RES, "axis_engagement.csv"))
enr = pd.read_csv(os.path.join(RES, "enrichment_gprofiler.csv"))
ppi = pd.read_csv(os.path.join(RES, "ppi_string_edges.csv"))
adme = pd.read_csv(os.path.join(RES, "adme_lipinski.csv"))
dock = pd.read_csv(os.path.join(RES, "docking_matrix.csv"))
panel = pd.read_csv(os.path.join(RES, "focus_metabolite_panel.csv"))
tpm = pd.read_csv(os.path.join(RES, "targets_per_metabolite.csv"))
reach = pd.read_csv(os.path.join(RES, "gmt_microbe_reach.csv"))
gnodes = pd.read_csv(os.path.join(RES, "gmt_network_nodes.csv"))

n_pdac = int(union["in_PDAC"].sum()); n_cp = int(union["in_CP"].sum())
n_both = int((union["in_PDAC"] & union["in_CP"]).sum())

# =====================================================================
# Figure 1 - study design + disease gene space
# =====================================================================
fig = plt.figure(figsize=(7.09, 4.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 0.95], hspace=0.52, wspace=0.55,
                      left=0.06, right=0.98, top=0.86, bottom=0.09)

ax = fig.add_subplot(gs[0, :])
ax.set_xlim(0, 100); ax.set_ylim(0, 30); ax.axis("off")
boxes = [
    (1, 15.5, 22, 12, "#EAE6F2", "Disease gene space",
     f"Open Targets platform\nPDAC (MONDO_0005184): {n_pdac} genes\nCP (MONDO_0005003): {n_cp} genes\nassociation score >= 0.2"),
    (27, 15.5, 22, 12, "#EAF1EA", "Gut microbial metabolites",
     "gutMGene v2.0 catalogue\n69 focus-class metabolites\n(indole/Trp, bile acids,\nSCFAs, polyamines)"),
    (53, 15.5, 22, 12, "#E8F0F7", "Metabolite-target mapping",
     "gutMGene curated records\n562 target genes\nChEMBL Tanimoto >= 0.85\n59 target genes"),
    (79, 15.5, 20, 12, "#FDF2E3", "Integration",
     "614 metabolite targets\n578 disease genes\n23 core targets"),
    (14, 1, 31, 9, "#F4F4F4", "Systems analysis",
     "GO / KEGG / Reactome enrichment\nSTRING PPI, hub genes (TP53, JUN, HIF1A)\nAHR / STAT3 / NF-kB axis engagement"),
    (55, 1, 42, 9, "#F4F4F4", "Validation",
     "TCGA-PAAD survival + tumour microenvironment\nGEO PDAC / CP cohorts\nAutoDock Vina docking + ADME"),
]
for x, y, w, h, fc, head, body in boxes:
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4",
                                fc=fc, ec="#555555", lw=0.7))
    ax.text(x + w / 2, y + h - 1.5, head, ha="center", va="top",
            fontsize=6.6, fontweight="bold")
    ax.text(x + w / 2, y + h - 3.8, body, ha="center", va="top", fontsize=5.7)
# horizontal arrows between top boxes
for x0, x1, y in ((23.5, 26.5, 21.5), (49.5, 52.5, 21.5), (75.5, 78.5, 21.5)):
    ax.annotate("", xy=(x1, y), xytext=(x0, y),
                arrowprops=dict(arrowstyle="-|>", color="#333333", lw=1.0))
# Integration feeds both Systems analysis and Validation (straight lines)
# Endpoints at y=10.8 are above the boxes' padded top edge (y=10.4);
# shrinkB pulls the arrowhead further back so it never touches the box.
ax.annotate("", xy=(33.5, 11.2), xytext=(80.5, 14.1),
            arrowprops=dict(arrowstyle="-|>", color="#333333", lw=1.0,
                            shrinkA=3, shrinkB=4))
ax.annotate("", xy=(74, 11.2), xytext=(86.5, 14.1),
            arrowprops=dict(arrowstyle="-|>", color="#333333", lw=1.0,
                            shrinkA=3, shrinkB=4))
ttl(ax, "a", "Study design")

ax = fig.add_subplot(gs[1, 0])
v = venn2(subsets=(n_pdac - n_both, n_cp - n_both, n_both),
          set_labels=(f"PDAC\n({n_pdac})", f"CP\n({n_cp})"), ax=ax)
for t, c in zip(v.set_labels, (C_PDAC, C_CP)):
    if t: t.set_fontsize(7)
for t in v.subset_labels:
    if t: t.set_fontsize(7.5)
if v.get_patch_by_id("10"): v.get_patch_by_id("10").set_color(C_PDAC); v.get_patch_by_id("10").set_alpha(0.5)
if v.get_patch_by_id("01"): v.get_patch_by_id("01").set_color(C_CP); v.get_patch_by_id("01").set_alpha(0.55)
if v.get_patch_by_id("11"): v.get_patch_by_id("11").set_color("#6B4E8E"); v.get_patch_by_id("11").set_alpha(0.8)
ttl(ax, "b", "PDAC-CP shared genes")

ax = fig.add_subplot(gs[1, 1])
x = np.arange(len(shared)); w = 0.38
ax.bar(x - w / 2, shared["pdac_score"], w, color=C_PDAC, label="PDAC score")
ax.bar(x + w / 2, shared["cp_score"], w, color=C_CP, label="CP score")
ax.set_xticks(x); ax.set_xticklabels(shared["symbol"], fontsize=7, fontweight="bold")
ax.set_ylabel("Open Targets association score")
ax.set_ylim(0, 0.68); ax.legend(loc="upper right", fontsize=6)
ttl(ax, "c", "Transition triad (PDAC & CP)")

ax = fig.add_subplot(gs[1, 2])
evc = ev.groupby("gene")["source"].nunique()
both = sorted(evc[evc == 2].index)
rows = []
for g in core:
    rows.append({"gene": g, "ChEMBL": "ChEMBL_similarity" in set(ev.loc[ev["gene"] == g, "source"]),
                 "gutMGene": "gutMGene_curated" in set(ev.loc[ev["gene"] == g, "source"])})
evm = pd.DataFrame(rows).set_index("gene")
for xoff, genes_blk in ((0.0, core[:12]), (2.55, core[12:])):
    y = np.arange(len(genes_blk))[::-1]
    for yi, g in zip(y, genes_blk):
        for xi, col in enumerate(("gutMGene", "ChEMBL")):
            on = evm.loc[g, col]
            ax.scatter(xoff + xi * 0.62, yi, s=46, marker="o",
                       facecolor="#2E6DA4" if on else "white",
                       edgecolor="#2E6DA4" if on else "#BBBBBB", lw=1.0, zorder=3)
        ax.text(xoff - 0.28, yi, g, ha="right", va="center", fontsize=5.8,
                fontweight="bold" if g in hub else "normal")
    ax.text(xoff, 12.0, "gm", ha="center", fontsize=5.4, style="italic")
    ax.text(xoff + 0.62, 12.0, "ch", ha="center", fontsize=5.4, style="italic")
ax.set_xticks([])
ax.set_xlim(-1.45, 3.6); ax.set_ylim(-1.9, 12.6)
ax.set_yticks([])
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.text(1.1, -1.55, f"gm = gutMGene curated, ch = ChEMBL similarity; dual: {', '.join(both)}",
        ha="center", fontsize=5.6, style="italic")
ttl(ax, "d", "Core-target evidence (bold = hub)")
save(fig, "Fig1_study_design_and_gene_overlap")

# =====================================================================
# Figure 2 - metabolite-target integration + microbiota network
# (merged from the former Fig 2 and Fig 5)
# =====================================================================
fig = plt.figure(figsize=(7.09, 7.6))
gs = fig.add_gridspec(3, 6, height_ratios=[1.0, 1.7, 0.95], hspace=0.55, wspace=1.9,
                 left=0.10, right=0.97, top=0.945, bottom=0.075)

ax = fig.add_subplot(gs[0, 0:2])
cc = panel["class"].value_counts()
ax.barh(range(len(cc)), cc.values, color=[CLS_COL.get(c, "#999") for c in cc.index])
ax.set_yticks(range(len(cc)))
ax.set_yticklabels([{"indole_trp": "Indole / Trp", "bile_acid": "Bile acids",
                     "scfa": "SCFAs", "polyamine": "Polyamines"}.get(c, c) for c in cc.index],
                   fontsize=7)
for i, v in enumerate(cc.values):
    ax.text(v + 0.4, i, str(v), va="center", fontsize=7)
ax.invert_yaxis()
ax.set_xlim(0, cc.max() * 1.18)
ax.set_xlabel("N metabolites (focus panel)")
ttl(ax, "a", f"Metabolite focus classes (n = {len(panel)})")

ax = fig.add_subplot(gs[0, 2:4])
gm_genes = set(ev.loc[ev["source"] == "gutMGene_curated", "gene"])
ch_genes = set(ev.loc[ev["source"] == "ChEMBL_similarity", "gene"])
v = venn2(subsets=(len(ch_genes - gm_genes), len(gm_genes - ch_genes), len(ch_genes & gm_genes)),
          set_labels=None, ax=ax)
for t in v.subset_labels:
    if t: t.set_fontsize(6.8)
if v.get_patch_by_id("10"): v.get_patch_by_id("10").set_color("#4F86A6"); v.get_patch_by_id("10").set_alpha(0.5)
if v.get_patch_by_id("01"): v.get_patch_by_id("01").set_color("#7FAF7A"); v.get_patch_by_id("01").set_alpha(0.5)
# set labels are stacked inside the axes so that they never spill into panel c
ax.text(0.0, -0.08, f"ChEMBL similarity ({len(ch_genes)})", transform=ax.transAxes,
        ha="left", va="top", fontsize=6.0, color="#2F5F7A")
ax.text(0.0, -0.21, f"gutMGene curated ({len(gm_genes)})", transform=ax.transAxes,
        ha="left", va="top", fontsize=6.0, color="#4A7C46")
ax.text(0.5, -0.44, "union of 614 metabolite target genes;\n23 overlap the PDAC/CP disease space",
        transform=ax.transAxes, ha="center", fontsize=5.8, style="italic")
ttl(ax, "b", "Target evidence sources")

ax = fig.add_subplot(gs[0, 4:6])
t = tpm[tpm["metabolite"].isin(set(panel["name"]))].sort_values("n_targets", ascending=False).head(12).iloc[::-1]
lab = t["metabolite"].str.replace(" acid", "", regex=False).str.replace(" Acid", "", regex=False)
ax.barh(range(len(t)), t["n_targets"], color="#C08A3E")
ax.set_yticks(range(len(t))); ax.set_yticklabels(lab, fontsize=5.6)
for i, v0 in enumerate(t["n_targets"]):
    ax.text(v0 + 1.5, i, str(v0), va="center", fontsize=6.3)
ax.set_xlim(0, t["n_targets"].max() * 1.20)
ax.set_xlabel("N target genes (gutMGene)")
ttl(ax, "c", "Focus-panel target breadth")

# --- 面板 d：微生物-代谢物-靶点网络（原 Figure 5 上部） ---
ax = fig.add_subplot(gs[1, :])
gm_edges = pd.read_csv(os.path.join(RES, "gmt_gm_edges.csv"))
mt_edges = pd.read_csv(os.path.join(RES, "gmt_mt_edges.csv"))
top_mic = reach.nlargest(14, "n_core_targets")["microbe"].tolist()
met_core = mt_edges.groupby("metabolite")["gene"].nunique().sort_values(ascending=False)
top_met = met_core.head(12).index.tolist()
tgt = sorted(set(mt_edges[mt_edges["metabolite"].isin(top_met)]["gene"]))
B = nx.Graph()
for m in top_mic: B.add_node(m, layer=0)
for m in top_met: B.add_node(m, layer=1)
for g in tgt: B.add_node(g, layer=2)
for _, r in gm_edges.iterrows():
    if r["microbe"] in top_mic and r["metabolite"] in top_met:
        B.add_edge(r["microbe"], r["metabolite"])
for _, r in mt_edges.iterrows():
    if r["metabolite"] in top_met and r["gene"] in tgt:
        B.add_edge(r["metabolite"], r["gene"])
layer_x = {0: 0.0, 1: 1.0, 2: 2.0}
spacing = {0: 1.0, 1: 1.45, 2: 1.0}
pos = {}
for L in (0, 1, 2):
    nodes = [n for n in B.nodes if B.nodes[n]["layer"] == L]
    nodes.sort(key=lambda n: -B.degree(n))
    sp = spacing[L]
    for i, n in enumerate(nodes):
        pos[n] = (layer_x[L], -i * sp + (len(nodes) - 1) * sp / 2.0)
def _nkey(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


# chemical class is assigned to every network metabolite, not only to the 69-panel
met_class = {_nkey(k): v for k, v in zip(panel["name"], panel["class"])}
# metabolites outside the 69-metabolite focus panel, classified by chemical family
_ext_class = {"Indoxyl sulfate": "indole_trp", "Succinate": "scfa",
              "Phenylacetylglutamine": "other", "p-Cresol sulfate": "other",
              "Genipin": "other", "Hydroquinone": "other",
              "10-Oxo-11-octadecenoic acid": "other",
              "10-keto-12Z-octadecenoic acid": "other"}
met_class.update({_nkey(k): v for k, v in _ext_class.items()})


def met_cls(name):
    k = _nkey(name)
    if k not in met_class:
        print("WARNING: unclassified metabolite ->", name)
    return met_class.get(k, "other")


for e0, e1 in B.edges():
    l0 = B.nodes[e0]["layer"]
    ax.plot([pos[e0][0], pos[e1][0]], [pos[e0][1], pos[e1][1]],
            color="#CCCCCC" if l0 == 0 else "#D8B96A", lw=0.5, zorder=1, alpha=0.7)
for n in B.nodes:
    L = B.nodes[n]["layer"]
    if L == 0:
        c, s = C_MIC, 30 + 16 * B.degree(n)
    elif L == 1:
        c = CLS_COL[met_cls(n)]; s = 40 + 22 * B.degree(n)
    else:
        c, s = C_TGT, 34 + 10 * B.degree(n)
    ax.scatter(*pos[n], s=s, c=c, edgecolor="white", lw=0.5, zorder=3)
    if L == 1:
        # metabolite labels: larger offset so they sit outside the big circles
        off = (16, 0); ha = "left"
    else:
        off = (-6, 0) if L == 0 else (6, 0)
        ha = "right" if L == 0 else "left"
    fw = "bold" if (L == 2 and n in hub) or L == 1 and B.degree(n) >= 4 else "normal"
    fs = 5.4 if L != 2 else 5.8
    ann = ax.annotate(n, pos[n], xytext=off, textcoords="offset points", ha=ha,
                va="center", fontsize=fs, fontweight=fw, zorder=4,
                fontstyle="italic" if L == 0 else "normal")
ax.set_xlim(-0.75, 2.75)
ys = [p[1] for p in pos.values()]
ax.set_ylim(min(ys) - 1.2, max(ys) + 3.4)
ax.axis("off")
_cls_lab = {"indole_trp": "Indole/Trp", "bile_acid": "Bile acid", "scfa": "SCFA",
            "polyamine": "Polyamine", "other": "Other"}
_handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=C_MIC,
                   markersize=6, label="Microbe (genus/species)")]
_handles += [Line2D([0], [0], marker="o", color="none", markerfacecolor=CLS_COL[k],
                    markersize=6, label=_cls_lab[k]) for k in _cls_lab]
_handles.append(Line2D([0], [0], marker="o", color="none", markerfacecolor=C_TGT,
                       markersize=6, label="Core target"))
ax.legend(handles=_handles, loc="upper left", fontsize=5.8, ncol=2,
          columnspacing=1.4, handletextpad=0.2, borderpad=0.2)
ttl(ax, "d", "Microbiota - metabolite - core-target network (top contributors)", pad=0.01)

# --- 面板 e、f：微生物覆盖与桥接代谢物（原 Figure 5 下部） ---
ax = fig.add_subplot(gs[2, 0:3])
r = reach.nlargest(12, "n_core_targets").iloc[::-1]
ax.barh(range(len(r)), r["n_core_targets"], color=C_MIC)
ax.set_yticks(range(len(r)))
ax.set_yticklabels(r["microbe"], fontsize=6.2, fontstyle="italic")
for i, (v0, m) in enumerate(zip(r["n_core_targets"], r["n_metabolites"])):
    ax.text(v0 + 0.15, i, f"{v0} targets / {m} metab.", va="center", fontsize=5.8)
ax.set_xlabel("N core targets reached")
ax.set_xlim(0, 13.5)
ttl(ax, "e", "Top microbes by core-target reach")

ax = fig.add_subplot(gs[2, 3:6])
gn = gnodes[gnodes["layer"] == "metabolite"].nlargest(10, "degree").iloc[::-1]
ax.barh(range(len(gn)), gn["degree"],
        color=[CLS_COL[met_cls(m)] for m in gn["node"]])
ax.set_yticks(range(len(gn)))
_short = {"Taurochenodeoxycholic acid": "TCDCA", "Tauroursodeoxycholic acid": "TUDCA",
          "Isoallolithocholic acid": "IsoalloLCA", "Taurocholic acid": "TCA",
          "Taurodeoxycholic acid": "TDCA"}
ax.set_yticklabels([textwrap.fill(_short.get(x, x.replace(" acid", " a.").replace(" Acid", " a.")), 16)
                    for x in gn["node"]], fontsize=6.0)
for i, v0 in enumerate(gn["degree"]):
    ax.text(v0 + 1, i, str(v0), va="center", fontsize=6.4)
ax.set_xlabel("Degree in microbiota-metabolite-target network")
ax.set_xlim(0, gn["degree"].max() * 1.18)
handles = [Line2D([0], [0], marker="s", color="none", markerfacecolor=c, markersize=6,
                  label={"indole_trp": "Indole/Trp", "bile_acid": "Bile acid",
                         "scfa": "SCFA", "polyamine": "Polyamine", "other": "Other"}[k])
           for k, c in CLS_COL.items()]
ax.legend(handles=handles, fontsize=5.8, loc="lower right")
ttl(ax, "f", "Bridging metabolites by degree")
save(fig, "Fig2_metabolite_integration_network")

# =====================================================================
# Figure 3 - functional enrichment
# =====================================================================
fig = plt.figure(figsize=(7.09, 4.9))
gs = fig.add_gridspec(1, 3, wspace=1.1, left=0.32, right=0.92, top=0.86, bottom=0.22)
specs = [("KEGG", "KEGG pathways"), ("REAC", "Reactome pathways"), ("GO:BP", "GO biological process")]
for i, (src, label) in enumerate(specs):
    ax = fig.add_subplot(gs[0, i])
    d = enr[enr["source"] == src].nsmallest(8, "p_value").iloc[::-1]
    d["mlogp"] = -np.log10(d["p_value"])
    sc = ax.scatter(d["mlogp"], range(len(d)), s=d["intersection_size"] * 14,
               c=d["mlogp"], cmap="viridis", edgecolor="#333333", lw=0.4, zorder=3)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([textwrap.fill(x, 17) for x in d["term_name"]], fontsize=5.9)
    ax.set_xlim(0, d["mlogp"].max() * 1.18)
    ax.set_xlabel("-log10(adj. p)")
    ax.grid(axis="x", color="#E5E5E5", lw=0.5, zorder=0)
    ttl(ax, "abc"[i], label)
cbar_ax = fig.add_axes([0.935, 0.22, 0.012, 0.55])
cb = fig.colorbar(sc, cax=cbar_ax)
cb.set_label("-log10(adj. p)", fontsize=6.6)
cb.ax.tick_params(labelsize=6.0)
fig.text(0.5, 0.01, "Bubble size = number of core targets in term; colour = -log10(adjusted p)",
         ha="center", fontsize=6.4, style="italic")
save(fig, "FigS7_functional_enrichment")

# =====================================================================
# Figure 4 - PPI network, hubs, axis engagement
# =====================================================================
fig = plt.figure(figsize=(7.09, 3.8))
gs = fig.add_gridspec(1, 3, wspace=0.62, left=0.055, right=0.985, top=0.70, bottom=0.17)

ax = fig.add_subplot(gs[0, 0])
G = nx.Graph()
G.add_nodes_from(core)
for _, r in ppi.drop_duplicates(subset=["source", "target"]).iterrows():
    G.add_edge(r["source"], r["target"], weight=r["combined_score"])
_iso = [n for n in list(G.nodes) if G.degree(n) == 0]
G.remove_nodes_from(_iso)
pos = nx.spring_layout(G, seed=7, k=1.75, iterations=400)
deg = dict(G.degree())
hub_list = [h for h in hub if h in G]
L = 1.15
_bb = ax.get_position()
_W = _bb.width * fig.get_figwidth() * 72.0
_H = _bb.height * fig.get_figheight() * 72.0
A, B = _W / (2 * L), _H / (2 * L)          # points per data unit
rad = {n: np.sqrt((10 + 13 * deg[n]) / np.pi) for n in G}
nodes_l = list(G.nodes)

# normalise the spring layout, then separate the discs in screen space so that
# no two circles can touch (label placement below relies on that clearance)
xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]
span = max(max(xs) - min(xs), max(ys) - min(ys)) / 2.0
for p in pos.values():
    p[0] = (p[0] - (min(xs) + max(xs)) / 2.0) / span * 0.93
    p[1] = (p[1] - (min(ys) + max(ys)) / 2.0) / span * 0.93
for _ in range(400):
    for i, a_ in enumerate(nodes_l):
        for b_ in nodes_l[i + 1:]:
            ddx = (pos[a_][0] - pos[b_][0]) * A
            ddy = (pos[a_][1] - pos[b_][1]) * B
            d = np.hypot(ddx, ddy) + 1e-9
            need = rad[a_] + rad[b_] + 2.6
            if d < need:
                push = (need - d) * 0.5
                pos[a_][0] += ddx / d * push / A; pos[a_][1] += ddy / d * push / B
                pos[b_][0] -= ddx / d * push / A; pos[b_][1] -= ddy / d * push / B
cx = np.mean([p[0] for p in pos.values()]); cy = np.mean([p[1] for p in pos.values()])

ws = [e[2]["weight"] for e in G.edges(data=True)]
wmin, wmax = min(ws), max(ws)
for e in G.edges(data=True):
    t = (e[2]["weight"] - wmin) / max(wmax - wmin, 1e-9)
    ax.plot([pos[e[0]][0], pos[e[1]][0]], [pos[e[0]][1], pos[e[1]][1]],
            color="#9A9A9A", lw=0.45 + 1.55 * t, alpha=0.55,
            solid_capstyle="round", zorder=1)

for n in sorted(G.nodes, key=lambda x: (x in hub_list, -deg[x])):
    c = "#B03A2E" if n in hub_list else "#5B7FA6"
    ax.scatter(*pos[n], s=10 + 13 * deg[n], c=c, edgecolor="white", lw=0.6,
               zorder=4 if n in hub_list else 3)
ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.axis("off")
_texts = []

# ---- labels: anchored on the node rim, pushed outward, never overlapping ---
node_pt = {n: ((pos[n][0] + L) / (2 * L) * _W, (pos[n][1] + L) / (2 * L) * _H) for n in G}
placed = []


def _box(anchor, ux, uy, text, fs, bold):
    w = len(text) * (0.44 if bold else 0.42) * fs
    h = 1.30 * fs
    ha = "left" if ux > 0.25 else "right" if ux < -0.25 else "center"
    va = "bottom" if uy > 0.25 else "top" if uy < -0.25 else "center"
    x0 = anchor[0] if ha == "left" else anchor[0] - w if ha == "right" else anchor[0] - w / 2
    y0 = anchor[1] if va == "bottom" else anchor[1] - h if va == "top" else anchor[1] - h / 2
    return (x0, y0, x0 + w, y0 + h), ha, va


def _cost(b, n):
    x0, y0, x1, y1 = b
    c = 0.0
    for ob in placed:
        dx = min(x1, ob[2]) - max(x0, ob[0]); dy = min(y1, ob[3]) - max(y0, ob[1])
        if dx > 0 and dy > 0:
            c += dx * dy + 60.0
    for m, (mx, my) in node_pt.items():
        if m == n:
            continue
        ddx = max(x0 - mx, 0.0, mx - x1); ddy = max(y0 - my, 0.0, my - y1)
        if np.hypot(ddx, ddy) < rad[m] + 0.9:
            c += 90.0
    return c


for n in sorted(G.nodes, key=lambda x: (x not in hub_list, -deg[x])):
    px, py = node_pt[n]
    a0 = np.arctan2(pos[n][1] - cy, pos[n][0] - cx)
    bold = n in hub_list
    fs = 5.8 if bold else 5.2
    best = None
    # stay on the outward hemisphere and as close to the rim as collisions allow
    for extra in (0.0, 2.5, 5.0, 7.5):
        for da in (0, 22, -22, 45, -45, 68, -68, 90, -90):
            a = a0 + np.radians(da)
            ux, uy = float(np.cos(a)), float(np.sin(a))
            gap = rad[n] + 2.2 + extra
            b, ha, va = _box((px + ux * gap, py + uy * gap), ux, uy, n, fs, bold)
            c = _cost(b, n) + 0.55 * abs(da) + 3.0 * extra
            if best is None or c < best[0]:
                best = (c, b, ux, uy, ha, va, gap)
    _, b, ux, uy, ha, va, gap = best
    placed.append(b)
    _t = ax.annotate(n, pos[n], fontsize=fs, ha=ha, va=va,
                     fontweight="bold" if bold else "normal",
                     color="#8B1A1A" if bold else "#222222",
                     xytext=(ux * gap, uy * gap), textcoords="offset points",
                     zorder=6 if bold else 5,
                     arrowprops=None if gap - rad[n] <= 3.2 else
                     dict(arrowstyle="-", lw=0.4, color="#B4B4B4",
                          shrinkA=0.0, shrinkB=1.2))
    _texts.append(_t)
_ms = lambda d: 2.0 * np.sqrt((10 + 13 * d) / np.pi) * 0.85
ax.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor="#B03A2E",
                          markeredgecolor="white", markersize=_ms(max(deg.values())),
                          label="Consensus hub"),
                   Line2D([0], [0], marker="o", color="none", markerfacecolor="#5B7FA6",
                          markeredgecolor="white", markersize=_ms(int(np.median(list(deg.values())))),
                          label="Other core target")],
          loc="upper center", fontsize=6, handletextpad=0.2, borderpad=0.1,
          bbox_to_anchor=(0.5, -0.06), frameon=False)
ax.text(0.5, -0.235, "Node area scales with STRING degree",
        transform=ax.transAxes, ha="center", fontsize=5.2, style="italic", color="#555555")
ax.text(0.5, -0.272, "ASPM: no edge at combined score >= 0.9 (not drawn)",
        transform=ax.transAxes, ha="center", fontsize=5.2, style="italic", color="#555555")
ttl(ax, "a", "Core-target PPI network\n(%d connected nodes, %d edges)"
    % (G.number_of_nodes(), G.number_of_edges()))

ax = fig.add_subplot(gs[0, 1])
v2 = votes.sort_values("methods_top5")
ax.barh(range(len(v2)), v2["methods_top5"],
        color=["#B03A2E" if g in hub else "#5B7FA6" for g in v2["gene"]])
ax.set_yticks(range(len(v2))); ax.set_yticklabels(v2["gene"], fontsize=6.8)
ax.set_xlabel("N cytoHubba methods (of 6) in top 5")
ax.set_xlim(0, 6.6)
for i, vv in enumerate(v2["methods_top5"]):
    ax.text(vv + 0.2, i, str(vv), va="center", fontsize=6.8)
ttl(ax, "b", "Hub consensus (cytoHubba)")

ax = fig.add_subplot(gs[0, 2])
top_ax = axe.sort_values("weight_sum", ascending=False).groupby("axis").head(5)
axes_order = ["AHR", "STAT3", "NFkB"]
genes_ax = []
for axx in axes_order:
    genes_ax += top_ax[top_ax["axis"] == axx].sort_values("weight_sum", ascending=False)["core_target"].tolist()
genes_ax = list(dict.fromkeys(genes_ax))[:12][::-1]
bottoms = np.zeros(len(genes_ax))
for axx in axes_order:
    vals = [axe[(axe["axis"] == axx) & (axe["core_target"] == g)]["weight_sum"].sum() for g in genes_ax]
    ax.barh(range(len(genes_ax)), vals, left=bottoms, color=C_AXIS[axx], label=axx + " axis", height=0.72)
    bottoms += np.array(vals)
ax.set_yticks(range(len(genes_ax))); ax.set_yticklabels(genes_ax, fontsize=6.6)
ax.set_xlabel("Cumulative STRING edge weight to axis genes")
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=6.2)
ttl(ax, "c", "AHR / STAT3 / NF-kB axis engagement")
save(fig, "Fig3_ppi_network_hubs")

# =====================================================================
# Figure 6 - TCGA-PAAD validation
# =====================================================================
surv = pd.read_csv(os.path.join(DATA, "TCGA", "TCGA_survival_pancanatlas.tsv"), sep="\t")
surv = surv[surv["cancer type abbreviation"] == "PAAD"]
ms = pd.read_csv(os.path.join(RES, "tcga_module_scores.csv"))
df = ms.merge(surv[["sample", "OS", "OS.time", "PFI", "PFI.time"]], on="sample", how="inner")
cox_os = pd.read_csv(os.path.join(RES, "tcga_cox_os.csv")).iloc[0]
cox_pfi = pd.read_csv(os.path.join(RES, "tcga_cox_pfi.csv")).iloc[0]
gcox = pd.read_csv(os.path.join(RES, "tcga_gene_cox_os.csv"))
ic = pd.read_csv(os.path.join(RES, "tcga_immune_corr.csv"))
ac = pd.read_csv(os.path.join(RES, "tcga_axis_corr.csv"))

fig = plt.figure(figsize=(7.09, 5.6))
gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.0, 0.85], hspace=0.78, wspace=0.85,
                      left=0.075, right=0.98, top=0.945, bottom=0.075)

for i, (tcol, ecol, label, cox) in enumerate(
        (("OS.time", "OS", "Overall survival", cox_os),
         ("PFI.time", "PFI", "Progression-free interval", cox_pfi))):
    ax = fig.add_subplot(gs[0, i])
    d = df.dropna(subset=[tcol, ecol]); d = d[d[tcol] > 0]
    hi = d["core_module"] >= d["core_module"].median()
    km = KaplanMeierFitter()
    for mask, lab, c in ((hi, f"High module (n={int(hi.sum())})", "#B03A2E"),
                         (~hi, f"Low module (n={int((~hi).sum())})", "#2E6DA4")):
        km.fit(d.loc[mask, tcol] / 30.44, d.loc[mask, ecol], label=lab)
        km.plot_survival_function(ax=ax, color=c, ci_show=True, ci_alpha=0.15, lw=1.3)
    p = cox["OS_median_logrank_p" if label.startswith("Overall") else "PFI_median_logrank_p"]
    ax.text(0.03, 0.06, f"log-rank p = {p:.4f}", transform=ax.transAxes, fontsize=7)
    ax.set_xlabel("Months"); ax.set_ylabel(label)
    ax.set_xlim(0, 120)
    ttl(ax, "ab"[i], f"{label} by core module")
    ax.legend(loc="upper right", fontsize=6.4)

ax = fig.add_subplot(gs[0, 2])
rows = [("OS, per SD", cox_os["OS_cont_HR"], cox_os["OS_cont_CI"]),
        ("OS, adjusted*", cox_os["OS_adj_HR"], cox_os["OS_adj_CI"]),
        ("PFI, per SD", cox_pfi["PFI_cont_HR"], cox_pfi["PFI_cont_CI"]),
        ("PFI, adjusted*", cox_pfi["PFI_adj_HR"], cox_pfi["PFI_adj_CI"])]
ylab = [r[0] for r in rows][::-1]
hrs = [r[1] for r in rows][::-1]
cis = [r[2].split("-") for r in rows][::-1]
lo = [float(c[0]) for c in cis]; hi_ = [float(c[1]) for c in cis]
ax.errorbar(hrs, range(len(rows)), xerr=[np.array(hrs) - lo, np.array(hi_) - hrs],
            fmt="o", color="#8E44AD", ecolor="#8E44AD", elinewidth=1.2, capsize=2.5, ms=4.5)
ax.axvline(1, color="#888888", lw=0.8, ls="--")
ax.set_yticks(range(len(rows))); ax.set_yticklabels(ylab, fontsize=6.8)
ax.set_xlabel("Cox hazard ratio (95% CI)")
ax.set_xlim(0.9, 3.0)
ax.set_ylim(-0.55, 3.9)
for i, h in enumerate(hrs):
    ax.text(h, i + 0.30, f"{h:.2f}", ha="center", fontsize=6.4)
ttl(ax, "c", "Module prognostic effect (*adj. age + stage)")

ax = fig.add_subplot(gs[1, 0])
g = gcox.nsmallest(10, "p").iloc[::-1]
ax.errorbar(g["HR"], range(len(g)),
            xerr=[g["HR"] - g["CI_low"], g["CI_high"] - g["HR"]],
            fmt="o", ms=4, color="#B03A2E", ecolor="#B03A2E", elinewidth=1.0, capsize=2)
ax.axvline(1, color="#888888", lw=0.8, ls="--")
ax.set_yticks(range(len(g)))
ax.set_yticklabels([f"{r.gene}  (q={r.q_bh:.2g})" for r in g.itertuples()], fontsize=6.2)
ax.set_xlabel("Univariate Cox OS hazard ratio")
ttl(ax, "d", "Per-gene prognostic effect (top 10)")

ax = fig.add_subplot(gs[1, 1])
ic2 = ic.sort_values("spearman_r")
col = ["#B03A2E" if q < 0.05 else "#BBBBBB" for q in ic2["q_bh"]]
ax.barh(range(len(ic2)), ic2["spearman_r"], color=col)
ax.set_yticks(range(len(ic2)))
ax.set_yticklabels(ic2["signature"], fontsize=5.6)
ax.axvline(0, color="#888888", lw=0.8)
ax.set_xlabel("Spearman r (core module vs signature)")
ttl(ax, "e", "Tumour-microenvironment\ncorrelations")
ax.text(0.97, 0.03, "red: BH q < 0.05", transform=ax.transAxes, fontsize=6,
        color="#B03A2E", ha="right")

ax = fig.add_subplot(gs[1, 2])
ac2 = ac.set_index("axis").loc[["AHR", "STAT3", "NFkB"]].reset_index()
ax.bar(range(3), ac2["spearman_r"], color=[C_AXIS[a] for a in ac2["axis"]], width=0.62)
ax.set_xticks(range(3)); ax.set_xticklabels(ac2["axis"], fontsize=7.2)
for i, (r0, p0) in enumerate(zip(ac2["spearman_r"], ac2["p"])):
    dx = -0.25 if ac2["axis"][i] == "AHR" else 0.25 if ac2["axis"][i] == "NFkB" else 0.0
    ytxt = max(r0, 0.0) + 0.02
    xoff_r = 9.0 if ac2["axis"][i] == "AHR" else -10.4 if ac2["axis"][i] == "NFkB" else 0
    xoff_p = 11.0 if ac2["axis"][i] == "AHR" else -10.2 if ac2["axis"][i] == "NFkB" else 0
    _trans_r = matplotlib.transforms.offset_copy(ax.transData, fig, x=xoff_r, units="points")
    _trans_p = matplotlib.transforms.offset_copy(ax.transData, fig, x=xoff_p, units="points")
    ax.text(i + dx, ytxt, f"r={r0:.2f}\n", ha="center", va="bottom", fontsize=5.6, transform=_trans_r)
    ax.text(i + dx, ytxt, f"p={p0:.2g}", ha="center", va="bottom", fontsize=5.6, transform=_trans_p)
ax.set_ylabel("Spearman r with core module", fontsize=6.5)
ax.set_ylim(0, 0.75)
ttl(ax, "f", "Axis module co-expression")

ax = fig.add_subplot(gs[2, :])
expr_dir = pd.read_csv(os.path.join(RES, "tcga_core_gene_expression.csv"))
expr_dir = expr_dir.sort_values("log2FC_tumor_vs_normal", ascending=False)
x = np.arange(len(expr_dir))
col = ["#B03A2E" if (l > 0 and q < 0.05) else "#2E6DA4" if (l < 0 and q < 0.05) else "#888888"
       for l, q in zip(expr_dir["log2FC_tumor_vs_normal"], expr_dir["q_bh"])]
ax.bar(x, expr_dir["log2FC_tumor_vs_normal"], color=col, width=0.75)
ax.set_xticks(x); ax.set_xticklabels(expr_dir["gene"], rotation=45, ha="right", fontsize=5.8)
ax.axhline(0, color="#888888", lw=0.8)
ax.set_ylabel("log2 ratio, tumour / normal\n(TCGA-PAAD, 178 vs 4)")
_gtrans = matplotlib.transforms.offset_copy(ax.transAxes, fig, y=38, units="points")
ax.text(0.995, 0.12, "grey: none reaches BH q < 0.05 (only 4 normal samples)",
        transform=_gtrans, fontsize=5.8, style="italic", ha="right")
ttl(ax, "g", "Core-target expression direction in TCGA-PAAD")
save(fig, "Fig4_tcga_validation")

# =====================================================================
# Figure 8 - docking + ADME
# =====================================================================
fig = plt.figure(figsize=(7.09, 4.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 0.9], hspace=0.5, wspace=0.6,
                      left=0.16, right=0.985, top=0.93, bottom=0.09)

ax = fig.add_subplot(gs[0, :])
d = dock.set_index("metabolite").loc[[m for m in adme["name"] if m in set(dock["metabolite"])]]
d = d.loc[d.min(axis=1).sort_values().index]
im = ax.imshow(d.values, aspect="auto", cmap="RdYlGn_r", vmin=-9, vmax=-2)
ax.set_xticks(range(d.shape[1])); ax.set_xticklabels(d.columns, fontsize=7.2)
ax.set_yticks(range(d.shape[0])); ax.set_yticklabels(d.index, fontsize=6.4)
cmap_obj = plt.get_cmap("RdYlGn_r")
for i in range(d.shape[0]):
    for j in range(d.shape[1]):
        v0 = d.values[i, j]
        rgba = cmap_obj((v0 + 9) / 7.0)
        lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
        ax.text(j, i, f"{v0:.1f}", ha="center", va="center", fontsize=5.6,
                color="white" if lum < 0.45 else "#333333")
cb = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.015)
cb.set_label("Vina affinity (kcal/mol)", fontsize=6.8)
cb.ax.tick_params(labelsize=6.4)
ttl(ax, "a", "Docking affinity heatmap (19 metabolites x AHR / STAT3)", pad=0.02)

ax = fig.add_subplot(gs[1, 0])
for cls, g in adme.groupby("class"):
    ax.scatter(g["MW"], g["cLogP"], s=14 + 40 * g["QED"],
               color=CLS_COL.get(cls, "#999"), label={"indole_trp": "Indole/Trp", "bile_acid": "Bile acid",
               "scfa": "SCFA", "polyamine": "Polyamine"}.get(cls, cls),
               edgecolor="#333333", lw=0.4, alpha=0.85, zorder=3)
ax.axvline(500, color="#BBBBBB", ls="--", lw=0.7)
ax.axhline(5, color="#BBBBBB", ls="--", lw=0.7)
ax.set_xlabel("Molecular weight (Da)"); ax.set_ylabel("cLogP")
ax.legend(fontsize=5.8, loc="upper left", handletextpad=0.1, borderpad=0.2, labelspacing=0.25)
ttl(ax, "b", "Physicochemical space (size = QED)")

ax = fig.add_subplot(gs[1, 1])
dm = dock.melt(id_vars="metabolite", var_name="receptor", value_name="aff")
best = dm.loc[dm.groupby("metabolite")["aff"].idxmin()].sort_values("aff").head(8).iloc[::-1]
ax.barh(range(len(best)), -best["aff"],
        color=["#B03A2E" if r == "AHR" else "#2E6DA4" for r in best["receptor"]])
ax.set_yticks(range(len(best)))
ax.set_yticklabels(best["metabolite"], fontsize=6.4)
for i, v0 in enumerate(-best["aff"]):
    ax.text(v0 + 0.08, i, f"{v0:.1f}", va="center", fontsize=6.4)
ax.set_xlabel("Best Vina |affinity| (kcal/mol)")
ax.set_xlim(0, (-best["aff"]).max() * 1.12)
ttl(ax, "c", "Top metabolite-receptor pairs")
save(fig, "FigS6_docking_adme")

print("ALL FIGURES DONE")
