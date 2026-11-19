"""Step 10b: GEO validation figures (Fig7 main, FigS1 supplementary QC).

Inputs (from 07c_geo_analysis.py):
  geo_module_scores.csv      gsm, module, group, dataset, contrast
  geo_core_validation.csv    gene, log2FC, p, q_bh, dataset, contrast, tcga_log2FC,
                             same_direction_as_tcga, module_p, n_a, n_b
  geo_fisher.csv             dataset, contrast, core_in_top500, core_total,
                             odds_ratio, p_fisher
  geo_samples_{GSE}.csv      gsm, group, text
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
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

C_GRP = {"PDAC": "#8E44AD", "CP": "#DF8A34", "normal": "#4F9A6D", "other": "#BBBBBB"}


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


ms = pd.read_csv(os.path.join(RES, "geo_module_scores.csv"))
cv = pd.read_csv(os.path.join(RES, "geo_core_validation.csv"))
fi = pd.read_csv(os.path.join(RES, "geo_fisher.csv"))

# per-contrast module p (first) and group means
mod_p = cv.groupby("contrast")["module_p"].first().to_dict()

# =====================================================================
# Figure 7 - GEO external validation
# =====================================================================
contrasts = cv["contrast"].unique().tolist()
# order: PDAC vs normal first, then CP contrasts
contrasts.sort(key=lambda c: ("CP" in c.split("_vs_")[0], c))

fig = plt.figure(figsize=(7.09, 6.7))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25], hspace=0.62, wspace=0.42,
                      left=0.09, right=0.985, top=0.945, bottom=0.09,
                      width_ratios=[1.0, 1.15])

ax = fig.add_subplot(gs[0, :])
mm = ms.drop_duplicates(subset=["gsm", "dataset", "contrast"])
xpos, xt, xl = [], [], []
tops = []
for i, tag in enumerate(contrasts):
    sub = mm[mm["contrast"] == tag]
    groups = [g for g in ("normal", "CP", "PDAC") if g in set(sub["group"])]
    rng = np.random.default_rng(7)
    top = -np.inf
    for j, g in enumerate(groups):
        v = sub.loc[sub["group"] == g, "module"].values
        x = i * 1.35 + j * 0.34
        bp = ax.boxplot([v], positions=[x], widths=0.28, showfliers=False,
                        patch_artist=True, medianprops={"color": "#222222", "lw": 1.0},
                        whiskerprops={"lw": 0.8}, capprops={"lw": 0.8})
        for patch in bp["boxes"]:
            patch.set_facecolor(C_GRP.get(g, "#BBBBBB"))
            patch.set_alpha(0.75)
        ax.scatter(x + rng.normal(0, 0.045, len(v)), v, s=3, color="#333333",
                   alpha=0.5, zorder=3, linewidths=0)
        top = max(top, np.nanmax(v))
    tops.append(top)
    xt.append(i * 1.35 + 0.17)
    ds, rest = tag.split("_", 1)
    xl.append(ds + "\n" + rest.replace("_vs_", " vs "))
for i, tag in enumerate(contrasts):
    p = mod_p.get(tag, np.nan)
    ax.text(xt[i], tops[i] + 0.06, f"p={p:.2g}", ha="center", fontsize=5.8)
ax.set_ylim(top=max(tops) + 0.42)
ax.set_xticks(xt)
ax.set_xticklabels(xl, fontsize=5.6)
ax.axhline(0, color="#888888", lw=0.6, ls=":")
ax.set_ylabel("Core-module score (mean z of 23 genes)")
handles = [plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=c,
                      markersize=6, label=g) for g, c in C_GRP.items() if g != "other"]
ax.legend(handles=handles, fontsize=6, loc="lower right", ncol=3)
ttl(ax, "a", "Core-module score by diagnostic group")

ax = fig.add_subplot(gs[1, 0])
fi2 = fi.set_index("contrast").loc[[c for c in contrasts if c in set(fi["contrast"])]].reset_index()
y = np.arange(len(fi2))[::-1]
col = ["#B03A2E" if p < 0.05 else "#BBBBBB" for p in fi2["p_fisher"]]
ax.barh(y, fi2["odds_ratio"], color=col, height=0.62)
for yi, (o, p, a, t) in enumerate(zip(fi2["odds_ratio"], fi2["p_fisher"],
                                      fi2["core_in_top500"], fi2["core_total"])):
    ax.text(o + 0.03, y[yi], f"OR={o:.2f}, p={p:.2g}  ({a}/{t})", va="center", fontsize=5.8)
ax.set_yticks(y)
ax.set_yticklabels([c.replace("_", " ").replace(" vs ", " vs\n") for c in fi2["contrast"]],
                   fontsize=5.6)
ax.axvline(1, color="#888888", lw=0.8, ls="--")
ax.set_xlabel("Fisher OR: core targets among top-500 DE genes")
ax.set_xlim(0, max(fi2["odds_ratio"].max() * 2.2, 2.2))
ttl(ax, "b", "Core-target enrichment (Fisher)")

ax = fig.add_subplot(gs[1, 1])
genes = sorted(cv["gene"].unique())
M = np.full((len(genes), len(contrasts)), np.nan)
S = np.zeros_like(M, dtype=bool)
for j, tag in enumerate(contrasts):
    sub = cv[cv["contrast"] == tag].set_index("gene")
    for i, g in enumerate(genes):
        if g in sub.index:
            M[i, j] = sub.loc[g, "log2FC"]
            S[i, j] = sub.loc[g, "q_bh"] < 0.05
# order rows by mean log2FC
order = np.argsort(np.nanmean(M, axis=1))
M, S, genes = M[order], S[order], [genes[i] for i in order]
vmax = np.nanmax(np.abs(M))
im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax.set_xticks(range(len(contrasts)))
ax.set_xticklabels([c.replace("_", " ").replace(" vs ", " vs ") for c in contrasts],
                   rotation=30, ha="right", fontsize=5.8)
ax.set_yticks(range(len(genes)))
ax.set_yticklabels(genes, fontsize=6.0)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if S[i, j]:
            ax.text(j, i, "*", ha="center", va="center", fontsize=7, color="#111111")
cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.01)
cb.set_label("log2FC (core targets)", fontsize=6.8)
cb.ax.tick_params(labelsize=6.4)
ttl(ax, "c", "Cross-dataset core-target differential expression (* BH q < 0.05)")
save(fig, "Fig5_geo_validation")

# =====================================================================
# Figure S1 - GEO QC / summary
# =====================================================================
fig = plt.figure(figsize=(7.09, 3.4))
gs = fig.add_gridspec(1, 3, wspace=0.6, left=0.08, right=0.98, top=0.86, bottom=0.18)

ax = fig.add_subplot(gs[0, 0])
ds_list = sorted(set(ms["dataset"]))
bottoms = np.zeros(len(ds_list := ds_list if False else ds_list))
bottoms = np.zeros(len(ds_list))
for g in ("normal", "CP", "PDAC", "other"):
    counts = []
    for ds in ds_list:
        sp = pd.read_csv(os.path.join(RES, f"geo_samples_{ds}.csv"))
        counts.append(int((sp["group"] == g).sum()))
    ax.bar(range(len(ds_list)), counts, bottom=bottoms, color=C_GRP[g], label=g, width=0.6)
    for i, (b, c) in enumerate(zip(bottoms, counts)):
        if c:
            ax.text(i, b + c / 2, str(c), ha="center", va="center", fontsize=6, color="white")
    bottoms += np.array(counts)
ax.set_xticks(range(len(ds_list))); ax.set_xticklabels(ds_list, fontsize=6.0, rotation=22, ha="right")
ax.set_ylabel("N samples")
ax.legend(fontsize=5.8)
ttl(ax, "a", "Sample classification per dataset")

ax = fig.add_subplot(gs[0, 1])
per_ds = cv.groupby("dataset")["gene"].nunique()
ax.bar(range(len(per_ds)), per_ds.values, color="#5B7FA6", width=0.6)
for i, v in enumerate(per_ds.values):
    ax.text(i, v + 0.2, str(v), ha="center", fontsize=6.8)
ax.set_xticks(range(len(per_ds))); ax.set_xticklabels(per_ds.index, fontsize=6.0, rotation=22, ha="right")
ax.set_ylim(0, 25)
ax.set_ylabel("N core targets mapped (of 23)")
ttl(ax, "b", "Core-target coverage per platform")

ax = fig.add_subplot(gs[0, 2])
mm2 = ms.drop_duplicates(subset=["gsm", "dataset", "contrast"])
for i, ds in enumerate(ds_list):
    v = mm2.loc[mm2["dataset"] == ds, "module"].values
    bp = ax.boxplot([v], positions=[i], widths=0.5, showfliers=False, patch_artist=True,
                    medianprops={"color": "#222222", "lw": 1.0},
                    whiskerprops={"lw": 0.8}, capprops={"lw": 0.8})
    bp["boxes"][0].set_facecolor("#9DC3E6"); bp["boxes"][0].set_alpha(0.8)
ax.set_xticks(range(len(ds_list))); ax.set_xticklabels(ds_list, fontsize=6.0, rotation=22, ha="right")
ax.axhline(0, color="#888888", lw=0.6, ls=":")
ax.set_ylabel("Core-module score (all samples)")
ttl(ax, "c", "Module-score distribution per dataset")
save(fig, "FigS1_geo_qc")

print("GEO FIGURES DONE")
