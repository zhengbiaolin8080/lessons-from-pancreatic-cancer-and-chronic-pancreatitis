"""Step 06: TCGA-PAAD validation.

- module score = mean z-scored expression of the 23 core targets (and AHR /
  STAT3 / NFkB axis modules) across primary tumours
- Cox PH (OS, PFI; PanCanAtlas Liu 2018 outcomes) continuous per-SD and
  median-dichotomised; multivariable adjustment for age + stage group
- per-gene univariate Cox OS
- Spearman correlation of core module with immune/stromal marker signatures
"""
import os
import gzip
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "TCGA")
RES = os.path.join(ROOT, "results")

core = pd.read_csv(os.path.join(RES, "core_targets.csv"))["core_target"].tolist()
axis_df = pd.read_csv(os.path.join(RES, "axis_membership.csv"))
axis_sets = {a: sorted(g["gene"].tolist()) for a, g in axis_df.groupby("axis")}

# ---- expression ----------------------------------------------------------------
print("loading expression ...")
expr = pd.read_csv(os.path.join(DATA, "PAAD_HiSeqV2.gz"), sep="\t", index_col=0)
expr = expr[~expr.index.duplicated(keep="first")]
tum = [c for c in expr.columns if c.split("-")[-1] == "01"]
norm = [c for c in expr.columns if c.split("-")[-1] == "11"]
print("samples: tumor", len(tum), "| normal", len(norm))

z = expr[tum].apply(lambda r: (r - r.mean()) / (r.std(ddof=0) or np.nan), axis=1)

core_in = [g for g in core if g in z.index]
print("core targets in matrix:", len(core_in), "| missing:", sorted(set(core) - set(core_in)))
mod = z.loc[core_in].mean()
scores = pd.DataFrame({"sample": mod.index, "core_module": mod.values})

for ax, genes in axis_sets.items():
    gi = [g for g in genes if g in z.index]
    scores[ax + "_module"] = z.loc[gi].mean().values
    print(ax, "module genes:", len(gi))

# per-sample tumour-vs-normal direction for core genes
dirn = pd.DataFrame({
    "gene": core_in,
    "mean_tumor": expr.loc[core_in, tum].mean(axis=1).values,
    "mean_normal": expr.loc[core_in, norm].mean(axis=1).values if norm else np.nan,
})
dirn["log2FC_tumor_vs_normal"] = dirn["mean_tumor"] - dirn["mean_normal"]
wil_p = []
for g in core_in:
    a = expr.loc[g, tum].astype(float)
    b = expr.loc[g, norm].astype(float)
    wil_p.append(stats.mannwhitneyu(a, b).pvalue if len(b) > 2 else np.nan)
dirn["p_mwu"] = wil_p
dirn["q_bh"] = stats.false_discovery_control(dirn["p_mwu"]) if dirn["p_mwu"].notna().all() else np.nan
dirn.to_csv(os.path.join(RES, "tcga_core_gene_expression.csv"), index=False)

# ---- survival -------------------------------------------------------------------
surv = pd.read_csv(os.path.join(DATA, "TCGA_survival_pancanatlas.tsv"), sep="\t")
surv = surv[surv["cancer type abbreviation"] == "PAAD"]
surv["sample01"] = surv["sample"]
df = scores.merge(surv[["sample01", "OS", "OS.time", "PFI", "PFI.time"]],
                  left_on="sample", right_on="sample01", how="left")
print("module scores:", len(df), "| with OS:", df["OS"].notna().sum(),
      "| with PFI:", df["PFI"].notna().sum())

clin = pd.read_csv(os.path.join(DATA, "PAAD_clinicalMatrix"), sep="\t", low_memory=False)
clin = clin.set_index("sampleID")
stage_map = {}
for s, v in clin["pathologic_stage"].items():
    v = str(v).upper()
    g = np.nan
    if v.startswith("STAGE"):
        v2 = v.replace("STAGE", "").strip()
        if v2.startswith("I") and not v2.startswith("IV"):
            g = 1 if v2 in ("I", "IA", "IB") else 2
        elif v2.startswith("IV"):
            g = 4
        elif v2.startswith("III"):
            g = 3
    stage_map[s] = g
df["stage_group"] = df["sample"].map(stage_map)
df["age"] = pd.to_numeric(
    clin["age_at_initial_pathologic_diagnosis"].reindex(df["sample"]).values,
    errors="coerce")


def cox_block(data, time, event, label):
    d = data[["sample", time, event, "core_module", "age", "stage_group"]].dropna(
        subset=[time, event, "core_module"])
    d = d[d[time] > 0]
    out = {}
    zd = d.copy()
    zd["core_z"] = (zd["core_module"] - zd["core_module"].mean()) / zd["core_module"].std()
    cph = CoxPHFitter()
    cph.fit(zd[[time, event, "core_z"]], time, event)
    s = cph.summary.loc["core_z"]
    out[f"{label}_cont_HR"] = round(np.exp(s["coef"]), 3)
    out[f"{label}_cont_CI"] = f"{np.exp(s['coef lower 95%']):.3f}-{np.exp(s['coef upper 95%']):.3f}"
    out[f"{label}_cont_p"] = s["p"]
    mv = zd.dropna(subset=["age", "stage_group"])
    if mv["stage_group"].nunique() > 1 and len(mv) > 30:
        cph2 = CoxPHFitter()
        cph2.fit(mv[[time, event, "core_z", "age", "stage_group"]], time, event)
        s2 = cph2.summary.loc["core_z"]
        out[f"{label}_adj_HR"] = round(np.exp(s2["coef"]), 3)
        out[f"{label}_adj_CI"] = f"{np.exp(s2['coef lower 95%']):.3f}-{np.exp(s2['coef upper 95%']):.3f}"
        out[f"{label}_adj_p"] = s2["p"]
        out[f"{label}_adj_n"] = len(mv)
    hi = d["core_module"] >= d["core_module"].median()
    lr = logrank_test(d.loc[hi, time], d.loc[~hi, time], d.loc[hi, event], d.loc[~hi, event])
    out[f"{label}_median_logrank_p"] = lr.p_value
    out[f"{label}_n"] = len(d)
    out[f"{label}_events"] = int(d[event].sum())
    return out, d


res_os, dos = cox_block(df, "OS.time", "OS", "OS")
res_pfi, dpfi = cox_block(df, "PFI.time", "PFI", "PFI")
pd.DataFrame([{**{"cohort": "TCGA-PAAD"}, **res_os}]).to_csv(
    os.path.join(RES, "tcga_cox_os.csv"), index=False)
pd.DataFrame([{**{"cohort": "TCGA-PAAD"}, **res_pfi}]).to_csv(
    os.path.join(RES, "tcga_cox_pfi.csv"), index=False)
print("\nOS:", res_os)
print("\nPFI:", res_pfi)
scores.to_csv(os.path.join(RES, "tcga_module_scores.csv"), index=False)

# ---- per-gene univariate Cox OS ---------------------------------------------------
rows = []
d = dos
for g in core_in:
    gx = z.loc[g, d["sample"]].values
    dd = pd.DataFrame({"t": d["OS.time"].values, "e": d["OS"].values, "x": gx})
    try:
        cph = CoxPHFitter()
        cph.fit(dd, "t", "e")
        s = cph.summary.loc["x"]
        rows.append({"gene": g, "HR": round(np.exp(s["coef"]), 3),
                     "CI_low": round(np.exp(s["coef lower 95%"]), 3),
                     "CI_high": round(np.exp(s["coef upper 95%"]), 3), "p": s["p"]})
    except Exception:
        pass
gcox = pd.DataFrame(rows).sort_values("p")
gcox["q_bh"] = stats.false_discovery_control(gcox["p"])
gcox.to_csv(os.path.join(RES, "tcga_gene_cox_os.csv"), index=False)
print("\nper-gene Cox OS (q<0.1):")
print(gcox[gcox["q_bh"] < 0.1].to_string(index=False))

# ---- immune/stromal signatures ----------------------------------------------------
SIG = {
    "CD8_T": ["CD8A", "CD8B"],
    "T_cell": ["CD3D", "CD3E", "CD3G"],
    "Treg": ["FOXP3", "IL2RA"],
    "B_cell": ["CD19", "MS4A1", "CD79A"],
    "NK": ["NCAM1", "NKG7", "KLRD1"],
    "Macrophage": ["CD68", "CD163", "CD14"],
    "M2_macrophage": ["CD163", "MRC1", "MSR1"],
    "Neutrophil": ["FCGR3B", "S100A8", "S100A9", "CXCR2"],
    "DC": ["ITGAX", "CD1C", "CLEC9A"],
    "CAF": ["ACTA2", "FAP", "COL1A1", "PDGFRB"],
    "Endothelial": ["PECAM1", "VWF", "KDR"],
    "Checkpoint": ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2"],
}
sig_scores = {}
for k, genes in SIG.items():
    gi = [g for g in genes if g in z.index]
    sig_scores[k] = z.loc[gi, tum].mean()
sig_df = pd.DataFrame(sig_scores)
sig_df.index.name = "sample"
sig_df.to_csv(os.path.join(RES, "tcga_immune_scores.csv"))

rows = []
for k in SIG:
    x = mod.reindex(sig_df.index)
    y = sig_df[k]
    ok = x.notna() & y.notna()
    r, p = stats.spearmanr(x[ok], y[ok])
    rows.append({"signature": k, "spearman_r": round(r, 3), "p": p, "n": int(ok.sum())})
ic = pd.DataFrame(rows).sort_values("p")
ic["q_bh"] = stats.false_discovery_control(ic["p"])
ic.to_csv(os.path.join(RES, "tcga_immune_corr.csv"), index=False)
print("\nimmune/stromal correlations with core module:")
print(ic.to_string(index=False))

# axis modules vs core module
rows = []
for ax in axis_sets:
    r, p = stats.spearmanr(scores["core_module"], scores[ax + "_module"])
    rows.append({"axis": ax, "spearman_r": round(r, 3), "p": p})
ac = pd.DataFrame(rows)
ac.to_csv(os.path.join(RES, "tcga_axis_corr.csv"), index=False)
print("\naxis module correlations with core module:")
print(ac.to_string(index=False))
print("DONE")
