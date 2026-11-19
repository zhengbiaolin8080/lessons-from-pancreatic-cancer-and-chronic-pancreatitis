"""Step 07c: GEO validation for PDAC and CP.

Parses series-matrix files, maps probes to gene symbols via GPL annotation,
classifies samples into groups by keywords, and computes:
  - per-gene DE for the 23 core targets (log2FC, MWU p, BH q) per contrast
  - core-module score per sample (mean z of core genes), group comparison
  - Fisher enrichment of core targets among significant DE genes
  - cross-dataset direction consistency vs TCGA tumour-vs-normal
"""
import os
import gzip
import re
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "GEO")
RES = os.path.join(ROOT, "results")

core = pd.read_csv(os.path.join(RES, "core_targets.csv"))["core_target"].tolist()
tcga_dir = pd.read_csv(os.path.join(RES, "tcga_core_gene_expression.csv"))
tcga_l2fc = dict(zip(tcga_dir["gene"], tcga_dir["log2FC_tumor_vs_normal"]))

DATASETS = {
    "GSE28735": {"gpl": "GPL6244",
                 "groups": {"PDAC": ["tumor", "tumour", "cancer", "carcinoma"],
                            "normal": ["non-tumor", "nontumor", "adjacent", "normal"]},
                 "contrasts": [("PDAC", "normal")]},
    "GSE62452": {"gpl": "GPL6244",
                 "groups": {"PDAC": ["tumor", "tumour", "cancer", "carcinoma"],
                            "normal": ["non-tumor", "nontumor", "adjacent", "normal"]},
                 "contrasts": [("PDAC", "normal")]},
    "GSE77858": {"gpl": "GPL7264",
                 "groups": {"PDAC": ["pdac: y"],
                            "CP": ["panreatitis", "pancreatitis"],
                            "normal": ["morphology: normal"]},
                 "contrasts": [("PDAC", "normal"), ("CP", "normal"), ("PDAC", "CP")]},
    "GSE101462": {"gpl": "GPL10558",
                  "groups": {"PDAC": ["tumor", "tumour", "cancer", "carcinoma", "pdac"],
                             "CP": ["pancreatitis"],
                             "normal": ["normal", "control", "healthy"]},
                  "contrasts": [("PDAC", "CP"), ("PDAC", "normal"), ("CP", "normal")]},
}


def read_series_matrix(path):
    meta = {}
    with gzip.open(path, "rt", errors="replace") as f:
        in_table = False
        header = None
        rows = []
        for line in f:
            line = line.rstrip("\n")
            if not in_table:
                if line.startswith("!series_matrix_table_begin"):
                    in_table = True
                elif line.startswith("!Sample_"):
                    k, _, v = line.partition("\t")
                    vals = [x.strip('"') for x in v.split("\t")]
                    if k in meta:  # repeated key (e.g. two characteristics rows)
                        meta[k] = [a + " | " + b for a, b in zip(meta[k], vals)]
                    else:
                        meta[k] = vals
                continue
            if header is None:
                header = [h.strip('"') for h in line.split("\t")[1:]]
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            parts = line.split("\t")
            rows.append(parts)
    ncol = len(header) + 1
    rows = [r for r in rows if len(r) == ncol]  # drop malformed lines
    df = pd.DataFrame([r[1:] for r in rows],
                      index=[r[0].strip('"') for r in rows],
                      columns=header)
    df = df.apply(pd.to_numeric, errors="coerce")  # '' -> NaN
    return df, meta


def _parse_platform_table(path):
    """Parse the !platform_table_begin..end block of a SOFT/annot file."""
    rows = []
    header = None
    with gzip.open(path, "rt", errors="replace") as f:
        in_pt = False
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_pt = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_pt:
                if header is None:
                    header = line.split("\t")
                    continue
                parts = line.split("\t")
                if len(parts) == len(header):
                    rows.append(parts)
    if header is None:
        raise ValueError(f"no platform table in {path}")
    return pd.DataFrame(rows, columns=header)


def gpl_map(gpl):
    p = os.path.join(DATA, gpl + ".annot.gz")
    if os.path.exists(p):
        annot = _parse_platform_table(p)
    else:
        # fallback: parse platform table from GSE family SOFT (GPL7264 case)
        soft = os.path.join(DATA, "GSE77858_family.soft.gz")
        annot = _parse_platform_table(soft)
    idcol = annot.columns[0]
    symcol = None
    for c in annot.columns:
        cl = c.strip().lower().replace("_", " ")
        if cl in ("gene symbol", "symbol", "ilmn gene", "genesymbol", "gene assignment"):
            symcol = c
            break
    if symcol is None:
        for c in annot.columns:
            if "symbol" in c.lower():
                symcol = c
                break
    print(f"  {gpl}: id={idcol!r} symbol={symcol!r}")
    m = {}
    for pid, sym in zip(annot[idcol], annot[symcol]):
        if pd.isna(sym):
            continue
        s = re.split(r"///|;", str(sym))[0].strip()
        if s and s != "---":
            m[str(pid)] = s
    return m


def classify(meta, rules):
    n = len(meta["!Sample_geo_accession"])
    texts = []
    # concatenate classification-relevant fields only (title / source / characteristics,
    # both channels for two-colour arrays); dataset-level description/protocol fields
    # would leak keywords (e.g. "non-tumor", "normalized") into every sample
    keys = [k for k in meta if k.startswith(("!Sample_title", "!Sample_source_name",
                                             "!Sample_characteristics"))]
    for i in range(n):
        t = " | ".join(meta.get(k, [""] * n)[i] for k in keys).lower()
        texts.append(t)
    grp = []
    for t in texts:
        g = "other"
        # CP first (a CP sample may also mention tumour absence), then PDAC, then normal
        for order_g in (["CP", "PDAC", "normal"] if "CP" in rules else list(rules)):
            kws = rules[order_g]
            if any(k in t for k in kws):
                if order_g == "PDAC" and ("non-tumor" in t or "non-tumour" in t
                                          or "adjacent" in t):
                    continue
                if order_g == "normal" and ("pancreatitis" in t or "panreatitis" in t):
                    continue
                g = order_g
                break
        grp.append(g)
    return pd.DataFrame({"gsm": meta["!Sample_geo_accession"], "group": grp,
                         "text": texts})


all_module = []
core_val_rows = []
fisher_rows = []

for gse, cfg in DATASETS.items():
    f = os.path.join(DATA, f"{gse}_series_matrix.txt.gz")
    if not os.path.exists(f):
        print("SKIP missing", gse)
        continue
    print("===", gse, "===")
    try:
        pmap = gpl_map(cfg["gpl"])
    except Exception as e:
        print(f"SKIP {gse}: platform map unavailable ({e})")
        continue
    expr, meta = read_series_matrix(f)
    expr.index = expr.index.map(lambda x: pmap.get(str(x), ""))
    expr = expr[expr.index != ""]
    expr = expr.groupby(expr.index).mean()  # collapse probes by mean
    if np.nanmax(expr.values) > 50:  # not log-transformed
        expr = np.log2(expr + 1.0)
    print("  genes:", expr.shape[0], "samples:", expr.shape[1])

    samples = classify(meta, cfg["groups"])
    print(samples.groupby("group").size().to_dict())
    samples.to_csv(os.path.join(RES, f"geo_samples_{gse}.csv"), index=False)

    z = expr.apply(lambda r: (r - r.mean()) / (r.std(ddof=0) or np.nan), axis=1)
    core_in = [g for g in core if g in z.index]
    print("  core genes mapped:", len(core_in))
    mod = z.loc[core_in].mean()

    for a, b in cfg["contrasts"]:
        sa = samples.loc[samples["group"] == a, "gsm"].tolist()
        sb = samples.loc[samples["group"] == b, "gsm"].tolist()
        if len(sa) < 3 or len(sb) < 3:
            print(f"  contrast {a}({len(sa)}) vs {b}({len(sb)}): too few, skip")
            continue
        tag = f"{gse}_{a}_vs_{b}"
        # module score
        ms = pd.DataFrame({"gsm": expr.columns,
                           "module": mod.reindex(expr.columns).values})
        ms["group"] = ms["gsm"].map(dict(zip(samples["gsm"], samples["group"])))
        ms["dataset"] = gse
        ms["contrast"] = tag
        all_module.append(ms[ms["group"].isin([a, b])])
        p_mod = stats.mannwhitneyu(ms.loc[ms["group"] == a, "module"],
                                   ms.loc[ms["group"] == b, "module"]).pvalue
        print(f"  {tag}: n={len(sa)}v{len(sb)} module MWU p={p_mod:.3g}")
        # genome-wide DE (core focus + fisher universe)
        rows = []
        for g in core_in:
            x = expr.loc[g, sa].dropna().astype(float)
            y = expr.loc[g, sb].dropna().astype(float)
            l2 = x.mean() - y.mean()
            p = stats.mannwhitneyu(x, y).pvalue
            rows.append({"gene": g, "log2FC": l2, "p": p})
        de = pd.DataFrame(rows)
        de["q_bh"] = stats.false_discovery_control(de["p"])
        de["dataset"] = gse
        de["contrast"] = tag
        de["tcga_log2FC"] = de["gene"].map(tcga_l2fc)
        de["same_direction_as_tcga"] = np.sign(de["log2FC"]) == np.sign(de["tcga_log2FC"])
        de.to_csv(os.path.join(RES, f"geo_de_{tag}.csv"), index=False)
        n_sig = ((de["q_bh"] < 0.05)).sum()
        n_dir = ((de["q_bh"] < 0.1) & (de["same_direction_as_tcga"])).sum()
        print(f"    core DE: {n_sig}/{len(de)} q<0.05 | {n_dir} q<0.1 & TCGA-consistent")
        for _, r in de.iterrows():
            core_val_rows.append({**r.to_dict(),
                                  "module_p": p_mod, "n_a": len(sa), "n_b": len(sb)})
        # Fisher: core among genome-wide significant
        lp = []
        for g in expr.index:
            x = expr.loc[g, sa].dropna().astype(float)
            y = expr.loc[g, sb].dropna().astype(float)
            if len(x) < 3 or len(y) < 3 or x.std() + y.std() == 0:
                continue
            tstat, pv = stats.ttest_ind(x, y, equal_var=False)
            lp.append((g, abs(tstat), pv))
        gw = pd.DataFrame(lp, columns=["gene", "abs_t", "p"]).dropna()
        top = set(gw.nsmallest(500, "p")["gene"])
        cs = set(core_in)
        a_ = len(top & cs); b_ = len(cs - top & cs)
        c_ = len(top - cs); d_ = len(gw) - a_ - b_ - c_
        or_, pf = stats.fisher_exact([[a_, b_], [c_, d_]])
        fisher_rows.append({"dataset": gse, "contrast": tag, "core_in_top500": a_,
                            "core_total": len(cs), "odds_ratio": round(a_ * d_ / max(b_ * c_, 1), 2),
                            "p_fisher": pf})

if all_module:
    pd.concat(all_module).to_csv(os.path.join(RES, "geo_module_scores.csv"), index=False)
pd.DataFrame(core_val_rows).to_csv(os.path.join(RES, "geo_core_validation.csv"), index=False)
pd.DataFrame(fisher_rows).to_csv(os.path.join(RES, "geo_fisher.csv"), index=False)
print("DONE")
