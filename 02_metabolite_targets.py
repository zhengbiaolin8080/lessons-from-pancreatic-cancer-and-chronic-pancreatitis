"""Step 02: metabolite panel (4 focus classes), gutMGene curated targets,
ChEMBL similarity-based target prediction, and core-target calling
(metabolite targets n PDAC+CP disease gene space).

Focus classes (per study design): indole/tryptophan derivatives, short-chain
fatty acids (SCFAs), bile acids, polyamines. Curated evidence: gutMGene v2.0
Table S2 (metabolite->host gene). Predicted evidence: ChEMBL structural
analogues (Tanimoto >= 70) with human binding activities (pChEMBL >= 6).
"""
import os
import re
import json
import time
import requests
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RES = os.path.join(ROOT, "results")

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# ---- metabolite universe & classes ------------------------------------------
cat = pd.read_csv(os.path.join(DATA, "metabolites.csv"))
cat = cat[cat["smiles"].notna()].copy()
s2 = pd.read_csv(os.path.join(DATA, "Table_S2.csv"))

CLASS_PATTERNS = [
    ("bile_acid", re.compile(r"cholic|cholanoic|deoxychol|lithochol|chenodeoxychol|ursodeoxychol|"
                             r"ursodiol|chenodiol|murichol|hyochol|hyodeoxychol|tauro.*chol|"
                             r"glyco.*chol|isoallolithochol|isodeoxychol|isolithochol|"
                             r"murideoxychol|ketolithochol|sulfochol", re.I)),
    ("polyamine", re.compile(r"putrescine|spermidine|spermine|cadaverine|agmatine|"
                             r"^L-ornithine$|^ornithine$", re.I)),
    ("indole_trp", re.compile(r"indol|tryptophan|kynuren|skatole|tryptamine|serotonin|"
                              r"melatonin|oxindole", re.I)),
    ("scfa", re.compile(r"^acetate$|^propionate$|^butyrate$|^valerate$|^valeric acid$|"
                        r"^isobutyric acid$|^isovaleric acid$|^isovalerate$|^caproic acid$|"
                        r"^2-hydroxybutyric acid$|^3-hydroxybutyric acid$|"
                        r"^l-2-aminobutyric acid$", re.I)),
]


def classify(name):
    for cls, rx in CLASS_PATTERNS:
        if rx.search(str(name)):
            return cls
    return "other"


cat["class"] = cat["name"].map(classify)
s2["class"] = s2["Metabolite"].astype(str).str.strip().map(classify)
cls_tab = pd.concat([
    cat[["name", "class"]].rename(columns={"name": "metabolite"}),
    s2[["Metabolite", "class"]].rename(columns={"Metabolite": "metabolite"}),
]).drop_duplicates("metabolite")
# focus-class membership takes priority if any source flags it
focus = cls_tab[cls_tab["class"] != "other"].groupby("metabolite")["class"].agg(
    lambda s: sorted(set(s))[0]).reset_index()
cls_tab = cls_tab.set_index("metabolite")
for _, r0 in focus.iterrows():
    cls_tab.loc[r0["metabolite"], "class"] = r0["class"]
cls_tab = cls_tab.reset_index()
cls_tab.to_csv(os.path.join(RES, "metabolite_classes.csv"), index=False)
print("class counts (gutMGene S2):", s2.groupby("class")["Metabolite"].nunique().to_dict())
print("class counts (catalog):", cat.groupby("class")["name"].nunique().to_dict())

FOCUS = {"indole_trp", "scfa", "bile_acid", "polyamine"}
panel = cat[cat["class"].isin(FOCUS)].copy().sort_values(["class", "name"])
panel.to_csv(os.path.join(RES, "focus_metabolite_panel.csv"), index=False)
print("focus panel (SMILES available):", panel.groupby("class")["name"].nunique().to_dict(),
      "| total:", len(panel))

# ---- ChEMBL similarity prediction (cached snapshots) -------------------------
SIM_F = os.path.join(DATA, "chembl_similar.json")
ACT_F = os.path.join(DATA, "chembl_activities.json")
TGT_F = os.path.join(DATA, "chembl_targets.json")


def chembl_get(path, params=None, tries=3):
    for i in range(tries):
        try:
            r = requests.get(CHEMBL + path, params=params, timeout=90)
            if r.status_code == 200:
                return r.json()
        except Exception:
            pass
        time.sleep(1.5 * (i + 1))
    return {}


if os.path.exists(SIM_F) and os.path.exists(ACT_F) and os.path.exists(TGT_F):
    sim_map = json.load(open(SIM_F))
    acts = json.load(open(ACT_F))
    tmap = json.load(open(TGT_F))
    print("loaded cached ChEMBL snapshots")
else:
    # 1) structural analogues (Tanimoto >= 70)
    sim_map = {}
    for _, r0 in panel.iterrows():
        name, smi = r0["name"], r0["smiles"]
        j = chembl_get(f"/similarity/{requests.utils.quote(smi, safe='')}/70.json",
                       params={"limit": 80})
        mols = [{"chembl_id": m["molecule_chembl_id"], "similarity": m["similarity"]}
                for m in j.get("molecules", [])]
        sim_map[name] = mols
        print(f"similarity {name}: {len(mols)}", flush=True)
        time.sleep(0.25)
    json.dump(sim_map, open(SIM_F, "w"))

    # 2) activities of all analogue molecules (batch __in)
    all_mols = sorted({m["chembl_id"] for v in sim_map.values() for m in v})
    print("unique analogue molecules:", len(all_mols))
    acts = []
    B = 40
    for i in range(0, len(all_mols), B):
        chunk = all_mols[i:i + B]
        params = {"molecule_chembl_id__in": ",".join(chunk), "limit": 1000,
                  "target_organism": "Homo sapiens"}
        j = chembl_get("/activity.json", params=params)
        for a in j.get("activities", []):
            try:
                pchembl = float(a.get("pchembl_value") or 0)
            except ValueError:
                pchembl = 0
            if pchembl < 6:
                continue
            acts.append({"molecule_chembl_id": a["molecule_chembl_id"],
                         "target_chembl_id": a["target_chembl_id"],
                         "target_pref_name": a.get("target_pref_name"),
                         "pchembl_value": a.get("pchembl_value"),
                         "standard_type": a.get("standard_type")})
        print(f"activities {i + len(chunk)}/{len(all_mols)} -> {len(acts)}", flush=True)
        time.sleep(0.25)
    json.dump(acts, open(ACT_F, "w"))

    # 3) target -> human gene symbol (single proteins)
    tids = sorted({a["target_chembl_id"] for a in acts})
    tmap = {}
    B = 40
    for i in range(0, len(tids), B):
        chunk = tids[i:i + B]
        j = chembl_get("/target.json", params={"target_chembl_id__in": ",".join(chunk),
                                               "limit": 200})
        for t in j.get("targets", []):
            gene = None
            if t.get("target_type") == "SINGLE PROTEIN":
                for c in t.get("target_components", []):
                    sy = (c.get("target_component_synonyms") or [])
                    gs = [x["component_synonym"] for x in sy
                          if x.get("syn_type") == "GENE_SYMBOL"]
                    if gs:
                        gene = gs[0]
                        break
            tmap[t["target_chembl_id"]] = {"pref_name": t.get("pref_name"),
                                           "type": t.get("target_type"),
                                           "gene": gene}
        print(f"targets {i + len(chunk)}/{len(tids)}", flush=True)
        time.sleep(0.25)
    json.dump(tmap, open(TGT_F, "w"))

# ---- merge predicted + curated ----------------------------------------------
mol2gene = {}
for a in acts:
    g = (tmap.get(a["target_chembl_id"]) or {}).get("gene")
    if g:
        mol2gene.setdefault(a["molecule_chembl_id"], set()).add(g)

rows = []
for name, hits in sim_map.items():
    genes = set()
    for h in hits:
        genes |= mol2gene.get(h["chembl_id"], set())
    for g in genes:
        rows.append({"metabolite": name, "gene": g, "source": "ChEMBL_similarity"})
pred = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["metabolite", "gene", "source"])

cur = s2[["Metabolite", "Gene"]].dropna()
cur.columns = ["metabolite", "gene"]
cur["metabolite_norm"] = cur["metabolite"].astype(str).str.strip().str.lower()
cat_names = {n.lower(): n for n in cat["name"]}
cur = cur[cur["metabolite_norm"].isin(cat_names.keys())]
cur["metabolite"] = cur["metabolite_norm"].map(cat_names)
cur = cur[["metabolite", "gene"]].drop_duplicates()
cur["source"] = "gutMGene_curated"

allmt = pd.concat([pred, cur], ignore_index=True).drop_duplicates(["metabolite", "gene", "source"])
allmt.to_csv(os.path.join(RES, "metabolite_targets.csv"), index=False)

uni_pred, uni_cur = set(pred["gene"]), set(cur["gene"])
print("ChEMBL predicted genes:", len(uni_pred), "| curated:", len(uni_cur),
      "| union:", len(uni_pred | uni_cur), "| dual:", len(uni_pred & uni_cur))

cnt = allmt.groupby("metabolite")["gene"].nunique().sort_values(ascending=False)
cnt.rename("n_targets").to_csv(os.path.join(RES, "targets_per_metabolite.csv"))

# ---- core targets: metabolite targets n disease gene space -------------------
dz = pd.read_csv(os.path.join(RES, "disease_genes_union.csv"))
dspace = set(dz["symbol"])
core = sorted((uni_pred | uni_cur) & dspace)
dual = sorted((uni_pred & uni_cur) & dspace)
pd.DataFrame({"core_target": core}).to_csv(os.path.join(RES, "core_targets.csv"), index=False)
pd.DataFrame({"dual_supported_core": dual}).to_csv(os.path.join(RES, "core_targets_dual.csv"),
                                                   index=False)
det = allmt[allmt["gene"].isin(core)].copy()
det.to_csv(os.path.join(RES, "core_target_evidence.csv"), index=False)
print("disease gene space:", len(dspace))
print("CORE TARGETS:", len(core), "| dual-supported:", len(dual))
print(core)

# ---- AHR / STAT3 / NF-kB axis annotation -------------------------------------
AXES = {
    "AHR": ["AHR", "ARNT", "AHRR", "AIP", "CYP1A1", "CYP1A2", "CYP1B1", "CYP2S1",
            "TIPARP", "NQO1", "ALDH3A1", "UGT1A1", "UGT1A6", "IL22"],
    "STAT3": ["STAT3", "JAK1", "JAK2", "JAK3", "TYK2", "IL6", "IL6R", "IL6ST",
              "IL11", "LIF", "OSM", "SOCS3", "PIM1", "MCL1", "BCL2L1", "CCND1",
              "MYC", "VEGFA", "HIF1A"],
    "NFkB": ["NFKB1", "NFKB2", "RELA", "RELB", "REL", "IKBKB", "IKBKG", "CHUK",
             "NFKBIA", "TNF", "IL1B", "TLR4", "MYD88"],
}
ax_rows = []
for ax, genes in AXES.items():
    for g in genes:
        ax_rows.append({"axis": ax, "gene": g,
                        "in_core": g in core,
                        "in_disease_space": g in dspace,
                        "is_metabolite_target": g in (uni_pred | uni_cur)})
axdf = pd.DataFrame(ax_rows)
axdf.to_csv(os.path.join(RES, "axis_membership.csv"), index=False)
print("\nAxis membership in core targets:")
print(axdf[axdf["in_core"]].to_string(index=False))