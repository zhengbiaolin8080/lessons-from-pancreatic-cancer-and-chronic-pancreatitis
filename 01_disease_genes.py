"""Step 01: disease-associated genes for PDAC and chronic pancreatitis (Open Targets).

Design: chronic pancreatitis (CP) is the strongest inflammatory risk factor for
pancreatic ductal adenocarcinoma (PDAC). Genes associated with BOTH diseases
constitute the fibroinflammatory-to-neoplastic transition substrate; genes
associated with either disease form the background disease-space.
Disease IDs are resolved by name search against the Open Targets GraphQL API.
"""
import os
import requests
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RES = os.path.join(ROOT, "results")

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"


def ot_query(query, variables=None):
    r = requests.post(OT_URL, json={"query": query, "variables": variables or {}}, timeout=120)
    r.raise_for_status()
    return r.json()


Q_SEARCH = """
query($q:String!){ search(queryString:$q){ hits{ id name entity description } } }
"""

Q_ASSOC = """
query($id:String!,$size:Int!){
  disease(efoId:$id){
    id name
    associatedTargets(page:{index:0,size:$size}){
      count
      rows{ target{ approvedSymbol approvedName } score }
    }
  }
}
"""

# ---- resolve disease IDs ---------------------------------------------------
for q in ["pancreatic ductal adenocarcinoma", "chronic pancreatitis"]:
    hits = ot_query(Q_SEARCH, {"q": q})["data"]["search"]["hits"]
    print("##", q)
    for h in hits[:8]:
        if h["entity"] == "disease":
            print("   ", h["id"], "|", h["name"], "|", (h.get("description") or "")[:80])

# ---- associated targets ------------------------------------------------------
DISEASES = {"PDAC": "MONDO_0005184", "CP": "MONDO_0005003"}
SCORE = 0.2

frames = []
for tag, did in DISEASES.items():
    snap = os.path.join(DATA, f"disease_genes_{tag}.csv")
    if os.path.exists(snap):
        d = pd.read_csv(snap)
    else:
        out = ot_query(Q_ASSOC, {"id": did, "size": 3000})["data"]["disease"]
        rows = [{"symbol": r["target"]["approvedSymbol"],
                 "name": r["target"]["approvedName"],
                 "score": r["score"], "disease": tag, "disease_id": out["id"],
                 "disease_name": out["name"]}
                for r in out["associatedTargets"]["rows"] if r["target"].get("approvedSymbol")]
        d = pd.DataFrame(rows)
        d.to_csv(snap, index=False)
    print(f"{tag} ({did}, {d['disease_name'].iloc[0]}): {len(d)} associated targets")
    frames.append(d)

both = pd.concat(frames, ignore_index=True)
sel = both[both["score"] >= SCORE].copy()
pdac = set(sel.loc[sel["disease"] == "PDAC", "symbol"])
cp = set(sel.loc[sel["disease"] == "CP", "symbol"])
shared = sorted(pdac & cp)
union = sorted(pdac | cp)
print(f"PDAC genes (score>={SCORE}): {len(pdac)}")
print(f"CP genes (score>={SCORE}): {len(cp)}")
print(f"shared PDAC∩CP: {len(shared)} | union: {len(union)}")

# overlap table with per-disease scores
p = sel[sel["disease"] == "PDAC"].set_index("symbol")
c = sel[sel["disease"] == "CP"].set_index("symbol")
rows = []
for g in shared:
    rows.append({"symbol": g,
                 "name": p.loc[g, "name"] if g in p.index else c.loc[g, "name"],
                 "pdac_score": round(float(p.loc[g, "score"]), 4) if g in p.index else None,
                 "cp_score": round(float(c.loc[g, "score"]), 4) if g in c.index else None})
ov = pd.DataFrame(rows)
ov["min_score"] = ov[["pdac_score", "cp_score"]].min(axis=1)
ov = ov.sort_values("min_score", ascending=False)
ov.to_csv(os.path.join(RES, "pdac_cp_shared_genes.csv"), index=False)

# full union with disease membership flags
urows = []
for g in union:
    urows.append({"symbol": g,
                  "name": p.loc[g, "name"] if g in p.index else c.loc[g, "name"],
                  "in_PDAC": g in pdac, "in_CP": g in cp,
                  "pdac_score": round(float(p.loc[g, "score"]), 4) if g in p.index else None,
                  "cp_score": round(float(c.loc[g, "score"]), 4) if g in c.index else None})
udf = pd.DataFrame(urows).sort_values(["in_PDAC", "in_CP", "symbol"], ascending=[False, False, True])
udf.to_csv(os.path.join(RES, "disease_genes_union.csv"), index=False)
sel.drop_duplicates(["symbol", "disease"]).to_csv(os.path.join(RES, "disease_genes_selected.csv"), index=False)

print("\nTop shared genes by min(disease scores):")
print(ov.head(30).to_string(index=False))
