"""Step 05: gut-microbiota -> metabolite -> core-target (G-M-T) tripartite network
from gutMGene v2.0 (Table S3 microbe-metabolite links + step-02 metabolite targets).
"""
import os
import pandas as pd
import networkx as nx

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RES = os.path.join(ROOT, "results")

s3 = pd.read_csv(os.path.join(DATA, "Table_S3.csv"))
mt = pd.read_csv(os.path.join(RES, "metabolite_targets.csv"))
core = set(pd.read_csv(os.path.join(RES, "core_targets.csv"))["core_target"])

gm = s3[["Gut Microbiota", "Metabolite", "Rank"]].dropna(subset=["Gut Microbiota", "Metabolite"])
gm["microbe"] = gm["Gut Microbiota"].astype(str).str.strip()
gm["metabolite"] = gm["Metabolite"].astype(str).str.strip()
gm["rank"] = gm["Rank"].astype(str).str.strip().str.lower()
gm = gm[["microbe", "metabolite", "rank"]].drop_duplicates()

mt2 = mt[mt["gene"].isin(core)].copy()
print("M-T edges to core targets:", len(mt2))
gm["met_norm"] = gm["metabolite"].str.lower()
mt2["met_norm"] = mt2["metabolite"].str.lower()
shared = set(gm["met_norm"]) & set(mt2["met_norm"])
print("metabolites bridging G-M and M-T:", len(shared))

gm_e = gm[gm["met_norm"].isin(shared)][["microbe", "metabolite", "met_norm", "rank"]].drop_duplicates()
mt_e = mt2[mt2["met_norm"].isin(shared)][["metabolite", "met_norm", "gene", "source"]].drop_duplicates()
canon = gm[gm["met_norm"].isin(shared)].drop_duplicates("met_norm").set_index("met_norm")["metabolite"].to_dict()
mt_e["metabolite"] = mt_e["met_norm"].map(canon)

G = nx.Graph()
for _, r in gm_e.iterrows():
    G.add_edge(("M", r["microbe"]), ("m", r["metabolite"]))
for _, r in mt_e.iterrows():
    G.add_edge(("m", r["metabolite"]), ("T", r["gene"]))

nodes = pd.DataFrame(
    [{"node": n[1], "layer": {"M": "microbiota", "m": "metabolite", "T": "target"}[n[0]],
      "degree": G.degree(n)} for n in G.nodes()]
).sort_values(["layer", "degree"], ascending=[True, False])
edges = pd.DataFrame(
    [{"layer_pair": f"{a[0]}-{b[0]}", "source": a[1], "target": b[1]}
     for a, b in G.edges()])
nodes.to_csv(os.path.join(RES, "gmt_network_nodes.csv"), index=False)
edges.to_csv(os.path.join(RES, "gmt_network_edges.csv"), index=False)
gm_e.to_csv(os.path.join(RES, "gmt_gm_edges.csv"), index=False)
mt_e.to_csv(os.path.join(RES, "gmt_mt_edges.csv"), index=False)

print("\nnetwork:", G.number_of_nodes(), "nodes,", G.number_of_edges(), "edges")
print("microbiota:", (nodes.layer == "microbiota").sum(),
      "| metabolites:", (nodes.layer == "metabolite").sum(),
      "| targets:", (nodes.layer == "target").sum())

# ---- summaries for the paper ---------------------------------------------------
# focus-class annotation
cls = pd.read_csv(os.path.join(RES, "metabolite_classes.csv")).set_index("metabolite")["class"].to_dict()
met_nodes = nodes[nodes.layer == "metabolite"].copy()
met_nodes["class"] = met_nodes["node"].map(lambda x: cls.get(x, "other"))
print("\nbridging metabolites by class:")
print(met_nodes.groupby("class")["node"].nunique().to_string())

mic_deg = nodes[nodes.layer == "microbiota"].nlargest(15, "degree")
print("\ntop microbiota by degree:")
print(mic_deg[["node", "degree"]].to_string(index=False))

met_deg = nodes[nodes.layer == "metabolite"].nlargest(15, "degree")
print("\ntop metabolites by degree:")
print(met_deg[["node", "degree"]].to_string(index=False))

# microbe -> core target reachability (2-hop via shared metabolite)
rows = []
for mic in gm_e["microbe"].unique():
    mets = gm_e.loc[gm_e["microbe"] == mic, "met_norm"]
    tg = mt_e[mt_e["met_norm"].isin(mets)]["gene"].unique()
    rows.append({"microbe": mic, "n_metabolites": len(mets), "n_core_targets": len(tg),
                 "core_targets": ";".join(sorted(tg))})
reach = pd.DataFrame(rows).sort_values("n_core_targets", ascending=False)
reach.to_csv(os.path.join(RES, "gmt_microbe_reach.csv"), index=False)
print("\ntop microbes reaching core targets:")
print(reach.head(12).to_string(index=False))
print("DONE")
