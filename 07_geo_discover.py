"""Step 07a: discover GEO expression datasets for chronic pancreatitis (CP)
and verify known PDAC-vs-normal datasets, via NCBI eutils (GDS database).

Prints a ranked candidate table and saves results/geo_candidates.csv for
manual selection before download/analysis (07_geo_analysis.py).
"""
import os
import time
import requests
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def get(path, params, tries=3):
    for i in range(tries):
        try:
            r = requests.get(EUTILS + path, params=params, timeout=60)
            if r.status_code == 200:
                return r
        except Exception:
            pass
        time.sleep(1.2 * (i + 1))
    return None


def esearch(term, retmax=300):
    r = get("/esearch.fcgi", {"db": "gds", "term": term, "retmax": retmax,
                              "retmode": "json"})
    if not r:
        return []
    return r.json()["esearchresult"].get("idlist", [])


def esummary(ids):
    out = []
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        r = get("/esummary.fcgi", {"db": "gds", "id": ",".join(chunk),
                                   "retmode": "json"})
        if not r:
            continue
        j = r.json()["result"]
        for uid in j["uids"]:
            e = j[uid]
            acc = e.get("accession", "")
            if not acc.startswith("GSE"):
                continue
            out.append({
                "gse": acc,
                "gds_type": e.get("gdstype", ""),
                "title": e.get("title", ""),
                "summary": e.get("summary", "")[:300],
                "n_samples": e.get("n_samples", ""),
                "gpl": e.get("gpl", ""),
                "platform": e.get("platform", ""),
                "taxon": e.get("taxon", ""),
                "pdat": e.get("pdat", ""),
                "suppfile": e.get("suppfile", ""),
            })
        time.sleep(0.4)
    return out


def score_cp(row):
    txt = (row["title"] + " " + row["summary"]).lower()
    s = 0
    s += 3 * txt.count("chronic pancreatitis")
    s += txt.count("pancreatitis")
    s += 2 * ("normal" in txt or "healthy" in txt or "control" in txt)
    s += 2 * ("expression profiling by array" in row["gds_type"].lower()
              or "expression profiling by high throughput" in row["gds_type"].lower())
    s -= 4 * ("cell line" in txt or "in vitro" in txt or "stellate" in txt)
    s -= 2 * ("mouse" in txt or "rat " in txt)
    try:
        n = int(row["n_samples"])
        s += 1 if n >= 20 else 0
    except Exception:
        pass
    return s


print("=== CP candidate search ===")
terms = [
    '"chronic pancreatitis"[All Fields] AND "Homo sapiens"[Organism] AND '
    '("expression profiling by array"[DataSet Type] OR "expression profiling by high throughput sequencing"[DataSet Type])',
    'pancreatitis[All Fields] AND "Homo sapiens"[Organism] AND "expression profiling by array"[DataSet Type]',
]
ids = []
for t in terms:
    got = esearch(t)
    print("term hits:", len(got), "|", t[:80])
    ids.extend(got)
ids = sorted(set(ids), key=ids.index)
rows = esummary(ids)
cp = pd.DataFrame(rows)
if len(cp):
    cp["score"] = cp.apply(score_cp, axis=1)
    cp = cp.sort_values("score", ascending=False)
    cp.to_csv(os.path.join(RES, "geo_cp_candidates.csv"), index=False)
    show = cp.head(20)
    for _, r in show.iterrows():
        print(f"\n{r['gse']}  score={r['score']}  n={r['n_samples']}  GPL{r['gpl']}  {r['pdat']}")
        print("  ", r["title"][:120])
        print("  ", r["summary"][:220].replace("\n", " "))

print("\n\n=== PDAC candidate verification ===")
pdac_acc = ["GSE62452", "GSE28735", "GSE15471", "GSE16515", "GSE71729"]
prows = []
for acc in pdac_acc:
    got = esearch(f"{acc}[Accession]")
    rows = esummary(got)
    prows.extend(rows)
    time.sleep(0.3)
pp = pd.DataFrame(prows).drop_duplicates("gse")
pp.to_csv(os.path.join(RES, "geo_pdac_candidates.csv"), index=False)
for _, r in pp.iterrows():
    print(f"\n{r['gse']}  n={r['n_samples']}  GPL{r['gpl']}  {r['gds_type']}")
    print("  ", r["title"][:120])
    print("  ", r["summary"][:220].replace("\n", " "))
print("\nDONE")
