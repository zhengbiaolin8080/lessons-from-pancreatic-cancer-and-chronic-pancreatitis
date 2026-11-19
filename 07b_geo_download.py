"""Step 07b: download GEO series matrices + platform annotations and inspect
sample composition for the validation datasets.

Datasets:
  GSE62452  PDAC vs adjacent normal (69T/61N), GPL6244 HuGene-1_0-st
  GSE28735  PDAC vs matched adjacent normal (45 pairs), GPL6244
  GSE77858  PDAC / normal / pancreatitis, GPL7264 (Agilent)
  GSE101462 PDAC vs chronic pancreatitis biomarker cohort, GPL10558 (Illumina)
"""
import os
import gzip
import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEO = os.path.join(ROOT, "data", "GEO")
os.makedirs(GEO, exist_ok=True)

GSES = ["GSE62452", "GSE28735", "GSE77858", "GSE101462"]
GPLS = ["GPL6244", "GPL570", "GPL7264", "GPL10558"]


def fetch(url, dst):
    if os.path.exists(dst) and os.path.getsize(dst) > 500:
        print("exists:", os.path.basename(dst), os.path.getsize(dst) // 1024, "KB")
        return True
    try:
        with requests.get(url, stream=True, timeout=600) as r:
            if r.status_code != 200:
                print("MISS", r.status_code, url)
                return False
            with open(dst + ".part", "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
        os.replace(dst + ".part", dst)
        print("saved:", os.path.basename(dst), os.path.getsize(dst) // 1024, "KB")
        return True
    except Exception as e:
        print("ERR", url, e)
        return False


def gse_url(gse):
    stub = gse[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/{gse}_series_matrix.txt.gz"


def gpl_url(gpl):
    stub = gpl[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{stub}/{gpl}/annot/{gpl}.annot.gz"


for gse in GSES:
    fetch(gse_url(gse), os.path.join(GEO, f"{gse}_series_matrix.txt.gz"))
for gpl in GPLS:
    fetch(gpl_url(gpl), os.path.join(GEO, f"{gpl}.annot.gz"))

# ---- inspect sample metadata -------------------------------------------------
for gse in GSES:
    f = os.path.join(GEO, f"{gse}_series_matrix.txt.gz")
    if not os.path.exists(f):
        continue
    print("\n" + "=" * 90 + f"\n{gse}")
    meta = {}
    with gzip.open(f, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                key, _, rest = line.partition("\t")
                vals = [v.strip().strip('"') for v in rest.rstrip("\n").split("\t")]
                meta[key] = vals
    ids = meta.get("!Sample_geo_accession", [])
    titles = meta.get("!Sample_title", [])
    src = meta.get("!Sample_source_name_ch1", [])
    chars = [v for k, v in meta.items() if k.startswith("!Sample_characteristics_ch1")]
    print("n_samples:", len(ids))
    for i in range(len(ids)):
        c = " | ".join(x[i] for x in chars) if chars else ""
        print(f"  {ids[i]:12s} {titles[i][:42]:42s} | {src[i][:38]:38s} | {c[:110]}")
print("\nDONE")
