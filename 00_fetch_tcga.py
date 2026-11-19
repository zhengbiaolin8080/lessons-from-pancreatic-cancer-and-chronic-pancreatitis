"""Step 00: fetch TCGA-PAAD expression + phenotype + PanCanAtlas survival (UCSC Xena).

HiSeqV2.gz is gene-level log2(norm_count+1) RSEM, rows = Hugo symbols.
Survival table = Liu et al. 2018 Cell PanCanAtlas clinical supplement
(OS/DSS/DFI/PFI per sample), the canonical TCGA outcome resource.
"""
import os
import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "TCGA")
os.makedirs(DATA, exist_ok=True)

FILES = {
    "PAAD_HiSeqV2.gz": "https://tcga.xenahubs.net/download/TCGA.PAAD.sampleMap/HiSeqV2.gz",
    "PAAD_clinicalMatrix": "https://tcga.xenahubs.net/download/TCGA.PAAD.sampleMap/PAAD_clinicalMatrix",
    "TCGA_survival_pancanatlas.tsv":
        "https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp",
}

for name, url in FILES.items():
    dst = os.path.join(DATA, name)
    if os.path.exists(dst) and os.path.getsize(dst) > 1000:
        print("exists:", name, os.path.getsize(dst) >> 20, "MB")
        continue
    print("downloading", url, flush=True)
    with requests.get(url, stream=True, timeout=900) as r:
        r.raise_for_status()
        with open(dst + ".part", "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    os.replace(dst + ".part", dst)
    print("saved", name, os.path.getsize(dst) >> 20, "MB", flush=True)

print("DONE")
