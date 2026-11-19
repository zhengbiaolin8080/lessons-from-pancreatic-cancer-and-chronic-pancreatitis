"""Step 08: ligand panel (18 focus-class metabolites) with isomeric SMILES
(catalog -> PubChem fallback) and RDKit Lipinski/ADME descriptors.
"""
import os
import time
import requests
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors, QED

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RES = os.path.join(ROOT, "results")

cat = pd.read_csv(os.path.join(DATA, "metabolites.csv"))
cat_map = {n.lower(): s for n, s in zip(cat["name"], cat["smiles"])}

# name, focus class, catalog name (or None -> PubChem)
PANEL = [
    ("Indole", "indole_trp", "Indole"),
    ("Indole-3-acetic acid", "indole_trp", "Indole-3-acetic acid"),
    ("Indole-3-propionic acid", "indole_trp", "3-Indolepropionic acid"),
    ("Indole-3-carboxaldehyde", "indole_trp", "Indole-3-carboxaldehyde"),
    ("Indole-3-acrylic acid", "indole_trp", "3-Indoleacrylic acid"),
    ("Indoxyl sulfate", "indole_trp", "Indoxyl sulfate"),
    ("Skatole", "indole_trp", "Skatole"),
    ("Tryptamine", "indole_trp", "Tryptamine"),
    ("Kynurenine", "indole_trp", "Kynurenine"),
    ("Acetic acid", "scfa", "Acetate"),
    ("Propionic acid", "scfa", "Propionate"),
    ("Butyric acid", "scfa", "Butyrate"),
    ("Deoxycholic acid", "bile_acid", "Deoxycholic Acid"),
    ("Lithocholic acid", "bile_acid", "Lithocholic acid"),
    ("Chenodeoxycholic acid", "bile_acid", "Chenodiol"),
    ("Ursodeoxycholic acid", "bile_acid", "Ursodiol"),
    ("Putrescine", "polyamine", "Putrescine"),
    ("Spermidine", "polyamine", "Spermidine"),
    ("Spermine", "polyamine", "Spermine"),
]


def pubchem_smiles(name):
    try:
        u = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
             f"{requests.utils.quote(name)}/property/IsomericSMILES,ConnectivitySMILES,SMILES/JSON")
        j = requests.get(u, timeout=60).json()
        p = j["PropertyTable"]["Properties"][0]
        return p.get("IsomericSMILES") or p.get("ConnectivitySMILES") or p.get("SMILES")
    except Exception as e:
        print("  PubChem FAIL", name, str(e)[:100])
        return None


rows = []
for name, cls, catname in PANEL:
    smi = cat_map.get((catname or "").lower()) if catname else None
    src = "gutMGene_catalog" if smi else "PubChem"
    if not smi:
        smi = pubchem_smiles(catname or name)
        time.sleep(0.3)
    if not smi:
        print("  !! no SMILES:", name)
        continue
    rows.append({"name": name, "class": cls, "isomeric_smiles": smi, "source": src})
    print(f"  {name:28s} {src:10s} {smi[:60]}")

lig = pd.DataFrame(rows)
lig.to_csv(os.path.join(DATA, "ligand_smiles.csv"), index=False)

# ---- ADME descriptors ----------------------------------------------------------
out = []
for _, r in lig.iterrows():
    m = Chem.MolFromSmiles(r["isomeric_smiles"])
    mw = Descriptors.MolWt(m)
    logp = Crippen.MolLogP(m)
    hbd = Lipinski.NumHDonors(m)
    hba = Lipinski.NumHAcceptors(m)
    tpsa = rdMolDescriptors.CalcTPSA(m)
    rotb = Lipinski.NumRotatableBonds(m)
    viol = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
    out.append({"name": r["name"], "class": r["class"], "MW": round(mw, 2),
                "cLogP": round(logp, 2), "HBD": hbd, "HBA": hba,
                "TPSA": round(tpsa, 2), "rotatable_bonds": rotb,
                "Lipinski_violations": viol, "QED": round(QED.qed(m), 3)})
adme = pd.DataFrame(out)
adme.to_csv(os.path.join(RES, "adme_lipinski.csv"), index=False)
print("\n", adme.to_string(index=False))
print("DONE")
