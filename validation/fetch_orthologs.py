"""Ground truth for validation B: Ensembl Compara orthologs of a random sample of cloned rice genes.

Writes data/known/orthologs.tsv (rap_id, species, target_gene, type). Cached; re-runnable.
"""

import csv
import json
import time
import urllib.request
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parents[1] / "data"
SPECIES = ["oryza_sativa_mh63", "oryza_sativa_zs97", "oryza_sativa_n22", "oryza_sativa_azucena"]
N_GENES = 300
URL = "https://rest.ensembl.org/homology/id/oryza_sativa/{}?compara=plants;type=orthologues;format=condensed"


def main():
    genes = pd.read_csv(DATA / "known/geneInfo.table.txt", sep="\t", quoting=csv.QUOTE_NONE,
                        encoding_errors="replace", dtype=str)
    ids = genes["RAPdb"].dropna()
    ids = ids[ids.str.match(r"^Os\d\dg\d{7}$")].drop_duplicates().sample(N_GENES, random_state=42)
    cache = DATA / "known/orthologs_cache.json"
    done = json.loads(cache.read_text()) if cache.exists() else {}
    for i, gid in enumerate(ids):
        if gid in done:
            continue
        req = urllib.request.Request(URL.format(gid), headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                done[gid] = json.load(r)["data"]
        except Exception as e:  # unknown IDs return 400
            done[gid] = {"error": str(e)}
        if i % 25 == 0:
            cache.write_text(json.dumps(done))
        time.sleep(0.08)
    cache.write_text(json.dumps(done))
    rows = []
    for gid, data in done.items():
        if not isinstance(data, list):
            continue
        for entry in data:
            for h in entry["homologies"]:
                if h["species"] in SPECIES:
                    rows.append((gid, h["species"], h["id"], h["type"]))
    out = pd.DataFrame(rows, columns=["rap_id", "species", "target_gene", "type"])
    out.to_csv(DATA / "known/orthologs.tsv", sep="\t", index=False)
    print(f"{out['rap_id'].nunique()} genes with orthologs; by species:\n{out.groupby(['species', 'type']).size()}")


if __name__ == "__main__":
    main()
