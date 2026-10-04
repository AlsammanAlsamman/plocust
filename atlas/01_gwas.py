"""Atlas step 1: mixed-model GWAS (REGENIE) for every usable trait of the 3K and RDP1 panels.

3K   : IRGCIS descriptors from the 3K phenotype workbook. Ordinal codes are analysed as quantitative;
       colour descriptors become binary (colourless vs coloured) using the workbook's own dictionary.
RDP1 : the 34 traits of Zhao et al. 2011, re-analysed for effect signs (MSU6 coordinates).
Outputs go to data/atlas/gwas/<panel>/<trait>.regenie and a trait table in atlas/results/traits.tsv.
"""

import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
WORK = DATA / "atlas" / "gwas"
RES = ROOT / "atlas" / "results"
REGENIE = DATA / "bin/regenie_v4.1.3.gz_x86_64_Linux"
PLINK = DATA / "bin/plink2"
MIN_N = 150
MIN_MINOR = 30  # binary traits: at least this many cases and controls

# 3K descriptor -> (short name, kind). kind: "q" ordinal/quantitative, "colour" binary colourless vs coloured,
# "present" binary absent (0) vs present (>0), "binary12" codes 1/2 -> 0/1
K3_TRAITS = {
    "CULT_CODE_REPRO": ("culm_length", "q"), "CUDI_CODE_REPRO": ("culm_diameter", "q"),
    "CUNO_CODE_REPRO": ("culm_number", "q"), "LLT_CODE": ("leaf_length", "q"), "PLT_CODE_POST": ("panicle_length", "q"),
    "SDHT_CODE": ("seedling_height", "q"), "SLLT_CODE": ("sterile_lemma_length", "q"), "LA": ("leaf_angle", "q"),
    "FLA_REPRO": ("flag_leaf_angle", "q"), "CUAN_REPRO": ("culm_angle", "q"), "PEX_REPRO": ("panicle_exsertion", "q"),
    "SECOND_BR_REPRO": ("secondary_branching", "q"), "LSEN": ("leaf_senescence", "q"), "SPKF": ("spikelet_fertility", "q"),
    "PTH": ("panicle_threshability", "q"), "PSH": ("panicle_shattering", "q"), "PTY": ("panicle_type", "q"),
    "LPPUB": ("lemma_pubescence", "q"), "BLPUB_VEG": ("blade_pubescence", "q"),
    "AWPR_REPRO": ("awn", "present"), "ENDO": ("glutinous_endosperm", "binary12"),
    "SCCO_REV": ("seed_coat_colour", "colour"), "APCO_REV_REPRO": ("apiculus_colour", "colour"),
    "LPCO_REV_POST": ("hull_colour", "colour"), "BLSCO_REV_VEG": ("leaf_sheath_colour", "colour"),
    "LIGCO_REV_VEG": ("ligule_colour", "colour"), "NOCO_REV": ("node_colour", "colour"),
    "INCO_REV_REPRO": ("internode_colour", "colour"), "AUCO_REV_VEG": ("auricle_colour", "colour"),
    "CCO_REV_VEG": ("collar_colour", "colour"), "SLCO_REV": ("sterile_lemma_colour", "colour"),
    "STCO_REV": ("stigma_colour", "colour"),
}
COLOURLESS = re.compile(r"WHITE|GREEN|STRAW|COLOURLESS|COLORLESS|NONE|ABSENT|LIGHT GREEN|PALE", re.I)


def run(cmd, log):
    with open(log, "w") as fh:
        subprocess.run([str(c) for c in cmd], check=True, stdout=fh, stderr=subprocess.STDOUT)


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def pheno_3k() -> tuple[pd.DataFrame, list[dict]]:
    book = DATA / "3k/3kRG_PhenotypeData_v20170411.xlsx"
    p = pd.read_excel(book, "Phenotype Data")
    dic = pd.read_excel(book, "<2007 Dictionary").ffill()
    sra = pd.read_csv(DATA / "3k/3K_list_sra_ids.txt", sep="\t")
    sra["key"] = sra["Genetic_Stock_varname"].map(norm)
    p["key"] = p["NAME"].map(norm)
    p = p.merge(sra[["key", "3K_DNA_IRIS_UNIQUE_ID"]], on="key").drop_duplicates("3K_DNA_IRIS_UNIQUE_ID")
    out = pd.DataFrame({"FID": p["3K_DNA_IRIS_UNIQUE_ID"], "IID": p["3K_DNA_IRIS_UNIQUE_ID"]})
    meta = []
    for col, (name, kind) in K3_TRAITS.items():
        if col not in p:
            continue
        v = pd.to_numeric(p[col], errors="coerce").where(lambda x: x != 999)
        if kind == "present":
            v = (v > 0).astype(float).where(v.notna())
        elif kind == "binary12":
            v = v.where(v.isin([1, 2])).map({1: 0.0, 2: 1.0})
        elif kind == "colour":
            d = dic[dic["Field Name"] == col]
            colourless = {float(c) for c, desc in zip(d["Value"], d["Value Description"])
                          if isinstance(desc, str) and COLOURLESS.search(desc) and str(c).replace(".", "").isdigit()}
            if not colourless:
                continue
            v = (~v.isin(colourless)).astype(float).where(v.notna())
        binary = kind != "q"
        if not binary and v.nunique() == 2:  # a two-code descriptor is binary
            v = (v == v.max()).astype(float).where(v.notna())
            binary = True
        n = int(v.notna().sum())
        minor = int(min((v == 0).sum(), (v == 1).sum())) if binary else None
        if n < MIN_N or (binary and minor < MIN_MINOR) or v.nunique() < 2:
            continue
        out[name] = v.values
        meta.append({"panel": "3K", "trait": name, "source_column": col, "binary": binary, "n": n, "minor_class": minor})
    return out, meta


def pheno_rdp1() -> tuple[pd.DataFrame, list[dict]]:
    raw = pd.read_csv(DATA / "rdp1/pheno_rdp1_raw.tsv", sep="\t")
    raw.columns = [c.strip() for c in raw.columns]
    out = pd.DataFrame({"FID": raw["HybID"], "IID": raw["NSFTVID"].astype(str)})
    meta = []
    for col in raw.columns[2:]:
        v = pd.to_numeric(raw[col], errors="coerce")
        name = re.sub(r"[^a-z0-9]+", "_", col.lower()).strip("_")
        if v.nunique() == 2:
            v = (v == v.max()).astype(float).where(v.notna())
        binary = set(v.dropna().unique()) <= {0, 1}
        n = int(v.notna().sum())
        minor = int(min((v == 0).sum(), (v == 1).sum())) if binary else None
        if n < MIN_N or (binary and minor < MIN_MINOR) or v.nunique() < 2:
            continue
        out[name] = v.values
        meta.append({"panel": "RDP1", "trait": name, "source_column": col, "binary": binary, "n": n, "minor_class": minor})
    return out, meta


def regenie(panel: str, bed: Path, step1: Path, covar: Path, pheno: pd.DataFrame, meta: list[dict]):
    d = WORK / panel
    d.mkdir(parents=True, exist_ok=True)
    pf = d / "pheno.tsv"
    pheno.to_csv(pf, sep="\t", index=False, na_rep="NA")
    for binary in (False, True):
        cols = [m["trait"] for m in meta if m["panel"] == panel and m["binary"] == binary]
        if not cols:
            continue
        tag = "bt" if binary else "qt"
        bt = ["--bt"] if binary else []
        common = ["--phenoFile", pf, "--phenoColList", ",".join(cols), "--covarFile", covar, "--threads", "8"]
        run([REGENIE, "--step", "1", "--bed", step1, "--bsize", "1000", "--lowmem", "--lowmem-prefix", d / f"tmp_{tag}",
             "--out", d / f"s1_{tag}", *bt, *common], d / f"s1_{tag}.log")
        firth = ["--firth", "--approx", "--pThresh", "0.01"] if binary else []
        run([REGENIE, "--step", "2", "--bed", bed, "--pred", d / f"s1_{tag}_pred.list", "--bsize", "400",
             "--minMAC", "20", "--out", d / "gwas", *bt, *firth, *common], d / f"s2_{tag}.log")
        print(f"{panel}: {tag} done ({len(cols)} traits)", flush=True)


def main():
    RES.mkdir(parents=True, exist_ok=True)
    p3, m3 = pheno_3k()
    pr, mr = pheno_rdp1()
    meta = pd.DataFrame(m3 + mr)
    meta.to_csv(RES / "traits.tsv", sep="\t", index=False)
    print(meta.groupby(["panel", "binary"]).size().to_string(), flush=True)

    rdp = DATA / "rdp1"
    regenie("3K", DATA / "3k/core3k", DATA / "3k/step1set", DATA / "3k/covar5.tsv", p3, m3)
    regenie("RDP1", rdp / "rdp1", rdp / "rstep1", rdp / "rcovar.tsv", pr, mr)

    # inflation per trait
    from scipy.stats import chi2

    rows = []
    for m in meta.itertuples():
        f = WORK / m.panel / f"gwas_{m.trait}.regenie"
        if not f.exists():
            continue
        lp = pd.read_csv(f, sep=" ", usecols=["LOG10P"])["LOG10P"].dropna().to_numpy()
        lam = np.median(chi2.isf(np.minimum(10.0 ** -lp, 1), 1)) / 0.4549
        rows.append({"panel": m.panel, "trait": m.trait, "lambda": round(lam, 3), "max_log10p": round(lp.max(), 1),
                     "n_bonferroni": int((lp > -np.log10(0.05 / len(lp))).sum())})
    q = pd.DataFrame(rows)
    q.to_csv(RES / "gwas_qc.tsv", sep="\t", index=False)
    print(q.to_string(index=False))


if __name__ == "__main__":
    main()
