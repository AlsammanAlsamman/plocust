"""The known-loci database: a versioned SQLite file, downloaded separately from the package.

Tables:
    meta(key, value)
    passports(passport_id, kind, species, trait, ontology_id, study, traits, json)
    placements(passport_id, build, chrom, start, end, lead_pos, method)
"""

from __future__ import annotations

import gzip
import hashlib
import sqlite3
import urllib.request
from datetime import date
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from . import __version__
from .io import Genome, normalize_chrom
from .passport import SCHEMA_VERSION, Anchor, AnchorRole, Gene, LocusKind, LocusPassport, Placement, Source, Trait

# Released databases. A release is a gzip-compressed SQLite file with its sha256.
REGISTRY: dict[str, dict] = {
    # "rice": {"url": "https://zenodo.org/records/<id>/files/plocust-db-rice-<version>.sqlite.gz", "sha256": "..."},
}

_SCHEMA = """
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE passports (
    passport_id TEXT PRIMARY KEY, kind TEXT, species TEXT, trait TEXT, ontology_id TEXT,
    study TEXT, traits TEXT, json TEXT NOT NULL);
CREATE TABLE placements (
    passport_id TEXT REFERENCES passports(passport_id), build TEXT, chrom TEXT,
    start INTEGER, "end" INTEGER, lead_pos INTEGER, method TEXT);
CREATE INDEX placements_region ON placements(build, chrom, start, "end");
CREATE INDEX passports_trait ON passports(trait);
"""


def _all_traits(p: LocusPassport) -> list[str]:
    return [p.trait.name] + [t.name for t in p.traits_other]


def create(path: str | Path, passports: Iterable[LocusPassport], name: str, version: str,
           extra_meta: Optional[dict] = None) -> int:
    """Write a new database file; returns the number of passports stored."""
    path = Path(path)
    if path.exists():
        path.unlink()
    con = sqlite3.connect(path)
    con.executescript(_SCHEMA)
    n = 0
    seen = set()
    for p in passports:
        if p.passport_id in seen:
            continue
        seen.add(p.passport_id)
        con.execute("INSERT INTO passports VALUES (?,?,?,?,?,?,?,?)",
                    (p.passport_id, p.kind.value, p.species, p.trait.name, p.trait.ontology_id, p.source.study,
                     "|".join(_all_traits(p)), p.model_dump_json(exclude_none=True)))
        for pl in p.placements:
            con.execute("INSERT INTO placements VALUES (?,?,?,?,?,?,?)",
                        (p.passport_id, pl.build, pl.chrom, pl.start, pl.end, pl.lead_pos, pl.method))
        n += 1
    meta = {"name": name, "version": version, "created": date.today().isoformat(), "schema_version": SCHEMA_VERSION,
            "plocust_version": __version__, "n_passports": str(n), **(extra_meta or {})}
    con.executemany("INSERT INTO meta VALUES (?,?)", meta.items())
    con.commit()
    con.close()
    return n


class LocusDB:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"database {path} not found (see `plocust db download` / `plocust db build-genes`)")
        self._con = sqlite3.connect(f"file:{self.path}?mode=ro", uri=True)

    def info(self) -> dict:
        meta = dict(self._con.execute("SELECT key, value FROM meta"))
        meta["builds"] = ",".join(r[0] for r in self._con.execute("SELECT DISTINCT build FROM placements"))
        meta["kinds"] = ",".join(f"{k}:{n}" for k, n in self._con.execute(
            "SELECT kind, count(*) FROM passports GROUP BY kind"))
        return meta

    def _load(self, ids) -> list[LocusPassport]:
        ids = list(dict.fromkeys(ids))
        if not ids:
            return []
        q = f"SELECT json FROM passports WHERE passport_id IN ({','.join('?' * len(ids))})"
        return [LocusPassport.model_validate_json(r[0]) for r in self._con.execute(q, ids)]

    def get(self, passport_id: str) -> Optional[LocusPassport]:
        found = self._load([passport_id])
        return found[0] if found else None

    def region(self, build: str, chrom, start: int, end: int) -> list[LocusPassport]:
        rows = self._con.execute('SELECT passport_id FROM placements WHERE build=? AND chrom=? AND "end">=? AND start<=?',
                                 (build, normalize_chrom(chrom), start, end))
        return self._load(r[0] for r in rows)

    def near(self, passports: list[LocusPassport], build: str, pad_kb: float = 1000) -> list[LocusPassport]:
        """Database records within pad_kb of any of the given passports on `build`."""
        pad = int(pad_kb * 1000)
        out = {}
        for p in passports:
            pl = p.placement(build)
            if pl is not None:
                for r in self.region(build, pl.chrom, pl.start - pad, pl.end + pad):
                    out[r.passport_id] = r
        return list(out.values())

    def by_trait(self, keyword: str) -> list[LocusPassport]:
        rows = self._con.execute("SELECT passport_id FROM passports WHERE lower(traits) LIKE ?", (f"%{keyword.lower()}%",))
        return self._load(r[0] for r in rows)


def download(name: str, dest_dir: str | Path) -> Path:
    if name not in REGISTRY:
        raise KeyError(f"no released database named {name!r}; available: {sorted(REGISTRY) or 'none yet'}. "
                       "Build one locally with `plocust db build-genes`.")
    entry = REGISTRY[name]
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    gz = dest_dir / Path(entry["url"]).name
    urllib.request.urlretrieve(entry["url"], gz)
    if hashlib.sha256(gz.read_bytes()).hexdigest() != entry["sha256"]:
        gz.unlink()
        raise ValueError("checksum mismatch; download removed")
    out = gz.with_suffix("")
    with gzip.open(gz) as src, open(out, "wb") as dst:
        dst.write(src.read())
    return out


# --------------------------------------------------------------- known genes


def gene_passports(genes: pd.DataFrame, gene_table: pd.DataFrame, genome: Genome, species: str, study: str,
                   keywords: Optional[pd.DataFrame] = None, flank: int = 100, id_col: str = "RAPdb",
                   symbol_col: str = "Symbol") -> list[LocusPassport]:
    """Passports for cloned / known genes.

    Anchors: `gene` (primary, at the gene start), `block_left` / `block_right`
    at the gene start and end. Traits come from the keyword table if given.
    """
    kw = {}
    if keywords is not None:
        for gid, k in keywords.groupby(id_col)["Keyword"]:
            kw[gid] = list(dict.fromkeys(str(x) for x in k if isinstance(x, str)))
    loc = genes.set_index("id")
    out = []
    for _, row in gene_table.iterrows():
        gid = row[id_col]
        if not isinstance(gid, str) or gid not in loc.index:
            continue
        g = loc.loc[gid]
        if isinstance(g, pd.DataFrame):
            g = g.iloc[0]
        chrom, start, end = g["chrom"], int(g["start"]), int(g["end"])
        five = start if g["strand"] == "+" else end
        # primary anchor just inside the gene's 5' end
        centre = five + flank if g["strand"] == "+" else five - flank
        anchors = []
        for role, pos in ((AnchorRole.gene, centre), (AnchorRole.block_left, start), (AnchorRole.block_right, end)):
            seq, off = genome.flank(chrom, pos, flank)
            anchors.append(Anchor(role=role, sequence=seq, variant_offset=off, label=gid))
        traits = kw.get(gid, [])
        symbol = str(row[symbol_col]).split("|")[0] if isinstance(row.get(symbol_col), str) else gid
        out.append(LocusPassport(
            kind=LocusKind.gene,
            species=species,
            trait=Trait(name=traits[0] if traits else "cloned gene"),
            traits_other=[Trait(name=t) for t in traits[1:]],
            source=Source(study=study),
            anchors=anchors,
            placements=[Placement(build=genome.build, chrom=chrom, start=start, end=end, lead_pos=centre,
                                  strand=g["strand"])],
            genes=[Gene(id=gid, name=symbol, start=start, end=end, distance_bp=0,
                        annotation=row[symbol_col] if isinstance(row.get(symbol_col), str) else None)],
        ))
    return out


def gene_symbol(p: LocusPassport) -> Optional[str]:
    if p.kind is not LocusKind.gene or not p.genes:
        return None
    return p.genes[0].name
