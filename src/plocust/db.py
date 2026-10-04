"""The known-loci database: a versioned SQLite file, downloaded separately from the package.

Tables:
    meta(key, value)
    passports(record_id, passport_id, kind, species, trait, ontology_id, study, traits, json)
    placements(record_id, passport_id, build, chrom, start, end, lead_pos, method)

A record is one locus in one study for one trait. `passport_id` identifies the locus by its DNA, so the
same SNP associated with several traits or genes gives several records sharing one passport_id:
looking a passport_id up across records is a pleiotropy query.
    markers(marker_id, sequence, panel)   flanking sequence of SNP-chip / KASP markers, variant as [A/G]
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
    record_id TEXT PRIMARY KEY, passport_id TEXT NOT NULL, kind TEXT, species TEXT, trait TEXT, ontology_id TEXT,
    study TEXT, traits TEXT, json TEXT NOT NULL);
CREATE INDEX passports_id ON passports(passport_id);
CREATE TABLE placements (
    record_id TEXT REFERENCES passports(record_id), passport_id TEXT, build TEXT, chrom TEXT,
    start INTEGER, "end" INTEGER, lead_pos INTEGER, method TEXT);
CREATE INDEX placements_region ON placements(build, chrom, start, "end");
CREATE INDEX passports_trait ON passports(trait);
CREATE TABLE IF NOT EXISTS markers (marker_id TEXT PRIMARY KEY, sequence TEXT NOT NULL, panel TEXT);
"""


def record_id(p: LocusPassport) -> str:
    """Key of one record: the locus (passport_id) in one study, for one trait (and target gene, if any)."""
    from .ids import sha512t24u

    target = p.genes[0].id if p.kind is LocusKind.eqtl and p.genes else ""
    return p.passport_id + "." + sha512t24u(f"{p.source.study}|{p.trait.name}|{target}".encode())[:12]


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
        rid = record_id(p)
        if rid in seen:
            continue
        seen.add(rid)
        con.execute("INSERT INTO passports VALUES (?,?,?,?,?,?,?,?,?)",
                    (rid, p.passport_id, p.kind.value, p.species, p.trait.name, p.trait.ontology_id, p.source.study,
                     "|".join(_all_traits(p)), p.model_dump_json(exclude_none=True)))
        for pl in p.placements:
            con.execute("INSERT INTO placements VALUES (?,?,?,?,?,?,?,?)",
                        (rid, p.passport_id, pl.build, pl.chrom, pl.start, pl.end, pl.lead_pos, pl.method))
        n += 1
    meta = {"name": name, "version": version, "created": date.today().isoformat(), "schema_version": SCHEMA_VERSION,
            "plocust_version": __version__, "n_passports": str(n), **(extra_meta or {})}
    con.executemany("INSERT INTO meta VALUES (?,?)", meta.items())
    con.commit()
    con.close()
    return n


def add_markers(path: str | Path, markers: pd.DataFrame, panel: str, id_col: str = "marker_id",
                seq_col: str = "sequence") -> int:
    """Add marker flanking sequences (variant written as [A/G]) to an existing database."""
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE IF NOT EXISTS markers (marker_id TEXT PRIMARY KEY, sequence TEXT NOT NULL, panel TEXT)")
    rows = [(str(m), str(q), panel) for m, q in zip(markers[id_col], markers[seq_col]) if isinstance(q, str) and "[" in q]
    con.executemany("INSERT OR REPLACE INTO markers VALUES (?,?,?)", rows)
    con.commit()
    con.close()
    return len(rows)


def attach_ld_panel(path: str | Path, bfile: str | Path, build: str, name: str) -> None:
    """Register a reference genotype panel (PLINK prefix) for LD on `build`.

    Stored relative to the database file when it sits next to it, so database and panel can be
    downloaded and moved together.
    """
    path, bfile = Path(path).resolve(), Path(str(bfile).removesuffix(".bed")).resolve()
    try:
        ref = str(bfile.relative_to(path.parent))
    except ValueError:
        ref = str(bfile)
    con = sqlite3.connect(path)
    con.executemany("INSERT OR REPLACE INTO meta VALUES (?,?)",
                    [(f"ld_panel:{build}", ref), (f"ld_panel_name:{build}", name)])
    con.commit()
    con.close()


class LocusDB:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"database {path} not found (see `plocust db download` / `plocust db build-genes`)")
        self._con = sqlite3.connect(f"file:{self.path}?mode=ro", uri=True)
        cols = {r[1] for r in self._con.execute("PRAGMA table_info(passports)")}
        self._key = "record_id" if "record_id" in cols else "passport_id"  # databases built before records existed

    def info(self) -> dict:
        meta = dict(self._con.execute("SELECT key, value FROM meta"))
        meta["builds"] = ",".join(r[0] for r in self._con.execute("SELECT DISTINCT build FROM placements"))
        meta["kinds"] = ",".join(f"{k}:{n}" for k, n in self._con.execute(
            "SELECT kind, count(*) FROM passports GROUP BY kind"))
        return meta

    def _load(self, record_ids) -> list[LocusPassport]:
        ids = list(dict.fromkeys(record_ids))
        out = []
        for i in range(0, len(ids), 900):  # SQLite parameter limit
            chunk = ids[i:i + 900]
            q = f"SELECT json FROM passports WHERE {self._key} IN ({','.join('?' * len(chunk))})"
            out += [LocusPassport.model_validate_json(r[0]) for r in self._con.execute(q, chunk)]
        return out

    def records(self, passport_id: str) -> list[LocusPassport]:
        """Every record of one locus (all studies and traits): the pleiotropy view."""
        rows = self._con.execute(f"SELECT {self._key} FROM passports WHERE passport_id = ?", (passport_id,))
        return self._load(r[0] for r in rows)

    def get(self, passport_id: str) -> Optional[LocusPassport]:
        found = self.records(passport_id)
        return found[0] if found else None

    def region(self, build: str, chrom, start: int, end: int) -> list[LocusPassport]:
        rows = self._con.execute(f'SELECT {self._key} FROM placements WHERE build=? AND chrom=? AND "end">=? AND start<=?',
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
                    out[record_id(r)] = r
        return list(out.values())

    def ld_panel(self, build: str):
        """The reference LD panel registered for `build`, or None."""
        from .io import Genotypes

        row = self._con.execute("SELECT value FROM meta WHERE key = ?", (f"ld_panel:{build}",)).fetchone()
        if row is None:
            return None
        prefix = Path(row[0])
        if not prefix.is_absolute():
            prefix = self.path.parent / prefix
        return Genotypes.open(prefix)

    def marker_flanks(self, panel: Optional[str] = None):
        """Marker flanking sequences as a FlankTable (empty when the database has none)."""
        from .io import FlankTable

        try:
            q = "SELECT marker_id, sequence FROM markers" + (" WHERE panel = ?" if panel else "")
            rows = self._con.execute(q, (panel,) if panel else ()).fetchall()
        except sqlite3.OperationalError:
            rows = []
        return FlankTable(pd.DataFrame(rows, columns=["marker_id", "sequence"]), "marker_id", "sequence")

    def by_trait(self, keyword: str) -> list[LocusPassport]:
        rows = self._con.execute(f"SELECT {self._key} FROM passports WHERE lower(traits) LIKE ?", (f"%{keyword.lower()}%",))
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
            kw[gid] = k.dropna().astype(str).value_counts().index.tolist()  # most-cited keyword first
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
