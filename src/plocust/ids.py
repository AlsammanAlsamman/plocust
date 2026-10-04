"""Build-independent identifiers derived from sequence.

The digest follows the GA4GH VRS convention (sha512, truncated to 24 bytes,
base64url): identical sequence gives an identical ID, whatever the genome
build or coordinates.
"""

import base64
import hashlib
import re

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
_VALID_SEQ = re.compile(r"^[ACGTN]+$")


def normalize_sequence(seq: str) -> str:
    seq = seq.strip().upper()
    if not _VALID_SEQ.match(seq):
        raise ValueError("sequence must contain only A, C, G, T, N")
    return seq


def reverse_complement(seq: str) -> str:
    return normalize_sequence(seq).translate(_COMPLEMENT)[::-1]


def canonical_sequence(seq: str) -> str:
    """Strand-independent form: the smaller of the sequence and its reverse complement."""
    seq = normalize_sequence(seq)
    return min(seq, reverse_complement(seq))


def sha512t24u(blob: bytes) -> str:
    digest = hashlib.sha512(blob).digest()[:24]
    return base64.urlsafe_b64encode(digest).decode("ascii")


def species_code(species: str) -> str:
    """'Oryza sativa' -> 'Os'."""
    parts = species.split()
    if len(parts) < 2:
        raise ValueError(f"expected a binomial species name, got {species!r}")
    return parts[0][0].upper() + parts[1][0].lower()


def passport_id(species: str, lead_anchor_sequence: str) -> str:
    """Stable passport ID, e.g. 'OsLP.<digest>'.

    Built from the species and the lead-SNP anchor only, so it never changes
    when the passport is placed on a new build or its description is updated.
    """
    blob = f"{species.strip()}|{canonical_sequence(lead_anchor_sequence)}".encode()
    return f"{species_code(species)}LP.{sha512t24u(blob)}"
