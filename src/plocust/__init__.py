"""LOCUST: build-independent locus passports for plant GWAS."""

from .passport import SCHEMA_VERSION, LocusPassport

__version__ = "0.0.1"
__all__ = ["LocusPassport", "SCHEMA_VERSION", "__version__"]
