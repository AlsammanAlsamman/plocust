"""Command-line interface: `plocust <command>`."""

import argparse
import json
import sys

from pydantic import ValidationError

from . import __version__
from .passport import LocusPassport, json_schema


def _cmd_schema(args: argparse.Namespace) -> int:
    print(json.dumps(json_schema(), indent=2))
    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    status = 0
    for path in args.passports:
        try:
            p = LocusPassport.from_json(path)
        except (ValidationError, ValueError, OSError) as e:
            print(f"FAIL  {path}\n{e}", file=sys.stderr)
            status = 1
        else:
            print(f"OK    {path}  {p.passport_id}")
    return status


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="plocust", description=__doc__)
    parser.add_argument("--version", action="version", version=f"plocust {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("schema", help="print the passport JSON Schema").set_defaults(func=_cmd_schema)

    p = sub.add_parser("validate", help="check passport JSON files against the schema")
    p.add_argument("passports", nargs="+")
    p.set_defaults(func=_cmd_validate)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
