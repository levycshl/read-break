#!/usr/bin/env python3
"""
Stream-parse MiSeq FASTQs with read-break and collect a small summary table.

Usage
-----
python test_driver.py \
    --run-dir  Z:/mnt/wigstore3/data/checkpointcharlie2/inessa/250707_M06142_0497_000000000-M2RFJ/Alignment_1/20250708_193955/Fastq \
    --config   parsers/a_deaminate_2025_07_10.yaml \
    --prefix   DB \
    --max-reads 10000
"""

from __future__ import annotations
import argparse, re, sys
from pathlib import Path
from collections import defaultdict

import yaml
import pandas as pd
from read_break.parser import ReadParser
from read_break.io import FastqReader
from read_break.logic import flatten_dot

# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

PAT = re.compile(
    r'^(?P<key>[^_]+_[^_]+)_L\d{3}_R(?P<read>[12])_\d{3}\.fastq\.gz$'
)

def discover_libraries(run_dir: Path, sample_prefix: str = "DB") -> dict[str, dict[str, Path]]:
    """Return {'DB1_S3': {'R1': path, 'R2': path}, …} limited by prefix."""
    libs: dict[str, dict[str, Path]] = defaultdict(lambda: {"R1": None, "R2": None})
    for f in run_dir.glob("*.fastq.gz"):
        m = PAT.match(f.name)
        if not m:
            continue
        key, read = m.group("key", "read")
        libs[key][f"R{read}"] = f.resolve()

    # keep only complete pairs that match the prefix
    return {k: v for k, v in libs.items()
            if k.startswith(sample_prefix) and None not in v.values()}


def build_parser(config_path: Path, r1: Path, r2: Path) -> ReadParser:
    """Instantiate ReadParser exactly like the CLI does."""
    return ReadParser.from_cli(
        config_path=config_path,
        r1_path=r1,
        r2_path=r2,
    )

# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main(argv: list[str] | None = None) -> None:
    parser_cli = argparse.ArgumentParser(description="run read-break on MiSeq libraries")
    parser_cli.add_argument("--run-dir",  type=Path, required=True, help="directory with *.fastq.gz")
    parser_cli.add_argument("--config",   type=Path, required=True, help="YAML parser spec")
    parser_cli.add_argument("--prefix",   default="DB", help="sample prefix to keep")
    parser_cli.add_argument("--max-reads",type=int, default=10_000, help="cap for quick tests")
    args = parser_cli.parse_args(argv)

    libs = discover_libraries(args.run_dir, args.prefix)
    if not libs:
        sys.exit("No matching libraries found.")

    # pick one library for the quick sanity run
    key = sorted(libs)[0]
    r1, r2 = libs[key]["R1"], libs[key]["R2"]
    print(f"⇢ Testing parser on {key}\n   R1={r1.name}\n   R2={r2.name}")

    parser = build_parser(args.config, r1, r2)

    results = []
    with FastqReader(r1, r2) as reader:
        for idx, (rec1, rec2) in enumerate(reader, 1):
            if idx > args.max_reads:
                break

            ctx = parser.parse(rec1, rec2)
            if ctx and ctx.get("status") == "ok":
                results.append(ctx)

            if idx % 1_000 == 0 or idx == 1:
                flat = flatten_dot(parser.get_parse_log())
                if idx == 1:
                    print("\t".join(flat))
                else:
                    print("\t".join(map(str, flat.values())))

    # ---- summary -----------------------------------------------------------
    df = pd.DataFrame(results)
    print(f"\nParsed {len(df):,} passing reads (of {args.max_reads:,} streamed).")
    # df.to_csv("parsed_subset.csv", index=False)  # optional

# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    main()          # <-- invoked only when run as a script
