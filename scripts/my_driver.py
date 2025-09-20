#!/usr/bin/env python3
from __future__ import annotations
import argparse, re, sys
from pathlib import Path
from collections import defaultdict
import yaml
from read_break.parser import ReadParser
from read_break.io import FastqReader
from read_break.logic import flatten_dot
import pandas as pd

PAT = re.compile(r'^(?P<key>[^_]+_[^_]+)_L\d{3}_R(?P<read>[12])_\d{3}\.fastq\.gz$')

def discover(run_dir: Path, prefix: str) -> dict[str, dict[str, Path]]:
    libs = defaultdict(lambda: {"R1": None, "R2": None})
    for f in run_dir.glob("*.fastq.gz"):
        m = PAT.match(f.name)
        if m:
            key, read = m.group("key", "read")
            libs[key][f"R{read}"] = f.resolve()
    return {k: v for k, v in libs.items()
            if k.startswith(prefix) and None not in v.values()}

def build_parser(cfg_path: Path) -> ReadParser:
    cfg = yaml.safe_load(cfg_path.read_text())
    return ReadParser(cfg, cfg["params"])

def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir",  type=Path, required=True)
    ap.add_argument("--config",   type=Path, required=True)
    ap.add_argument("--prefix",   default="DB")
    ap.add_argument("--max-reads",type=int, default=10_000)
    args = ap.parse_args(argv)

    libs = discover(args.run_dir, args.prefix)
    if not libs:
        sys.exit("No libraries found.")

    key = sorted(libs)[0]
    r1, r2 = libs[key]["R1"], libs[key]["R2"]
    print(f"⇢ Streaming {key}")

    parser = build_parser(args.config)

    results = []
    with FastqReader(r1, r2) as reader:
        for i, pair in enumerate(reader, 1):
            if i > args.max_reads:
                break
            ctx = parser.parse(*pair)
            if ctx and ctx.get("status") == "ok":
                results.append(ctx)
            if i % 1_000 == 0 or i == 1:
                flat = flatten_dot(parser.get_parse_log())
                if i == 1:
                    print(*flat.keys(), sep="\t")
                else:
                    print(*flat.values(), sep="\t")

    df = pd.DataFrame(results)
    print(f"\n✓ {len(df):,} passing reads of {args.max_reads:,}")

if __name__ == "__main__":   # <-- now import-safe
    main()
