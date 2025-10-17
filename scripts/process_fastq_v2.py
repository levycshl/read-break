#!/usr/bin/env python3
"""
Process FASTQ files using ReadParser and FastqReader/Writer.
"""

import os
import re
from pathlib import Path
from collections import defaultdict

import yaml
import pandas as pd
import numpy as np

from read_break.parser import ReadParser, read_clip_and_write
from read_break.logic import flatten_dot
from read_break.io import FastqReader, FastqWriter


def collect_libraries(run_dir: Path) -> dict[str, dict[str, Path]]:
    """Collect paired FASTQ files into a dict keyed by library."""
    pattern = re.compile(
        r'^(?P<key>[^_]+_[^_]+)_L\d{3}_R(?P<read>[12])_\d{3}\.fastq\.gz$'
    )
    libs: dict[str, dict[str, Path]] = defaultdict(lambda: {"R1": None, "R2": None})

    for f in run_dir.glob("*.fastq.gz"):
        match = pattern.match(f.name)
        if not match:
            continue
        key, read = match.group("key", "read")
        libs[key][f"R{read}"] = f.resolve()

    # sanity-check: complain if any library is missing a mate
    missing = {k: v for k, v in libs.items() if None in v.values()}
    if missing:
        raise ValueError(f"Mate not found for: {missing}")

    return libs


def main():
    root_dir = Path(os.getcwd()).parent

    # data_dir = Path(
    #     r"Z:/mnt/wigstore3/data/checkpointcharlie2/Junyi/"
    #     "250825_M06142_0500_000000000-M5NRD/Alignment_1/"
    #     "20250826_160009/Fastq"
    # )

    data_dir = Path(r"Z:/mnt/wigstore3/data/checkpointcharlie/Debmalya/250918_M00164_0701_000000000-M6HLD/Alignment_1/20250919_152010/Fastq")
    output_dir = Path(r"Z:/data/safe/levy/read_break/2025_09_19_alu_capture_deaminate_second")

    run_dir = Path(data_dir)

    libs = collect_libraries(run_dir)

    # restrict to keys that start with 'D'
    libs = {k: v for k, v in libs.items()}
    key_list = sorted(libs.keys())
    
    print(f"\nKeys: {key_list}")

    # load parser configuration
    parsers_dir = root_dir / "parsers"
    parse_config = parsers_dir / "alu_bag_2025_09_03.yaml"
    parser_cfg = yaml.safe_load(parse_config.read_text())
    print(parser_cfg)

    for key in key_list:
        print(f"\nProcessing library: {key}")
        read1_filename, read2_filename = libs[key]["R1"], libs[key]["R2"]
        print(f"Read 1: {read1_filename}\nRead 2: {read2_filename}")

        reader = FastqReader(read1_filename, read2_filename)
        parser = ReadParser(parser_cfg, parser_cfg["params"])
        writer = FastqWriter(output_dir, key)

        read_clip_and_write(
            reader,
            parser,
            writer,
            verbose=True,
            verbose_interval=10000,
        )


if __name__ == "__main__":
    main()
