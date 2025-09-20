
import os
import sys
from pathlib import Path
import re
from collections import defaultdict
import yaml
from read_break.parser import ReadParser
from read_break.logic import flatten_dot
from read_break.io import FastqReader
import pandas as pd


root_dir = Path(os.getcwd()).parent

data_dir = Path( r"Z:\mnt\wigstore3\data\checkpointcharlie2\inessa\250707_M06142_0497_000000000-M2RFJ\Alignment_1\20250708_193955\Fastq" )

# directory that contains the MiSeq run
RUN_DIR = Path(data_dir)          
# regex: capture key = <field1>_<field2>, read = 1 or 2
PAT = re.compile(r'^(?P<key>[^_]+_[^_]+)_L\d{3}_R(?P<read>[12])_\d{3}\.fastq\.gz$')
libs: dict[str, dict[str, Path]] = defaultdict(lambda: {"R1": None, "R2": None})
for f in RUN_DIR.glob("*.fastq.gz"):
    m = PAT.match(f.name)
    if not m:
        continue                       # skip unexpected names
    key, read = m.group("key", "read")
    libs[key][f"R{read}"] = f.resolve()

# sanity-check: complain if any library is missing a mate
missing = {k: v for k, v in libs.items() if None in v.values()}
if missing:
    raise ValueError(f"Mate not found for: {missing}")

## restrict to keys that start with 'DB'
libs = {k: v for k, v in libs.items() if k.startswith("DB")}
key_list = sorted(libs.keys())
print(f"\nKeys: {key_list}")
key = key_list[0]
read1_filename, read2_filename = libs[key]["R1"], libs[key]["R2"]
print(f"\nRead 1: {read1_filename}\nRead 2: {read2_filename}")


## get the parser configuration

parsers_dir = root_dir / "parsers"
parse_config = parsers_dir / "a_deaminate_2025_07_10.yaml"

parser_cfg = yaml.safe_load( parse_config.read_text() )

parser = ReadParser(parser_cfg,
                    parser_cfg['params'])

print(parser)

## collect results
results = []

MAX_READS = 10000  # limit for testing, set to None for no limit
with FastqReader(read1_filename, read2_filename) as reader:
    for ind, read_pair in enumerate( reader ):
        if ind >= MAX_READS:
            break
        if ind % 1000 == 0:        
            parse_log = flatten_dot( parser.get_parse_log() )        
            if ind == 0:
                print(*list( parse_log.keys() ), sep="\t")
            else:
                print(*list( parse_log.values() ), sep="\t")                            
        ctx = parser.parse(*read_pair)
        if not ctx or ctx.get("status") != "ok":
            continue
        else:
            results.append(ctx)

parser.get_parse_log()

df = pd.DataFrame.from_dict(results)