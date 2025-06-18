import yaml
from read_break.parser import ReadParser

# -----------------------------------------------------------------------------
# Load pipeline YAML file
# -----------------------------------------------------------------------------
with open("parsers/test_pipeline.yaml", "r") as f:
    pipeline_cfg = yaml.safe_load(f)

# No external globals required for this test
parser = ReadParser(pipeline_cfg)

# -----------------------------------------------------------------------------
# Synthetic read pair (only read1 is used in this pipeline)
# -----------------------------------------------------------------------------
read1 = "NNNGGGTACCTAG"  # GGGTAC match at index 3, then CTAG starts at 9
read2 = "CCCCCCCCCCCC"   # unused

# -----------------------------------------------------------------------------
# Run parser
# -----------------------------------------------------------------------------
result = parser.parse("test_read", read1, "~~~~", read2, "~~~~")

# -----------------------------------------------------------------------------
# Output + Assertions
# -----------------------------------------------------------------------------
print("Result:", result)

assert result["status"] == "ok"
assert result["s1_start"] == 3
assert result["tag"] == "CTAG"
assert result["check_flank_r1"] is True

print("✅ YAML pipeline test passed.")


