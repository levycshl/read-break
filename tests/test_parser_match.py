from read_break.parser import ReadParser

# Synthetic read where target appears after 3 bases of wobble
read1 = "NNNGGGTACCTAG"     # match at offset=3
read2 = "AAAACCCCGGGG"

# Define minimal pipeline to just match something in read1
pipeline_cfg = {
    "pipeline": [
        {
            "id": "test_match",
            "read": 1,
            "op": "match",
            "ref": "GGGTAC",
            "hamming_fn": "hamming",
            "max_wobble": 3,
            "max_mismatch": 0,
            "store_pos_as": "s1_start"
        }
    ]
}

parser = ReadParser(pipeline_cfg)
result = parser.parse("test_read", read1, "~~~~", read2, "~~~~")
print(result)
assert result["s1_start"] == 3
assert result["read_id"] == "test_read"

print("✅ match-only parser test passed.")
