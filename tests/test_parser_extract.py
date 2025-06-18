from read_break.parser import ReadParser

# Read with offset match and tag downstream
read1 = "NNNGGGTACCTAG"
read2 = "AAAACCCCGGGG"

pipeline_cfg = {
    "pipeline": [
        {
            "id": "match_s1",
            "read": 1,
            "op": "match",
            "ref": "GGGTAC",
            "hamming_fn": "hamming",
            "max_wobble": 5,
            "max_mismatch": 0,
            "store_pos_as": "s1_start"
        },
        {
            "id": "extract_tag",
            "read": 1,
            "op": "extract",
            "start": 3 + 6,  # match starts at 3, sequence is 6bp long
            "length": 4,
            "store_seq_as": "tag"
        }
    ]
}

parser = ReadParser(pipeline_cfg)
result = parser.parse("ex_read", read1, "~~~~", read2, "~~~~")
assert result["s1_start"] == 3
assert result["tag"] == "CTAG"
print("✅ extract step test passed.")
