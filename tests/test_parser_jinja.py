from read_break.parser import ReadParser

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
            "start": "{{ s1_start + 6 }}",  # match at 3, tag at 9
            "length": 4,
            "store_seq_as": "tag"
        }
    ]
}

# Inject constants (could be loaded from YAML later)
cfg = {"LT_LEN": 15}

parser = ReadParser(pipeline_cfg, globals_cfg={"cfg": cfg})
result = parser.parse("jinja_read", read1, "~~~~", read2, "~~~~")
assert result["s1_start"] == 3
assert result["tag"] == "CTAG"
print("✅ jinja-templated extract test passed.")