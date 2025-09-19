from read_break.parser import ReadParser

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

# Read 1 is not used in this test. Read 2 contains a target for match + flank check.
read1 = "NNNGGGTACCTAG"
read2 = "AAAGGGTTTTCCGGG"  # 'TTTTCC' is expected to be checked with hamming
read3 = "AAAGGGTTTTTTGGG"  # 'TTTTTT' has too many mismatches vs. 'TTTTCC' but not with hammingTC

# Define a 2-step pipeline:
# 1. Match "AAAGGG" in read2 (at position 0)
# 2. Then hamming_test against "TTTTCC" at offset +6 (i.e., positions 6–11)
pipeline_cfg = {
    "pipeline": [
        {
            "id": "match_s2",
            "read": 2,
            "op": "match",
            "ref": "AAAGGG",
            "hamming_fn": "hamming",
            "max_wobble": 3,
            "max_mismatch": 0,
            "store_pos_as": "s2_start"
        },
        {
            "id": "check_flank",
            "read": 2,
            "op": "hamming_test",
            "ref": "TTTTCC",
            "start": "{{ s2_start + 6 }}",  # This evaluates to 6
            "length": 6,
            "hamming_fn": "hamming",
            "max_mismatch": 1,
            "must_pass": False
        },
        {
            "id": "check_flank_TC",
            "read": 2,
            "op": "match",
            "ref": "TTTTCC",
            "hamming_fn": "hammingTC",
            "start": "{{ s2_start + 6 }}",  # This evaluates to 6
            "length": 6,
            "max_mismatch": 1,
            "must_pass": True
        }
    ]
}

parser = ReadParser(pipeline_cfg)

# -----------------------------------------------------------------------------
# ✅ Test Case 1: Successful flank match (within mismatch tolerance)
# -----------------------------------------------------------------------------

result = parser.parse("hamming_pass", read1, "~~~~", read2, "~~~~")

assert result is not None, "Expected result dictionary, got None"
assert result["status"] == "ok"
assert result["read_id"] == "hamming_pass"
assert result["check_flank"] is True  # Hamming test result stored under its step ID
assert result["s2_start"] == 0

print("✅ PASS: Hamming test succeeded when mismatches were within tolerance.")

# -----------------------------------------------------------------------------
# ❌ Test Case 2: Failing flank match (too many mismatches)
# -----------------------------------------------------------------------------

bad_read2 = "AAAGGGAAAAAACCCC"  # mismatches vs. TTTTCC at positions 6–11

result = parser.parse("hamming_fail", read1, "~~~~", bad_read2, "~~~~")

assert result is not None, "Expected structured fail dictionary"
assert result["status"] == "fail"
assert result["read_id"] == "hamming_fail"
assert result["failed_step"] == "check_flank"
print(result["message"])  # Should contain details about the failure
assert "Hamming" in result["message"]

print("✅ PASS: Hamming test correctly failed with structured fail response.")
