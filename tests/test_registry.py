from read_break.registry import HAMMING_FUNCS

# --- vanilla Hamming ---
assert HAMMING_FUNCS["hamming"]("ACGT", "ACGT") == 0
assert HAMMING_FUNCS["hamming"]("ACGT", "CTAG") == 4

# --- T→C–tolerant Hamming ---
s1, s2 = "TTTTT", "CCCCG"          # 1 true mismatch (pos 4: A vs A) + 3 T→C conversions
assert HAMMING_FUNCS["hammingTC"](s1, s2) == 1

# --- A→G–tolerant Hamming ---
s1, s2 = "GGGGG", "GGGGT"          # 1 true mismatch (last pos) + 4 A→G conversions
assert HAMMING_FUNCS["hammingAG"](s1, s2) == 1

print("✅ registry stand-alone tests passed.")
