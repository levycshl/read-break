# read_break/registry.py
"""
Registry of helper functions that the YAML pipeline can reference by name.
Feel free to extend this dict later with new mismatch-tolerant metrics.
"""

from read_break.logic import hamming, hamming35

HAMMING_FUNCS = {
    "hamming": hamming,
    "hammingTC": lambda x, y: hamming35(x, y, "T", "C"),
    "hammingAG": lambda x, y: hamming35(x, y, "A", "G"),
}
