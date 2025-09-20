from operator import ne
from collections.abc import Mapping


"""
Low-level sequence logic for read-break.

Includes base-to-integer conversion, hamming distance functions (standard and asymmetric),
and fuzzy matching with wobble offsets.

These functions are stateless and pure.
"""

from typing import List, Callable


BASES    = ['A', 'C', 'G', 'T', 'N']
BASE2INT = {base: idx for idx, base in enumerate(BASES)}
# Dictionary used to invert read sequences
INDICT = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', ' ':' ',
          'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'}

def reverse_complement(seq: str) -> str:
    """
    Returns the reverse complement of the given sequence.
    
    Args:
        seq (str): A nucleotide sequence
        
    Returns:
        str: The reverse complement of the input sequence
    """
    return "".join([INDICT[base] for base in seq[::-1]])

def seq_to_int(seq: str) -> List[int]:
    """
    Converts a nucleotide sequence into a list of integer indices.

    Args:
        seq (str): A string of nucleotide bases (A, C, G, T, N)

    Returns:
        List[int]: List of integer-encoded bases, using BASE2INT mapping.
    """
    return [BASE2INT[base] for base in seq]





def hamming(x: str, y: str) -> int:
    """
    Computes the Hamming distance between two equal-length strings.

    Args:
        x (str): First string.
        y (str): Second string.

    Returns:
        int: Number of mismatched positions.
    """
    return sum(map(ne, x, y))


def hamming35(
    from_string: str,
    to_string: str,
    from_letter: str,
    to_letter: str
) -> int:
    """
    Variant of Hamming distance that ignores a specific base conversion.
    Designed for detecting asymmetric enzyme-induced base changes (e.g. A→G or T→C).

    Args:
        from_string (str): Original sequence.
        to_string (str): Converted sequence.
        from_letter (str): Base to ignore in conversions.
        to_letter (str): Target base that from_letter is allowed to change into.

    Returns:
        int: Modified Hamming distance (ignores specified conversions).
    """
    distance = 0
    for a, b in zip(from_string, to_string):
        if a != b and not (a == from_letter and b == to_letter):
            distance += 1
    return distance


def wobble_match(
    test: str,
    target: str,
    max_wobble: int,
    max_hamming: int,
    base_offset: int = 0,
    hamming_func: Callable[[str, str], int] = hamming
) -> int:
    """
    Attempts to align a target sequence within a wobble window of the test sequence.

    Args:
        test (str): Read sequence to search within.
        target (str): Reference sequence to match.
        max_wobble (int): Number of bases to shift forward for fuzzy alignment.
        max_hamming (int): Maximum allowable mismatches for a valid match.
        base_offset (int): Starting position in test to begin offset from.
        hamming_func (Callable): Hamming function (possibly asymmetric).

    Returns:
        int: Offset relative to base_offset if match is found, else -1.
    """
    for offset in range(base_offset, base_offset + max_wobble + 1):
        test_substring = test[offset:(offset + len(target))]
        hval = hamming_func(target, test_substring)
        if len(test_substring) != len(target):
            break
        if hamming_func(target, test_substring) <= max_hamming:
            return offset - base_offset
    return -1


def flatten_dot(d: Mapping, prefix: str = "", sep: str = ".") -> dict[str, object]:
    """Return a flat dict: {'a.b.c': value, ...}"""
    flat = {}
    for k, v in d.items():
        path = f"{prefix}{sep}{k}" if prefix else k
        if isinstance(v, Mapping):
            flat.update(flatten_dot(v, path, sep=sep))
        else:
            flat[path] = v
    return flat