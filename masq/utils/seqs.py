import re

_COMPL_NUC = {
    'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
    'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'
}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a sequence."""
    return "".join([_COMPL_NUC[base] for base in seq[::-1]])


def complement(seq: str) -> str:
    """Return the complement of a sequence."""
    return "".join([_COMPL_NUC[base] for base in seq])


def check_sequence_for_cut_site(sequence: str, pattern: str) -> bool:
    results = re.findall(pattern, sequence)
    if len(results) > 0:
        return True  # pattern is in seuqence

    return False


BASES = ["A", "C", "G", "T", "N"]
BASE2INT = {base: index for index, base in enumerate(BASES)}


def base2int(base: str) -> int:
    return BASE2INT[base]


def int2base(intbase: int) -> str:
    return BASES[intbase]
