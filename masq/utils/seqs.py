

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
