import operator
import re

import editdistance

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


def hamming(
    str1: str, str2: str, use_edit_distance: bool = False
) -> int:
    """Compute hamming distance."""
    # only compares up until end of shortest string
    ne = operator.ne
    if use_edit_distance:
        min_len = min(len(str1), len(str2))
        return int(editdistance.eval(str1[:min_len],str2[:min_len]))

    return int(sum(map(ne, str1, str2)))


def convert_quality_score(qual_ascii: str) -> list[int]:
    qual_numeric = [ord(x)-33 for x in qual_ascii]
    return qual_numeric


def check_tag_structure(vt: str, structure: str) -> bool:
    structureregex = structure
    hits = re.findall(
        '([ACGTN]{3}[AT]{1}[ACGTN]{3}[AT]{1}[ACGTN]{3}[AT]{1}[ACGTN]{3}[AT]{1}[ACGT]{3})',
        vt)
    return len(hits)>0
    # vt1="ACTTGGTACCGTTTTAAAG" # perfect
    # vt2="ACTGGGTACCGTTTTCAAG" # not
    # structure="NNNWNNNWNNNWNNNWNNN"
