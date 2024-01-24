# pylint: disable=W0621,C0114,C0116,W0212,W0613
from typing import Callable
from contextlib import closing

import pytest

from masq.utils.reference_genome import ReferenceGenome


@pytest.fixture(scope='session')
def acgt_genome(fixtures_filename: Callable[[str], str]) -> ReferenceGenome:
    filename = fixtures_filename('reference_genome/allChr.fa')
    index_filename = fixtures_filename('reference_genome/allChr.fa.fai')
    return ReferenceGenome(filename, index_filename)


def test_reference_genome_simple(acgt_genome: ReferenceGenome) -> None:
    assert acgt_genome is not None
    assert not acgt_genome.is_open()

    with closing(acgt_genome.open()) as genome:
        assert genome.is_open()
        assert genome.get_chromosomes() == ['chr1', 'chr2', 'chr3']


@pytest.mark.parametrize('chrom, start, stop, expected', [
    ('chr1', 1, 4, 'ACGT'),
    ('chr1', 101, 104, 'ACGT'),
    ('chr2', 101, 104, 'ACGT'),
    ('chr2', 201, 204, 'ACGT'),
    ('chr3', 197, 204, 'ACGTACGT'),
])
def test_reference_genome_get_sequence(
    acgt_genome: ReferenceGenome,
    chrom: str,
    start: int,
    stop: int,
    expected: str
) -> None:
    with closing(acgt_genome.open()) as genome:
        seq = genome.get_sequence(chrom, start, stop)
        assert seq == expected


@pytest.mark.parametrize('chrom, start, stop, expected', [
    ('chr1', 1, 4, 'ACGT'),
    ('chr1', 101, 104, 'ACGT'),
    ('chr2', 101, 104, 'ACGT'),
    ('chr2', 201, 204, 'ACGT'),
    ('chr3', 197, 204, 'ACGTACGT'),
])
def test_reference_genome_contextmanager_get_sequence(
    acgt_genome: ReferenceGenome,
    chrom: str,
    start: int,
    stop: int,
    expected: str
) -> None:
    with acgt_genome as genome:
        seq = genome.get_sequence(chrom, start, stop)
        assert seq == expected
