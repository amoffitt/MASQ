# pylint: disable=W0621,C0114,C0116,W0212,W0613

from typing import Callable

import pytest

from masq.utils.paired_reads import PairedReads


def test_paired_reads_simple(fixtures_filename: Callable[[str], str]) -> None:
    paired_reads = PairedReads(
        fixtures_filename('paired_reads/r1.fastq'),
        fixtures_filename('paired_reads/r2.fastq'))
    assert paired_reads is not None

    count = 0
    for entry1, entry2 in paired_reads.fetch():
        count += 1
        assert entry1.name[: -1] == entry2.name[: -1]

    assert count == 10


@pytest.mark.parametrize('r1_filename, r2_filename', [
    ('r1.fastq', 'r2_missing_last_1_line.fastq'),
    ('r1_missing_last_1_line.fastq', 'r2.fastq'),
    ('r1_missing_last_1_line.fastq', 'r2_missing_last_1_line.fastq'),
    ('r1.fastq', 'r2_missing_last_2_line.fastq'),
    ('r1_missing_last_2_line.fastq', 'r2.fastq'),
    ('r1_missing_last_2_line.fastq', 'r2_missing_last_2_line.fastq'),
    ('r1.fastq', 'r2_missing_last_3_line.fastq'),
    ('r1_missing_last_3_line.fastq', 'r2.fastq'),
    ('r1_missing_last_3_line.fastq', 'r2_missing_last_3_line.fastq'),
    ('r1.fastq', 'r2_missing_last_4_line.fastq'),
    ('r1_missing_last_4_line.fastq', 'r2.fastq'),
    ('r1_missing_last_4_line.fastq', 'r2_missing_last_4_line.fastq'),
])
def test_paired_reads_bad_r1_or_r2_reads(
    fixtures_filename: Callable[[str], str],
    r1_filename: str,
    r2_filename: str
) -> None:
    paired_reads = PairedReads(
        fixtures_filename(f'paired_reads/{r1_filename}'),
        fixtures_filename(f'paired_reads/{r2_filename}'))
    assert paired_reads is not None

    count = 0
    for entry1, entry2 in paired_reads.fetch():
        count += 1
        assert entry1.name[: -1] == entry2.name[: -1]

    assert count == 9


def test_paired_reads_context_manager(
    fixtures_filename: Callable[[str], str]
) -> None:
    with PairedReads(
            fixtures_filename('paired_reads/r1.fastq'),
            fixtures_filename('paired_reads/r2.fastq')) as paired_reads:
        assert paired_reads is not None

        count = 0
        for entry1, entry2 in paired_reads:
            count += 1
            assert entry1.name[: -1] == entry2.name[: -1]

    assert count == 10
