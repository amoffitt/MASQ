from __future__ import annotations

import os
import logging
from typing import Optional, TextIO
from types import TracebackType

logger = logging.getLogger(__name__)


class ReferenceGenome:

    def __init__(
        self, filename: str, index_filename: Optional[str] = None
    ) -> None:
        if not os.path.exists(filename):
            raise IOError(
                f"reference genome file is missing: {filename}")
        if index_filename is None:
            index_filename = f"{filename}.fai"
        if not os.path.exists(index_filename):
            raise IOError(
                f"reference genome index file is missing: {index_filename}")
        self.filename = filename
        self.index_filename = index_filename
        self.chromosomes: dict[str, int] = {}
        self._genome_index: dict[str, dict[str, int]] = {}
        self._genome_file: Optional[TextIO] = None

    def _load_index(self) -> None:
        self._genome_index = {}
        with open(self.index_filename, "r") as infile:
            while True:
                line = infile.readline()
                if not line:
                    break
                parts = [str(p.strip()) for p in line.split()]
                self._genome_index[parts[0]] = {
                    "length": int(parts[1]),
                    "start": int(parts[2]),
                    "seq_line_length": int(parts[3]),
                    "line_length": int(parts[4]),
                }
                self.chromosomes[parts[0]] = int(parts[1])

    def is_open(self) -> bool:
        return self._genome_file is not None

    def open(self) -> ReferenceGenome:
        # pylint: disable=consider-using-with
        self._genome_file = open(self.filename, "r")
        self._load_index()
        return self

    def close(self) -> None:
        if self.is_open():
            assert self._genome_file is not None
            self._genome_file.close()
            self._genome_file = None

        self._genome_index = {}

    def get_chrom_length(self, chrom: str) -> int:
        if not self.is_open():
            raise IOError(f"reference genome is not open: {self.filename}")
        if chrom not in self._genome_index:
            raise ValueError(f"contig {chrom} not found in reference genome")
        return self._genome_index[chrom]["length"]

    def get_chromosomes(self) -> list[str]:
        if not self.is_open():
            raise IOError(f"reference genome is not open: {self.filename}")
        return list(self._genome_index.keys())

    def get_sequence(
        self, chrom: str, start: int, stop: int
    ) -> str:
        """Returns a sequence from the reference genome.

        Returns a subsequence from the reference genome at chromosome `chrom`
        in the interval [start, stop). The coordinates are 0-based.
        Args:
            chrom: chromosome name
            start: 0-based start position
            stop: 0-based stop position
        """
        if not self.is_open():
            raise IOError(f"reference genome is not open: {self.filename}")
        assert self._genome_file is not None

        if chrom not in self._genome_index:
            raise ValueError(f"unknown chromosome: {chrom}")

        seq_line_length = self._genome_index[chrom]["seq_line_length"]
        seq_start = self._genome_index[chrom]["start"]
        start_index: int = seq_start \
            + start \
            + start // seq_line_length
        self._genome_file.seek(start_index)

        ll = stop - start
        x = 1 + ll // seq_line_length

        w = self._genome_file.read(ll + x)
        w = w.replace("\n", "")[:ll]

        return w.upper()

    def __enter__(
        self
    ) -> ReferenceGenome:
        return self.open()

    def __exit__(
            self,
            exc_type: type[BaseException] | None,
            exc_value: Optional[BaseException],
            exc_tb: TracebackType | None) -> None:
        self.close()
