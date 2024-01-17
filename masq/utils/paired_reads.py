from __future__ import annotations

import logging
import os
from typing import Optional, Generator
from types import TracebackType

import pysam

logger = logging.getLogger(__name__)


class PairedReads:
    '''Class for walking through a paired-end data of fastq.gz files.

    The reads are organized in pairs in the file.
    Once initialized with two filenames (for read1 and read2 fastq.gz files)
    it returns an iterator which returns one read pair at a time.
    '''

    def __init__(
        self, read1filename: str, read2filename: str
    ) -> None:
        if not os.path.exists(read1filename) or \
                not os.path.exists(read2filename):
            raise IOError(
                f"File not found: {read1filename} or {read2filename}")
        self.r1filename = read1filename
        self.r2filename = read2filename
        self.read1file: Optional[pysam.FastxFile] = None
        self.read2file: Optional[pysam.FastxFile] = None

    def open(self) -> PairedReads:
        self.read1file = pysam.FastxFile(self.r1filename)
        self.read2file = pysam.FastxFile(self.r2filename)
        return self

    def is_open(self) -> bool:
        return self.read1file is not None and self.read2file is not None

    def fetch(
        self
    ) -> Generator[tuple[pysam.FastqRecord, pysam.FastqRecord], None, None]:
        if self.read1file is None or self.read2file is None:
            raise IOError("Files not opened")

        try:
            for entry1, entry2 in zip(self.read1file, self.read2file):
                if not entry1.name or not entry2.name:
                    raise ValueError(
                        "missing name line(s) in fastq")
                if entry1.quality is None or entry2.quality is None:
                    raise ValueError(
                        "quality is None; missing quality line(s) in fastq")
                if not entry1.sequence or not entry2.sequence:
                    raise ValueError(
                        "missing sequence line(s) in fastq")
                yield entry1, entry2
        except ValueError:
            logger.error("value error", exc_info=True)

    def close(self) -> None:
        if not self.is_open():
            return
        assert self.read1file is not None
        assert self.read2file is not None
        self.read1file.close()
        self.read2file.close()

    def __enter__(
        self
    ) -> Generator[tuple[pysam.FastqRecord, pysam.FastqRecord], None, None]:
        self.open()
        return self.fetch()

    def __exit__(
            self,
            exc_type: type[BaseException] | None,
            exc_value: Optional[BaseException],
            exc_tb: TracebackType | None) -> None:
        self.close()
