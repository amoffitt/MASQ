import time
from dataclasses import dataclass
from collections import defaultdict
from typing import Any, Optional

import pandas as pd


@dataclass
class EnzymeDescriptor:

    def __init__(
        self, name: str, motif: str, recognition_site: str, cut_offset: int
    ):
        self.name = name
        self.motif = motif
        self.recognition_site = recognition_site
        self.cut_offset = cut_offset


def load_enzyme_descriptors(
    fn: str,
    enzymes: Optional[list[str]] = None
) -> dict[str, EnzymeDescriptor]:
    df = pd.read_csv(
        fn, sep="\t", header=None,
        names=["name", "motif", "recognition_site", "cut_offset"])

    result = {}
    for _, row in df.iterrows():
        if enzymes is not None and row["name"] not in enzymes:
            continue
        result[row["name"]] = EnzymeDescriptor(
            row["name"], row["motif"],
            row["recognition_site"], row["cut_offset"]
        )
    return result


def process_enzyme_cut_sites(
    enzymes: list[EnzymeDescriptor],
    enzymes_cut_sites_folder: str,
    genome_build: str = "hg19"
) -> tuple[dict[str, Any], dict[str, Any]]:
    start = time.time()
    cut_sites_top: dict[str, Any] = defaultdict(dict)
    cut_sites_btm: dict[str, Any] = defaultdict(dict)
    # Loop over enzyme files and store cut site information
    for edesc in enzymes:
        ename = edesc.name
        fname = f"{enzymes_cut_sites_folder}/{ename}.{genome_build}.sort.gz"

        off_top = edesc.cut_offset
        off_btm = len(edesc.motif) - off_top

        df = pd.read_csv(fname, sep="\t", header=None, names=["chrom", "pos"])
        for chrom, index in df.groupby("chrom").groups.items():
            positions = df.loc[index, "pos"]
            pos_top = positions.values + off_top
            pos_btm = positions.values + off_btm

            cut_sites_top[ename][chrom] = list(pos_top)
            cut_sites_btm[ename][chrom] = list(pos_btm)
        elapsed = time.time() - start
        print(f"{ename} processed in {elapsed:.2f} seconds")

    elapsed = time.time() - start
    print(f"enzimes cut sites processed in {elapsed:.2f} seconds")
    return cut_sites_top, cut_sites_btm
