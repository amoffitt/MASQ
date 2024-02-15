"""Converts primer table output from primer design functions to input table
needed for MASQ analysis pipeline
"""

import sys
import argparse
from typing import Optional

import pandas as pd

from masq.utils.regions import Region
from masq.utils.seqs import reverse_complement
from masq.utils.reference_genome import ReferenceGenome


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Converts primer table from primer design "
        "functions to input table needed for MASQ analysis pipeline"
    )
    parser.add_argument(
        "primer_table",
        help="input primer table",
    )
    parser.add_argument(
        "output_table",
        help="output table",
    )
    parser.add_argument(
        "ref_genome",
        help="reference genome",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]
    args = parse_args(argv)

    primer_df = pd.read_csv(args.primer_table, sep="\t")
    result = []

    with ReferenceGenome(args.ref_genome) as ref_genome:
        for index, rec in primer_df.iterrows():
            chrom = rec["chrom"]
            pos = rec["pos"]
            primer1 = rec["downstream_primerseq"]
            primer2 = rec["cutadj_primerseq"]
            strand = rec["strand"]
            refbase = rec["ref_trinuc"][1]
            altbase = rec["alt_trinuc"][1]
            ref_alt = f"{refbase}_{altbase}"

            # get target sequence coordinates based on primer Coordinates
            # pull sequence from dictionary
            # reverse complement depending on strand
            # target position is relative to target seq
            primer_region1 = Region.from_string(
                rec["cutadj_primer_coordinates"])
            primer_region2 = Region.from_string(
                rec["downstream_primer_coordinates"])
            assert primer_region1.chrom == primer_region2.chrom == chrom

            primercoords = [
                primer_region1.start, primer_region1.stop,
                primer_region2.start, primer_region2.stop
            ]
            primercoords.sort()

            target_start = primercoords[1]
            target_end = primercoords[2] - 1

            target_seq_ref = ref_genome.get_sequence(
                chrom, target_start, target_end)

            if strand == "bottom":
                target_seq_stranded = reverse_complement(target_seq_ref)
                relative_pos = target_end - pos
            elif strand == "top":
                target_seq_stranded = target_seq_ref
                relative_pos = pos - target_start - 1

            result.append({
                "loc": index,
                "chr": chrom,
                "posi": pos,
                "specific-primer-1": primer1,
                "specific-primer-2": primer2,
                "trimmed-target-seq": target_seq_stranded,
                "target_locs": relative_pos,
                "ref-alt_allele": ref_alt,
            })

    result_df = pd.DataFrame.from_records(result)
    result_df.to_csv(args.output_table, sep="\t", index=False)


if __name__ == "__main__":
    main()
