import sys
from typing import Any

import pandas as pd

from loguru import logger

from masq.utils.reference_genome import ReferenceGenome
from masq.utils.seqs import reverse_complement


def initial_snp_dict(
    snpdict: dict[str, dict[str, Any]],
    snpid: str,
    chrom: str,
    pos: int,
    strand: str,
    ref: str,
    alt: str,
    reftrinuc: str,
    alttrinuc: str
) -> None:
    snpdict[snpid]['chrom'] = chrom
    snpdict[snpid]['pos'] = pos
    snpdict[snpid]['strand'] = strand
    snpdict[snpid]['ref'] = ref
    snpdict[snpid]['alt'] = alt
    snpdict[snpid]['ref_trinuc'] = reftrinuc
    snpdict[snpid]['alt_trinuc'] = alttrinuc
    snpdict[snpid]['status'] = 'pass'


def print_snp_dict(snpdict: dict[str, dict[str, Any]], passonly: bool) -> None:
    for snpid, info in snpdict.items():
        if passonly and snpdict[snpid]['status'] != 'pass':
            continue
        print(snpid)
        for key, val in info.items():
            print(f"{key}: {val}")
        print("\n######################\n")
    sys.stdout.flush()


def process_snps(
    variant_file: str, ref_genome: ReferenceGenome
) -> dict[str, Any]:
    snpdict: dict[str, dict] = {}

    df = pd.read_csv(variant_file, sep='\t', header=None)
    if len(df.columns) < 4 or len(df.columns) > 5:
        logger.error(
            "wrong number of columns in the variant file %s", variant_file)
        raise ValueError("wrong number of columns in the variant file")
    if len(df.columns) == 4:
        df[4] = ''

    for _index, rec in df.iterrows():
        chrom = rec[0]
        pos = int(rec[1])
        ref = rec[2]
        alt = rec[3]
        strand = rec[4]

        # check reference base
        seqref = ref_genome.get_sequence(chrom, pos-1, pos)
        if strand == 'bottom':
            seqref = reverse_complement(seqref)

        assert seqref == ref

        if seqref != ref:
            logger.error(
                "reference base (and strand) provided do not match "
                "reference genome: chrom=%s, pos=%d, strand=%s",
                chrom, pos, strand)
            raise ValueError(
                "reference base (and strand) provided do not match "
                "reference genome")

        if strand == 'top':
            ref_trinuc = ref_genome.get_sequence(chrom, pos-2, pos+1)
            alt_trinuc = ref_trinuc[0] + alt + ref_trinuc[2]

            snpid = '_'.join([chrom, str(pos), strand])
            snpdict[snpid] = {}
            initial_snp_dict(
                snpdict, snpid,
                chrom, pos, strand,
                ref, alt,
                ref_trinuc, alt_trinuc)

        elif strand == 'bottom':
            ref_trinuc = reverse_complement(
                ref_genome.get_sequence(chrom, pos-2, pos+1))
            alt_trinuc = ref_trinuc[0] + alt + ref_trinuc[2]

            snpid = '_'.join([chrom, str(pos), strand])
            snpdict[snpid] = {}
            initial_snp_dict(
                snpdict, snpid,
                chrom, pos, strand,
                ref, alt,
                ref_trinuc, alt_trinuc)
        else:
            ref_trinuc_fwd = ref_genome.get_sequence(chrom, pos-2, pos+1)
            alt_trinuc_fwd = ref_trinuc_fwd[0] + alt + ref_trinuc_fwd[2]

            ref_trinuc_rev = reverse_complement(ref_trinuc_fwd)
            alt_trinuc_rev = ref_trinuc_rev[0] + \
                reverse_complement(alt) + \
                ref_trinuc_rev[2]

            strand = 'top'
            snpid = '_'.join([chrom, str(pos), strand])
            snpdict[snpid] = {}
            initial_snp_dict(
                snpdict, snpid,
                chrom, pos, strand,
                ref, alt,
                ref_trinuc_fwd, alt_trinuc_fwd)

            strand = 'bottom'
            ref_rev = reverse_complement(ref)
            alt_rev = reverse_complement(alt)
            snpid = '_'.join([chrom, str(pos), strand])
            snpdict[snpid] = {}
            initial_snp_dict(
                snpdict, snpid,
                chrom, pos, strand,
                ref_rev, alt_rev,
                ref_trinuc_rev, alt_trinuc_rev)
    return snpdict
