import sys
import logging
from typing import Any, Union
from collections import defaultdict

import pandas as pd
import numpy as np

import pysam

from masq.utils.reference_genome import ReferenceGenome
from masq.utils.regions import Region
from masq.utils.seqs import reverse_complement, check_sequence_for_cut_site, \
    base2int, int2base
from masq.primer_design.enzymes import EnzymeDescriptor


logger = logging.getLogger(__name__)


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
    """Load SNPs into dictionary with initial information.

    Chrom, pos, ref, alt, strand (optional: top or bottom)
    If strand is included - ref and alt are expected to be flipped for bottom
    """
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


def filter_high_error_trinucleotides(
    trinucleotide_file: str,
    snpdict: dict[str, dict[str, Any]]
) -> dict[str, dict[str, Any]]:
    """Filter trinucleotides with high error rate."""
    high_error_trinucs = []
    with open(trinucleotide_file, 'rt') as infile:
        for line in infile:
            parts = [p.strip() for p in line.strip().split()]
            assert len(parts) == 2
            reftrinuc = parts[0]
            alttrinuc = parts[1]
            high_error_trinucs.append((reftrinuc, alttrinuc))

    for _snpid, snpdesc in snpdict.items():
        if snpdesc["status"] != "pass":
            continue
        trinucs = (snpdesc['ref_trinuc'], snpdesc['alt_trinuc'])
        if trinucs in high_error_trinucs:
            snpdesc["status"] = "drop"
            snpdesc["drop_reason"] = "high_error_trinuc"
    return snpdict


def check_snps_for_enzyme_cut_sites(
    snpdict: dict[str, dict[str, Any]],
    enzyme_descriptors: dict[str, EnzymeDescriptor],
    ref_genome: ReferenceGenome
) -> dict[str, list[str]]:
    """Check mutation changes against enzyme cut sites.

    Return a list of bad enzymes per mutation.
    Recognition sites only need checking in top strand - if they are there,
    they are in bottom too.
    """
    logger.info("Checking mutations for cut sites")
    bad_enzyme_choices = defaultdict(list)

    for snpid, snpdesc in snpdict.items():
        if snpdesc['status'] != 'pass':
            continue

        chrom = snpdesc['chrom']
        pos = snpdesc['pos']
        alt = snpdesc['alt']
        strand = snpdesc['strand']
        ref_context = ref_genome.get_sequence(chrom, pos-7, pos+6)
        if strand == 'top':
            alt_context = ref_genome.get_sequence(chrom, pos-7, pos-1) + \
                alt + \
                ref_genome.get_sequence(chrom, pos, pos+6)

        else:
            alt_context = ref_genome.get_sequence(chrom, pos-7, pos-1) + \
                reverse_complement(alt) + \
                ref_genome.get_sequence(chrom, pos, pos+6)

        for ename, edesc in enzyme_descriptors.items():
            ref_hit = check_sequence_for_cut_site(ref_context, edesc.motif)
            alt_hit = check_sequence_for_cut_site(alt_context, edesc.motif)
            if (alt_hit and not ref_hit):
                logger.info(
                    "%s: mutation introduces cut site for %s",
                    snpid, ename)
                bad_enzyme_choices[snpid].append(ename)
            if (ref_hit and not alt_hit):
                logger.info(
                    "%s: mutation removes cut site for %s",
                    snpid, ename)
                bad_enzyme_choices[snpid].append(ename)

    logger.info(
        "Bad enzyme choices due to cut sites introduced by mutations: %s",
        bad_enzyme_choices)
    return bad_enzyme_choices


def strip_chrom_prefix(chrom: str) -> str:
    """Strip 'chr' prefix from chromosome name."""
    if chrom.startswith('chr'):
        return chrom[3:]
    return chrom


def snps_or_indels_in_region(
    bamfiles: Union[str, list[str]],
    regionstring: str,
    ref_genome: ReferenceGenome,
    basequal_cutoff: int = 28,
    vaf_cutoff: float = 0.05,
    indelaf_cutoff: float = 0.05,
    var_count_cutoff: int = 2,
    indel_count_cutoff: int = 2
) -> Any:
    # Coordinates should be 0-based: (VCF/IGV position)-1 = Python Position
    # Regions are defined as [start,stop), including start but not stop
    # positions
    # bamfiles here can be one string for one bam, or a list of strings for
    # multiple bams
    if isinstance(bamfiles, str):
        bamlist = [bamfiles]  # convert to iterable list
    else:
        bamlist = bamfiles

    # # Extract coordinates
    print(regionstring)
    # [chrom,start,end] = split_region_string(regionstring)
    region = Region.from_string(regionstring)
    if region.chrom not in ref_genome.chromosomes:
        region.chrom = strip_chrom_prefix(region.chrom)
    assert region.chrom in ref_genome.chromosomes

    # Get reference sequence
    # region_ref_seq = seq_dic[chrom][(start-1):end]
    region_ref_seq = ref_genome.get_sequence(
        region.chrom, region.start-1, region.stop)

    print(region_ref_seq)

    # Initialize
    snp_positions: list[int] = []
    snp_alt_bases: list[str] = []

    for bam in bamlist:
        # Load BAM file
        snps_or_indels_in_bamfile_region(
            bam, region, ref_genome,
            snp_positions, snp_alt_bases,
            basequal_cutoff,
            vaf_cutoff,
            indelaf_cutoff,
            var_count_cutoff,
            indel_count_cutoff)

    # Return 2 lists
    # Positions of alterations in region
    # Altered base/indel at each position
    seq_positions = [int(x) - region.start for x in snp_positions]

    return [
        snp_positions,
        seq_positions,
        snp_alt_bases,
        region_ref_seq
    ]


def snps_or_indels_in_bamfile_region(
    bam: str, region: Region, ref_genome: ReferenceGenome,
    snp_positions: list[int], snp_alt_bases: list[str],
    basequal_cutoff: int = 28,
    vaf_cutoff: float = 0.05,
    indelaf_cutoff: float = 0.05,
    var_count_cutoff: int = 2,
    indel_count_cutoff: int = 2
) -> None:
    # Load BAM file
    with pysam.AlignmentFile(bam, "r") as bamfile:
        print(bamfile)
        # Pileup columns
        for pileupcolumn in bamfile.pileup(
                region.chrom,
                region.start - 1,
                region.stop + 1,
                truncate=True,
                stepper='nofilter',
                max_depth=10000000):
            # stepper=all filters out pcr duplicates,
            # set stepper=nofilter to not filter pcr duplicates

            # What position are we at?
            currpos = pileupcolumn.reference_pos  # type: ignore
            currbase = ref_genome.get_sequence(
                region.chrom, currpos, currpos+1)
            refbaseint = base2int(currbase)

            basecounts = np.array([0, 0, 0, 0, 0])
            totalcount = 0
            indelcount = 0
            has_something = False

            for pileupread in pileupcolumn.pileups:  # type: ignore
                if pileupread.indel != 0:
                    indelcount += 1
                    totalcount += 1
                    continue
                if ((pileupread.alignment.qual is None) or
                        (pileupread.query_position is None)):
                    continue
                basequal = ord(pileupread.alignment.qual[
                    pileupread.query_position]) - 33
                if basequal < basequal_cutoff:
                    # skip reads with quality less than filter value
                    continue
                mapqual = pileupread.alignment.mapping_quality
                if mapqual < 10:
                    # skip reads with quality less than filter value
                    continue
                # Extract base from read at this position
                totalcount += 1
                base = pileupread.alignment.seq[pileupread.query_position]
                baseint = base2int(base)
                basecounts[baseint] += 1

            # Drop N's
            basecounts = basecounts[1:]
            # Check if this pileup column has snp
            if np.sum(basecounts) == 0:
                base_ratios = basecounts / 1.0
            else:
                base_ratios = basecounts / np.sum(basecounts)
            # Drop reference base from ratios
            ref_zero_base_ratios = base_ratios
            ref_zero_base_ratios[refbaseint-1] = 0
            # how many non-ref bases are here
            non_ref_base_count = np.sum(ref_zero_base_ratios > vaf_cutoff)
            if non_ref_base_count > 0:
                ref_zero_base_counts = basecounts
                ref_zero_base_counts[refbaseint-1] = 0
                if max(ref_zero_base_counts) > var_count_cutoff:
                    altbase = int2base(int(np.argmax(ref_zero_base_ratios)))
                    has_something = True
            if indelcount > indel_count_cutoff:
                if float(indelcount)/totalcount > indelaf_cutoff:
                    altbase = 'I'  # indel
                    has_something = True
            if has_something:
                # put back into igv 1-based cooridinates
                snp_positions.append(currpos + 1)
                snp_alt_bases.append(altbase)


def check_snps_in_target_region_for_cut_sites(
    snpdict: dict[str, dict[str, Any]],
    enzyme_descriptors: dict[str, EnzymeDescriptor],
    ref_genome: ReferenceGenome,
    bamfiles: Union[str, list[str]],
    config: dict[str, Any],
    bad_enzyme_choices: dict[str, list[str]]
) -> dict[str, list[str]]:
    """Check SNPs in target region for cut sites."""
    for s in snpdict:
        if snpdict[s]['status'] != 'pass':
            continue

        target_pos = snpdict[s]['pos']
        chrom = snpdict[s]['chrom']
        strand = snpdict[s]['strand']
        if strand == 'top':
            target_start = max(target_pos + config['frag_end_range'][0], 1)
            target_end = target_pos + config['good_cut_range'][1]
        else:
            target_start = max(target_pos - config['good_cut_range'][1], 1)
            target_end = target_pos - config['frag_end_range'][0]
        target_region = f"{chrom}:{target_start}-{target_end}"

        [snp_positions, _seq_positions, snp_alt_bases, _region_ref_seq] = \
            snps_or_indels_in_region(
                bamfiles, target_region,
                ref_genome,
                basequal_cutoff=config['basequal_cutoff'],
                vaf_cutoff=config['vaf_cutoff'],
                indelaf_cutoff=config['indelaf_cutoff'],
                var_count_cutoff=config['var_count_cutoff'],
                indel_count_cutoff=config['indel_count_cutoff'])
        non_target_snps = [x for x in snp_positions if x != target_pos]

        # Save indel positions to check and warn later
        snpdict[s]['indel_positions'] = [
            x for x, y in zip(snp_positions, snp_alt_bases) if y == 'I']
        if len(snpdict[s]['indel_positions']) > 0:
            logger.info(
                "Indel positions: %s; %s",  s, snpdict[s]['indel_positions'])

        if config['drop_snps_in_full_target']:
            if len(non_target_snps) > 0:
                logger.info("%s: found SNPs in nearby region", s)
                snpdict[s]['status'] = 'drop'
                snpdict[s]['drop_reason'] = 'snps_in_target_region'
                continue

        if config['drop_indel_in_full_target']:
            if 'I' in snp_alt_bases:
                snpdict[s]['status'] = 'drop'
                snpdict[s]['drop_reason'] = 'indel_in_target_region'
                continue
        else:
            for x, y in zip(non_target_snps, snp_alt_bases):
                ref_context = ref_genome.get_sequence(chrom, x-7, x+6)
                alt_context = ref_genome.get_sequence(chrom, x-7, x-1) + \
                    y + \
                    ref_genome.get_sequence(chrom, x, x+6)
                for ename, edesc in enzyme_descriptors.items():
                    ref_hit = check_sequence_for_cut_site(
                        ref_context, edesc.motif)
                    alt_hit = check_sequence_for_cut_site(
                        alt_context, edesc.motif)
                    if (alt_hit and not ref_hit):
                        logger.info(
                            "%s: SNP in target region introduces "
                            "cut site for %s, %s", s, ename, edesc.motif)
                        bad_enzyme_choices[s].append(ename)
                    if (ref_hit and not alt_hit):
                        logger.info(
                            "%s: SNP in target region removes "
                            "cut site for %s, %s", s, ename, edesc.motif)
                        bad_enzyme_choices[s].append(ename)

    print("Bad enzyme choices due to cut sites introduced by SNPs")
    print(bad_enzyme_choices)
    return bad_enzyme_choices
