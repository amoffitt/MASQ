import logging
from typing import TextIO, Any

from masq.utils.seqs import reverse_complement
from masq.utils.reference_genome import ReferenceGenome

logger = logging.getLogger(__name__)


def tabprint(line: list) -> str:
    return "\t".join(map(str, line))


SNV_HEADER = [
    'loc', 'chr', 'posi', 'specific-primer-1',
    'specific-primer-2', 'trimmed-target-seq', 'target_locs',
    'ref-alt_allele', 'indel_start', 'indel_length', 'indel_seq',
    'add-targets', 'strand', 'fragment-start', 'fragment-end'
]


def load_snv_table(infile: TextIO) -> dict[str, Any]:
    """Read locus/SNV table from file."""
    snv_table: dict[str, Any] = {}
    headings = infile.readline().strip().split("\t")
    for h in headings:
        if h not in SNV_HEADER:
            raise IOError(f'Unknown column header: {h}')
        snv_table[h] = []

    for line in infile:
        entry = line.strip("\n").split("\t")
        for i, e in enumerate(entry):
            snv_table[headings[i]].append(e)
    infile.close()
    return snv_table


def write_snv_table(table: dict[str, Any], outfile: TextIO) -> None:
    """Write locus/SNV table to file."""
    outfile.write(tabprint(list(table.keys()))+"\n")
    for i in range(len(table['loc'])):
        line = []
        for _, val in table.items():
            line.append(val[i])
        outfile.write(tabprint(line)+"\n")


def process_target_info(
    snv_info: dict[str, Any],
    ref_genome: ReferenceGenome,
) -> list[list]:
    """Process input file to extend it with strand and coordinates of seq."""
    logger.info('Parsing specific SNV info fields')

    target_info = []
    target_locs_array = [
        list(map(int, x.split(";") ) ) if len(x)>0 else []
        for x in snv_info['target_locs']
    ]

    for i,loc in enumerate(snv_info['loc']):
        chrom=snv_info['chr'][i]
        pos = int(snv_info['posi'][i])
        target_locs = target_locs_array[i]
        aml_loc = target_locs[0]
        logger.debug('Position: %d', pos)
        logger.debug('AML Loc: %d', aml_loc)
        int_seq = snv_info['trimmed-target-seq'][i]
        length = len(int_seq)

        if chrom=="0":
            strand="+"
            start=1
            end=100
            targets=target_locs
        elif ('strand' in snv_info.keys()) and \
                ('fragment-start' in snv_info.keys()) and \
                ('fragment-end' in snv_info.keys()) and \
                ('add-targets' in snv_info.keys()):
            strand=snv_info['strand'][i]
            start=snv_info['fragment-start'][i]
            end=snv_info['fragment-end'][i]
            end=snv_info['add-targets'][i]
        else:
            logger.info(
                'Inferring strand, start, and end from match to '
                'reference genome')
            # identifying start and end genomic positions of the targeted region
            # and getting target sequence from ref genome
            top_start = pos - 1 - aml_loc
            top_end = top_start + length
            top_match = ref_genome.get_sequence(chrom, top_start, top_end)
            # if the target sequence was on the bottom strand...
            bottom_start = pos - (length - aml_loc)
            bottom_end = bottom_start + length
            bottom_match = reverse_complement(
                ref_genome.get_sequence(chrom, bottom_start, bottom_end))

            logger.debug('Target seq : %s', int_seq)
            logger.debug('Top match  : %s', top_match)
            logger.debug('Btm match  : %s', bottom_match)

            # check which strand the sequence came from and set coordinates
            if top_match == int_seq:
                logger.debug('Locus %s: positive strand', loc)
                strand = "+"
                start = top_start
                end = top_end
                targets = target_locs
            else:
                logger.debug('Locus %s: negative strand', loc)
                strand = "-"
                start = bottom_start
                end = bottom_end
                targets = [length - x - 1 for x in target_locs]

        target_info.append([chrom, strand, start, end, targets])
    return target_info


def extend_snv_info_with_target_info(
    snv_info: dict[str, Any],
    target_info: list[list],
    all_targets: list[list[int]],
) -> dict[str, Any]:
    snv_info['add-targets']=[]
    snv_info['strand']=[]
    snv_info['fragment-start']=[]
    snv_info['fragment-end']=[]
    print("Start loop")

    for i, (new_targets, more_info) in enumerate(zip(all_targets, target_info)):
        print(i)
        prev_targets = list(map(int, snv_info['target_locs'][i].split(";")))
        add_targets = [x for x in new_targets if x not in prev_targets]
        snv_info['add-targets'].append(";".join(list(map(str, add_targets))))
        snv_info['strand'].append(more_info[1])
        snv_info['fragment-start'].append(more_info[2])
        snv_info['fragment-end'].append(more_info[3])
    print("End loop")
    print(snv_info)

    return snv_info
