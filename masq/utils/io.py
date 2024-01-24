from typing import TextIO, Any


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
