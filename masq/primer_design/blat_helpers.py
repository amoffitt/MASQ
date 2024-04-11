import subprocess
from typing import Any
from collections import Counter

from masq.primer_design.primer3_helpers import get_primer_coordinates
from masq.primer_design.snp import snps_or_indels_in_region
from masq.utils.reference_genome import ReferenceGenome


def run_blat(
    inputfile: str,
    outputfile: str,
    config: dict[str, Any],
    mode: str = 'primers'
) -> None:
    ref = config['ref_fa']
    blat = "blat"  # if installed via conda or on path
    if mode == 'primers':
        cmd = \
            f"{blat} {ref} {inputfile} {outputfile} " \
            f"-tileSize={config['tileSize']} " \
            f"-stepSize={config['stepSize']} " \
            f"-minIdentity={config['minIdentity']} " \
            f"-minScore={config['minScore']} " \
            f"-maxIntron={config['maxIntron']} " \
            f"-noHead"
    else:
        # full length blat
        cmd = f"{blat} {ref} {inputfile} {outputfile} " \
            f"-tileSize={config['tileSize']} " \
            f"-stepSize={config['stepSize']} " \
            f"-minIdentity={config['minIdentity_full']} " \
            f"-minScore={config['minScore_full']} " \
            f"-maxIntron={config['maxIntron_full']} " \
            f"-noHead"
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) as p:
        _out, _err = p.communicate()


def process_blat_results(
    blatresultfile: str,
    config: dict[str, Any]
) -> dict[str, Any]:
    # Process BLAT results
    # Counter for number of hits per sequence
    blat_hits: dict[str, Any] = {}

    with open(blatresultfile, 'r') as blatr:
        for line in blatr:
            s = line.split()[9]
            s_split = s.split("_")
            snpid = "_".join(s_split[0:3])
            leftright = s_split[4]
            primerid = s_split[5]
            gaps = int(line.split()[6])
            # Only count entry if score is X away from length of primer(query)
            plen = int(line.split()[10])
            score = int(line.split()[0])
            if score >= (plen - config['blat_num_mismatches']):
                print(f"{s} - {plen} - {score} - {gaps}")
                if snpid in blat_hits:
                    if leftright in blat_hits[snpid]:
                        blat_hits[snpid][leftright].update([primerid])
                    else:
                        blat_hits[snpid][leftright] = Counter()
                        blat_hits[snpid][leftright].update([primerid])
                else:
                    blat_hits[snpid] = {}
                    blat_hits[snpid][leftright] = Counter()
                    blat_hits[snpid][leftright].update([primerid])
    return blat_hits


def find_valid_pairs_of_primers_based_on_blat_hits(
    snpdict: dict[str, dict[str, Any]],
    primer3results: dict[str, Any],
    blat_hits: dict[str, Any],
    ref_genome: ReferenceGenome,
    config: dict[str, Any],
) -> dict[str, Any]:
    """Look for valid pairs of primers based on blat hits.

    Also check for SNPs in primer pairs
    """
    valid_primer_pairs = {}
    best_primer_pair = {}
    bam = config["wgs_bam"]

    for snpid, counts in blat_hits.items():
        perfectfound = False
        perfectleft = []
        perfectright = []
        okleft = []
        okright = []

        print(snpid)
        print(counts)
        print("")

        if 'LEFT' in counts:
            for pid, ct in counts['LEFT'].items():
                if ct == config['blat_perfect_num_hits']:
                    perfectleft.append(pid)
                    print(f"Perfect left: {snpid} - {pid}")
                elif ct < config['blat_ok_num_hits']:
                    okleft.append(pid)
                    print(f"OK left: {snpid} - {pid}")
        if 'RIGHT' in counts:
            for pid, ct in counts['RIGHT'].items():
                if ct == config['blat_perfect_num_hits']:
                    perfectright.append(pid)
                    print(f"Perfect right: {snpid} - {pid}")
                elif ct < config['blat_ok_num_hits']:
                    okright.append(pid)
                    print(f"OK right: {snpid} - {pid}")

        # check for perfect pair
        perfectpairs = list(set(perfectleft).intersection(perfectright))
        if len(perfectpairs) > 0:
            valid_primer_pairs[snpid] = perfectpairs
            ps = list(valid_primer_pairs[snpid])
            # check perfect pairs for snps in primers before skipping next step
            for p in ps:
                (primerstring_l, primerstring_r) = get_primer_coordinates(
                    p, snpid, primer3results, snpdict)
                [snp_positions_l, _seq_positions,
                 _snp_alt_bases, _region_ref_seq] = \
                    snps_or_indels_in_region(
                        bam, primerstring_l,
                        ref_genome,
                        basequal_cutoff=config['basequal_cutoff'],
                        vaf_cutoff=config['vaf_cutoff'],
                        indelaf_cutoff=config['indelaf_cutoff'],
                        var_count_cutoff=config['var_count_cutoff'],
                        indel_count_cutoff=config['indel_count_cutoff'])
                [snp_positions_r, _seq_positions,
                 _snp_alt_bases, _region_ref_seq] = \
                    snps_or_indels_in_region(
                        bam, primerstring_r,
                        ref_genome,
                        basequal_cutoff=config['basequal_cutoff'],
                        vaf_cutoff=config['vaf_cutoff'],
                        indelaf_cutoff=config['indelaf_cutoff'],
                        var_count_cutoff=config['var_count_cutoff'],
                        indel_count_cutoff=config['indel_count_cutoff'])
                if (len(snp_positions_l) > 0) or (len(snp_positions_r) > 0):
                    valid_primer_pairs[snpid].remove(p)
                    print(f"Found SNP in primer pair: {p}")
                    print(f"Left primer: {primerstring_l}")
                    print(f"Right primer: {primerstring_r}")
            if len(valid_primer_pairs[snpid]) > 0:
                perfectfound = True

        if not perfectfound:
            # check for one perfect and one ok
            ok_perf_pairs = list(set(perfectleft).intersection(okright))
            ok_perf_pairs.extend(list(set(perfectright).intersection(okleft)))

            best_pairs = []
            # min hits combined across 2 primers
            m = config['blat_perfect_num_hits'] + config['blat_ok_num_hits']
            for p in ok_perf_pairs:
                m_obs = blat_hits[snpid]['LEFT'][p] + \
                    blat_hits[snpid]['RIGHT'][p]
                if m_obs < m:
                    m = m_obs
                    best_pairs = [p]
                elif m_obs == m:
                    best_pairs.append(p)

            valid_primer_pairs[snpid] = best_pairs

        print(valid_primer_pairs[snpid])
        print("Has valid primer pairs before SNP checking")
        # Further selection based on product size (larger is better)
        if len(valid_primer_pairs[snpid]) > 0:
            m = 0
            ps = valid_primer_pairs[snpid]
            # sort so ties are broken by lowest number, which has best score
            # from primer3
            ps.sort(key=float)
            print(ps)

            ps_no_snps = list(ps)
            # Check valid primer pairs for snps in primer, drop if SNP in
            # primer
            for p in ps:
                (primerstring_l, primerstring_r) = \
                    get_primer_coordinates(p, snpid, primer3results, snpdict)
                [snp_positions_l, _seq_positions,
                 _snp_alt_bases, _region_ref_seq] = \
                    snps_or_indels_in_region(
                        bam, primerstring_l,
                        ref_genome,
                        basequal_cutoff=config['basequal_cutoff'],
                        vaf_cutoff=config['vaf_cutoff'],
                        indelaf_cutoff=config['indelaf_cutoff'],
                        var_count_cutoff=config['var_count_cutoff'],
                        indel_count_cutoff=config['indel_count_cutoff'])
                [snp_positions_r, _seq_positions,
                 _snp_alt_bases, _region_ref_seq] = \
                    snps_or_indels_in_region(
                        bam, primerstring_r,
                        ref_genome,
                        basequal_cutoff=config['basequal_cutoff'],
                        vaf_cutoff=config['vaf_cutoff'],
                        indelaf_cutoff=config['indelaf_cutoff'],
                        var_count_cutoff=config['var_count_cutoff'],
                        indel_count_cutoff=config['indel_count_cutoff'])
                if (len(snp_positions_l) > 0) or (len(snp_positions_r) > 0):
                    print(f"Found SNP in primer pair: {p}")
                    print(f"Left primer: {primerstring_l}")
                    print(f"Right primer: {primerstring_r}")
                    ps_no_snps.remove(p)

            # Still has valid primer options
            if len(ps_no_snps) > 0:
                m = 0
                for p in ps_no_snps:
                    prodsize = int(
                        primer3results[snpid][f"PRIMER_PAIR_{p}_PRODUCT_SIZE"])
                    print(f"{prodsize}: {p}")
                    if prodsize > m:
                        m = prodsize
                        bestprimer = p
                        print("bigger product")
                best_primer_pair[snpid] = bestprimer
                print(bestprimer)
        else:
            snpdict[snpid]['status'] = 'drop'
            snpdict[snpid]['drop_reason'] = 'blat_hits'

    return best_primer_pair
