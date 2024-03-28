import sys
import subprocess
from typing import Any
from textwrap import dedent


def write_primer3_input_file(
    fn: str,
    snpid: str,
    templateseq: str,
    strand: str,
    dist: int,
    config: dict[str, Any],
) -> None:
    # Prepare PRIMER3 input file and run PRIMER3 on each SNP
    # SEQUENCE_TARGET should be position, followed by length!
    primer3text = dedent(
        f"""SEQUENCE_ID={snpid}
SEQUENCE_TEMPLATE={templateseq}
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE={config['PRIMER_OPT_SIZE']}
PRIMER_MIN_SIZE={config['PRIMER_MIN_SIZE']}
PRIMER_MAX_SIZE={config['PRIMER_MAX_SIZE']}
PRIMER_PRODUCT_SIZE_RANGE={config['PRIMER_PRODUCT_SIZE_RANGE'][0]}-{config['PRIMER_PRODUCT_SIZE_RANGE'][1]}
PRIMER_PRODUCT_OPT_SIZE={config['PRIMER_PRODUCT_OPT_SIZE']}
PRIMER_MIN_TM={config['PRIMER_MIN_TM']}
PRIMER_MAX_TM={config['PRIMER_MAX_TM']}
PRIMER_OPT_TM={config['PRIMER_OPT_TM']}
PRIMER_PAIR_MAX_DIFF_TM={config['PRIMER_PAIR_MAX_DIFF_TM']}
PRIMER_MIN_GC={config['PRIMER_MIN_GC']}
PRIMER_MAX_GC={config['PRIMER_MAX_GC']}
PRIMER_MAX_HAIRPIN_TH={config['PRIMER_MAX_HAIRPIN_TH']}
PRIMER_MAX_POLY_X={config['PRIMER_MAX_POLY_X']}
PRIMER_NUM_RETURN={config['PRIMER_NUM_RETURN']}
PRIMER_TM_FORMULA=0
PRIMER_SALT_CORRECTIONS=0
PRIMER_EXPLAIN_FLAG=1
PRIMER_THERMODYNAMIC_PARAMETERS_PATH={config['primer3_thermo_param_folder']}
=""")

    if strand == 'top':
        forcetext = f"SEQUENCE_FORCE_RIGHT_START={len(templateseq) - 1}\n"
        targettext = f"SEQUENCE_TARGET={len(templateseq)-dist-1},2\n"
    else:
        forcetext = "SEQUENCE_FORCE_LEFT_START=0\n"
        targettext = f"SEQUENCE_TARGET={dist-1},2\n"
    primer3text_plusforce = targettext + forcetext + primer3text

    with open(fn, 'w') as f:
        f.write(primer3text_plusforce)
    print("writing file")
    print(snpid)


def run_primer3(
    snpdict: dict[str, dict[str, Any]],
    blatqueryfile: str,
    primer3file: str,
    primer3cmd: str,
    config: dict[str, Any],

) -> dict[str, Any]:
    primer3results: dict[str, Any] = {}
    with open(blatqueryfile, 'w') as blatf:
        for snpid, snpinfo in snpdict.items():
            if snpinfo['status'] != 'pass':
                continue

            print("snp")
            write_primer3_input_file(
                primer3file,
                snpid,
                snpdict[snpid]['target_seq_for_primer_search'],
                snpdict[snpid]['strand'],
                snpdict[snpid]['dist_from_mut_to_upstream_cut'],
                config)
            print("")
            with subprocess.Popen(
                    f"{primer3cmd} {primer3file}",
                    shell=True,
                    stdout=subprocess.PIPE) as p:
                primer3out, _err = p.communicate()
                print(primer3out.decode('ascii'))
                print("")
                sys.stdout.flush()

                # Store all the primer3 results
                primer3results[snpid] = {}
                for line in primer3out.decode('ascii').split('\n'):
                    if line.startswith('PRIMER'):
                        t, val = line.split("=")
                        primer3results[snpid][t] = val

                if "PRIMER_PAIR_NUM_RETURNED=0" in primer3out.decode('ascii'):
                    snpdict[snpid]['status'] = 'drop'
                    snpdict[snpid]['drop_reason'] = 'primer3_nonefound'

                    snpdict[snpid]['left_primer_explanation'] = \
                        primer3results[snpid]['PRIMER_LEFT_EXPLAIN']
                    snpdict[snpid]['right_primer_explanation'] = \
                        primer3results[snpid]['PRIMER_RIGHT_EXPLAIN']
                elif "PRIMER_ERROR" in primer3out.decode('ascii'):
                    snpdict[snpid]['status'] = 'drop'
                    snpdict[snpid]['drop_reason'] = 'primer3_error_seelog'

                else:
                    # primer pairs found!
                    for i in range(config['PRIMER_NUM_RETURN']):
                        t = f"PRIMER_LEFT_{i}_SEQUENCE"
                        if t in primer3results[snpid].keys():
                            seq = primer3results[snpid][t]
                            blatf.write(">"+snpid+"_"+t+"\n")
                            blatf.write(seq+"\n")

                        t = f"PRIMER_RIGHT_{i}_SEQUENCE"
                        if t in primer3results[snpid].keys():
                            seq = primer3results[snpid][t]
                            blatf.write(">"+snpid+"_"+t+"\n")
                            blatf.write(seq+"\n")
    return primer3results
