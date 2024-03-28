import subprocess
from typing import Any
from collections import Counter


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
            if (score >= (plen - config['blat_num_mismatches'])):
                print("%s - %d - %d - %d" % (s, plen, score, gaps))
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