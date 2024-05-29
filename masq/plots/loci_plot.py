from collections import defaultdict
import logging
from typing import Any

import numpy as np

import matplotlib
import matplotlib.pyplot as plt

import pysam

from masq.utils.seqs import BASE2INT, BASES
from masq.utils.reference_genome import ReferenceGenome


matplotlib.use('Agg')

logger = logging.getLogger(__name__)

REGION_BUFFER = 30
BUFF = 30
MIN_QUAL = 20
MIN_MAP = 40

# Plotting colors
BASE_COLORS = ["#00ABBA",
               "#9ACA3C",
               "#F26421",
               "#672D8F"]


def _collect_reads(
    sample_bam_fn: str,
    chrom: str,
    start: int,
    end: int,
) -> dict[str, list]:
    # open the bam file
    read_pair_dict = defaultdict(list)
    with pysam.AlignmentFile(sample_bam_fn, "rb") as sample_bam:

        # for each position of interest
        # get an iterator that walks over the reads that include our event
        read_iterator = sample_bam.fetch(
            chrom, start - BUFF, end + BUFF)  # ref has chr1,chr2,etc.

        # store as pairs
        for read in read_iterator:
            qname = read.query_name
            assert qname is not None

            read_pair_dict[qname].append(read)  # if pair is there
    return read_pair_dict


def _process_reads(
    read_pair_dict: dict[str, list],
    start: int,
    end: int,
) -> Any:
    # store an integer at every position for each read(pair) in interval
    # matrix (number of reads (N) by number of positions (end-start))
    read_stack = np.zeros(shape=(len(read_pair_dict), end - start), dtype=int)
    # get all read pairs
    for rcounter, reads in enumerate(read_pair_dict.values()):
        for read in reads:
            read_sequence = read.query_sequence
            read_qualities = read.query_qualities
            if read_sequence is None or read_qualities is None:
                continue
            # get_aligned_pairs:
            # returns paired list of ref seq position and read seq position
            for rindex, pindex in read.get_aligned_pairs():
                # rindex and pindex return ALL positions, including soft-clipped
                # rindex: position in read (1-150 for example)
                # pindex: position in reference genome
                # if there's a base in the read
                # and we are in range
                if rindex is None or pindex is None:
                    continue
                # Separated these lines, added pindex None check
                if not start <= pindex < end:
                    continue

                # compute the base score
                base_score = read_qualities[rindex]
                # if it's good enough
                if not base_score >= MIN_QUAL:
                    continue
                # check if there is already a value in place
                current = read_stack[rcounter, pindex - start]
                base = read_sequence[rindex]
                if base == "N":
                    baseint = 0
                else:
                    # A-1, C-2, G-3, T-4
                    baseint = BASE2INT[base] + 1
                # if there is no value, store the value
                if current == 0:
                    read_stack[rcounter, pindex - start] = baseint
                else:
                    # if there is a mismatch between the two reads,
                    # set value back to 0
                    # this value is just for the 2 paired reads
                    if current != baseint:
                        read_stack[rcounter, pindex - start] = 0
    return read_stack


def _build_summary(
    read_stack: np.ndarray,
) -> np.ndarray:
    summary_list = []
    # iterating over numpy array - by rows
    # transpose first to iterate over positions
    for x in read_stack.transpose():
        # gets counts of N,A,C,G,T per position as array
        # append to summary list
        summary_list.append(np.bincount(x, minlength=5))
    # convert summary to array, fills by row
    # drop the N count, and transpose again
    # .T tranposes, only if num dimensions>1
    return np.array(summary_list)[:, 1:].T

def _plot_bases_vs_freq(
    ax2: Any,
    base_ratio: np.ndarray,
    local_int: np.ndarray,
    seq_len: int,
) -> None:
    # Plot colored circle for each base at position vs. frequency
    for base_color, base, yvals in zip(BASE_COLORS, BASES, base_ratio):
        ax2.plot(yvals, 'o', markersize=13, label=base, color=base_color)
        base_filter = local_int == BASE2INT[base]  # same as ref base
        x = np.arange(seq_len)
        # label ref bases with  black circle
        ax2.plot(x[base_filter], yvals[base_filter], 'o',
                    markersize=13, mfc="None", mec='black', mew=2)

def _label_number_of_reads(
    ax2: Any,
    has_something: np.ndarray,
    summary: np.ndarray,
    base_cover: np.ndarray,
    base_ratio: np.ndarray,
    local_int: np.ndarray,
) -> None:
    # for non-ref sites...
    # label number of reads for non-ref base in red
    for hs_ind in np.where(has_something)[0]:
        base_counts, _ = summary[:, hs_ind], base_cover[hs_ind]
        ref_base = local_int[hs_ind]
        for bind in range(4):
            # loop through A/C/G/T, find the non-zero counts
            if bind != ref_base and base_counts[bind] > 0:
                # text_string = r"$\frac{%d}{%d}$" % (
                    # base_counts[bind],total)
                # text_string = r"%d/%d" % (base_counts[bind], total)
                text_string = f"{base_counts[bind]}"
                # add the count to the plot
                ax2.text(hs_ind + 1, base_ratio[bind, hs_ind], text_string,
                            fontsize=20, color='red', weight='semibold')



def loci_plot(
    loc_ind: int,
    image_filename: str,
    chrom: str,
    true_start: int,
    true_end: int,
    true_targets: list[int],
    strand: str,
    ref_genome: ReferenceGenome,
    sample_bams: dict[str, str],
) -> Any:
    """Plot loci for a given region."""
    # add buffer to target sequence on either side and update values
    start = true_start - REGION_BUFFER
    end = true_end + REGION_BUFFER
    targets = [x + REGION_BUFFER for x in true_targets]
    # get local sequence
    local_seq = ref_genome.get_sequence(chrom, start, end)
    seq_len = len(local_seq)

    logger.debug('Chrom: %s', chrom)
    logger.debug('start: %d', start)
    logger.debug('end: %d', end)
    logger.debug('Local_seq: %s', local_seq)

    # also in integer form
    local_int = np.array([BASE2INT[x] for x in local_seq])
    # intialize for keeping track of variants at this position
    has_something = np.zeros(seq_len, dtype=bool)

    # for each sample
    fig = plt.figure(figsize=(20, 18))
    fig.suptitle(f"{chrom}:{start}-{end} {strand}", fontsize=20)

    for ind, (sample_name, sample_bam_fn) in enumerate(sample_bams.items()):
        logger.info('Plotting locus %d', loc_ind)

        # open the bam file
        read_pair_dict = _collect_reads(sample_bam_fn, chrom, start, end)
        read_stack = _process_reads(read_pair_dict, start, end)

        summary = _build_summary(read_stack)
        # now we have base (4) by position array as summary
        # base_cover: coverage at each position (sum of A/C/G/T counts)
        base_cover = np.sum(summary, axis=0)
        # base_ratio: A/C/G/T counts over total coverage
        # aka frequency of each base
        base_ratio = summary.astype(float) / np.maximum(1, base_cover)
        # update has something
        # EDIT - 0 coverage is not has_something??
        has_something += (
            (base_ratio[local_int, np.arange(seq_len)] < 0.9)
            & (base_cover>3)
        ) # reference ratio is less than 0.9

        # Set up the plot
        ax1 = fig.add_subplot(len(sample_bams), 1, ind + 1)
        ax2 = ax1.twinx()
        ax1.tick_params(axis='both', labelsize=15)
        ax2.tick_params(axis='both', labelsize=15)

        ########################################################################
        # Plot variants in each region

        # plot the coverage first
        ax1.plot(base_cover, color='k', label='coverage', alpha=0.5)
        # draw lines for boundaries of event
        # and targets
        if strand == "+":
            ax2.axvline(REGION_BUFFER-0.5, color='g', lw=2)
            ax2.axvline(
                seq_len-REGION_BUFFER-0.5,
                color='g', linestyle='--', lw=2)
        else:
            ax2.axvline(REGION_BUFFER-0.5, color='g', linestyle='--', lw=2)
            ax2.axvline(seq_len-REGION_BUFFER-0.5, color='g', lw=2)
        for pos in targets:
            ax1.axvline(pos, color='r', lw=2)

        ax1.set_ylabel(sample_name, fontsize=25)

        _plot_bases_vs_freq(ax2, base_ratio, local_int, seq_len)
        _label_number_of_reads(
            ax2,
            has_something,
            summary,
            base_cover,
            base_ratio,
            local_int,
        )

        ax2.set_xlim(-10, 10+seq_len)
        ax2.set_ylim(0, 1)
        ax2.legend(loc="upper center", numpoints=1, ncol=4,
                    fontsize=18, columnspacing=0.8, handletextpad=-0.2)
        ax1.legend(loc="lower center", numpoints=2,
                    fontsize=18, columnspacing=0.5)

    plt.tight_layout(rect=(0, 0, 1, 0.98))
    plt.savefig(image_filename, dpi=200, facecolor='w', edgecolor='w',
                papertype=None, format=None,
                transparent=False)
    plt.close()

    ########################################################################
    # Add non-ref het/homo sites to "target" list
    rlen = true_end - true_start
    new_targets = np.where(has_something)[0] - REGION_BUFFER
    if strand == "-":
        new_targets = true_end - true_start - new_targets - 1

    # new_targets>=0 (in original target region)
    new_targets = new_targets[(new_targets >= 0) * (new_targets < rlen)]
    return list(new_targets)
