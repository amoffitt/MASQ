from collections import Counter, defaultdict
from typing import Any

import numpy as np


## simple class for counting pairs
def double_counter() -> list[int]:
    return [0, 0]

########################################################################
# takes as input a:
# 1) counter with tags and tag counts
# 2) dictionary with tags as keys, and counter as values ((read1,read2) as key, count as value)
# returns updated versions of these

def cluster_rollup2(
    vt_counter: Counter,
    vt_seq_counter: dict[str, Counter],
    flip_counter: dict[tuple[str, str, str], list[int]],
    show_progress: bool = False
) -> Any:

    # Unique list and counts sorted from most to least common
    unique_list = [x[0] for x in vt_counter.most_common()]
    unique_count = [x[1] for x in vt_counter.most_common()]

    # IF READS ARE IN COUNTER:
    if len(unique_list) == 0:
        return vt_counter, vt_seq_counter, [] , [], {}, flip_counter

    taglen = len(unique_list[0])

    # match_dict stores the first hamming-radius match for each entry
    # (initialized to itself)
    match_dict = np.arange(len(unique_list))
    # omit dictionary is a dictionary for each position 0 <= x < TAGLEN
    # such that when index == index,
    # omit_dictionary[i] has keys of length TAGLEN -1
    # derived from the previous entries tag[k] for k < index
    # by omitting the i-th letter,
    # and the value is the first index with that key.
    omit_dictionary: list[dict] = [{} for _ in range(taglen)]

    counter = 0

    for index, tag in enumerate(unique_list):
        counter += 1
        if show_progress and counter % 100000 == 0:
            print(counter)
        ## keep a set of values where the i-th omit matches a previous key.
        match       = set()
        for i in range(taglen):
            # print(i)
            newtag = tag[:i] + tag[(i+1):]
            try:
                ## check if it's already in the dictionary
                has_match = omit_dictionary[i][newtag]
                ## if it is, add it to the list of matches
                match.add(has_match)
            except KeyError:
                ## if not, there is a error raised
                ## and we now assign the present index to this entry
                omit_dictionary[i][newtag] = index

        # You are assigned to your top match by a sequence of hamming-1
        # distances that connect lower counts to higher counts
        matches = list(np.sort(list(match))) + [index]
        # select lowest index matching tag - then go look up what value it
        # has stored in the dictionary
        match_dict[index]  = match_dict[matches[0]]

    # Now re-process the original counters using match_dict to combine tags
    new_vt_counter: dict[str, int]     = Counter()
    new_vt_seq_counter: dict[str, Counter] = defaultdict(Counter)
    new_flip_counter: dict[tuple[str, str, str], list[int]] = \
        defaultdict(double_counter)

    for i, match_item in enumerate(match_dict):
        # go through each unique tag once

        # get the assignment for this tag
        clusterindex = match_item

        # get the tag seq that it was assigned to
        clustertag = unique_list[clusterindex]

        origtag = unique_list[i]

        # count of originally unique tag
        origcount = unique_count[i]

        # update new counter with old count, new tag
        new_vt_counter[clustertag] += origcount
        seqs_for_vt = vt_seq_counter[origtag]

        # add reads to new counter with new tag
        for (r1, r2), readcount in seqs_for_vt.most_common():
            new_vt_seq_counter[clustertag][(r1,r2)] += readcount

            new_flip_counter[(clustertag,r1,r2)][0] += \
                flip_counter[(origtag,r1,r2)][0]
            new_flip_counter[(clustertag,r1,r2)][1] += \
                flip_counter[(origtag,r1,r2)][1]

    return (
        new_vt_counter,
        new_vt_seq_counter,
        unique_list,
        unique_count,
        match_dict,
        new_flip_counter,
    )
