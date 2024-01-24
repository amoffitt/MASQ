'''
Collapses tag to combine tags with one base error
Updates all counters to reflect rolled-up tags
'''

import os
import numpy as np
import gzip
from collections import Counter, defaultdict
import fileinput
import operator
import pickle
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from masq_helper_functions import cluster_rollup2

from masq.utils.io import tabprint

########################################################################

# Start timer
t0_all = time.time()

########################################################################

# Redirect prints to log file
old_stdout = sys.stdout
log_file = open(str(snakemake.log),"w")
sys.stdout = log_file
sys.stderr = log_file

########################################################################

# Filenames and parameters
vt_counter_filename = snakemake.input.vt_counter
vt_seq_counter_filename = snakemake.input.vt_seq_counter
flip_counter_filename = snakemake.input.flip_counter

DNA_INPUT_NANOGRAMS = snakemake.params.dna_input_ng

# WHICH REGION ARE WE PROCESSING
REGION = snakemake.params.region
# Sample name
sample = snakemake.params.sample

# Output report file
outfilename = snakemake.output.tagcounts
outfile = open(outfilename,"w")

# New counters
new_vt_counter_filename = snakemake.output.vt_counter
new_vt_seq_counter_filename = snakemake.output.vt_seq_counter
new_flip_counter_filename = snakemake.output.flip_counter

########################################################################

# Load sequence data
t0 = time.time()
print("loading sequence data...")
vt_counter = pickle.load(open(vt_counter_filename, 'rb')) # only one region
vt_seq_counter = pickle.load(open(vt_seq_counter_filename, 'rb'))
flip_counter = pickle.load(open(flip_counter_filename, 'rb'))
print("data loaded in %0.2f seconds" % (time.time() - t0))

########################################################################

# Do cluster rollup 
new_vt_counter, new_vt_seq_counter, unique_list, unique_count, match_dict, new_flip_counter = cluster_rollup2(vt_counter,vt_seq_counter, flip_counter, show_progress=False)

########################################################################

# Save new counters
pickle.dump(new_vt_counter, open(new_vt_counter_filename, 'wb'), pickle.HIGHEST_PROTOCOL)
pickle.dump(new_vt_seq_counter, open(new_vt_seq_counter_filename, 'wb'), pickle.HIGHEST_PROTOCOL)
pickle.dump(new_flip_counter, open(new_flip_counter_filename, 'wb'), pickle.HIGHEST_PROTOCOL)

########################################################################

# Calculate tag counts for report
num_obs_tags = sum(vt_counter.values())
num_uniq_tags = len(vt_counter)
num_collapsed_tags = len(new_vt_counter)
########################################################################

# Report on original tags, unique tags, and rolled up tags
outfile.write(tabprint(["Region","Observed Tags", "Unique Tags", "Rolled-Up Tags",
                        "Avg Reads Per Unique Tag","Avg Reads Per Rolled-Up Tag",
                        "Fraction of Unique Tags that are Rolled-Up","Yield: Rolled-Up Tags"])+"\n")
rolledupyield=float(num_collapsed_tags)/( DNA_INPUT_NANOGRAMS/(3.59*0.001) ) 
counts = [int(REGION),num_obs_tags,num_uniq_tags,num_collapsed_tags,
          float(num_obs_tags)/max(1,num_uniq_tags), float(num_obs_tags)/max(1,num_collapsed_tags),
          float(num_collapsed_tags)/max(1,num_uniq_tags),rolledupyield]
outfile.write(tabprint(counts)+"\n")

########################################################################

# End timer
t1 = time.time()
td = (t1 - t0) / 60
print("done in %0.2f minutes" % td)

########################################################################
# Put standard out back...

sys.stdout = old_stdout
log_file.close()

########################################################################
