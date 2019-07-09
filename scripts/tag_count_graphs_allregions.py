'''
For all regions, plot number of reads per tag/template
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


from masq_helper_functions import tabprint
from masq_helper_functions import plot_number_of_reads_per_tag
from masq_helper_functions import plot_at_least_x_reads_per_tag
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
vt_counter_filenames = snakemake.input.vt_counters
NREG = len(vt_counter_filenames)

# Sample name
sample = snakemake.params.sample

# Output report file
outfilename = snakemake.output.tagcounts
outfile = open(outfilename,"w")
########################################################################

# Load sequence data
t0 = time.time()
print("loading sequence data...")
vt_counters=[]
vt_counters = [Counter() for _ in range(NREG)]
for r in range(NREG):
    vt_counters[r]=pickle.load(open(vt_counter_filenames[r], 'rb'))
print("data loaded in %0.2f seconds" % (time.time() - t0))

########################################################################

# tabulate distribution of reads per tag
reads_per_tag = Counter()
Nreads_total = 0
Ntags_total = 0

for r in range(NREG):
    vt_counter = vt_counters[r]
    Ntags = len(vt_counter)
    Ntags_total += Ntags
    for tag, tag_count in vt_counter.items():
        reads_per_tag[tag_count]+=1
        Nreads_total += tag_count

numreads=np.array(list(reads_per_tag.keys()))
numtags=np.array(list(reads_per_tag.values()))

########################################################################

# Plot reads per tag
plot_number_of_reads_per_tag(
                                numreads,
                                numtags,
                                Nreads_total,
                                Ntags_total,
                                sample=sample,
                                filename=snakemake.output.plot1,
                                logscale=True)
plot_at_least_x_reads_per_tag(
                                numreads,
                                numtags,
                                Nreads_total,
                                Ntags_total,
                                sample=sample,
                                filename=snakemake.output.plot2,
                                logscale=True,
                                maxcount=100)
########################################################################

# Write bar graph counts out to file
header = ["Region","Reads Per Tag","Number of Tags"]
outfile.write(tabprint(header)+"\n")
for i in range(len(numreads)):
    outfile.write(tabprint(["Combined",numreads[i],numtags[i]])+"\n")

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
