'''
Plot results of tag roll-up / collapsing for all loci
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

########################################################################

# Start timer
t0 = time.time()

########################################################################

# Redirect prints to log file
old_stdout = sys.stdout
log_file = open(str(snakemake.log),"w")
sys.stdout = log_file
sys.stderr = log_file

########################################################################
input_file = snakemake.input.combined_report
output_file = snakemake.output.plot1
sample = snakemake.params.sample
########################################################################

# Process report file to get counts
num_obs_tags = []
num_uniq_tags = []
num_collapsed_tags = []

infile = open(input_file,"r")
linecount=1
regioncount = 0
for x in infile:
    line = x.strip().split()
    if linecount>1:
        region = line[0]
        regioncount+=1
        num_obs_tags.append(int(line[1]))
        num_uniq_tags.append(int(line[2]))
        num_collapsed_tags.append(int(line[3]))
    linecount+=1
infile.close()
print(num_obs_tags)
print(num_uniq_tags)
print(num_collapsed_tags)

########################################################################

# Graph for each region, total tags, unique tags, collapsed tags
N = regioncount
fig = plt.figure(figsize=(50,10))
fig.suptitle(sample+"\n"+
             "Results of Tag Rollup - 1 error allowed", fontsize=16)
x=range(N)
width=0.25
ax = plt.subplot(111)
ax.bar(x,num_obs_tags,width,alpha=0.5,color='purple',label='Total')
ax.bar([p + width for p in x],num_uniq_tags,width,alpha=0.5,color='green',label='Unique')
ax.bar([p + width*2 for p in x],num_collapsed_tags,width,alpha=0.5,color='blue',label='Collapsed')
ax.legend(['Total','Unique','Collapsed'], loc='upper left')
ax.ticklabel_format(style='plain')
ax.get_yaxis().set_major_formatter(
    matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

plt.xticks([p + width for p in x],range(N))
plt.xlim([min(x)-width,max(x)+width*4])
plt.xlabel("Region",fontsize=14)
plt.ylabel("Tag Counts",fontsize=14)
plt.savefig(output_file, dpi=200, facecolor='w', edgecolor='w',
            papertype=None, format=None,
            transparent=False)
plt.close()


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
