'''
Extracts only the target variant bases from the full base report
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

# Extract only the variant bases to a new report
# Could add additional calculations here if necessary

input_file = snakemake.input.combined_report
output_file = snakemake.output.variant_report
outfile = open(output_file,"w")

infile=open(input_file,"r")
linecounter=1
for x in infile:
    line = x.strip().split()
    if linecounter==1: # header
        outfile.write(x)
    elif int(line[0])== 2:
        outfile.write(x)
    linecounter+=1
infile.close()
outfile.close()


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
