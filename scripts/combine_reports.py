'''
Combines reports from individual loci into one report
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

# Input is list of region specific report files which all have a header
# Output is combined report files with one header

input_file_list = snakemake.input.region_reports
output_file = snakemake.output.combined_report
outfile = open(output_file,"w")

counter=0
for f in input_file_list:
    counter +=1
    infile=open(f,"r")
    linecounter=1
    for line in infile:
        if linecounter==1:
            if counter==1: # first file header
                outfile.write(line)
        else:
            outfile.write(line)
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
