'''
Combines within-template error statistics across all loci
'''

import os
import numpy as np
import gzip
from collections import Counter, defaultdict
import fileinput
import operator
import dill as pickle
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

# Input is all the pickle files with error tables
# Output is 2 table files - counts, and converted to fractions

input_file_list = snakemake.input.within_tag_error_pickles
output_file1 = snakemake.output.within_tag_table1
output_file2 = snakemake.output.within_tag_table2

outfile1 = open(output_file1,"w")
outfile2 = open(output_file2,"w")
########################################################################
NUCS = ["A", "C", "G", "T"]
REF_TRINUCS = [x+y+z for x in NUCS for y in NUCS for z in NUCS]
PAIRED_TRINUCS=[ (a , a[0]+x+a[2]) for a in REF_TRINUCS for x in NUCS if x!=a[1] ]
ERR_RANGE1 = np.arange(-0.1,1.0,0.1)

# Intiialize table
WITHIN_TAG_ERRS_ALL = np.zeros(shape=(len(PAIRED_TRINUCS),len(ERR_RANGE1),2), dtype=int)

for f in input_file_list:
    # Add to current table
    WITHIN_TAG_ERRS_REGION = pickle.load(open(f,'rb'))
    WITHIN_TAG_ERRS_ALL = WITHIN_TAG_ERRS_ALL  +  WITHIN_TAG_ERRS_REGION

# Original table
for i,pt in enumerate(PAIRED_TRINUCS):
    outfile1.write(tabprint(['R1']+list(pt)+list(WITHIN_TAG_ERRS_ALL[i,:,0]))+'\n')
    outfile1.write(tabprint(['R2']+list(pt)+list(WITHIN_TAG_ERRS_ALL[i,:,1]))+'\n')

# Convert to fractions, skip the 0 error bin
for i,pt in enumerate(PAIRED_TRINUCS):
    fracR1=WITHIN_TAG_ERRS_ALL[i,1:,0]/WITHIN_TAG_ERRS_ALL[i,1:,0].sum(keepdims=True)
    fracR2=WITHIN_TAG_ERRS_ALL[i,1:,1]/WITHIN_TAG_ERRS_ALL[i,1:,1].sum(keepdims=True)
    outfile2.write(tabprint(['R1']+list(pt)+list(fracR1))+'\n')
    outfile2.write(tabprint(['R2']+list(pt)+list(fracR2))+'\n')


########################################################################

# Close output files
outfile1.close()
outfile2.close()

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
