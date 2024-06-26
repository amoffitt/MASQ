'''
Filters base tables to remove loci that have been filtered for QC metrics
QC fail loci defined in masq_QC_plots.R
'''

import os
import numpy as np
import gzip
from collections import Counter, defaultdict
import fileinput
import operator
import pickle
import time
import logging

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
########################################################################

logger = logging.getLogger('qcfilter_report')

logger.info('Starting process')

########################################################################

# Input is combined base report file 
# Output is same file, with QC filtered loci removed
logger.info('Filtering report')

qc_fail_loci=[]
with open(snakemake.input.qcfail,'r') as f:
    for line in f:
        qc_fail_loci.append(line.strip())

c=0
with open(snakemake.input.base_report,'r') as f:
    with open(snakemake.output.filtered_base_report,'w') as fout:
        for line in f:
            if (c==0):
                fout.write(line)
            else:
                x=line.split()
                if (x[1] not in qc_fail_loci):
                    fout.write(line)
            c = c+1
