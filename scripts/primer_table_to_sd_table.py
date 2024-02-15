'''
Converts primer table output from primer design functions to input table needed for MASQ analysis pipeline
'''
import sys
import pickle

from masq.utils.io import tabprint
from masq.utils.seqs import reverse_complement
from masq.utils.reference_genome import ReferenceGenome


primer_table = snakemake.input.oldtable

# sequence dictionary
ref_genome_filename = snakemake.params.ref_genome
ref_genome = ReferenceGenome(ref_genome_filename).open()


c=0
loci=0
with open(snakemake.output.newtable,'w') as fout:
    with open(primer_table,'r') as f:
        for line in f:
            if c==0:
                header=line.strip().split("\t")
                fout.write(tabprint(["loc","chr","posi","specific-primer-1","specific-primer-2","trimmed-target-seq","target_locs","ref-alt_allele"])+"\n")
            else:
                X=line.strip().split("\t")

                #chrom=X[header.index("chrom")][3:]
                chrom=X[header.index("chrom")]
                pos=int(X[header.index("pos")])
                primer1=X[header.index("downstream_primerseq")]
                primer2=X[header.index("cutadj_primerseq")]
                strand=X[header.index("strand")]
                refbase=X[header.index("ref_trinuc")][1]
                altbase=X[header.index("alt_trinuc")][1]
                ref_alt=refbase+"_"+altbase

                # get target sequence coordinates based on primer Coordinates
                # pull sequence from dictionary
                # reverse complement depending on strand
                # target position is relative to target seq
                primercoords1=X[header.index("cutadj_primer_coordinates")].split(':')[1].split('-')
                primercoords1.extend(X[header.index("downstream_primer_coordinates")].split(':')[1].split('-'))
                primercoords = [int(x) for x in primercoords1]
                primercoords.sort()

                target_start = primercoords[1]
                target_end = primercoords[2]-1

                target_seq_ref = ref_genome.get_sequence(
                    chrom, target_start, target_end)

                if strand=="bottom":
                    target_seq_stranded = reverse_complement(target_seq_ref)
                    relative_pos = target_end - pos
                elif strand=="top":
                    target_seq_stranded = target_seq_ref
                    relative_pos = pos - target_start - 1

                fout.write(tabprint([loci,chrom,pos,primer1,primer2,target_seq_stranded,relative_pos,ref_alt])+"\n")

                loci=loci+1
            c=c+1
