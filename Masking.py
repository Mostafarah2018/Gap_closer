import os
import subprocess
from Bio import SeqIO
import pandas as pd
from itertools import chain
import sys as syst



seq1 = "AAAACCCTTAGCAAATAAGCTTAGAATATAATAAAGCGCGAATTAAAA"
seq2 = "TTTTAATTCGCGCTTTATTATATTCTAAGCTTATTTGCTAAGGGTTTT"

# Reading fastq files in the directory
files = []
files_Revco_list = []
files_Revco_set = set()
for file in os.listdir(os.getcwd()):
    if file.endswith(".fastq"):
        files.append(file)

# Counting the number of reads in each raw reads files


def count_seq(file):
    count = 0
    suffix = "fasta" if file.endswith("fasta") else "fastq"
    parsed_iter = SeqIO.parse(file, suffix)
    count = sum(1 for i in parsed_iter)
    return count

files_CC_list=[]
# Grep reads with CC strand
for file in files:

    name = file[:-6]
    fastq_parser = SeqIO.parse(file, "fastq")
    wanted = [rec for rec in fastq_parser if seq1 in rec]
    SeqIO.write(wanted, name+"_CC.fastq", "fastq")  # to be deleted

    # first out put fasta
    SeqIO.convert(name+"_CC.fastq", 'fastq', name+"_CC.fasta", 'fasta')
    os.remove(name+"_CC.fastq")
    files_CC_list.append(name+"_CC.fasta")
    
    
    # and Reverse complementing the results as fasta file
    # name_GG=name+"_GG.fastq"
    fastq_parser = SeqIO.parse(file, "fastq")
    wanted = (rec.reverse_complement(id=rec.id, description="Reverse complement")
              for rec in fastq_parser if seq2 in rec)
    
    SeqIO.write(wanted, name+"_GG_RevCo.fasta", "fasta")
    files_Revco_list.append(name+"_GG_RevCo.fasta")
    if name[:11] not in files_Revco_set:
        files_Revco_set.add(name[:11])
        

print("Combining telomeric reads ....")
# Combine all 4 tel files for two R1 and R2 raw reads
Comb_files = []
for name in files_Revco_set:
    temp_list = [s for s in files_CC_list+files_Revco_list if name in s]
    parsed_iter = [SeqIO.parse(s, "fasta") for s in temp_list]
    #print("name=", name, "list length=", len(temp_list))
    #print(temp_list)
    output = chain()
    for gen in parsed_iter:
        output = chain(output, gen)
    SeqIO.write(output, "Combined_"+name+".fasta", "fasta")
    Comb_files.append("Combined_"+name+".fasta")

##############################
# Remove all duplicate files based on their sequences
# read all Comnined files

# find duplicates based on their sequences
combs_uniq_mask = []
for fil in Comb_files:

    name = fil[:-6]
    uniq_seqrecs = []
    seen_seqs = []
    for rec in SeqIO.parse(fil, 'fasta'):
        if rec.seq not in seen_seqs:
            seen_seqs.append(rec.seq)
            uniq_seqrecs.append(rec)

    SeqIO.write(uniq_seqrecs, name+"_uniq.fasta", 'fasta')
    combs_uniq_mask.append(name+"_uniq.fasta")

# use replace instead of  runing sed command TT
####### FOAD >>>>>>>>>> please remove all the sequence until the end of the marker 




combs_uniq = []
for J in combs_uniq_mask:

    with open(J, "r") as inp:
        data = inp.read().replace("AAAACCCTTAGCAAATAAGCTTAGAATATAATAAAGCGCGAATTAAAA", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")

    output = J[:-6]+"_Masked.fasta"
    with open(output, "w") as out:
        out.write(data)
    combs_uniq.append(output)
