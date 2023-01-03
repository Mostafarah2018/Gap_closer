import pysam
from Bio import SeqIO
import subprocess 
import csv


#patterns
seq1 = "CCCTAACCCTAA"
seq2 = "CCTAACCCTAAC"
seq3 = "CTAACCCTAACC"
seq4 = "TAACCCTAACCC"
seq5 = "AACCCTAACCCT"
seq6 = "ACCCTAACCCTA"

seq7 = "TTAGGGTTAGGG"
seq8 = "TAGGGTTAGGGT"
seq9 = "AGGGTTAGGGTT"
seq10 = "GGGTTAGGGTTA"
seq11 = "GGTTAGGGTTAG"
seq12 = "GTTAGGGTTAGG"



def get_read_positions(read):
  """Get read positions while allowing for inserted and soft-clipped bases.
  Args:
      read: A pysam.AlignedRead instance.
  Returns:
      A list of read positions equal in length to read.seq. 
      The result is identical to read.positions if the read does not contain any insertions or soft-clips. Read-positions that are insertions or soft-clips have None as the corresponding element in the returned list.
  """
  #print(read.qname)
  #print(read.cigar)
  # Check read actually has CIGAR
  if read.cigar is None or read.cigar==[]:
    # No CIGAR string so positions must be [] because there is no alignment.
    read_positions = []
  else:
    # From the SAM spec (http://samtools.github.io/hts-specs/SAMv1.pdf), "S may only have H operations between them and the ends of the CIGAR string".
    n = len(read.cigar)
    # If first CIGAR operation is H (5), check whether second is S (4).
    if read.cigar[0][0] == 5:
      if n > 1:
        if read.cigar[1][0] == 4:
          read_positions = [None] * read.cigar[1][1]
        else:
          read_positions = []
    # Check if first CIGAR operation is S (4).
    elif read.cigar[0][0] == 4:
      read_positions = [None] * read.cigar[0][1]
    # Otherwise there can't be any leftmost soft-clipping.
    else:
      read_positions = []
    # Add "internal" read-positions, which do not contain S/H operations and so can be extracted from the aligned_pairs property
    read_positions = read_positions + [y[1] for y in read.aligned_pairs if not y[0] is None]
    # If last CIGAR operation is H (5), check whether second-last is S (4).
    if read.cigar[n - 1][0] == 5:
      if n > 1:
        # If second-last positions is S (4), then need to pad but otherwise nothing to do (and also no need for "closing" else).
        if read.cigar[n - 2][0] == 4:
          read_positions = read_positions + [None] * read.cigar[n - 2][1]
    # Check if last CIGAR operation is S (4).
    elif read.cigar[n - 1][0] == 4:
      read_positions = read_positions + [None] * read.cigar[n - 1][1]
  return read_positions
  
############ a function to extract soft clipping lenght from .cigartuples function (belong to pysam library)
def get_S_clipping_len(read):
    S_len_start = 0
    S_len_end = 0
    if read.cigartuples[0][0] ==4:
        S_len_start = read.cigartuples[0][1]
    if read.cigartuples[-1][0] ==4:
        S_len_end = read.cigartuples[-1][1]
    
    return S_len_start, S_len_end 
    
  ##################### extracting start and end position of reference genome that reads match to it and 
                     ###  soft clipping length
bf = pysam.AlignmentFile("Test_alignment.bam", "rb")

cig_pos = []
for aligned_read in bf.fetch():
    # using function to get soft clipping length
    Soft_clipping = get_S_clipping_len(aligned_read)
    # a loop to extract start and end position of ref that aligned to the reads using get_read_position function that are numbers (int) and not None
    pos = [int(x) for x in get_read_positions(aligned_read) if x]
    # reference location of aligned reads
    rname = bf.getrname(aligned_read.tid)
    
    cig_pos.append((aligned_read.qname, rname, pos[0], pos[-1], Soft_clipping[0],Soft_clipping[1]))


with open('Test_output.csv', mode='w',newline='') as res:
    writer = csv.writer(res, delimiter=',')
    writer.writerow(["qname","rname","pos_start","pos_end","S_start","S_end"])
    writer.writerows(cig_pos)

