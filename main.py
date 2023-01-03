import pysam
from Bio import SeqIO

class Reads:
    def __init__(self, bam_file_name, chromosome_file_name,reads_file_name):
        self.ref={}
        with open(chromosome_file_name) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                self.ref[record.id]=record.seq
            print(self.ref.keys())
        self.bam_file=pysam.AlignmentFile(bam_file_name, "rb")
        self.reads={}
       # with open(reads_file_name) as handle:
        #    for record in SeqIO.parse(handle, "fasta"):
         #       self.reads[record.id]=record.seq
    
        
    def get_next(self,i):
        if i == 0:
            return 1
        elif i==-1:
            return -2
        else:
            return None


    def get_clip_len(self,cigar, pos):
        if pos.lower()=="start":
            p=0
        elif pos.lower()=="end":
            p=-1
        if cigar[p][0]==4 :
            clip=cigar[p][1]
        elif cigar[p][0]==5 and cigar[get_next(p)][0]==4:
            clip=cigar[get_next(p)][1]
        else:
            clip=0
        return clip


    def checker(self,read ,pos,thr=50):
        if pos.lower()=="end":
            length=len(self.ref[read.reference_name])-thr
        else:
            length=0
        diff=read.reference_start-length
        if diff<=thr and self.get_clip_len(read.cigar, pos)>diff:
            return True
        else:
            return False
    def get_top_reads(self, pos, chrom):
        clips={}
        for aligned_read in bf.fetch():
            if aligned_read.reference_name!=chrom:
                continue
            if self.checker(aligned_read,pos):        
                cigar=aligned_read.cigar
                s=sum([i[1] for i in cigar])
                clips[aligned_read.qname]={"start_clip_len":self.get_clip_len(cigar,'start'),
                                           "start_clip_rate":self.get_clip_len(cigar,'start')/s,
                                           "end_clip_len":self.get_clip_len(cigar,'end'),
                                           "end_clip_rate":self.get_clip_len(cigar,'end')/s}
                
        return clips
