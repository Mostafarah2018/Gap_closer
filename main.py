import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord




class Reads:
    def __init__(self, bam_file_name, chromosome_file_name,reads_file_name):
        self.ref={}
        with open(chromosome_file_name) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                self.ref[record.id]=record.seq
            print(self.ref.keys())
        self.bam_file=pysam.AlignmentFile(bam_file_name, "rb")
        self.reads={}
        main_list = ["CCCTAACCCTAA","CCTAACCCTAAC","CTAACCCTAACC", "TAACCCTAACCC", "AACCCTAACCCT", "ACCCTAACCCTA"]
        
        rev_list=[ "TTAGGGTTAGGG", "TAGGGTTAGGGT", "AGGGTTAGGGTT", "GGGTTAGGGTTA", "GGTTAGGGTTAG", "GTTAGGGTTAGG"]
        self.main_list=rev_list+main_list
        with open(reads_file_name) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq=str(record.seq)
                self.reads[record.id]=seq
    def pattern_matcher(self,read_name):
        read=self.reads[read_name]
        print(read)
        for start in range(len(read)-len(self.main_list[0])):
                for rep in self.main_list:
                    if rep==read[start:start+len(rep)]:
                            return True
        return False
        
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
        elif cigar[p][0]==5 and cigar[self.get_next(p)][0]==4:
            clip=cigar[self.get_next(p)][1]
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
    def return_seq(self,name,start, end ):
        return self.read[name][start:end]
    
    def get_top_reads(self, pos, chrom):
        clips={}
        for aligned_read in self.bam_file.fetch():
            if self.reads.get(aligned_read.qname,None) is None:
                continue
            if aligned_read.reference_name!=chrom:
                continue
            if self.checker(aligned_read,pos):
                #clips[aligned_read.qname]=self.reads[aligned_read.qname]
                if self.pattern_matcher(aligned_read.qname):
                    cigar=aligned_read.cigar
                    s=sum([i[1] for i in cigar])
                    clips[aligned_read.qname]={"start_clip_len":self.get_clip_len(cigar,'start'),
                                               "start_clip_rate":self.get_clip_len(cigar,'start')/s,
                                               "end_clip_len":self.get_clip_len(cigar,'end'),
                                               "end_clip_rate":self.get_clip_len(cigar,'end')/s,
                                              "read_size":s,
                                              "alignment_size":aligned_read.reference_end -aligned_read.reference_start+1}

        return clips
    def save_seq(self,seq,out):
        records =[ SeqRecord(
                        Seq(seq[i]),
                        id=i,
                    ) for i in seq]
        SeqIO.write(records, out, "fasta")
    def ref_repair(self, chrom,pos,out_name):
        data=self.get_top_reads(pos=pos, chrom=chrom)
        if len(data)==1:
            read_name= list(data.keys())[0]
            chrm=self.ref[chrom]
            cliped_seq=self.reads[read_name][:data[read_name]["start_clip_len"]]
            
            if pos == "start":
                self.ref[chrom]=cliped_seq+chrm
            elif pos == "end":
                self.ref[chrom]=chrm+cliped_seq
            self.save_seq(self.ref,out_name)
            
            

r=Reads(bam_file_name="./Test_alignment.bam",chromosome_file_name="ref.fasta",reads_file_name="./Test_Reads.fasta")
