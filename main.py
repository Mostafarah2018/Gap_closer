import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os

#B i 1 solved

class Reads:
    def __init__(self, bam_file_name, chromosome_file_name,reads_file_name, intrv_file_name=None):
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
        for z,i in enumerate(self.main_list):
            self.main_list[z]=i+i[:6]
        print(self.main_list)
        with open(reads_file_name) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq=str(record.seq)
                self.reads[record.id]=seq
        if intrv_file_name is not None:
            self.set_intervals(intrv_file_name)
        else:
            self.intervals=None
            
    def set_intervals(self,f_neme):
        
        df=pd.read_csv(f_name,header=None)
        df.columns=["Unk","chromosome","start","end"]
        df=df[df.chromosome==chrm]
        self.intervals=[[i["start"],i["end"]] for  _,i in df.iterrows()]
        self.intervals.sort()
        
    def pattern_matcher(self,read_name):
        read=self.reads[read_name]
       # print(read)
        l=len(read)
        for start in range(len(read)-len(self.main_list[0])):
            if start<100 or start>l-100:
                for rep in self.main_list:
                    if rep==read[start:start+len(rep)]:
                            return True#,start,rep
        return False#,None,None
    
    def seq_pattern_extractor(self):
        seqs={}
        for read in self.reads:
            a=self.pattern_matcher(read)
            if a:
                #print(b,c,read)
                seqs[read]=self.reads[read]
        return seqs
        
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


    def checker(self,read ,pos,start,thr_seq=500, thr_tail=50):
        if pos.lower()=="end":
            length=len(self.ref[read.reference_name])-thr_tail
            diff=length-read.reference_end
           # print(length,diff,read.reference_end)
        else:
            ls=[i for i in read.positions if i>=start]
            if len(ls)>0:
                length=min(ls)
                diff=length-start
            else:
                return False
        
        if diff<=thr_seq and self.get_clip_len(read.cigar, pos)>thr_tail:
            return True
        else:
            return False
    def return_seq(self,name,start, end ):
        return self.read[name][start:end]
    
    def  check_merge_start_Intervals(self,chrm,diff=500):
        # Sort the array on the basis of start values of intervals.
        if self.intervals is None:
            return None
        if self.intervals[0][0]>diff:
            return None
        
        stack = []
        stack=self.intervals[0]
        for i in self.intervals[1:]:
            if stack[0] <= i[0] and i[0] <= stack[-1]+diff:
                stack[-1] = max(stack[-1], i[-1])
            else:
                break

        return stack
    
    def get_top_reads(self, pos, chrom):
        clips={"marker":{},"no_marker":{}}
        
        inter=self.check_merge_start_Intervals(chrom,)
        if inter is None:
            if pos=="start":
                start=0
            else:
                start= len(self.ref[chrom])
        else:
            start=inter[-1]
        for aligned_read in self.bam_file.fetch():
            if self.reads.get(aligned_read.qname,None) is None:
                continue
            if aligned_read.reference_name!=chrom:
                continue
            if self.checker(aligned_read,pos,start):
                #clips[aligned_read.qname]=self.reads[aligned_read.qname]
                if self.pattern_matcher(aligned_read.qname):
                    key="marker"
                else:
                    key="no_marker"
                cigar=aligned_read.cigar
                s=sum([i[1] for i in cigar])
                clips[key][aligned_read.qname]={"start_clip_len":self.get_clip_len(cigar,'start'),
                                               "start_clip_rate":self.get_clip_len(cigar,'start')/s,
                                               "end_clip_len":self.get_clip_len(cigar,'end'),
                                               "end_clip_rate":self.get_clip_len(cigar,'end')/s,
                                              "read_size":s,
                                              "alignment_size":aligned_read.reference_end -aligned_read.reference_start+1}

        return clips
    def save_chrom(self, out="updated_chromosome.fasta"):
        self.save_seq(self.reads,out)
    def save_seq(self,seq,out):
        records =[ SeqRecord(
                        Seq(seq[i]),
                        id=i,
                    ) for i in seq]
        SeqIO.write(records, out, "fasta")
    def chrom_repair_unibase(self, read_name, chrom, clip_start,clip_end):
        chrm=self.ref[chrom]
        cliped_seq=self.reads[read_name][s]
    def no_marker_repair(self,pos,reads_stat,chrmm,cut=2000):#check_ratio=0.2):
        # selesct maximum clipings with alignment size > X defult 2000
        target_clip_len_key=pos+"_clip_len"
        if pos == "start":
            other_side_clip_key="end_clip_len"
        elif pos == "end":
            other_side_clip_key="start_clip_len"
        print(target_clip_len_key,other_side_clip_key)
        
       # check_size=int(len(reads_stat)*check_ratio)
        read=dict(sorted([(i,(reads_stat[i][target_clip_len_key],reads_stat[i][other_side_clip_key])) for i in reads_stat if reads_stat[i]["alignment_size"]>cut], key=lambda x: x[1][0]))#[-check_size:] )
        print(read)
        #align=dict(sorted([(i,reads_stat[i]["alignment_size"]) for i in reads_stat], key=lambda x: x[1])[-check_size:]) 
        #intersection=np.intersect1d(list(read.keys()),list(align.keys()))
        check_data=read#{key:reads_stat[key][target_clip_len_key] for key in intersection}
        res={}
        if pos == "start":
            for rec in read:
                res[rec]=self.reads[rec][:read[rec][0]]
        else:
            for rec in read:
                res[rec]=self.reads[rec][read[rec][0]:]
        self.save_seq(res,"clip.fasta")
        os.system("blastn -subject clip.fasta -query Read_test4.fa -outfmt 6 -out out.txt")
        
        
    def ref_repair(self, chrom,pos,out_name):
        data=self.get_top_reads(pos=pos, chrom=chrom)
        if len(data["marker"])==0 and len(data["no_marker"])>0:
            if len(data["no_marker"])==1:
                cliped_seq=self.reads[read_name][:data[read_name]["start_clip_len"]["marker"]]
                
                
                
        if len(data["marker"])==1:
            read_name= list(data.keys())[0]
            chrm=self.ref[chrom]
            cliped_seq=self.reads[read_name][:data[read_name]["start_clip_len"]["marker"]]
            
            if pos == "start":
                self.ref[chrom]=cliped_seq+chrm
            elif pos == "end":
                self.ref[chrom]=chrm+cliped_seq
            #self.save_seq(self.ref,out_name)
        
        elif len(data)>1:
            rs={}
            for marker in data:
                clips={}
                for read in data[marker]:
                    if pos == "start":
                        clips[read]=self.reads[read]#[:data[read]["start_clip_len"]]
                    elif pos == "end":
                        clips[read]=self.reads[read]#[-data[read]["end_clip_len"]:]
                rs[marker]=clips
            return rs
            
                
                            
                
              #  print(str(aligned_read))
             #   if abs( aligned_read.reference_start- aligned_read.reference_end)+1/s  >0.50:
r=Reads(bam_file_name="./Test4.bam",chromosome_file_name="./Ref_test4.fa",reads_file_name="./Read_test4.fa")
data=r.get_top_reads("end",'tig00018')
r.no_marker_repair("end",data["no_marker"],"tig00018",cut=2000)
