import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
import argparse
import json
import os

#B i 1 solved

class Log:
    def __init__(self,path="./log.txt"):
        self.file=open(path,"w")
    def write_dic(self,dic):
        self.file.write(json.dumps(dic))
        self.file.write("\n")
    def write_df(self,df):
        dfAsString = df.to_string()
        self.file.write(dfAsString)
        self.file.write("\n")
    def write_txt(self,txt):
        self.file.write(str(txt))
        self.file.write("\n")
    def top_reads_logger(self,data):
        print(len(data["marker"]))
        print(len(data["no_marker"]))
        if len(data["marker"])==0:
            self.write_txt("* Non of reads has marker")
        else:
            self.write_txt("*  marker read stat : ")
            self.write_df(pd.DataFrame(data["marker"]).T)
        
        if len(data["no_marker"])==0:
            self.write_txt("* all of the reads has marker")
        else:
            self.write_txt("* no marker read stat : ")
            self.write_df(pd.DataFrame(data["no_marker"]).T)
#B i 1 solved

class Reads:
    def __init__(self, bam_file_name, chromosome_file_name,reads_file_name, intrv_file_name=None,log_path="./log.txt"):
        self.log=Log()
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

        self.log.top_reads_logger(clips)
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
    def no_marker_repair(self,pos,reads_stat,chrmm,cut=2000,check_ratio=0.2):
        # selest maximum clipings with alignment size > X defult 2000
        target_clip_len_key=pos+"_clip_len"
        if pos == "start":
            other_side_clip_key="end_clip_len"
        elif pos == "end":
            other_side_clip_key="start_clip_len"        
        self.log.write_txt("proccessing reads with no marker")
        read = sorted([(i,(reads_stat[i][target_clip_len_key],reads_stat[i][other_side_clip_key])) for i in reads_stat if reads_stat[i]["alignment_size"]>cut], key=lambda x: x[1][0])
        check_size=int(len(reads_stat)*check_ratio)
        read=dict(read[-check_size:] )
        
        self.log.write_txt("Selected reads with %f percent max length : ".format(check_ratio))
        self.log.write_dic(read)
        #align=dict(sorted([(i,reads_stat[i]["alignment_size"]) for i in reads_stat], key=lambda x: x[1])[-check_size:]) 
        #intersection=np.intersect1d(list(read.keys()),list(align.keys()))
        check_data=read#{key:reads_stat[key][target_clip_len_key] for key in intersection}
        res={}
        if pos == "start":
            for rec in read:
                res[rec]=self.reads[rec][:read[rec][0]]
        else:
            for rec in read:
                res[rec]=self.reads[rec][-read[rec][0]:]
                print(rec,len(self.reads[rec][-read[rec][0]:]))
        self.log.write_txt("Aligning clipping regions of selected reads")
        self.save_seq(res,"clip.fasta")
        os.system("blastn -subject clip.fasta -query Read_test4.fa -outfmt 6 -out out.txt")
        self.log.write_txt("Alignment done")

        seqs=''
        df=pd.read_csv("./out.txt",sep="\t",header=None)
        
        longest_clip_read=sorted([(reads_stat[i][target_clip_len_key],i) for i in reads_stat])[-1][1]
        self.log.write_txt("Longest read selected : "+str(longest_clip_read))
        seqs+=res[longest_clip_read]
        df=df[df[1]==longest_clip_read].sort_values(3,ascending=False)
        df["pattern"]=df[0].apply(lambda x: r.pattern_matcher(x))
        df["len"]=df[0].apply(lambda x: len(r.reads[x]))
        df["loc"]=df[0].apply(lambda x: self.pattern_loc_matcher(x))
        if pos=="end":
            df["rem"]=df["loc"]-df[7]
        else:
            df["rem"]=df[6]
            
        f=df[df["pattern"]==True]
        
        self.log.write_txt("Alignemt dataFrame : ")
        self.log.write_df(f)
        
        read_len=self.read_select(list(f.rem))
        f=f[f["rem"]<=read_len]
        self.log.write_txt("Marker read with median size %d selected".format(read_len))
        self.log.write_df(df)
        final_rec=f[f["rem"]==max(f.rem)].iloc[0]
        read_name=final_rec[0]
        rem_len=final_rec["rem"]
        self.log.write_txt("chromosome Updated")
        if pos=="end":
            seqs+=self.reads[read_name][-rem_len:]
        else:
            seqs=self.reads[read_name][:rem_len]+seqs
        return seqs
    def read_select(self,lenght):
        data=sorted([(i,z) for z,i in enumerate(lenght)])
        diff=1000
        min_rec=0
        clust=[]
        cl=[data[0]]
        for z,rec in enumerate(data):
            if z==0:
                continue
            dif=rec[0]-data[z-1][0]
            if dif>diff:
                print(cl,len(cl),rec)
                if len(cl)==0:
                    print(rec)
                    cl=[rec]
                clust.append(cl)
                cl=[]

            cl.append(rec)
        clust.append(cl)
        read_len=np.median([j[0] for j in clust[np.argmax([len(i) for  i in clust])]])
        return read_len
        
        
        
    def pattern_loc_matcher(self,read_name):
                read=r.reads[read_name]
               # print(read)
                l=len(read)
                for start in range(len(read)-len(r.main_list[0])):
                    if start<100 or start>l-100:
                        for rep in r.main_list:
                            if rep==read[start:start+len(rep)]:
                                    return start#,start,rep
                return -1#
    
        
    def ref_repair(self, chrom,pos,out_name="updated_chomome.fasta"):
        self.log.write_txt("* Extracting data : ")
        chrm=self.ref[chrom]
        data=self.get_top_reads(pos=pos, chrom=chrom)
        
        
        
        if len(data["marker"])==0 and len(data["no_marker"])>0:
            
            self.log.write_txt("data has no marker reads pocessing reads with no marker")
            cliped_seq=self.no_marker_repair(pos,data["no_marker"],chrom)
                
                
        if len(data["marker"])==1:
            
            read_name= list(data["marker"].keys())[0]
            self.log.write_txt("data has reads with marker, adding read wtih "+str(read_name)+" id")
            
            cliped_seq=self.reads[read_name][:data[read_name]["start_clip_len"]["marker"]]
            
        if pos == "start":
                self.ref[chrom]=cliped_seq+chrm
        elif pos == "end":
                self.ref[chrom]=chrm+cliped_seq
        self.save_chrom("updated_chomome.fasta")
        self.log.write_txt("chromosome updated!")
        

def main():
    parser = argparse.ArgumentParser(description='A test program.')
    parser.add_argument("-r", "--read", help="Read name", type=str,required=True)
    parser.add_argument("-p", "--position", help="position of read", required=True)

    args = parser.parse_args()
    r=Reads(bam_file_name="./Test4.bam",chromosome_file_name="./Ref_test4.fa",reads_file_name="./Read_test4.fa")
    r.ref_repair(args.read,args.position)
    r.log.file.close()

if __name__ == "__main__":
    main()


