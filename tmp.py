import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
import json
from intervaltree import Interval, IntervalTree


import os

def complement_reverse(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement
#B i 1 solved

class Reads:
    def __init__(self, bam_file_name, chromosome_file_name,reads_file_name, rep_file=None ,intrv_file_name=None,
                 marker_read="CCCTAACCCTAA", marker_check_range=100):
        self.marker_check_range=marker_check_range
        self.ref={}
        with open(chromosome_file_name) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                self.ref[record.id]=record.seq
            print(self.ref.keys())
        self.logs=[]
        self.bam_file=pysam.AlignmentFile(bam_file_name, "rb")
        self.reads={}
        r_marker_read=complement_reverse(marker_read)
        self.main_list=[marker_read[r:] + marker_read[:r] for r in range(len(marker_read))]\
                        +[r_marker_read[r:] + r_marker_read[:r] for r in range(len(r_marker_read))]
        self.main_list=list(set(self.main_list))

        for z,i in enumerate(self.main_list):
            self.main_list[z]=i+i[:6]

        with open(reads_file_name) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq=str(record.seq)
                self.reads[record.id]=seq
        
        if intrv_file_name is not None:
            self.set_intervals(intrv_file_name)
        else:
            self.intervals=None
        if rep_file is not None:
            self.reps=records = list(SeqIO.parse("./Repeats.fa",format="fasta"))
        else:
            self.reps=None
            
    def set_intervals(self,f_name):
        df=pd.read_csv(f_name,header=None)
        df=df[[1,2,3]]
        df.columns=["chromosome","start","end"]
        intervals={chrm:IntervalTree() for chrm in set(df["chromosome"])}
        for z,rec in df.iterrows():
            intervals[rec["chromosome"]][rec["start"]:rec["end"]]=["repeat"]
        for chrm in intervals:
            intervals[chrm].merge_overlaps()
        self.intervals=intervals

    def pattern_matcher(self,read_name):
        read=self.reads[read_name]
        ranges=[range(self.marker_check_range) , range(len(read)-len(self.main_list[0])-self.marker_check_range)]
        for rang in ranges:
            for start in rang:
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
    
 #   def repeat_region(self, chrom, pos, gap=500):
  #      self.intervals[chrom].
    def pattern_select(self,df,pos,seqs=""):
            df["len"]=df[0].apply(lambda x: len(self.reads[x]))
            df["loc"]=df[0].apply(lambda x: self.pattern_loc_matcher(x))
            if pos=="end":
                df["rem"]=df["loc"]-df[7]
            else:
                df["rem"]=df[6]

            f=df[df["pattern"]==True]

            self.log.write_txt("Alignemt dataFrame : ")
            self.log.write_df(f)
            f["align_len"]=f[7]-f[6]

            final_rec=f[f["align_len"]==max(f["align_len"])].iloc[0]
            read_name=final_rec[0]
            print("Marker read with id "+str(read_name)+" selected")
            self.log.write_txt("Marker read with id "+str(read_name)+" selected")

            rem_len=final_rec["rem"]
            self.log.write_txt("chromosome Updated")
            if pos=="end":
                seqs+=self.reads[read_name][-rem_len:]
            else:
                seqs=self.reads[read_name][:rem_len]+seqs
            return seqs
    
    def checker(self,read ,pos, inter ,thr_seq=500, thr_tail=50):
        
        pos_read=[i for i in read.positions if len(inter[i])>0]
        cover=len(pos_read)
        if cover==0:
            return False
        diff=0
        if pos=="start":
            read_start=min(pos_read)
            for interval in inter:
                if interval.begin<read_start:
                    diff+=min(interval.end,read_start)-interval.begin
        if pos=="end":
            read_end=max(pos_read)
            for interval in inter:
                if interval.end>read_end:
                    diff=interval.end-max(interval.begin,read_end)
        
       # if pos.lower()=="end":
            
      #      length=len(self.ref[read.reference_name])-thr_tail
      #      diff=length-read.reference_end
      #  else:
       #     ls=[i for i in read.positions if i>=start]
        #    if len(ls)>0:
         #       length=min(ls)
          #      diff=length-start
          #  else:
           #     return False
        
        if diff<=thr_seq and self.get_clip_len(read.cigar, pos)>thr_tail:
            return True
        else:
            return False
        
    def save_seq(self,seq,out):
        records =[ SeqRecord(
                        Seq(seq[i]),
                        id=i,
                    ) for i in seq]
        SeqIO.write(records, out, "fasta")
    def interval_region_extractor(self,chrm,pos,start,end,thr_seq=500):
            res=IntervalTree()
            reps=self.intervals[chrm][start:end]
            if len(reps)==0:
                res[start:end]=["seq"]
            else:
                sum_inter=0
                while sum_inter<thr_seq:

                    if pos.lower()=="start":
                        rep=sorted([(rep.begin,rep) for rep in reps])[0][1]
                        lenght=min(rep.begin-start,thr_seq-sum_inter)
                        res[start:start+lenght]=["seq"]
                        sum_inter+=lenght
                        if sum_inter>=thr_seq:
                            break
                        end=rep.end+end-(start+lenght)
                        start=rep.end
                        reps=self.intervals[chrm][start:end]
                        print(reps,start,end)
                        if len(reps)==0:
                            print(start,end)
                            res[start:end]=["seq"]
                            break

                    elif pos.lower()=="end":
                        rep=sorted([(rep.end,rep) for rep in reps])[-1][1]
                        print("*",rep)
                        lenght=min(end-rep.end,thr_seq-sum_inter)
                        res[end-lenght:end]=["seq"]

                        sum_inter+=lenght
                        if sum_inter>=thr_seq:
                            break
                        end=rep.begin#+end-(start+lenght)
                        start=rep.begin-(thr_seq-sum_inter)
                        reps=self.intervals[chrm][start:end]
                        if len(reps)==0:
                            print(start,end)
                            res[start:end]=["seq"]
                            break
            return res
    def region_finder(self,chrm,pos,start,end,thr_seq=500,similarity=95):
        res=IntervalTree()
        reps=IntervalTree()
        lenght=0
        while(lenght<thr_seq):
            if thr_seq-lenght<500:
                end=start+500
            else:
                end=start+thr_seq-lenght
            
            seq=self.ref[chrm][start:end]
            self.save_seq({"seq":seq},"rf_ref_search.fasta")
            os.system("blastn -subject out.fasta -query Repeats.fa -outfmt 6 -out rf_ref_search.txt")
            df=pd.read_csv("rf_ref_search.txt",sep="\t",header=None)
            df=df[df[2]>=similarity]
            if len(df[2])==0:
                res[start:end]=["seq"]
                lenght+=end-start

            for row,data in df.iterrows():
                inter_start=min(data[8],data[9])
                inter_end=max(data[8],data[9])
                reps[inter_start:inter_end]=["rep"]
            reps.merge_overlaps()
            rep_pos=[(inter.begin,inter.end) for inter in reps]
            rep_pos.sort()
            tree_start=0
            for inter in rep_pos:
                end_tree=min(inter[0]-tree_start,thr_seq-lenght)
                if lenght==thr_seq:
                    break
                res[tree_start:tree_start+end_tree]=["seq"]
                lenght+=end_tree
                tree_start=inter[1]
            start=end
        
        return res
    
    def region_extractor(self,chrm,pos,thr_seq=500,similarity=95):
        
        if pos.lower()=="start":
            start=0
            end=thr_seq
        elif pos.lower()=="end":
                end=len(self.ref[chrm])
                start=end-thr_seq
        
        if self.intervals is not None:
            return self.interval_region_extractor(chrm,pos,start,end,thr_seq)
        if self.intervals is None:
            return self.region_finder(chrm,pos,start,end, thr_seq,similarity)
    
    def  check_merge_start_Intervals(self,chrm,pos,diff,similarity):
        return self.region_extractor(chrm,pos,thr_seq=diff,similarity=similarity)
    
    def get_top_reads(self, pos, chrom,diff=500,similarity=500):
        clips={"marker":{},"no_marker":{}}
        
        inter=self.check_merge_start_Intervals(chrom,pos,diff,similarity)
        
        
        
        
        for aligned_read in self.bam_file.fetch():
            if self.reads.get(aligned_read.qname,None) is None:
                continue
            if aligned_read.reference_name!=chrom:
                continue
            if self.checker(aligned_read,pos,inter):
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
        self.save_seq(self.ref,out)
    
    
    def no_marker_repair(self,pos,reads_stat,chrmm,cut=2000,check_ratio=0.2,max_len=True,align_cnt=5,n_iter=3):
        # selest maximum clipings with alignment size > X defult 2000
        target_clip_len_key=pos+"_clip_len"
        if pos == "start":
            other_side_clip_key="end_clip_len"
        elif pos == "end":
            other_side_clip_key="start_clip_len"        
        read = sorted([(i,(reads_stat[i][target_clip_len_key],reads_stat[i][other_side_clip_key])) for i in reads_stat if reads_stat[i]["alignment_size"]>cut], key=lambda x: x[1][0])
        check_size=int(len(reads_stat)*check_ratio)
        read=dict(read[-check_size:] )
        
        #self.log.write_txt("Selected reads with %f percent max length : ".format(check_ratio))
        #self.log.write_dic(read)
        check_data=read#{key:reads_stat[key][target_clip_len_key] for key in intersection}
        res={}
        if pos == "start":
            for rec in read:
                res[rec]=self.reads[rec][:read[rec][0]]
        else:
            for rec in read:
                res[rec]=self.reads[rec][-read[rec][0]:]
                print(rec,len(self.reads[rec][-read[rec][0]:]))
        #self.log.write_txt("Aligning clipping regions of selected reads")
        if max_len is True:
            max_len_id=np.argmax([len(res[i]) for i in res])
            key=list(res.keys())[max_len_id]
            res={key:res[key]}
        self.save_seq(res,"clip.fasta")
        print("Runing blast on "+key+" read")
        os.system("blastn -subject clip.fasta -query Read_test4.fa -outfmt 6 -out outI.txt")
        #self.log.write_txt("Alignment done")

        seqs=''
        df=pd.read_csv("./outI.txt",sep="\t",header=None)
        
        longest_clip_read=sorted([(reads_stat[i][target_clip_len_key],i) for i in reads_stat])[-1][1]
        print(sorted([(reads_stat[i][target_clip_len_key],i) for i in reads_stat]))
        #self.log.write_txt("Longest read selected : "+str(longest_clip_read))
        self.logs.append(res[longest_clip_read])
        seqs+=res[longest_clip_read]
        
        print(1)
        print(df)
        df=df[df[1]==longest_clip_read]
        print(2)
        print(df)
        df=df[df[0]!=longest_clip_read].sort_values(3,ascending=False)
        print(longest_clip_read)
        print(df)
        df["pattern"]=df[0].apply(lambda x: self.pattern_matcher(x))
        if True in df["pattern"]:
            print("read ")
            print(df[df["pattern"]==True][0])
            seqs=self.pattern_select(df,pos,seqs)
        else:
            select_reads=[longest_clip_read]
            for align_iter in range(n_iter):
                df["align_len"]=df[7]-df[6]
                df=df[df["align_len"]>1000]
                if pos=="end":
                    df["len"]=df[0].apply(lambda x:len(self.reads[x]))
                    df["rem"]=df["len"]-df[7]
                    start=df[7]
                    end=df["len"]
                elif pos=="start":
                    df["rem"]=df[6]
                    start=0
                    end=df[6]
                df=df[df.rem==max(df.rem)]
                max_len_read=df.iloc[0][0]
                print("unmarked read "+max_len_read+" with alignment "+str(df.iloc[0]["align_len"])+" clip size "++str(df.iloc[0]["rem"])+" selected")
                select_reads.append(max_len_read)
                self.logs.append(len(seq))
                if pos=="start":
                    seq=self.reads[max_len_read][start:end+1]+seqs
                elif pos=="end":
                    seqs=seq+self.reads[max_len_read][start:end+1]
                
                self.save_seq({max_len_read:self.reads[max_len_read][start:end+1]},"clip_"+str(align_iter)+".fasta")
                print(" At step "+str(align_iter+1)+" Runining Blast")
                os.system("blastn -subject clip_"+str(align_iter)+".fasta -query Read_test4.fa -outfmt 6 -out out_"+str(align_iter)+".txt")
                df=pd.read_csv("./out_"+str(align_iter)+".txt",sep="\t",header=None)
                df=df[[i not in select_reads for i in df[0]]]
                df["pattern"]=df[0].apply(lambda x: self.pattern_matcher(x))
                if True in df["pattern"]:
                    seqs=self.pattern_select(df,pos,seqs)
                    break
                else:
                    print("no marker Found ")
                if align_iter==n_iter-1:
                    raise Exception("No marker  read found after "+str(n_iter)+" iterrations")
                
                    
      #      for align_iter in range(align_cnt):
                
        return seqs
    
    
    def ref_repair(self, chrom,pos,out_name="updated_chomome.fasta"):
        chrm=self.ref[chrom]
        data=self.get_top_reads(pos=pos, chrom=chrom)
                
        if len(data["marker"])==0 and len(data["no_marker"])>0:
            print("No Read Found with marker")
            cliped_seq=self.no_marker_repair(pos,data["no_marker"],chrom)
                
        if len(data["marker"])>0:
            #Issue : Which Record shoulb be selected
            #2000 match and biggest
            read_name= list(data["marker"].keys())[0]
            print("read "+read_name+" found with marker")
            self.log,append(len(cliped_seq))
            cliped_seq=self.reads[read_name][:data["marker"][read_name]["start_clip_len"]]
        if pos == "start":
                    
                self.ref[chrom]=cliped_seq+chrm
        elif pos == "end":
                self.ref[chrom]=chrm+cliped_seq
        print("Process Finished")
        self.save_chrom("updated_chomome.fasta")
    

r=Reads(bam_file_name="./Test5.bam",chromosome_file_name="./Test5_Ref.fasta",reads_file_name="./Read_test4.fa",
       intrv_file_name="./Test5_Repeat_Location.csv")
r.ref_repair( "Chr4","start",out_name="updated_chomome.fasta")
