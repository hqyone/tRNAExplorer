#  tRNAExplorer v.1.0 - tRNA profiles and regulation networks
#  Copyright (C) 2020  Dr. Quanyuan He
#  School of Medicine, Hunan Normal University
#  Email: hqyone@hotmail.com
#  Freely distributed under the GNU General Public License (GPLv3)
#

import lib_code.share as share
import re

class tRNA():
    def __init__(self):
        self.name = ""
        self.chrom=""
        self.start=""
        self.end=""
        self.strand=""
        self.full_seq=""   # May contain introns
        self.seq=""        # without_intron
        self.intron_infor = ""
        self.anticodon = ""
        self.anticodon_start = 0
        self.anticodon_end = 0
        self.acceptor = ""
        self.family = ""
        # Fasta file
        self.seq_utr = ""
        self.utr_len = 0
        # tRNAScan SE results
        self.map_start = ""
        self.map_end = ""
        self.map_len = 0
        self.map_structure_str = ""
        self.map_seq = ""
        self.anticodon_start_in_map = 0
        self.anticodon_end_in_map = 0
        self.map_scan_score = 0
        self.possible_type = ""
        self.d_loop={'start':-1,'l_start':-1,'l_end':-1,'end':-1, 'struct_str':"",'for_str':"",'rev_str':"", 'loop_str':""}
        self.a_loop={'start':-1,'l_start':-1,'l_end':-1,'end':-1, 'struct_str':"",'for_str':"",'rev_str':"", 'loop_str':""}
        self.v_loop={'start':-1,'l_start':-1,'l_end':-1,'end':-1, 'struct_str':"",'for_str':"",'rev_str':"", 'loop_str':""}
        self.t_loop={'start':-1,'l_start':-1,'l_end':-1,'end':-1, 'struct_str':"",'for_str':"",'rev_str':"", 'loop_str':""}
        self.stem_for = {'start':-1,'end':-1, 'str':""}
        self.stem_rev = {'start':-1, 'end':-1, 'str': ""}

    def GetTabTitle(self):
        return "#chrom\tstart\tend\tname\tscore\tstrand\tfull_seq\tseq\tintron_infor\tanticodon\tanticodon_start"+ \
               "\tanticodon_end\tacceptor"+ \
               "\tfamily\tseq_utr\tutr_len\tmap_start\tmap_end\tmap_len\tmap_structure_str\tmap_seq"+\
               "\tanticodon_start_in_map\tanticodon_end_in_map\tmap_scan_score\tpossible_type"+\
               "\td_loop_start\td_loop_lstart\td_loop_lend\td_loop_end" + \
               "\ta_loop_start\ta_loop_lstart\ta_loop_lend\ta_loop_end" + \
               "\tv_loop_start\tv_loop_lstart\tv_loop_lend\tv_loop_end" + \
               "\tt_loop_start\tt_loop_lstart\tt_loop_lend\tt_loop_end" + \
               "\tstem_for_start\tstem_for_end\tstem_rev_start\tstem_rev_end"

    # Print a bed like file
    def GetTabStr(self):
        a = str(self.chrom)+"\t"+str(self.start)+"\t"+str(self.end)+"\t"+self.name+"\t0\t"+str(self.strand)+"\t"+str(self.full_seq)+"\t"+str(self.seq)+"\t"+ \
            self.intron_infor + "\t" + self.anticodon+"\t"+str(self.anticodon_start) +"\t"+str(self.anticodon_end) + "\t" + \
            str(self.acceptor)+"\t"+str(self.family) +"\t"+ str(self.seq_utr) + "\t" +str(self.utr_len) + "\t" +\
            str(self.map_start)+"\t"+str(self.map_end)+ "\t"+str(self.map_len) +"\t"+self.map_structure_str +"\t"+self.map_seq +"\t"+ \
            str(self.anticodon_start_in_map)+"\t"+str(self.anticodon_end_in_map)+"\t"+str(self.map_scan_score)+"\t"+str(self.possible_type)+"\t"+ \
            str(self.d_loop["start"])+"\t"+str(self.d_loop["l_start"])+"\t"+str(self.d_loop["l_end"])+"\t"+str(self.d_loop["end"])+"\t"+ \
            str(self.a_loop["start"]) + "\t" + str(self.a_loop["l_start"]) + "\t" + str(self.a_loop["l_end"]) + "\t" + str(self.a_loop["end"]) + "\t" + \
            str(self.v_loop["start"]) + "\t" + str(self.v_loop["l_start"]) + "\t" + str(self.v_loop["l_end"]) + "\t" + str(self.v_loop["end"]) + "\t" + \
            str(self.t_loop["start"]) + "\t" + str(self.t_loop["l_start"]) + "\t" + str(self.t_loop["l_end"]) + "\t" + str(self.t_loop["end"]) + "\t" + \
            str(self.stem_for["start"]) + "\t" + str(self.stem_for["end"])+ "\t"  + \
            str(self.stem_rev["start"]) + "\t" + str(self.stem_rev["end"])
        return a

    def LoadStr(self, str):
        contents = str.strip().split("\t")
        if len(contents)>42:
            self.chrom = contents[0]
            self.start = int(contents[1])
            self.end = int(contents[2])
            self.name = contents[3]
            self.strand = contents[5]
            self.full_seq = contents[6]
            self.seq = contents[7]
            self.intron_infor = contents[8]
            self.anticodon = contents[9]
            self.anticodon_start =int(contents[10])
            self.anticodon_end = int(contents[11])
            self.acceptor = contents[12]
            self.family = contents[13]
            # Fasta file
            self.seq_utr = contents[14]
            self.utr_len = int(contents[15])
            # tRNAScan SE results
            self.map_start = int(contents[16])
            self.map_end = int(contents[17])
            self.map_len = int(contents[18])
            self.map_structure_str = contents[19]
            self.map_seq = contents[20]
            self.anticodon_start_in_map = int(contents[21])
            self.anticodon_end_in_map = int(contents[22])
            self.map_scan_score = float(contents[23])
            self.possible_type = contents[24]
            self.d_loop = {'start': int(contents[25]), 'l_start': int(contents[26]), 'l_end': int(contents[27]),
                           'end': int(contents[28]), 'struct_str': "", 'for_str': "",
                           'rev_str': "", 'loop_str': ""}
            self.a_loop = {'start': int(contents[29]), 'l_start':int(contents[30]), 'l_end': int(contents[31]), 'end': int(contents[32]), 'struct_str': "", 'for_str': "",
                           'rev_str': "", 'loop_str': ""}
            self.v_loop = {'start': int(contents[33]), 'l_start': int(contents[34]), 'l_end': int(contents[35]),
                           'end': int(contents[36]), 'struct_str': "", 'for_str': "",
                           'rev_str': "", 'loop_str': ""}
            self.t_loop = {'start': int(contents[37]), 'l_start': int(contents[38]), 'l_end': int(contents[39]),
                           'end': int(contents[40]), 'struct_str': "", 'for_str': "",
                           'rev_str': "", 'loop_str': ""}
            self.stem_for = {'start': int(contents[41]), 'end': int(contents[42]), 'str': ""}
            self.stem_rev = {'start': int(contents[43]), 'end': int(contents[44]), 'str': ""}

    def GetMatureSeq(self):
        if self.intron_infor!="":
            s = re.search(r'intron:\s+(\d+)-(\d+)\s+\((\d+)-(\d+)\)', self.intron_infor)
            s_start = 0
            s_end = 0
            if s:
                s_start = int(s.group(1))
                s_end = int(s.group(2))
            return self.full_seq[0:s_start-1]+self.full_seq[s_end:]
        else:
            return self.full_seq

    def GetPreMatureNoIntronSeq(self):
        if self.intron_infor!="":
            s = re.search(r'intron:\s+(\d+)-(\d+)\s+\((\d+)-(\d+)\)', self.intron_infor)
            s_start = 0
            s_end = 0
            if s:
                s_start = int(s.group(3))
                s_end = int(s.group(4))
            return self.seq_utr[0:s_start-1]+self.seq_utr[s_end:]
        else:
            return self.seq_utr

    def GetKeySitesInfor(self, type):
        if type=="I":
            return [int(self.map_start), int(self.map_start+self.anticodon_start_in_map),int(self.map_start+self.anticodon_end_in_map),int(self.map_end)]
        elif type =="P":
            return [int(self.map_start), int(self.map_start+self.anticodon_start_in_map),int(self.map_start+self.anticodon_end_in_map),int(self.map_end)]
        elif type =="M":
            return [1, int(self.anticodon_start_in_map),int(self.anticodon_end_in_map),int(self.map_end)-int(self.map_start)+1]
        elif type =="C":
            return [1, int(self.anticodon_start_in_map),int(self.anticodon_end_in_map),int(self.map_end)-int(self.map_start)+4]  #Because CCA
        else:
            return []

    def IsQualifiedtRNA(self, no_mit_tRNA=True, no_pseudogenes=True, min_qscore=30):
        if no_mit_tRNA and self.name.startswith("nm"):
            return False
        if no_pseudogenes and self.possible_type == "pseudogene":
            return False
        if self.map_scan_score<min_qscore:
            return False
        return True

    # Calculate the position of reads in the tRNA with intron and UTR regions
    # read_class : four styles of tRNA:
    #   I (tRNA with UTRs and Intron) ,
    #   P (tRNA with UTRs),
    #   M (mature),
    #   C (meature tRNA with CCA)
    # start: start position of read in object position
    # end: end position of read in object position
    # I,1,75
    def CalculateAlignmentLocInI(self, read_class, start, end):
        offset = self.utr_len
        if read_class!="A" or read_class!="B":
            start+=offset
            end+=offset
        # Remove the effect of CCA
        if read_class=="I":
            end-=3
        if self.intron_infor!="":
            s = re.search(r'intron:\s+(\d+)-(\d+)\s+\((\d+)-(\d+)\)', self.intron_infor)
            s_start = 0
            s_end = 0
            if s:
                s_start = int(s.group(3))
                s_end = int(s.group(4))
                if self.intron_infor!="" and (read_class=="E" or read_class=="G" or read_class=="H" or read_class=="I"):
                    if start>s_start:
                        start+=s_end-s_start+1
                    if end>s_start:
                        end+=s_end-s_start+1
            return [start, end]
        else:
            return [start, end]

    # The position should be index of URL_Seq_INTRON
    # pos is 1 basded
    def FindPosType(self, pos):
        ptype="Other"
        if pos<-1:  # Upsteam of cleavage
            ptype="5-UTR"
        elif pos == -1:
            ptype = "Start"
        elif pos>len(self.full_seq)-1:
            ptype = "3-UTR"
        elif pos == len(self.full_seq)-1:
            ptype = "End"
        elif pos >=self.anticodon_start-self.utr_len and pos <=self.anticodon_end-self.utr_len:
            ptype = "Anticodon"
        elif self.d_loop["start"]!=-1 and  pos >= self.d_loop["start"] and pos <= self.d_loop["end"]:
            ptype = "D-loop"
        elif self.a_loop["start"]!=-1 and  pos >=  self.a_loop["start"] and pos <= self.a_loop["end"]:
            ptype = "A-loop"
        elif self.t_loop["start"]!=-1 and  pos >=  self.t_loop["start"] and pos <= self.t_loop["end"]:
            ptype = "T-loop"
        elif self.stem_for["start"]!=-1 and  pos >=  self.stem_for["start"] and pos <= self.stem_for["end"]:
            ptype = "stem_for"
        elif self.stem_rev["start"]!=-1 and  pos >=  self.stem_rev["start"] and pos <= self.stem_rev["end"]:
            ptype = "stem_rev"
        return ptype

    def FindCleavageSites(self, start_pos_profile, end_pos_profile, min_intensity=100, min_sn_ratio=10,context_seq_len=8):
        cleavageSite_Dic={}
        combined_profile = start_pos_profile[:]
        if len(start_pos_profile)==len(end_pos_profile):
            for i in range(len(end_pos_profile)):
                combined_profile[i]+=end_pos_profile[i]
        newList = sorted(combined_profile)
        median = newList[int(len(newList)/2)]
        start_peak_dic = self.DetectPeak(start_pos_profile, "startpos", context_seq_len)
        end_peak_dic = self.DetectPeak(end_pos_profile, "endpos", context_seq_len)
        # Find the peak list by combined start and end profiles
        start_peak_ls = list(start_peak_dic.keys())
        end_peak_ls = list(end_peak_dic.keys())
        combined_peak_id_ls = start_peak_ls + list(set(start_peak_ls) - set(end_peak_ls))
        for id in combined_peak_id_ls:
            combined_peak={
                "pos": 0,
                "int": 0,
                "int5": 0,
                "int3": 0,
                "sn_ratio":0,
                "c_5_seq": "",
                "c_3_seq": ""
            }
            if id in start_peak_dic:
                start_peak = start_peak_dic[id]
                combined_peak["pos"] = start_peak["pos"]
                combined_peak["int"]+=start_peak["int"]
                combined_peak["ptype"] = start_peak["ptype"]
                combined_peak["int3"] = start_peak["int"]
                combined_peak["c_5_seq"] = start_peak["c_5_seq"]
                combined_peak["c_3_seq"] = start_peak["c_3_seq"]
            if id in end_peak_dic:
                end_peak = end_peak_dic[id]
                combined_peak["pos"] = end_peak["pos"]
                combined_peak["int"]+=end_peak["int"]
                combined_peak["ptype"] = end_peak["ptype"]
                combined_peak["int5"] = end_peak["int"]
                combined_peak["c_5_seq"] = end_peak["c_5_seq"]
                combined_peak["c_3_seq"] = end_peak["c_3_seq"]
            if median==0:
                median=1
            if combined_peak["int"]>min_intensity and (float(combined_peak["int"])/median)>min_sn_ratio:
                combined_peak['sn_ratio'] = int(float(combined_peak["int"])/median)
                cleavageSite_Dic[id]=combined_peak
        return cleavageSite_Dic

    def DetectPeak(self, profile, type, context_seq_len=5):
        potential_peak={}
        for i in range(1, len(profile)-1):
            if profile[i-1]<profile[i] and profile[i]>profile[i+1]:
                pos = i
                # pos is 0 based on
                context_5_seq = ""
                context_3_seq =""
                if type == "startpos":
                    context_5_seq= self.seq_utr[max(0,pos-context_seq_len):min(len(self.seq_utr),pos)]
                    context_3_seq = self.seq_utr[max(0,pos): min(len(self.seq_utr),pos+context_seq_len)]
                elif type == "endpos":
                    context_5_seq = self.seq_utr[max(0, pos - context_seq_len+1):min(len(self.seq_utr), pos+1)]
                    context_3_seq = self.seq_utr[max(0, pos+1): min(len(self.seq_utr), pos + context_seq_len+1)]
                if type == "startpos":
                    pos-=1
                potential_peak["P"+str(pos)]={
                    "pos":pos-self.utr_len,
                    "int":profile[i],
                    "ptype":self.FindPosType(pos-self.utr_len),
                    #"sn_ratio":round(float(profile[i])/median,3),
                    "c_5_seq" :context_5_seq,
                    "c_3_seq": context_3_seq
                }
        return potential_peak

    #I,1,75
    def CreateTrfProfiles(self, brief_mapping_infor_ls, offset):
        profile_length = len(self.seq_utr)
        profiles = {
            "total_reads_num": 0,
            "pileup_height": 0,
            "start_pos": [0] * profile_length,
            "end_pos": [0] * profile_length,
            "total": [0] * profile_length,
            "cleavage_site_dic":{}
        }
        # Both "start_pos" and "end_pos" are 1 based
        for bi in brief_mapping_infor_ls:
            i=bi.split(",")
            c = i[0]  #I
            trna_start = int(i[1]) #1
            trna_end = int(i[2])   #75
            read_num = float(i[3])
            loc = self.CalculateAlignmentLocInI(c, trna_start, trna_end)
            profiles["total"]= share.plusNumList(profiles["total"], loc[0] - 1, loc[1] - 1, read_num)
            loc = self.CalculateAlignmentLocInI(c, trna_start, trna_start)
            profiles["start_pos"] = share.plusNumList(profiles["start_pos"], loc[0] - 1, loc[0] - 1, read_num)
            loc = self.CalculateAlignmentLocInI(c, trna_end, trna_end)
            profiles["end_pos"] = share.plusNumList(profiles["end_pos"], loc[1] - 1, loc[1] - 1, read_num)
            profiles["total_reads_num"] += read_num
        profiles["pileup_height"]=max(profiles["total"])
        profiles["cleavage_site_dic"] = self.FindCleavageSites(profiles["start_pos"],profiles["end_pos"])
        return profiles


def getTRFType(tc, s, e, c_offset=2, t_offset=3):
    t_s = tc[0]+1
    c_s = tc[1]+1
    c_e = tc[2]+1
    t_e = tc[3]+1
    if (s >= t_s - t_offset and s <= t_s + t_offset) and (e >= t_e - t_offset and e <= t_e + t_offset):
        return "full_tRNA"
    elif (s < t_s - t_offset) and (e > t_e + t_offset):
        return "full_U_tRNA"
    elif s<t_s-2 and (e<=c_e+c_offset and e>=c_s-c_offset):
        return "5_U_tRNA_halve"
    elif s < t_s-t_offset and e < c_s-c_offset-1:
        return "5_U_tRF"
    elif (s>=t_s-t_offset and s<=t_s+t_offset) and (e>=c_s-c_offset and e<=c_e+c_offset):
        return "5_tRNA_halve"
    elif (s >= t_s-t_offset and s <= t_s+t_offset) and (e<c_s-c_offset or e>c_s+c_offset) and (e<t_e+t_offset):
        return "5_tRF"
    elif (s>=c_s-c_offset and s<=c_e+c_offset) and (e >= t_e - t_offset and e <= t_e + t_offset):
        return "3_tRNA_halve"
    elif (s > t_s + t_offset) and (s > c_e + c_offset or s < c_s - c_offset) and (
            e >= t_e - t_offset and e <= t_e + t_offset):
        return "3_tRF"
    elif s >c_e+c_offset and e>t_e+t_offset:
        return "3_U_tRF"
    elif (s>=c_s-c_offset and s<=c_e+c_offset) and e > t_e + t_offset:
        return "3_U_tRNA_halve"
    elif (s >= t_s + t_offset) and e <= t_e - t_offset:
        return "i-tRF"
    else:
        return "other"

#type = getTRFType([60,93,95,132],98,129)  #3_tRF