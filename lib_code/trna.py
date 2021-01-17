#  tRNAExplorer v.1.0 - tRNA profiles and regulation networks
#  Copyright (C) 2020  Dr. Quanyuan He
#  School of Medicine, Hunan Normal University
#  Email: hqyone@hotmail.com
#  Freely distributed under the GNU General Public License (GPLv3)
#

import lib_code.share as share
import re

def CheckTabStr(func):
    def check(obj,tabstr):
        contents = tabstr.strip().split("\t")
        min_len = len(obj.GetTabTitle().split("\t"))
        if len(contents)>=min_len:
            func(obj, tabstr)
        else:
            print (f"{tabstr} doesn't contain enough {min_len} fields")   
    return check

class Seq():
    __slots__=["name","seq"]
    def __init__(self, **kwargs):
        self.name=""
        self.seq=""
        for k, v in kwargs.items():
            setattr(self, k, v)

    def GetTabTitle(self):
        return "name\tseq"

    def GetTabStr(self):
        return f"{self.name}\t{self.seq}"
    
    @CheckTabStr
    def LoadStr(self, tabstr):
        contents = tabstr.strip().split("\t")
        self.name = contents[0]
        self.seq = contents[1]        

class tRNA_Loop(Seq):
    __slots__ = ["start", "l_start", "l_end", "end", "struct_str","for_str","rev_str","loop_str","type"]
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.start=-1
        self.l_start=-1
        self.l_end=-1
        self.end=-1
        self.struct_str=""
        self.for_str=""
        self.rev_str=""
        self.loop_str=""
        self.type=""
        for k, v in kwargs.items():
            setattr(self, k, v)
                    
    def GetTabTitle(self):
        return f"{self.type}_loop_start\t{self.type}_loop_lstart\t{self.type}_loop_lend\t{self.type}_loop_end"
    
    def GetTabStr(self):
        return f"{self.start}\t{self.l_start}\t{self.l_end}\t{self.end}"

    @CheckTabStr
    def LoadStr(self, tabstr):
        contents = tabstr.strip().split("\t")
        self.start = int(contents[0])
        self.l_start = int(contents[1])
        self.l_end = int(contents[2])
        self.end = contents[3] 

class tRNA(Seq):
    __slots__ = ["name", "chrom", "start", "end", "strand",
                 "full_seq", "seq", "intron_infor", "anticodon", "anticodon_start",
                 "anticodon_end", "acceptor", "family", "seq_utr", "utr_len", "map_start",
                 "map_end", "map_len", "map_structure_str", "map_seq", "anticodon_start_in_map",
                 "anticodon_end_in_map", "map_scan_score", "possible_type", "d_loop", "a_loop",
                 "v_loop", "t_loop", "stem_for", "stem_rev"]

    def __init__(self, **kwargs):
        self.name = ""
        self.chrom = ""
        self.start = ""
        self.end = ""
        self.strand = ""
        self.full_seq = ""   # May contain introns
        self.seq = ""        # without_intron
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
        self.d_loop = tRNA_Loop(type="d")
        self.a_loop = tRNA_Loop(type="a")
        self.v_loop = tRNA_Loop(type="v")
        self.t_loop = tRNA_Loop(type="t")
        self.stem_for={"start":-1, "end":-1}
        self.stem_rev={"start":-1, "end":-1}
        for k, v in kwargs.items():
            setattr(self, k, v)

    def GetTabTitle(self):
        return "#chrom\tstart\tend\tname\tscore\tstrand\tfull_seq\tseq\tintron_infor\tanticodon\tanticodon_start" + \
               "\tanticodon_end\tacceptor" + \
               "\tfamily\tseq_utr\tutr_len\tmap_start\tmap_end\tmap_len\tmap_structure_str\tmap_seq" +\
               "\tanticodon_start_in_map\tanticodon_end_in_map\tmap_scan_score\tpossible_type" +\
               "\t"+self.d_loop.GetTabTitle() + \
               "\t"+self.a_loop.GetTabTitle() + \
               "\t"+self.v_loop.GetTabTitle() + \
               "\t"+self.t_loop.GetTabTitle() + \
               "\tstem_for_start\tstem_for_end\tstem_rev_start\tstem_rev_end"

    # Print a bed like file
    def GetTabStr(self):
        a = str(self.chrom)+"\t"+str(self.start)+"\t"+str(self.end)+"\t"+self.name+"\t0\t"+str(self.strand)+"\t"+str(self.full_seq)+"\t"+str(self.seq)+"\t" + \
            self.intron_infor + "\t" + self.anticodon+"\t"+str(self.anticodon_start) + "\t"+str(self.anticodon_end) + "\t" + \
            str(self.acceptor)+"\t"+str(self.family) + "\t" + str(self.seq_utr) + "\t" + str(self.utr_len) + "\t" +\
            str(self.map_start)+"\t"+str(self.map_end) + "\t"+str(self.map_len) + "\t"+self.map_structure_str + "\t"+self.map_seq + "\t" + \
            str(self.anticodon_start_in_map)+"\t"+str(self.anticodon_end_in_map)+"\t"+str(self.map_scan_score)+"\t"+str(self.possible_type)+"\t" + \
            self.d_loop.GetTabStr() +"\t"+\
            self.a_loop.GetTabStr() +"\t"+\
            self.v_loop.GetTabStr() +"\t"+\
            self.t_loop.GetTabStr() +"\t"+\
            str(self.stem_for["start"]) + "\t" + str(self.stem_for["end"]) + "\t" + \
            str(self.stem_rev["start"]) + "\t" + str(self.stem_rev["end"])
        return a

    @CheckTabStr
    def LoadStr(self, str):
        contents = str.strip().split("\t")
        self.chrom = contents[0]
        self.start = int(contents[1])
        self.end = int(contents[2])
        self.name = contents[3]
        self.strand = contents[5]
        self.full_seq = contents[6]
        self.seq = contents[7]
        self.intron_infor = contents[8]
        self.anticodon = contents[9]
        self.anticodon_start = int(contents[10])
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
        self.d_loop = tRNA_Loop(type="d", start=int(contents[25]), l_start=int(contents[26]), l_end=int(contents[27]), end=int(contents[28]))
        self.a_loop = tRNA_Loop(type="d", start=int(contents[29]), l_start=int(contents[30]), l_end=int(contents[31]), end=int(contents[32]))
        self.v_loop = tRNA_Loop(type="d", start=int(contents[33]), l_start=int(contents[34]), l_end=int(contents[35]), end=int(contents[36]))
        self.t_loop = tRNA_Loop(type="d", start=int(contents[37]), l_start=int(contents[38]), l_end=int(contents[39]), end=int(contents[40]))
        self.stem_for = {'start': int(contents[41]), 'end': int(contents[42]), 'str': ""}
        self.stem_rev = {'start': int(contents[43]), 'end': int(contents[44]), 'str': ""}

    def GetMatureSeq(self):
        if self.intron_infor != "":
            s = re.search(
                r'intron:\s+(\d+)-(\d+)\s+\((\d+)-(\d+)\)', self.intron_infor)
            s_start = 0
            s_end = 0
            if s:
                s_start = int(s.group(1))
                s_end = int(s.group(2))
            return self.full_seq[0:s_start-1]+self.full_seq[s_end:]
        else:
            return self.full_seq

    def GetPreMatureNoIntronSeq(self):
        if self.intron_infor != "":
            s = re.search(
                r'intron:\s+(\d+)-(\d+)\s+\((\d+)-(\d+)\)', self.intron_infor)
            s_start = 0
            s_end = 0
            if s:
                s_start = int(s.group(3))
                s_end = int(s.group(4))
            return self.seq_utr[0:s_start-1]+self.seq_utr[s_end:]
        else:
            return self.seq_utr

    def GetIntronLocationInI(self):
        if self.intron_infor != "":
            s = re.search(
                r'intron:\s+(\d+)-(\d+)\s+\((\d+)-(\d+)\)', self.intron_infor)
            s_start = 0
            s_end = 0
            if s:
                s_start = int(s.group(3))
                s_end = int(s.group(4))
            return [s_start, s_end]
        else:
            return [-1, -1]

    def GetKeySitesInfor(self, trans_type):
        [intron_start, intron_end] = self.GetIntronLocationInI()
        intron_len = intron_end-intron_start+1
        anti_codon_start = int(self.map_start + self.anticodon_start_in_map)
        anti_codon_end = int(self.map_start + self.anticodon_end_in_map)
        trna_start_pos = int(self.map_start)
        trna_end_pos = int(self.map_end)
        if trans_type == "I":
            return [int(self.map_start), anti_codon_start,
                    anti_codon_end, int(self.map_end)]
        elif trans_type == "P":
            if intron_start >= 0:  # If there are is intron
                if anti_codon_start > intron_start:
                    anti_codon_start = anti_codon_start - intron_len
                if anti_codon_end > intron_start:
                    anti_codon_end = anti_codon_end - intron_len
                if trna_end_pos > intron_start:
                    trna_end_pos -= intron_len
        elif trans_type == "M":
            if intron_start >= 0:
                if anti_codon_start > intron_start:
                    anti_codon_start = anti_codon_start - intron_len
                if anti_codon_end > intron_start:
                    anti_codon_end = anti_codon_end - intron_len
                if trna_end_pos > intron_start:
                    trna_end_pos -= intron_len
            anti_codon_start -= (trna_start_pos-1)
            anti_codon_end -= (trna_start_pos-1)
            trna_end_pos -= (trna_start_pos-1)
            trna_start_pos = 1
        elif trans_type == "C":
            if intron_start >= 0:
                if anti_codon_start > intron_start:
                    anti_codon_start = anti_codon_start - intron_len
                if anti_codon_end > intron_start:
                    anti_codon_end = anti_codon_end - intron_len
                if trna_end_pos > intron_start:
                    trna_end_pos -= intron_len
            anti_codon_start -= (trna_start_pos - 1)
            anti_codon_end -= (trna_start_pos - 1)
            trna_end_pos -= (trna_start_pos - 1)
            trna_end_pos += 3
            trna_start_pos = 1
        else:
            return []
        return [trna_start_pos, anti_codon_start, anti_codon_end, trna_end_pos]

    def IsQualifiedtRNA(self, no_mit_tRNA=True, no_pseudogenes=True, min_qscore=30):
        if no_mit_tRNA and self.name.startswith("nm"):
            return False
        if no_pseudogenes and self.possible_type == "pseudogene":
            return False
        if self.map_scan_score < min_qscore:
            return False
        return True

    # Calculate the position of reads in the tRNA with intron and UTR regions
    # read_class : nine styles of tRNA:
    # start: start position of read in object position
    # end: end position of read in object position
    # I,1,75
    def CalculateAlignmentLocInI(self, read_class, start, end):
        offset = self.utr_len
        if not (read_class == "A" or read_class == "B" or read_class == "E" or read_class == "G"):
            start += offset
            end += offset
        # Remove the effect of CCA
        if read_class == "I":
            if end > start-3:
                end -= 3
        if self.intron_infor != "":
            s = re.search(
                r'intron:\s+(\d+)-(\d+)\s+\((\d+)-(\d+)\)', self.intron_infor)
            s_start = 0
            s_end = 0
            if s:
                s_start = int(s.group(3))
                s_end = int(s.group(4))
                if self.intron_infor != "" and (read_class == "G" or read_class == "H" or read_class == "I"):
                    if start > s_start:
                        start += s_end-s_start+1
                    if end > s_start:
                        end += s_end-s_start+1
            return [start, end]
        else:
            return [start, end]

    # Transform Mutation location (the tRNA with intron and UTR regions) with to the location in mature tRNA
    # read_class : nine styles of tRNA:
    # loc: start position of read in object position 1based
    # I,1,75
    def CalculateMutLocItoM(self, loc):
        offset = self.utr_len
        if loc <= offset or loc > len(self.full_seq) + offset:
            return -1
        else:
            if self.intron_infor != "":
                s = re.search(
                    r'intron:\s+(\d+)-(\d+)\s+\((\d+)-(\d+)\)', self.intron_infor)
                s_start = 0
                s_end = 0
                if s:
                    s_start = int(s.group(3))
                    s_end = int(s.group(4))
                    if loc >= s_start and loc <= s_end:
                        return -1
                    elif loc < s_start:
                        return loc-offset
                    elif loc > s_end:
                        return loc-(s_end-s_start+1)-offset
                else:
                    return -1
            else:
                return loc-offset

    # The position should be pos of URL_Seq_INTRON
    # pos is 1 based

    def FindPosType(self, pos):
        ptype = "Other"
        if pos < -1:  # Upsteam of cleavage
            ptype = "5-UTR"
        elif pos == -1:
            ptype = "Start"
        elif pos > len(self.full_seq)-1:
            ptype = "3-UTR"
        elif pos == len(self.full_seq)-1:
            ptype = "End"
        elif pos >= self.anticodon_start-self.utr_len and pos <= self.anticodon_end-self.utr_len:
            ptype = "Anticodon"
        elif self.d_loop.start != -1 and pos >= self.d_loop.start and pos <= self.d_loop.end:
            ptype = "D-loop"
        elif self.a_loop.start != -1 and pos >= self.a_loop.start and pos <= self.a_loop.end:
            ptype = "A-loop"
        elif self.v_loop.start != -1 and pos >= self.v_loop.start and pos <= self.v_loop.end:
            ptype = "V-loop"
        elif self.t_loop.start != -1 and pos >= self.t_loop.start and pos <= self.t_loop.end:
            ptype = "T-loop"
        elif self.stem_for["start"] != -1 and pos >= self.stem_for["start"] and pos <= self.stem_for["end"]:
            ptype = "stem_for"
        elif self.stem_rev["start"] != -1 and pos >= self.stem_rev["start"] and pos <= self.stem_rev["end"]:
            ptype = "stem_rev"
        return ptype

    def FindCleavageSites(self, start_pos_profile, end_pos_profile, min_intensity=100, min_sn_ratio=30, context_seq_len=8):
        cleavageSite_Dic = {}
        combined_profile = start_pos_profile[:]
        if len(start_pos_profile) == len(end_pos_profile):
            #for i in range(len(end_pos_profile)):
            for i, v in enumerate(end_pos_profile):
                combined_profile[i] += v
        newList = sorted(combined_profile)
        median = newList[int(len(newList)/2)]
        start_peak_dic = self.DetectPeak(
            start_pos_profile, "startpos", context_seq_len)
        end_peak_dic = self.DetectPeak(
            end_pos_profile, "endpos", context_seq_len)
        # Find the peak list by combined start and end profiles
        start_peak_ls = list(start_peak_dic.keys())
        end_peak_ls = list(end_peak_dic.keys())
        combined_peak_id_ls = list(set(start_peak_ls + end_peak_ls))
        for id in combined_peak_id_ls:
            combined_peak = {
                "pos": 0,
                "int": 0,
                "int5": 0,
                "int3": 0,
                "sn_ratio": 0,
                "c_5_seq": "",
                "c_3_seq": ""
            }
            if id in start_peak_dic:
                start_peak = start_peak_dic[id]
                combined_peak["pos"] = start_peak["pos"]
                combined_peak["int"] += start_peak["int"]
                combined_peak["ptype"] = start_peak["ptype"]
                combined_peak["int3"] = start_peak["int"]
                combined_peak["c_5_seq"] = start_peak["c_5_seq"]
                combined_peak["c_3_seq"] = start_peak["c_3_seq"]
            if id in end_peak_dic:
                end_peak = end_peak_dic[id]
                combined_peak["pos"] = end_peak["pos"]
                combined_peak["int"] += end_peak["int"]
                combined_peak["ptype"] = end_peak["ptype"]
                combined_peak["int5"] = end_peak["int"]
                combined_peak["c_5_seq"] = end_peak["c_5_seq"]
                combined_peak["c_3_seq"] = end_peak["c_3_seq"]
            if median == 0:
                median = 1
            if combined_peak["int"] > min_intensity and (float(combined_peak["int"])/median) > min_sn_ratio:
                combined_peak['sn_ratio'] = int(
                    float(combined_peak["int"])/median)
                cleavageSite_Dic[id] = combined_peak
        return cleavageSite_Dic

    def DetectPeak(self, profile, type, context_seq_len=5):
        potential_peak = {}
        for i in range(1, len(profile)-1):
            if profile[i-1] < profile[i] and profile[i] > profile[i+1]:
                pos = i
                # pos is 0 based on
                context_5_seq = ""
                context_3_seq = ""
                if type == "startpos":
                    context_5_seq = self.seq_utr[max(
                        0, pos-context_seq_len):min(len(self.seq_utr), pos)]
                    context_3_seq = self.seq_utr[max(0, pos): min(
                        len(self.seq_utr), pos+context_seq_len)]
                elif type == "endpos":
                    context_5_seq = self.seq_utr[max(
                        0, pos - context_seq_len+1):min(len(self.seq_utr), pos+1)]
                    context_3_seq = self.seq_utr[max(
                        0, pos+1): min(len(self.seq_utr), pos + context_seq_len+1)]
                if type == "startpos":
                    pos -= 1
                potential_peak["P"+str(pos)] = {
                    "pos": pos-self.utr_len,
                    "int": profile[i],
                    "ptype": self.FindPosType(pos-self.utr_len),
                    # "sn_ratio":round(float(profile[i])/median,3),
                    "c_5_seq": context_5_seq,
                    "c_3_seq": context_3_seq
                }
        return potential_peak

    # I,1,75
    def CreateTrfProfiles(self, brief_mapping_infor_ls, offset):
        profile_length = len(self.seq_utr)
        profiles = {
            "total_reads_num": 0,
            "isequence": self.seq_utr,
            "pileup_height": 0,
            "start_pos": [0] * profile_length,
            "end_pos": [0] * profile_length,
            "total": [0] * profile_length,
            "cleavage_site_dic": {},
            "mutation_dic_str": ""
        }
        # Both "start_pos" and "end_pos" are 1 based
        for bi in brief_mapping_infor_ls:
            i = bi.split(",")
            c = i[0]  # I
            trna_start = int(i[1])  # 1
            trna_end = int(i[2])  # 75
            read_num = float(i[3])
            loc = self.CalculateAlignmentLocInI(c, trna_start, trna_end)
            profiles["total"] = share.plusNumList(
                profiles["total"], loc[0] - 1, loc[1] - 1, read_num)
            loc = self.CalculateAlignmentLocInI(c, trna_start, trna_start)
            profiles["start_pos"] = share.plusNumList(
                profiles["start_pos"], loc[0] - 1, loc[0] - 1, read_num)
            loc = self.CalculateAlignmentLocInI(c, trna_end, trna_end)
            profiles["end_pos"] = share.plusNumList(
                profiles["end_pos"], loc[1] - 1, loc[1] - 1, read_num)
            profiles["total_reads_num"] += read_num
        profiles["pileup_height"] = max(profiles["total"])
        profiles["cleavage_site_dic"] = self.FindCleavageSites(
            profiles["start_pos"], profiles["end_pos"])
        profiles["mutation_dic_str"] = self.CreateMutationDicStr(
            brief_mapping_infor_ls, profiles["total"])
        return profiles

    # The brief_mapping_infor can derived from transcript with intron or not
    # The final position conside intron.
    def BriefMappingInfor2MutationDic(self, brief_mapping_infor, total_profile):
        # Mutation IDs
        # Mutation: M:B>A:Location->Reads
        # Deletion: D:B:Location->Reads
        # Insertion: I:A:Location->Reads
        i = brief_mapping_infor.split(",")
        c = i[0]  # I
        trna_start = int(i[1])  # 1
        trna_end = int(i[2])  # 75
        loc = self.CalculateAlignmentLocInI(c, trna_start, trna_start)
        read_num = round(float(i[3]), 3)
        read_seq = i[4]
        ref_seq = i[5]
        mut_dic = {}
        if read_seq == ref_seq:
            return mut_dic
        if len(read_seq) == len(ref_seq):
            #for i in range(0, len(read_seq)):
            for index, A in enumerate(read_seq):
                B = ref_seq[index]
                Mutaion_Str = ""
                location = loc[0]+index  # 1 based
                if (location-1) < len(total_profile):
                    total_intensity = round(
                        float(total_profile[location-1]), 3)
                    if A == "-":
                        # Deletion
                        Mutaion_Str = f"D:{B}:{location}:{total_intensity}"
                    elif B == "-":
                        # Insertion
                        Mutaion_Str = f"I:{A}:{location}:{total_intensity}"
                    elif A != B:
                        # Mutaion
                        Mutaion_Str = f"M:{B}>{A}:{location}:{total_intensity}"
                    else:
                        # Unexpect
                        continue
                else:
                    print("Some thing unexpected mapping infor:" +
                          brief_mapping_infor)
                    print(f"{location-1}:{len(total_profile)-1}: {self.name}: {self.full_seq}")
                if Mutaion_Str != "":
                    if Mutaion_Str not in mut_dic:
                        mut_dic[Mutaion_Str] = 0
                    mut_dic[Mutaion_Str] += read_num
        return mut_dic

    def CreateMutationDicStr(self, brief_mapping_infor_ls, total_profile):
        # Both "start_pos" and "end_pos" are 1 based
        mut_dic = {}
        for bi in brief_mapping_infor_ls:
            cur_dic = self.BriefMappingInfor2MutationDic(bi, total_profile)
            for k in cur_dic:
                if k not in mut_dic:
                    mut_dic[k] = cur_dic[k]
                else:
                    mut_dic[k] += cur_dic[k]
        mut_str_ls = []
        for key in mut_dic:
            value = mut_dic[key]
            mut_str_ls.append(key+"="+str(value))
        if len(mut_str_ls) == 0:
            return '-'
        else:
            return ",".join(mut_str_ls)


def getTRFType(tc, s, e, c_offset=2, t_offset=3):
    if len(tc) == 4:
        t_s = tc[0]+1
        c_s = tc[1]+1
        c_e = tc[2]+1
        t_e = tc[3]+1
        if (s >= t_s - t_offset and s <= t_s + t_offset) and (e >= t_e - t_offset and e <= t_e + t_offset):
            return "full_tRNA"
        elif (s < t_s - t_offset) and (e > t_e + t_offset):
            return "full_U_tRNA"
        elif s < t_s-2 and (e <= c_e+c_offset and e >= c_s-c_offset):
            return "5_U_tRNA_halve"
        elif s < t_s-t_offset and e < c_s-c_offset-1:
            return "5_U_tRF"
        elif (s >= t_s-t_offset and s <= t_s+t_offset) and (e >= c_s-c_offset and e <= c_e+c_offset):
            return "5_tRNA_halve"
        elif (s >= t_s-t_offset and s <= t_s+t_offset) and (e < c_s-c_offset or e > c_s+c_offset) and (e < t_e+t_offset):
            return "5_tRF"
        elif (s >= c_s-c_offset and s <= c_e+c_offset) and (e >= t_e - t_offset and e <= t_e + t_offset):
            return "3_tRNA_halve"
        elif (s > t_s + t_offset) and (s > c_e + c_offset or s < c_s - c_offset) and (
                e >= t_e - t_offset and e <= t_e + t_offset):
            return "3_tRF"
        elif s > c_e+c_offset and e > t_e+t_offset:
            return "3_U_tRF"
        elif (s >= c_s-c_offset and s <= c_e+c_offset) and e > t_e + t_offset:
            return "3_U_tRNA_halve"
        elif (s >= t_s + t_offset) and e <= t_e - t_offset:
            return "i-tRF"
        else:
            return "other"
    else:
        return "unknow"

# type = getTRFType([60,93,95,132],98,129)  #3_tRF
