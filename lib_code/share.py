# coding=utf-8

#  tRNAExplorer v.1.0 - tRNA profiles and regulation networks
#  Copyright (C) 2020  Dr. Quanyuan He
#  School of Medicine, Hunan Normal University
#  Email: hqyone@hotmail.com
#  Freely distributed under the GNU General Public License (GPLv3)
#

import os
import re
from pyfaidx import Fasta

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-':'-'}
def reverse_complement(seq):
    return "".join(complement.get(base, base) for base in reversed(seq))

def getFasta(ref_fasta, chr, start, end):
    genes = Fasta(ref_fasta)
    return genes[chr][start:end].__str__().upper()

def getExtFileList(dir,ext):
    file_ls = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(ext):
                file_ls.append(os.path.join(root, file))
    return file_ls

def plusNumList(ls, start, end, value):
    try:
        for i in range(start,end+1):
            if i <len(ls) and value:
                ls[i] = float(ls[i])+float(value)
        return ls
    except Exception as a:
        print(str(start)+":"+str(end))
        print(a)

def stringfyNumList(ls):
    str_ls = []
    for i in ls:
        str_ls.append(str(int(i)))
    return str_ls

def trimseq(rawseq,  f_patterm="", r_patterm=""):
    seq = rawseq
    if f_patterm != "":
        p = "^" + f_patterm.replace('N', '\w')
        b = re.search(p, seq)
        if b:
            seq = seq[b.span()[1]:]
    if r_patterm != "":
        p = "CCA" + r_patterm.replace('N', '\w')
        b = re.search(p, seq)
        if b:
            seq = seq[:b.span()[0] + 3]
        else:
            p = r_patterm.replace('N', '\w') + "$"
            b = re.search(p, seq)
            if b:
                seq = seq[:b.span()[0]]
    return seq
# filter out reads with quality lower than 28 at any position
def filterFastQ2FastA(fastq, filtered_fasta, num_dic_txt, f_patterm="", r_patterm="", qcutoff=0, num_cutoff=50):
    big_file = False
    if os.path.getsize(fastq)>3259059200:
        big_file=True
    dic = {}
    name_dic = {}
    total_num = 0
    low_quality_num=0
    high_quality_num = 0
    line_index=0
    seq=""
    title=""
    title2=""
    quanlity=""
    with open(fastq,'r') as FASTQ, open(filtered_fasta, 'w') as Filterecd_FASTA, open(num_dic_txt,"w") as NumDicFile:
        for line in FASTQ:
            line  = line.strip()
            line_index+=1
            if line_index%4==2:
                seq = trimseq(line, f_patterm,r_patterm)
            elif line_index%4==3:
                continue
                # title2 = line.replace(" ","_")
            elif line_index%4==0:
                quanlity = line
            elif line.startswith('@') and line_index%4==1 and seq!="":
                low_quanlity = False
                title=line.replace(" ","_")
                if total_num % 1000000==0:
                    print("reads ("+str(total_num)+")")
                    seq_ls= list(dic.keys())
                    if big_file:
                        for seq in seq_ls:
                            if (dic[seq]["count"]==1):
                                del dic[seq]
                                del name_dic[seq]
                    #print("dic ("+str(len(dic.keys()))+")")
                    # gc.collect()
                for i in quanlity:
                    if ord(i)<=qcutoff+33:  #Based on Phred 33
                        low_quanlity = True
                        break
                if not low_quanlity :
                    high_quality_num+=1
                    if seq not in dic:
                        dic[seq] = {
                            "title":title,
                            #"seq":seq,
                            "count":1
                        }
                        name_dic[seq] = title
                    else:
                        dic[seq]["count"]+=1
                total_num+=1
        total_num += 1
        non_redundent_num = 0
        for seq in dic:
            title = name_dic[seq]
            if dic[seq]["count"]>=num_cutoff:
                Filterecd_FASTA.write(">" + dic[seq]["title"] + "\n")
                Filterecd_FASTA.write(seq + "\n")
                non_redundent_num += 1
                NumDicFile.write(title+"\t"+str(dic[seq]["count"])+"\t"+seq+"\n")
    return({"total_num":int(total_num),
        "removed_num":int(total_num-high_quality_num),
        "survived_num":int(high_quality_num),
        "non_redundent_num":int(non_redundent_num)})


# filter out reads with quality lower than 28 at any position and have small numbers
def filterFastQ(fastq, filtered_fastq, num_dic_txt, cutoff=28):
    FASTQ  = open(fastq,'r')
    iFASTQ = iter(FASTQ)
    dic = {}
    name_dic = {}
    Filterecd_FASTQ = open(filtered_fastq, 'w')
    NumDicFile = open(num_dic_txt,"w")
    total_num = 0
    low_quality_num=0
    high_quality_num = 0
    non_redundent_num=0
    for line in iFASTQ:
        line  = line.strip()
        if line.startswith('@'):
            total_num+=1
            title = line
            seq = next(iFASTQ).strip()
            title2 = next(iFASTQ).strip()
            quanlity = next(iFASTQ).strip()
            low_quanlity = False

            for i in quanlity:
                if ord(i)<=cutoff+33:  #Based on Phred 33
                    low_quanlity = True
                    break
            if not low_quanlity:
                high_quality_num+=1
                if seq not in dic:
                    Filterecd_FASTQ.write(title + "\n")
                    Filterecd_FASTQ.write(seq + "\n")
                    Filterecd_FASTQ.write(title2 + "\n")
                    Filterecd_FASTQ.write(quanlity + "\n")
                    name_dic[seq] = title
                    non_redundent_num+=1
                else:
                    dic[seq]+=1
            else:
                low_quality_num+=1
    for seq in dic:
        title = name_dic[seq]
        NumDicFile.write(title+"\t"+str(dic[seq])+"\n")
    NumDicFile.close()
    Filterecd_FASTQ.close()
    return({"trimed_num":total_num,
          "removed_num":low_quality_num,
          "survived_num":high_quality_num,
          "non_redundent_num":non_redundent_num})


# test
# F="/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/fastq/SRR1836126_1.fastq"
# FF="/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/fastq/SRR1836126_1_filtered.fastq"
# T="/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/fastq/SRR1836126_1.num_dic.tab"
# filterFastQ(F,FF,T,25)