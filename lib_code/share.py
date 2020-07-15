# coding=utf-8

#  tRNAExplorer v.1.0 - tRNA profiles and regulation networks
#  Copyright (C) 2020  Dr. Quanyuan He
#  School of Medicine, Hunan Normal University
#  Email: hqyone@hotmail.com
#  Freely distributed under the GNU General Public License (GPLv3)
#

import os
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

# filter out reads with quality lower than 28 at any position
def filterFastQ2FastA(fastq, filtered_fasta, num_dic_txt, qcutoff=28, num_cutoff=50):
    FASTQ  = open(fastq,'r')
    iFASTQ = iter(FASTQ)
    dic = {}
    name_dic = {}
    Filterecd_FASTA = open(filtered_fasta, 'w')
    NumDicFile = open(num_dic_txt,"w")
    total_num = 0
    low_quality_num=0
    high_quality_num = 0
    for line in iFASTQ:
        line  = line.strip()
        if line.startswith('@'):
            total_num+=1
            title = line.replace(" ","_")
            seq = next(iFASTQ).strip()
            title2 = next(iFASTQ).strip()
            quanlity = next(iFASTQ).strip()
            low_quanlity = False

            for i in quanlity:
                if ord(i)<=qcutoff+33:  #Based on Phred 33
                    low_quanlity = True
                    break
            if not low_quanlity:
                high_quality_num+=1
                if seq not in dic:
                    dic[seq] = {
                        "title":title,
                        "seq":seq,
                        "count":1
                    }
                    name_dic[seq] = title
                else:
                    dic[seq]["count"]+=1
            else:
                low_quality_num+=1
    non_redundent_num = 0
    for seq in dic:
        title = name_dic[seq]
        if dic[seq]["count"]>=num_cutoff:
            Filterecd_FASTA.write(">" + dic[seq]["title"] + "\n")
            Filterecd_FASTA.write(dic[seq]["seq"] + "\n")
            non_redundent_num += 1
            NumDicFile.write(title+"\t"+str(dic[seq]["count"])+"\t"+dic[seq]["seq"]+"\n")
    NumDicFile.close()
    Filterecd_FASTA.close()
    return({"total_num":int(total_num),
          "removed_num":int(low_quality_num),
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