import os, subprocess
import lib_code.share as share
import lib_code.trna as trna
from lib_code.cmd_tools import BLASTN

# Indexing the tRNA FASTA file
def CreateBLASTdb(fasta):
    if os.path.isfile(fasta):
        trna_db_dir= os.path.dirname(os.path.abspath(fasta))
        base = os.path.basename(fasta)
        out_name = os.path.splitext(base)[0]
        cmd_file = trna_db_dir+"/cmd.sh"
        CMD_FILE = open(cmd_file, "w")
        b = BLASTN()
        # The second step do blast
        # The manual of BLAST+ is https://www.ncbi.nlm.nih.gov/books/NBK279690/
        BLAST_CMD = b.getCreateBLASTdbCMD(fasta,out_name)+"\n"
        CMD_FILE.write(BLAST_CMD)
        CMD_FILE.close()
        # os.popen("bash " + blast_cmd_file)
        process = subprocess.Popen("bash " + cmd_file, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print (process.returncode)
        return 0
    else:
        return -1

# Run BLAST for a bam file
def RunBLASTN(blastn, mkdb, id, db_fasta, query_fasta, blast_out_dir, eval=0.01, hit_number=30):
    try:
        print("Begin process "+id+" ...")
        cmd_bash = blast_out_dir+"/"+id+"_blast.sh"
        CMD_FILE = open(cmd_bash, "w")

        blast_out_file = blast_out_dir + "/"+id+"_tRNA_blast_out.tab"
        #qseqid sseqid nident sstart send
        blastn = BLASTN(blastn,mkdb)
        if not os.path.isfile(db_fasta+".nhr"):
            CMD_FILE.write(blastn.getCreateBLASTdbCMD(db_fasta, id)+"\n")
        blastn_cmd = blastn.getAlignmentCMD(db_fasta,query_fasta,eval,blast_out_file,hit_number)
        print(blastn_cmd)
        CMD_FILE.write(blastn_cmd)
        CMD_FILE.close()
        #os.popen("bash " + cmd_bash)
        process = subprocess.Popen("bash " + cmd_bash, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print (process.returncode)
        os.remove(cmd_bash)
        print("Finish processing " + id)
        return blast_out_file
    except():
        return ""

# Run BLAST for a bam file

# id="ENCFF715RQU"
# ds_dir = "/home/hqyone/python_code/testrepo/szf"
# db_fasta=ds_dir+"/hg38_tRNA_nointron_60.fasta"
# query_fasta="/home/hqyone/python_code/data/ENCFF715RQU_tRNA_60.fasta"
# blast_out_dir="/home/hqyone/python_code/testrepo/szf/blast_out"
# CreateBlastOutTab(id, db_fasta,query_fasta,blast_out_dir)


# Analysis BLAST output
# The output format is as following
# 1. read
# 2. tRNAs
# 3. percentage of identity
# 4. match length
# 5. start of read
# 6. end of read
# 7. start of tRNAs
# 8. end of tRNAs
# 9. e-value
# 10. score
# 12. qseq
# 13. sseq

# trna_read_dic is unique map reads number
def AnalysisBlastOut2(blast_out_file, read_num_dic_file, tRNA_dic, tRNA_reads_count_file, tRNA_reads_hit_file, url_len, tRF_Min_Length=19, max_mismatch=2):
    dist_dic = {}
    trna_unique_read_dic={}  #For quantifiation

    BLAST_OUT = open(blast_out_file, 'r')
    # Load tRNA structure annotations
    tRNA_anno_dic={}
    tRNA_id_family_dic = {}
    for t_id in tRNA_dic:
        t = tRNA_dic[t_id]
        tRNA_anno_dic[t.name]=t
        tRNA_id_family_dic[t.name] = t.family

    # Read the read numnber map
    READ_NUM = open(read_num_dic_file, 'r')
    read_num_dic ={}
    read_seq_dic={}
    for line in READ_NUM:
        contents = line.strip().split("\t")
        if len(contents)>1:
            read_num_dic[contents[0]] = int(contents[1])
        if len(contents)>2:
            read_seq_dic[contents[0]] = contents[2]
    READ_NUM.close()

    read_trna_dic = {}   # percent,
    # -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send qseq sseq qlen slen evalue"
    for line in BLAST_OUT:
        contents = line.strip().split("\t")
        cur_read_id = contents[0] #qseqid
        cur_tRNA_id = contents[1] #sseqid
        cur_percent = float(contents[2]) #pident
        cur_length = int(contents[3]) #length
        cur_mismatch = int(contents[4])  #mismatch
        cur_gap = int(contents[5]) #gaps
        read_start = float(contents[6]) #qstart
        read_end = int(contents[7])  #qend
        trna_start = float(contents[8]) #sstart
        trna_end = int(contents[9])  #send
        qseq = contents[10]  # qseq
        sseq = contents[11]  # sseq
        qlen = int(contents[12])  # qlen
        slen = int(contents[13])  # slen
        evalue = float(contents[14])  # evalue
        direction = "+"
        # As we use DNA to BLAST to mRNA so the sequences should be reversed complemental so the tRF ID will not redundent
        if trna_start>trna_end:
            a = trna_end
            trna_end = trna_start
            trna_start = a
            b = read_end
            read_end = read_start
            read_start = b
            qseq = share.reverse_complement(qseq)
            sseq = share.reverse_complement(sseq)
            direction = "-"
        #Cutoff here 98 means only reads length longer than 60 bp were allowed for two missmatch
        if cur_percent>=98 and cur_mismatch+cur_gap<=max_mismatch and len(sseq)>tRF_Min_Length:
            # Get the best hits for each read and accumulate matched tRNA id
            if cur_read_id not in read_trna_dic:
                read_trna_dic[cur_read_id] = {
                    "percent": cur_percent,
                    "length": cur_length,
                    "direction":direction,
                    "sseq":sseq,
                    "hit_type":"", # multiple or unique
                    "tRNAs": {}
                }
                c = cur_tRNA_id.split("::")
                if len(c)==2:
                    trna_seq_classes = c[0].strip()
                    trna_id = c[1].strip()
                    read_trna_dic[cur_read_id]["tRNAs"][trna_id]= {}
                for c in trna_seq_classes:
                    read_trna_dic[cur_read_id]["tRNAs"][trna_id][c]={
                            "read_start": read_start,
                            "read_end": read_end,
                            "trna_start": trna_start,
                            "trna_end": trna_end,
                            "qseq": qseq,
                            "sseq": sseq,
                            "qlen": qlen,  # read length
                            "slen": slen,
                            "evalue": evalue
                        }
            else:
                obj = read_trna_dic[cur_read_id]
                if cur_percent >= obj["percent"] and cur_length >= obj["length"] and cur_tRNA_id not in obj["tRNAs"] and (sseq == obj["sseq"] or sseq == share.reverse_complement(obj["sseq"])):
                    obj["percent"] = cur_percent
                    obj["length"] = cur_length
                    obj["sseq"] = sseq
                    c = cur_tRNA_id.split("::")
                    if len(c) == 2:
                        trna_seq_classes = c[0].strip()
                        trna_id = c[1].strip()
                        obj = read_trna_dic[cur_read_id]["tRNAs"]
                        if not trna_id in obj:
                            obj[trna_id] = {}
                    for c in trna_seq_classes:
                        obj[trna_id][c] = {
                            "read_start": read_start,
                            "read_end": read_end,
                            "trna_start": trna_start,
                            "trna_end": trna_end,
                            "qseq": qseq,
                            "sseq": sseq,
                            "qlen": qlen,
                            "slen": slen,
                            "evalue": evalue
                        }

    # Create mean reads number for every tRNA
    tRNA_mean_read_dic = {} #Seperate
    pretRNA_mean_read_dic = {}
    tRNA_total_obj_dic = {}
    for cur_read_id in read_trna_dic:
        # Get all unique read for tRNA
        if "tRNAs" in read_trna_dic[cur_read_id]:
            read_trna_dic[cur_read_id]["hit_tRNAs_num"] = len(read_trna_dic[cur_read_id]["tRNAs"].keys())
            for tRNA_id in read_trna_dic[cur_read_id]["tRNAs"]:
                number = read_num_dic[cur_read_id]
                mean_exp_number = float(number)/len(read_trna_dic[cur_read_id]["tRNAs"].keys())
                if not tRNA_id in tRNA_mean_read_dic:
                    tRNA_mean_read_dic[tRNA_id] = 0
                tRNA_mean_read_dic[tRNA_id]+= mean_exp_number

                if tRNA_id not in tRNA_total_obj_dic:
                    tRNA_total_obj_dic[tRNA_id] = {}
                if cur_read_id not in tRNA_total_obj_dic[tRNA_id]:
                    tRNA_total_obj_dic[tRNA_id][cur_read_id] = {}
                tRNA_total_obj_dic[tRNA_id][cur_read_id]["mean_exp"] = mean_exp_number
                tRNA_total_obj_dic[tRNA_id][cur_read_id]["strand"] = mean_exp_number
                tRNA_total_obj_dic[tRNA_id][cur_read_id]["hit_tRNAs_num"] =read_trna_dic[cur_read_id]["hit_tRNAs_num"]
                tRNA_total_obj_dic[tRNA_id][cur_read_id]["direction"] = read_trna_dic[cur_read_id]["direction"]
                if "classes" not in tRNA_total_obj_dic[tRNA_id][cur_read_id]:
                    tRNA_total_obj_dic[tRNA_id][cur_read_id]["classes"] = {}
                for c in read_trna_dic[cur_read_id]["tRNAs"][tRNA_id]:
                    tRNA_total_obj_dic[tRNA_id][cur_read_id]["classes"][c]={}
                    read_start = read_trna_dic[cur_read_id]["tRNAs"][tRNA_id][c]["read_start"]
                    read_end = read_trna_dic[cur_read_id]["tRNAs"][tRNA_id][c]["read_end"]
                    qseq = read_trna_dic[cur_read_id]["tRNAs"][tRNA_id][c]["qseq"]
                    sseq = read_trna_dic[cur_read_id]["tRNAs"][tRNA_id][c]["sseq"]
                    qlen = read_trna_dic[cur_read_id]["tRNAs"][tRNA_id][c]["qlen"]
                    slen = read_trna_dic[cur_read_id]["tRNAs"][tRNA_id][c]["slen"]
                    trna_start = read_trna_dic[cur_read_id]["tRNAs"][tRNA_id][c]["trna_start"]
                    trna_end = read_trna_dic[cur_read_id]["tRNAs"][tRNA_id][c]["trna_end"]
                    tRNA_total_obj_dic[tRNA_id][cur_read_id]["classes"][c] = {
                        "read_start": int(read_start),
                        "read_end": int(read_end),
                        "trna_start": int(trna_start),
                        "trna_end": int(trna_end),
                        "qseq": qseq,
                        "sseq": sseq,
                        "qlen": qlen,
                        "slen": slen
                        # "hit_seq":read_seq_dic[cur_read_id][start_read-1:end_read]
                    }

    # Write the result to files
    tRNA_READS_COUNT = open(tRNA_reads_count_file,'w')
    tRNA_READS_COUNT.write("#tRNA_family\ttRNA_id\ttotal_reads\n")
    for tRNA_id in tRNA_mean_read_dic:
        family_id = tRNA_id
        if tRNA_id in tRNA_id_family_dic:
            family_id = tRNA_id_family_dic[tRNA_id]
        total_reads = tRNA_mean_read_dic[tRNA_id]
        tRNA_READS_COUNT.write(family_id+"\t"+tRNA_id+"\t"+str(total_reads)+"\n")
    tRNA_READS_COUNT.close()

    tRNA_READS_HIT = open(tRNA_reads_hit_file,'w')
    tRNA_READS_HIT.write("#tRNA_family\ttRNA_id\tread_id\tdirection"+  #4
                         "\tI\tI_read_start\tI_read_end\tI_trna_start\tI_trna_end"+ #5
                         "\tP\tP_read_start\tP_read_end\tP_trna_start\tP_trna_end"+
                         "\tM\tM_read_start\tM_read_end\tM_trna_start\tM_trna_end"+
                         "\tC\tC_read_start\tC_read_end\tC_trna_start\tC_trna_end"+
                         "\t5T\t3T\t5C\t3C\tmean_number\thit_tRNAs_num\tTRF_type\tbrief_mapping_infor\tRead_type"+
                         "\tmapping_ratio\tread_5_fragment\tread_fragment\tread_3_fragment\n")
    for tRNA_id in tRNA_total_obj_dic:
        family_id = tRNA_id
        if tRNA_id in tRNA_id_family_dic:
            family_id = tRNA_id_family_dic[tRNA_id]
        for read_id in tRNA_total_obj_dic[tRNA_id]:
            mean_exp = tRNA_total_obj_dic[tRNA_id][read_id]["mean_exp"]
            hit_tRNAs_num = tRNA_total_obj_dic[tRNA_id][read_id]["hit_tRNAs_num"]
            class_obj = tRNA_total_obj_dic[tRNA_id][read_id]["classes"]
            direction = tRNA_total_obj_dic[tRNA_id][read_id]["direction"]

            # Get mapping location information
            M_5T = 0
            M_3T = 0
            M_5C = 0
            M_3C = 0
            I, I_read_start, I_read_end, I_trna_start, I_trna_end, I_qseq, I_sseq = 0, -1, -1, -1, -1,"",""
            P, P_read_start, P_read_end, P_trna_start, P_trna_end, P_qseq, P_sseq  = 0, -1, -1, -1, -1,"",""
            M, M_read_start, M_read_end, M_trna_start, M_trna_end, M_qseq, M_sseq  = 0, -1, -1, -1, -1,"",""
            C, C_read_start, C_read_end, C_trna_start, C_trna_end, C_qseq, C_sseq  = 0, -1, -1, -1, -1,"",""
            read_start = 0
            read_end = 0
            brief_mapping_infor = ""
            trna_start = 0
            trna_end = 0
            qseq = ""
            sseq = ""
            for c in class_obj:
                obj = class_obj[c]
                read_start = obj["read_start"]
                read_end = obj["read_end"]
                if c=="I":
                    I, I_read_start, I_read_end, I_trna_start, I_trna_end,I_qseq, I_sseq = 1, obj["read_start"], obj["read_end"],obj[
                        "trna_start"],obj["trna_end"],obj["qseq"], obj["sseq"]
                    trna_start, trna_end,qseq, sseq  = obj["trna_start"],obj["trna_end"],obj["qseq"], obj["sseq"]
                if c=="P":
                    P, P_read_start, P_read_end, P_trna_start, P_trna_end, P_qseq, P_sseq = 1, obj["read_start"], obj["read_end"], obj[
                        "trna_start"], obj["trna_end"],obj["qseq"], obj["sseq"]
                    trna_start, trna_end, qseq, sseq = obj["trna_start"], obj["trna_end"], obj["qseq"], obj["sseq"]
                if c=="M":
                    M, M_read_start, M_read_end, M_trna_start, M_trna_end, M_qseq, M_sseq= 1, obj["read_start"], obj["read_end"], obj[
                        "trna_start"], obj["trna_end"],obj["qseq"], obj["sseq"]
                    trna_start, trna_end, qseq, sseq = obj["trna_start"], obj["trna_end"], obj["qseq"], obj["sseq"]
                if c=="C":
                    C, C_read_start, C_read_end, C_trna_start, C_trna_end, C_qseq, C_sseq= 1, obj["read_start"], obj["read_end"], obj[
                        "trna_start"], obj["trna_end"],obj["qseq"], obj["sseq"]
                    trna_start, trna_end, qseq, sseq = obj["trna_start"], obj["trna_end"], obj["qseq"], obj["sseq"]
                if c=="I" or c=="P":
                    if obj["trna_start"]<url_len-2:
                        M_5C = 1
                    if obj["trna_end"]>obj["slen"]-url_len+1:
                        M_3C = 1
                else: # M or C
                    if obj["trna_start"]<=1:
                        M_5T = 1
                    if obj["trna_end"]>=obj["slen"]-1:
                        M_3T = 1
            Read_type = ""
            # Get read type based on (I, P, M, C) and location information
            # Detailed can be found figure 1 in manuscript
            if I==1 and P==1 and M==0 and C==0 and M_5T==0 and M_3T==0 and M_5C==0 and M_3C==0:
                Read_type += "A"
            if I==1 and P==1 and M==0 and C==0 and M_5C==1:
                Read_type += "B"
            if I==1 and P==1 and M==1 and C==1 and M_5T==1:
                Read_type += "C"
            if I==1 and P==1 and M==1 and C==1 and M_5T==0 and M_3T==0 and M_5C==0 and M_3C==0:
                Read_type += "D"
            if I==1 and P==0 and M==0 and C==0 and M_5T==0 and M_3T==0 and M_5C==0 and M_3C==0:
                Read_type += "E"
            if I==0 and P==1 and M==1 and C==1 and M_5T==0 and M_3T==0 and M_5C==0 and M_3C==0:
                Read_type += "F"
            if I==1 and P==1 and M==0 and C==0 and M_3C==1:
                Read_type += "G"
            if I==1 and P==1 and M==1 and C==1 and M_3T==1:
                Read_type += "H"
            if I==0 and P==0 and M==0 and C==1 and M_3T==1:
                Read_type += "I"
            if Read_type=="":
                Read_type = "U"

            if brief_mapping_infor == "":
                brief_mapping_infor = Read_type + "," + str(trna_start) + "," + str(trna_end) + "," + str(
                    round(mean_exp, 3)) + "," + qseq.replace("-", "") + "," + sseq.replace("-", "")
            read_seq = read_seq_dic[read_id]
            if direction == "-":
                read_seq = share.reverse_complement(read_seq)
            read_5_fragment = read_seq[0:read_start - 1]
            read_fragment = read_seq[read_start - 1:read_end]
            read_3_fragment = read_seq[read_end:]

            trf_type="Unknown"
            mapping_ratio = 0.0
            if tRNA_id in tRNA_anno_dic:
                if "P" in class_obj:
                    trf_type = trna.getTRFType(tRNA_anno_dic[tRNA_id].GetKeySitesInfor("P"), class_obj["P"]["trna_start"], class_obj["P"]["trna_end"])
                    mapping_ratio = (class_obj["P"]["trna_end"] - class_obj["P"]["trna_start"] + 1) / class_obj["P"]["qlen"]
                else:
                    for c in class_obj:
                        trf_type = trna.getTRFType(tRNA_anno_dic[tRNA_id].GetKeySitesInfor(c),
                                                   class_obj[c]["trna_start"], class_obj[c]["trna_end"])
                        mapping_ratio=(class_obj[c]["trna_end"]-class_obj[c]["trna_start"]+1)/class_obj[c]["qlen"]
                        break
            tRNA_READS_HIT.write(family_id+"\t"+tRNA_id + "\t" + read_id +"\t" +str(direction)+
                                 "\t" +str(I)+"\t" +str(I_read_start)+"\t" +str(I_read_end)+"\t" +str(I_trna_start)+"\t" +str(I_trna_end)+
                                 "\t" +str(P) + "\t" +str(P_read_start) + "\t" + str(P_read_end) + "\t" + str(P_trna_start) + "\t" + str(P_trna_end) +
                                 "\t" +str(M)+"\t" +str(M_read_start)+"\t" +str(M_read_end)+"\t" +str(M_trna_start)+"\t" +str(M_trna_end)+
                                 "\t" +str(C)+"\t" +str(C_read_start)+"\t" +str(C_read_end)+"\t" +str(C_trna_start)+"\t" +str(C_trna_end)+
                                 "\t"+str(M_5T)+"\t"+str(M_3T)+"\t"+str(M_5C)+"\t"+str(M_3C)+ "\t" + str(round(mean_exp,3)) +"\t"+
                                 str(hit_tRNAs_num) + "\t"+trf_type+"\t"+brief_mapping_infor+ "\t"+Read_type+"\t"+str(round(mapping_ratio,3))+"\t"+read_5_fragment+"\t"+read_fragment+"\t"+read_3_fragment+"\n")
    tRNA_READS_HIT.close()
    return dist_dic
