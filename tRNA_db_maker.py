import re, os, subprocess
from lib_code import share
from lib_code.trna import tRNA
import sys, getopt

#Generate four tRNA fasta files
# offset is the URL length
def CreatetRNAFastas(bed,name,ref_fasta,trna_db_dir,offset=60):
    #f0 = open(trna_db_dir+"/"+name+"_"+str(offset)+".bed",'w')  # BED file
    fasta_file = trna_db_dir + "/"+name+"_"+str(offset)+".fasta"
    print("Create tRNA FASTA file: " + fasta_file)
    f1 = open(fasta_file,'w') # fasta without introns but UTR

    tRNA_Dic = {}
    Tab = open(bed,'r')
    for line in Tab:
        if (not line.startswith("#")):
            contents = line.strip().split("\t")
            chr = str(contents[0])
            start = int(contents[1])
            end = int(contents[2])
            name = contents[3]
            strand = contents[5].strip()
            intron= ""
            if len(contents)>=10:
                intron = contents[9]
            u_start = start - offset
            u_end = end + offset

            #f0.write(chr+"\t"+str(u_start)+"\t"+str(u_end)+"\t"+name+"\t0\t"+strand+"\n")
            #full_seq= share.getFasta(ref_fasta, chr, start, end)
            seq_utr = share.getFasta(ref_fasta, chr, u_start, u_end).strip()
            #seq_sel = seq[offset:offset+end-start]
            if strand=="-":
                seq_utr = share.reverse_complement(seq_utr)

            # <BLOCKQUOTE>tRNA intron at pos 38-49<BR>GAGGCTGAAGGC(12 bp)</BLOCKQUOTE>
            s = re.search(r'pos\s(\d+)-(\d+)<BR>(\w+)\(', intron)

            s_start = 0
            s_end = 0
            s_seq = ""
            if s:
                s_start = int(s.group(1))
                s_end = int(s.group(2))
                s_seq = s.group(3)

            full_seq = seq_utr[offset:offset+end-start]
            seq = full_seq.replace(s_seq,"")
            t = tRNA()
            t.name = name
            t.chrom = chr
            t.start = start
            t.end = end
            t.strand = strand
            t.seq = seq
            t.full_seq = seq
            t.seq_utr = seq_utr
            t.utr_len = offset
            t.intron_infor = intron

            tRNA_Dic[t.name] = t

            f1.write(">"+name+"\n"+seq_utr.upper()+'\n')
            #f2.write(">" + name + "\n" + seq2.upper() + '\n')
            #f3.write(">" + name + "\n" + seq3.upper() + '\n')
            #f4.write(">" + name + "\n" + seq4.upper() + '\n')

    f1.close()
    Tab.close()
    return {"fastafile":fasta_file, "dic":tRNA_Dic}

#Run tRNAScan-SE
#tRNAscan-SE -E -f hg38_tRNA_nointron_60_structure.txt hg38_tRNA_nointron_60.fasta
def RuntRNAScan(wdir, tRNAscanSE, trna_fasta, out_file):
    try:
        if os.path.isfile(out_file):
            os.remove(out_file)
        cmd_file = wdir+"/cmd.sh"
        if os.path.isfile(tRNAscanSE):
            CMD_FILE = open(cmd_file, "w")
            print("Begin to run tRNAScan-SE, which may take several minutes ....")
            print("tRNAscanSE command : "+tRNAscanSE+" -E -f " + out_file + " " + trna_fasta + "\n")
            CMD_FILE.write(tRNAscanSE+" -E -f " + out_file + " " + trna_fasta + "\n")
        CMD_FILE.close()
        process = subprocess.Popen("bash " + cmd_file, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)
        return process.returncode
    except:
        print("Exceptions happen during running tRNAScanSE")
        return -1

# Generate annotation tab file and FASTA file for BLASTN
def parsetRNAScanFile(tRNAscan, tRNA_Dir, tabFile, faFile, no_mit_tRNA=True, no_pseudogenes=True, min_qscore=30):
    try:
        tRNA_ls = []
        FILE = open(tRNAscan, 'r')
        OUT_TAB = open(tabFile, 'w')
        OUT_FASTA = open(faFile, 'w')
        cur_trna = tRNA()
        rna_name = ""
        for line in FILE:
            if line.strip().endswith("bp"):
                if rna_name!="":
                    tRNA_ls.append(cur_trna)
                rna_name = line.strip().split(' ')[0].replace(".trna1","")
                if rna_name in tRNA_Dir:
                    cur_trna = tRNA_Dir[rna_name]
                else:
                    cur_trna = tRNA()
                    cur_trna.name = rna_name

                m = re.search(r'\((\d+)-(\d+)\)', line)
                if m:
                    cur_trna.map_start = int(m.group(1))
                    cur_trna.map_end = int(m.group(2))
                    cur_trna.map_len = abs(int(cur_trna.map_end-cur_trna.map_start)+1)

                m = re.search(r'(\d+)\s+bp$', line)
                if m:
                    cur_trna.tRNA_length = int(m.group(1))

                #print(line)
            elif line.strip().startswith("Type"):
                m = re.search(r'Type:\s(\w+)', line)
                if m:
                    cur_trna.acceptor = m.group(1)

                m = re.search(r'Anticodon:\s(\w+)', line)
                if m:
                    cur_trna.anticodon = m.group(1)

                m = re.search(r'(\d+)-(\d+)', line)
                if m:
                    cur_trna.anticodon_start_in_map = int(m.group(1))
                    cur_trna.anticodon_end_in_map = int(m.group(2))

                m = re.search(r'\((\d+)-(\d+)\)', line)
                if m:
                    cur_trna.anticodon_start = int(m.group(1))
                    cur_trna.anticodon_end = int(m.group(2))

                m = re.search(r'Score:\s+([\d\.]+)', line)
                if m:
                    cur_trna.map_scan_score = float(m.group(1))
            elif line.strip().startswith("Possible"):
                poss_str= line.strip().replace("Possible ","")
                if "intron" in poss_str:
                    cur_trna.intron_infor = poss_str
                else:
                    cur_trna.possible_type = poss_str
                cur_trna.seq = cur_trna.GetMatureSeq()
            elif line.strip().startswith("Seq:"):
                cur_trna.map_seq = line.strip().replace("Seq: ","")
            elif line.strip().startswith("Str:"):
                cur_trna.map_structure_str = line.strip().replace("Str: ", "")
                # Analysis structure string.
                # m = re.findall(r'([>]{2,}[\.]{3,}[<]{2,})', cur_trna.tRNA_structure_str)
                m = re.match(r'^[\.]{0,1}([>]{3,})', cur_trna.map_structure_str)
                if (m):
                    cur_trna.stem_for['start'] = m.regs[1][0]
                    cur_trna.stem_for['end'] = m.regs[1][1]
                    cur_trna.stem_for['str'] = cur_trna.map_seq[cur_trna.stem_for['start']:cur_trna.stem_for['end']]
                stem_len = len(cur_trna.stem_for['str'])
                pattern = "([<]{"+str(stem_len)+"})[\.]{0,1}$"
                #pattern = "([<]{7})[\.]{0,1}$"
                m1 = re.search(pattern, cur_trna.map_structure_str)
                if (m1):
                    cur_trna.stem_rev['start'] = m1.regs[1][0]
                    cur_trna.stem_rev['end'] = m1.regs[1][1]
                    cur_trna.stem_rev['str'] = cur_trna.map_seq[cur_trna.stem_rev['start']:cur_trna.stem_rev['end']]
                for m in re.finditer(r'([>]{2,}[\.]{4,}[<]{2,})', cur_trna.map_structure_str):
                    s = m.start()
                    e = m.end()
                    struct_str=m.group(0)

                    sm = re.match(r'([>]{2,})([\.]{4,})([<]{2,})', struct_str)
                    for_struct_str = sm.group(1)
                    loop_struct_str = sm.group(2)
                    rev_struct_str = sm.group(3)
                    segment = cur_trna.map_seq[s:e]
                    l_s = s+sm.regs[2][0]
                    l_e = s+sm.regs[2][1]
                    for_str = segment[sm.regs[1][0]:sm.regs[1][1]]
                    loop_str = segment[sm.regs[2][0]: sm.regs[2][1]]
                    rev_str = segment[sm.regs[3][0]: sm.regs[3][1]]
                    for_str_regs = sm.regs[1]
                    loop_str_regs = sm.regs[2]
                    rev_str_regs = sm.regs[3]
                    if s<12:
                        cur_trna.d_loop['start'] = s
                        cur_trna.d_loop['end'] = e
                        cur_trna.d_loop['l_start'] = l_s
                        cur_trna.d_loop['l_end'] = l_e
                        cur_trna.d_loop['for_str'] = for_str
                        cur_trna.d_loop['rev_str'] = rev_str
                        cur_trna.d_loop['loop_str'] = loop_str
                        cur_trna.d_loop['struct_str'] = struct_str
                    elif cur_trna.anticodon in loop_str and len(for_str)==len(rev_str):
                        cur_trna.a_loop['start'] = s
                        cur_trna.a_loop['end'] = e
                        cur_trna.a_loop['l_start'] = l_s
                        cur_trna.a_loop['l_end'] = l_e
                        cur_trna.a_loop['for_str'] = for_str
                        cur_trna.a_loop['rev_str'] = rev_str
                        cur_trna.a_loop['loop_str'] = loop_str
                        cur_trna.a_loop['struct_str'] = struct_str
                    elif len(rev_str)>len(for_str) or s>45:
                        e=s+loop_str_regs[1]+len(for_str)
                        struct_str = for_struct_str+loop_struct_str+'<'*len(for_str)
                        cur_trna.t_loop['start'] = s
                        cur_trna.t_loop['end'] = e
                        cur_trna.t_loop['l_start'] = l_s
                        cur_trna.t_loop['l_end'] = l_e
                        cur_trna.t_loop['for_str'] = for_str
                        cur_trna.t_loop['rev_str'] = rev_str[0:len(for_str)]
                        cur_trna.t_loop['struct_str'] = struct_str[0:len(for_str)]

        # OUT_TAB.write(tRNA_ls[0].GetTabTitle()+"\n")
        # Sort list
        sorted_RNA_ls = sorted(tRNA_ls, key=lambda x: x.seq+"_"+x.name, reverse=False)
        # Get family id
        seq = ""
        cur_family=""
        OUT_TAB.write(tRNA().GetTabTitle() + "\n")
        for t in sorted_RNA_ls:
            if seq == "":
                t.family = "TF_"+t.name
                cur_family = t.family
                seq = t.seq
            if seq != t.seq:
                t.family = "TF_"+t.name
                cur_family = t.family
                seq = t.seq
            else:
                t.family = cur_family
            OUT_TAB.write(t.GetTabStr() + "\n")
            if t.IsQualifiedtRNA(no_mit_tRNA=no_mit_tRNA, no_pseudogenes=no_pseudogenes,min_qscore=min_qscore):
                if t.intron_infor=="":
                    # Because the is no intron in tRNA
                    OUT_FASTA.write(">PI::" + t.name + "\n")  # Pre tRNA and Pre Intron tRNA
                    OUT_FASTA.write(t.seq_utr + "\n")
                else:
                    OUT_FASTA.write(">I::" + t.name + "\n")  # Pre tRNA and Pre Intron tRNA
                    OUT_FASTA.write(t.seq_utr + "\n")
                    OUT_FASTA.write(">P::" + t.name + "\n")  # Pre tRNA and Pre Intron tRNA
                    OUT_FASTA.write(t.GetPreMatureNoIntronSeq() + "\n")
                OUT_FASTA.write(">C::"+t.name + "\n") # Mature CCA tRNA
                OUT_FASTA.write(t.GetMatureSeq() + "CCA\n")
                OUT_FASTA.write(">M::" + t.name + "\n") # Mature tRNA
                OUT_FASTA.write(t.GetMatureSeq() + "\n")
        OUT_TAB.close()
        OUT_FASTA.close()
        FILE.close()
        return 0
    except:
        print("An exception occurred")
        return -1

# Function for generate tRNA datanbase for following analysis
# Filter and extract structure information for tRNAs
# Input:
#   1) tRNA_bed: bed file for tRNA which can be downloaded from UCSC table browser
#   2) db_name: The name of you tRNA databse such as tRNA_hg38_60 recommend patterm tRNA_[genome]_offset
#   3) tRNAscanSE: The absolute path of tRNAScan-SE
#   4) offset is the length of the 5' abd 3' UTRs
# Options:
#   1) no_mit <1/else> : 1 means removing all nucleic mitochondrial tRNA genes'
#   2) no_pseu <1/else> : 1 means removing all pseudogene tRNA genes'
#   3) minq <number> : The minimum of quality scores of tRNAs provided by tRNAScan-SE
# Output: (these files will be created at the same folder of the tRNA_bed file):
#   1) A FASTA file of tRNAs : trna_db_dir+"/"+db_name+"_" + offset +".fa"
#   2) A bed file with all tRNA annotations: trna_db_dir + "/" + db_name + "_" + offset + ".bed"


def tRNA_DB_Preparing(db_name, tRNA_bed, ref_fasta, tRNAscanSE, offset=60, no_mit_tRNA=True, no_pseudogenes=True, min_qscore=30):
    print("######---------- tRNA Explorer : tRNA database preparing  -----------------#")
    print("###### Begining tRNA database preparing ... ")
    trna_db_dir = os.path.dirname(os.path.abspath(tRNA_bed))
    print("###### Create FASTA file for tRNAscan SE scanning ")
    result = CreatetRNAFastas(tRNA_bed, db_name, ref_fasta, trna_db_dir, offset)
    trna_utr_fasta= result['fastafile']
    trna_scan_out = trna_db_dir+"/"+db_name+"_scan_"+str(offset)+".out"
    trna_detail_bed = trna_db_dir + "/" + db_name + "_" + str(offset) + ".bed"
    output_code = RuntRNAScan(trna_db_dir, tRNAscanSE, trna_utr_fasta, trna_scan_out)
    if output_code!=-1:
        print("###### tRNAScan-SE running is down, Generate bed file with annotations ...")
        fastfile = trna_db_dir+"/"+db_name+"_" + str(offset) +".fa"
        print("######  Parsing tRNAScan-SE output file")
        if parsetRNAScanFile(trna_scan_out,
                             result["dic"],
                             trna_detail_bed,
                             fastfile,
                             no_mit_tRNA=no_mit_tRNA,
                             no_pseudogenes=no_pseudogenes,
                             min_qscore=min_qscore)==0:
            print("###### tRNA database named \'"+db_name+"\' was generated successfully #######\n")
            print("###### Output 1 :FASTA of mixed tRNAs fragment :"+fastfile)
            print("###### Output 2 :Annotated  BED FILE:" + trna_detail_bed)
            print("######----------- Finished -------------------")
            return [fastfile,trna_detail_bed]
        else:
            print("Something wrong when parsing tRNAScan out file: "+trna_scan_out)
    else:
        print("Something wrong during running tRNA ScanSE with FASTA file :"+trna_utr_fasta)
    return []


def main(argv):
    # Default settings
    db_name= "hg38_tRNA"
    tRNA_bed= "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer_Old/test_data/ChIP/ucsc_hg38_tRNA.bed"
    ref_fasta ="/Users/hqyone/Desktop/hg38/hg38.fa"
    tRNAscanSE="/usr/local/bin/tRNAscan-SE"
    offset=60
    no_mit_tRNA = True
    no_pseudogenes = True
    min_qscore = 30
    try:
        opts, args = getopt.getopt(argv, "n:b:r:s:u:", ["name=","bed=","ref=","scan=","offset=","no_mit=", "no_pseu=","minq="])
    except getopt.GetoptError:
        print('# Usage: python tRNA_db_maker.py -n <name> -b <bed> -r <ref> -s <tRNAScanSE> -u <offset> --no_mit <1/else> --no_pseu <1/else> --minq <number>')
        print('# -n <name> : The pref name of database')
        print('# -b <bed> : Absolute path of bed file for tRNAs')
        print('# -r <ref> : Absolute path of Genome FASTA file coordinating with bed file')
        print('# -s <tRNAScanSE> : Absolute path of tRNAscan SE program')
        print('# -u <offset> : The length of UTRs')
        print('# --no_mit <1/else> : 1 means removing all nucleic mitochondrial tRNA genes')
        print('# --no_pseu <1/else> : 1 means removing all pseudogene tRNA genes')
        print('# --minq <number> : The minimum of quality scores of tRNAs, default 30')
        print('# Output 1: annotation bed file of tRNAs <Directory of tRNA bed file>/<name>.bed')
        print('# Output 2: FASTA file of hybrid tRNA library <Directory of tRNA bed file>/<name>.fa')
        sys.exit(2)
    print('The configs are as following"')
    for opt, arg in opts:
        if opt == '-h':
            print('# Usage: python tRNA_db_maker.py -n <name> -b <bed> -r <ref> -s <tRNAScanSE> -o <offset> --no_mit <1/else> --no_pseu <1/else> --minq <number>')
            print('# -n <name> : The name of database')
            print('# -b <bed> : Absolute path of bed file for tRNAs')
            print('# -r <ref> : Absolute path of Genome FASTA file coordinating with bed file')
            print('# -s <tRNAScanSE> : Absolute path of tRNAscan SE program')
            print('# -u <offset> : The length of UTRs')
            print('# --no_mit <1/else> : 1 means removing all nucleic mitochondrial tRNA genes')
            print('# --no_pseu <1/else> : 1 means removing all pseudogene tRNA genes')
            print('# --minq <number> : The minimum of quality scores of tRNAs')
            print('# Output 1: annotation bed file of tRNAs <Directory of tRNA bed file>/<name>.bed')
            print('# Output 2: FASTA file of hybrid tRNA library <Directory of tRNA bed file>/<name>.fa')
            sys.exit(0)
        elif opt in ("-n", "--name"):
            db_name = arg
            print("database name: " + db_name)
        elif opt in ("-b", "--bed"):
            tRNA_bed = arg
            print("tRNA bed file: " + tRNA_bed)
        elif opt in ("-r", "--ref"):
            ref_fasta = arg
            print("genome references: " + ref_fasta)
        elif opt in ("-s", "--scan"):
            tRNAscanSE = arg
            print("tRNAscanSE: " + tRNAscanSE)
        elif opt in ("-u", "--offset"):
            offset = int(arg)
            print("offset: " + str(offset))
        elif opt in ("--no_mit"):
            if arg.strip() != "1":
                no_mit_tRNA = False
        elif opt in ("--no_pseu"):
            if arg.strip() != "1":
                no_pseudogenes = False
        elif opt in ("--minq"):
            min_qscore = int(arg)
    if os.path.isfile(tRNAscanSE):
        tRNA_DB_Preparing(db_name, tRNA_bed, ref_fasta, tRNAscanSE, offset=offset, no_mit_tRNA=no_mit_tRNA, no_pseudogenes=no_pseudogenes, min_qscore=min_qscore)
        exit(0)
    else:
        print("The tRNAscanSE execute file can not be found. Abort!")
        exit(-1)

if __name__ == "__main__":
    main(sys.argv[1:])

# Local Test
# bed = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer_Old/test_data/ChIP/ucsc_hg38_tRNA.bed"
# name = "hg38_tRNA"
# ref_fasta="/Users/hqyone/Desktop/hg38/hg38.fa"
# trna_db_dir = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools"
# tRNAscanSE="/usr/local/bin/tRNAscan-SE"
# print(tRNA_DB_Preparing(name, bed, ref_fasta, tRNAscanSE))

# Server Test
# bed = "/home/hqyone/python_code/testrepo/szf/ucsc_hg38_tRNA.bed"
# name = "hg38_tRNA"
# ref_fasta="/home/hqyone/python_code/testrepo/szf/hg38.fa"
# tRNAscanSE="/usr/local/bin/tRNAscan-SE"
# print(tRNA_DB_Preparing(bed, name, ref_fasta, tRNAscanSE))
# exit(0)
