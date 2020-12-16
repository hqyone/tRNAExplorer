# coding=utf-8

#  tRNAExplorer v.1.0 - tRNA profiles and regulation networks
#  Copyright (C) 2020  Dr. Quanyuan He
#  School of Medicine, Hunan Normal University
#  Email: hqyone@hotmail.com
#  Freely distributed under the GNU General Public License (GPLv3)
#

import os
from lib_code.trna import tRNA
from lib_code import blast_tools, make_report, share
from lib_code.cmd_tools import Trimmomatic
import subprocess
import sys
import getopt
import time
import pathlib
import shutil


version = 1.0
# Function for running BLASTN and classify/quantify tRFs for samples
# Create multiple matrix and tab file to summarize the results


def rseq_blastn_pipeline(proj_name,
                         trna_fa,
                         trna_anno_bed,
                         sample_tsv,
                         fastq_dir,
                         out_dir,
                         url_len=60,
                         no_indexing=False,
                         no_alignment=False,
                         t_path="",
                         t_adapter="",
                         t_phred=33,
                         t_leading=3,
                         t_tailing=3,
                         t_slidingwindow="4:15",
                         t_minlen=18,
                         t_threads=2,
                         min_read_qscore=20,
                         min_reads_count=50,
                         blastn="",
                         mkdb="",
                         blastn_e_cutoff=0.001,
                         blastn_max_mismatch=2,
                         blastn_max_hit_num=40,
                         blastn_pident=98,
                         min_trf_len=18,
                         trim_seq=1,
                         ):
    try:
        print("Begin the tRNA processing pipeline with fastq....")
        # Read sample_tsv
        # Load information about tRNAs        if
        if not os.path.isfile(trna_anno_bed):
            print("The tRNA annotation bed file :" +
                  trna_anno_bed+" is not exist. Abort!")
            return -1
        tRNAFile = open(trna_anno_bed, 'r')
        tRNA_dic = {}
        for line in tRNAFile:
            if not line.startswith("#"):
                t = tRNA()
                t.LoadStr(line.strip())
                tRNA_dic[t.name] = t
        tRNAFile.close()

        # Load sample information
        sample_dic = {}
        if not os.path.isfile(sample_tsv):
            print("The sample file :"+sample_tsv+" is not exist. Abort!")
            return -1
        # Copy sample tsv for visualization
        cp_sample_tsv = out_dir+"/samples"
        shutil.copy(sample_tsv, cp_sample_tsv)
        SAMPLES = open(sample_tsv, 'r')
        for line in SAMPLES:
            if not line.startswith("#"):
                contents = line.strip().split("\t")
                if len(contents) > 1:
                    sample_dic[contents[0]] = {
                        "des": contents[1], "adapters": ""}
                if len(contents) > 2:
                    sample_dic[contents[0]]['adapters'] = contents[2]

        SAMPLES.close()

        loginfor = {}
        # Run fastq file one by one
        ext = ".fastq"
        fastq_files = share.getExtFileList(fastq_dir, ".fastq")
        if len(fastq_files) == 0:
            fastq_files = share.getExtFileList(fastq_dir, ".fq")
            ext = ".fq"

        # Create visual_config file
        # trna vs sample: pileup max matrix
        visual_config = out_dir + "/visual_config.tsv"
        VISUAL_CONFIG = open(visual_config, 'w')
        VISUAL_CONFIG.write("out_dir="+out_dir+"\n")
        VISUAL_CONFIG.write("sample_tsv=" + sample_tsv+"\n")
        VISUAL_CONFIG.write("trna_anno_bed=" + trna_anno_bed)
        VISUAL_CONFIG.close()
        if not no_alignment:
            for f in fastq_files:
                s_id = os.path.basename(f).replace(
                    ".fastq", "").replace(".fq", "")
                fastq_dir = os.path.dirname(os.path.abspath(f))
                if s_id in sample_dic:
                    des = sample_dic[s_id]["des"]
                    adapters = sample_dic[s_id]["adapters"].strip().split(',')
                    f_adapter = str(adapters[0])
                    r_adapter = ""
                    if len(adapters) > 1:
                        r_adapter = adapters[1]

                    start_time = time.time()
                    end_time = time.time()
                    processing_time = 0
                    raw_fastq = fastq_dir + "/" + s_id + ext
                    cmd_bash = fastq_dir + "/" + s_id + "_blast.sh"
                    CMD_FILE = open(cmd_bash, "w")
                    trimmed_fastq = raw_fastq

                    trim_start_time = time.time()
                    if trim_seq != 0:
                        # Trimed fastq
                        trimmed_fastq = fastq_dir + "/" + s_id + "_trimmed"+ext
                        adapter_fasta = t_adapter
                        if not os.path.isfile(adapter_fasta):
                            adapter_fasta = ""
                        T = Trimmomatic(t_path)
                        trimmomatics_cmd = T.TrimSE(raw_fastq, trimmed_fastq, adapter_fa=adapter_fasta, phred=t_phred, LEADING=t_leading,
                                                    TRAILING=t_tailing, SLIDINGWINDOW=t_slidingwindow, MINLEN=t_minlen, threads=t_threads)

                        CMD_FILE.write(trimmomatics_cmd + "\n")
                        CMD_FILE.close()
                        process = subprocess.Popen(
                            "bash " + cmd_bash, shell=True, stdout=subprocess.PIPE)
                        process.wait()
                    trim_end_time = time.time()

                    filter_start_time = time.time()
                    # Removed redundant read Filter and get the read number file
                    filtered_fasta = fastq_dir + "/" + s_id + "_filtered.fa"
                    num_dic_txt = fastq_dir + "/" + s_id + "_num_dic.txt"
                    static_infor = {}
                    if not no_indexing:
                        static_infor = share.filterFastQ2FastA(trimmed_fastq, filtered_fasta, num_dic_txt,
                                                               f_patterm=f_adapter, r_patterm=r_adapter, qcutoff=min_read_qscore, num_cutoff=min_reads_count)
                        print(static_infor)
                    else:
                        print("Skip filtering and indexing FASTQ flie : "+f)
                    loginfor[s_id] = static_infor
                    # delete trmmed fastq to save space
                    if trim_seq != 0 and os.path.exists(trimmed_fastq):
                        os.remove(trimmed_fastq)
                    filter_end_time = time.time()

                    # Do BLASTN
                    blastn_start_time = time.time()
                    blast_out_file = blast_tools.RunBLASTN(
                        blastn, mkdb, s_id, trna_fa, filtered_fasta, out_dir, eval=blastn_e_cutoff, hit_number=blastn_max_hit_num)
                    #blast_out_file = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/output/SRR1836126_1_tRNA_blast_out.tab"
                    # Analysis BLASTN result
                    if blast_out_file != "":
                        print(
                            "BLASTN Successfully, processing BLAST out tab file: " + blast_out_file)
                        tRNA_reads_count_file = out_dir + "/" + s_id + "_" + proj_name + "_count.tab"
                        tRNA_reads_hit_file = out_dir + "/" + s_id + "_" + proj_name + "_hit.tab"
                        blast_tools.AnalysisBlastOut2(blast_out_file, num_dic_txt, tRNA_dic,
                                                      tRNA_reads_count_file, tRNA_reads_hit_file, url_len, max_mismatch=blastn_max_mismatch, pident=blastn_pident)
                        print("Hit file :"+tRNA_reads_hit_file)
                        end_time = time.time()
                        processing_time = end_time - start_time
                        print("Processing time :" +
                              str(round(processing_time, 3))+" Secs")
                    else:
                        print("Something wrong while blastn " + s_id)
                    blastn_end_time = time.time()

                    loginfor[s_id]["start_time"] = time.asctime(
                        time.localtime(start_time))
                    loginfor[s_id]["end_time"] = time.asctime(
                        time.localtime(end_time))
                    loginfor[s_id]["processing_time"] = str(
                        int(processing_time))
                    loginfor[s_id]["trim_time"] = str(
                        int(trim_end_time-trim_start_time))
                    loginfor[s_id]["filter_time"] = str(
                        int(filter_end_time-filter_start_time))
                    loginfor[s_id]["blastn_time"] = str(
                        int(blastn_end_time - blastn_start_time))
        else:
            print("Skip indexing and alignment for FASTQ files")

        static_dic = make_report.getTrfReportFile(
            proj_name, out_dir, tRNA_dic, sample_dic, out_dir)
        if no_indexing or no_alignment:
            print("No indexing and alignment so skip update logfile.")
        else:
            # Write loginfor into a file
            for s_id in static_dic:
                if s_id in loginfor:
                    for key in static_dic[s_id]:
                        loginfor[s_id][key] = static_dic[s_id][key]
            logfile = out_dir + "/static.log"
            LOG = open(logfile, 'w')
            # LOG.write("#FASTQ_STATISTICS\n")
            title = ""
            key_ls = ['total_num', 'removed_num', 'survived_num',
                      'non_redundent_num', 'start_time', 'end_time', 'processing_time',
                      'trim_time', 'filter_time', 'blastn_time', 'A', 'B', 'C', 'D', 'E',
                      'F', 'G', 'H', 'I', 'total', 'intro_cl_ratio', 'u5_cl_ratio', 'u3_cl_ratio', 'cca_add_ratio']
            for s_id in loginfor:
                statis_obj = loginfor[s_id]
                if title == "":
                    title = "#SampleID"+"\tDescription\t" + "\t".join(key_ls)
                    LOG.write(title + "\n")

                sample_des = s_id
                if s_id in sample_dic:
                    sample_des = sample_dic[s_id]["des"]
                line = s_id+"\t"+sample_des
                for key in key_ls:
                    if key in statis_obj:
                        line += "\t"+str(statis_obj[key])
                    else:
                        line += "\t" + ""
                LOG.write(line+"\n")
            LOG.close()
        return 0
    except Exception as inst:
        print(type(inst))  # the exception instance
        print(inst)
        print('Some things wrong during pipeline running!')
        return -1
        # sys.exit(2)


class Config:
    def __init__(self):
        # Load The Config File
        self.config = {}

    def loadConfig(self, config_file):
        if os.path.isfile(config_file):
            with open(config_file, 'r') as FILE:
                for line in FILE:
                    if line.startswith("#"):
                        continue
                    contents = line.split("=")
                    if len(contents) == 2:
                        KEY = contents[0].strip().strip('"')
                        VAL = contents[1].strip()
                        if VAL.startswith('"'):
                            VAL = VAL.strip('"')
                        elif VAL.lower().strip() == 'true':
                            VAL = True
                        elif VAL.lower().strip() == 'false':
                            VAL = False
                        else:
                            try:
                                if "." in VAL:
                                    VAL = float(VAL)
                                else:
                                    VAL = int(VAL)
                            except Exception as e:
                                print(e)
                                print("The option parser fails " + line)
                                print(
                                    "Please check the config file. The String options should be enbraced by \"")
                                return False
                        self.config[KEY] = VAL
        return True


def rseq_blastn_pipeline2(config):
    try:
        rseq_blastn_pipeline(
            config["proj_name"],
            config["trna_fa"],
            config["trna_anno_file"],
            config["sample_tsv"],
            config["fastq_dir"],
            config["out_dir"],
            config["url_len"],
            config["no_indexing"],
            config['no_alignment'],
            t_path=config["t_path"],
            t_adapter=config["t_adapter"],
            t_phred=config["t_phred"],
            t_leading=config["t_leading"],
            t_tailing=config["t_tailing"],
            t_slidingwindow=config["t_slidingwindow"],
            t_minlen=config["t_minlen"],
            t_threads=config["t_threads"],
            min_read_qscore=config["min_read_qscore"],
            min_reads_count=config["min_reads_count"],
            blastn=config["blastn"],
            mkdb=config["mkdb"],
            blastn_e_cutoff=config["blastn_e_cutoff"],
            blastn_max_mismatch=config["blastn_max_mismatch"],
            blastn_max_hit_num=config["blastn_max_hit_num"],
            blastn_pident=config["blastn_pident"],
            min_trf_len=config["min_trf_len"],
            trim_seq=config["t_do"])
        print('Project ['+config["proj_name"]+'] running completely!')
        return 0
        # sys.exit(0)
    except Exception as inst:
        print(type(inst))  # the exception instance
        print(inst)
        print('Some things wrong during pipeline running!')
        return -1
        # sys.exit(2)

# Main function


def main(argv):
    config_key_ls = [
        "proj_name",
        "trna_fa",
        "trna_anno_file",
        "sample_tsv",
        "fastq_dir",
        "out_dir",
        "url_len",
        "no_indexing",
        "no_alignment",
        "t_do",
        "t_path",
        "t_adapter",
        "t_phred",
        "t_leading",
        "t_tailing",
        "t_slidingwindow",
        "t_minlen",
        "t_threads",
        "min_read_qscore",
        "min_read_length",
        "blastn",
        "mkdb",
        "blastn_e_cutoff",
        "blastn_pident",
        "blastn_max_mismatch",
        "blastn_max_hit_num",
        "min_trf_len"
    ]
    currentDirectory = os.getcwd()
    wdir = pathlib.Path(__file__).parent.absolute()
    config = {
        # Default settings
        #########################################################
        # Project Settings
        #########################################################
        "proj_name": "test",
        "trna_fa": str(wdir)+"/test/trna_db/hg38_tRNA_60.fa",
        "trna_anno_file": str(wdir)+"/test/trna_db/hg38_tRNA_60.bed",
        "sample_tsv": str(wdir)+"/test/samples",
        "fastq_dir": str(wdir)+"/test/fastq",
        "out_dir": str(wdir)+"/test/output",
        "url_len": 60,
        "no_indexing": False,
        "no_alignment": False,

        #########################################################
        # Trimmomatic
        #########################################################
        "t_do": 0,
        "t_adapter": "",
        "t_path": "/Users/hqyone/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar",
        "t_phred": 33,
        "t_leading": 3,
        "t_tailing": 3,
        "t_slidingwindow": "4:15",
        "t_minlen": 18,
        "t_threads": 2,

        #########################################################
        # Read Filtering
        #########################################################
        "min_read_qscore": 5,
        "min_read_length": 18,
        "min_reads_count": 10,

        #########################################################
        # BLASTN settings
        #########################################################
        "blastn": "/Users/hqyone/Downloads/ncbi-blast-2.10.0+/bin/blastn",
        "mkdb": "/Users/hqyone/Downloads/ncbi-blast-2.10.0+/bin/makeblastdb",
        "blastn_e_cutoff": 0.001,
        "blastn_max_mismatch": 2,
        "blastn_max_hit_num": 40,
        "blastn_pident": 98,

        #########################################################
        # TRF analyiss settinss
        #########################################################
        "min_trf_len": 18
    }

    c = Config()
    config_file = ""
    try:
        opts, args = getopt.getopt(sys.argv[1:], "c:n:f:a:s:i:o:h:v", [
                                   "config", 'no-indexing', 'no-alignment'])
    except getopt.GetoptError as err:
        print(err)
        print('# Usage: python tRNAExplorer.py -c <configfile> ')
        print('# -c config file> : The absolute path of config file')
        print('# -n <proj_name> : The name of project')
        print('# -f <trna_fa> : Absolute path of FASTA file for tRNAs which was created by tRNA_db_maker')
        print('# -a <trna_anno_file> : Absolute path of bed file for tRNA annotations which was created by tRNA_db_maker')
        print('# -s <sample tsv> : Absolute path of sample information')
        print('# -i <fastq_dir> : The directory storing fastq files (Input Directory)')
        print('# -h : Show the help information')
        print('# -o <out_dir> : The directory of output files (Out Directory)')
        print('# --no-indexing : Skip sequence indexing step')
        print('# --no-alignmnent : Skip sequence indexing and alignment step')
        print('# Output 1: <out_dir>/static.log , Reads numbers, processing time for each sample')
        print('# Output 2: <out_dir>/trf_sample_matrix.tsv , Read number matrix of tRFs across tRFs and samples')
        print('# Output 3: <out_dir>/trna_sample_readcount_matrix.tsv , Read number matrix across tRNAs and samples')
        print('# Output 4: <out_dir>/trna_sample_pileup_matrix.tsv , Pileup depth matrix across tRNAs and samples')
        print('# Output 5: <out_dir>/trna_trftype_matrix.tsv , The read number matrix across samples/tRNAs and tRF types')
        print('# Output 6: <out_dir>/sample_trftype_matrix.tsv , The read number matrix across samples and tRF types')
        print('# Output 7: <out_dir>/cleavage_sites.tsv , Cleavage sites information for tRNAs in different samples')
        print('# Output 8: <out_dir>/profiles.tsv , Pileup information for tRNAs in different samples')
        print(
            '# Output 9: <out_dir>/variants.tsv , A tsv file to summarize mismatched sites in each gene for all samples .')
        print(
            '# Output 10: <out_dir>/visual_config.tsv , A tsv file including paths of files required for the visualization module')
        sys.exit(2)
    print('###########################################################################################')
    print('############                      tRNAExplorer ['+str(
        version)+"]                          #############")
    print('############      A tool for tRNA and tRF quantification and annotation       #############')
    print('############      Author:  Quanyuan He Ph.D  Contact: hqyone@hotmail.com      #############')
    print('############      Institution :  School of Medicine, Hunan Normal University  #############')
    print('############  Freely distributed under the GNU General Public License (GPLv3) #############')
    print('############                             2020/06/15                           #############')
    print('###########################################################################################')
    for opt, arg in opts:
        if opt == '-h':
            print('# Usage: python tRNAExplorer.py -c <configfile> ')
            print('# -c <path to config file> : The absolute path of config file')
            print('# -h : Show the help information')
            print(
                '# -f <trna_fa> : Absolute path of FASTA file for tRNAs which was created by tRNA_db_maker')
            print(
                '# -a <trna_anno_file> : Absolute path of bed file for tRNA annotations which was created by tRNA_db_maker')
            print('# -s <sample tsv> : Absolute path of sample information')
            print(
                '# -i <fastq_dir> : The directory storing fastq files (Input Directory)')
            print('# -h : Show the help information')
            print('# -o <out_dir> : The directory of output files (Out Directory)')
            print('# --no-indexing : Skip sequence indexing step')
            print('# --no-alignmnent : Skip sequence indexing and alignment step')
            print(
                '# Output 1: <out_dir>/static.log , Reads numbers, processing time for each sample')
            print('# Output 2: <out_dir>/trf_sample_matrix.tsv , Read number matrix of tRFs across tRFs and samples')
            print(
                '# Output 3: <out_dir>/trna_sample_readcount_matrix.tsv , Read number matrix across tRNAs and samples')
            print('# Output 4: <out_dir>/trna_sample_pileup_matrix.tsv , Pileup depth matrix across tRNAs and samples')
            print(
                '# Output 5: <out_dir>/trna_trftype_matrix.tsv , The read number matrix across samples/tRNAs and tRF types')
            print(
                '# Output 6: <out_dir>/sample_trftype_matrix.tsv , The read number matrix across samples and tRF types')
            print(
                '# Output 7: <out_dir>/cleavage_sites.tsv , Cleavage sites information for tRNAs in different samples')
            print(
                '# Output 8: <out_dir>/profile.tsv , Pileup information for tRNAs in different samples')
            print(
                '# Output 9: <out_dir>/variants.tsv , A tsv file to summarize mismatched sites in each gene for all samples .')
            print('# Output 10: <out_dir>/visual_config.tsv , A tsv file including paths of files required for the visualization module')
            sys.exit(0)
        elif opt == '-n':
            config["proj_name"] = arg
        elif opt == '-f':
            config["trna_fa"] = arg
        elif opt == '-a':
            config["trna_anno_file"] = arg
        elif opt == '-s':
            config["sample_tsv"] = arg
        elif opt == '-i':
            config["fastq_dir"] = arg
        elif opt == '-o':
            config["out_dir"] = arg
        elif opt == '-h':
            config["url_len"] = int(arg)
        elif opt == '--no-indexing':
            config["no_indexing"] = True
        elif opt == '--no-alignment':
            config["no_alignment"] = True
        elif opt in ("-c"):
            print("Load "+arg)
            config_file = arg
    print('##  ------------------ Initialization ....  -----------------')
    init_file = str(wdir)+"/init"
    if os.path.isfile(init_file):
        ic = Config()
        if ic.loadConfig(init_file):
            for key in ic.config:
                if key in config:
                    config[key] = ic.config[key]
                    print(key + ":"+ic.config[key])
        print("Loading init file successfully! ")
    else:
        print("No 'init' file was given, Abort!")
    print('##  ------------------ The loading config ....  -----------------')
    if config_file != "":
        input_cfile = config_file
        if not os.path.isfile(config_file):
            config_file = currentDirectory+"/"+config_file
        if os.path.isfile(config_file):
            if c.loadConfig(config_file):
                for key in c.config:
                    if key in config:
                        config[key] = c.config[key]
            print("Load config file: " + config_file)
        else:
            print("Can't find config file: " + input_cfile)
            print("Use default settings")
    else:
        print("No Config file was given, run for test data!")
    print('##  ------------------ The configs are as following ....  -----------------')
    for key in config_key_ls:
        if key in config:
            print("## "+key+"="+str(config[key])+"\n")
    print('##  ------------------            End Settings             ----------------')
    print('###########################################################################################')
    print('## Setting testing .... ')
    if not os.path.isfile(config["blastn"]):
        print("The BLASTN path "+config["blastn"]+" is not valid. Abort!!!\n")
        print("Please modify the file <wdir>/init \n")
        exit(-1)
    if not os.path.isfile(config["mkdb"]):
        print("The BLASTN mkdb path " +
              config["mkdb"]+" is not valid. Abort!!!\n")
        print("Please modify the file <wdir>/init \n")
        exit(-1)
    print('## Run pipeline ..... ')
    print(rseq_blastn_pipeline2(config))
    exit(0)
    # Run pipeline


if __name__ == "__main__":
    main(sys.argv[1:])
# Test
# db_dir = "/home/hqyone/python_code/testrepo/szf"
# tRNA_bed = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/ChIP/hg38_tRNA_60.bed"
# genome_fa = db_dir + "/hg38.fa"
# # bam_dir = "/home/hqyone/python_code/fastq/trna_method/bam"
# bam_dir = "/media/server2/e88c6737-2e83-4500-a2f4-9f4b28c1f3d2/sunzf/bam"
# # blast_out_dir="/home/hqyone/python_code/fastq/trna_method/out"
# blast_out_dir = "/media/server2/e88c6737-2e83-4500-a2f4-9f4b28c1f3d2/sunzf/tRNA_analysis_version_2/blast_out"
# #out_dir = blast_out_dir
# tRNA_anno_file = "/home/hqyone/python_code/testrepo/new_tools/tRNA_anticodon_position.txt"
# tRNA_class_file = "/home/hqyone/python_code/testrepo/new_tools/tRNA_families.txt"

# proj_name= "test"
# config_file="/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorerII/config.txt"
#
# c = Config()
# c.loadConfig(config_file)
# rseq_blastn_pipeline2(c)
# #rseq_fastq_pipeline(proj_name, genome_fa,tRNA_bed, sample_tsv, fastq_dir, out_dir, doBLAST=True)
# print("tRNA pipeline finish successfully!")
# sys.exit(0)
