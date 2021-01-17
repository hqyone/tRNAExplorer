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
def singleton(class_):
    instances = {}
    def getinstance(*args, **kwargs):
        if class_ not in instances:
            instances[class_] = class_(*args, **kwargs)
        return instances[class_]
    return getinstance

@singleton
class Config():
    def __init__(self, wdir):
        # Default configuration list
        self.config_key_ls = [
            "proj_name","trna_fa","trna_anno_bed",
            "sample_tsv","fastq_dir","out_dir",
            "url_len","no_indexing","no_alignment",
            "t_do","t_path","t_adapter","t_phred",
            "t_leading","t_tailing","t_slidingwindow",
            "t_minlen","t_threads","min_read_qscore",
            "min_read_length","blastn",
            "mkdb","blastn_e_cutoff","blastn_pident",
            "blastn_max_mismatch","blastn_max_hit_num",
            "min_trf_len"
        ]
        # Default configuration
        self.config = {
            # Default settings
            #########################################################
            # Project Settings
            #########################################################
            "proj_name": "test",
            "trna_fa": str(wdir)+"/test/trna_db/hg38_tRNA_60.fa",
            "trna_anno_bed": str(wdir)+"/test/trna_db/hg38_tRNA_60.bed",
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

    # Load The Config File  
    def loadConfigFile(self, config_file):
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
    
    # Print configure settings
    def printConfig(self):
        print('##  ------------------ The configs are as following ....  -----------------')
        for key in self.config_key_ls:
            if key in self.config:
                print("## "+key+"="+str(self.config[key])+"\n")
        print('##  ------------------            End Settings             ----------------')
        print('###########################################################################################')
        print('## Setting testing .... ')
        if not os.path.isfile(self.config["blastn"]):
            print("The BLASTN path "+self.config["blastn"]+" is not valid. Abort!!!\n")
            print("Please modify the file <wdir>/init \n")
            exit(-1)
        if not os.path.isfile(self.config["mkdb"]):
            print("The BLASTN mkdb path " +
                self.config["mkdb"]+" is not valid. Abort!!!\n")
            print("Please modify the file <wdir>/init \n")
            exit(-1)
    
    # Create the visual setting file
    def createVisualSettingFile(self):
        visual_config = self.config["out_dir"] + "/visual_config.tsv"
        with open(visual_config, 'w') as VISUAL_CONFIG:
            VISUAL_CONFIG.write("out_dir="+self.config["out_dir"]+"\n")
            VISUAL_CONFIG.write("sample_tsv=" + self.config["sample_tsv"]+"\n")
            VISUAL_CONFIG.write("trna_anno_bed=" + self.config["trna_anno_bed"])

class pipeline():
    def __init__(self, config):
      self.config = config
    
    def run(self):
        pass

# Function for running BLASTN and classify/quantify tRFs for samples
# Create multiple matrix and tab file to summarize the results
class rseqBlastnPipeline(pipeline):
    def __init__(self, config):
      self.config = config
    
    def run(self):
        cfg = self.config
        loginfor = {}
        try:
            proj_name = cfg["proj_name"]
            out_dir = cfg["out_dir"]

            print("Begin the tRNA processing pipeline with fastq....")
            # Read sample_tsv
            trna_anno_bed = cfg["trna_anno_bed"]
            tRNA_dic = self.load_tRNA_database(trna_anno_bed)

            # Load sample information
            sample_dic = self.loadSampleInfor(cfg["sample_tsv"], cfg["out_dir"])

            # Get list of fastq files
            ext,fastq_files =  self.getFastqList(cfg["fastq_dir"])

            # Process FASTQ file one by one
            if not cfg["no_alignment"]:
                for f in fastq_files:
                    s_id = os.path.basename(f).replace(
                        ".fastq", "").replace(".fq", "")
                    fastq_dir = os.path.dirname(os.path.abspath(f))
                    if s_id in sample_dic:
                        start_time = time.time()
                        
                        # Trimed fastq
                        trim_time, trimmed_fastq = self.runTrimmomatics(s_id, ext)

                        # Transform FASTQ to FASTQ file and filtering
                        filter_time, static_infor, filtered_fasta, num_dic_txt = \
                            self.fastq2Fasta(sample_dic, s_id, fastq_dir, trimmed_fastq)
                        
                        # Mapping with BLASTN
                        blast_time, tRNA_reads_count_file, tRNA_reads_hit_file = \
                            self.doMapping(s_id, num_dic_txt, tRNA_dic, filtered_fasta)

                        end_time = time.time()
                        processing_time = end_time - start_time
                        print("Processing time :" + str(round(processing_time, 3))+" Secs")
                        
                        # recoding time stamps in loginfor
                        loginfor[s_id] = static_infor
                        loginfor[s_id]["start_time"] = time.asctime(
                            time.localtime(start_time))
                        loginfor[s_id]["end_time"] = time.asctime(
                            time.localtime(end_time))
                        loginfor[s_id]["processing_time"] = str(int(processing_time))
                        loginfor[s_id]["trim_time"] = str(int(trim_time))
                        loginfor[s_id]["filter_time"] = str(int(filter_time))
                        loginfor[s_id]["blastn_time"] = str(int(blast_time))
            else:
                print("Skip indexing and alignment for FASTQ files")
            
            # Summarize all information from all fastq files
            static_dic = make_report.getTrfReportFile(
                proj_name, out_dir, tRNA_dic, sample_dic, out_dir)

            # Write down all login information
            self.writeLogFile(static_dic, loginfor, sample_dic)
            return 0
        except Exception as inst:
            print(type(inst))  # the exception instance
            print(inst)
            print('Some things wrong during pipeline running!')
            return -1
            # sys.exit(2)

    def load_tRNA_database(self,trna_anno_bed):
        # Load information about tRNAs if
        if not os.path.isfile(trna_anno_bed):
            print("The tRNA annotation bed file :" +
                    trna_anno_bed+" is not exist. Abort!")
            return -1
        tRNA_dic = {}
        with open(trna_anno_bed, 'r') as tRNAFile:
            for line in tRNAFile:
                if not line.startswith("#"):
                    t = tRNA()
                    t.LoadStr(line.strip())
                    tRNA_dic[t.name] = t
        return tRNA_dic

    def loadSampleInfor(self, sample_tsv, out_dir):
        # Load sample information
        sample_dic = {}
        if not os.path.isfile(sample_tsv):
            print("The sample file :"+sample_tsv+" is not exist. Abort!")
            return -1
        # Copy sample tsv for visualization
        if not os.path.isdir(out_dir):
            try:
                print(f"Output dir {out_dir} is not a directory. Try to create it ...append() ")
                os.mkdir(out_dir)
            except Exception as inst:
                print("Create output directory :"+out_dir+" fail. Abort!")
                print(inst)
                return -1
        cp_sample_tsv = out_dir+"/samples"
        shutil.copy(sample_tsv, cp_sample_tsv)
        with open(sample_tsv, 'r') as SAMPLES:
            for line in SAMPLES:
                if not line.startswith("#"):
                    contents = line.strip().split("\t")
                    if len(contents) > 1:
                        sample_dic[contents[0]] = {
                            "des": contents[1], "adapters": ""}
                    if len(contents) > 2:
                        sample_dic[contents[0]]['adapters'] = contents[2]
        return sample_dic

    def getFastqList(self, fastq_dir):
        ext = ".fastq"
        fastq_files = share.getExtFileList(fastq_dir, ".fastq")
        if len(fastq_files) == 0:
            fastq_files = share.getExtFileList(fastq_dir, ".fq")
            ext = ".fq"
        
        return [ext, fastq_files]

    def writeLogFile(self, static_dic, loginfor, sample_dic):
        cfg = self.config
        out_dir = cfg["out_dir"]
        # Write login information
        if cfg["no_indexing"] or cfg["no_alignment"]:
            print("No indexing and alignment so skip updating logfile.")
        else:
            # Write loginfor into a file
            for s_id in static_dic:
                if s_id in loginfor:
                    for key in static_dic[s_id]:
                        loginfor[s_id][key] = static_dic[s_id][key]
            logfile = out_dir + "/static.log"
            with open(logfile, 'w') as LOG:
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

    def fastq2Fasta(self, sample_dic, s_id, fastq_dir, trimmed_fastq):
        # Process adapter information
        cfg = self.config
        des = sample_dic[s_id]["des"]
        adapters = sample_dic[s_id]["adapters"].strip().split(',')
        f_adapter = str(adapters[0])
        r_adapter = ""
        if len(adapters) > 1:
            r_adapter = adapters[1]

        filter_start_time = time.time()
        # Removed redundant read Filter and get the read number file
        filtered_fasta = fastq_dir + "/" + s_id + "_filtered.fa"
        num_dic_txt = fastq_dir + "/" + s_id + "_num_dic.txt"
        static_infor = {}
        if not cfg["no_indexing"]:
            min_read_qscore = cfg["min_read_qscore"]
            min_reads_count = cfg["min_reads_count"]
            static_infor = share.filterFastQ2FastA(trimmed_fastq, filtered_fasta, num_dic_txt,
                                        f_patterm=f_adapter, r_patterm=r_adapter, qcutoff=min_read_qscore, num_cutoff=min_reads_count)
            print(static_infor)
        else:
            print("Skip filtering and indexing FASTQ flie : "+ trimmed_fastq)
            # delete trmmed fastq to save space
        if cfg["t_do"] != 0 and os.path.exists(trimmed_fastq):
            os.remove(trimmed_fastq)
        filter_end_time = time.time()
        filter_time = int(filter_end_time-filter_start_time)
        return [filter_time, static_infor, filtered_fasta, num_dic_txt]
    
    def doMapping(self, s_id, num_dic_txt, tRNA_dic, fasta):
        cfg = self.config
        blastn = cfg["blastn"]
        mkdb = cfg["mkdb"]
        trna_fa = cfg["trna_fa"]
        out_dir = cfg["out_dir"]
        blastn_e_cutoff = cfg["blastn_e_cutoff"]
        blastn_max_hit_num = cfg["blastn_max_hit_num"]
        proj_name = cfg["proj_name"]
        url_len = cfg["url_len"]
        blastn_max_mismatch = cfg["blastn_max_mismatch"]
        blastn_pident = cfg['blastn_pident']
        # Do BLASTN
        lastn_start_time = time.time()
        blast_out_file = blast_tools.RunBLASTN(
                blastn, mkdb, s_id, trna_fa, fasta, out_dir, eval=blastn_e_cutoff, hit_number=blastn_max_hit_num)
        #blast_out_file = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/output/SRR1836126_1_tRNA_blast_out.tab"
        # Analysis BLASTN result
        tRNA_reads_count_file = ""
        tRNA_reads_hit_file = ""
        if blast_out_file != "":
            print("BLASTN Successfully, processing BLAST out tab file: " + blast_out_file)
            tRNA_reads_count_file = out_dir + "/" + s_id + "_" + proj_name + "_count.tab"
            tRNA_reads_hit_file = out_dir + "/" + s_id + "_" + proj_name + "_hit.tab"
            blast_tools.AnalysisBlastOut2(blast_out_file, num_dic_txt, tRNA_dic,
                                tRNA_reads_count_file, tRNA_reads_hit_file, url_len, max_mismatch=blastn_max_mismatch, pident=blastn_pident)
            print("Hit file :"+tRNA_reads_hit_file)
        else:
            print("Something wrong while blastn " + s_id)
        blastn_end_time = time.time()
        blast_time = int(blastn_end_time-lastn_start_time)
        return [blast_time, tRNA_reads_count_file, tRNA_reads_hit_file]


    def runTrimmomatics(self, s_id, ext):
        cfg= self.config
        trim_time =0
        out_file = ""
        fastq_dir = cfg["fastq_dir"]
        raw_fastq = fastq_dir + "/" + s_id + ext
        if cfg["t_do"] != 0:
            trim_start_time = time.time()           
            t_adapter = cfg["t_adapter"]
            t_path = cfg["t_path"]
            t_phred = cfg["t_phred"]
            t_leading = cfg["t_leading"]
            t_tailing = cfg["t_tailing"]
            t_slidingwindow = cfg["t_slidingwindow"]
            t_minlen = cfg["t_minlen"]
            t_threads = cfg["t_threads"]
            # Trimed fastq
            raw_fastq = fastq_dir + "/" + s_id + ext
            cmd_bash = fastq_dir + "/" + s_id + "_blast.sh"
            trimmed_fastq = fastq_dir + "/" + s_id + "_trimmed"+ext
            with open(cmd_bash, "w") as CMD_FILE:
                adapter_fasta = t_adapter
                if not os.path.isfile(adapter_fasta):
                    adapter_fasta = ""
                    T = Trimmomatic(t_path)
                    trimmomatics_cmd = T.TrimSE(raw_fastq, trimmed_fastq, adapter_fa=adapter_fasta, phred=t_phred, LEADING=t_leading,
                                    TRAILING=t_tailing, SLIDINGWINDOW=t_slidingwindow, MINLEN=t_minlen, threads=t_threads)
                    CMD_FILE.write(trimmomatics_cmd + "\n")
                process = subprocess.Popen("bash " + cmd_bash, shell=True, stdout=subprocess.PIPE)
                process.wait()
            out_file = trimmed_fastq
            trim_end_time = time.time()
            trim_time = int(trim_end_time-trim_start_time)
        else:
            out_file = raw_fastq
        return [trim_time, out_file] 

def printHelpInfor():
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

# Main function
def main(argv):
    currentDirectory = os.getcwd()
    wdir = pathlib.Path(__file__).parent.absolute()
    c = Config(wdir)
    config_file = ""
    try:
        opts, args = getopt.getopt(sys.argv[1:], "c:n:f:a:s:i:o:h:v", [
                                   "config", 'no-indexing', 'no-alignment'])
    except getopt.GetoptError as err:
        print(err)
        printHelpInfor()
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
            printHelpInfor()
            sys.exit(0)
        elif opt == '-n':
            c.config["proj_name"] = arg
        elif opt == '-f':
            c.config["trna_fa"] = arg
        elif opt == '-a':
            c.config["trna_anno_file"] = arg
        elif opt == '-s':
            c.config["sample_tsv"] = arg
        elif opt == '-i':
            c.config["fastq_dir"] = arg
        elif opt == '-o':
            c.config["out_dir"] = arg
        elif opt == '-h':
            c.config["url_len"] = int(arg)
        elif opt == '--no-indexing':
            c.config["no_indexing"] = True
        elif opt == '--no-alignment':
            c.config["no_alignment"] = True
        elif opt in ("-c"):
            print("Load "+arg)
            config_file = arg
    print('##  ------------------ Initialization ....  -----------------')
    init_file = str(wdir)+"/init"
    if os.path.isfile(init_file):
        c.loadConfigFile(init_file)
        print("Loading init file successfully! ")
    else:
        print("No 'init' file was given, Using differnt!")
    print('##  ------------------ The loading config ....  -----------------')
    if config_file != "":
        if not os.path.isfile(config_file):
            config_file = currentDirectory+"/"+config_file
        if os.path.isfile(config_file):
            c.loadConfigFile(config_file)
            print("Loaded config file: " + config_file)
            c.printConfig()
        else:
            print("Can't find config file: " + config_file)
            print("Use default settings")
    else:
        print("No Config file was given, run for test data!")
    
    print('## Run tRNAExplorer pipeline ..... ')
    if rseqBlastnPipeline(c.config).run()==0:
        c.createVisualSettingFile()
        print('## The Pipeline processing done successfully ... ')
        exit(0)
    else:
        print('## Errors ocurr during processing  ... ')
        exit(-1)

if __name__ == "__main__":
    main(sys.argv[1:])

# Test
# proj_name= "test"
# config_file="/home/hqyone/mnt/sdc/trna/software/tRNAExplorer/config.txt"
#
# c = Config()
# c.loadConfig(config_file)
#c = Config()
#p = rseqBlastnPipeline(c.config)
#p.run()
# print("tRNA pipeline finish successfully!")
# sys.exit(0)
