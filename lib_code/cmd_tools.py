# coding=utf-8

#  tRNAExplorer v.1.0 - tRNA profiles and regulation networks
#  Copyright (C) 2020  Dr. Quanyuan He
#  School of Medicine, Hunan Normal University
#  Email: hqyone@hotmail.com
#  Freely distributed under the GNU General Public License (GPLv3)
#

import io, os
import subprocess

# An abstract class for commands tools
class CMDTool():
    def __init__(self, exe_path):
        self.exe_path = exe_path

    def getCMD(self, *argv, **kwargs):
        kwargs_ls= []
        for key, value  in kwargs.items():
            kwargs_ls.append("{}{}{}".format(key," ",value))
        return "{} {} {}".format(self.exe_path, " ".join(argv), " ".join(kwargs_ls))

    def executeCMD(self, cmd):
        try:
            print(cmd)
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            process.wait()
            if process.returncode==0:
                print('The command executes sucessfully!')
                return 0
            else:
                print('Command excutation failed with return {}'.format(process.returncode))
                return -1
        except:
            print("Errors happened druing execute the cmd:{}".format(cmd))
            return -1

# The details about output format can be found at
# https://sites.google.com/site/wiki4metagenomics/tools/blast/blastn-output-format-6
class BLASTN(CMDTool):
    def __init__(self, exe_path, output_format_str='-outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send qseq sseq qlen slen evalue"'):
        super().__init__(exe_path)
        self.output_format_str = output_format_str

    def alignment(self, db_fasta, query_fasta, eval, blast_out_file, hit_number=30, num_threads = 4, *argv, **kwrgs):
        if os.path.isfile(db_fasta) and os.path.isfile(query_fasta):
            cmd = self.getCMD('-task blastn -db {} -query {} -evalue {} -num_threads {} {} -num_alignments {} > {}'.format(db_fasta,query_fasta,eval,num_threads, self.output_format_str,hit_number,blast_out_file),  *argv, **kwrgs)
            #print(cmd)
            return self.executeCMD(cmd)
        else:
            print('BLASTN error: Some files can not be found.')
            return -1

class BLASTDBMaker(CMDTool):
    def __init__(self, exe_path="/Users/hqyone/Downloads/ncbi-blast-2.10.0+/bin/makeblastdb"):
        super().__init__(exe_path)
    
    def formatDB(self,fasta,out_name, dbtype='nucl', *argv, **kwrgs):
        if os.path.isfile(fasta):
            cmd = self.getCMD(" -in {} -dbtype {} -title {}".format(fasta,dbtype,out_name))
            return self.executeCMD(cmd)
        else:
            print("The referece fasta file ({}) is not exist!".format(fasta))
            return -1

class BWA(CMDTool):
    def __init__(self, exe_path):
        super().__init__(exe_path)
    
    def indexRefDB(self,genome_fa):
        if os.path.isfile(genome_fa):
            return self.executCMD(self.getCMD( index= genome_fa))
        else:
            print("The referece fasta file ({}) is not exist!".format(genome_fa))
            return -1
    
    def alignment(self, algorithm, fastq1, fastq2, genome_fa, out_sam,*argv, **kwrgs):
        cmd = self.getCMD(algorithm, fastq1, fastq2, genome_fa, *argv, **kwrgs)+" > "+out_sam
        self.executeCMD(cmd)

class SAMTOOLS(CMDTool):
    def __init__(self, exe_path):
        super().__init__(exe_path)

    def sort(self, bam, sorted_bam, process_num):
        if os.path.isfile(bam):
            return self.executeCMD(self.getCMD("sort",bam,"-@",process_num,'-o',sorted_bam))
        else:
            print("SAMTools(BAMSort Error): The bam file ({}) is not exist!".format(bam))
            return -1

    def indexBAM(self, sorted_bam):
        if os.path.isfile(sorted_bam):
            return self.executeCMD(self.getCMD("index",sorted_bam))
        else:
            print("SAMTools(Indexing BAM Error): The bam file ({}) is not exist!".format(sorted_bam))
            return -1

    def indexRefFASTA(self, ref_fasta):
        if os.path.isfile(ref_fasta):
            return self.executeCMD(self.getCMD("faidx",ref_fasta))
        else:
            print("SAMTools(Index Ref Error): The FASTA file ({}) is not exist!".format(ref_fasta))
            return -1

    def filterBAM(self,bam,bed,sam):
        if os.path.isfile(bam) and os.path.isfile(bed):
            return self.executeCMD(self.getCMD("view -L",bed,bam)+">"+sam)
        else:
            print("SAMTools(Filter Ref Error): BAM {} or BED {} file are missing!".format(bam, bed))
            return -1

class Trimmomatic(CMDTool):
    def __init__(self, exe_path):
        super().__init__(exe_path)


    def trimSE(self, input_fastq, out_fastq, adapter_fa="", phred="-phred33", LEADING=3,
               TRAILING=3, SLIDINGWINDOW="4:15", MINLEN=18, threads=2):
        ILLUMINACLIP = ""
        if os.path.isfile(adapter_fa):
            ILLUMINACLIP = "ILLUMINACLIP:{}:2:30:10".format(adapter_fa)
        else:
            ILLUMINACLIP = "ILLUMINACLIP:TruSeq3-SE:2:30:10"
            
        if os.path.isfile(input_fastq):
            return self.executeCMD("java -jar "+self.getCMD("SE",phred,'-threads',threads,input_fastq, out_fastq,\
                                                ILLUMINACLIP,\
                                                'LEADING:{}'.format(LEADING),\
                                                'TRAILING:{}'.format(TRAILING),\
                                                'SLIDINGWINDOW:{}'.format(SLIDINGWINDOW),\
                                                'MINLEN:{}'.format(MINLEN)
                                                ))
        else:
            print("Trimmomatic (trimSE Error): input_fastq {} does't exist!".format(input_fastq))
            return -1

class BWA_SAMTOOLS:
    def __init__(self, bwa_path="", samtools_path=""):
        self.bwa_path = "/Users/hqyone/Downloads/bwa-0.7.17/bwa"
        if bwa_path!="":
            self.bwa_path = bwa_path
        self.samtools_path = "/Users/hqyone/Downloads/samtools-1.10/samtools"
        if samtools_path!="":
            self.samtools_path = samtools_path
        self.p = 2

    def getAlignmentCMD(self, fastq1, fastq2, genome_fa, sorted_bam):
        cmd = ""
        if os.path.isfile(fastq1) and os.path.isfile(genome_fa):
            cmd = self.bwa_path + " mem -t " + str(self.p) + " " + genome_fa + " " + fastq1
            if os.path.isfile(fastq2):
                cmd += " " + fastq2
            cmd += " | "
            cmd += self.samtools_path+ " sort -@ "+str(self.p)+" -o " + sorted_bam+ " - "
        print(cmd)
        return cmd

class PICARD:
    def __init__(self,path=""):
        self.path = "/Users/hqyone/Downloads/picard.jar"
        if path!="":
            self.path = path
        self.ASSUME_SORTED = True
        self.REMOVE_DUPLICATE = True

    def getMarkDuplicatesCMD(self, bam, filterecd_bam, matrix):
        cmd = "java -jar "+self.path+" MarkDuplicates I="+bam+" O="+filterecd_bam + " M=" +matrix
        if self.ASSUME_SORTED:
            cmd += " ASSUME_SORTED=true "
        if self.REMOVE_DUPLICATE:
            cmd += " REMOVE_DUPLICATES=true "
        print(cmd)
        return cmd

class MACS2:
    def __init__(self, path=""):
        self.path = "MACS2"
        if path !="":
            self.path = path
        self.name = ""
        self.out_dir = ""
        self.q_value = 0
    def CallPeakBAMsCMD(self, sbam, cbam, out_dir, genome="hs"):
        cmd = self.path+" callpeak "
        cmd += " -t "+sbam
        cmd += " -c "+cbam
        if self.name != "":
            cmd += " -n " + self.name
        if self.q_value!=0:
            cmd+=" -q "+str(self.q_value)
        cmd += " -g "+genome
        cmd += " --outdir "+ out_dir
        print(cmd)
        return cmd


class FASTQ_DUMP:
    def __init__(self, path=""):
        self.path = "/Users/hqyone/OneDrive/Software_Tools/software/sratoolkit.2.9.2-mac64/bin/fastq-dump"
        if path!="":
            self.path = path
    def SRA2FASTQ(self, sra, out_dir="."):
        cmd = self.path + " --split-files "+sra + " "
        if out_dir!=".":
            cmd += " --outdir "+out_dir
        return cmd


"bamCoverage" #From deeptools
# using "sudo pip3 install deepTools==2.5.3"
# other version will result in "does not have BAM or CRAM format" error
class BAMCOVERAGE:
    def __init__(self, path=""):
        self.path = "/anaconda3/bin/bamCoverage"
        if path!="":
            self.path = path
        self.binsize = 30
        self.p = 1
    def getBAM2BiWigCMD(self, bam, bw):
        cmd = self.path+" -b "+bam+" -o "+bw+" -bs "+str(self.binsize) + " -p "+str(self.p)
        return cmd

#Deeptools
# https://deeptools.readthedocs.io/en/develop/content/tools/bigwigCompare.html
class BIGWIGCOMPARE:
    def __init__(self,path=""):
        self.path = "/anaconda3/bin/bigwigcompare"
        if path!="":
            self.path = path
        self.binSize = 50
        self.p = "max/2"
        self.scaleFactors = "1:1"
        self.region = ""  #CHR:START:END
        self.operation = "log2"  #log2,ratio,subtract,add,mean,reciprocal_ratio,first,second

    def getCommand(self, biwig1, biwig2, out_biwig):
        cmd =  self.path+" --binSize "+str(self.binSize)
        if self.region!="":
            cmd += " --region "+self.region
        cmd += " -o "+out_biwig
        cmd += " --bigwig1 "+ biwig1
        cmd += " --bigwig2 " + biwig2
        cmd += " --scaleFactor " + self.scaleFactors
        cmd += " --ratio " + self.operation
        return cmd

#computeMatrix reference-point -S <biwig file(s)> -R <bed file> -a 3000 -b 3000
# https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html
class COMPUTEMATRIX:
    def __init__(self, path=""):
        self.path= "/anaconda3/bin/computeMatrix"
        if path!="":
            self.path = path
        self.type="reference-point"
        self.binSize = 10
        self.sortRegions = "keep"

    def getCommand(self, biwig_ls, bed, out_gz, type="reference-point", a=1000, b=1000):
        cmd =  self.path+" "+type
        cmd += " -S " + " ".join(biwig_ls)
        cmd += " -R " + bed
        cmd += " -a " + str(a)
        cmd += " -b " + str(b)
        cmd += " -o " + out_gz
        return cmd

#https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html
class PLOTHEATMAP:
    def __init__(self, path=""):
        self.path = "plotheatmap"
        if path!="":
            self.path = path
        self.sortUsing = "mean"
        self.sortRegions = "descend"
        self.outFileSortedRegions=""
        self.kmeans = 3
        self.refPointLabel = "tRNA"
        # ‘Blues’, ‘BrBG’, ‘BuGn’, ‘BuPu’, ‘CMRmap’, ‘GnBu’, ‘Greens’, ‘Greys’, ‘OrRd’, ‘Oranges’, ‘PRGn’, ‘PiYG’, ‘PuBu’, ‘PuBuGn’, ‘PuOr’,
        self.colorMap = "CMRmap"
        self.alpha = 1
        self.colorList = ""  #black,yellow,blue
        self.zMax = -1
        self.heatmapWidth = 4
        self.plotTitle = ""
        self.legendLocation = "best" #best,upper-right,upper-left,upper-center,lower-left,lower-right
        self.ylabel = ""
        self.xlabel = ""
        self.dpi =  200

    def getComamnd(self, matrixFile, outFileName):
        cmd = self.path
        cmd += " -m " + matrixFile
        cmd += " -o " + outFileName
        cmd += " --sortUsing " + self.sortUsing
        if self.outFileSortedRegions!="" and os.path.isfile(self.outFileSortedRegions):
            cmd +=" --outFileSortedRegions "+self.outFileSortedRegions
        cmd += " --kmeans " + str(self.kmeans)
        cmd += " --sortRegions "+ self.sortRegions
        if self.refPointLabel!="":
            cmd += " --refPointLabel "+self.refPointLabel
        cmd += ' --alpha ' + str(self.alpha)
        cmd += ' --dpi '+str(self.dpi)
        cmd += " --colorMap " + self.colorMap
        if self.colorList!="":
            cmd += " --colorMap "+self.colorList
        if self.zMax >0:
            cmd += '--zMax '+str(self.zMax)
        cmd += " --heatmapWidth "+ str(self.heatmapWidth)
        if self.plotTitle!="":
            cmd += " --plotTitle " + self.plotTitle
        cmd += " --legendLocation " + self.legendLocation
        if self.xlabel!="":
            cmd += " --xAxisLabel "+self.xlabel
        if self.ylabel!="":
            cmd += " --yAxisLabel "+self.ylabel
        return cmd

#########################################################
# Files and Directories
#########################################################
def LoadConfig(config_file):
    config = {}
    if os.path.isfile(config_file):
        FILE = open(config_file,'r')
        for line in FILE:
            if line.startswith("#"):
                continue
            contents = line.split("=")
            if len(contents)==2:
                KEY = contents[0].strip().strip('"')
                VAL = contents[1].strip().strip('"')
                config[KEY] = VAL
    return config