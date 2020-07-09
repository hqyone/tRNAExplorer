# coding=utf-8
import io, os
#########################################################
# Programs
#########################################################
class BWA:
    def __init__(self, path=""):
        self.path = "/Users/hqyone/Downloads/bwa-0.7.17/bwa"
        if path != "":
            self.path = path
        self.p = 2
    def getIndexingCMD(self, ref_fasta):
        return self.path+ ' index ' +ref_fasta

    def getAlignmentCMD(self, fastq1, fastq2, genome_fa, sam):
        cmd = ""
        if os.path.isfile(fastq1) and os.path.isfile(genome_fa):
            cmd = self.path + " mem -t " + str(self.p)+" "+genome_fa+" "+fastq1
            if os.path.isfile(fastq2):
                cmd += " "+fastq2
            cmd += ">"+sam
        print(cmd)
        return cmd

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

# The details about output format can be found at
# https://sites.google.com/site/wiki4metagenomics/tools/blast/blastn-output-format-6
class BLASTN:
    def __init__(self, cmdpath = "", mkdbpath=""):
        self.cmdpath = "/Users/hqyone/Downloads/ncbi-blast-2.10.0+/bin/blastn"
        if cmdpath!="":
            self.cmdpath = cmdpath
        self.mkdbpath = "/Users/hqyone/Downloads/ncbi-blast-2.10.0+/bin/makeblastdb"
        if mkdbpath!="":
            self.mkdbpath = mkdbpath
    def getCreateBLASTdbCMD(self,fasta,out_name):
        cmd = self.mkdbpath+" -in " + fasta + " -dbtype nucl -title "+out_name+"\n"
        return cmd

    def getAlignmentCMD(self, db_fasta, query_fasta, eval, blast_out_file, hit_number=30):
        cmd = ""
        output_format_str = '-outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send qseq sseq qlen slen evalue"'
        if os.path.isfile(db_fasta) and os.path.isfile(query_fasta):
            cmd = self.cmdpath + " -task blastn -db " + db_fasta + " -query " + \
                query_fasta + " -evalue " + str(eval) + " -num_threads 8 " + output_format_str + \
                "  -num_alignments " + \
                str(hit_number) + " > " + blast_out_file
        return cmd

class SAMTOOLS:
    def __init__(self, path=""):
        self.path = "/Users/hqyone/Downloads/samtools-1.10/samtools"
        if path!="":
            self.path = path
        self.p = 1

    def getSortCMD(self, bam, sorted_bam):
        cmd =  self.path+ ' sort '+bam+' -@ '+str(self.p)+' -o '+ sorted_bam
        print(cmd)
        return cmd

    def getIndexingCMD(self, sorted_bam):
        cmd = self.path+ ' index ' + sorted_bam
        print(cmd)
        return cmd

    def getIndexingGenomeCMD(self, ref_fasta):
        cmd = self.path+ ' faidx ' + ref_fasta
        print(cmd)
        return cmd

    def FilterBAMCMD(self,bam,bed,sam):
        cmd = self.path+" view -L " + bed + " " + bam + " > " + sam + "\n"
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

class Trimmomatic:
    def __init__(self, path=""):
        self.path = "/Users/hqyone/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar"
        if path !="":
            self.path = path

    def TrimSE(self, input_fastq, out_fastq, adapter_fa="", phred=33, LEADING=3,
               TRAILING=3, SLIDINGWINDOW="4:15", MINLEN=18, threads=2):
        cmd = "java -jar "+ self.path+" SE "
        if phred==33:
            cmd += " -phred33 "
        else:
            cmd += " -phred64 "
        cmd += " -threads " + str(threads)
        cmd += " " +input_fastq+" " +out_fastq
        if adapter_fa=="" or os.path.isfile(adapter_fa):
            trimmomatic_dir = os.path.dirname(os.path.abspath(self.path))
            default_adapter_fa = trimmomatic_dir+"/adapters/TruSeq3-SE.fa"
            if os.path.isfile(default_adapter_fa):
                cmd += " ILLUMINACLIP:"+default_adapter_fa+":2:30:10 "  # Adapter for illuminate
            else:
                cmd += " ILLUMINACLIP:TruSeq3-SE:2:30:10 "  # Adapter for illuminate
        cmd += " LEADING:"+str(LEADING)
        cmd += " TRAILING:" + str(TRAILING)
        cmd += " SLIDINGWINDOW:" + str(SLIDINGWINDOW)
        cmd += " MINLEN:" + str(MINLEN)
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