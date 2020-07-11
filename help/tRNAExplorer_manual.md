# Manual for tRNAExplorer.py
## 1. Usage
Usage: python tRNAExplorer.py -c <configfile>
## 2. Settings

### 2.1 Command options
* User can use command options to launch pipeline in 

| Option  | Description  |
| :------------ |:--------------------------------| 
| -c     | The absolute path of config file | 
| -n     | The name of project  |  
| -f     | Absolute path of FASTA file for tRNAs which was created by [tRNA_db_maker.py](./tRNA_db_maker_manual.md) |  
| -a     | Absolute path of bed file for tRNA annotations which was created by  [tRNA_db_maker.py](./tRNA_db_maker_manual.md)  |  
| -s     | Absolute path of sample information  |  
| -i     | The directory storing fastq files (Input Directory)  |  
| -h     | Show the help information  |  
| -o     | The directory of output files (Out Directory)| 

### 2.2 Configure file [(see a sample)](config.txt)
* Configure file is a place to set all options available for the pipeline. 
* User can use [config.txt](config.txt) as a template to customize it.
* Removing/Adding the "#" at the start of line will enable/disable following option
* Notice: The option in configure file will overwrite the settings in command line

| Name  | Section  | Default  |Description
| :---------|:---------|:-------------|:--------------------------------| 
| proj_name  | General  | NA | Project name | 
| trna_fa    | General  | NA | Absolute path of fasta file for tRNAs which was created by tRNA_db_maker |   
| trna_anno_file    | General | NA  | Absolute path of bed file for tRNA annotations which was created by tRNA_db_maker |
| sample_tsv  | General | NA  | Absolute path of sample information | 
| fastq_dir    | General | NA | The directory storing FASTQ files (Input Directory) |   
| out_dir    | General | NA | The directory of output files (Out Directory) |
| url_len    | General | 60 | The length of UTRs, default = 60', it should match utr length in trna_anno bed |
| t_do  | Trimmomatic  | 1 | Whether trimming sequences 0 means "no" others means "yes"| 
| t_adapter    | Trimmomatic | Space | Absolute path for adapter FASTA file, When the file is not available using <trimmomatic_dir>+"/adapters/TruSeq3-SE.fa |   
| t_phred    | Trimmomatic  | 33 | phred cutoff to trim sequence |
| t_leading  | Trimmomatic  | 3 | Cut bases off the start of a read, if below a threshold quality | 
| t_tailing    | Trimmomatic  | 3 | Cut bases off the end of a read, if below a threshold quality |   
| t_slidingwindow   | Trimmomatic |"4:15" | Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold |
| t_minlen    | Trimmomatic  | 18 | Drop the read if it is below a specified length |
| t_threads   | Trimmomatic  | 2 | The number of thread to run trimmomatic |
| min_read_qscore  | Filtering  | 20 | Remove reads, if it contains any base below a threshold quality| 
| min_reads_count    | Filtering | 50 | Remove reads, if the read appears less than certain times in the FASTQ file |   
| blastn_e_cutoff    | BLASTN  | 0.001 | The evalue cutoff to report the match|
| blastn_max_hit_num  | BLASTN  | 40 | The maximum number of matched hits | 
| min_trf_len    | TRF  | 20 | The minimum of tRF length |   
| blastn_max_mismatch   | TRF | 2 | The max mismatch allow for valid matches |

## 3. Output
### 3.1 static.log [(see a sample)](../test/output/static.log)
The file includes reads numbers, processing time for each sample.

| Column  | Description  |
| :------------ |:--------------------------------| 
| #SampleID      | sample ID | 
| total_num     | Total number of reads in input FASTQ file  |  
| removed_num     | Total number of reads removed by trimming and filtering |  
| survived_num     | Total number of reads in the trimmed and filtered FASTQ file  |  
| non_redundant_num     | The number of non redundant reads in FASTA file  |  
| start_time     |  Time stamp to starting processing  |  
| end_time     | Time stamp to ending processing  |  
| processing_time     | Total time used for processing the sample  |  
| trim_time     | Time used by trimming (seconds)  |  
| filter_time     | Time used by filtering and removing redundant reads (seconds)  |  
| blastn_time     | Time used by mapping using BLASTN (seconds)  |  
| A - I    | The numbers of nine types of reads  |  
| intro_cl_ratio     | Intron cleavage ratio = F/(E+F)  |  
| u5_cl_ratio     | 5’ UTR cleavage ratio = (C)/(C+B)  | 
| u3_cl_ratio     | 3’ UTR cleavage ratio = (H+I)/(G+H+I)  | 
| cca_add_ratio     | CCA addition ratio = (I)/(I+H)| 


### 3.2 trf_sample_matrix.tsv 
* Read number matrix of tRFs across tRFs and samples
* See a [sample](../test/output/trf_sample_matrix.tsv)
| Column  | Description  |
| :------------ |:--------------------------------| 
| tRF_ID      | Sequence based ID of tRF | 
| tRNA_Families     | tRNA Family  |  
| tRNA_IDs     | ID list of tRNA which may generated the tRF |  
| seq     | Sequence of tRF  |  
| sample IDs ....     | List of sample IDs  |  

### 3.3 trna_sample_readcount_matrix.tsv 
* Read number matrix across tRNAs and samples
* See a [sample](../test/output/trna_sample_readcount_matrix.tsv)

| Column  | Description  |
| :------------ |:--------------------------------| 
| #RNA_family      | tRNA Family ID | 
| tRNA_ID     | ID of tRNA  |  
| sample IDs ....     | List of sample IDs  |    

### 3.4 trna_sample_pileup_matrix.tsv
* Pileup depth matrix across tRNAs and samples 
* See a [sample](../test/output/trna_sample_pileup_matrix.tsv)

| Column  | Description  |
| :------------ |:--------------------------------| 
| #RNA_family      | tRNA Family ID | 
| tRNA_ID     | ID of tRNA  |  
| sample IDs ....     | List of sample IDs  | 
  
### 3.5 trna_trftype_matrix.tsv
* The read number matrix across samples/tRNAs and tRF types
* See a [sample](../test/output/trna_trftype_matrix.tsv)

| Column  | Description  |
| :------------ |:--------------------------------| 
| #SampleID      | Sample ID | 
| tRNA_ID     | ID of tRNA  |  
| Sample IDs ....     | List of sample IDs  | 

### 3.6 sample_trftype_matrix.tsv
* The read number matrix across samples and tRF types
* See a [sample](../test/output/sample_trftype_matrix.tsv)

| Column  | Description  |
| :------------ |:--------------------------------| 
| #SampleID      | Sample ID | 
| tRF type ...     | The type of tRFs includes full_tRNA,   |  
| Sample IDs ....     | List of sample IDs  | 

### 3.7 cleavage_sites.tsv
* Cleavage sites information for tRNAs in different samples
* See a [sample](../test/output/cleavage_sites.tsv)

| Column  | Description  |
| :------------ |:--------------------------------| 
| #SampleID      | Sample ID | 
| tRF type ...     | The type of tRFs includes full_tRNA,   |  
| Sample IDs ....     | List of sample IDs  | 


### 3.8 profiles.tsv
* Pileup information for tRNAs in different samples
* See a [sample](../test/output/profiles.tsv)

| Column  | Description  |
| :------------ |:--------------------------------| 
| #SampleID      | Sample ID | 
| tRF type ...     | The type of tRFs includes full_tRNA,   |  
| Sample IDs ....     | List of sample IDs  | 
