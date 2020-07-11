# Manual for tRNAExplorer.py
## 1. Usage
Usage: python tRNAExplorer.py -c <configfile>
## 2. Settings
### 2.1 -c : path to config file
Config file include all options available to 

## 3. Output
### 3.1 static.log
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
Read number matrix of tRFs across tRFs and samples

| Column  | Description  |
| :------------ |:--------------------------------| 
| tRF_ID      | Sequence based ID of tRF | 
| tRNA_Families     | tRNA Family  |  
| tRNA_IDs     | ID list of tRNA which may generated the tRF |  
| seq     | Sequence of tRF  |  
| sample IDs ....     | List of sample IDs  |  

### 3.3 trna_sample_readcount_matrix.tsv
Read number matrix across tRNAs and samples

| Column  | Description  |
| :------------ |:--------------------------------| 
| #RNA_family      | tRNA Family ID | 
| tRNA_ID     | ID of tRNA  |  
| sample IDs ....     | List of sample IDs  |    

### 3.4 trna_sample_pileup_matrix.tsv
Pileup depth matrix across tRNAs and samples 

| Column  | Description  |
| :------------ |:--------------------------------| 
| #RNA_family      | tRNA Family ID | 
| tRNA_ID     | ID of tRNA  |  
| sample IDs ....     | List of sample IDs  | 
  
### 3.5 trna_trftype_matrix.tsv
The read number matrix across samples/tRNAs and tRF types

| Column  | Description  |
| :------------ |:--------------------------------| 
| #SampleID      | Sample ID | 
| tRNA_ID     | ID of tRNA  |  
| Sample IDs ....     | List of sample IDs  | 

### 3.6 sample_trftype_matrix.tsv
The read number matrix across samples and tRF types

| Column  | Description  |
| :------------ |:--------------------------------| 
| #SampleID      | Sample ID | 
| tRF type ...     | The type of tRFs includes full_tRNA,   |  
| Sample IDs ....     | List of sample IDs  | 

### 3.7 cleavage_sites.tsv
Cleavage sites information for tRNAs in different samples
### 3.8 profiles.tsv
Pileup information for tRNAs in different samples
