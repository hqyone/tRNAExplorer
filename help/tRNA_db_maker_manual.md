# Manual for tRNA_db_maker.py
The program can generate two file as a tRNA database which will be used by downstream [tRNAExplorer.py](../tRNAExplorer.py) processing
## 1. Requirements
1. tRNA bed file which can be download from UCSC table browser ([Sample](../test/trna_db/hg38_tRNA.bed)) (https://genome.ucsc.edu/cgi-bin/hgTables)
2. FASTA file of the genome with same version using by tRNA bed file.
3. tRNAScan-SE should be download and installed 

## 1. Usage
Usage: `python tRNA_db_maker.py -n <name> -b <bed> -r <ref> -s <tRNAScanSE> -o <offset> --no_mit <1/else> --no_pseu <1/else> --minq <number>`
## 2. Settings

| Option  | Default | Description  |
| :------------| :------ |:--------------------------------| 
| -n    |NA| The name of tRNA database  |  
| -b     |NA | BED file of tRNA genes which can be downloaded from UCSC table browser (https://genome.ucsc.edu/cgi-bin/hgTables) |  
| -r      |NA| Absolute path of the FASTA file of genome  |  
| -s      |NA| Absolute path of tRNAScan-SE  |  
| -o     |60 | The length of UTR |  
| --no_mit    |1| 1 means removing all nucleic mitochondrial tRNA genes |  
| --no_pseu   |1 | 1 means removing all pseudogene tRNA genes  | 
| --minq | 30| The minimum of quality scores of tRNAs  | 


## 3. Output

tRNA_db_maker will create two files:

1. FASTA file of hybrid tRNA library ([Example](../test/trna_db/hg38_tRNA_60.fa))
   *   pre-tRNA with intron(s) (I)
   *   pre-tRNA without intron(s) (P)
   *   mature tRNA (M)
   *   mature tRNA with CCA (C)
   
2. A extended bed file with tRNAs sequence and structure annotations ([Example](../test/trna_db/hg38_tRNA_60.bed))
* All position are 1 based

| Column  | Description  |
| :------------ |:--------------------------------| 
| #chrom    | Chromosome name  |  
| start     | Start position in chromosome | 
| end     | End position in chromosome |   
| name     | Name of tRNA  |  
| score     | No used , always 0  |  
| strand     | The strand of tRNA gene (-/+) |  
| full_seq     | tRNA sequence with intron |  
| seq     | Mature tRNA sequence without intron  | 
| intron_infor | A string describe intron : intron: 38-58 (98-118) | 
| anticodon     | Anticodon of the tRNA |  
| anticodon_start  | start position at UTR+tRNA with intron | 
| anticodon_end | end position at UTR+tRNA with intron |
| acceptor | Acceptor (Amino acid carried by the tRNA)|
| family | Family ID|
| seq_utr | tRNA seq |
| utr_len | The length of URT |
| map_start | The start position of region predicted by tRNAScan-SE in tRNA|
| map_end | The end position of region predicted by tRNAScan-SE in tRNA|
| map_len | The length of region predicted by tRNAScan-SE in tRNA|
| map_structure_str| String for structure prediction predicted by tRNAScan-SE|
| map_seq|The sequenced predicted by tRNAScan-SE in tRNA. It should be the same of full_seq |
| anticodon_start_in_map|The start position of anticodon in map_seq |
| anticodon_end_in_map|The end position of anticodon in map_seq|
| map_scan_score|tRNAScan-SE score indicating possibility to be tRNA |
| possible_type|pseudogene or not|
| d_loop_start| The start position of D loop stem |
| d_loop_lstart|The start position of D loop (not including the stem of D loop)|
| d_loop_lend|The end position of D loop (not including the stem of D loop)|
| d_loop_end|The end position of A loop stem|
| a_loop_start| The start position of A loop stem |
| a_loop_lstart|The start position of A loop (not including the stem of A loop)|
| a_loop_lend|The end position of A loop (not including the stem of A loop)|
| a_loop_end|The end position of A loop stem|
| v_loop_start| The start position of Variant loop stem |
| v_loop_lstart|The start position of Variant loop (not including the stem of Variant loop)|
| v_loop_lend|The end position of Variant loop (not including the stem of Variant loop)|
| v_loop_end|The end position of Variant loop stem|
| t_loop_start| The start position of T loop stem |
| t_loop_lstart|The start position of T loop (not including the stem of T loop)|
| t_loop_lend|The end position of T loop (not including the stem of T loop)|
| t_loop_end|The end position of T loop stem|
| stem_for_start|The start position of 5' stem of tRNA |
| stem_for_end|The end position of 5' stem of tRNA |
| stem_rev_start|The start position of 3' stem of tRNA |
| stem_rev_end|The end position of 3' stem of tRNA |
   

