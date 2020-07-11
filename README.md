# tRNAExplorer (v1.0)
## 1. Introduction
tRNAExplorer is a Python pipeline optimized for analysing tRF (tRNA-derived fragments) profiles of multiple samples using small-RNA-seq data.
* Major functions:
    * Categorize and quantify tRNA/tRFs.
    * Discover novel tRNA cleavage sites.
    * Discover new base additions on tRFs 
    * Comparison and visualization of tRF profile across multiple samples

## 2. Requirements
1. Trimmomatics 0.39 http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
2. BLASTN (ncbi-blast-2.10.0+) https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
3. tRNAScan-SE  http://lowelab.ucsc.edu/tRNAscan-SE/
4. Python (version 2.7 or 3.x) and packages
   * jupyter>=1.0.0
   * scipy==1.1.0
   * numpy>=1.16.1
   * pandas==0.23.0
   * matplotlib==2.2.5
   * seaborn>=0.9.0
   * pysam >= 0.8
   * HTSeq>=0.11.0
   * pathlib

## 3. Installation
* Step 1: Download and unzip the pipeline  https://github.com/hqyone/tRNAExplorer/archive/master.zip
* Step 2: Install requirements ...
* Step 3: Initialization settings:
    * Find and modify the file named [init](init) at the root directory of tRNAExplorer 
    * Change the absolute path for Trimmomatics, BLASTN and makeblastdb in the file.

## 4. Running
* Run test data by      `python <root path>/tRNAExploer.py`
* Get help information  `python <root path>/tRNAExploer.py -h`
* Create a sample tsv file with two columns ([example](./test/samples))

    | Column  | Description  |
    | :------------ |:--------------------------------| 
    | ID    | ID of sample, match the name of FASTQ file in fastq_dir. For example: ID.fastq or ID.fq | 
    | Description     | The short string descript the sample |  

* Run with customized data: `python <root path>/tRNAExploer.py -n <proj_name> -f <trna_fa> -a <trna_anno_file> -s <sample tsv> -i <fastq_dir> -o <out_dir>`
* Run with advanced settings using config file: 
    *   Modify config.txt  (taken [config.txt](config.txt) in root directory as template)
    *   Run `python <root path>/tRNAExplorer.py -c config.txt`
    *   Details about tRNAExploer.py can be found [here](./help/tRNAExplorer_manual.md)
* Make Customized tRNA Databases
    *   Run `python tRNA_db_maker.py -n <name> -b <bed> -r <ref> -s <tRNAScanSE> -o <offset> --no_mit <1/else> --no_pseu <1/else> --minq <number>`
    *   Details about tRNA_db_maker.py can be found [here](./help/tRNA_db_maker_manual.md)

## 5. Mapping Strategy
tRNAExplorer uses BLASTN as the engine to map small-RNA-Seq reads to a hybrid tRNA sequences database containing four major style of tRNA gene transcripts:
   *   pre-tRNA with intron(s) (I)
   *   pre-tRNA without intron(s) (P)
   *   mature tRNA (M)
   *   mature tRNA with CCA (C)
 
 Only best alignment hits will be keep for further analysis. 
 Based on the mapped transcript types and mapping locations, all mapped reads can be categorized into 9 types (A-I). 
 Using this information we can elucidate the relative abundance of four transcript types for each tRNA.
 We also can monitor the efficiency of intron, 5'-term & 3'-term cleavage and CCA addition. (Figure 1)
 
![alt text](./images/read_classfication.png)

Figure 1: The definitions of four tRNA transcript types and nine read types. 

## 6. Architecture
The tRNAExplorer contains three modules:
*   tRNAExplorer.py : main program to mapping and quantify tRNA/tRF
*   tRNA_db_maker.py : program to format the tRNA database which will be used by tRNAExplorer.py
*   Report_create.ipynb : program to analysis data and draw chart to visualization data.

Processing steps are summarized in Figure 2 

![alt text](./images/architecture.png)

Figure 2. The architecture of tRNAExplorer. 

## 7. tRF types
Actually, tRFs are any RNA fragments or transcripts that divided from tRNA genes.
tRNAExplorer defines 12 types of tRFs based on their mapping place on tRNA genes.

![alt text](./images/trf_types.png)

Figure 3. tRF types in tRNAExplorer 

## 8. Visualization
tRNAExplorer will generate several tsv files. 
A jupyter notebook [Report_creator.ipynb](./visualization/Report_creater.ipynb) to analysis and visualization them.
This notebook is self-explained.
The user can also do their own analysis based on these tsv files.

## 9. License
Copyright (c) 2020 Quanyuan He Ph.D, 
School of Medicine, Hunan Normal University.
Released under GPLv3. See
[license](LICENSE.txt) for details.