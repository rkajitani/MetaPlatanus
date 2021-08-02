# MetaPlatanus README.md

## Description
MetaPlatanus is a de novo assembler for metagenome (microbiome). The features of this tool are as follows:
(1) It can utilize various types of long-range information such as Oxford-Nanopore/PacBio long reads, mate-pairs (jumping libraries), and 10x linked reads (barcoded reads).
(2) Coverage depths, k-mer frequencies and results of the binning tool are also employed to extend sequences and correct mis-assemblies, reducing inter-species misassemblies.
(3) Contig-assembly, scaffolding, gap-closing and binning are automatically executed at once.
(4) MetaPlatanus requires at least one short-read paired-end library.

## Version
v1.3.0

## Web site
<https://github.com/rkajitani/MetaPlatanus>
<http://platanus.bio.titech.ac.jp>  

## Author
Rei Kajitani at Tokyo Institute of Technology wrote key source codes.  
Address for this tool: <platanus@bio.titech.ac.jp>


## Installation
Currently MetaPlatanus can be executed in Linux. There are two ways.

1. Using Miniconda or Anaconda  
MetaPlatanus is registered in Bioconda channel; Linux only.
```sh
conda install -c conda-forge -c bioconda mataplatanus
```

2. Using Docker  
For any OS if Docker is available.
```sh
docker pull rkajitani/metaplatanus
# run metaplatanus
docker run rkajitani/metaplatanus metaplatanus 
```

3. Build from source  
Install the dependencies above.
If the following commands are available, you will be able to run metaplatanus.
  - minimap2
  - samtools
  - seqkit
  - metabat2
  - jgi_summarize_bam_contig_depths
  - bwa
  - tgsgapcloser
  - racon
  - nextPolish

Compile (make), and copy metaplatanus and sub_bin to a directory listed in PATH (e.g. $HOME/bin).
```sh
make
cp -r sub_bin metaplatanus $HOME/bin
```
The main program is "metaplatanus".
Note that the directory "sub_bin", which consists of Perl-scripts and other tools, should be specified (-sub_bin option) or put in the same directory of metaplatanus.
There are two ways to install metaplatanus.


## Synopsis
### Inputs
* Illumina paired-end: PE_1.fq PE_2.fq (mandatory)
* Oxford Nanopor long-reads: ONT.fq (optional)

### Commands
```sh
metaplatanus -IP1 PE_1.fq PE_2.fq -ont ONT.fq >log.txt 2>&1
```

### Output
The files below in out_result (directory). The prefix "out" can be changed using "-o" option.
* out_final.fa  
  Final scaffolds as one FASTA file (gap-closed and polished).
* out_finalClusters  
  Final scaffolds as separated files of bins.
* out_finalClusters.tsv  
  TSV file of scaffold-bin relations.
* out_preClose.fa  
  A scaffold FASTA file before processes of TGS-GapCloser and NextPolish.
* out_preCloseClusters  
  Separated FASTA files before processes of TGS-GapCloser and NextPolish.
* out_preCloseClusters.tsv  
  TSV file of scaffold-bin relations before processes of TGS-GapCloser and NextPolish.


## Dependency (the tools below are included in this package) 
* Minimap2
    - <https://github.com/lh3/minimap2>
    - Only required to use Oxford-Nanopore/PacBio long reads.
	- v2.17-r943 or newer

* BWA
    - <http://bio-bwa.sourceforge.net/>
    - v0.7.17-r1194 or newer

* MetaBAT2
    - <https://bitbucket.org/berkeleylab/metabat>
	- v2.12.1 or newer

* SAMtools
    - <http://www.htslib.org/>
	- v1.3.1 or newer

* SeqKit
    - <https://bioinf.shenwei.me/seqkit/>
	- v1.16.1 or newer

* TGS-GapCloser (optional)
    - <https://github.com/BGI-Qingdao/TGS-GapCloser>
	- v1.0.1 or newer

* Racon (optional)
    - <https://github.com/isovic/racon>
	- v1.4.20 or newer

* NextPolish (optional)
    - <https://github.com/Nextomics/NextPolish>
	- v1.3.1 or newer


---
## Usage
### Command
```sh
metaplatanus -IP1 short_R1.fastq(a) short_R2.fastq(a) [Options] ...
```

### Options
```
-IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq; at least one library required)
-OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq; aka mate-pairs or jumping-library)
-binning_IP{INT} FWD1 REV1 ...     : lib_id inward_pair_files for binning process. (reads in 2 files, fasta or fastq; the data are usually from another sample)
-p FILE1 [FILE2 ...]               : PacBio long-read file (fasta or fastq)
-ont FILE1 [FILE2 ...]             : Oxford Nanopore long-read file (fasta or fastq)
-x PAIR1 [PAIR2 ...]               : barcoded_pair_files (10x Genomics) (reads in 1 file, interleaved, fasta or fastq)
-X FWD1 REV1 [FWD2 REV2 ...]       : barcoded_pair_files (10x Genomics) (reads in 2 files, fasta or fastq)
-t INT                             : number of threads (<= 1; default, 1)
-m INT                             : memory limit for making kmer distribution (unit, GB; default, 0.75 * available_memory))
-o STR                             : prefix of output files (default "out")
-tmp DIR                           : directory for temporary files (default, ".")
-sub_bin DIR                       : directory for sub-executables, such as mata_plantaus and minimap2 (default, directory-of-this-script/sub_bin)
-min_cov_contig                    : k-mer coverage cutoff for contig-assembly of MetaPlatanus (default, 2)
-overwrite                         : overwrite the previous results, not re-start (default, off)
-no_tgsgapcloser                   : do not use TGS-GapCloser and NextPolish (default, off)
-no_nextpolish                     : do not use NextPolish (default, off)
-h, -help                          : display usage
-v, -version                       : display version 
```


### Outputs:
    PREFIX_result (directory)

PREFIX is specified by -o
  
  
---
## Notes
* Options related to run time
Although -t (number of threads) of all commands and -m (memory amount) of the "assemble" command are 
not mandatory to run, it is recommended to set the values adjusting your machine-environment.
These options may severely effect the run time.  
e.g.,  
Available number of threads and memory amount are 4 and 16GB, respectively.  
->  -t 4 -m 16

* Paired-end (mate-pair) input  
The "phase" and "consensus" accept paired-end and/or mate-pair libraries. Paired libraries are 
classified into "inward-pair" and "outward-pair" according to the sequence direction. 
For file formats, separate and interleaved files can be input through -IP (-OP) and -ip (-op) 
options, respectively.

Inward-pair (usually called "paired-end", accepted in options "-IP" or "-ip"):

    FWD --->
        5' -------------------- 3'
        3' -------------------- 5'
                        <--- REV 

Outward-pair (usually called "mate-pair", accepted in options "-OP" or "-op"):

                        ---> REV 
        5' -------------------- 3'
        3' -------------------- 5'
    FWD <---

Example inputs:

    Inward-pair (separate, insert=300)   : PE300_1.fq PE300_2.fq
    Outward-pair (separate, insert=2k)   : MP2k_1.fa MP2k_2.fq

Corresponding options:

    -IP1 PE300_1_pair.fq PE300_2.fq \
    -OP2 MP2k_1.fq MP2k_2.fq


To utilize multiple-samples data, MetaPlatanus can accept the short-reads
that exclusivelly used for the binning process through -binning-IP# options.
e.g.,
```
metaplatanus \
	-IP1 sample1_R1.fq sample1_R2.fq \
	-binning_IP1 sample1_R1.fq sample1_R2.fq \
	-binning_IP2 sample2_R1.fq sample2_R2.fq \
	-binning_IP3 sample3_R1.fq sample3_R2.fq \
...
```
