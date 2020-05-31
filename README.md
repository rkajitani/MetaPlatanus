# MetaPlatanus README.md

## Description
MetaPlatanus is a de novo assembler for metagenome (microbiome). The features of this tool are as follows:
(1) It can utilize various types of long-range information such as mate-pairs (jumping libraries), Oxford-Nanopore/PacBio long reads and 10x linked reads (barcoded reads).
(2) Coverage depths, k-mer frequencies and results of the binning tool are also employed to extend sequences and correct mis-assemblies.
(3) Contig-assembly, scaffolding, gap-closing and binning are automatically executed at once.
(4) As the other mode, haplotypes (strain genomes) in a metagenome can be distinguished and constructed.
(5) MetaPlatanus requires at least one short-read paired-end library.


## Version
v1.2.2

## Web site
<http://platanus.bio.titech.ac.jp/>

## Author
Rei Kajitani at Tokyo Institute of Technology wrote key source codes.  
Address for this tool: <platanus@bio.titech.ac.jp>


## Dependency (the tools below are included in this package) 
* Minimap2
    - <https://github.com/lh3/minimap2>
    - Only required to use Oxford-Nanopore/PacBio long reads.
	- v2.17-r943

* BWA
    - <http://bio-bwa.sourceforge.net/>
    - v0.7.17-r1194

* MEGAHIT
    - <https://github.com/voutcn/megahit>
	- v1.1.3

* MetaBAT2
    - <https://bitbucket.org/berkeleylab/metabat>
	- v2.12.1

* SAMtools
    - <http://www.htslib.org/>
	- v1.3.1


## Installation
Currently MetaPlatanus can be executed in Linux.
To compile (build), just type "make".
```
make
```

The main program is "meta_platanus.pl".
Note that the directory "sub_bin", which consists of Perl-scripts and other tools, should be specified (-sub_bin option) or put in the same directory of meta_platanus.pl.
There are four ways to execute meta_platanus.pl.
Below, "build-directory" means the directory where you typed "make" and it consists of meta_platanus.pl and sub_bin.

(1) Directly type the path to build meta_platanus.pl.
```
/build-directory/mata_platanus.pl ...
```

(2) Add the build-directory to PATH.
```
export PATH=/build-directory/:$PATH
meta_platanus.pl ...
```

(3) Copy meta_platanus.pl and sub_bin to a directory listed in PATH (e.g. $HOME/bin).
```
cp -r sub_bin meta_platanus.pl $HOME/bin
meta_platanus.pl ...
```

(3) Copy meta_platanus.pl to a directory listed in PATH, and use the "-sub_bin" option.
```
cp meta_platanus.pl $HOME/bin
meta_platanus.pl [command] -sub_bin /build-directory/sub_bin ...
```


## Synopsis of consensus assembly
### Inputs
* Illumina paired-end: PE_1.fq PE_2.fq (mandatory)
* Illumina mate-pair: MP_1.fq MP_2.fq (optional)
* Oxford Nanopor long-reads: ONT.fq (optional)

### Commands
```
meta_platanus.pl cons_asm -IP1 PE_1.fq PE_2.fq -OP2 MP_1.fq MP_2.fq -ont ONT.fq >log.txt 2>&1
```

### Final output
    out_finalClusters_all.fa
    out_finalClusters (directory)


## Synopsis of haplotype-phasing assembly
### Inputs
* Illumina paired-end: PE_1.fq PE_2.fq (mandatory)
* Illumina mate-pair: MP_1.fq MP_2.fq (optional)
* Oxford Nanopor long-reads: ONT.fq (optional)

### Commands
```
meta_platanus.pl phase_asm -IP1 PE_1.fq PE_2.fq -OP2 MP_1.fq MP_2.fq -ont ONT.fq >log.txt 2>&1
```

### Final output
    out_allPhasedBlock.fa (including sequences below)
    out_primaryBubble.fa
    out_secondaryBubble.fa
    out_nonBubbleOther.fa
	...



---
## Consensus assembly usage
### Command
```
meta_platanus.pl cons_asm -IP1 short_R1.fastq(a) short_R2.fastq(a) [Options] ...
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
    -megahit_min_len                   : minimum length of contigs of MEGAHIT (default, 500)
    -overwrite                         : overwrite the previous results, not re-start (default, off)
    -h, -help                          : display usage
```


### Outputs:
    PREFIX_finalClusters_all.fa
    PREFIX_finalClusters (directory

PREFIX is specified by -o
  
  
## Assemble each haplotype in a metanogme (phasing)
### Command
```sh
meta_platanus.pl phase_asm -IP1 short_R1.fastq(a) short_R2.fastq(a) [Options] ...
```

### Options
```
    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq; at least one library required)
    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq; aka mate-pairs or jumping-library)
    -p FILE1 [FILE2 ...]               : PacBio long-read file (fasta or fastq)
    -ont FILE1 [FILE2 ...]             : Oxford Nanopore long-read file (fasta or fastq)
    -x PAIR1 [PAIR2 ...]               : barcoded_pair_files (10x Genomics) (reads in 1 file, interleaved, fasta or fastq)
    -X FWD1 REV1 [FWD2 REV2 ...]       : barcoded_pair_files (10x Genomics) (reads in 2 files, fasta or fastq)
    -t INT                             : number of threads (<= 1; default, 1)
    -m INT                             : memory limit for making kmer distribution (unit, GB; default, 0.75 * available_memory))
    -o STR                             : prefix of output files (default "out")
    -tmp DIR                           : directory for temporary files (default, ".")
    -sub_bin DIR                       : directory for sub-executables, such as mata_plantaus and minimap2 (default, directory-of-this-script/sub_bin)
    -overwrite                         : overwrite the previous results, not re-start (default, off)
    -h, -help                          : display usage
```


### Outputs:
```
    PREFIX_allPhasedBlock.fa (including sequences below)
    PREFIX_primaryBubble.fa
    PREFIX_secondaryBubble.fa
    PREFIX_nonBubbleOther.fa
	...
```

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
meta_platanus.pl cons_asm \
	 -IP1 sample1_R1.fq  sample1_R2.fq \
	-binning_IP1 sample1_R1.fq  sample1_R2.fq \
	-binning_IP2 sample2_R1.fq  sample2_R2.fq \
	-binning_IP3 sample3_R1.fq  sample3_R2.fq \
...
```