# Platanus_B README.md


## Version
v1.2.3

## Web site
<http://platanus.bio.titech.ac.jp/>

## Author
Rei Kajitani at Tokyo Institute of Technology wrote key source codes.  
Address for this tool: <platanus@bio.titech.ac.jp>


## Requirements
* GCC 
    - <https://gcc.gnu.org/>
    - version >= 4.4, with OpenMP
    - To compile the source code.

* Minimap2
    - <https://github.com/lh3/minimap2>
    - Only required to use PacBio/Oxford-Nanopore long reads.
   
   
## Installation
```sh
make
cp meta_platanus <installation_path>
```


## Synopsis
### Inputs
* Illumina paired-end: PE_1.fq PE_2.fq (mandatory)

### Commands
```
platanus_b assemble -f PE_1.fq PE_2.fq 2>assemble.log

platanus_b iterate -c out_contig.fa out_junctionKmer.fa -IP1 PE_1.fq PE_2.fq 2>iterate.log
```

### Final output
    out_iterativeAssembly.fa
