/*
Copyright (C) 2014 Itoh Laboratory, Tokyo Institute of Technology

This file is part of Platanus.

Platanus is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Platanus is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with Platanus; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "consensus.h"
#include <sys/stat.h>
#include <sys/types.h>

using std::vector;
using std::string;


//////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////
Consensus::Consensus()
{
    optionSingleArgs["-k"] = "";
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-e"] = "";
    optionSingleArgs["-v"] = "32";
    optionSingleArgs["-l"] = "3";
    optionSingleArgs["-u"] = "0";
    optionSingleArgs["-L"] = "200000";
    optionSingleArgs["-t"] = "1";
    optionMultiArgs["-c"] = vector<string>();
    optionMultiArgs["-b"] = vector<string>();
    optionMultiArgs["-p"] = vector<string>();
    optionMultiArgs["-x"] = vector<string>();
    optionMultiArgs["-X"] = vector<string>();

    optionMultiArgs["-s"] = vector<string>(3);
    optionMultiArgs["-s"][0] = "32";
    optionMultiArgs["-s"][1] = "64";
    optionMultiArgs["-s"][2] = "96";

    optionMultiArgs["-S"] = vector<string>(1);
    optionMultiArgs["-S"][0] = "20";


    optionBool["-no_scaffold"] = false;
    optionBool["-unphase"] = false;

    pairedEndSingleFileType.push_back("-ip");
    pairedEndSingleFileType.push_back("-op");
    pairedEndPairFileType.push_back("-IP");
    pairedEndPairFileType.push_back("-OP");
    optionSingleArgs["-tmp"] = ".";

    optionSingleArgs["-mapper"] = "";
    optionBool["-minimap2"] = false;
    optionBool["-minialign"] = false;
    optionBool["-kmer_align"] = false;
}




//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void Consensus::usage(void) const
{

    std::cerr << "\nUsage: meta_platanus consensus [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : contig_file (fasta format)\n"
              << "    -k FILE                            : k-mer-occurrence file (binary format)\n"
//              << "    -b FILE1 [FILE2 ...]               : bubble_seq_file (fasta format)\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (interleaved file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (separate forward and reverse files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (interleaved, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (separate forward and reverse files, fasta or fastq)\n"
              << "    -p PAIR1 [PAIR2 ...]               : long-read file (PacBio, Nanopore) (fasta or fastq)\n"
              << "    -x PAIR1 [PAIR2 ...]               : linked-reads files (paired-ends, 10x Genomics) (interleaved, fasta or fastq)\n"
              << "    -X FWD1 REV1 [FWD2 REV2 ...]       : linked-reads files (paired-ends, 10x Genomics) (separate forward and reverse files, fasta or fastq)\n"
//              << "    -n{INT} INT                        : lib_id minimum_insert_size\n"
//              << "    -a{INT} INT                        : lib_id average_insert_size\n"
//              << "    -d{INT} INT                        : lib_id SD_insert_size\n"
              << "    -e FLOAT                           : coverage depth of homozygous region (default auto)\n"
              << "    -L INT                             : maximum fragment length of tag (10x Genomics) (default " << optionSingleArgs.at("-L") << ")\n"
              << "    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default " << optionMultiArgs.at("-s")[0] << " " << optionMultiArgs.at("-s")[1]  << ")\n"
              << "    -S INT1 [INT2 ...]                 : mapping seed length for long reads (default " << optionMultiArgs.at("-S")[0] << ", only effective with -kmer_align option)\n"
              << "    -l INT                             : minimum number of links to scaffold (default " << optionSingleArgs.at("-l") << ")\n"
              << "    -t INT                             : number of threads (<= " << optionSingleArgs.at("-t") << ", default 1)\n"
              << "    -mapper FILE                       : path of mapper executable file (default, minimap, only effective with -p option)\n"
              << "    -minimap2                          : use minimap2 insterd of minimap (default) for alignment of long reads\n"
//              << "    -minialign                         : use minialign insterd of minimap (default) for alignment of long reads\n"
//              << "    -kmer_align                        : use built-in alignment based on k-mer match insterd of minimap (default) for long reads\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n"

              << "\n\n"
              << "Outputs:\n"
              << "    PREFIX_consensusScaffold.fa\n"
              << "    PREFIX_consensusScaffoldComponent.bed\n"
              << "\n"
              << "Note, PREFIX is specified by -o\n"
              << std::endl;
}


void Consensus::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());

	std::cerr << "ececuting meta_platanus solve_DBG internally ..." << std::endl;
	this->execSolveDBG();
	std::cerr << "phase completed!" << std::endl;
}

void Consensus::execSolveDBG(void)
{
	vector<string> singleArgOptionList = {"-k", "-o", "-t", "-tmp", "-e", "-l", "-L", "-mapper"};
	vector<string> multiArgOptionList = {"-c", "-b", "-p", "-x", "-X", "-s", "-S"};

    std::ostringstream oss;
    oss << this->platanusBinary << " solve_DBG -unphase";

	for (auto optItr = singleArgOptionList.begin(); optItr != singleArgOptionList.end(); ++optItr) {
		if (!(optionSingleArgs[*optItr].empty())) {
			oss << " " << *optItr << " " << optionSingleArgs[*optItr];
		}
	}

	for (auto optItr = multiArgOptionList.begin(); optItr != multiArgOptionList.end(); ++optItr) {
		if (!(optionMultiArgs[*optItr].empty())) {
			oss << " " << *optItr;
			for (auto argItr = optionMultiArgs[*optItr].begin(); argItr != optionMultiArgs[*optItr].end(); ++argItr) {
				oss << " " << *argItr;
			}
		}
	}

    for (auto itr = optionPairFile.begin(); itr != optionPairFile.end(); ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (!(itr->fileSecond.empty())) {
            oss << " " << itr->fileSecond;
        }
    }

    if (system(oss.str().c_str()) != 0) {
        throw platanus::SolveDBGError();
    }
}
