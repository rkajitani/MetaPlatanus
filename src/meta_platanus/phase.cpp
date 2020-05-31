/*
Copyright (C) 2017 Itoh Laboratory, Tokyo Institute of Technology

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

#include "phase.h"
#include <sys/stat.h>
#include <sys/types.h>

using std::vector;
using std::string;


//////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////
Phase::Phase()
{
    optionSingleArgs["-k"] = "";
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-e"] = "";
    optionSingleArgs["-v"] = "32";
    optionSingleArgs["-j"] = "1";
    optionSingleArgs["-l"] = "3";
    optionSingleArgs["-L"] = "200000";
    optionSingleArgs["-t"] = "1";
    optionMultiArgs["-c"] = vector<string>();
    optionSingleArgs["-i"] = "2";
    optionSingleArgs["-tmp"] = ".";

    pairedEndSingleFileType.push_back("-ip");
    pairedEndSingleFileType.push_back("-op");
    pairedEndPairFileType.push_back("-IP");
    pairedEndPairFileType.push_back("-OP");
    optionMultiArgs["-p"] = vector<string>();
    optionMultiArgs["-x"] = vector<string>();
    optionMultiArgs["-X"] = vector<string>();

    optionSingleArgs["-mapper"] = "";
    optionBool["-minimap2"] = false;
    optionBool["-minialign"] = false;
    optionBool["-kmer_align"] = false;
}




//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void Phase::usage(void) const
{

    std::cerr << "\nUsage: platanus2 phase [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file and directory (do not use \"/\", default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : contig (or scaffold) file (fasta format)\n"
              << "    -k FILE                            : k-mer-occurrence file (binary format)\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (interleaved file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (separate forward and reverse files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (interleaved, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (separate forward and reverse files, fasta or fastq)\n"
              << "    -p PAIR1 [PAIR2 ...]               : long-read file (PacBio, Nanopore) (fasta or fastq)\n"
              << "    -x PAIR1 [PAIR2 ...]               : linked-reads files (paired-ends, 10x Genomics) (interleaved, fasta or fastq)\n"
              << "    -X FWD1 REV1 [FWD2 REV2 ...]       : linked-reads files (paired-ends, 10x Genomics) (separate forward and reverse files, fasta or fastq)\n"
              << "    -i INT                             : number of iterations (default " << optionSingleArgs.at("-i") << ")\n"
              << "    -l INT                             : minimum number of links to scaffold (default " << optionSingleArgs.at("-l") << ")\n"
              << "    -j INT                             : minimum number of links to phase variants (default " << optionSingleArgs.at("-k") << ")\n"
              << "    -t INT                             : number of threads (default " << optionSingleArgs.at("-t") << ")\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n"
              << "    -mapper FILE                       : path of mapper executable file (default minimap2, only effective with -p option)\n"
//              << "    -minialign                         : use minialign insterd of minimap2 (default) for alignment of long reads\n"
//              << "    -kmer_align                        : use built-in alignment based on k-mer match insterd of minimap (default) for long reads\n"
              << "\n\n"
              << "Outputs:\n"
              << "    PREFIX_allPhaseBlock.fa\n"
              << "    PREFIX_primaryBubble.fa\n"
              << "    PREFIX_secondaryBubble.fa\n"
              << "    PREFIX_nonBubbleHomoCandidate.fa\n"
              << "    PREFIX_nonBubbleHetero.fa\n"
              << "    PREFIX_intermediateResults (directory)\n"
              << "\n"
              << "Note, PREFIX is specified by -o\n"
              << std::endl;
}


void Phase::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());

	string intermediateDirectoryName = optionSingleArgs["-o"];
	intermediateDirectoryName += "_intermediateResults";
	this->setDirectoryName(intermediateDirectoryName);
	this->createDirectory();

    unsigned long numIterate = atoi(optionSingleArgs["-i"].c_str());
    if (!(optionMultiArgs["-p"].empty()) || !(optionMultiArgs["-x"].empty()) || !(optionMultiArgs["-X"].empty())) {
		++numIterate;
	}
		
    for (unsigned long times = 1; times <= numIterate; ++times) {
		std::ostringstream oss;
		oss << intermediateDirectoryName << "/round" << times;
		this->setDirectoryName(oss.str());

        this->createDirectory();
		this->execSolveDBG(times, numIterate);
		this->execGapClose();
    }

	this->moveAndConcatenateFinalRoundResult(intermediateDirectoryName, numIterate);

	std::cerr << "phase completed!" << std::endl;
}


void Phase::setDirectoryName(const string newDirectoryName)
{
	this->previousDirectoryName = this->directoryName;
	this->directoryName = newDirectoryName;
}


void Phase::createDirectory(void) const
{
    struct stat st;

    if (stat(directoryName.c_str(), &st) != 0) {
        int returnValue = mkdir(directoryName.c_str(), 0755);
        if (returnValue != 0) {
            throw platanus::CreateDirError(directoryName);
        }
    }
}


void Phase::execSolveDBG(const unsigned long times, const unsigned long numIterate)
{
	vector<string> singleArgOptionList = {"-k", "-t", "-tmp", "-l", "-j", "-mapper"};
	vector<string> multiArgOptionList = {"-x", "-X"};
	if (times > 1) {
		multiArgOptionList.push_back("-p");
	}

	string inputPrefix = this->previousDirectoryName;
	inputPrefix += "/";
	inputPrefix += optionSingleArgs["-o"];

	string outputPrefix = this->directoryName;
	outputPrefix += "/";
	outputPrefix += optionSingleArgs["-o"];


    std::ostringstream oss;
    oss << this->platanusBinary << " solve_DBG";

	if (times == 1) {
		multiArgOptionList.push_back("-c");
		multiArgOptionList.push_back("-b");
	}
	else {
        oss << " -c " << inputPrefix << "_gapClosed_nonBubbleOther.fa"
			<< " -b " << inputPrefix << "_gapClosed_primaryBubble.fa" << " " << inputPrefix << "_gapClosed_secondaryBubble.fa";
	}

	if (times != numIterate) {
		oss << " -no_barcode_scaffold";
	}

    for (auto itr = optionPairFile.begin(); itr != optionPairFile.end(); ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (!(itr->fileSecond.empty())) {
            oss << " " << itr->fileSecond;
        }
    }

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

	oss << " -o " << outputPrefix;
    oss << " 2>" << outputPrefix << ".solveDBGLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::SolveDBGError();
    }
}


void Phase::execGapClose(void)
{
	string inputPrefix = this->directoryName;
	inputPrefix += "/";
	inputPrefix += optionSingleArgs["-o"];

	string outputPrefix = this->directoryName;
	outputPrefix += "/";
	outputPrefix += optionSingleArgs["-o"];

	vector<string> singleArgOptionList = {"-t", "-tmp"};


    std::ostringstream oss;
    oss << this->platanusBinary << " gap_close";

	oss << " -c"
		<< " " << inputPrefix << "_primaryBubble.fa"
		<< " " << inputPrefix << "_secondaryBubble.fa"
		<< " " << inputPrefix << "_nonBubbleOther.fa"
		<< " " << inputPrefix << "_nonBubbleHetero.fa";

	for (auto optItr = singleArgOptionList.begin(); optItr != singleArgOptionList.end(); ++optItr) {
		if (!(optionSingleArgs[*optItr].empty())) {
			oss << " " << *optItr << " " << optionSingleArgs[*optItr];
		}
	}

    for (auto itr = optionPairFile.begin(); itr != optionPairFile.end(); ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (!(itr->fileSecond.empty())) {
            oss << " " << itr->fileSecond;
        }
    }

	oss << " -o " << outputPrefix;
    oss << " 2>" << outputPrefix << ".gapCloseLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::GapError();
    }
}

void Phase::moveAndConcatenateFinalRoundResult(const string intermediateDirectoryName, const unsigned long times)
{
    std::ostringstream prefixOss;
	prefixOss << intermediateDirectoryName << "/round" << times << "/" << optionSingleArgs["-o"];
	string finalRoundPrefix = prefixOss.str();

    std::ostringstream cmdOss;
    cmdOss << "mv " << finalRoundPrefix << "_gapClosed_primaryBubble.fa" << " " << optionSingleArgs["-o"] << "_primaryBubble.fa" << " && "
		<< "mv " << finalRoundPrefix << "_gapClosed_secondaryBubble.fa" << " " << optionSingleArgs["-o"] << "_secondaryBubble.fa" << " && "
		<< "mv " << finalRoundPrefix << "_gapClosed_nonBubbleHetero.fa" << " " << optionSingleArgs["-o"] << "_nonBubbleHetero.fa" << " && "
		<< "mv " << finalRoundPrefix << "_gapClosed_nonBubbleOther.fa" << " " << optionSingleArgs["-o"] << "_nonBubbleHomoCandidate.fa" << " && "

		<< "cat " << optionSingleArgs["-o"] << "_primaryBubble.fa" << " " 
		<< optionSingleArgs["-o"] << "_secondaryBubble.fa" << " "
		<< optionSingleArgs["-o"] << "_nonBubbleHetero.fa" << " "
		<< optionSingleArgs["-o"] << "_nonBubbleHomoCandidate.fa" << " "
		<< ">" << optionSingleArgs["-o"] << "_allPhaseBlock.fa";

    if (system(cmdOss.str().c_str()) != 0) {
        throw platanus::CreateLinkError();
    }
}
