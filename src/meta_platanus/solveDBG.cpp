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

#include "solveDBG.h"
#include "seqlib.h"
#include "kmer.h"
#include <cmath>
#include <algorithm>
#include <sstream>
#include <climits>
#include <cfloat>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>

using std::vector;
using std::string;
using std::unordered_map;
using std::cerr;
using std::endl;


//////////////////////////////////////////////////////////////////////////////////////
// const parameter define
//////////////////////////////////////////////////////////////////////////////////////
const double SolveDBG::MIN_LONG_READ_LENGTH_CUTOFF_FACTOR = 1;
const double SolveDBG::MAX_LONG_READ_LENGTH_CUTOFF_FACTOR = 8;
const double SolveDBG::LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR = 0.1;
const long SolveDBG::LONG_READ_MIN_ALIGNMENT_LENGTH = 1000;
const double SolveDBG::LONG_READ_MIN_ALIGNMENT_COVERAGE = 0.8;
const double SolveDBG::LONG_READ_MIN_ALIGNMENT_IDENTITY = 0.8;

const long SolveDBG::MIN_TAG_SCAFFOLD_LENGTH_CUTOFF_FACTOR = 1;
const long SolveDBG::MAX_TAG_SCAFFOLD_LENGTH_CUTOFF_FACTOR = 2;
const double SolveDBG::TAG_SCAFFOLD_LENGTH_CUTOFF_UNIT_FACTOR = 0.5;
const double SolveDBG::TAG_SCAFFOLD_LENGTH_LOWER_CUTOFF_UNIT_FACTOR = 0.25;

const long SolveDBG::MIN_TAG_SCAFFOLD_ISLAND_LENGTH_CUTOFF_FACTOR = 1;
const long SolveDBG::MAX_TAG_SCAFFOLD_ISLAND_LENGTH_CUTOFF_FACTOR = 2;
const long SolveDBG::MIN_NUM_FOR_TAG_SCAFFOLD_ISLAND = 3;
const double SolveDBG::TAG_SCAFFOLD_ISLAND_LENGTH_CUTOFF_UNIT_FACTOR = 0.1;

//////////////////////////////////////////////////////////////////////////////////////
// default constructor
//////////////////////////////////////////////////////////////////////////////////////
SolveDBG::SolveDBG()
: Scaffold()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-e"] = "";
    optionSingleArgs["-v"] = "32";
    optionSingleArgs["-r"] = "0.05";
    optionSingleArgs["-l"] = "2";
    optionSingleArgs["-j"] = "1";
    optionSingleArgs["-J"] = "2";
    optionSingleArgs["-u"] = "0";
    optionSingleArgs["-L"] = "200000";
    optionSingleArgs["-t"] = "1";
    optionSingleArgs["-k"] = "";
    optionSingleArgs["-g"] = "";
    optionMultiArgs["-c"] = vector<string>();
    optionMultiArgs["-b"] = vector<string>();
    optionMultiArgs["-p"] = vector<string>();
    optionMultiArgs["-ont"] = vector<string>();
    optionMultiArgs["-gc"] = vector<string>();
    optionMultiArgs["-x"] = vector<string>();
    optionMultiArgs["-X"] = vector<string>();

//    optionMultiArgs["-s"] = vector<string>(3);
//    optionMultiArgs["-s"][0] = "32";
//    optionMultiArgs["-s"][1] = "64";
//    optionMultiArgs["-s"][2] = "96";
    optionMultiArgs["-s"] = vector<string>(1);
    optionMultiArgs["-s"][0] = "32";

    optionMultiArgs["-S"] = vector<string>(1);
    optionMultiArgs["-S"][0] = "20";


    optionBool["-no_scaffold"] = false;
    optionBool["-no_barcode_scaffold"] = false;
    optionBool["-unphase"] = false;
    optionBool["-trim_overlap"] = false;
    optionBool["-trim_cluster"] = false;

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
void SolveDBG::usage(void) const
{
    std::cerr << "\nUsage: platanus2 solveDBG [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : contig_file (fasta format)\n"
              << "    -k FILE                            : k-mer-occurrence file (binary format)\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -p FILE1 [FILE2 ...]               : PacBio long-read file (fasta or fastq)\n"
              << "    -ont FILE1 [FILE2 ...]             : Oxford Nanopore long-read file (fasta or fastq)\n"
              << "    -gc FILE1 [FILE2 ...]              : Guiding contig file; i.e. other assemblies, synthetic long-reads or corrected reads (fasta or fastq)\n"
              << "    -x PAIR1 [PAIR2 ...]               : barcoded_pair_files (10x Genomics) (reads in 1 file, fasta or fastq)\n"
              << "    -X FWD1 REV1 [FWD2 REV2 ...]       : barcoded_pair_files (10x Genomics) (reads in 2 files, fasta or fastq)\n"
              << "    -n{INT} INT                        : lib_id minimum_insert_size\n"
              << "    -a{INT} INT                        : lib_id average_insert_size\n"
              << "    -d{INT} INT                        : lib_id SD_insert_size\n"
              << "    -r FLOAT                           : threshold of conditional probability in edge-cut (default " << optionSingleArgs.at("-r") << ")\n"
              << "    -L INT                             : maximum fragment length of tag (10x Genomics) (default " << optionSingleArgs.at("-L") << ")\n"
              << "    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default " << optionMultiArgs.at("-s")[0] << ")\n"
              << "    -j INT                             : minimum number of links to phase variants (default " << optionSingleArgs.at("-j") << ")\n"
              << "    -J INT                             : minimum number of links to barcode scaffolding (default " << optionSingleArgs.at("-J") << ")\n"
              << "    -l INT                             : minimum number of links to scaffold (default " << optionSingleArgs.at("-l") << ")\n"
              << "    -u FLOAT                           : maximum difference for bubble crush (identity, default " << optionSingleArgs.at("-u") << ")\n"
              << "    -t INT                             : number of threads (<= " << optionSingleArgs.at("-t") << ", default 1)\n"
              << "    -mapper FILE                       : path of mapper executable file (default, minimap, only effective with -p option)\n"
              << "    -unphase                           : not phase heterozygous regions and construct consensus scaffolds (default false)\n"
              << "    -trim_overlap                      : trim overlapping edges of scaffolds (default, off)\n"
              << "    -trim_cluster                      : only trim overlapping edges of scaffolds within each cluster (default, off)\n"
              << "    -minimap2                          : use minimap2 insterd of minimap (default) for alignment of long reads\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n"
			  << "\n\n"

              << "Outputs:\n"
			  << "    PREFIX_*.fa\n"
             << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// initialize parameters
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::initializeParameters(void)
{
	seedLength = platanus::ConstParam::SCAFFOLD_HASH_OVERLAP;
    multiSeedLengthForShortRead.clear();
	if (!(optionPairFile.empty())) {
		for (auto itr = optionMultiArgs["-s"].begin(); itr != optionMultiArgs["-s"].end(); ++itr) {
			unsigned length = std::stoi(*itr);
			multiSeedLengthForShortRead.push_back(length);
			if (seedLength > length)
				seedLength = length;
		}
	}

    keyLength = std::min(seedLength, platanus::ConstParam::SCAFFOLD_HASH_OVERLAP);
    bubbleThreshold = atof(optionSingleArgs["-u"].c_str());
    minLink = atoi(optionSingleArgs["-l"].c_str());
    minLinkToPhase = atoi(optionSingleArgs["-j"].c_str());
    minLinkTagScaffold = atoi(optionSingleArgs["-J"].c_str());
    minOverlapForScaffolding = atoi(optionSingleArgs["-v"].c_str());
    numThread = atoi(optionSingleArgs["-t"].c_str());
    pairedDBG.setSeedLength(seedLength);
    pairedDBG.setMinTolerenceFactor(MIN_TOL_FACTOR);
    pairedDBG.setMaxFragmentLengthOfTag(atoi(optionSingleArgs["-L"].c_str()));
	pairedDBG.setUondProbThreshold(atof(optionSingleArgs["-r"].c_str()));

    sort(optionPairFile.begin(), optionPairFile.end());
    numFilePerLibraryID.resize(optionPairFile.size());
    libraryIDList.resize(optionPairFile.size());
    unsigned numLibrary = 0;
    for (unsigned i = 0; i < optionPairFile.size(); ++i) {
        ++(numFilePerLibraryID[numLibrary]);
        libraryIDList[numLibrary] = optionPairFile[i].libraryID;
        if (optionPairFile[i].libraryID != optionPairFile[i + 1].libraryID) {
            ++numLibrary;
        }
    }

    libraryMT.resize(numLibrary);
    omp_set_num_threads(numThread);

	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
}


//////////////////////////////////////////////////////////////////////////////////////
// exec scaffold
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::exec(void)
{
    initializeParameters();
    mapLibraryAndInitGraph(numThread);


    if (optionBool["-trim_cluster"]) {
		pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, true, true);
		pairedDBG.markRedundantResultSeq(numThread, true);

		vector<vector<string > > outputSeq;
		vector<vector<string> > outputName;
		pairedDBG.splitResultSeqForClusters(outputSeq, outputName);
		printClusteredSeq(outputSeq, outputName, optionSingleArgs["-o"] + "_trimmedClusters");
		pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_trimmedClustersComponent.bed");

		cerr << "solve_DBG completed!" << endl;
		return;
	}


//	unsigned correctionMode;
	if (!clusterFlag) {
//		correctionMode = PairedDBG::INTRA_OTU_MODE;
		pairedDBG.setMgcEmpDistr();
		pairedDBG.setCoverageEmpDistr(optionSingleArgs["-k"]);
	}
	else {
//		correctionMode = PairedDBG::CLUSTERING_MODE;
	}

	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.makeGraph(numThread);


	if (optionBool["-unphase"]) {
		pairedDBG.calculateHeteroAndAverageCoverageUnphase();
		pairedDBG.clearEdges();

		extendConsensus(false, true, true);

		if (!(libraryMT.empty()))
			pairedDBG.setTolerence(MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
		else
			pairedDBG.setTolerence(this->contigMaxK);
		pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, optionBool["-trim_overlap"], false);
		pairedDBG.markRedundantResultSeq(numThread);

		if (!clusterFlag) {
			outputGraph("_consensusScaffold.fa");
			pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_consensusScaffoldComponent.bed");
		}
		else {
			vector<vector<string > > outputSeq;
			vector<vector<string> > outputName;
			pairedDBG.splitResultSeqForClusters(outputSeq, outputName);
			printClusteredSeq(outputSeq, outputName, optionSingleArgs["-o"] + "_extendedClusters");

			outputGraph("_postClusteringScaffold.fa");
			pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_extendedClustersComponent.bed");
		}

		cerr << "solve_DBG completed!" << endl;
		return;
	}


	pairedDBG.extractDBGBubbleInformation();
	pairedDBG.clearEdges();

	extendConsensus(true, true, false);
	pairedDBG.resetGraph();


	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.setCutoffLength(0);
	pairedDBG.makeGraph(numThread);
	pairedDBG.extractDBGBubbleInformation();
	pairedDBG.setOppositeBubbleContigIDByEndMatch();
	pairedDBG.setOppositeBubbleContigIDOverlapped(numThread);
	pairedDBG.setBubbleJunctionContigIDOverlapped();
	pairedDBG.clearEdges();


	for (long outerIteration = 0; outerIteration < 2; ++outerIteration) {
		for (long iteration = 0; iteration < 2; ++iteration) {
			if (iteration == 0)
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE);
			else
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
			pairedDBG.setCutoffLength(0);

			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

				pairedDBG.setTargetLibraryIndex(i);
				unsigned tolerenceFactor = MAX_TOL_FACTOR;
				pairedDBG.setTolerence(tolerenceFactor * libraryMT[i][0].getSDInsSize());
				cerr << "TOLERENCE_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);
			}
			if (longReadLibraryMT.size() > 0) {
				if (iteration == 0)
					pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE);
				else
					pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
				cerr << "[LONG_READ_LIBRARY]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageInsSize() << ", SD_READ_LENGTH = " << longReadLibraryMT[0].getSDInsSize() << endl;
				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::SCORE, clusterFlag, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);
			}

			if (iteration == 0)
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE);
			else
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				pairedDBG.setTolerenceFactor(MAX_TOL_FACTOR);
				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.makeGraphAllLibraries(numThread);
				if (clusterFlag)
					pairedDBG.deleteChimericEdge();
				else
					pairedDBG.deleteInterOTUEdge();
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);
			}
			if (optionBool["-no_scaffold"]) {
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
				pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
			}
		}

//		pairedDBG.solveSimpleCrossStructureCoverageIterative(clusterFlag, numThread);

		if (optionBool["-no_scaffold"]) {
			pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE);
			pairedDBG.solveUniquePathBetweenLinkedNodePairAllLibrariesIterative(numThread);

			pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
			pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
			pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
			continue;
		}

		pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
		pairedDBG.setCutoffLength(0);
		pairedDBG.setMinLink(minLink);
		pairedDBG.makeGraph(numThread);
		pairedDBG.clearEdges();


		for (long iteration = 0; iteration < 2; ++iteration) {
			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				if (iteration == 0)
					pairedDBG.setMinLink(minLinkToPhase);
				else
					pairedDBG.setMinLink(minLink);

				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

				pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE);

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;
					cerr << "LENGTH_CUTOFF = " << pairedDBG.getCutoffLength() << endl;

					pairedDBG.trimSparseEnd();
					if (iteration > 0) {
						if (clusterFlag)
							pairedDBG.deleteChimericEdge();
						else
							pairedDBG.deleteInterOTUEdge();

						pairedDBG.joinUnambiguousNodePairGappedIterative(numThread);
					}
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);
				}
			}
			if (longReadLibraryMT.size() > 0) {
				pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE);
				cerr << "[LONG_READ_LIBRARY]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageLength() << endl;
				pairedDBG.setTolerence(2 * this->contigMaxK);
				if (iteration == 0)
					pairedDBG.setMinLink(minLinkToPhase);
				else
					pairedDBG.setMinLink(minLink);

				for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
					pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR * longReadLibraryMT[0].getAverageInsSize());
					pairedDBG.trimSparseEnd();
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::SCORE, clusterFlag, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);
				}
			}
			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE);

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setTolerenceFactor(tolerenceFactor);
					pairedDBG.setCutoffLengthFactor(MIN_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					if (iteration == 0)
						pairedDBG.setMinLink(minLinkToPhase);
					else {
						pairedDBG.setMinLink(minLink);
						if (clusterFlag)
							pairedDBG.deleteChimericEdge();
						else
							pairedDBG.deleteInterOTUEdge();
						pairedDBG.joinUnambiguousNodePairGappedIterativeAllLibraries(numThread);
					}
					pairedDBG.trimSparseEnd();
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);
				}
			}
		}


		pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE);
		pairedDBG.solveUniquePathBetweenLinkedNodePairAllLibrariesIterative(numThread);


		pairedDBG.trimSparseEnd();
		pairedDBG.divideNodeUsingBubbleContigPair(numThread);

		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
		pairedDBG.divideBubbleContigInNonHeteroNode();
		pairedDBG.divideBubbleJunctionNode(false);

		pairedDBG.setMinOverlap(minOverlapForScaffolding);
		for (long iteration = 0; iteration < 2; ++iteration) {
			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				long linkThreshold;
				if (iteration % 2 == 0)
					linkThreshold = minLink;
				else
					linkThreshold = std::max(minLink, pairedDBG.estimateLink());

				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;
					cerr << "LENGTH_CUTOFF = " << pairedDBG.getCutoffLength() << endl;

					pairedDBG.trimSparseEnd();

					pairedDBG.setMinLink(linkThreshold);
					pairedDBG.makeGraph(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					if (clusterFlag)
						pairedDBG.deleteChimericEdge();
					else
						pairedDBG.deleteInterOTUEdge();
					pairedDBG.makeScaffold();

					pairedDBG.joinUnambiguousNodePairGappedIterative(numThread);
					pairedDBG.setMinLink(minLink);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);

					pairedDBG.setMinLink(linkThreshold);
					pairedDBG.makeGraph(numThread);
					pairedDBG.deleteErroneousEdgeNumTagRateIterative(false, numThread);
					pairedDBG.deleteDifferentBubbleEdgeIterative(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					if (clusterFlag)
						pairedDBG.deleteChimericEdge();
					else
						pairedDBG.deleteInterOTUEdge();
					pairedDBG.makeScaffold();
				}
			}
			if (outerIteration < 1)
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::GAP, numThread, this->contigMaxK);
			else
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		}
		pairedDBG.setMinOverlap(this->contigMaxK - 1);


		pairedDBG.trimSparseEnd();
		pairedDBG.divideNodeUsingBubbleContigPair(numThread);
		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		pairedDBG.divideNodeBasedOnBubblesIterative(false, numThread);


		pairedDBG.setMinOverlap(minOverlapForScaffolding);
		for (long iteration = 0; iteration < 2; ++iteration) {
			long linkThreashold = (iteration + 1) * minLink;
			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setTolerenceFactor(tolerenceFactor);
					pairedDBG.setCutoffLengthFactor(MAX_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);

					pairedDBG.trimSparseEnd();

					pairedDBG.setMinLink(linkThreashold);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					pairedDBG.makeGraphAllLibraries(numThread);
					pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					if (clusterFlag)
						pairedDBG.deleteChimericEdge();
					else
						pairedDBG.deleteInterOTUEdge();
					pairedDBG.makeScaffold();

					pairedDBG.setMinLink(minLink);
					pairedDBG.joinUnambiguousNodePairGappedIterativeAllLibraries(numThread);
					pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);

					pairedDBG.setMinLink(linkThreashold);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					pairedDBG.makeGraphAllLibraries(numThread);
					pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
					pairedDBG.deleteErroneousEdgeNumTagRateIterative(false, numThread);
					pairedDBG.deleteDifferentBubbleEdgeIterative(numThread);
					pairedDBG.deleteConflictingBubbleEdge(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					if (clusterFlag)
						pairedDBG.deleteChimericEdge();
					else
						pairedDBG.deleteInterOTUEdge();
					pairedDBG.makeScaffold();
				}
			}
			if (outerIteration < 1)
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::GAP, numThread, this->contigMaxK);
			else
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		}
		pairedDBG.setMinOverlap(this->contigMaxK - 1);


		if (longReadLibraryMT.size() > 0) {
			pairedDBG.setMinOverlap(minOverlapForScaffolding);

			pairedDBG.divideNodeUsingBubbleContigPair(numThread);

			for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
				if (outerIteration == 0)
					pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);
				else
					pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

				pairedDBG.setTolerence(2 * this->contigMaxK);
				pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR *  longReadLibraryMT[0].getAverageInsSize());

				pairedDBG.trimSparseEnd();

				pairedDBG.setMinLink(minLink);
				pairedDBG.makeGraph(numThread);
				pairedDBG.deleteErroneousEdgeNumLinkRateIterative(numThread);
				pairedDBG.deleteRepeatEdge();
				pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
				if (clusterFlag)
					pairedDBG.deleteChimericEdge();
				else
					pairedDBG.deleteInterOTUEdge();
				pairedDBG.makeScaffold();

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);
				pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.makeGraph(numThread);
				pairedDBG.deleteErroneousEdgeNumTagRateIterative(false, numThread);
				pairedDBG.deleteDifferentBubbleEdgeIterative(numThread);
				pairedDBG.deleteErroneousEdgeNumLinkRateIterative(numThread);
				pairedDBG.deleteRepeatEdge();
				pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
				if (clusterFlag)
					pairedDBG.deleteChimericEdge();
				else
					pairedDBG.deleteInterOTUEdge();
				pairedDBG.makeScaffold();

				pairedDBG.divideNodeUsingBubbleContigPair(numThread);
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE);
				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.setCutoffLength(0);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);
			}

			pairedDBG.setMinOverlap(this->contigMaxK - 1);
		}


		if (tagLibraryMT.size() > 0 && !(optionBool["-no_scaffold"]) && !(optionBool["-no_barcode_scaffold"])) {
			pairedDBG.setMinOverlap(minOverlapForScaffolding);
			pairedDBG.setMode(PairedDBG::TAG_SCAFFOLD_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

			for (unsigned cutoffFactor = MIN_TAG_SCAFFOLD_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_TAG_SCAFFOLD_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
//				pairedDBG.setTolerence(TAG_SCAFFOLD_LENGTH_CUTOFF_UNIT_FACTOR * tagLibraryMT[0].getAverageInsSize());
				pairedDBG.setTolerence(tagLibraryMT[0].getAverageInsSize());
				pairedDBG.setCutoffLength(cutoffFactor * TAG_SCAFFOLD_LENGTH_CUTOFF_UNIT_FACTOR *  tagLibraryMT[0].getAverageInsSize());

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);
				pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.makeGraph(numThread);
				pairedDBG.deleteConflictingEdgeToSameNode(numThread);
				pairedDBG.deleteErroneousEdgeNumTagRateIterative(false, numThread);
				pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
				if (clusterFlag)
					pairedDBG.deleteChimericEdge();
				else
					pairedDBG.deleteInterOTUEdge();
				pairedDBG.makeScaffold();

				if (outerIteration < 1)
					pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::GAP, numThread, this->contigMaxK);
				else
					pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
			}

			pairedDBG.setMinOverlap(this->contigMaxK - 1);
		}


		pairedDBG.trimSparseEnd();
		pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
		pairedDBG.divideBubbleJunctionNode(true);
		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);

		if (outerIteration < 1) {
			pairedDBG.joinUnambiguousNodePairIterative(numThread);
			pairedDBG.divideGappedNode((outerIteration + 1) * this->contigMaxK);
		}
	}

	pairedDBG.divideBubbleContigInNonHeteroNode();

	pairedDBG.setMode(0);
	pairedDBG.makeGraph(numThread);
	pairedDBG.adjustOppositeBubbleNodeIDDirection();
	pairedDBG.copyAllNodes(phasedGraph);

	pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
	pairedDBG.clearEdges();
	pairedDBG.makeScaffold();

	extendConsensus(false, false, true);

	if (!(libraryMT.empty()))
		pairedDBG.setTolerence(MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
	else
		pairedDBG.setTolerence(this->contigMaxK);
    pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, false, false);
	outputGraph("_preliminaryConsensusScaffold.fa");

	pairedDBG.setMode(0);
	pairedDBG.remakeGraphRecoveringSecondaryBubble(phasedGraph);
	pairedDBG.makeGraph(numThread);
	pairedDBG.divideNodeBasedOnBubblesIterative(false, numThread);
	pairedDBG.makeGraph(numThread);
	pairedDBG.adjustOppositeBubbleNodeIDDirection();
	if (!(libraryMT.empty()))
		pairedDBG.setTolerence(MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
	else
		pairedDBG.setTolerence(this->contigMaxK);
    pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, optionBool["-trim_overlap"], false);
	pairedDBG.markRedundantResultSeq(numThread);
	pairedDBG.outputResultSeqWithBubble(optionSingleArgs["-o"], "_primaryBubble.fa", "_secondaryBubble.fa", "_nonBubbleHetero.fa", "_nonBubbleOther.fa", "_bubbleRelation.tsv", this->contigMaxK, this->contigReadLength);

	pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_phasedScaffoldComponent.bed");
	
	cerr << "solve_DBG completed!" << endl;
	return;
}



//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::mapLibraryAndInitGraph(const int numThread)
{
    platanus::Contig contig;

    std::unique_ptr<HeteroMapper> mapper(new HeteroMapper(seedLength, keyLength));

    readLibrary(mapper, contig, numThread);
    cerr << "CONTIG_AVERAGE_COVERAGE = " << averageCoverage << endl;

	mapper->setMultiSeedLength(multiSeedLengthForShortRead);
    unsigned nowFileNumber = 0;
    for (unsigned i = 0; i < libraryMT.size(); ++i) {
		libraryMT[i][0].setAverageInsSize(0);
        int nowLibraryID = optionPairFile[nowFileNumber].libraryID;
        cerr << "[LIBRARY " << libraryIDList[i] << "]" << endl;
        // set average length and minimum insert size
        // estimate insert size
        libraryMT[i][0].setAverageLength((long)((double)(libraryMT[i][0].getTotalLength()) / (2 * libraryMT[i][0].getNumPair()) + 0.5));
        long minInsertion = optionMinIns.find(nowLibraryID) == optionMinIns.end() ? 0 : optionMinIns[nowLibraryID];

        mapper->contigMap.mapPairAndSaveReadLink(libraryMT[i], minInsertion, this->contigMaxK + 1, numThread);

        libraryMT[i][0].setInsCutoffRate(optionInsCutoffRate.find(nowLibraryID) == optionInsCutoffRate.end() ? DEFAULT_INS_CUTOFF_RATE : optionInsCutoffRate[nowLibraryID]);
        if (optionAveIns.find(nowLibraryID) != optionAveIns.end() || optionSDIns.find(nowLibraryID) != optionSDIns.end()) {
            if (optionAveIns.find(nowLibraryID) != optionAveIns.end()) {
                libraryMT[i][0].setAverageInsSize(optionAveIns[nowLibraryID]);
                std::cerr << "Average insert size specified: AVE = " << libraryMT[i][0].getAverageInsSize() << std::endl;
            }
            if (optionSDIns.find(nowLibraryID) != optionSDIns.end()) {
                libraryMT[i][0].setSDInsSize(optionSDIns[nowLibraryID]);
            } else {
                libraryMT[i][0].setSDInsSize(static_cast<long>(static_cast<double>(libraryMT[i][0].getAverageInsSize()) / 10.0 + 0.5));
            }
        }
        nowFileNumber += numFilePerLibraryID[i];
    }

	if (longReadLibraryMT.size() > 0) {
		cerr << "[LONG_READ LIBRARY]" << endl;

		string alignerOutFilename(optionSingleArgs["-o"]);
		alignerOutFilename += "_longReadAlignment.tsv";

		string aligner;

		if (optionSingleArgs["-mapper"].empty())
			aligner = "minimap2";

		vector<string> targetFilename = optionMultiArgs["-c"];
		execMinimap2(targetFilename, optionMultiArgs["-b"], pacBioLongReadFilename, alignerOutFilename, numThread, aligner, "-p 1 -x map-pb");
		execMinimap2(targetFilename, optionMultiArgs["-b"], nanoporeLongReadFilename, alignerOutFilename, numThread, aligner, "-p 1 -x map-ont");
		execMinimap2(targetFilename, optionMultiArgs["-b"], guideContigLongReadFilename, alignerOutFilename, numThread, aligner, "-p 1 -x asm10");

		longReadLibraryMT[0].setAverageLength((long)((double)(longReadLibraryMT[0].getTotalLength()) / (2 * longReadLibraryMT[0].getNumPair()) + 0.5));
		mapper->contigMap.readLongReadPAFfileAndSaveLink(alignerOutFilename, longReadLibraryMT, LONG_READ_MIN_ALIGNMENT_LENGTH, LONG_READ_MIN_ALIGNMENT_COVERAGE, LONG_READ_MIN_ALIGNMENT_IDENTITY, contigMaxK, numThread);

		vector<long> insSizeDistribution;
		longReadLibraryMT[0].readInsertSizeFile(insSizeDistribution);
		longReadLibraryMT[0].estimateInsSize(insSizeDistribution, 0, platanus::ConstParam::LONG_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::LONG_READ_INS_SIZE_UPPER_BOUND_FACTOR);

		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_longReadLibrary" << "_readDistribution.tsv";
		longReadLibraryMT[0].printInsertSizeFreq(outStream.str());
		cerr << "[LONG_READ_LIBRARY " << 1 << "]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageInsSize()
			 << ", SD_READ_LENGTH = " << longReadLibraryMT[0].getSDInsSize() << endl;
	}	

	if (tagLibraryMT.size() > 0) {
		mapper->setMultiSeedLength(multiSeedLengthForShortRead);

		cerr << "[TAG LIBRARY]" << endl;
		mapper->contigMap.mapTagPairMT(tagLibraryMT, numThread);
	}	

	if (libraryMT.size() > 0) {
		pairedDBG.setAllLibraryMT(&libraryMT);
		pairedDBG.setTargetLibraryIndex(0);
		pairedDBG.initScaffolding(contig.coverage, mapper->contigMap, averageCoverage, bubbleThreshold);
		pairedDBG.setContigName(contig.name);

		vector<long> insSizeDistribution;
		libraryMT[0][0].readInsertSizeFile(insSizeDistribution);
		pairedDBG.insertSizeDistribution(libraryMT[0], insSizeDistribution, numThread);

		vector<long> seqLengths;
		pairedDBG.scaffoldLengthList(seqLengths);

		if (libraryMT[0][0].getAverageInsSize() == 0)
			libraryMT[0][0].estimateInsSize(insSizeDistribution, 0, platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);

		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_lib" << 1 << "_insFreq.tsv";
		libraryMT[0][0].printInsertSizeFreq(outStream.str());
		cerr << "[LIBRARY " << 1 << "]\nAVE_INS = " << libraryMT[0][0].getAverageInsSize()
			 << ", SD_INS = " << libraryMT[0][0].getSDInsSize() << endl;
	}
	else {
		pairedDBG.initScaffolding(contig.coverage, mapper->contigMap, averageCoverage, bubbleThreshold);
		pairedDBG.setContigName(contig.name);
	}

	if (optionSingleArgs["-g"] != "")
		clusterFlag = pairedDBG.setClusterID(contig, optionSingleArgs["-g"]);
	else if (optionBool["-trim_cluster"])
		pairedDBG.setClusterIDfromSeqName(contig);

    pairedDBG.setContigMaxK(this->contigMaxK);
    pairedDBG.setMinOverlap(this->contigMaxK - 1);
//pairedDBG.setMinOverlap(32);
    pairedDBG.saveOverlap(mapper->contigMap, this->contigMaxK - 1, this->contigMaxK, numThread);
    pairedDBG.classifyNode();

	if (longReadLibraryMT.size() > 0) {
		pairedDBG.setLongReadLibraryMT(&longReadLibraryMT);
	}

	if (tagLibraryMT.size() > 0) {
		pairedDBG.setTagLibraryMT(&tagLibraryMT);
		pairedDBG.countMappedTagForEachContig(numThread);
	}

    cerr << "destructing mapper objects..." << std::endl;
}


void SolveDBG::readLibrary(std::unique_ptr<HeteroMapper> &mapper, platanus::Contig &contig, const int numThread)
{
    omp_set_num_threads(numThread);
    # pragma omp parallel for schedule(static, 1)
    for (int i = -3; i < static_cast<int>(libraryMT.size()); ++i) {
        try {
            // load contig file
            if (i == -3) {
                for (unsigned i = 0; i < optionMultiArgs["-c"].size(); ++i)
                    contig.readFastaCoverage(optionMultiArgs["-c"][i]);

				long j = contig.numSeq;
                for (unsigned i = 0; i < optionMultiArgs["-b"].size(); ++i)
                    contig.readFastaCoverage(optionMultiArgs["-b"][i]);

				pairedDBG.setNumInputBubbleContig(contig.numSeq - j);
				contig.setNameIndex();

                this->contigMaxK = contig.getMaxKFromFastaHeader(optionMultiArgs["-c"][0]);
                this->contigReadLength = contig.getReadLengthFromFastaHeader(optionMultiArgs["-c"][0]);
                if (this->contigReadLength == 0)
                    this->contigReadLength = platanus::ConstParam::DEFAULT_CONTIG_READ_LEN;

				if (optionSingleArgs["-e"] == "")
					averageCoverage = contig.calculateAverageCoverageExcludingOutlier(contig.getSeqLengthMedian());
				else
					averageCoverage = atof(optionSingleArgs["-e"].c_str());

                mapper->setContigMap(contig);
                mapper->makeKmerTableContigMap();
            }
			else if (i == -2) {
                if (optionMultiArgs["-p"].empty() && optionMultiArgs["-ont"].empty() && optionMultiArgs["-gc"].empty())
					continue;

				vector<string> filenames(optionMultiArgs["-p"]);
				pacBioLongReadFilename = optionMultiArgs["-p"];

				std::copy(optionMultiArgs["-ont"].begin(), optionMultiArgs["-ont"].end(), std::back_inserter(filenames));
				nanoporeLongReadFilename = optionMultiArgs["-ont"];

				std::copy(optionMultiArgs["-gc"].begin(), optionMultiArgs["-gc"].end(), std::back_inserter(filenames));
				guideContigLongReadFilename = optionMultiArgs["-gc"];

				longReadLibraryMT.resize(numThread);
				for (long j = 0; j < numThread; ++j)
					longReadLibraryMT[j].pairFP = platanus::makeTemporaryFile();

				for (unsigned j = 0; j < static_cast<long>(filenames.size()); ++j) {
					platanus::FILETYPE fileFormat = checkFileFormat(filenames[j]);
					bool isFastq;
					switch (fileFormat) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaSingleMT(longReadLibraryMT, filenames[j], numThread, false, isFastq, true);
				}
            }
			else if (i == -1) {
                if (optionMultiArgs["-x"].size() == 0 && optionMultiArgs["-X"].size() == 0)
					continue;

				vector<string> filenames(optionMultiArgs["-x"]);
				filenames.insert(filenames.end(), optionMultiArgs["-X"].begin(), optionMultiArgs["-X"].end());

				std::unordered_map<string, int> tagStringConverter;
				setTagStringConverter(filenames, tagStringConverter);

				tagLibraryMT.resize(numThread);
				for (long j = 0; j < numThread; ++j)
					tagLibraryMT[j].pairFP = platanus::makeTemporaryFile();

				for (unsigned j = 0; j < static_cast<long>(optionMultiArgs["-x"].size()); ++j) {
					platanus::FILETYPE fileFormat = checkFileFormat(optionMultiArgs["-x"][j]);
					bool isFastq;
					switch (fileFormat) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaSingleTaggedMT(tagLibraryMT, optionMultiArgs["-x"][j], numThread, false, isFastq, true, tagStringConverter);
				}

				for (unsigned j = 0; j < static_cast<long>(optionMultiArgs["-X"].size()); j += 2) {
					platanus::FILETYPE fileFormat1 = checkFileFormat(optionMultiArgs["-X"][j]);
					platanus::FILETYPE fileFormat2 = checkFileFormat(optionMultiArgs["-X"][j + 1]);
					if (fileFormat1 != fileFormat2) {
						throw platanus::FormatError("Different file type in paired-file (-X).");
					}
					bool isFastq;
					switch (fileFormat1) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaPairTaggedMT(tagLibraryMT, optionMultiArgs["-X"][j] , optionMultiArgs["-X"][j + 1], numThread, false, isFastq, tagStringConverter);
				}
            }
			else {
                unsigned nowFileNumber = 0;
                libraryMT[i].resize(numThread);
                for (int j = 0; j < i; ++j)
                    nowFileNumber += numFilePerLibraryID[j];

                for (int j = 0; j < numThread; ++j) {
                    libraryMT[i][j].makeTempPairFP();
                }
                for (int j = 0; j < numFilePerLibraryID[i]; ++j) {
                    bool isFastq, isMate;
                    isMate = optionPairFile[nowFileNumber+j].libraryType == "-op"
                           || optionPairFile[nowFileNumber+j].libraryType == "-OP" ? true : false;
                    if (optionPairFile[nowFileNumber+j].libraryType == "-ip" || optionPairFile[nowFileNumber+j].libraryType == "-op") {
                        platanus::FILETYPE fileFormat = checkFileFormat(optionPairFile[nowFileNumber+j].fileFirst);
                        switch (fileFormat) {
                        case platanus::UNKNOWN:
                            throw platanus::FormatError();
                            break;
                        case platanus::FASTA:
                            isFastq = false;
                            break;
                        default:
                            isFastq = true;
                            break;
                        }
                        ReadFastaSingleMT(libraryMT[i], optionPairFile[nowFileNumber+j].fileFirst, numThread, isMate, isFastq);
                    }
					else {
                        platanus::FILETYPE fileFormat1 = checkFileFormat(optionPairFile[nowFileNumber+j].fileFirst);
                        platanus::FILETYPE fileFormat2 = checkFileFormat(optionPairFile[nowFileNumber+j].fileSecond);
                        if (fileFormat1 != fileFormat2) {
                            throw platanus::FormatError("Different file type in paired-file.");
                        }
                        switch (fileFormat1) {
                        case platanus::UNKNOWN:
                            throw platanus::FormatError();
                            break;
                        case platanus::FASTA:
                            isFastq = false;
                            break;
                        default:
                            isFastq = true;
                            break;
                        }
                        ReadFastaPairMT(libraryMT[i], optionPairFile[nowFileNumber+j].fileFirst, optionPairFile[nowFileNumber+j].fileSecond, numThread, isMate, isFastq);
                    }
                }
            }
        } catch (platanus::ErrorBase &e) {
            e.showErrorMessage();
            exit(e.getID());
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// destroy graph
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::outputAndAfterTreatment(void)
{
    cerr << "writing scaffold files..." << endl;
    string outFilename(optionSingleArgs["-o"]);
    string componentFilename(outFilename);
    outFilename += "_solvedContig.fa";
    componentFilename += "_solvedContigComponent.tsv";
    pairedDBG.cutAndPrintSeq(this->contigMaxK, this->contigReadLength, outFilename, componentFilename);
}

void SolveDBG::outputGraph(const char *suffix)
{
    cerr << "writing scaffold files..." << endl;
    string outFilename(optionSingleArgs["-o"]);
    outFilename += suffix;
    pairedDBG.printResultSeq(outFilename);
}

void SolveDBG::updateAndWiteInsertSize(const long libraryIndex)
{
	pairedDBG.updateInsertLengthFP(libraryMT[libraryIndex], numThread);
	vector<long> insSizeDistribution;
	libraryMT[libraryIndex][0].readInsertSizeFile(insSizeDistribution);
	pairedDBG.insertSizeDistribution(libraryMT[libraryIndex], insSizeDistribution, numThread);
	if (libraryIndex > 0)
		libraryMT[libraryIndex][0].estimateInsSize(insSizeDistribution, libraryMT[libraryIndex - 1][0].getAverageInsSize(), platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);
	else
		libraryMT[libraryIndex][0].estimateInsSize(insSizeDistribution, 0, platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);

	std::ostringstream outStream;
	outStream << optionSingleArgs["-o"] << "_lib" << (libraryIndex + 1) << "_insFreq.tsv";
	printInsertSizeFreq(outStream.str(), insSizeDistribution);
}


void SolveDBG::execMinialign(const string targetFilename, const string &readFilename, const string &outFilename, const long numThread, const string minialignExecutable)
{
	std::ostringstream oss;

	oss << minialignExecutable <<  " -x pacbio -m 0 -O paf" << " -t " << numThread << " " << targetFilename << " " << readFilename << " >" << outFilename;

	std::cerr << "Executing minialign ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	cerr << "minialign fineshed." << endl;
}

void SolveDBG::execMinimap(const vector<string> &contigFilenames, const vector<string> &bubbleFilenames, const vector<string> &readFilenames, const string &outFilename, const long minAlignmentLength, const long numThread, const string minimapExecutable)
{
	std::ostringstream oss;

	oss << "cat ";
	for (auto itr = contigFilenames.begin(); itr != contigFilenames.end(); ++itr)
		oss << " " << *itr;
	for (auto itr = bubbleFilenames.begin(); itr != bubbleFilenames.end(); ++itr)
		oss << " " << *itr;

	oss << " | " << minimapExecutable <<  " -t " << numThread <<  " -L " << minAlignmentLength << " - ";

	for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr)
		oss << " " << *itr;
	oss << " | cut -f1-11 ";
	oss << " >" << outFilename;

	std::cerr << "Executing minimap ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	cerr << "minimap fineshed." << endl;
}

/*
void SolveDBG::execMinimap2(const vector<string> &contigFilenames, const vector<string> &bubbleFilenames, const vector<string> &readFilenames, const string &outFilename, const long numThread, const string minimap2Executable)
{
	std::ostringstream oss;

	oss << "cat ";
	for (auto itr = contigFilenames.begin(); itr != contigFilenames.end(); ++itr)
		oss << " " << *itr;
	for (auto itr = bubbleFilenames.begin(); itr != bubbleFilenames.end(); ++itr)
		oss << " " << *itr;

	oss << " | " << minimap2Executable << " -c " << " -t " << numThread << " - ";

	for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr)
		oss << " " << *itr;
	oss << " | cut -f1-11 ";
	oss << " >" << outFilename;

	std::cerr << "Executing minimap2 ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	cerr << "minimap2 fineshed." << endl;
}
*/
void SolveDBG::execMinimap2(const vector<string> &contigFilenames, const vector<string> &bubbleFilenames, const vector<string> &readFilenames, const string &outFilename, const long numThread, const string minimap2Executable, const string minimap2Option)
{
	if (readFilenames.empty())
		return;

	std::ostringstream oss;

	oss << "cat ";
	for (auto itr = contigFilenames.begin(); itr != contigFilenames.end(); ++itr)
		oss << " " << *itr;
	for (auto itr = bubbleFilenames.begin(); itr != bubbleFilenames.end(); ++itr)
		oss << " " << *itr;

	oss << " | " << minimap2Executable << " -c " << " -t " << numThread;
	if (optionBool["-minimap2_sensitive"])
		oss << " -p 0";
	
	oss << " " << minimap2Option << " - ";


	bool bzip2Flag = false;
	char bzip2TempFileName[platanus::globalTmpFileDir.size() + 7 + 1];
	platanus::FILECOMPRESSION format;

	for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr) {
		format = platanus::checkFileCompression(*itr);
		if (format == platanus::FILECOMPRESSION::BZIP2) {
			bzip2Flag = true;
			break;
		}
	}

	if (bzip2Flag) {
		strcpy(bzip2TempFileName, platanus::globalTmpFileDir.c_str());
		strcat(bzip2TempFileName, "/XXXXXX"); 

        int fd = mkstemp(bzip2TempFileName);
        if (fd == -1) {
            throw platanus::TMPError();
		}
        FILE *fp = fdopen(fd, "w+");
		fclose(fp);

		std::ostringstream bzip2Oss;
		bzip2Oss << "bzip2 -cd ";

		for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr) {
			format = platanus::checkFileCompression(*itr);
			if (format == platanus::FILECOMPRESSION::BZIP2)
				bzip2Oss << " " << *itr;
			else
				oss << " " << *itr;
		}

		bzip2Oss << " >"  << bzip2TempFileName;
		if (system(bzip2Oss.str().c_str()) != 0) {
			throw platanus::AlignerError();
		}

		oss << " " << bzip2TempFileName;
	}
	else {
		for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr)
			oss << " " << *itr;
	}

	oss << " | perl -pne \'s/cg:Z:\\S+//\' ";
	oss << " >>" << outFilename;

	std::cerr << "Executing minimap2 ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	if (bzip2Flag) {
        unlink(bzip2TempFileName);
	}

	cerr << "minimap2 finished." << endl;
}


void SolveDBG::outputForException(void)
{
    string outFilename(optionSingleArgs["-o"]);
    string componentFilename(outFilename);
    outFilename += "_scaffold.fa";
    this->printCatContigAsScaffold(outFilename, optionSingleArgs["-c"]);
    if (optionMultiArgs["-b"].size() != 0) {
        outFilename = optionSingleArgs["-o"] + "_scaffoldBubble.fa";
        this->printCatContigAsScaffold(outFilename, optionSingleArgs["-b"]);
    }
}


void SolveDBG::printCatContigAsScaffold(const std::string &outputFilename, const string &inputFilename)
{
    std::ofstream out(outputFilename.c_str());
	std::ifstream ifs(inputFilename);

	while (ifs.good()) {
		int a = ifs.get();
		if (!ifs.eof()) {
			out.put(a);
		}
	}
	ifs.close();
    out.close();
}


void SolveDBG::printClusteredSeq(vector<vector<string > > &seq,  vector<vector<string> > &name, const string &outPrefix)
{
	string dirName = outPrefix;
	int returnValue = mkdir(dirName.c_str(), 0755);
	if (returnValue != 0) {
		cerr << "Warning: directory " << dirName << " was not created successfully." << endl;
	}

	for (unsigned long clusterID = 0; clusterID < seq.size(); ++clusterID) {
		std::ostringstream oss;
		if (clusterID > 0)
			oss << dirName << "/cluster" << clusterID << ".fa";
		else
			oss << dirName << "/unclassified.fa";
		std::ofstream ofs(oss.str());

		for (unsigned long seqID = 0; seqID < seq[clusterID].size(); ++seqID) {
			if (seq[clusterID][seqID].empty())
				continue;

			ofs << '>' << name[clusterID][seqID] << '\n';

			for (unsigned long i = 0; i < seq[clusterID][seqID].size(); ++i) {
				ofs.put(platanus::Bin2Char(seq[clusterID][seqID][i]));
				if ((i + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
					ofs.put('\n');
			}
			if (seq[clusterID][seqID].size() % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
				ofs.put('\n');
		}
		ofs.close();
	}
}


void SolveDBG::extendConsensus(const bool bubbleRemovalFlag, const bool insertEstimationFlag, const bool tagFlag)
{
	pairedDBG.clearContigPreviousParentNodeID();

	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.setMinLink(minLink);
	pairedDBG.makeGraph(numThread);
	pairedDBG.joinUnambiguousNodePairIterative(numThread);
	if (bubbleRemovalFlag) {
		pairedDBG.makeGraph(numThread);
		pairedDBG.setOppositeBubbleContigIDOverlapped(numThread);
		pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
		pairedDBG.clearEdges();
		pairedDBG.makeScaffold();
		pairedDBG.joinUnambiguousNodePairIterative(numThread);
	}


	for (long iteration = 0; iteration < 2; ++iteration) {
		pairedDBG.setMinTagReadToLoad(1);

		for (unsigned i = 0; i < libraryMT.size(); ++i) {
			cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

			pairedDBG.setTargetLibraryIndex(i);
			unsigned tolerenceFactor = MAX_TOL_FACTOR;
			pairedDBG.setTolerence(tolerenceFactor * libraryMT[i][0].getSDInsSize());
			cerr << "TOLERENCE_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

			pairedDBG.setMinLink(minLinkToPhase);
			pairedDBG.joinUnambiguousNodePairIterative(numThread);
			pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);

			pairedDBG.setMinLink(minLink);
			pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);
		}
		if (longReadLibraryMT.size() > 0) {
			if (iteration == 0)
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE);
			else
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
			cerr << "[LONG_READ_LIBRARY]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageInsSize() << ", SD_READ_LENGTH = " << longReadLibraryMT[0].getSDInsSize() << endl;
			pairedDBG.setMinLink(minLinkToPhase);
			pairedDBG.joinUnambiguousNodePairIterative(numThread);
			pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::SCORE, clusterFlag, numThread);

			pairedDBG.setMinLink(minLink);
			pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);
		}

		if (iteration == 0)
			pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE);
		else
			pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
		for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
			pairedDBG.setTargetLibraryIndex(i);
			pairedDBG.setTolerenceFactor(MAX_TOL_FACTOR);
			pairedDBG.setMinLink(minLinkToPhase);
			pairedDBG.makeGraphAllLibraries(numThread);
			pairedDBG.joinUnambiguousNodePairIterative(numThread);
			pairedDBG.solveSimpleCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);

			pairedDBG.setMinLink(minLink);
			pairedDBG.solveSimpleCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);
		}

	}

/*
	if (tagFlag) {
		pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
		pairedDBG.makeGraph(numThread);
		pairedDBG.deleteErroneousEdgeNumTagRateIterative(true, numThread);
		pairedDBG.joinUnambiguousNodePair(numThread);
	}
*/

	if (optionBool["-no_scaffold"])
		return;


	for (long iteration = 0; iteration < 2; ++iteration) {
		pairedDBG.setMinTagReadToLoad(1);

		for (unsigned i = 0; i < libraryMT.size(); ++i) {
			pairedDBG.setTargetLibraryIndex(i);
			if (iteration == 0)
				pairedDBG.setMinLink(minLinkToPhase);
			else
				pairedDBG.setMinLink(minLink);

			cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE);

			for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
				pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
				pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
//				pairedDBG.setCutoffLength(tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
//				pairedDBG.setTolerence(pairedDBG.getCutoffLength());
				cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;
				cerr << "LENGTH_CUTOFF = " << pairedDBG.getCutoffLength() << endl;

				pairedDBG.trimSparseEnd();
//					if (iteration > 0)
//						pairedDBG.joinUnambiguousNodePairGappedIterative(numThread);
				pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);
				pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);
			}
		}
		if (longReadLibraryMT.size() > 0) {
			pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE);
			cerr << "[LONG_READ_LIBRARY]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageLength() << endl;
			pairedDBG.setTolerence(2 * this->contigMaxK);
			if (iteration == 0)
				pairedDBG.setMinLink(minLinkToPhase);
			else
				pairedDBG.setMinLink(minLink);

			for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
				pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR * longReadLibraryMT[0].getAverageInsSize());
				pairedDBG.trimSparseEnd();
				pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::SCORE, clusterFlag, numThread);
				pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);
			}
		}
		for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
			pairedDBG.setTargetLibraryIndex(i);
			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE);

			for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
				pairedDBG.setTolerenceFactor(tolerenceFactor);
				pairedDBG.setCutoffLengthFactor(MIN_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);
				pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
				pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
				if (iteration == 0)
					pairedDBG.setMinLink(minLinkToPhase);
				else {
					pairedDBG.setMinLink(minLink);
//						pairedDBG.joinUnambiguousNodePairGappedIterativeAllLibraries(numThread);
				}
				pairedDBG.trimSparseEnd();
				pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, clusterFlag, numThread);
				pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, clusterFlag, numThread);
			}
		}
	}

	pairedDBG.trimSparseEnd();


	pairedDBG.setMinOverlap(minOverlapForScaffolding);

	unsigned long long previousNumNode = pairedDBG.getNumNode();

    while(1) {
		pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);
		pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE);

		for (unsigned libraryIndex = 0; libraryIndex < libraryMT.size(); ++libraryIndex) {
			if (insertEstimationFlag)
				updateAndWiteInsertSize(libraryIndex);

			pairedDBG.setTargetLibraryIndex(libraryIndex);
			pairedDBG.setMinLink(minLink);

			cerr << "[LIBRARY " << libraryIndex + 1 << "]\nAVE_INS = " << libraryMT[libraryIndex][0].getAverageInsSize() << ", SD_INS = " << libraryMT[libraryIndex][0].getSDInsSize() << endl;
			for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
				pairedDBG.setCutoffLength(2.0 * tolerenceFactor * libraryMT[libraryIndex][0].getSDInsSize());
				pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
				cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

				pairedDBG.makeGraph(numThread);
				if (bubbleRemovalFlag) {
					pairedDBG.setOppositeBubbleContigIDGapped(numThread);
					pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
				}

				if (libraryIndex > 0) {
					pairedDBG.deleteRepeatEdge();
					pairedDBG.split();
					pairedDBG.makeGraph(numThread);
				}

				if (tagFlag)
					pairedDBG.deleteErroneousEdgeNumTagRateIterative(false, numThread);

				if (clusterFlag) {
					pairedDBG.deleteChimericEdge();
					pairedDBG.detectRepeat(0.0);
				}
				else {
					pairedDBG.detectRepeat(0.0);
					pairedDBG.deleteInterOTUEdge();
				}

				pairedDBG.makeScaffold();
			}
		}


		if (libraryMT.size() > 1) {
			long libraryIndex = libraryMT.size() - 1;
			pairedDBG.setTargetLibraryIndex(libraryIndex);
			pairedDBG.setMinLink(minLink);
			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE);

			for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
				pairedDBG.setTolerenceFactor(tolerenceFactor);
				pairedDBG.setCutoffLengthFactor(MAX_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);
				pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[libraryIndex][0].getSDInsSize(), 0.1 * libraryMT[libraryIndex][0].getAverageInsSize()));
				pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);

				pairedDBG.makeGraphAllLibraries(numThread);
				pairedDBG.deleteLongEdge(libraryMT[libraryIndex][0].getAverageInsSize());
				pairedDBG.deleteRepeatEdge();

				if (clusterFlag) {
					pairedDBG.deleteChimericEdge();
					pairedDBG.detectRepeat(0.0);
				}
				else {
					pairedDBG.detectRepeat(0.0);
					pairedDBG.deleteInterOTUEdge();
				}

				pairedDBG.makeScaffold();
			}
		}


        if (previousNumNode <= pairedDBG.getNumNode())
			break;

        previousNumNode = pairedDBG.getNumNode();
	}


	if (longReadLibraryMT.size() > 0) {
		pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE);

		for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
			pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR *  longReadLibraryMT[0].getAverageInsSize());
			pairedDBG.setTolerence(std::min(MIN_LONG_READ_LENGTH_CUTOFF_FACTOR * longReadLibraryMT[0].getAverageInsSize(), 0.5 * pairedDBG.getCutoffLength()));
			pairedDBG.setMinLink(minLink);
			pairedDBG.makeGraph(numThread);
			if (bubbleRemovalFlag) {
				pairedDBG.setOppositeBubbleContigIDGapped(numThread);
				pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
			}
			if (tagFlag)
				pairedDBG.deleteErroneousEdgeNumTagRateIterative(false, numThread);

			pairedDBG.deleteErroneousEdgeNumLinkRateIterative(numThread);
			pairedDBG.deleteRepeatEdge();

			if (clusterFlag) {
				pairedDBG.deleteChimericEdge();
				pairedDBG.detectRepeat(0.0);
			}
			else {
				pairedDBG.detectRepeat(0.0);
				pairedDBG.deleteInterOTUEdge();
			}

			pairedDBG.makeScaffold();

			pairedDBG.trimSparseEnd();
			pairedDBG.trimRepeatEnd();
		}
	}


	if (optionBool["-no_barcode_scaffold"] || tagLibraryMT.size() == 0 || !(tagFlag || insertEstimationFlag)) {
		pairedDBG.setMinOverlap(this->contigMaxK - 1);
		return;
	}


	pairedDBG.setMinLink(minLinkTagScaffold);
	pairedDBG.setMinTagReadToLoad(MIN_NUM_FOR_TAG_SCAFFOLD_ISLAND);

	if (insertEstimationFlag) {
		pairedDBG.setMode(PairedDBG::TAG_SCAFFOLD_MODE);
		pairedDBG.estimateTaggedMoleculeLength(numThread);

		if (bubbleRemovalFlag) {
			pairedDBG.setMinOverlap(this->contigMaxK - 1);
			return;
		}
	}


	pairedDBG.setMode(PairedDBG::TAG_SCAFFOLD_MODE | PairedDBG::LENGTH_CUTOFF_MODE);

	for (long cutoffFactor = MIN_TAG_SCAFFOLD_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_TAG_SCAFFOLD_LENGTH_CUTOFF_FACTOR + 1; ++cutoffFactor) {
		if (cutoffFactor <= MAX_TAG_SCAFFOLD_LENGTH_CUTOFF_FACTOR)
			pairedDBG.setCutoffLength(cutoffFactor * TAG_SCAFFOLD_LENGTH_CUTOFF_UNIT_FACTOR * tagLibraryMT[0].getAverageInsSize());
		else
			pairedDBG.setCutoffLength(TAG_SCAFFOLD_LENGTH_LOWER_CUTOFF_UNIT_FACTOR * tagLibraryMT[0].getAverageInsSize());

		pairedDBG.setTolerence(tagLibraryMT[0].getAverageInsSize());

		for (long islandCutoffFactor = MIN_TAG_SCAFFOLD_ISLAND_LENGTH_CUTOFF_FACTOR; islandCutoffFactor <= MAX_TAG_SCAFFOLD_ISLAND_LENGTH_CUTOFF_FACTOR; ++islandCutoffFactor) {
			pairedDBG.setMinTagIslandLengthToScaffold(islandCutoffFactor * TAG_SCAFFOLD_ISLAND_LENGTH_CUTOFF_UNIT_FACTOR *  tagLibraryMT[0].getAverageInsSize());

			pairedDBG.makeGraph(numThread);
			pairedDBG.deleteConflictingEdgeToSameNode(numThread);

			pairedDBG.setTolerence(1);
			pairedDBG.deleteErroneousEdgeNumLinkRateIterative(numThread);

			pairedDBG.setTolerence(tagLibraryMT[0].getAverageInsSize());

			if (clusterFlag) {
				pairedDBG.deleteChimericEdge();
				pairedDBG.detectRepeat(0.0);
			}
			else {
				pairedDBG.detectRepeat(0.0);
				pairedDBG.deleteInterOTUEdge();
			}

			pairedDBG.makeScaffold();
		}
	}


	pairedDBG.setMinOverlap(this->contigMaxK - 1);

}
