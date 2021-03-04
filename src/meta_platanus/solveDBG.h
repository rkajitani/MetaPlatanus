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

#ifndef SOLVE_DBG_H
#define SOLVE_DBG_H

#include "baseCommand.h"
#include "mapper.h"
#include "scaffold.h"
#include "scaffoldGraph.h"
#include "pairedDBG.h"
#include <cmath>


//class SolveDBG : public BaseCommand
class SolveDBG : public Scaffold
{
private:

	static const double MIN_LONG_READ_LENGTH_CUTOFF_FACTOR;
	static const double MAX_LONG_READ_LENGTH_CUTOFF_FACTOR;
    static const double LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR;
	static const long LONG_READ_MIN_ALIGNMENT_LENGTH;
	static const double LONG_READ_MIN_ALIGNMENT_COVERAGE;
	static const double LONG_READ_MIN_ALIGNMENT_IDENTITY;

	static const long MIN_TAG_SCAFFOLD_LENGTH_CUTOFF_FACTOR;
	static const long MAX_TAG_SCAFFOLD_LENGTH_CUTOFF_FACTOR;
	static const double TAG_SCAFFOLD_LENGTH_CUTOFF_UNIT_FACTOR;
	static const double TAG_SCAFFOLD_LENGTH_LOWER_CUTOFF_UNIT_FACTOR;

	static const long MIN_TAG_SCAFFOLD_ISLAND_LENGTH_CUTOFF_FACTOR;
	static const long MAX_TAG_SCAFFOLD_ISLAND_LENGTH_CUTOFF_FACTOR;
	static const long MIN_NUM_FOR_TAG_SCAFFOLD_ISLAND;
	static const double TAG_SCAFFOLD_ISLAND_LENGTH_CUTOFF_UNIT_FACTOR;

    PairedDBG pairedDBG;
    PairedDBG unphasedGraph;
    PairedDBG phasedGraph;
	std::vector<int> multiSeedLengthForShortRead;
	std::vector<int> multiSeedLengthForLongRead;
	long minLinkToPhase;
	long minLinkTagScaffold;
	long minOverlapForScaffolding;
	bool clusterFlag;
    std::vector<std::string> pacBioLongReadFilename;
    std::vector<std::string> nanoporeLongReadFilename;
    std::vector<std::string> guideContigLongReadFilename;

    virtual bool parseArgs(int argc, char **argv)
    {
        if (BaseCommand::parseArgs(argc, argv)) {
            int optInd = 2;
            while (optInd < argc) {
                if (strstr(argv[optInd], "-n") == argv[optInd]) {
                    int pairNum = divideArgvInt(argv[optInd]);
                    if (pairNum > 0) {
						int minInsertSize = atoi(argv[optInd + 1]);
						optionMinIns[pairNum] = minInsertSize;
						optInd += 2;
                    }
					else {
						++optInd;
					}
                } else if (strstr(argv[optInd], "-a") == argv[optInd]) {
                    int pairNum = divideArgvInt(argv[optInd]);
                    if (pairNum > 0) {
						int minInsertSize = atoi(argv[optInd + 1]);
						optionAveIns[pairNum] = minInsertSize;
						optInd += 2;
					}
					else {
						++optInd;
					}
                } else if (strstr(argv[optInd], "-d") == argv[optInd]) {
                    int pairNum = divideArgvInt(argv[optInd]);
                    if (pairNum > 0) {
						int minInsertSize = atoi(argv[optInd + 1]);
						optionSDIns[pairNum] = minInsertSize;
						optInd += 2;
                    }
					else {
						++optInd;
					}
                } else if (strstr(argv[optInd], "-z") == argv[optInd]) {
                    int pairNum = divideArgvInt(argv[optInd]);
                    if (pairNum > 0) {
						double insCutoffRate = atof(argv[optInd + 1]);
						optionInsCutoffRate[pairNum] = insCutoffRate;
						optInd += 2;
					}
					else {
						++optInd;
					}
                } else {
                    ++optInd;
                }
            }
        } else return false;
        return true;
    }

    virtual bool checkFileEnough(void)
    {
		if (optionMultiArgs["-c"].size() == 0) {
			std::cerr << "Error: not specified contig file!!" << std::endl;
			return false;
		}
		if (!(optionBool["-trim_cluster"]) && optionSingleArgs["-k"] == "" && optionSingleArgs["-g"] == "") {
			std::cerr << "Error: both -k(k-mer-occurrence-file) and -g(clustering-result-file) is not specified!!" << std::endl;
			return false;
		}
		return true;
    }

    virtual int checkOtherOption(char *argv) const
    {
        if (strstr(argv, "-n") == argv)
            return 2;
        if (strstr(argv, "-a") == argv)
            return 2;
        if (strstr(argv, "-d") == argv)
            return 2;
        if (strstr(argv, "-z") == argv)
            return 2;
        else
            return 0;
    }

public:
    SolveDBG();
    SolveDBG(const SolveDBG &) = delete;
    SolveDBG &operator=(const SolveDBG &) = delete;
    ~SolveDBG() = default;

    virtual void usage(void) const;
    virtual void exec(void);
    virtual void initializeParameters(void);
    virtual void readLibrary(std::unique_ptr<HeteroMapper> &mapper, platanus::Contig &contig, const int numThread);
    virtual void mapLibraryAndInitGraph(const int numThread);
    virtual void outputAndAfterTreatment(void);
	void outputGraph(const char *suffix);
	void updateAndWiteInsertSize(const long libraryIndex);
	void extendConsensus(const bool bubbleRemovalFlag, const bool insertEstimationFlag, const bool tagFlag);
	void execMinialign(const std::string targetFilename, const std::string &readFilename, const std::string &outFilename, const long numThread, const std::string minialignExecutable);
	void execMinimap(const std::vector<std::string> &contigFilenames, const std::vector<std::string> &bubbleFilenames, const std::vector<std::string> &readFilenames, const std::string &outFilename, const long minAlignmentLength, const long numThread, const std::string minimapExecutable);
	void execMinimap2(const std::vector<std::string> &contigFilenames, const std::vector<std::string> &bubbleFilenames, const std::vector<std::string> &readFilenames, const std::string &outFilename, const long numThread, const std::string minimap2Executable, const std::string minimap2Option);
    void outputForException(void);
	void printCatContigAsScaffold(const std::string &outputFilename, const std::string &inputFilename);
	void printClusteredSeq(std::vector<std::vector<std::string > > &seq,  std::vector<std::vector<std::string> > &name, const std::string &outPrefix);
};




#endif
