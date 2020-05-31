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

#include "assemble.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <climits>
#include <cfloat>

using std::string;
using std::vector;
using std::unordered_map;
using std::cerr;
using std::endl;
using std::ifstream;
using platanus::SEQ;

//////////////////////////////////////////////////////////////////////////////////////
// static constant
//////////////////////////////////////////////////////////////////////////////////////
const unsigned Assemble::MIN_KMER_LENGTH = 25;
const unsigned Assemble::SMOOTHING_WINDOW = 1;
const double Assemble::MAX_COVERAGE_CUT_DIFF_RATE = 0.25;


//////////////////////////////////////////////////////////////////////////////////////
// default constructor
//////////////////////////////////////////////////////////////////////////////////////
Assemble::Assemble()
: numExtension(0), averageLength(0), averageCoverage(0), numThread(0), lengthCutoff(0), numKmer(0), memory(0), kmer31Graph(), kmer63Graph(), kmer95Graph(), kmer127Graph(), kmer159Graph(), kmerNGraph(), kmer31Counter(), kmer63Counter(), kmer95Counter(), kmer127Counter(), kmer159Counter(), kmerNCounter(), doubleHashSize()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-e"] = "";
    optionSingleArgs["-k"] = "0.25";
    optionSingleArgs["-K"] = "0.9";
    optionSingleArgs["-s"] = "2";
    optionSingleArgs["-n"] = "0";

    optionSingleArgs["-c"] = "4";
    optionSingleArgs["-C"] = "8";
    optionSingleArgs["-l"] = "2.0";
    optionSingleArgs["-u"] = "0.02";

    optionSingleArgs["-d"] = "0.5";
    optionSingleArgs["-t"] = "1";
    optionSingleArgs["-m"] = "16";
    optionMultiArgs["-f"] = vector<string>();
    optionSingleArgs["-tmp"] = ".";
    optionBool["-repeat"] = false;
}



//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void Assemble::usage(void) const
{
    std::cerr << "\nUsage meta_platanus2 assemble [Options]\n"
              << "Options:\n"
              << "    -o STR               : prefix of output files (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -f FILE1 [FILE2 ...] : reads file (fasta or fastq, number <= "<< platanus::ConstParam::MAX_FILE_NUM << ")\n"
              << "    -k INT               : initial k-mer size (default " << optionSingleArgs.at("-k") << ")\n"
              << "    -K FLOAT             : maximum-k-mer factor (maximum-k = FLOAT*read-length, default  " << optionSingleArgs.at("-K") << ")\n"
              << "    -s INT               : step size of k-mer extension (>= 1, default " << optionSingleArgs.at("-s") << ")\n"
              << "    -n INT               : initial k-mer coverage cutoff (default " << optionSingleArgs.at("-n") << ", 0 means auto)\n"
              << "    -c INT               : minimun k-mer coverage (default " << optionSingleArgs.at("-c") << ")\n"
              << "    -u FLOAT             : maximum difference for bubble crush (identity, default " << optionSingleArgs.at("-u") << ")\n"
              << "    -d FLOAT             : maximum difference for branch cutting (coverage ratio, default " << optionSingleArgs.at("-d") << ")\n"
              << "    -e FLOAT             : k-mer coverage depth (k = initial k-mer size specified by -k) of homozygous region (default auto)\n"
              << "    -t INT               : number of threads (<= " << platanus::ConstParam::MAX_THREAD << ", default " << optionSingleArgs.at("-t") << ")\n"
              << "    -m INT               : memory limit for making kmer distribution (GB, >=1, default " << optionSingleArgs.at("-m") << ")\n"
              << "    -tmp DIR             : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n"
              << "\n\n"
              << "Outputs:\n"
              << "    PREFIX_contig.fa\n"
              << "    PREFIX_junctionKmer.fa\n"
              << "    PREFIX_kmerFrq.tsv\n"
              << "\n"
              << "Note, PREFIX is specified by -o\n"
              << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// initilaize parametors
//////////////////////////////////////////////////////////////////////////////////////
void Assemble::initializeParameters(void)
{
    this->numThread = atoi(optionSingleArgs["-t"].c_str());
    this->memory = atoi(optionSingleArgs["-m"].c_str()) * static_cast<unsigned long long>(1000000000);
    this->numExtension = atoi(optionSingleArgs["-s"].c_str());
    this->upperCoverageCutoff = atoi(optionSingleArgs["-C"].c_str());
    this->lowerCoverageCutoff = atoi(optionSingleArgs["-c"].c_str());
    this->minKmerLengthRatio = atof(optionSingleArgs["-k"].c_str());
    this->maxKmerLengthRatio = atof(optionSingleArgs["-K"].c_str());
    this->tableBinaryName = optionSingleArgs["-o"] + "_kmer_occ.bin";
    const double bubble = atof(optionSingleArgs["-u"].c_str());
    const double branch = atof(optionSingleArgs["-d"].c_str());
    initializeGraph(bubble, branch);
    platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
}

//////////////////////////////////////////////////////////////////////////////////////
// exec assemble
//////////////////////////////////////////////////////////////////////////////////////
void Assemble::exec(void)
{
    initializeParameters();
    vector<unsigned> k;
    FILE *readFP[numThread];
    FILE *allReadFP[numThread];

    string outputFilename;


    k.push_back(atoi(optionSingleArgs["-k"].c_str()));

    omp_set_num_threads(numThread);
    // read file
    for (unsigned long long i = 0; i < numThread; ++i)
        readFP[i] = platanus::makeTemporaryFile();
    for (unsigned i = 0; i < optionMultiArgs["-f"].size(); ++i) {
        readInputFile(optionMultiArgs["-f"][i], readFP, numThread);
    }
    for (unsigned long long i = 0; i < numThread; ++i)
        allReadFP[i] = readFP[i];


    this->extendKmer(readFP, k);

    // execute initial assemble
    if (k[0] <= 32)
        initialKmerAssemble(k, readFP, allReadFP, numThread, kmer31Graph, kmer31Counter);
    else if (k[0] <= 64)
        initialKmerAssemble(k, readFP, allReadFP, numThread, kmer63Graph, kmer63Counter);
    else if (k[0] <= 96)
        initialKmerAssemble(k, readFP, allReadFP, numThread, kmer95Graph, kmer95Counter);
    else if (k[0] <= 128)
        initialKmerAssemble(k, readFP, allReadFP, numThread, kmer127Graph, kmer127Counter);
    else if (k[0] <= 160)
        initialKmerAssemble(k, readFP, allReadFP, numThread, kmer159Graph, kmer159Counter);
    else
        initialKmerAssemble(k, readFP, allReadFP, numThread, kmerNGraph, kmerNCounter);

    // kmer extension and assemble iterative
    for (unsigned i = 1; i < numKmer; ++i) {
        if (k[i - 1] <= 32)
            saveAndRedoAssemble(k[i], k[i - 1], readFP, allReadFP, numThread, kmer31Graph, kmer31Counter, i);
        else if (k[i - 1] <= 64)
            saveAndRedoAssemble(k[i], k[i - 1], readFP, allReadFP, numThread, kmer63Graph, kmer63Counter, i);
        else if (k[i - 1] <= 96)
            saveAndRedoAssemble(k[i], k[i - 1], readFP, allReadFP, numThread, kmer95Graph, kmer95Counter, i);
        else if (k[i - 1] <= 128)
            saveAndRedoAssemble(k[i], k[i - 1], readFP, allReadFP, numThread, kmer127Graph, kmer127Counter, i);
        else if (k[i - 1] <= 160)
            saveAndRedoAssemble(k[i], k[i - 1], readFP, allReadFP, numThread, kmer159Graph, kmer159Counter, i);
        else
            saveAndRedoAssemble(k[i], k[i - 1], readFP, allReadFP, numThread, kmerNGraph, kmerNCounter, i);
    }


    if (k[numKmer - 1] <= 32) {
        outputAndAfterTreatment(k[numKmer - 1], k[0], allReadFP, kmer31Graph, numThread);
    } else if (k[numKmer - 1] <= 64) {
        outputAndAfterTreatment(k[numKmer - 1], k[0], allReadFP, kmer63Graph, numThread);
    } else if (k[numKmer - 1] <= 96) {
        outputAndAfterTreatment(k[numKmer - 1], k[0], allReadFP, kmer95Graph, numThread);
    } else if (k[numKmer - 1] <= 128) {
        outputAndAfterTreatment(k[numKmer - 1], k[0], allReadFP, kmer127Graph, numThread);
    } else if (k[numKmer - 1] <= 160) {
        outputAndAfterTreatment(k[numKmer - 1], k[0], allReadFP, kmer159Graph, numThread);
    } else {
        outputAndAfterTreatment(k[numKmer - 1], k[0], allReadFP, kmerNGraph, numThread);
    }

    cerr << "assemble completed!" << endl;
}


template <typename KMER>
void Assemble::mergeContig(platanus::Contig &contig, Counter<KMER> &counter, BruijnGraph<KMER> &graph, const unsigned long long kmerLength, FILE **mergedContigFP)
{
    graph.setBubbleAndBranch(atof(optionSingleArgs["-u"].c_str()), atof(optionSingleArgs["-d"].c_str()));
    counter.setKmerLength(kmerLength);
    graph.setKmerLength(kmerLength);
	double averageCoverage = contig.calculateAverageCoverageExcludingOutlier(contig.getSeqLengthMedian());
    unsigned long long doubleHashSize = counter.makeKmerReadDistributionFromContig(contig, kmerLength, 1, memory);

    FILE *sortedKeyFP = counter.sortedKeyFromKmerFile(0);
    counter.loadKmer(0, doubleHashSize);
    graph.makeInitialBruijnGraph(counter, sortedKeyFP);
    fclose(sortedKeyFP);
    counter.deleteAllTable();

    graph.cutBranchIterative(numThread);
    graph.crushBubbleIterative();

    graph.saveContigSimple(*mergedContigFP, 1.0);
}




//////////////////////////////////////////////////////////////////////////////////////
// initial assemble
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Assemble::initialKmerAssemble(vector<unsigned> &k, FILE **readFP, FILE **allReadFP, const unsigned long long numThread, BruijnGraph<KMER> &graph, Counter<KMER> &counter)
{
    cerr << "K = " << k[0] << ", saving kmers from reads..." << endl;
    FILE *sortedKeyFP;
    std::ostringstream oss;

    counter.setKmerLength(k[0]);
    graph.setKmerLength(k[0]);

    // make kmer distribution
    this->doubleHashSize = counter.makeKmerReadDistributionMT(k[0], readFP, memory, numThread);

    // calculate various value and decide kmer extension
    averageCoverage = counter.calcOccurrenceDistributionAverage(1, counter.getMaxOccurrence());
    averageCoverage = averageCoverage * averageLength / (averageLength - k[0] + 1.0);

    // output kmer frequency distribution
    oss << optionSingleArgs["-o"] << '_' << k[0] << "merFrq.tsv";
    string kmerFrqFilename = oss.str();
    counter.outputOccurrenceDistribution(kmerFrqFilename);

    // load kmer from temporary file
    cerr << "COVERAGE_CUTOFF = " << lowerCoverageCutoff << endl;
    sortedKeyFP = counter.sortedKeyFromKmerFile(this->lowerCoverageCutoff);
    this->doubleHashSize = counter.loadKmer(this->lowerCoverageCutoff, this->doubleHashSize);
    counter.outputOccurrenceTableBinary(this->tableBinaryName);

    // make bruijn graph
    graph.makeInitialBruijnGraph(counter, sortedKeyFP);
    fclose(sortedKeyFP);
    counter.deleteAllTable();

    graph.cutBranchIterative(numThread);
    graph.deleteErroneousStraightNodeIterative(static_cast<long>(atof(optionSingleArgs["-l"].c_str()) * averageLength + 0.5) - k[0] + 1, this->upperCoverageCutoff, numThread);
    graph.cutBranchIterative(numThread);
    graph.crushBubbleIterative();
    //graph.divideStraightNode(allReadFP, this->lowerCoverageCutoff, this->doubleHashSize, numThread);

}


//////////////////////////////////////////////////////////////////////////////////////
// save contig and kmer extension assemble
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Assemble::saveAndRedoAssemble(const unsigned k, const unsigned previousk, FILE **readFP, FILE **allReadFP, const unsigned long long numThread, BruijnGraph<KMER> &lowKmerGraph, Counter<KMER> &lowKmerCounter, const unsigned long long position)
{

//    lowKmerGraph.divideStraightNode(allReadFP, coverageCutoff[position - 1], this->doubleHashSize, numThread);
    saveGraph(k, readFP, numThread, lowKmerGraph, lowKmerCounter);
    if (k <= 32) {
        updateDoubleHashSize<Kmer31::keyType>(lowKmerGraph, k, previousk);
        DoubleHash<Kmer31::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        lowKmerGraph.template saveContig<Kmer31>(k, (averageLength - k + 1.0) / (averageLength - previousk + 1.0), this->upperCoverageCutoff, tmpOccurrenceTable);
        lowKmerGraph.deleteAllTable();
        kmer31Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
        redoAssemble(k, previousk, readFP, allReadFP, numThread, kmer31Graph, kmer31Counter, position);
    } else if (k <= 64) {
        updateDoubleHashSize<KmerN<Binstr63>::keyType>(lowKmerGraph, k, previousk);
        DoubleHash<KmerN<Binstr63>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        lowKmerGraph.template saveContig<KmerN<Binstr63> >(k, (averageLength - k + 1.0) / (averageLength - previousk + 1.0), this->upperCoverageCutoff, tmpOccurrenceTable);
        lowKmerGraph.deleteAllTable();
        kmer63Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
        redoAssemble(k, previousk, readFP, allReadFP, numThread, kmer63Graph, kmer63Counter, position);
    } else if (k <= 96) {
        updateDoubleHashSize<KmerN<Binstr95>::keyType>(lowKmerGraph, k, previousk);
        DoubleHash<KmerN<Binstr95>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        lowKmerGraph.template saveContig<KmerN<Binstr95> >(k, (averageLength - k + 1.0) / (averageLength - previousk + 1.0), this->upperCoverageCutoff, tmpOccurrenceTable);
        lowKmerGraph.deleteAllTable();
        kmer95Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
        redoAssemble(k, previousk, readFP, allReadFP, numThread, kmer95Graph, kmer95Counter, position);
    } else if (k <= 128) {
        updateDoubleHashSize<KmerN<Binstr127>::keyType>(lowKmerGraph, k, previousk);
        DoubleHash<KmerN<Binstr127>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        lowKmerGraph.template saveContig<KmerN<Binstr127> >(k, (averageLength - k + 1.0) / (averageLength - previousk + 1.0), this->upperCoverageCutoff, tmpOccurrenceTable);
        lowKmerGraph.deleteAllTable();
        kmer127Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
        redoAssemble(k, previousk, readFP, allReadFP, numThread, kmer127Graph, kmer127Counter, position);
    } else if (k <= 160) {
        updateDoubleHashSize<KmerN<Binstr159>::keyType>(lowKmerGraph, k, previousk);
        DoubleHash<KmerN<Binstr159>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        lowKmerGraph.template saveContig<KmerN<Binstr159> >(k, (averageLength - k + 1.0) / (averageLength - previousk + 1.0), this->upperCoverageCutoff, tmpOccurrenceTable);
        lowKmerGraph.deleteAllTable();
        kmer159Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
        redoAssemble(k, previousk, readFP, allReadFP, numThread, kmer159Graph, kmer159Counter, position);
    } else {
        updateDoubleHashSize<KmerN<binstr_t>::keyType>(lowKmerGraph, k, previousk);
        DoubleHash<KmerN<binstr_t>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        lowKmerGraph.template saveContig<KmerN<binstr_t> >(k, (averageLength - k + 1.0) / (averageLength - previousk + 1.0), this->upperCoverageCutoff, tmpOccurrenceTable);
        lowKmerGraph.deleteAllTable();
        kmerNCounter.swapOccurrenceTable(k, tmpOccurrenceTable);
        redoAssemble(k, previousk, readFP, allReadFP, numThread, kmerNGraph, kmerNCounter, position);
    }
}



//////////////////////////////////////////////////////////////////////////////////////
// save contig made brujin graph and prepare next kmer step
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Assemble::saveGraph(const unsigned k, FILE **readFP, const unsigned long long numThread, BruijnGraph<KMER> &graph, Counter<KMER> &counter)
{
    unsigned long long i;
    counter.setOccurrenceTableSize(this->doubleHashSize);
    graph.saveEdgeKmer(counter, k);
    cerr << "extracting reads (containing kmer used in contig assemble)..." << endl;
    # pragma omp parallel for schedule(static, 1)
    for (i = 0; i < numThread; ++i)
        counter.pickupReadMatchedEdgeKmer(&readFP[i]);
    counter.deleteAllTable();
}



//////////////////////////////////////////////////////////////////////////////////////
// load new kmer and make de brujin graph
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Assemble::redoAssemble(const unsigned k, const unsigned previousk, FILE **readFP, FILE **allReadFP, unsigned long long numThread, BruijnGraph<KMER> &graph, Counter<KMER> &counter, const unsigned long long position)
{
    FILE *extractFP[numThread];
    counter.makeKmerReadDistributionConsideringPreviousGraph(k, readFP, memory, numThread);
    cerr << "COVERAGE_CUTOFF = " << upperCoverageCutoff << endl;


    // load kmer from temporary file
    FILE *sortedKeyFP = counter.sortedKeyFromKmerFile(this->upperCoverageCutoff);
    this->doubleHashSize = counter.loadKmer(this->upperCoverageCutoff, this->doubleHashSize);

    // make bruijn graph
    graph.setKmerLength(k);
    graph.makeInitialBruijnGraph(counter, sortedKeyFP);
    fclose(sortedKeyFP);
    counter.deleteAllTable();

    graph.cutBranchIterative(numThread);
    graph.crushBubbleIterative();
    // graph.divideStraightNode(allReadFP, this->upperCoverageCutoff, this->doubleHashSize, numThread);
    # pragma omp parallel for schedule(static, 1)
    for (unsigned long long i = 0; i < numThread; ++i) {
        graph.extractRead(readFP[i], &extractFP[i]);
    }

    DoubleHash<typename KMER::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
    graph.template saveContig<KMER>(k, 1.0, this->lowerCoverageCutoff / 2, tmpOccurrenceTable);
    graph.deleteAllTable();
    counter.swapOccurrenceTable(k, tmpOccurrenceTable);


    counter.makeKmerReadDistributionConsideringPreviousGraph(k, extractFP, memory, numThread);
    for (unsigned long long i = 0; i < numThread; ++i) {
        fclose(extractFP[i]);
    }

    cerr << "COVERAGE_CUTOFF = " << lowerCoverageCutoff / 2 << endl;

    sortedKeyFP = counter.sortedKeyFromKmerFile(this->lowerCoverageCutoff / 2);
    this->doubleHashSize = counter.loadKmer(this->lowerCoverageCutoff / 2, this->doubleHashSize);
    graph.makeInitialBruijnGraph(counter, sortedKeyFP);
    fclose(sortedKeyFP);
    counter.deleteAllTable();
 
    graph.cutBranchIterative(numThread);
    graph.deleteErroneousStraightNodeIterative(static_cast<long>(atof(optionSingleArgs["-l"].c_str()) * averageLength + 0.5) - k + 1, this->upperCoverageCutoff, numThread);
    graph.cutBranchIterative(numThread);
    graph.crushBubbleIterative();
    // graph.divideStraightNode(allReadFP, this->lowerCoverageCutoff, this->doubleHashSize, numThread);

}


template <typename KMER>
void Assemble::outputAndAfterTreatment(const unsigned k, const unsigned initialK, FILE **readFP, BruijnGraph<KMER> &graph, const unsigned long long numThread)
{
    string outputFilename;

    // graph.cutBranchIterative(numThread);

    graph.crushBubbleIterative();
    graph.deleteErroneousStraightNodeIterative(k, this->upperCoverageCutoff * 2, numThread);
    graph.crushBubbleIterative();
    for (unsigned long long i = 0; i < numThread; ++i)
        fclose(readFP[i]);

    FILE *contigFP = platanus::makeTemporaryFile();
    graph.saveContigSimple(contigFP);

    FILE *junctionFP = platanus::makeTemporaryFile();
    graph.saveJunction(junctionFP);

    graph.deleteAllTable();

    if (initialK <= 32) {
        kmer31Counter.readOccurrenceTableBinary(this->tableBinaryName);
        kmer31Counter.reCalculateCoverage(&contigFP);
        kmer31Counter.reCalculateCoverage(&junctionFP);
    } else if (initialK <= 64) {
        kmer63Counter.readOccurrenceTableBinary(this->tableBinaryName);
        kmer63Counter.reCalculateCoverage(&contigFP);
        kmer63Counter.reCalculateCoverage(&junctionFP);
    } else if (initialK <= 96) {
        kmer95Counter.readOccurrenceTableBinary(this->tableBinaryName);
        kmer95Counter.reCalculateCoverage(&contigFP);
        kmer95Counter.reCalculateCoverage(&junctionFP);
    } else if (initialK <= 128) {
        kmer127Counter.readOccurrenceTableBinary(this->tableBinaryName);
        kmer127Counter.reCalculateCoverage(&contigFP);
        kmer127Counter.reCalculateCoverage(&junctionFP);
    } else if (initialK <= 160) {
        kmer159Counter.readOccurrenceTableBinary(this->tableBinaryName);
        kmer159Counter.reCalculateCoverage(&contigFP);
        kmer159Counter.reCalculateCoverage(&junctionFP);
    } else {
        kmerNCounter.readOccurrenceTableBinary(this->tableBinaryName);
        kmerNCounter.reCalculateCoverage(&contigFP);
        kmerNCounter.reCalculateCoverage(&junctionFP);
    }

    outputFilename = optionSingleArgs["-o"];
    outputFilename += "_contig.fa";
    platanus::printContig(outputFilename, contigFP, averageLength / (averageLength - initialK + 1.0) , averageLength, k, "seq");
	fclose(contigFP);

    outputFilename = optionSingleArgs["-o"];
    outputFilename += "_junctionKmer.fa";
    platanus::printContig(outputFilename, junctionFP, averageLength / (averageLength - initialK + 1.0), averageLength, k, "junction");
	fclose(junctionFP);
}


//////////////////////////////////////////////////////////////////////////////////////
// decide how extend kmer
//////////////////////////////////////////////////////////////////////////////////////
void Assemble::extendKmer(FILE **readFP, vector<unsigned> &kmer)
{
    this->averageLength = this->getAverageLengthFromReadFile(readFP);

    cerr << "AVE_READ_LEN=" << this->averageLength << endl;
//    kmer[0] = static_cast<long>(this->averageLength * this->minKmerLengthRatio + 0.5);
    kmer[0] = 25;

	if (kmer[0] < MIN_KMER_LENGTH && this->averageLength * this->maxKmerLengthRatio > MIN_KMER_LENGTH) {
		kmer[0] = MIN_KMER_LENGTH;
	}

    const double diffMaxMinKmerLength = this->averageLength * this->maxKmerLengthRatio - kmer[0];
    long lengthStep = static_cast<long>(diffMaxMinKmerLength / static_cast<double>(this->numExtension) + 0.5);
    numKmer = static_cast<long>(diffMaxMinKmerLength / lengthStep) + 1;
    kmer.resize(numKmer);
    for (unsigned long long i = 1; i < numKmer; ++i) {
        kmer[i] = kmer[i - 1] + lengthStep;
    }
    if (static_cast<long>(diffMaxMinKmerLength) % lengthStep != 0) {
        kmer.emplace_back(diffMaxMinKmerLength + kmer[0]);
		++numKmer;
    }
}


double Assemble::getAverageLengthFromReadFile(FILE **readFP)
{
    unsigned long long numRead = 0;
    unsigned long long sumBase = 0;

    # pragma omp parallel for schedule(static, 1) reduction(+: numRead, sumBase)
    for (unsigned long long i = 0; i < this->numThread; ++i) {
        rewind(readFP[i]);
        platanus::SEQ seq;
        while (seq.readTemporaryFile(readFP[i])) {
            ++numRead;
            sumBase += seq.length;
        }
    }
    return static_cast<double>(sumBase) / numRead;
}


//////////////////////////////////////////////////////////////////////////////////////
// calculate something value (I don't know detail...)
//////////////////////////////////////////////////////////////////////////////////////
double Assemble::calcLogProbabilityJoin(const unsigned long long coverageCutoff, const double averageCoverage, const double averageLength, const unsigned long long largeKmer, const unsigned long long smallKmer) const
{
    double p;
    double s = 0;
    double tmpAverageCoverage = averageCoverage * (averageLength - largeKmer + 1.0) / averageLength;

    for (unsigned long long i = 0; i < coverageCutoff; ++i) {
        p = 0;
        for (unsigned long long j = 1; j <= i; ++j)
            p += log(tmpAverageCoverage) - log(static_cast<double>(j));
        s += exp(p);
    }
    s = exp(-tmpAverageCoverage + log(s));

    return ((largeKmer - smallKmer) + 1.0) * (-s);
}




//////////////////////////////////////////////////////////////////////////////////////
// decrease coverage cutoff value
//////////////////////////////////////////////////////////////////////////////////////
unsigned long long Assemble::decreaseCoverageCutoff(const unsigned long long coverageCutoff, const double averageCoverage, const double averageLength, const double minLogPJoin, const unsigned long long largeKmer, const unsigned long long smallKmer) const
{

    unsigned long long i;

    if (coverageCutoff <= 1)
        return 1;

    for (i = coverageCutoff; i > 1; --i) {
        if (calcLogProbabilityJoin(i, averageCoverage, averageLength, largeKmer, smallKmer) > minLogPJoin)
            break;
    }
    return i;
}


template <typename KEY, typename KMER>
void Assemble::updateDoubleHashSize(const BruijnGraph<KMER> &graph, const unsigned k, const unsigned previousk)
{
    unsigned long long estimate = graph.estimateNumKmerOnStraight();
    unsigned long long size = log(estimate / platanus::ConstParam::DOUBLE_HASH_MAX_LOAD_FACTOR) / log(2);
    size = std::pow(2, size + 1);

    long base = sizeof(std::pair<KEY, unsigned short>);
    unsigned long long tmpMemory = this->memory / base;
    base = log(tmpMemory) / log(2);
    this->doubleHashSize = std::pow(2, base);
    if (this->doubleHashSize > this->memory) {
        this->doubleHashSize >>= 1;
    }

    if (size > this->doubleHashSize) {
        platanus::MemoryAlert();
    }
    this->doubleHashSize = std::max(size, this->doubleHashSize);

}


// below this line there are only read function

//////////////////////////////////////////////////////////////////////////////////////
// read functions
//////////////////////////////////////////////////////////////////////////////////////
void Assemble::readInputFile(const string &inputFilename, FILE **outputMT, const unsigned long long numThread) const
{
    platanus::FILETYPE fileType = checkFileFormat(inputFilename);
    ifstream ifs(inputFilename);
    ifs.close();
    if (fileType == platanus::FILETYPE::FASTA)
        readFasta(inputFilename, outputMT, numThread);
    else if (fileType > platanus::FILETYPE::FASTA)
        readFastq(inputFilename, outputMT, numThread);
    else
        throw platanus::ReadError("Read file exception!!\nRead file is not FASTA/FASTQ format.");
}


//////////////////////////////////////////////////////////////////////////////////////
// read FASTA
//////////////////////////////////////////////////////////////////////////////////////
void Assemble::readFasta(const string &inputFilename, FILE **outputMT, const unsigned long long numThread) const
{
    ifstream ifs(inputFilename);
    SEQ seq;
    string oneLine;
    string read = "";
    unsigned i = 0;

    for (unsigned long long i = 0; i < numThread; ++i)
        fseek(outputMT[i], 0, SEEK_END);
    while (ifs && getline(ifs, oneLine)) {
        if (oneLine[0] == '>') {
            break;
        }
    }

    while (ifs && getline(ifs, oneLine)) {
        if (oneLine[0] != '>') {
            read += oneLine;
        } else {
            if (read.length() != 0) {
                seq.convertFromString(read);
                seq.writeTemporaryFile(outputMT[i]);
                read = "";
                i = (i + 1) % numThread;
            }
        }
    }
    seq.convertFromString(read);
    seq.writeTemporaryFile(outputMT[i]);

    ifs.close();
}


//////////////////////////////////////////////////////////////////////////////////////
// read FASTQ
//////////////////////////////////////////////////////////////////////////////////////
void Assemble::readFastq(const string &inputFilename, FILE **outputMT, const unsigned long long numThread) const
{
    ifstream ifs(inputFilename);
    SEQ seq;
    string oneLine;
    string read = "";
    unsigned i = 0;
    bool flag = true;

    for (unsigned long long i = 0; i < numThread; ++i)
        fseek(outputMT[i], 0, SEEK_END);
    while (ifs && getline(ifs, oneLine)) {
        if (oneLine[0] == '@') {
            break;
        }
    }

    while (ifs && getline(ifs, oneLine)) {
        if (oneLine.length() != 0) {
            if (oneLine[0] != '@') {
                if (flag && oneLine[0] != '+')
                    read += oneLine;
                else
                    flag = false;
            } else {
                if (read.length() != 0) {
                    seq.convertFromString(read);
                    seq.writeTemporaryFile(outputMT[i]);
                    read = "";
                    i = (i + 1) % numThread;
                }
                flag = true;
            }
        }
    }
    seq.convertFromString(read);
    seq.writeTemporaryFile(outputMT[i]);

    ifs.close();

}


