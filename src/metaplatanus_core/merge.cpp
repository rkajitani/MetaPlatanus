#include "merge.h"




//////////////////////////////////////////////////////////////////////////////////////
// default constructor
//////////////////////////////////////////////////////////////////////////////////////
ContigMerger::ContigMerger()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-c"] = "1";
    optionSingleArgs["-k"] = "1.5";
    optionSingleArgs["-l"] = "2.0";
    optionSingleArgs["-u"] = "0.02";
    optionSingleArgs["-d"] = "0.5";
    optionSingleArgs["-m"] = "16";
    optionMultiArgs["-f"] = std::vector<std::string>();
    optionSingleArgs["-tmp"] = ".";
    optionBool["-fixed_k_len"] = false;
}


//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void ContigMerger::usage(void) const
{
    std::cerr << "\nUsage: platanus_b merge [Options]\n"
              << "Options:\n"
              << "    -o STR               : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -f FILE1 [FILE2 ...] : contig or scaffold_file (fasta format)\n"
              << "    -c INT               : minimum coverage (default " << optionSingleArgs.at("-c") << ")\n"
              << "    -k FLOAT             : k-mer size factor (k = FLOAT * read_length) (default read from fasta-header)\n"
              << "    -l FLOAT             : minimum length factor (minimum_length = FLOAT * read_length) (default " << optionSingleArgs.at("-l") << ")\n"
              << "    -u FLOAT             : maximum difference for bubble crush (identity, default " << optionSingleArgs.at("-u") << ")\n"
              << "    -d FLOAT             : maximum difference for branch cutting (coverage ratio, default " << optionSingleArgs.at("-d") << ")\n"
              << "    -m INT               : memory limit for making kmer distribution (GB, >=1, default " << optionSingleArgs.at("-m") << ")\n"
			  << "    -fixed_k_len         : use k-mer size in input fasta file (\"maxK\" in the first line) and do not filter short contigs\n"
              << "    -tmp DIR             : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n\n\n"

              << "Output:\n"
              << "    PREFIX_merged.fa\n"
              << std::endl;
}




//////////////////////////////////////////////////////////////////////////////////////
// exec divide
//////////////////////////////////////////////////////////////////////////////////////
void ContigMerger::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
    this->readLength = this->contig.getReadLengthFromFastaHeader(optionMultiArgs["-f"][0]);

    unsigned long long kmerLength;
	if (!optionBool["-fixed_k_len"]) {
		kmerLength = atof(optionSingleArgs["-k"].c_str()) * readLength + 0.5;
	}
	else {
		kmerLength = contig.getMaxKFromFastaHeader(optionMultiArgs["-f"][0]);
	}

    if (kmerLength <= 32) {
        Counter<Kmer31> counter;
        BruijnGraph<Kmer31> graph;
        this->exec2_ForIntegrateKmer(counter, graph, kmerLength);
    } else if (kmerLength <= 64) {
        Counter<KmerN<Binstr63> > counter;
        BruijnGraph<KmerN<Binstr63> > graph;
        this->exec2_ForIntegrateKmer(counter, graph, kmerLength);
    } else if (kmerLength <= 96) {
        Counter<KmerN<Binstr95> > counter;
        BruijnGraph<KmerN<Binstr95> > graph;
        this->exec2_ForIntegrateKmer(counter, graph, kmerLength);
    } else if (kmerLength <= 128) {
        Counter<KmerN<Binstr127> > counter;
        BruijnGraph<KmerN<Binstr127> > graph;
        this->exec2_ForIntegrateKmer(counter, graph, kmerLength);
    } else if (kmerLength <= 160) {
        Counter<KmerN<Binstr159> > counter;
        BruijnGraph<KmerN<Binstr159> > graph;
        this->exec2_ForIntegrateKmer(counter, graph, kmerLength);
    } else {
        Counter<KmerN<binstr_t> > counter;
        BruijnGraph<KmerN<binstr_t> > graph;
        this->exec2_ForIntegrateKmer(counter, graph, kmerLength);
    }
}



template <typename KMER>
void ContigMerger::exec2_ForIntegrateKmer(Counter<KMER> &counter, BruijnGraph<KMER> &graph, const unsigned long long kmerLength)
{
    unsigned long long memory = atoi(optionSingleArgs["-m"].c_str()) * static_cast<unsigned long long>(1000000000);
    unsigned long long contigReadLength = platanus::Contig::getReadLengthFromFastaHeader(optionMultiArgs["-f"][0]);

    unsigned long long minContigLength;
	if (!optionBool["-fixed_k_len"]) {
		minContigLength = atof(optionSingleArgs["-l"].c_str()) * contigReadLength + 0.5;
	}
	else {
		minContigLength = kmerLength;
	}

    graph.setBubbleAndBranch(atof(optionSingleArgs["-u"].c_str()), atof(optionSingleArgs["-d"].c_str()));
    counter.setKmerLength(kmerLength);
    graph.setKmerLength(kmerLength);
    for (auto itr = optionMultiArgs["-f"].begin(), end = optionMultiArgs["-f"].end(); itr != end; ++itr) {
        contig.readFastaCoverageCutN(*itr, minContigLength);
        //contig.readFastaCoverage(*itr);
    }
    unsigned long long doubleHashSize = counter.makeKmerReadDistributionFromContig(contig, minContigLength, atoi(optionSingleArgs["-c"].c_str()), memory);
    contig.clear();

    FILE *sortedKeyFP = counter.sortedKeyFromKmerFile(0);
    counter.loadKmer(0, doubleHashSize);
    graph.makeInitialBruijnGraph(counter, sortedKeyFP);
    fclose(sortedKeyFP);
    counter.deleteAllTable();

    graph.cutBranchIterative(1);
    graph.crushBubbleIterative();

    FILE *contigFP = platanus::makeTemporaryFile();
    graph.saveContigSimple(contigFP);

    FILE *junctionFP = platanus::makeTemporaryFile();
    graph.saveJunction(junctionFP);


    std::string outputFilename = optionSingleArgs["-o"];
    outputFilename += "_merged.fa";
    platanus::printContig(outputFilename, contigFP, 1.0, readLength, kmerLength, "seq");

    outputFilename = optionSingleArgs["-o"];
    outputFilename += "_mergedJunctionKmer.fa";
    platanus::printContig(outputFilename, junctionFP, 1.0, this->readLength, kmerLength, "junction");
	fclose(junctionFP);
}
