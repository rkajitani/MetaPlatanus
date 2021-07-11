#ifndef MERGE_H
#define MERGE_H

#include "common.h"
#include "counter.h"
#include "graph.h"


class ContigMerger : public BaseCommand
{
private:
    platanus::Contig contig;
    unsigned long long readLength;
    unsigned long long contigmaxK;

public:
    ContigMerger();
    ContigMerger(const ContigMerger &) = delete;
    ~ContigMerger() = default;

    void usage(void) const;
    void exec(void);
    template <typename KMER> void exec2_ForIntegrateKmer(Counter<KMER> &counter, BruijnGraph<KMER> &graph, const unsigned long long kmerLength);


    bool checkFileEnough(void)
    {
        if (optionMultiArgs["-f"].size() == 0) {
            std::cerr << "Error: not specified contig file!!" << std::endl;
            return false;
        }
        return true;
    }

    int checkOtherOption(char *argv) const
    {
        return 0;
    }

};





#endif

