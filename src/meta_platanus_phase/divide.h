#ifndef DIVIDE_H
#define DIVIDE_H

#include "common.h"
#include "counter.h"
#include <algorithm>
#include <memory>
#include <vector>
#include <deque> 
#include <float.h>


class ContigDivider : public BaseCommand
{
private:

    class OccurrenceArray
    {
    private:
        unsigned long long kmerLength;
        unsigned long long numSeq;
        std::vector<std::vector<unsigned long long> > occurrenceArray;

    public:
        OccurrenceArray(const unsigned long long kmer, const platanus::Contig &contig): kmerLength(kmer), numSeq(contig.numSeq), occurrenceArray()
        {
            occurrenceArray.resize(numSeq);
            for (unsigned long long i = 0; i < numSeq; ++i) {
                occurrenceArray[i].resize(contig.seq[i].length - kmerLength + 1);
            }
        }

        unsigned long long getKmerLength(void) const { return kmerLength; }
        unsigned long long getNumSeq(void) const { return numSeq; }
        unsigned long long getNumKmer(const unsigned long long id) const {return occurrenceArray[id].size(); }
        std::vector<unsigned long long> getOccurrenceArray(const unsigned long long id) const { return occurrenceArray[id]; }
        unsigned long long getOccurrence(const unsigned long long id, const unsigned long long pos) const { return occurrenceArray[id][pos]; }

        void setOccurrenceArray(const unsigned long long id, const unsigned long long pos, const unsigned long long value)
        {
            occurrenceArray[id][pos] = value;
        }
    };

    platanus::Contig contig;
    unsigned long long readLength;
    unsigned long long contigMaxK;
    double coverageRatioThreshold;
    std::unique_ptr<OccurrenceArray> contigItself;
    std::vector<std::deque<unsigned long long> > breakPoint;

public:
    ContigDivider();
    ContigDivider(const ContigDivider &) = delete;
    ~ContigDivider() = default;

    void usage(void) const;
    void exec(void);
    unsigned long long getKmerLengthFromBinary(const std::string &filename);
    template <typename KMER> void getOccurrenceArray(const Counter<KMER> &counter);
    void decideContigBreakPoint(void);
    unsigned long long findMedian(const std::vector<unsigned long long> &vec) const;
    void divideAndPrintContig(const std::string &filename);
    unsigned long long calcAverageCoverage(const unsigned long long contigID, const unsigned long long start, const unsigned long long end);


    bool checkFileEnough(void)
    {
        if (optionMultiArgs["-c"].size() == 0) {
            std::cerr << "Error: not specified contig file!!" << std::endl;
            return false;
        }
        if (optionSingleArgs["-k"] == "") {
            std::cerr << "Error: not specified kmer occurrence file!!!" << std::endl;
            return false;
        }
        return true;
    }

    int checkOtherOption(char *argv) const
    {
        return 0;
    }


    inline bool withinRatioThreshold(const unsigned long long val1, const unsigned long long val2) const
    {
        return (val1 * this->coverageRatioThreshold >= val2) && (val2 * this->coverageRatioThreshold >= val1);
    }

};





#endif

