#ifndef CONSENSUS_H
#define CONSENSUS_H

#include "baseCommand.h"
#include <memory>
#include <vector>
#include <string>
#include <sstream>

class Consensus : public BaseCommand
{
private:
    std::string directoryName;
    std::string previousDirectoryName;
    std::string platanusBinary;


    virtual bool checkFileEnough(void)
    {
        if (optionMultiArgs["-c"].empty()) {
            std::cerr << "Error: No contig fasta files (-c) were specified!" << std::endl;
            return false;
        }
        if (optionSingleArgs["-k"] == "") {
            std::cerr << "Error: both -k(k-mer-occurrence-file) is not specified!!" << std::endl;
            return false;
        }
        return true;
    }

    virtual int checkOtherOption(char *argv) const
    {
        return 0;
    }


public:
    Consensus();
    Consensus(const Consensus &) = delete;
    Consensus &operator=(const Consensus &) = delete;
    ~Consensus() = default;

    void usage(void) const;
    void exec(void);

    virtual bool parseArgs(int argc, char **argv)
    {
        if (BaseCommand::parseArgs(argc, argv)) {
            this->platanusBinary = argv[0];
            return true;
        }
        return false;
    }

	void execSolveDBG(void);
};



#endif
