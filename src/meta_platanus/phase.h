#ifndef PHASE_H
#define PHASE_H

#include "baseCommand.h"
#include <memory>
#include <vector>
#include <string>
#include <sstream>

class Phase : public BaseCommand
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
    Phase();
    Phase(const Phase &) = delete;
    Phase &operator=(const Phase &) = delete;
    ~Phase() = default;

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


	void setDirectoryName(const std::string newDirectoryName);
    void createDirectory(void) const;
	void execSolveDBG(const unsigned long times, const unsigned long numIterate);
    void execGapClose(void);
	void moveAndConcatenateFinalRoundResult(const std::string intermediateDirectoryName, const unsigned long times);
};



#endif
