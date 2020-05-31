#ifndef ITERATE_H
#define ITERATE_H

#include "baseCommand.h"
#include "scaffold.h"
#include "gapClose.h"
#include "divide.h"
#include "merge.h"
#include <memory>
#include <sstream>

class IterateScaffold : public BaseCommand
{
private:
    static const std::string CONTIG_FOOTER;
    static const std::string SCAF_FOOTER;
    static const std::string GAP_FOOTER;
    static const std::string EX_FOOTER;
    static const std::string DIV_FOOTER;
    static const std::string MERGE_FOOTER;
    static const std::string ITERATION_FOOTER;
    std::string directoryName;
    std::string previousDirectoryName;
    std::string platanusBinary;


    virtual bool checkFileEnough(void)
    {
        if (optionMultiArgs["-c"].empty()) {
            std::cerr << "Error: not specified contig file!!" << std::endl;
            return false;
        }
        if (optionSingleArgs["-k"] == "") {
            std::cerr << "Error: not specified kmer occurrence file!!" << std::endl;
            return false;
        }
        return true;
    }

    virtual int checkOtherOption(char *argv) const
    {
        return 0;
    }


public:
    IterateScaffold();
    IterateScaffold(const IterateScaffold &) = delete;
    IterateScaffold &operator=(const IterateScaffold &) = delete;
    ~IterateScaffold() = default;

    void usage(void) const;
    void exec(void);

    void setDirectoryName(const unsigned long long times)
    {
        std::ostringstream oss;
        oss << optionSingleArgs["-o"] << times;
        previousDirectoryName = directoryName;
        directoryName = oss.str();
    }

    virtual bool parseArgs(int argc, char **argv)
    {
        if (BaseCommand::parseArgs(argc, argv)) {
            this->platanusBinary = argv[0];
            return true;
        }
        return false;
    }


    void createDirectory(void) const;
    void createContig(const unsigned long long times);
    void execDivideMode(void);
    void execMergeMode(void);
    void execScaffoldMode(const unsigned long times, const std::string &mode);
    void execGapCloseMode(const std::string &mode);
	void execClusterMode(void);
	void execPostClusteringScaffoldMode(void);

};



#endif
