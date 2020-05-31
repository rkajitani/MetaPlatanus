#ifndef CLUSTER_FILL_H
#define CLUSTER_FILL_H

#include "baseCommand.h"
#include "scaffold.h"
#include <memory>
#include <sstream>

class ClusterFill : public BaseCommand
{
private:
    std::string directoryName;
    std::string platanusBinary;


    virtual bool checkFileEnough(void)
    {
        if (optionMultiArgs["-c"].empty()) {
            std::cerr << "Error: No contig fasta files (-c) were specified!" << std::endl;
            return false;
        }
        return true;
    }

    virtual int checkOtherOption(char *argv) const
    {
        return 0;
    }


public:
    ClusterFill();
    ClusterFill(const ClusterFill &) = delete;
    ClusterFill &operator=(const ClusterFill &) = delete;
    ~ClusterFill() = default;

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

	void setDirectoryName(const std::string newDirectoryName) { this->directoryName = newDirectoryName; }
    void createDirectory(void) const;
	void execGapClose(void);
	void execSolveDBG(void);
	void processResultFile(void);
};



#endif
