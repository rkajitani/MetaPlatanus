#ifndef CLUSTER_SCAFFOLD_H
#define CLUSTER_SCAFFOLD_H

#include "baseCommand.h"
#include "scaffold.h"
#include <memory>
#include <sstream>

class ClusterScaffold : public BaseCommand
{
private:
    std::string platanusBinary;


    virtual bool checkFileEnough(void)
    {
        if (optionSingleArgs["-c"] == "") {
            std::cerr << "Error: not specified contig file!!" << std::endl;
            return false;
        }
        return true;
    }

    virtual int checkOtherOption(char *argv) const
    {
        return 0;
    }


public:
    ClusterScaffold();
    ClusterScaffold(const ClusterScaffold &) = delete;
    ClusterScaffold &operator=(const ClusterScaffold &) = delete;
    ~ClusterScaffold() = default;

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

	void execClusterMode(void);
	void execPostClusteringScaffoldMode(void);
};



#endif
