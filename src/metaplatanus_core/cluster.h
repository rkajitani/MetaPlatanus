#ifndef CLUSTER_H
#define CLUSTER_H

#include "baseCommand.h"
#include "common.h"
#include "mgc.h"
#include <string>


class Cluster : public BaseCommand
{
private:
	static const unsigned KMER_LENGTH_FOR_CLUSTERING;

public:
    Cluster();
    Cluster(const Cluster &) = delete;
    ~Cluster() = default;

    void usage(void) const;
    void exec(void);

    bool checkFileEnough(void)
    {
        if (optionSingleArgs["-c"] == "") {
            std::cerr << "Error: not specified scaffold(contig) file!!" << std::endl;
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

