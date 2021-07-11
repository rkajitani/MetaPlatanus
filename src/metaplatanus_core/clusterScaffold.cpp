#include "clusterScaffold.h"
#include <sys/stat.h>
#include <sys/types.h>


using std::vector;
using std::string;


//////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////
ClusterScaffold::ClusterScaffold()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-sl"] = "";
    optionSingleArgs["-c"] = "";
    optionSingleArgs["-t"] = "1";
    pairedEndSingleFileType.emplace_back("-ip");
    pairedEndSingleFileType.emplace_back("-op");
    pairedEndPairFileType.emplace_back("-IP");
    pairedEndPairFileType.emplace_back("-OP");
    optionMultiArgs["-p"] = vector<string>();
    optionMultiArgs["-x"] = vector<string>();
    optionMultiArgs["-X"] = vector<string>();
    optionSingleArgs["-cl"] = "";
    optionSingleArgs["-n"] = "";
    optionSingleArgs["-tmp"] = ".";
    optionSingleArgs["-mga"] = "mga";
}




//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void ClusterScaffold::usage(void) const
{

    std::cerr << "\nUsage: metaplatanus_core cluster_scaffold [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file and directory (do not use \"/\", default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1                           : contig (or scaffold) file (fasta format)\n"
              << "    -sl INT                            : -l value of \"scaffold\" step\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -cl INT                            : -l value of \"cluster\" step\n"
              << "    -n INT                             : -n value of \"cluster\"\n"
              << "    -t INT                             : number of threads (default " << optionSingleArgs.at("-t") << ")\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n"
              << "    -mga FILE                          : name of mga executable file (default " << optionSingleArgs.at("-mga") << ")\n"
              << "    -keep_file                         : keep all intermediate files (default, off)\n"
			  << "\n\n"
              << "Outputs:\n"
              << "    [contig_file(-c)].clusters.tsv\n"
              << "    PREFIX_finalClusters\n"
              << "    PREFIX_finalClusters_all.fa\n"
              << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////
// exec gap close
//////////////////////////////////////////////////////////////////////////////////////
void ClusterScaffold::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());

	this->execClusterMode();
	this->execPostClusteringScaffoldMode();

    std::cerr << "clustering and scaffolding completed!" << std::endl;
}


void ClusterScaffold::execClusterMode(void)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " cluster"
        << " -t " << optionSingleArgs["-t"]
        << " -c " << optionSingleArgs["-c"];

    if (optionSingleArgs["-cl"] != "") {
        oss << " -l " << optionSingleArgs["-cl"];
    }
    if (optionSingleArgs["-n"] != "") {
        oss << " -n " << optionSingleArgs["-n"];
    }
    if (optionSingleArgs["-mga"] != "") {
        oss << " -mga " << optionSingleArgs["-mga"];
    }

    oss << " 2>" << optionSingleArgs["-o"] << ".clusterLog";

    if (system(oss.str().c_str()) != 0) {
//        throw platanus::ClusterError();
    }
}


void ClusterScaffold::execPostClusteringScaffoldMode(void)
{
	vector<string> multiArgOptionList = {"-p", "-x", "-X"};

    std::ostringstream oss;
    oss << this->platanusBinary << " solve_DBG -unphase"
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -c " << optionSingleArgs["-c"]
        << " -g " << optionSingleArgs["-c"] << ".clusters.tsv"
        << " -o " << optionSingleArgs["-o"];
    if (optionSingleArgs["-sl"] != "") {
        oss << " -l " << optionSingleArgs["-sl"];
    }

    for (auto itr = optionPairFile.begin(), end = optionPairFile.end(); itr != end; ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (itr->fileSecond != "") {
            oss << " " << itr->fileSecond;
        }
    }

	for (auto optItr = multiArgOptionList.begin(); optItr != multiArgOptionList.end(); ++optItr) {
		if (!(optionMultiArgs[*optItr].empty())) {
			oss << " " << *optItr;
			for (auto argItr = optionMultiArgs[*optItr].begin(); argItr != optionMultiArgs[*optItr].end(); ++argItr) {
				oss << " " << *argItr;
			}
		}
	}

    oss << " 2>" << optionSingleArgs["-o"] << ".postClustScafLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::ScaffoldError();
    }

    std::ostringstream ossForMv;
    ossForMv << "rm -rf " << optionSingleArgs["-o"] << "_finalClusters >/dev/null 2>&1;";
    ossForMv << "mv " << optionSingleArgs["-o"] << "_extendedClusters " << optionSingleArgs["-o"] << "_finalClusters";
	system(ossForMv.str().c_str());

    std::ostringstream ossForCat;
    ossForCat << "cat " << optionSingleArgs["-o"] << "_finalClusters/*.fa >" << optionSingleArgs["-o"] << "_finalClusters_all.fa";
	system(ossForCat.str().c_str());

	if (!optionBool["-keep_file"]) {
		std::ostringstream rmOss;
		rmOss << "rm -f " << optionSingleArgs["-o"] << "_longReadAlignment.tsv >/dev/null 2>&1";
		system(rmOss.str().c_str());
	}
}
