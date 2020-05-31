#include "clusterFill.h"
#include <sys/stat.h>
#include <sys/types.h>


using std::vector;
using std::string;


//////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////
ClusterFill::ClusterFill()
{
    optionSingleArgs["-o"] = "out";
    optionMultiArgs["-c"] = vector<string>();
    optionSingleArgs["-t"] = "1";
    pairedEndSingleFileType.emplace_back("-ip");
    pairedEndSingleFileType.emplace_back("-op");
    pairedEndPairFileType.emplace_back("-IP");
    pairedEndPairFileType.emplace_back("-OP");
    optionSingleArgs["-tmp"] = ".";

    optionMultiArgs["-s"] = vector<string>(1);
    optionMultiArgs["-s"][0] = "32";
}




//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void ClusterFill::usage(void) const
{

    std::cerr << "\nUsage: meta_platanus cluster_fill [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file and directory (do not use \"/\", default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : clustered (binned) contig files (fasta format)\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default " << optionMultiArgs.at("-s")[0] << ")\n"
              << "    -t INT                             : number of threads (default " << optionSingleArgs.at("-t") << ")\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n"
			  << "\n\n"
              << "Outputs:\n"
              << "    PREFIX_filledClusters\n"
              << "    PREFIX_filledClusters_all.fa\n"
              << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////
// exec gap close
//////////////////////////////////////////////////////////////////////////////////////
void ClusterFill::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());

	string intermediateDirectoryName = optionSingleArgs["-o"];
	intermediateDirectoryName += "_clusterFillIntermediateResults";
	this->setDirectoryName(intermediateDirectoryName);
	this->createDirectory();

	this->execGapClose();
	this->execSolveDBG();
	this->processResultFile();

    std::cerr << "cluster_fill completed!" << std::endl;
}


void ClusterFill::createDirectory(void) const
{
    struct stat st;

    std::ostringstream oss;
    oss << "rm -rf " << directoryName << " >/dev/null 2>&1;";
	system(oss.str().c_str());

    if (stat(directoryName.c_str(), &st) != 0) {
        int returnValue = mkdir(directoryName.c_str(), 0755);
        if (returnValue != 0) {
            throw platanus::CreateDirError(directoryName);
        }
    }
}


void ClusterFill::execGapClose(void)
{
	vector<string> singleArgOptionList = {"-t", "-tmp"};
	vector<string> multiArgOptionList = {"-c", "-s"};

	string outputPrefix = this->directoryName;
	outputPrefix += "/";
	outputPrefix += optionSingleArgs["-o"];

    std::ostringstream oss;
    oss << this->platanusBinary << " gap_close -cluster -extend";

    for (auto itr = optionPairFile.begin(); itr != optionPairFile.end(); ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (!(itr->fileSecond.empty())) {
            oss << " " << itr->fileSecond;
        }
    }

	for (auto optItr = singleArgOptionList.begin(); optItr != singleArgOptionList.end(); ++optItr) {
		if (!(optionSingleArgs[*optItr].empty())) {
			oss << " " << *optItr << " " << optionSingleArgs[*optItr];
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

	oss << " -o " << outputPrefix;
    oss << " 2>" << outputPrefix << ".gapCloseLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::GapError();
    }
}


void ClusterFill::execSolveDBG(void)
{
	vector<string> singleArgOptionList = {"-o", "-t", "-tmp"};

	string outputPrefix = this->directoryName;
	outputPrefix += "/";
	outputPrefix += optionSingleArgs["-o"];

    std::ostringstream oss;
    oss << this->platanusBinary << " solve_DBG -trim_cluster";

	oss << " -c"
		<< " " << outputPrefix << "*_gapClosed_*.fa"
		<< " " << outputPrefix << "_extraContig.fa";

	for (auto optItr = singleArgOptionList.begin(); optItr != singleArgOptionList.end(); ++optItr) {
		if (!(optionSingleArgs[*optItr].empty())) {
			oss << " " << *optItr << " " << optionSingleArgs[*optItr];
		}
	}
	
	oss << " -o " << outputPrefix
    	<< " 2>" << outputPrefix << ".solveDBGLog";


    if (system(oss.str().c_str()) != 0) {
        throw platanus::SolveDBGError();
    }
}


void ClusterFill::processResultFile(void)
{
    std::ostringstream oss;

    oss << "rm -rf " << optionSingleArgs["-o"] << "_filledClusters >/dev/null 2>&1;"
    	<< "mv " << this->directoryName << "/" << optionSingleArgs["-o"] << "_trimmedClusters " << optionSingleArgs["-o"] << "_filledClusters;"
    	<< "cat " << optionSingleArgs["-o"] << "_filledClusters/*.fa >" << optionSingleArgs["-o"] << "_filledClusters_all.fa";

	system(oss.str().c_str());
}
