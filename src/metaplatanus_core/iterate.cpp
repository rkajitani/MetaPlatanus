#include "iterate.h"
#include <sys/stat.h>
#include <sys/types.h>


const std::string IterateScaffold::CONTIG_FOOTER = "_contig.fa";
const std::string IterateScaffold::SCAF_FOOTER = "_consensusScaffold.fa";
const std::string IterateScaffold::GAP_FOOTER = "_gapClosed_consensusScaffold.fa";
const std::string IterateScaffold::EX_FOOTER = "_extraContig.fa";
const std::string IterateScaffold::COV_TRIM_FOOTER = "_cov_trimmed.fa";
const std::string IterateScaffold::DIV_FOOTER = "_divided.fa";
const std::string IterateScaffold::MERGE_FOOTER = "_merged.fa";
const std::string IterateScaffold::ITERATION_FOOTER = "_iterativeAssembly.fa";


//////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////
IterateScaffold::IterateScaffold()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-i"] = "5";
    optionSingleArgs["-f"] = "1";
    optionSingleArgs["-m"] = "16";
    optionSingleArgs["-l"] = "";
    optionSingleArgs["-r"] = "";
    optionSingleArgs["-k"] = "";
    optionSingleArgs["-t"] = "1";
    optionMultiArgs["-c"] = std::vector<std::string>();
    optionMultiArgs["-p"] = std::vector<std::string>();
    optionMultiArgs["-ont"] = std::vector<std::string>();
    optionMultiArgs["-gc"] = std::vector<std::string>();
    optionMultiArgs["-x"] = std::vector<std::string>();
    optionMultiArgs["-X"] = std::vector<std::string>();
    pairedEndSingleFileType.emplace_back("-ip");
    pairedEndSingleFileType.emplace_back("-op");
    pairedEndPairFileType.emplace_back("-IP");
    pairedEndPairFileType.emplace_back("-OP");
    optionSingleArgs["-tmp"] = ".";
    optionSingleArgs["-mga"] = "mga";
	optionBool["-keep_file"] = false;
}




//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void IterateScaffold::usage(void) const
{

    std::cerr << "\nUsage: metaplatanus_core iterate [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file and directory (do not use \"/\", default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : contig (or scaffold) file (fasta format)\n"
              << "    -k FILE                            : k-mer-occurrence file (binary format)\n"
              << "    -i INT                             : number of iterations (default " << optionSingleArgs.at("-i") << ")\n"
              << "    -f INT                             : continue from the INT-th iteration (default " << optionSingleArgs.at("-f") << ")\n"
              << "    -l INT                             : -l value of \"scaffold\" step\n"
              << "    -l INT                             : -l value of \"scaffold\" step\n"
              << "    -r FLOAT                           : -r value of \"scaffold\" and \"cov_trim\" step\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -ont FILE1 [FILE2 ...]             : Oxford Nanopore long-read file (fasta or fastq)\n"
              << "    -p FILE1 [FILE2 ...]               : PacBio long-read file (fasta or fastq)\n"
              << "    -gc FILE1 [FILE2 ...]              : Guiding contig file; i.e. other assemblies, synthetic long-reads or corrected reads (fasta or fastq)\n"
              << "    -x PAIR1 [PAIR2 ...]               : linked-reads files (paired-ends, 10x Genomics) (interleaved, fasta or fastq)\n"
              << "    -X FWD1 REV1 [FWD2 REV2 ...]       : linked-reads files (paired-ends, 10x Genomics) (separate forward and reverse files, fasta or fastq)\n"
              << "    -t INT                             : number of threads (default " << optionSingleArgs.at("-t") << ")\n"
              << "    -m INT                             : memory limit for making kmer distribution (GB, >=1, default " << optionSingleArgs.at("-m") << ")\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n"
              << "    -keep_file                         : keep all intermediate files (default, off)\n\n\n"

              << "Input format:\n"
              << "    Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -ip, -IP, -op, -OP, -ont, -p and -gc options.\n"
              << "\n\n"

              << "Outputs:\n"
              << "    PREFIX_iterativeAssembly.fa\n"
              << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////
// exec gap close
//////////////////////////////////////////////////////////////////////////////////////
void IterateScaffold::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());

    const unsigned long long iterateTimes = atoi(optionSingleArgs["-i"].c_str());

    for (unsigned long long times = 1; times <= iterateTimes; ++times) {
        this->setDirectoryName(times);
		if (times < atoi(optionSingleArgs["-f"].c_str()))
			continue;

        this->createDirectory();
        this->createContig(times);
        this->execDivideMode();
        this->execCovTrimMode();

		std::string scaffoldMode;
		if (times <= 3) {
			scaffoldMode = "-no_scaffold";
		}
		else if (times <= 4) {
			scaffoldMode = "-no_barcode_scaffold";
		}
		else {
			scaffoldMode = "";
		}
		this->execScaffoldMode(times, scaffoldMode);

        if (times == iterateTimes) {
            this->execGapCloseMode("");
        } else {
            this->execGapCloseMode("-extend");
        }
    }
	this->execFinalDivideMode();
//	this->execFinalGapTrimMode();

    std::ostringstream oss;
//    oss << "mv " << optionSingleArgs["-o"] << "_final_gapTrimmed.fa " << optionSingleArgs["-o"] << ITERATION_FOOTER;
    oss << "mv " << optionSingleArgs["-o"] << "_final_divided.fa " << optionSingleArgs["-o"] << ITERATION_FOOTER;
    if (system(oss.str().c_str()) != 0) {
        throw platanus::CreateLinkError();
    }

	if (!optionBool["-keep_file"]) {
		std::ostringstream rmOss;
		rmOss << "rm -f " << optionSingleArgs["-o"] << "_final_longReadAlignment.tsv >/dev/null 2>&1";
		system(rmOss.str().c_str());
	}
}


void IterateScaffold::createDirectory(void) const
{
    struct stat st;

    if (stat(directoryName.c_str(), &st) != 0) {
        int returnValue = mkdir(directoryName.c_str(), 0755);
        if (returnValue != 0) {
//            throw platanus::CreateDirError(directoryName);
        }
    }
}


void IterateScaffold::createContig(const unsigned long long times)
{
    if (times == 1) {
        std::ostringstream oss;
        oss << "cat";
		for (auto itr = optionMultiArgs["-c"].begin(); itr != optionMultiArgs["-c"].end(); ++itr) {
			oss << " " << *itr;
		}
        oss << " >" << this->directoryName << "/" << optionSingleArgs["-o"] << CONTIG_FOOTER;
        if (system(oss.str().c_str()) != 0) {
            throw platanus::CreateLinkError();
        }
    } else {
        this->execMergeMode();
        std::ostringstream oss;
        oss << "cat " << this->directoryName << "/" << optionSingleArgs["-o"] << "_merged.fa"
            << " " << this->directoryName << "/" << optionSingleArgs["-o"] << "_mergedJunctionKmer.fa"
            << " >" << this->directoryName << "/" << optionSingleArgs["-o"] << CONTIG_FOOTER
            << "; rm " <<  this->directoryName << "/" << optionSingleArgs["-o"] << "_merged.fa"
            << " " << this->directoryName << "/" << optionSingleArgs["-o"] << "_mergedJunctionKmer.fa";
 
        if (system(oss.str().c_str()) != 0) {
            throw platanus::CreateLinkError();
        }
    }
}


void IterateScaffold::execDivideMode(void)
{
	std::vector<std::string> multiArgOptionList = {};
	multiArgOptionList.push_back("-p");
	multiArgOptionList.push_back("-ont");
	multiArgOptionList.push_back("-gc");

    std::ostringstream oss;

	if (optionMultiArgs["-p"].empty() && optionMultiArgs["-ont"].empty() && optionMultiArgs["-ont"].empty()) {
		oss << "cp " << this->directoryName << "/"<< optionSingleArgs["-o"] << CONTIG_FOOTER << " "
			<< this->directoryName << "/" << optionSingleArgs["-o"] << DIV_FOOTER;
	}
	else {
		oss << this->platanusBinary << " divide"
			<< " -t " << optionSingleArgs["-t"]
			<< " -tmp " << optionSingleArgs["-tmp"]
			<< " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << CONTIG_FOOTER
			<< " -o " << this->directoryName << "/" << optionSingleArgs["-o"];

		for (auto optItr = multiArgOptionList.begin(); optItr != multiArgOptionList.end(); ++optItr) {
			if (!(optionMultiArgs[*optItr].empty())) {
				oss << " " << *optItr;
				for (auto argItr = optionMultiArgs[*optItr].begin(); argItr != optionMultiArgs[*optItr].end(); ++argItr) {
					oss << " " << *argItr;
				}
			}
		}
		oss << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".divLog";
	}

    if (system(oss.str().c_str()) != 0) {
        throw platanus::DivideError();
    }

	if (!optionBool["-keep_file"]) {
		std::ostringstream rmOss;
		rmOss << "rm -f " << this->directoryName << "/" << optionSingleArgs["-o"] << CONTIG_FOOTER << ">/dev/null 2>&1";
		system(rmOss.str().c_str());
	}
}


void IterateScaffold::execFinalDivideMode(void)
{
	std::vector<std::string> multiArgOptionList = {};
	multiArgOptionList.push_back("-p");
	multiArgOptionList.push_back("-ont");
	multiArgOptionList.push_back("-gc");

    std::ostringstream oss;

	if (optionMultiArgs["-p"].empty() && optionMultiArgs["-ont"].empty() && optionMultiArgs["-ont"].empty()) {
		oss << "mv " << this->directoryName << "/"  << optionSingleArgs["-o"] << GAP_FOOTER << " "
			<< optionSingleArgs["-o"] << "_final_divided.fa";
	}
	else {
		oss << this->platanusBinary << " divide -long"
			<< " -t " << optionSingleArgs["-t"]
			<< " -tmp " << optionSingleArgs["-tmp"]
			<< " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << GAP_FOOTER
			<< " -o " << optionSingleArgs["-o"] << "_final";

		for (auto optItr = multiArgOptionList.begin(); optItr != multiArgOptionList.end(); ++optItr) {
			if (!(optionMultiArgs[*optItr].empty())) {
				oss << " " << *optItr;
				for (auto argItr = optionMultiArgs[*optItr].begin(); argItr != optionMultiArgs[*optItr].end(); ++argItr) {
					oss << " " << *argItr;
				}
			}
		}
		oss << " 2>"  << optionSingleArgs["-o"] << ".finalDivLog";
	}

    if (system(oss.str().c_str()) != 0) {
        throw platanus::DivideError();
    }

	if (!optionBool["-keep_file"]) {
		std::ostringstream rmOss;
		rmOss << "rm -f " << this->directoryName << "/" << optionSingleArgs["-o"] << GAP_FOOTER << ">/dev/null 2>&1";
		system(rmOss.str().c_str());
	}
}


void IterateScaffold::execFinalGapTrimMode(void)
{
	std::vector<std::string> multiArgOptionList = {};
	multiArgOptionList.push_back("-p");
	multiArgOptionList.push_back("-ont");
	multiArgOptionList.push_back("-gc");

    std::ostringstream oss;

	if (optionMultiArgs["-p"].empty() && optionMultiArgs["-ont"].empty() && optionMultiArgs["-ont"].empty()) {
		oss << "mv " << this->directoryName << "/"  << optionSingleArgs["-o"] << GAP_FOOTER << " "
			<< optionSingleArgs["-o"] << "_final_divided.fa";
	}
	else {
		oss << this->platanusBinary << " solve_DBG -trim_gap_side"
			<< " -t " << optionSingleArgs["-t"]
			<< " -tmp " << optionSingleArgs["-tmp"]
			<< " -c "  << optionSingleArgs["-o"] << "_final_divided.fa"
			<< " -o " << optionSingleArgs["-o"] << "_final";

		for (auto optItr = multiArgOptionList.begin(); optItr != multiArgOptionList.end(); ++optItr) {
			if (!(optionMultiArgs[*optItr].empty())) {
				oss << " " << *optItr;
				for (auto argItr = optionMultiArgs[*optItr].begin(); argItr != optionMultiArgs[*optItr].end(); ++argItr) {
					oss << " " << *argItr;
				}
			}
		}
		oss << " 2>"  << optionSingleArgs["-o"] << ".finalGapTrimLog";
	}

    if (system(oss.str().c_str()) != 0) {
        throw platanus::SolveDBGError();
    }

	if (!optionBool["-keep_file"]) {
		std::ostringstream rmOss;
		rmOss << "rm -f " << optionSingleArgs["-o"] << "_final_divided.fa >/dev/null 2>&1";
		system(rmOss.str().c_str());
	}
}


void IterateScaffold::execCovTrimMode(void)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " cov_trim"
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << DIV_FOOTER
        << " -k " << optionSingleArgs["-k"]
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"];
    if (optionSingleArgs["-r"] != "") {
        oss << " -r " << optionSingleArgs["-r"];
    }
    oss << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".covTrimLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::CovTrimError();
    }

	if (!optionBool["-keep_file"]) {
		std::ostringstream rmOss;
		rmOss << "rm -f " << this->directoryName << "/" << optionSingleArgs["-o"] << DIV_FOOTER << ">/dev/null 2>&1";
		system(rmOss.str().c_str());
	}
}



void IterateScaffold::execMergeMode(void)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " merge"
        << " -m " << optionSingleArgs["-m"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -f " << this->previousDirectoryName << "/" << optionSingleArgs["-o"] << GAP_FOOTER
        << " " << this->previousDirectoryName << "/" << optionSingleArgs["-o"] << EX_FOOTER
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"]
        << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".mergeLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::MergeError();
    }

	if (!optionBool["-keep_file"]) {
		std::ostringstream rmOss;
		rmOss << "rm -f"
			<< " " << this->previousDirectoryName << "/" << optionSingleArgs["-o"] << GAP_FOOTER
			<< " " << this->previousDirectoryName << "/" << optionSingleArgs["-o"] << EX_FOOTER
			<< ">/dev/null 2>&1";
		system(rmOss.str().c_str());
	}
}


void IterateScaffold::execScaffoldMode(const unsigned long times, const std::string &mode)
{
	std::vector<std::string> multiArgOptionList = {"-x", "-X"};
	if (times > 1) {
		multiArgOptionList.push_back("-p");
		multiArgOptionList.push_back("-ont");
		multiArgOptionList.push_back("-gc");
	}

    std::ostringstream oss;
    oss << this->platanusBinary << " solve_DBG -unphase"
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -k " << optionSingleArgs["-k"]
        << " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << COV_TRIM_FOOTER
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"]
        << " " << mode;

    if (optionSingleArgs["-l"] != "") {
        oss << " -l " << optionSingleArgs["-l"];
    }
    if (optionSingleArgs["-r"] != "") {
        oss << " -r " << optionSingleArgs["-r"];
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

    oss << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".scafLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::ScaffoldError();
    }

	if (!optionBool["-keep_file"]) {
		std::ostringstream rmOss;
		rmOss << "rm -f"
			<< " " << this->directoryName << "/" << optionSingleArgs["-o"] << COV_TRIM_FOOTER
			<< " " << this->directoryName << "/" << optionSingleArgs["-o"] << "_longReadAlignment.tsv"
			<< ">/dev/null 2>&1";
		system(rmOss.str().c_str());
	}
}


void IterateScaffold::execGapCloseMode(const std::string &mode)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " gap_close"
		<< " -no_partial"
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << SCAF_FOOTER
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"]
        << " " << mode;
    for (auto itr = optionPairFile.begin(), end = optionPairFile.end(); itr != end; ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (itr->fileSecond != "") {
            oss << " " << itr->fileSecond;
        }
    }
    oss << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".gapLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::GapError();
    }

	if (!optionBool["-keep_file"]) {
		std::ostringstream rmOss;
		rmOss << "rm -f " << this->directoryName << "/" << optionSingleArgs["-o"] << SCAF_FOOTER << ">/dev/null 2>&1";
		system(rmOss.str().c_str());
	}
}


void IterateScaffold::execClusterMode(void)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " cluster"
        << " -t " << optionSingleArgs["-t"]
        << " -c " << optionSingleArgs["-o"] << ITERATION_FOOTER;
    if (optionSingleArgs["-n"] != "") {
        oss << " -n " << optionSingleArgs["-n"];
    }
    oss << " 2>" << optionSingleArgs["-o"] << ".clusterLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::ClusterError();
    }
}
