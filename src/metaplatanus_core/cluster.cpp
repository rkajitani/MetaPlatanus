#include "cluster.h"


using std::cerr;
using std::endl;
using std::string;


const unsigned Cluster::KMER_LENGTH_FOR_CLUSTERING = 4;

//////////////////////////////////////////////////////////////////////////////////////
// default constructor
//////////////////////////////////////////////////////////////////////////////////////
Cluster::Cluster()
{
    optionSingleArgs["-c"] = "";
    optionSingleArgs["-n"] = "0";
    optionSingleArgs["-i"] = "100";
    optionSingleArgs["-l"] = "1200";
    optionSingleArgs["-t"] = "1";
    optionSingleArgs["-mga"] = "mga";
    optionBool["-hcl"] = false;
//    optionSingleArgs["-p"] = "0.1";
}


//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void Cluster::usage(void) const
{
    std::cerr << "\nUsage: metaplatanus_core cluster [Options]\n"
              << "Options:\n"
              << "    -c FILE     : contig or scaffold_file (fasta format)\n"
              << "    -i INT      : maximu iteration count (default " << optionSingleArgs.at("-i") << ")\n"
              << "    -l INT      : minimum length for clustering (non-N bases, default " << optionSingleArgs.at("-l") << ")\n"
              << "    -t INT      : number of threads (default " << optionSingleArgs.at("-t") << ")\n"
              << "    -hcl        : hierarchical clustering (default k-means)\n"
              << "    -n INT      : number of clusters for k-means clustering (default 0, 0 means auto)\n"
//              << "    -p FLOAT    : probability threthold for hierarchical clustering\n"
              << "    -mga FILE   : name of mga executable file (default " << optionSingleArgs.at("-mga") << ")\n"

              << "Output:\n"
              << "    [contig_file(-c)].mga.txt\n"
              << "    [contig_file(-c)].clusters.tsv\n"
              << "    (-hcl) [contig_file(-c)].nwk\n"
              << std::endl;
}


void Cluster::exec(void)
{
    omp_set_num_threads(atoi(optionSingleArgs["-t"].c_str()));

    platanus::Contig contig;
	contig.readFastaCoverage(optionSingleArgs["-c"]);
	contig.setNameIndex();

	Mgc mgc(KMER_LENGTH_FOR_CLUSTERING);

	if (!(optionBool["-hcl"])) {
		string mgaOutFilename(optionSingleArgs["-c"]);
		mgaOutFilename += ".mga.txt";
		contig.execMga(optionSingleArgs["-c"], mgaOutFilename, optionSingleArgs["-mga"]);
		contig.readMgaResult(mgaOutFilename);

		mgc.setContigCDS(contig, atoi(optionSingleArgs["-l"].c_str()));
		mgc.setTotalContigLength(contig, atoi(optionSingleArgs["-l"].c_str()));
		mgc.dicodonStat();
		mgc.kMeansClustering(atoi(optionSingleArgs["-n"].c_str()), atoi(optionSingleArgs["-i"].c_str()));
//		mgc.setContig(contig, atoi(optionSingleArgs["-l"].c_str()));
//		mgc.kmerStat();
//		mgc.kmerKMeansClustering(atoi(optionSingleArgs["-n"].c_str()), atoi(optionSingleArgs["-i"].c_str()));

		string clusterOutFilename(optionSingleArgs["-c"]);
		clusterOutFilename += ".clusters.tsv";
		mgc.printClusters(clusterOutFilename);
	}
	else {
		mgc.setContig(contig, atoi(optionSingleArgs["-l"].c_str()));
		mgc.setTotalContigLength(contig, 0);
		mgc.kmerStat();
		mgc.kmerKMeansClustering(atoi(optionSingleArgs["-n"].c_str()), atoi(optionSingleArgs["-i"].c_str()));
		mgc.setMatrix(false);
		mgc.initUPGMAtree();
		mgc.UPGMA();
		Mgc (KMER_LENGTH_FOR_CLUSTERING).calcDistributionFromContig(contig.seq, atoi(optionSingleArgs["-l"].c_str()), mgc.empDistr);
//		mgc.clusterUPGMAnodes(atof(optionSingleArgs["-p"].c_str()));

		string newickOutFilename(optionSingleArgs["-c"]);
		newickOutFilename += ".nwk";
		mgc.printNewick(newickOutFilename);
	}

    cerr << "clustering completed!" << endl;
}
