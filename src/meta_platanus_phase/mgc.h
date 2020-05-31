#ifndef MGC_H
#define MGC_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <float.h>
#include <omp.h>
#include "common.h"


class Mgc
{
private:
	struct MgcSeq
	{
		unsigned int nc;
		int cl;
		double cov;
		double score;
		std::vector<char> seq;
		std::vector<unsigned> kmer;
		std::vector<unsigned> pre_kmer;
		std::vector<double> frq;
		std::vector<unsigned> mono;
		std::vector<double> f_mono;
		std::string name;

		MgcSeq(): nc(0) {}
	};

	struct CL {
		unsigned len;
		double cov;
		std::vector<int> member;
		std::vector<double> kmer;
		std::vector<double> f_kmer;
		std::vector<double> di;
		std::vector<double> f_di;
		std::vector<int> pre_kmer;
		std::vector<int> ndi;
		std::vector<int> mono;
		std::vector<double> f_mono;

		CL () {
			mono.resize (64);
			di.resize (4096);
			f_mono.resize (64);
			f_di.resize (4096);
			ndi.resize (64);
		}
	};

	struct DistanceMatrix {
		std::vector<unsigned> id;
		std::vector<std::vector<double> > dist;
		std::vector<std::string> name;
		std::vector<double> d;
		std::vector<unsigned int> num;
	};

	struct UPGMAnode {
		double height;
		unsigned id;
		long child1;
		long child2;
	};


public:
	unsigned order;
	unsigned num_kmer;
	unsigned long totalContigLength;
    std::vector<std::string> name;
    std::vector<MgcSeq> seq;
	std::vector<CL> cl;
    DistanceMatrix mat;
	std::vector<UPGMAnode> UPGMAtree;
	platanus::empiricalDistribution empDistr;

	Mgc() {}
	Mgc(const int k): order(k - 1), num_kmer(1<<(k<<1)) {}
	double getDistance(const int i, const int j) { return mat.dist[i][j-i-1]; }

	void clear();
	void kmerStat();
	void dicodonStat();
	void setMatrix(const bool seqWeightFlag);
	void calcDistributionFromContig(std::vector<platanus::SEQ> &contig, const unsigned minLength, platanus::empiricalDistribution &empDistr);
	void setTotalContigLength(const platanus::Contig &contig, const unsigned minCDSLen);
	void setContig(const platanus::Contig &contig, const unsigned minLen);
	void setContigCDS(const platanus::Contig &contig, const unsigned minLen);
	void usage (std::vector<CL> &cl);
	void kmerUsage (std::vector<CL> &cl);
	double cl_score (std::vector<CL> &cl);
	double cl_score (CL &cl);
	double cl_kmer_score (std::vector<CL> &cl);
	double cl_kmer_score (CL &cl);
	std::vector<int> pivot (CL &cl, std::vector<double> &fm, double fmav);
	std::vector<int> kmerPivot (CL &cl, std::vector<double> &fm, double fmav);
	void init_fmap (std::vector<CL> &wcl, CL &cl, unsigned pn);
	void kmer_init_fmap (std::vector<CL> &wcl, CL &cl, unsigned pn);
	void init_rand (std::vector<CL> &wcl, CL &cl);
	void kmer_init_rand (std::vector<CL> &wcl, CL &cl);
	int recalc (CL &par, std::vector<CL> &cl);
	int kmerRecalc (CL &par, std::vector<CL> &cl);
	double calcCoverageScore();
	void kMeansClustering(unsigned long num_cluster, const unsigned maxitr);
	void kmerKMeansClustering(unsigned long num_cluster, const unsigned maxitr);
	void printClusters(const std::string &outFile);
	void initUPGMAtree();
	void UPGMA();
	long UPGMAmergePair();
	std::string rename (const std::string &name, const double d);
	void clusterUPGMAnodes(const double probThreshold);
	void UPGMAleafNodeIDs(long roodID, std::vector<long> &ret);
	void printNewick(const std::string &outFile);

	template <typename T>
	void addSeq(const T &inSeq)
	{
		seq.resize(seq.size() + 1);
		seq.rbegin()->seq.resize(inSeq.size());
		for (unsigned i = 0; i < inSeq.size(); ++i)
			seq.rbegin()->seq[i] = inSeq[i];
	}

	template <typename T>
	double calcPairDistance(const T &seq1, const T &seq2)
	{
		this->clear();
		this->addSeq(seq1);
		this->addSeq(seq2);
		this->kmerStat();
		this->setMatrix(false);
		return this->getDistance(0, 1);
	}
};


#endif
