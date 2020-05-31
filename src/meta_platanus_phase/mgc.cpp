#include "mgc.h"
#include <sstream>
#include <stack>

using std::vector;
using std::string;
using std::cerr;
using std::endl;


void Mgc::kmerStat()
{
    unsigned i, j, k, dof, mask, dame;
	unsigned key = 0;
	vector<int> tbl;

    mask = num_kmer - 1;

    tbl.assign (num_kmer, -1);
    for (i = 0, dof = 0; i < num_kmer; ++i) {
        if (tbl[i] != -1) continue;
        tbl[i] = dof;
        for (k = key = 0; k <= order; ++k) key = (key<<2)|((i>>k*2)&3);
		key ^= (0xffffffff&mask);
        tbl[key] = dof;
        ++dof;
    }

	#pragma omp parallel for private (j, key, dame)
    for (i = 0; i < seq.size (); ++i) {
		key = dame = 0;
//        seq[i].kmer.assign (dof, 1);
        seq[i].kmer.resize (num_kmer, 1);
		seq[i].nc = dof;

		for (j = 0; j < order; ++j) {
			key <<= 2;
			key |= seq[i].seq[j]&3;
			dame <<= 2;
			dame |= seq[i].seq[j]>>2;
		}
		for (; j < seq[i].seq.size (); ++j) {
			key <<= 2;
			key |= seq[i].seq[j]&3;
			dame <<= 2;
			dame |= seq[i].seq[j]>>2;
			dame &= mask;
			if (dame) continue;
			key &= mask;
			++seq[i].kmer[tbl[key]];
			++seq[i].nc;
		}

		seq[i].frq.assign (num_kmer, 1);
        seq[i].pre_kmer.assign (num_kmer >> 2, 4);

//		vector<double> pre (num_kmer/4, 4);

		for (j = 0; j < num_kmer; ++j) {
			seq[i].frq[j] = seq[i].kmer[tbl[j]];
			seq[i].pre_kmer[j>>2] += seq[i].kmer[tbl[j]];
		}
		for (j = 0; j < num_kmer; ++j)
			seq[i].frq[j] /= seq[i].pre_kmer[j>>2];
    }
}

void Mgc::dicodonStat()
{
    unsigned i, j, k, key, key2;

	#pragma omp parallel for private (j, k, key, key2)
    for (i = 0; i < seq.size (); i++) {
		if (seq[i].seq[1] == 3 && seq[i].seq[2] == 2) k = 3;
		else k = 0;

        seq[i].kmer.assign(4096, 0);
		seq[i].pre_kmer.assign(64, 0);
        seq[i].frq.resize(4096);
		seq[i].mono.assign(64, 0);
		seq[i].f_mono.resize(64);
		seq[i].nc = 0;

/*
		for (j = 0; j < 64; j++) {
			if (j == 48 || j == 50 || j == 56 || j >= 64 ) continue;
			seq[i].mono[j]++;
			seq[i].nc++;
			for (k = 0; k < 64; k++) {
				if (k == 48 || k == 50 || k == 56 || k >= 64 ) continue;
				seq[i].pre_kmer[j]++;
				key = (j<<6)|k;
				seq[i].kmer[key]++;
			}
		}
*/

		for (j = k; j < seq[i].seq.size ()-6; j += 3) {
			key = (seq[i].seq[j]<<4)|(seq[i].seq[j+1]<<2)|seq[i].seq[j+2];
			if (key == 48 || key == 50 || key == 56 || key >= 64 ) continue;
			seq[i].mono[key]++;
			seq[i].nc++;

			key2 = (seq[i].seq[j+3]<<4)|(seq[i].seq[j+4]<<2)|seq[i].seq[j+5];
			if (key2 == 48 || key2 == 50 || key2 == 56 || key2 >= 64 ) continue;

			seq[i].pre_kmer[key]++;
			key = (key<<6)|key2;
			seq[i].kmer[key]++;
		}
		key = (seq[i].seq[j]<<4)|(seq[i].seq[j+1]<<2)|seq[i].seq[j+2];
		if (!(key == 48 || key == 50 || key == 56 || key >= 64)) {
			seq[i].mono[key]++;
			seq[i].nc++;
		}

		for (j = 0; j < 64; j++) {
			seq[i].f_mono[j] = seq[i].nc > 0 ? (double)seq[i].mono[j]/seq[i].nc : 0.0;
			if (seq[i].pre_kmer[j] == 0) continue;
			for (k = 0; k < 64; k++) {
				key = (j<<6)|k;
				seq[i].frq[key] = (double)seq[i].kmer[key]/seq[i].pre_kmer[j];
			}
		}
    }
}




void Mgc::clear()
{
	mat.dist.clear();
//	name.clear();
	mat.d.clear();
	mat.num.clear();
	seq.clear();
}


void Mgc::setMatrix(const bool seqWeightFlag)
{
    unsigned i, j, k;
    double f, w;

    mat.id.resize(seq.size());
    mat.dist.resize(seq.size());

	#pragma omp parallel for private (j, k, f, w) schedule (dynamic)
    for (i = 0; i < seq.size (); ++i) {
		mat.id[i] = i;
		mat.dist[i].resize (seq.size ()-i-1);
		for (j = i+1; j < seq.size (); ++j) {
			f = 0.0;
			for (k = 0; k < seq[i].frq.size (); ++k) {
				w = seq[i].frq[k]-seq[j].frq[k];
				f += w*w;
			}
			mat.dist[i][j-i-1] = sqrt (f);
		}
    }

    mat.d.assign (seq.size (), 0.0);
    if (seqWeightFlag) {
		mat.num.assign (seq.size (), 1);
    }else {
		mat.num.resize (seq.size ());
		#pragma omp parallel for
		for (i = 0; i < seq.size (); ++i) mat.num[i] = seq[i].nc;
    }

//    vector<MgcSeq> ().swap(seq);
}

void Mgc::initUPGMAtree()
{
	UPGMAtree.resize(seq.size());

    for (unsigned i = 0; i < seq.size (); ++i) {
		UPGMAtree[i].id = i;
		UPGMAtree[i].height = 0.0;
		UPGMAtree[i].child1 = -1;
		UPGMAtree[i].child2 = -1;
	}
}

void Mgc::UPGMA()
{
    mat.name.resize(seq.size());
    for (unsigned i = 0; i < seq.size (); ++i)
		mat.name[i] = name[i];

	while (UPGMAmergePair() > 1);
}


long Mgc::UPGMAmergePair()
{
    DistanceMatrix work;
    long i, j, k, mini, minj, wi, wj;
    double sc, min; 

    min = DBL_MAX;
	mini = minj = 0;
    for (i = 0; static_cast<unsigned>(i) < mat.name.size ()-1; ++i) {
		for (j = i+1; static_cast<unsigned>(j) < mat.name.size (); ++j) {
			wj = j-i-1;
			if (mat.dist[i][wj] < min) {
				min = mat.dist[i][wj];
				mini = i;
				minj = j;
			}
		}
    }

    sc = mat.dist[mini][minj-mini-1]/2.0;
    mat.name[mini] = rename(mat.name[mini], sc-mat.d[mini]);
    mat.name[minj] = rename(mat.name[minj], sc-mat.d[minj]);

	UPGMAnode wnode;
	wnode.height = sc;
	wnode.child1 = mat.id[mini];
	wnode.child2 = mat.id[minj];
	UPGMAtree.push_back(wnode);

    work.name.resize (mat.name.size ()-1);
    work.dist.resize (work.name.size ());
    work.d.resize (work.name.size ());
    work.num.resize (work.name.size ());
    work.id.resize (work.name.size ());

    work.name[0] += "("+mat.name[mini]+","+mat.name[minj]+")";
    work.d[0] = sc;
    work.num[0] = mat.num[mini]+mat.num[minj];

    work.dist[0].resize (mat.name.size ()-2);
	work.id[0] = UPGMAtree.size() - 1;
    for (i = k = 0; static_cast<unsigned>(i) < mat.name.size (); ++i) {
		if (i == mini || i == minj) continue;
		if (i < mini) {
			wi = i;
			wj = mini-i-1;
		}else {
			wi = mini;
			wj = i-mini-1;
		}
		work.dist[0][k] = mat.num[mini]*mat.dist[wi][wj];

		if (i < minj) {
			wi = i;
			wj = minj-i-1;
		}else {
			wi = minj;
			wj = i-minj-1;
		}
		work.dist[0][k] += mat.num[minj]*mat.dist[wi][wj];
		work.dist[0][k] /= mat.num[mini]+mat.num[minj];
		++k;
    }
    
    for (i = 0, wi = 1; static_cast<unsigned>(i) < mat.name.size (); ++i) {
		if (i == mini || i == minj) continue;
		work.id[wi] = mat.id[i];
		work.name[wi] = mat.name[i];
		work.d[wi] = mat.d[i];
		work.num[wi] = mat.num[i];
		work.dist[wi].resize (mat.name.size ()-2-wi);
		for (j = i+1, wj = 0; static_cast<unsigned>(j) < mat.name.size (); ++j) {
			if (j == mini || j == minj) continue;
			work.dist[wi][wj++] = mat.dist[i][j-i-1];
		}
		++wi;
    }
    
    mat = work;

    return mat.name.size ();
}


string Mgc::rename (const string &name, const double d)
{
    std::ostringstream os;
    os << name << ":" << d;
    return os.str();
}


void Mgc::clusterUPGMAnodes(const double probThreshold)
{
	long numClusters = 0;
	cl.clear();

	for (auto nodeIt = UPGMAtree.begin(); nodeIt != UPGMAtree.end(); ++nodeIt) {
		if (nodeIt->child1 >= 0 && empDistr.probGreater(2.0 * nodeIt->height) < probThreshold) {
			if (empDistr.probGreater(2.0 * UPGMAtree[nodeIt->child1].height) >= probThreshold) {
				cl.resize(numClusters + 1);
				vector<long> leafNodeIDs;
				UPGMAleafNodeIDs(nodeIt->child1, leafNodeIDs);
				cl[numClusters].member.clear();
				for (auto leafIt = leafNodeIDs.begin(); leafIt != leafNodeIDs.end(); ++leafIt)
					cl[numClusters].member.push_back(UPGMAtree[*leafIt].id);
				++numClusters;
			}

			if (empDistr.probGreater(2.0 * UPGMAtree[nodeIt->child2].height) >= probThreshold) {
				cl.resize(numClusters + 1);
				vector<long> leafNodeIDs;
				UPGMAleafNodeIDs(nodeIt->child2, leafNodeIDs);
				cl[numClusters].member.clear();
				for (auto leafIt = leafNodeIDs.begin(); leafIt != leafNodeIDs.end(); ++leafIt)
					cl[numClusters].member.push_back(UPGMAtree[*leafIt].id);
				++numClusters;
			}
		}
	}

	if (numClusters == 0) {
		if (UPGMAtree.size() > 1) {
			double max= 0.0;
			unsigned maxi = 0;
			for (unsigned i = 0; i <  UPGMAtree.size(); ++i) {
				if (UPGMAtree[i].height > max) {
					max = UPGMAtree[i].height;
					maxi = i;
				}
			}
			vector<long> leafNodeIDs;
			UPGMAleafNodeIDs(maxi, leafNodeIDs);
			cl[numClusters].member.clear();
			for (auto leafIt = leafNodeIDs.begin(); leafIt != leafNodeIDs.end(); ++leafIt)
				cl[numClusters].member.push_back(UPGMAtree[*leafIt].id);
		}
		else {
			cl.resize(1);
			cl[0].member.clear();
			cl[0].member.push_back(0);
		}
	}
}

void Mgc::UPGMAleafNodeIDs(long rootID, vector<long> &ret)
{
	ret.clear();

	std::stack<long> buf;
	buf.push(rootID);

	while (!buf.empty()) {
		long id = buf.top();
		buf.pop();
		if (UPGMAtree[id].child1 >= 0) {
			buf.push(UPGMAtree[id].child1);
			buf.push(UPGMAtree[id].child2);
		}
		else {
			ret.push_back(id);
		}
	}
}
	
void Mgc::printNewick(const string &outFile)
{
    std::ofstream ofs(outFile.c_str());
    if (!ofs) throw platanus::FILEError(outFile);

//	ofs << "(" << mat.name[0] 
//	     << "," << mat.name[1] << ");" << endl;
	ofs << mat.name[0] << endl;

	ofs.close();
}


void Mgc::calcDistributionFromContig(vector<platanus::SEQ> &contig, const unsigned minLength, platanus::empiricalDistribution &empDistr)
{
	std::cerr << "calculating the distribution of kmer-frequency-distance..." << std::endl;

	for (auto seqIt = contig.begin(); seqIt != contig.end(); ++seqIt) {
		unsigned leftLen = seqIt->base.size() / 2;
		unsigned rightLen = seqIt->base.size() - leftLen;

		if (leftLen < minLength) continue;

		this->addSeq(seqIt->base.substr(0, leftLen));
		this->addSeq(seqIt->base.substr(leftLen, rightLen));
		empDistr.count(calcPairDistance(seqIt->base.substr(0, leftLen), seqIt->base.substr(leftLen, rightLen)));
	}
	empDistr.setDistribution();
}

void Mgc::setTotalContigLength(const platanus::Contig &contig, const unsigned minCDSLen)
{
	totalContigLength = 0;
	for (auto it = contig.seq.begin(); it != contig.seq.end(); ++it) {
		if (it->CDS.size() >= minCDSLen) {
			totalContigLength += it->base.size();
		}
	}
}


void Mgc::setContig(const platanus::Contig &contig, const unsigned minLen)
{
	name.clear();

	for (auto it = contig.seq.begin(); it != contig.seq.end(); ++it) {
		if (platanus::getNonNLength(it->base) >= minLen) {
			name.push_back(contig.name[it - contig.seq.begin()]);
			addSeq(it->base);
		}
	}
}
		

void Mgc::setContigCDS(const platanus::Contig &contig, const unsigned minLen)
{
	name.clear();

	for (auto it = contig.seq.begin(); it != contig.seq.end(); ++it) {

//		if (platanus::getNonNLength(it->CDS) >= minLen) {
		if (it->CDS.size() >= minLen) {
			name.push_back(contig.name[it - contig.seq.begin()]);
			seq.resize(seq.size() + 1);
			seq.rbegin()->seq.resize(it->CDS.size());
			seq.rbegin()->cov = static_cast<double>(contig.coverage[it - contig.seq.begin()]);
			for (unsigned i = 0; i < it->CDS.size(); ++i)
				seq.rbegin()->seq[i] = it->CDS[i] != 4 ? it->CDS[i] : 64;
		}
	}
}


void Mgc::usage (vector<CL> &cl)
{
    unsigned i, j, k, l, key;
    unsigned long len;

#pragma omp parallel for private (j, k, l, key, len)
    for (i = 0; i < cl.size (); i++) {
		cl[i].cov = 0.0;
		len = 0;
		for (k = 0; k < 64; k++) {
			cl[i].ndi[k] = 0;
			for (l = 0; l < 64; l++) cl[i].di[((k<<6)|l)] = 0;
		}
		for (j = 0; j < cl[i].member.size (); j++) {
			cl[i].cov += static_cast<double>(seq[cl[i].member[j]].cov * seq[cl[i].member[j]].seq.size());
			len += seq[cl[i].member[j]].seq.size();
			for (k = 0; k < 64; k++) {
				for (l = 0; l < 64; l++) {
					key = (k<<6)|l;
					cl[i].di[key] += seq[cl[i].member[j]].kmer[key];
					cl[i].ndi[k] += seq[cl[i].member[j]].kmer[key];
				}
			}
		}
		cl[i].len = len;
		cl[i].cov = cl[i].cov / len;

		for (k = 0; k < 64; k++) {
			for (l = 0; l < 64; l++) {
				key = (k<<6)|l;
				if (cl[i].ndi[k] == 0) {
					cl[i].di[key] = -10;
					cl[i].f_di[key] = 0;
				}else {
					if (cl[i].di[key] == 0) cl[i].di[key] = 0.01;
					cl[i].f_di[key] = cl[i].di[key]/cl[i].ndi[k];
					cl[i].di[key] = log2 (cl[i].f_di[key]*61);
				}
			}
		}
    }
}


void Mgc::kmerUsage (vector<CL> &cl)
{
    unsigned i, j, k;
    unsigned long len;

#pragma omp parallel for private (j, k, len)
    for (i = 0; i < cl.size (); i++) {
		cl[i].cov = 0.0;
		len = 0;
		cl[i].kmer.assign(num_kmer, 0.0);
		cl[i].pre_kmer.assign(num_kmer>>2, 0);
		cl[i].f_kmer.assign(num_kmer, 0.0);

		for (j = 0; j < cl[i].member.size (); j++) {
			cl[i].cov += static_cast<double>(seq[cl[i].member[j]].cov * seq[cl[i].member[j]].seq.size());
			len += seq[cl[i].member[j]].seq.size();

			for (k = 0; k < num_kmer; ++k) {
				cl[i].kmer[k] += seq[cl[i].member[j]].kmer[k];
				cl[i].pre_kmer[k>>2] += seq[cl[i].member[j]].kmer[k];
			}
		}
		cl[i].len = len;
		cl[i].cov = cl[i].cov / len;

		for (k = 0; k < num_kmer; ++k) {
			if (cl[i].pre_kmer[k>>2] == 0) {
				cl[i].kmer[k] = -10;
				cl[i].f_kmer[k] = 0;
			}else {
				if (cl[i].kmer[k] == 0) cl[i].kmer[k] = 1.0 / cl[i].len;
				cl[i].f_kmer[k] = cl[i].kmer[k] / cl[i].pre_kmer[k>>2];
//				cl[i].kmer[k] = log2(cl[i].f_kmer[k]*num_kmer);
				cl[i].kmer[k] = log2(cl[i].f_kmer[k]);
			}
		}
    }
}


double Mgc::cl_score (vector<CL> &cl)
{
    unsigned i, j, k;
    double sc;
    
    sc = 0;
    for (i = 0; i < cl.size (); i++) {
		for (j = 0; j < cl[i].member.size (); j++) {
			for (k = 0; k < 4096; k++)
				sc += cl[i].di[k]*seq[cl[i].member[j]].kmer[k];
		}
	}
    return sc;
}


double Mgc::cl_kmer_score (vector<CL> &cl)
{
    unsigned i, j, k;
    double sc;
    
    sc = 0;
    for (i = 0; i < cl.size (); i++) {
		for (j = 0; j < cl[i].member.size (); j++) {
			for (k = 0; k < num_kmer; k++)
				sc += cl[i].kmer[k]*seq[cl[i].member[j]].kmer[k];
		}
	}
    return sc;
}


double Mgc::cl_score (CL &cl)
{
    unsigned j, k;
    double sc;
    
    sc = 0;
    for (j = 0; j < cl.member.size (); j++) {
		for (k = 0; k < 4096; k++)
			sc += cl.di[k]*seq[cl.member[j]].kmer[k];
    }
    return sc;
}


double Mgc::cl_kmer_score (CL &cl)
{
    unsigned j, k;
    double sc;
    
    sc = 0;
    for (j = 0; j < cl.member.size (); j++) {
		for (k = 0; k < num_kmer; k++) {
			sc += cl.kmer[k]*seq[cl.member[j]].kmer[k];
		}
    }
    return sc;
}


vector<int> Mgc::pivot (CL &cl, vector<double> &fm, double fmav)
{
    unsigned i, j, k, wp;
    double d, wd, max;
    vector<int> p, cln;

    p.resize (2); cln.resize (2);

    p[0] = cl.member[0];
    cln[0] = 0;
    for (i = 0; i < 6; i++) {
		wp = i%2;
		max = 0;
		for (j = 0; j < cl.member.size (); j++) {
			d = 0;
			for (k = 0; k < 64; k++) {
				wd = seq[cl.member[j]].f_mono[k]-seq[p[wp]].f_mono[k];
				d += wd*wd;
				//if (cl.f_mono[k] > 0) wmono = cl.f_mono[k];
				//else wmono = 0.0001;
				//d += wd*wd/(wmono*wmono);
			}
			if (fm[cln[wp]]-fm[j] != 0) {
				wd = fm[cln[wp]]-fm[j];
				d -= wd*wd;
				//if (fmav != 0) wmono = fmav;
				//else wmono = 0.0001;
				//if (wmono > 0) d -= wd*wd/(wmono*wmono);
			}

			if (d > max) {
				max = d;
				p[wp^1] = cl.member[j];
				cln[wp^1] = j;
			}
		}
    }
    return p;
}


vector<int> Mgc::kmerPivot (CL &cl, vector<double> &fm, double fmav)
{
    unsigned i, j, k, wp;
    double d, wd, max;
    vector<int> p, cln;

    p.resize (2); cln.resize (2);

    p[0] = cl.member[0];
    cln[0] = 0;
    for (i = 0; i < 6; i++) {
		wp = i%2;
		max = 0;
		for (j = 0; j < cl.member.size (); j++) {
			d = 0;
			for (k = 0; k < num_kmer; k++) {
				wd = seq[cl.member[j]].frq[k] - seq[p[wp]].frq[k];
				d += wd*wd;
			}
			if (fm[cln[wp]]-fm[j] != 0) {
				wd = fm[cln[wp]]-fm[j];
				d -= wd*wd;
			}

			if (d > max) {
				max = d;
				p[wp^1] = cl.member[j];
				cln[wp^1] = j;
			}
		}
    }
    return p;
}


void Mgc::init_fmap (vector<CL> &wcl, CL &cl, unsigned pn)
{
    unsigned i, j, k, dim;
    vector<double> fm (cl.member.size ());
    double d, wd, max;
    vector<int> p (2);
    
    pn++; wd = 1;
    for (dim = 0; dim < pn; dim++) {
	p = pivot (cl, fm, wd);

	d = 0;
	for (k = 0; k < 64; k++) {
	    wd = seq[p[0]].f_mono[k]-seq[p[1]].f_mono[k];
	    d += wd*wd;
	}
	max = 2*sqrt (d);

	wd = 0;
	for (i = 0; i < cl.member.size (); i++) {
	    fm[i] = 0;
	    for (k = 0; k < 64; k++) {
			fm[i] += (seq[p[0]].f_mono[k]-seq[cl.member[i]].f_mono[k])*(seq[p[0]].f_mono[k]-seq[cl.member[i]].f_mono[k])-(seq[p[1]].f_mono[k]-seq[cl.member[i]].f_mono[k])*(seq[p[1]].f_mono[k]-seq[cl.member[i]].f_mono[k]);
	    }
	    fm[i] = (fm[i]+d)/max;
	    wd += fm[i];
	}
	wd /= cl.member.size ();
    }

    for (i = 0; i < wcl.size (); i++) wcl[i].member.clear ();    
    for (i = 0; i < cl.member.size (); i++) {
	if (fm[i] < wd) j = 0;
	else j = 1;
	seq[cl.member[i]].cl = j;
	wcl[j].member.push_back (cl.member[i]);
    }

    usage (wcl);
}


void Mgc::kmer_init_fmap (vector<CL> &wcl, CL &cl, unsigned pn)
{
    unsigned i, j, k, dim;
    vector<double> fm (cl.member.size ());
    double d, wd, max;
    vector<int> p (2);
    
    pn++; wd = 1;
    for (dim = 0; dim < pn; dim++) {
	p = kmerPivot (cl, fm, wd);

	d = 0;
	for (k = 0; k < num_kmer; k++) {
	    wd = seq[p[0]].frq[k]-seq[p[1]].frq[k];
	    d += wd*wd;
	}

	max = 2*sqrt (d);
	wd = 0;
	for (i = 0; i < cl.member.size (); i++) {
	    fm[i] = 0;
	    for (k = 0; k < num_kmer; k++) {
			fm[i] += (seq[p[0]].frq[k]-seq[cl.member[i]].frq[k])*(seq[p[0]].frq[k]-seq[cl.member[i]].frq[k]) - (seq[p[1]].frq[k]-seq[cl.member[i]].frq[k])*(seq[p[1]].frq[k]-seq[cl.member[i]].frq[k]);
	    }
	    fm[i] = (fm[i]+d)/max;
	    wd += fm[i];
	}
	wd /= cl.member.size ();
    }

    for (i = 0; i < wcl.size (); i++) wcl[i].member.clear ();    
    for (i = 0; i < cl.member.size (); i++) {
	if (fm[i] < wd) j = 0;
	else j = 1;
	seq[cl.member[i]].cl = j;
	wcl[j].member.push_back (cl.member[i]);
    }

    kmerUsage (wcl);
}


void Mgc::init_rand (vector<CL> &wcl, CL &cl)
{
    unsigned i, j;

    for (i = 0; i < wcl.size (); i++) wcl[i].member.clear ();
    
    srand ((unsigned)time (NULL));
//    srand (1);
    for (i = 0; i < cl.member.size (); i++) {
		j = rand ()%wcl.size ();
		seq[cl.member[i]].cl = j;
		wcl[j].member.push_back (cl.member[i]);
    }

    usage (wcl);
}


void Mgc::kmer_init_rand (vector<CL> &wcl, CL &cl)
{
    unsigned i, j;

    for (i = 0; i < wcl.size (); i++) wcl[i].member.clear ();
    
//    srand ((unsigned)time (NULL));
    srand (1);
    for (i = 0; i < cl.member.size (); i++) {
		j = rand ()%wcl.size ();
		seq[cl.member[i]].cl = j;
		wcl[j].member.push_back (cl.member[i]);
    }

    kmerUsage (wcl);
}



int Mgc::recalc (CL &par, vector<CL> &cl)
{
    int maxcl, miss;
    unsigned i, j, k;
    double score, maxsc;

    for (j = 0; j < cl.size (); j++) cl[j].member.clear ();
    miss = 0;

#pragma omp parallel for private (j, k, score, maxcl, maxsc) reduction (+: miss)
    for (i = 0; i < par.member.size (); i++) {
		maxsc = -1e-30;
		maxcl = 0;
		//maxsc = 999;
		for (j = 0; j < cl.size (); j++) {
			score = 0;
			for (k = 0; k < 4096; k++) {
				score += cl[j].di[k]*seq[par.member[i]].kmer[k];
				//score += (cl[j].f_di[k]-seq[par.member[i]].frq[k])*(cl[j].f_di[k]-seq[par.member[i]].frq[k]);
			}
			score /= seq[par.member[i]].nc;
			if (maxsc < score) {
			//if (maxsc > score) {
				maxsc = score;
				maxcl = j;
			}
		}

		if (maxcl != seq[par.member[i]].cl) {
			miss++;
			seq[par.member[i]].cl = maxcl;
		}
		seq[par.member[i]].score = maxsc;
		//cl[maxcl].member.push_back (par.member[i]);
    }
    for (i = 0; i < par.member.size (); i++)
		cl[seq[par.member[i]].cl].member.push_back (par.member[i]);
    
    usage (cl);

    return miss;
}


int Mgc::kmerRecalc (CL &par, vector<CL> &cl)
{
    int maxcl, miss;
    unsigned i, j, k;
    double score, maxsc;

    for (j = 0; j < cl.size (); j++) cl[j].member.clear ();
    miss = 0;

#pragma omp parallel for private (j, k, score, maxcl, maxsc) reduction (+: miss)
    for (i = 0; i < par.member.size (); i++) {
//		maxsc = -1e-30;
		maxsc = -DBL_MAX;
		maxcl = 0;
		//maxsc = 999;
		for (j = 0; j < cl.size (); j++) {
			score = 0;
			for (k = 0; k < num_kmer; k++) {
				score += cl[j].kmer[k]*seq[par.member[i]].kmer[k];
			}
			score /= seq[par.member[i]].nc;
			if (maxsc < score) {
				maxsc = score;
				maxcl = j;
			}
		}

		if (maxcl != seq[par.member[i]].cl) {
			miss++;
			seq[par.member[i]].cl = maxcl;
		}
		seq[par.member[i]].score = maxsc;
		//cl[maxcl].member.push_back (par.member[i]);
    }
    for (i = 0; i < par.member.size (); i++)
		cl[seq[par.member[i]].cl].member.push_back (par.member[i]);
    
    kmerUsage (cl);

    return miss;
}


double Mgc::calcCoverageScore()
{
    unsigned i, j, k;
	unsigned long numInconsistent;
    double selfDist;

	numInconsistent = 0;
    for (i = 0; i < cl.size(); i++) {
		for (j = 0; j < cl[i].member.size(); ++j) {
			selfDist = std::abs(cl[i].cov - seq[cl[i].member[j]].cov);
			for (k = 0; k < cl.size(); ++k) {
				if (selfDist > std::abs(cl[k].cov - seq[cl[i].member[j]].cov)) {
//					numInconsistent += seq[cl[i].member[j]].seq.size();
					++numInconsistent;
					break;
				}
			}
		}
	}
    return -static_cast<double>(numInconsistent);
}


void Mgc::kMeansClustering(unsigned long num_cluster, const unsigned maxitr)
{
	const unsigned long EXPECTED_BACTERIAL_GENOME_SIZE = 3000000;

    vector<CL> maxcl, wcl, max_cl;
    CL par;
    unsigned i, j, k, size, miss, maxsize;
    int maxi = 0;
    double clsc, maxsc, check;

	maxsize = num_cluster;
	if (maxsize == 0)
		maxsize = std::max(totalContigLength /  EXPECTED_BACTERIAL_GENOME_SIZE, 2ul);

	cerr << "\nNumber of clusters: " << maxsize << endl;

	cl.clear();
    for (i = 0; i < seq.size (); i++) par.member.push_back (i);
    cl.push_back (par);
    usage (cl);
    check = cl_score (cl);

    wcl.resize (2);
    for (size = 2; size <= maxsize; size++) {
		maxsc = 0;
		for (i = 0; i < cl.size (); i++) {
			for (k = 0; k < 3; k++) {
				init_fmap (wcl, cl[i], k);
//				init_rand(wcl, cl[i]);

				j = 0;
				while (1) {
					miss = recalc (cl[i], wcl);
					j++;
					if (miss == 0 || j == maxitr) break;
				}
				clsc = cl_score(wcl) - cl_score(cl[i]);
				if (maxsc < clsc) {
					maxsc = clsc;
					maxcl = wcl;
					maxi = i;
				}
			}
		}
		if (maxsc == 0) break;

		cl[maxi] = maxcl[0];
		cl.push_back (maxcl[1]);
		for (i = 0; i < cl.size (); i++) {
			for (j = 0; j < cl[i].member.size (); j++)
				seq[cl[i].member[j]].cl = i;
		}

		j = 0;
		while (1) {
			miss = recalc (par, cl);
			j++;
			if (miss == 0 || j == maxitr) break;
		}
		clsc = cl_score(cl);
		check = clsc;

    }
}


void Mgc::kmerKMeansClustering(unsigned long num_cluster, const unsigned maxitr)
{
	const unsigned long EXPECTED_BACTERIAL_GENOME_SIZE = 3000000;

    vector<CL> maxcl, wcl, max_cl;
    CL par;
    unsigned i, j, k, size, miss, maxsize;
    int maxi = 0;
    double clsc, maxsc, check;

	maxsize = num_cluster;
	if (maxsize == 0)
		maxsize = std::max(totalContigLength /  EXPECTED_BACTERIAL_GENOME_SIZE, 2ul);

	cerr << "\nNumber of clusters: " << maxsize << endl;

	cl.clear();
    for (i = 0; i < seq.size (); i++) par.member.push_back (i);
    cl.push_back (par);
    kmerUsage (cl);
    check = cl_kmer_score (cl);

    wcl.resize (2);
    for (size = 2; size <= maxsize; size++) {
		cerr << "current cluster size: " << size << endl;
		maxsc = 0;
		for (i = 0; i < cl.size (); i++) {
			for (k = 0; k < 3; k++) {
				kmer_init_fmap (wcl, cl[i], k);
//				kmer_init_rand(wcl, cl[i]);

				j = 0;
				while (1) {
					miss = kmerRecalc (cl[i], wcl);
					j++;
					if (miss == 0 || j == maxitr) break;
				}
				clsc = cl_kmer_score (wcl)-cl_kmer_score (cl[i]);
//				clsc = cl[i].len * (cl_kmer_score(wcl) - cl_kmer_score(cl[i]));
				if (maxsc < clsc) {
					maxsc = clsc;
					maxcl = wcl;
					maxi = i;
				}
			}
		}
		if (maxsc == 0) break;

		cl[maxi] = maxcl[0];
		cl.push_back (maxcl[1]);
		for (i = 0; i < cl.size (); i++) {
			for (j = 0; j < cl[i].member.size (); j++)
				seq[cl[i].member[j]].cl = i;
		}

		j = 0;
		while (1) {
			miss = kmerRecalc (par, cl);
			j++;
			if (miss == 0 || j == maxitr) break;
		}
		clsc = cl_kmer_score (cl);
		check = clsc;

    }
}


void Mgc::printClusters(const std::string &outFile)
{
    std::ofstream ofs(outFile.c_str());
    if (!ofs) throw platanus::FILEError(outFile);

    for (unsigned j = 0; j < cl.size (); j++) {
		for (unsigned k = 0; k < cl[j].member.size (); k++)
			ofs << name[cl[j].member[k]] <<'\t' << j + 1 << endl;
	}

	ofs.close();
}
