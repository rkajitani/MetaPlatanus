#include "cov_trim.h"




//////////////////////////////////////////////////////////////////////////////////////
// default constructor
//////////////////////////////////////////////////////////////////////////////////////
ContigDivider::ContigDivider(): contig(), readLength(100), contigMaxK(100), coverageRatioThreshold(), contigItself(), breakPoint()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-k"] = "";
    optionSingleArgs["-r"] = "2.0";
    optionMultiArgs["-c"] = std::vector<std::string>();
    optionSingleArgs["-tmp"] = ".";
    optionBool["-recalc_cov"] = false;
}


//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void ContigDivider::usage(void) const
{
    std::cerr << "\nUsage: metaplatanus_core cov_trim [Options]\n"
              << "Options:\n"
              << "    -o STR               : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...] : contig or scaffold_file (fasta format)\n"
              << "    -k STR               : k-mer-occurrence file (binary format)\n"
              << "    -r FLOAT             : coverage ratio threshold (default " << optionSingleArgs.at("-r") << ")\n"
              << "    -recalc_cov          : just re-calculate coverage depth and not trim sequences\n"
              << "    -tmp DIR             : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n\n\n"

              << "Output:\n"
              << "    PREFIX_cov_trimmed.fa (when \"-recalc_cov\" not specified)\n"
              << "    PREFIX_recalc.fa (when \"-recalc_cov\" specified)\n"
              << std::endl;
}




void ContigDivider::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
    unsigned long long kmerLength = platanus::getKmerLengthFromBinary(optionSingleArgs["-k"]);

    this->readLength = platanus::Contig::getReadLengthFromFastaHeader(optionMultiArgs["-c"][0]);
	this->contigMaxK = contig.getMaxKFromFastaHeader(optionMultiArgs["-c"][0]);

	if (!optionBool["-recalc_cov"]) {
		this->coverageRatioThreshold = atof(optionSingleArgs["-r"].c_str());
	}
	else {
		this->coverageRatioThreshold = DBL_MAX;
	}

    for (auto itr = optionMultiArgs["-c"].begin(), end = optionMultiArgs["-c"].end(); itr != end; ++itr) {
        contig.readFastaCoverage(*itr);
    }
	contig.setNameIndex();

    if (kmerLength <= 32) {
        Counter<Kmer31> counter(kmerLength);
        counter.readOccurrenceTableBinary(optionSingleArgs["-k"]);
        this->getOccurrenceArray(counter);
    } else if (kmerLength <= 64) {
        Counter<KmerN<Binstr63> > counter(kmerLength);
        counter.readOccurrenceTableBinary(optionSingleArgs["-k"]);
        this->getOccurrenceArray(counter);
    } else if (kmerLength <= 96) {
        Counter<KmerN<Binstr95> > counter(kmerLength);
        counter.readOccurrenceTableBinary(optionSingleArgs["-k"]);
        this->getOccurrenceArray(counter);
    } else if (kmerLength <= 128) {
        Counter<KmerN<Binstr127> > counter(kmerLength);
        counter.readOccurrenceTableBinary(optionSingleArgs["-k"]);
        this->getOccurrenceArray(counter);
    } else if (kmerLength <= 160) {
        Counter<KmerN<Binstr159> > counter(kmerLength);
        counter.readOccurrenceTableBinary(optionSingleArgs["-k"]);
        this->getOccurrenceArray(counter);
    } else {
        Counter<KmerN<binstr_t> > counter(kmerLength);
        counter.readOccurrenceTableBinary(optionSingleArgs["-k"]);
        this->getOccurrenceArray(counter);
    }

    this->decideContigBreakPoint();
    std::string outputFilename = optionSingleArgs["-o"]; 
	if (!optionBool["-recalc_cov"]) {
		outputFilename += "_cov_trimmed.fa";
	}
	else {
		outputFilename += "_recalc.fa";
	}

    this->divideAndPrintContig(outputFilename);
}



/*
//////////////////////////////////////////////////////////////////////////////////////
// get kmer length from binary file
//////////////////////////////////////////////////////////////////////////////////////
unsigned long long ContigDivider::getKmerLengthFromBinary(const std::string &filename)
{
    unsigned long long kmerLength;
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    if (!ifs) throw platanus::FILEError(filename);

    ifs.read(reinterpret_cast<char *>(&kmerLength), sizeof(unsigned long long));
    ifs.close();

    return kmerLength;
}
*/



//////////////////////////////////////////////////////////////////////////////////////
// exec divide
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void ContigDivider::getOccurrenceArray(const Counter<KMER> &counter)
{
    unsigned long long kmerLength = counter.getKmerLength();
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));
    KMER kmer(kmerLength);
    contigItself.reset(new OccurrenceArray(kmerLength, this->contig));

    for (unsigned long long i = 0; i < contigItself->getNumSeq(); ++i) {
        unsigned long long numKmer = contigItself->getNumKmer(i);
        if (numKmer == 0) continue;

        const platanus::SEQ &seq = contig.seq[i];
        unsigned long long start = 0;
        bool isInit = true;
        while (start < numKmer) {
            if (isInit) {
                unsigned long long j = 0;
                for (; j < kmerLength - 1; ++j) {
                    if (seq.base[start + j] == 4) break;

                    kmer.setForward(kmerLength - 2 - j, seq.base[start + j]);
                    kmer.setReverse(j + 1, 0x3 ^ seq.base[start + j]);
                }

                if (j == kmerLength - 1)
                    isInit = false;
                else {
                    start += j + 1;
                    continue;
                }
            }

            if (seq.base[start + kmerLength - 1] == 4) {
                start += kmerLength;
                isInit = true;
                continue;
            }
            kmer.forward <<= 2;
            kmer.maskForward(mask);
            kmer.reverse >>= 2;
            kmer.setForward(0, seq.base[start + kmerLength - 1]);
            kmer.setReverse(kmerLength - 1, 0x3 ^ seq.base[start + kmerLength - 1]);
            contigItself->setOccurrenceArray(i, start, std::max(counter.findValue(std::min(kmer.forward, kmer.reverse)), static_cast<unsigned short>(1)));
            ++start;
        }
    }
}


void ContigDivider::decideContigBreakPoint(void)
{
	const double SD_THRESHOLD = 2.0;
    const unsigned long long kmerLength = contigItself->getKmerLength();
    const unsigned long long windowSize = 1000 - kmerLength + 1;
    this->breakPoint.resize(this->contig.numSeq);

	for (unsigned long long i = 0; i < this->contig.numSeq; ++i) {
		if (contig.seq[i].length < 2 * (windowSize + kmerLength - 1))
			continue;

		std::vector<unsigned long long> kmerCoverage = contigItself->getOccurrenceArray(i);
		std::vector<double> windowMean(contig.seq[i].length - kmerLength - windowSize + 2, 0.0);
		std::vector<double> windowSD(contig.seq[i].length - kmerLength - windowSize + 2, 0.0);

		for (unsigned long long j = 0; j < windowSize; ++j) {
			windowMean[0] += kmerCoverage[j];
			windowSD[0] += std::pow(kmerCoverage[j], 2.0);
		}

		for (unsigned long long j = 1; j < kmerCoverage.size() - windowSize + 1; ++j) {
			windowMean[j] = windowMean[j - 1] - kmerCoverage[j - 1] + kmerCoverage[j + windowSize - 1];
			windowSD[j] = windowSD[j - 1] - std::pow(kmerCoverage[j - 1], 2.0) + std::pow(kmerCoverage[j + windowSize - 1], 2.0);
		}
		for (unsigned long long j = 0; j < kmerCoverage.size() - windowSize + 1; ++j) {
			windowMean[j] /= (double)windowSize;
			windowSD[j] = std::sqrt(windowSD[j] / windowSize - std::pow(windowMean[j], 2.0));
		}

		for (unsigned long long j = 0; j < kmerCoverage.size() - (2 * windowSize + kmerLength) + 1; ++j) {
            if (!(this->withinRatioThreshold(windowMean[j], windowMean[j + windowSize + kmerLength])) &&
				std::abs(windowMean[j] - windowMean[j + windowSize + kmerLength]) > SD_THRESHOLD * std::max(windowSD[j], windowSD[j + windowSize + kmerLength])) {

                this->breakPoint[i].emplace_back(j + windowSize + kmerLength - 1);
			}
		}
	}

    for (unsigned long long i = 0; i < this->contig.numSeq; ++i) {
        std::vector<unsigned long long> kmerCoverage = contigItself->getOccurrenceArray(i);
        const unsigned long long median = this->findMedian(kmerCoverage);
        this->breakPoint[i].emplace_back(0);

        for (unsigned long long j = 0; j < contig.seq[i].length - kmerLength + 1; ++j) {
            if (this->withinRatioThreshold(kmerCoverage[j], median)) {
                this->breakPoint[i].emplace_back(j);
                break;
            }
        }

        for (unsigned long long j = contig.seq[i].length - kmerLength; j >= 0; --j) {
            if (this->withinRatioThreshold(kmerCoverage[j], median)) {
                this->breakPoint[i].emplace_back(j + kmerLength);
                break;
            }
        }
        this->breakPoint[i].emplace_back(this->contig.seq[i].length);
    }
}


unsigned long long ContigDivider::findMedian(const std::vector<unsigned long long> &vec) const
{
    std::vector<unsigned long long> tmp = vec;
    std::sort(tmp.begin(), tmp.end());
    return tmp[tmp.size() / 2];
}


void ContigDivider::divideAndPrintContig(const std::string &filename)
{
    unsigned long long minLength = this->contigMaxK;
    unsigned long long seqID = 1;
    std::ofstream ofs(filename.c_str());

    for (unsigned long long contigID = 0; contigID < this->contig.numSeq; ++contigID) {
        std::sort(this->breakPoint[contigID].begin(), this->breakPoint[contigID].end());

        for (auto itr = breakPoint[contigID].begin(), end = breakPoint[contigID].end(); itr != end; ++itr) {
            auto next = itr + 1;
            if (next == end)
				break;

            if (*next - *itr < minLength) continue;
			if (*next - *itr < minLength) {
				while (next + 1 != end && *(next + 1) - *next < minLength)
					++next;
			}

            if (*next - *itr < minLength)
				continue;

            unsigned long long coverage = this->calcAverageCoverage(contigID, *itr, *next);
            ofs << ">seq" << seqID << "_len" << *next - *itr << "_cov" << coverage << "_read" << this->readLength << "_maxK" << this->contigMaxK << "\n";
            unsigned long long j = 0;
            for (unsigned long long pos = *itr; pos < *next; ++pos, ++j) {
                ofs.put(platanus::Bin2Char(this->contig.seq[contigID].base[pos]));
                if ((j + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                    ofs.put('\n');
            }
            if (j % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
                ofs.put('\n');
        }
        ++seqID;
    }
    ofs.close();
}


unsigned long long ContigDivider::calcAverageCoverage(const unsigned long long contigID, const unsigned long long start, const unsigned long long end)
{
    double coverage = 0;
    unsigned long long kmerLength = contigItself->getKmerLength();
    for (unsigned long long i = start; i < end - kmerLength + 1; ++i) {
        coverage += contigItself->getOccurrence(contigID, i);
    }

    coverage /= (end - start + kmerLength + 1);
    return static_cast<unsigned long long>((coverage * this->readLength) / (this->readLength - kmerLength + 1) + 0.5);
}








