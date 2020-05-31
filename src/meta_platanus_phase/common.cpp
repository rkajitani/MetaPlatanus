/*
Copyright (C) 2014 Itoh Laboratory, Tokyo Institute of Technology

This file is part of Platanus.

Platanus is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Platanus is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with Platanus; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "common.h"

//////////////////////////////////////////////////////////////////////////////////////
// static constant
//////////////////////////////////////////////////////////////////////////////////////
const unsigned long long platanus::ConstParam::MAX_READ_LEN = 500000;
const unsigned platanus::ConstParam::SCAFFOLD_HASH_OVERLAP = 32;
const unsigned platanus::ConstParam::OUTPUT_LINE_LENGTH = 80;
const unsigned platanus::ConstParam::MAX_FILE_NUM = 100;
const unsigned platanus::ConstParam::MAX_FILE_LEN = 200;
const unsigned platanus::ConstParam::MAX_THREAD = 100;
const std::string platanus::ConstParam::VERSION = "1.2.2";
const double platanus::ConstParam::DOUBLE_HASH_MAX_LOAD_FACTOR = 0.9;
const unsigned long long platanus::ConstParam::DEFAULT_CONTIG_READ_LEN = 100;
const double platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR = 0.25;
const double platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR = 1.75;
const double platanus::ConstParam::LONG_READ_INS_SIZE_LOWER_BOUND_FACTOR = 0.25;
const double platanus::ConstParam::LONG_READ_INS_SIZE_UPPER_BOUND_FACTOR = 10.0;
const double platanus::ConstParam::DEFAULT_EMP_DISTR_BIN_SIZE = 0.0001;
std::string platanus::globalTmpFileDir = ".";


void platanus::setLimitNumFileOpenMax()
{
	struct rlimit rl;

	getrlimit(RLIMIT_NOFILE, &rl);

	rl.rlim_cur = rl.rlim_max;
	setrlimit(RLIMIT_NOFILE, &rl);
}



void platanus::setFastaFileNameAndNumber(const vector<string> fastaFileName, std::vector<std::string> &scaffoldFileName, std::vector<long> &numScaffoldInFile)
{
	scaffoldFileName.assign(fastaFileName.size(), "");
	numScaffoldInFile.assign(fastaFileName.size(), 0);

	for (unsigned fileIndex = 0; fileIndex < fastaFileName.size(); ++fileIndex) {
		scaffoldFileName[fileIndex] = fastaFileName[fileIndex];
		std::ifstream ifs(fastaFileName[fileIndex]);
		if (!ifs)
			continue;

		string line;
		while (ifs && getline(ifs, line)) {
			if (line[0] == '>')
				++numScaffoldInFile[fileIndex];
		}
		ifs.close();
	}
}


unsigned long long platanus::getKmerLengthFromBinary(const std::string &filename)
{
	unsigned long long kmerLength;

	std::ifstream ifs(filename.c_str(), std::ios::binary);
	if (!ifs) throw FILEError(filename);

	ifs.read(reinterpret_cast<char *>(&kmerLength), sizeof(unsigned long long));
	ifs.close();

	return kmerLength;
}


void platanus::Contig::execMga(const string &inFilename, const string &outFilename, const string &mgaExecutable)
{
	std::ostringstream oss;

	oss << mgaExecutable << " " << inFilename << " -m  >" << outFilename;

	std::cerr << "Executing MeteGeneAnnotator..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl;

	system(oss.str().c_str());
	if (system(oss.str().c_str()) != 0) {
		throw platanus::MgaError();
	}
}


void platanus::Contig::readMgaResult(const string &mgaOutFilename)
{
	std::string oneLine, buf, seqName, strand, startStopFlags;
	long start, end, frame;
	platanus::SEQ *seqP = NULL;

	std::ifstream ifs(mgaOutFilename);

	while (ifs && getline(ifs, oneLine)) {
		std::istringstream ss(oneLine);
		if (oneLine[0] == '#') {
			ss >> buf >> seqName;
			getline(ifs, oneLine);
			getline(ifs, oneLine);

			auto tmpIt = nameIndex.find(seqName);
			if (tmpIt != nameIndex.end()) {
				seqP = &(seq[tmpIt->second]);
			}
			else {
				seqP = NULL;
			}
		}
		else if (seqP != NULL) {
			ss >> buf >> start >> end >> strand >> frame >> startStopFlags;
			--start;

			if (seqP->CDS.size() > 0)
				seqP->CDS.append("\4\4\4");

			if (strand == "+") {
				if (startStopFlags[0] == '1')
					start += 3;
				else
					start += frame;

				if (startStopFlags[1] == '1')
					end -= 3;
				else
					end -= (end - start) % 3;

				seqP->CDS.append(seqP->base.substr(start, end - start));
			}
			else {
				if (startStopFlags[0] == '1')
					end -= 3;
				else
					end -= frame;

				if (startStopFlags[1] == '1')
					start += 3;
				else
					start += (end - start) % 3;

				for (long i = end - 1; i >= start; --i)
					seqP->CDS.push_back(seqP->base[i] != 4 ? seqP->base[i] ^ 0x3 : 4);
			}
		}

	}
		
	ifs.close();
}


void platanus::empiricalDistribution::count(double value)
{
	unsigned index = static_cast<unsigned>(value / binSize);

	if (histo.size() < index + 1)
		histo.resize(index + 1, 0);
	
	++histo[index];
}


void platanus::empiricalDistribution::setDistribution()
{
	if (histo.empty())
		return;

	unsigned sum = 0;
	for (auto it = histo.begin(); it != histo.end(); ++it)
		sum += *it;

	distribution.clear();
	distribution.resize(histo.size());
	distribution[0] = histo[0];
	for (unsigned i = 1; i < histo.size(); ++i)
		distribution[i] = distribution[i - 1] + (static_cast<double>(histo[i]) / sum);
}

double platanus::empiricalDistribution::probGreater(const double value)
{
	if (static_cast<unsigned>(value / binSize) < distribution.size())
		return (1.0 - distribution[static_cast<unsigned>(value / binSize)]);
	else 
		return 0.0;
}
