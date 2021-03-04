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

#include "pairedDBG.h"
#include <omp.h>
#include <climits>
#include <cfloat>
#include <string>
#include <sstream>
#include <queue>
#include <array>

using std::vector;
using std::string;
using std::pair;
using std::cerr;
using std::endl;


const unsigned PairedDBG::MIN_LENGTH_FOR_MGC = 1000;
const unsigned PairedDBG::MIN_LENGTH_FOR_COV_ESTIMATE = 1000;
//const unsigned PairedDBG::MIN_LENGTH_FOR_MGC = 200;
//const unsigned PairedDBG::MIN_LENGTH_FOR_COV_ESTIMATE = 200;
const double PairedDBG::EMP_DISTR_BIN_SIZE = 0.0001;

// for node-state
const unsigned PairedDBG::DBG_HETERO = 0x8;
const unsigned PairedDBG::DBG_PRIMARY_BUBBLE = 0x10;
const unsigned PairedDBG::DBG_SECONDARY_BUBBLE = 0x20;

// for edge-state
const char PairedDBG::DBG_OVERLAP = 0x1;

// for contig-state
const char PairedDBG::DBG_CONTIG_BUBBLE_JUNCTION = 0x1;
const char PairedDBG::DBG_CONTIG_PRIMARY_BUBBLE = 0x2;
const char PairedDBG::DBG_CONTIG_SECONDARY_BUBBLE = 0x4;

// related to scaffolding modes
const unsigned PairedDBG::OVERLAP_MODE = 0x1;
const unsigned PairedDBG::PAIRED_END_LINK_MODE = 0x2;
const unsigned PairedDBG::LENGTH_CUTOFF_MODE = 0x4;
const unsigned PairedDBG::NON_DBG_OVERLAP_MODE = 0x8;
const unsigned PairedDBG::LONG_READ_LINK_MODE = 0x10;
const unsigned PairedDBG::BUBBLE_AWARE_MODE = 0x20;
const unsigned PairedDBG::SECONDARY_BUBBLE_REMOVAL_MODE = 0x40;
const unsigned PairedDBG::PREVIOUS_DIVISION_AWARE_MODE = 0x80;
const unsigned PairedDBG::TAG_SCAFFOLD_MODE = 0x100;

const long PairedDBG::MAX_ITERATION_OF_CROSS_SOLUTION = 5;

const double PairedDBG::NON_PAIRED_END_TOLERENCE_FACTOR = 0.1;
//const double PairedDBG::HETERO_COVERAGE_THRESHOLD_FACTOR = 1.75;
const double PairedDBG::HETERO_COVERAGE_THRESHOLD_FACTOR = USHRT_MAX;
const double PairedDBG::CROSS_LINK_RATE_THRESHOLD = 0.25;
const double PairedDBG::CROSS_SCORE_RATE_THRESHOLD = 0.5;
const double PairedDBG::MIN_BUBBLE_COUNT_FACTOR = 10000.0;


PairedDBG::PairedDBG(): ScaffoldGraph(), minTagIslandLengthToScaffold(0), heteroCoverage(0.0), tolerenceFactor(1.0), cutoffLengthFactor(0.0), contigPreviousParentNodeID(0), oppositeBubbleContigID(0), mode(0x1 | 0x2), numInputBubbleContig(0), mgcEmpDistr(EMP_DISTR_BIN_SIZE), covEmpDistr(EMP_DISTR_BIN_SIZE) {}


void PairedDBG::storeGraphLinkFromOverlap(vector<GraphLinkWithFlag> &graphLinkPool)
{
    # pragma omp parallel for schedule(dynamic)
    for (unsigned d = 0; d < TABLE_DIVID; ++d) {
		auto overlapIterator = overlapTable[d].begin();
		auto overlapEnd = overlapTable[d].end();
		for (; overlapIterator != overlapEnd; ++overlapIterator) {
			GraphLinkWithFlag overlapLink;
			overlapLink.overlapFlag = true;
			overlapLink.id1 = overlapIterator->second.id1;
			overlapLink.id2 = overlapIterator->second.id2;

            long i = id2Index(overlapLink.id1);
            overlapLink.id1 = overlapLink.id1 > 0 ? contigPositionInScaffold[i].id : -(contigPositionInScaffold[i].id);

            long j = id2Index(overlapLink.id2);
            overlapLink.id2 = overlapLink.id2 > 0 ? contigPositionInScaffold[j].id : -(contigPositionInScaffold[j].id);

			if (overlapLink.id1 * overlapLink.id2 == 0 || abs(overlapLink.id1) == abs(overlapLink.id2))
				continue;

			overlapLink.gap = -(getScaffoldOverlap(overlapLink.id1, overlapLink.id2));

			if (overlapLink.gap == -(this->minOverlap)) {
				if (abs(overlapLink.id1) > abs(overlapLink.id2)) {
					std::swap(overlapLink.id1, overlapLink.id2);
					overlapLink.id1 *= -1;
					overlapLink.id2 *= -1;
				}
				#pragma omp critical
				{
					graphLinkPool.push_back(overlapLink);
				}
			}
		}
    }
}

void PairedDBG::storeGraphLinkFromMappedPair(vector<GraphLinkWithFlag> &graphLinkPool, long numThread)
{
	if (!(this->mode & PAIRED_END_LINK_MODE) || allLibraryMT == NULL)
		return;

	# pragma omp parallel for schedule(static, 1)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		platanus::Position forwardResult;
		platanus::Position reverseResult;
		unsigned numPairLink;
		unsigned numReadLink;

		rewind((*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
		while (fread(&numPairLink, sizeof(unsigned), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);

				long i = id2Index(forwardResult.id);
				if (contigPositionInScaffold[i].id == 0) continue;
				forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[i].id : -(contigPositionInScaffold[i].id);
				forwardResult.offset = contigPositionInScaffold[i].id > 0 ? forwardResult.offset : contig[i].length - forwardResult.offset - 1;
				forwardResult.offset += node[id2Index(forwardResult.id)].contig[contigPositionInScaffold[i].offset].start;

				long j = id2Index(reverseResult.id);
				if (contigPositionInScaffold[j].id == 0) continue;
				reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[j].id : -(contigPositionInScaffold[j].id);
				reverseResult.offset = contigPositionInScaffold[j].id > 0 ? reverseResult.offset : contig[j].length - reverseResult.offset - 1;
				reverseResult.offset += node[id2Index(reverseResult.id)].contig[contigPositionInScaffold[j].offset].start;

				if (abs(forwardResult.id) == abs(reverseResult.id))
					continue;

				GraphLinkWithFlag graphLink;
				graphLink.overlapFlag = false;
				if (abs(forwardResult.id) < abs(reverseResult.id)) {
					graphLink.id1 = forwardResult.id;
					graphLink.offset1 = contigPositionInScaffold[i].offset;
					graphLink.id2 = -(reverseResult.id);
					graphLink.offset2 = contigPositionInScaffold[j].offset;
				}
				else {
					graphLink.id1 = reverseResult.id;
					graphLink.offset1 = contigPositionInScaffold[j].offset;
					graphLink.id2 = -(forwardResult.id);
					graphLink.offset2 = contigPositionInScaffold[i].offset;
				}

				std::pair<int, int> redundancyCheckKey = std::make_pair(graphLink.id1, graphLink.id2);
				if (redundancyCheckSet.find(redundancyCheckKey) != redundancyCheckSet.end())
					continue;

				redundancyCheckSet.insert(redundancyCheckKey);


				// calc gap length (average length minus own node length plus offset)
				// image
				//                    ---------average gap----------
				// -----------------|node length |NNNNNNNNNNNNNNNNNN----------------
				//                  -- <- offset
				graphLink.gap = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();
				long gapMinus;
				if (forwardResult.id > 0) {
					gapMinus = node[forwardResult.id - 1].length - forwardResult.offset;
				} else {
					gapMinus = forwardResult.offset + 1;
				}
//				if (node[id2Index(forwardResult.id)].length < cutoffLength) continue;
				graphLink.gap -= gapMinus;


				if (reverseResult.id > 0) {
	//				if (node[reverseResult.id - 1].length < cutoffLength || abs(forwardResult.id) == abs(reverseResult.id)) continue;
					if (abs(forwardResult.id) == abs(reverseResult.id)) continue;
					graphLink.gap -= node[reverseResult.id-1].length - reverseResult.offset;
				} else {
	//				if (node[(-1 * reverseResult.id) - 1].length < cutoffLength || abs(forwardResult.id) == abs(reverseResult.id)) continue;
					if (abs(forwardResult.id) == abs(reverseResult.id)) continue;
					graphLink.gap -= reverseResult.offset + 1;
				}

				if (-graphLink.gap > std::min(tolerence, std::min(node[id2Index(graphLink.id1)].length, node[id2Index(graphLink.id2)].length)) + this->getScaffoldOverlap(graphLink.id1, graphLink.id2))
					continue;

				# pragma omp critical (push)
				{
					graphLinkPool.push_back(graphLink);
				}
			}


			for (unsigned pairIndex = 0; pairIndex < numReadLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);

				long i = id2Index(forwardResult.id);
				if (contigPositionInScaffold[i].id == 0) continue;
				forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[i].id : -(contigPositionInScaffold[i].id);
				forwardResult.offset = contigPositionInScaffold[i].id > 0 ? forwardResult.offset : contig[i].length - forwardResult.offset - 1;
				forwardResult.offset += node[id2Index(forwardResult.id)].contig[contigPositionInScaffold[i].offset].start;

				long j = id2Index(reverseResult.id);
				if (contigPositionInScaffold[j].id == 0) continue;
				reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[j].id : -(contigPositionInScaffold[j].id);
				reverseResult.offset = contigPositionInScaffold[j].id > 0 ? reverseResult.offset : contig[j].length - reverseResult.offset - 1;
				reverseResult.offset += node[id2Index(reverseResult.id)].contig[contigPositionInScaffold[j].offset].start;

				if (abs(forwardResult.id) == abs(reverseResult.id))
					continue;

				GraphLinkWithFlag graphLink;
				graphLink.overlapFlag = false;

				long forwardLeft;
				long forwardRight;
				if (forwardResult.id > 0) {
					forwardLeft = -(forwardResult.offset);
					forwardRight = node[forwardResult.id - 1].length - forwardResult.offset - 1;
				}
				else {
					forwardLeft = -(node[-(forwardResult.id) - 1].length - forwardResult.offset - 1);
					forwardRight = forwardResult.offset;
				}

				long reverseLeft;
				long reverseRight;
				if (reverseResult.id > 0) {
					reverseLeft = -(reverseResult.offset);
					reverseRight = node[reverseResult.id - 1].length - reverseResult.offset - 1;
				}
				else {
					reverseLeft = -(node[-(reverseResult.id) - 1].length - reverseResult.offset - 1);
					reverseRight = reverseResult.offset;
				}

				if (forwardLeft <= reverseLeft) {
					if (forwardRight > reverseRight)
						continue;
					graphLink.gap = -(forwardRight - reverseLeft + 1);
				}
				else {
					if (reverseRight > forwardRight)
						continue;
					graphLink.gap = -(reverseRight - forwardLeft + 1);
				}

				if (abs(forwardResult.id) < abs(reverseResult.id)) {
					if (forwardRight < reverseRight) {
						graphLink.id1 = forwardResult.id;
						graphLink.offset1 = contigPositionInScaffold[i].offset;
						graphLink.id2 = reverseResult.id;
						graphLink.offset2 = contigPositionInScaffold[j].offset;
					}
					else {
						graphLink.id1 = -(forwardResult.id);
						graphLink.offset1 = contigPositionInScaffold[i].offset;
						graphLink.id2 = -(reverseResult.id);
						graphLink.offset2 = contigPositionInScaffold[j].offset;
					}
				}
				else {
					if (forwardRight < reverseRight) {
						graphLink.id1 = -(reverseResult.id);
						graphLink.offset1 = contigPositionInScaffold[j].offset;
						graphLink.id2 = -(forwardResult.id);
						graphLink.offset2 = contigPositionInScaffold[i].offset;
					}
					else {
						graphLink.id1 = reverseResult.id;
						graphLink.offset1 = contigPositionInScaffold[j].offset;
						graphLink.id2 = forwardResult.id;
						graphLink.offset2 = contigPositionInScaffold[i].offset;
					}
				}

				std::pair<int, int> redundancyCheckKey = std::make_pair(graphLink.id1, graphLink.id2);
				if (redundancyCheckSet.find(redundancyCheckKey) != redundancyCheckSet.end())
					continue;

				redundancyCheckSet.insert(redundancyCheckKey);


/*
				// check contigs aren't overlapeed too large area
				if (-(graphLink.gap) > this->minOverlap) continue;
				# pragma omp critical (push)
				{
					graphLinkPool.push_back(graphLink);
				}
*/
			}
		}
	}
}


void PairedDBG::storeGraphLinkFromMappedLongRead(vector<GraphLinkWithFlag> &graphLinkPool, long numThread)
{
	if (!(this->mode & LONG_READ_LINK_MODE) || longReadLibraryMT == NULL)
		return;
	
	# pragma omp parallel for schedule(static, 1)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		unsigned numResult;
		platanus::Position mapResult;
		platanus::Position forwardResult;
		platanus::Position reverseResult;
		long score;

		rewind((*longReadLibraryMT)[threadID].mappedReadFP);
		while (fread(&numResult, sizeof(unsigned), 1, (*longReadLibraryMT)[threadID].mappedReadFP)) {

			vector<MapInfoForGraph> mapInfoBuffer;
			for (unsigned long i = 0; i < numResult; ++i) {
				fread(&mapResult, sizeof(platanus::Position), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
				fread(&score, sizeof(long), 1, (*longReadLibraryMT)[threadID].mappedReadFP);

				long contigIndex = id2Index(mapResult.id);
				if (contigPositionInScaffold[contigIndex].id == 0) continue;
				mapResult.id = mapResult.id > 0 ? contigPositionInScaffold[contigIndex].id : -(contigPositionInScaffold[contigIndex].id);
				mapResult.offset = contigPositionInScaffold[contigIndex].id > 0 ? mapResult.offset : contig[contigIndex].length - mapResult.offset - 1;
				mapResult.offset += node[id2Index(mapResult.id)].contig[contigPositionInScaffold[contigIndex].offset].start;

				mapInfoBuffer.emplace_back(mapResult, contigIndex, score);
			}
			
			if (mapInfoBuffer.empty())
				continue;

			std::stable_sort(mapInfoBuffer.begin(), mapInfoBuffer.end(), MapInfoForGraph::PositionIDLessScoreGreater());

			vector<MapInfoForGraph> mergedMapResultBuffer;
			auto maxScoreItr = mapInfoBuffer.begin();
			while (maxScoreItr != mapInfoBuffer.end()) {
				auto boundItr = std::upper_bound(maxScoreItr, mapInfoBuffer.end(), *(maxScoreItr), MapInfoForGraph::PositionIDLess());

				score = 0;
				for (auto itr = maxScoreItr; itr != boundItr; ++itr)
					score += itr->score;

				mergedMapResultBuffer.emplace_back(maxScoreItr->position, maxScoreItr->contigIndex, score);

				maxScoreItr = boundItr;
			}

			for (auto mapItr1 = mergedMapResultBuffer.begin(); mapItr1 != mergedMapResultBuffer.end() - 1; ++mapItr1) {
				for (auto mapItr2 = mapItr1 + 1; mapItr2 != mergedMapResultBuffer.end(); ++mapItr2) {
					if (abs(mapItr1->position.id) == abs(mapItr2->position.id))
						continue;

					GraphLinkWithFlag graphLink;
					graphLink.overlapFlag = false;

					long forwardLeft;
					long forwardRight;
					if (mapItr1->position.id > 0) {
						forwardLeft = -(mapItr1->position.offset);
						forwardRight = node[mapItr1->position.id - 1].length - mapItr1->position.offset - 1;
					}
					else {
						forwardLeft = -(node[-(mapItr1->position.id) - 1].length - mapItr1->position.offset - 1);
						forwardRight = mapItr1->position.offset;
					}

					long reverseLeft;
					long reverseRight;
					if (mapItr2->position.id > 0) {
						reverseLeft = -(mapItr2->position.offset);
						reverseRight = node[mapItr2->position.id - 1].length - mapItr2->position.offset - 1;
					}
					else {
						reverseLeft = -(node[-(mapItr2->position.id) - 1].length - mapItr2->position.offset - 1);
						reverseRight = mapItr2->position.offset;
					}

					if (forwardLeft <= reverseLeft) {
						if (forwardRight > reverseRight)
							continue;
						graphLink.gap = -(forwardRight - reverseLeft + 1);
					}
					else {
						if (reverseRight > forwardRight)
							continue;
						graphLink.gap = -(reverseRight - forwardLeft + 1);
					}

					if (abs(mapItr1->position.id) < abs(mapItr2->position.id)) {
						if (forwardRight < reverseRight) {
							graphLink.id1 = mapItr1->position.id;
							graphLink.offset1 = contigPositionInScaffold[mapItr1->contigIndex].offset;
							graphLink.id2 = mapItr2->position.id;
							graphLink.offset2 = contigPositionInScaffold[mapItr2->contigIndex].offset;
						}
						else {
							graphLink.id1 = -(mapItr1->position.id);
							graphLink.offset1 = contigPositionInScaffold[mapItr1->contigIndex].offset;
							graphLink.id2 = -(mapItr2->position.id);
							graphLink.offset2 = contigPositionInScaffold[mapItr2->contigIndex].offset;
						}
					}
					else {
						if (forwardRight < reverseRight) {
							graphLink.id1 = -(mapItr2->position.id);
							graphLink.offset1 = contigPositionInScaffold[mapItr2->contigIndex].offset;
							graphLink.id2 = -(mapItr1->position.id);
							graphLink.offset2 = contigPositionInScaffold[mapItr1->contigIndex].offset;
						}
						else {
							graphLink.id1 = mapItr2->position.id;
							graphLink.offset1 = contigPositionInScaffold[mapItr2->contigIndex].offset;
							graphLink.id2 = mapItr1->position.id;
							graphLink.offset2 = contigPositionInScaffold[mapItr1->contigIndex].offset;
						}
					}

					if (-(graphLink.gap) > this->tolerence)
						continue;
					
					graphLink.score = mapItr1->score + mapItr2->score;

					# pragma omp critical (push)
					{
						graphLinkPool.push_back(graphLink);
					}
				}
			}
		}
	}

}


void PairedDBG::calcLink(const long libraryIndex, const long linkThreshold, const long numThread)
{
    if (graphLinkFP[libraryIndex] != NULL)
        fclose(graphLinkFP[libraryIndex]);
    graphLinkFP[libraryIndex] = platanus::makeTemporaryFile();

//    cerr << "linking scaffolds (MIN_LINK = " << minLink << ")" << endl;
    cerr << "linking scaffolds..." << endl;

    vector<GraphLinkWithFlag> graphLinkPool;

	if (this->mode & OVERLAP_MODE)
		storeGraphLinkFromOverlap(graphLinkPool);

	if (this->mode & PAIRED_END_LINK_MODE)
		storeGraphLinkFromMappedPair(graphLinkPool, numThread);

	if (this->mode & LONG_READ_LINK_MODE)
		storeGraphLinkFromMappedLongRead(graphLinkPool, numThread);

	if (this->mode & TAG_SCAFFOLD_MODE)
		storeGraphLinkFromTagReadPair(graphLinkPool, true, numThread);

    long totalLink = graphLinkPool.size();

	if (totalLink == 0)
		return;

    cerr << "sorting links in contigID order..." << endl;
    std::stable_sort(graphLinkPool.begin(), graphLinkPool.end());
    graphLinkPool.emplace_back();
    graphLinkPool[totalLink].id1 = 0;

    cerr << "estimating contig distances..." << endl;

    vector<GraphLinkWithFlagPoolIndex> indexes(1, 0);
	++indexes.back().numLink;
	if (graphLinkPool[0].overlapFlag)
		indexes.back().overlapFlag = true;
    for (long idx = 1; idx < totalLink; ++idx) {
        if (graphLinkPool[idx - 1].id1 != graphLinkPool[idx].id1 || graphLinkPool[idx - 1].id2 != graphLinkPool[idx].id2) {
            if (!(indexes.back().overlapFlag) && indexes.back().numLink < linkThreshold) {
                indexes.pop_back();
            }
            indexes.emplace_back(idx);
        }

        ++indexes.back().numLink;
		if (graphLinkPool[idx].overlapFlag) {
			indexes.back().overlapFlag = true;
			indexes.back().gap = graphLinkPool[idx].gap;
		}
    }
    if (!(indexes.back().overlapFlag) && indexes.back().numLink < linkThreshold) {
        indexes.pop_back();
    }
    sort(indexes.begin(), indexes.end(), GraphLinkWithFlagPoolIndexGreater());

	if (numThread == 1) {
		for (unsigned long idxIndex = 0; idxIndex < indexes.size(); ++idxIndex) {
			calcLinkAndWriteGraphLinkWithFlagFile(graphLinkPool, indexes[idxIndex], libraryIndex, false);
		}
	}
	else {
		// pragma omp parallel for schedule(dynamic) num_threads(numThread)
		for (unsigned long idxIndex = 0; idxIndex < indexes.size(); ++idxIndex) {
			calcLinkAndWriteGraphLinkWithFlagFile(graphLinkPool, indexes[idxIndex], libraryIndex, true);
		}
	}
}

void PairedDBG::calcLinkAndWriteGraphLinkWithFlagFile(const std::vector<GraphLinkWithFlag>& graphLinkPool, const GraphLinkWithFlagPoolIndex& index, const long libraryIndex, const bool multiThreadFlag)
{
    unsigned long indexBegin = index.index;
    long numLink = index.numLink;
	long score = 0;

    GraphLinkWithFlag graphLink = graphLinkPool[indexBegin];
	graphLink.overlapFlag = index.overlapFlag;
    vector<long> breakdown1(node[id2Index(graphLink.id1)].numContig, 0);
    vector<long> breakdown2(node[id2Index(graphLink.id2)].numContig, 0);
    for (long idxLink = indexBegin, endLink = idxLink + numLink; idxLink < endLink; ++idxLink) {
        ++(breakdown1[graphLinkPool[idxLink].offset1]);
        ++(breakdown2[graphLinkPool[idxLink].offset2]);
    }

	if (index.overlapFlag) {
		graphLink.gap = index.gap;
	}
	else {
		vector<GraphLink> targetLinkPool(numLink);
		for (long i = 0; i < numLink; ++i) {
			targetLinkPool[i] = static_cast<GraphLink>(graphLinkPool[i + indexBegin]);
			score += graphLinkPool[i + indexBegin].score;
		}
//		graphLink.gap = estimateGapSize(targetLinkPool, 0, numLink);
		graphLink.gap = estimateGapSizeAverage(targetLinkPool, 0, numLink);
	}
    if (graphLink.gap == LONG_MIN) return;

    long breakdownSize = breakdown1.size() + breakdown2.size();
    std::unique_ptr<long[]> tmp(new long[breakdownSize]());
    std::copy(breakdown1.begin(), breakdown1.end(), tmp.get());
    std::copy(breakdown2.begin(), breakdown2.end(), &(tmp[breakdown1.size()]));

	if (multiThreadFlag) {
		#pragma omp critical (write_link)
		{
			fwrite(&numLink, sizeof(long), 1, graphLinkFP[libraryIndex]);
			fwrite(&score, sizeof(long), 1, graphLinkFP[libraryIndex]);
			fwrite(&graphLink, sizeof(GraphLinkWithFlag), 1, graphLinkFP[libraryIndex]);
			fwrite(tmp.get(), sizeof(long), breakdownSize, graphLinkFP[libraryIndex]);
		}
	}
	else {
		fwrite(&numLink, sizeof(long), 1, graphLinkFP[libraryIndex]);
		fwrite(&score, sizeof(long), 1, graphLinkFP[libraryIndex]);
		fwrite(&graphLink, sizeof(GraphLinkWithFlag), 1, graphLinkFP[libraryIndex]);
		fwrite(tmp.get(), sizeof(long), breakdownSize, graphLinkFP[libraryIndex]);
	}
}

void PairedDBG::makeGraph(const long numThread)
{
	unsigned numPairedEndLibrary = 0;
	if (allLibraryMT != NULL)
		numPairedEndLibrary = (*allLibraryMT).size();

	unsigned numLongReadLibrary = 0;
	if (longReadLibraryMT != NULL)
		numLongReadLibrary = 1;

	if (graphLinkFP.empty())
		graphLinkFP.resize(numPairedEndLibrary + numLongReadLibrary, NULL);

    calcLink(this->targetLibraryIndex, this->minLink, numThread);

    cerr << "constructing paired de Bruijn graph" << endl;

    for (int i = 0; i < numNode; ++i) {
		node[i].edge.clear();
		node[i].numEdge = 0;
	}

    rewind(graphLinkFP[this->targetLibraryIndex]);
    long numLink;
    long score;
    while (fread(&numLink, sizeof(long), 1, graphLinkFP[this->targetLibraryIndex])) {
		fread(&score, sizeof(long), 1, graphLinkFP[this->targetLibraryIndex]);
        GraphLinkWithFlag link;
        fread(&link, sizeof(GraphLinkWithFlag), 1, graphLinkFP[this->targetLibraryIndex]);

        long i = id2Index(link.id1);
        long j = id2Index(link.id2);
        long m = node[i].numEdge;
        long n = node[j].numEdge;

		node[i].edge.resize(m + 1);
		node[j].edge.resize(n + 1);
		

		if (link.overlapFlag) {
			node[i].edge[m].state = DBG_OVERLAP;
			node[j].edge[n].state = DBG_OVERLAP;
			node[i].edge[m].numLink = numLink - 1;
			node[j].edge[n].numLink = numLink - 1;
		}
		else {
			node[i].edge[m].state = 0x0;
			node[j].edge[n].state = 0x0;
			node[i].edge[m].numLink = numLink;
			node[j].edge[n].numLink = numLink;
		}

		node[i].edge[m].score = score;
		node[j].edge[n].score = score;

		node[i].edge[m].breakdown.resize(node[i].numContig, 0);
        for (unsigned int k = 0; k < node[i].numContig; ++k) {
            fread(&(node[i].edge[m].breakdown[k]), sizeof(long), 1, graphLinkFP[this->targetLibraryIndex]);
        }
		node[j].edge[n].breakdown.resize(node[j].numContig, 0);
        for (unsigned int k = 0; k < node[j].numContig; ++k) {
            fread(&(node[j].edge[n].breakdown[k]), sizeof(long), 1, graphLinkFP[this->targetLibraryIndex]);
        }

        node[i].edge[m].length = link.gap;
        node[j].edge[n].length = link.gap;
        node[i].edge[m].direction = static_cast<char>(link.id1 / (i + 1));
        node[j].edge[n].direction = static_cast<char>(-(link.id2) / (j + 1));
        if (link.id1 * link.id2 > 0) {
            node[i].edge[m].end = j + 1;
            node[j].edge[n].end = i + 1;
        }
        else {
            node[i].edge[m].end = -(j + 1);
            node[j].edge[n].end = -(i + 1);
        }
        ++(node[i].numEdge);
        ++(node[j].numEdge);
    }
    fclose(graphLinkFP[this->targetLibraryIndex]);
    graphLinkFP[this->targetLibraryIndex] = NULL;

    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        GraphNode &tmp = node[nodeID];
        std::sort(tmp.edge.begin(), tmp.edge.begin() + tmp.numEdge);
    }


	setOppositeBubbleNodeIDAndStateForEachNode();

	if (getMode() & LENGTH_CUTOFF_MODE) {
		if (getMode() & BUBBLE_AWARE_MODE)
			deleteEdgeFromShortNodeKeepingBubble(this->cutoffLength);
		else
			deleteEdgeFromShortNode(this->cutoffLength);
	}

	if (getMode() & PREVIOUS_DIVISION_AWARE_MODE)
		deleteEdgeFromDifferentPreviousParent(numThread);


//	if (getMode() & SECONDARY_BUBBLE_REMOVAL_MODE)
//		deleteEdgeFromSecondaryBubble();

	this->updateNodeClusterID();
}

void PairedDBG::makeGraphAllLibraries(const long numThread)
{
	unsigned numPairedEndLibrary = 0;
	if (allLibraryMT != NULL)
		numPairedEndLibrary = this->targetLibraryIndex + 1;
//		numPairedEndLibrary = (*allLibraryMT).size();

	unsigned numLongReadLibrary = 0;
	if (longReadLibraryMT != NULL)
		numLongReadLibrary = 1;

	if (graphLinkFP.empty())
		graphLinkFP.resize(numPairedEndLibrary + numLongReadLibrary, NULL);

	unsigned currentMode = getMode();
	for (unsigned libraryIndex = 0; libraryIndex < numPairedEndLibrary; ++libraryIndex) {
		this->setTargetLibraryIndex(libraryIndex);

//		if (getMode() & LENGTH_CUTOFF_MODE)
//			this->setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * (*allLibraryMT)[libraryIndex][0].getSDInsSize(), 0.1 * (*allLibraryMT)[libraryIndex][0].getAverageInsSize()));
//		else
//			this->setCutoffLength(0);

//		this->setTolerence(this->cutoffLength / 2);

		if (libraryIndex > 0 && (getMode() & OVERLAP_MODE))
			unsetMode(OVERLAP_MODE);

		calcLink(libraryIndex, 1, numThread);
	}

	for (unsigned libraryIndex = numPairedEndLibrary; libraryIndex < numPairedEndLibrary + numLongReadLibrary; ++libraryIndex) {
		this->setTargetLibraryIndex(libraryIndex);
//		this->setTolerence(2 * (this->minOverlap + 1));
//		if (getMode() & LENGTH_CUTOFF_MODE)
//			this->setCutoffLength(this->cutoffLengthFactor * (*longReadLibraryMT)[0].getAverageInsSize());
//		else
//			this->setCutoffLength(0);
			
		if (libraryIndex > 0 && (getMode() & OVERLAP_MODE))
			unsetMode(OVERLAP_MODE);
		
		unsetMode(PAIRED_END_LINK_MODE);
		calcLink(libraryIndex, 1, numThread);
	}

	setMode(currentMode);

    cerr << "constructing paired de Bruijn graph using all libraries simultaneously" << endl;

    for (int i = 0; i < numNode; ++i) {
		node[i].edge.clear();
		node[i].numEdge = 0;
	}

	for (unsigned libraryIndex = 0; libraryIndex < numPairedEndLibrary + numLongReadLibrary; ++libraryIndex) {
		rewind(graphLinkFP[libraryIndex]);
		long numLink;
		long score;
		while (fread(&numLink, sizeof(long), 1, graphLinkFP[libraryIndex])) {
			fread(&score, sizeof(long), 1, graphLinkFP[libraryIndex]);
			GraphLinkWithFlag link;
			fread(&link, sizeof(GraphLinkWithFlag), 1, graphLinkFP[libraryIndex]);

			long i = id2Index(link.id1);
			long m = 0 ;
			for (; m < node[i].numEdge; ++m) {
				if (node[i].edge[m].end == sign(link.id1) * link.id2 && node[i].edge[m].direction == sign(link.id1))
					break;
			}
			if (m == node[i].numEdge) {
				node[i].edge.resize(m + 1);
				++(node[i].numEdge);
			}

			long j = id2Index(link.id2);
			long n = 0;
			for (; n < node[j].numEdge; ++n) {
				if (node[j].edge[n].end == sign(link.id2) * link.id1 && node[j].edge[n].direction == -sign(link.id2))
					break;
			}
			if (n == node[j].numEdge) {
				node[j].edge.resize(n + 1);
				++(node[j].numEdge);
			}
			

			if (link.overlapFlag) {
				node[i].edge[m].state = DBG_OVERLAP;
				node[j].edge[n].state = DBG_OVERLAP;
				node[i].edge[m].numLink += (numLink - 1);
				node[j].edge[n].numLink += (numLink - 1);
			}
			else {
				node[i].edge[m].numLink += numLink;
				node[j].edge[n].numLink += numLink;
			}

			node[i].edge[m].score = score;
			node[j].edge[n].score = score;

			long breakdownNum = 0;
			node[i].edge[m].breakdown.resize(node[i].numContig, 0);
			for (unsigned int k = 0; k < node[i].numContig; ++k) {
				fread(&breakdownNum, sizeof(long), 1, graphLinkFP[libraryIndex]);
				node[i].edge[m].breakdown[k] += breakdownNum;
			}
			node[j].edge[n].breakdown.resize(node[j].numContig, 0);
			for (unsigned int k = 0; k < node[j].numContig; ++k) {
				fread(&breakdownNum, sizeof(long), 1, graphLinkFP[libraryIndex]);
				node[j].edge[n].breakdown[k] += breakdownNum;
			}

			if (!(link.overlapFlag)) {
				node[i].edge[m].length += link.gap * numLink;
				node[j].edge[n].length += link.gap * numLink;
			}
			else {
				node[i].edge[m].length = link.gap;
				node[j].edge[n].length = link.gap;
			}

			node[i].edge[m].direction = sign(link.id1);
			node[j].edge[n].direction = -sign(link.id2);

			if (link.id1 * link.id2 > 0) {
				node[i].edge[m].end = j + 1;
				node[j].edge[n].end = i + 1;
			}
			else {
				node[i].edge[m].end = -(j + 1);
				node[j].edge[n].end = -(i + 1);
			}
		}
		fclose(graphLinkFP[libraryIndex]);
		graphLinkFP[libraryIndex] = NULL;
	}

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
        GraphNode &tmp = node[nodeIndex];
        std::sort(tmp.edge.begin(), tmp.edge.begin() + tmp.numEdge);
		
		for (long edgeIndex = 0; edgeIndex < node[nodeIndex].numEdge; ++edgeIndex) {
			if (!(node[nodeIndex].edge[edgeIndex].state & DBG_OVERLAP))
				node[nodeIndex].edge[edgeIndex].length /= node[nodeIndex].edge[edgeIndex].numLink;
		}
    }

	deleteThinEdgeCostantKeepingOverlap(this->minLink);

	if (this->mode & PAIRED_END_LINK_MODE && allLibraryMT != NULL)
		this->setTargetLibraryIndex(numPairedEndLibrary - 1);

	setOppositeBubbleNodeIDAndStateForEachNode();

	if (getMode() & LENGTH_CUTOFF_MODE) {
		if (getMode() & BUBBLE_AWARE_MODE)
			deleteEdgeFromShortNodeKeepingBubble(this->cutoffLength);
		else
			deleteEdgeFromShortNode(this->cutoffLength);
	}

//	if (getMode() & SECONDARY_BUBBLE_REMOVAL_MODE)
//		deleteEdgeFromSecondaryBubble();

	this->updateNodeClusterID();
}

void PairedDBG::resetGraph()
{
    this->numNode = this->numContig;
    node.clear();
    node.resize(this->numNode);

    contigPositionInScaffold.clear();
    contigPositionInScaffold.resize(this->numContig);

    numBubble.clear();
    numBubble.resize(this->numContig, 0);

    for (long i = 0; i < numNode; ++i) {
        node[i].state = 0;
        node[i].numContig = 1;
        node[i].contig.clear();
        node[i].contig.resize(1);
        node[i].contig[0].id = i + 1;
        node[i].length = node[i].contig[0].end = this->contig[i].length;
        contigPositionInScaffold[i].id = i + 1;
        contigPositionInScaffold[i].offset = 0;
    }

	for (unsigned long i = 0; i < this->contigState.size(); ++i) {
		contigState[i] = 0;
		oppositeBubbleContigID[i] = 0;
	}
}

void PairedDBG::getOverlappedBubbleNodeIndex(vector<long> &bubbleNodeIndex)
{
	vector<long> nodeIDBuffer;
	vector<vector<long> > sinkNodeIDBuffer(2);
	vector<char> bubbleNodeFlag(this->numNode, 0);

    for (long sourceNodeIndex = 0; sourceNodeIndex < numNode; ++sourceNodeIndex) {
		for (char sourceDirection = -1; sourceDirection <= 1; sourceDirection += 2) {
			getOverlappedNode(sourceNodeIndex, sourceDirection, nodeIDBuffer);
			if (nodeIDBuffer.size() != 2)
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				if (sourceDirection == 1) {
					if (nodeIDBuffer[i] >= 0)
						getOverlappedNode(nodeIDBuffer[i] - 1, 1, sinkNodeIDBuffer[i]);
					else
						getOverlappedNode(-(nodeIDBuffer[i]) - 1, -1, sinkNodeIDBuffer[i]);
				}
				else {
					if (nodeIDBuffer[i] >= 0)
						getOverlappedNode(nodeIDBuffer[i] - 1, -1, sinkNodeIDBuffer[i]);
					else
						getOverlappedNode(-(nodeIDBuffer[i]) - 1, 1, sinkNodeIDBuffer[i]);
				}

				if (sinkNodeIDBuffer[i].size() != 1)
					break;

				sinkNodeIDBuffer[i][0] *= sign(nodeIDBuffer[i]);
			}

			if (i == 2 && sinkNodeIDBuffer[0][0] == sinkNodeIDBuffer[1][0]) {
				bubbleNodeFlag[id2Index(nodeIDBuffer[0])] = 1;
				bubbleNodeFlag[id2Index(nodeIDBuffer[1])] = 1;
			}
		}
    }

	bubbleNodeIndex.clear();
    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (bubbleNodeFlag[nodeIndex]) {
			bubbleNodeIndex.push_back(nodeIndex);
		}
	}
}

void PairedDBG::getOverlappedBubbleNodePairID(vector<pair<long, long> > &bubbleNodePairID)
{
//	const double HOMO_COVERAGE_THRESHOLD = 2.0 * HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	vector<long> nodeIDBuffer;
	vector<vector<long> > sinkNodeIDBuffer(2);
	vector<char> bubbleNodeFlag(this->numNode, 0);

	bubbleNodePairID.clear();

    for (long sourceNodeIndex = 0; sourceNodeIndex < numNode; ++sourceNodeIndex) {
//		if (calcNodeCoverage(node[sourceNodeIndex]) > HOMO_COVERAGE_THRESHOLD)
//			continue;

		for (char sourceDirection = -1; sourceDirection <= 1; sourceDirection += 2) {
			getOverlappedNode(sourceNodeIndex, sourceDirection, nodeIDBuffer);
			if (nodeIDBuffer.size() != 2)
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				getOverlappedNode(id2Index(nodeIDBuffer[i]), sign(nodeIDBuffer[i]) * sourceDirection, sinkNodeIDBuffer[i]);
				if (sinkNodeIDBuffer[i].size() != 1)
					break;

				sinkNodeIDBuffer[i][0] *= sign(nodeIDBuffer[i]);
			}

			if (i != 2 || sinkNodeIDBuffer[0][0] != sinkNodeIDBuffer[1][0])
				continue;

			vector<long> nodeIDBufferForCheck;
			getOverlappedNode(id2Index(sinkNodeIDBuffer[0][0]), -(sign(sinkNodeIDBuffer[0][0]) * sourceDirection), nodeIDBufferForCheck);
			if (nodeIDBufferForCheck.size() != 2)
				continue;

			bubbleNodePairID.emplace_back(nodeIDBuffer[0], nodeIDBuffer[1]);
		}
    }
}

void PairedDBG::getGappedBubbleNodePairID(vector<pair<long, long> > &bubbleNodePairID)
{
//	const double HOMO_COVERAGE_THRESHOLD = 2.0 * HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	vector<NodeIDWithGap> nodeIDBuffer;
	vector<vector<NodeIDWithGap> > sinkNodeIDBuffer(2);

	bubbleNodePairID.clear();

    for (long sourceNodeIndex = 0; sourceNodeIndex < numNode; ++sourceNodeIndex) {
//		if (calcNodeCoverage(node[sourceNodeIndex]) > HOMO_COVERAGE_THRESHOLD)
//			continue;

		for (char sourceDirection = -1; sourceDirection <= 1; sourceDirection += 2) {
			getUniqueConflictingNode(sourceNodeIndex, sourceDirection, nodeIDBuffer);
			if (nodeIDBuffer.size() != 2 || abs(nodeIDBuffer[0].id) == abs(nodeIDBuffer[1].id))
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				getLinkedNode(id2Index(nodeIDBuffer[i].id), sign(nodeIDBuffer[i].id) * (-sourceDirection), sinkNodeIDBuffer[i]);
				if (sinkNodeIDBuffer[i].size() != 1)
					break;

				getLinkedNode(id2Index(nodeIDBuffer[i].id), sign(nodeIDBuffer[i].id) * sourceDirection, sinkNodeIDBuffer[i]);
				if (sinkNodeIDBuffer[i].size() != 1)
					break;

				sinkNodeIDBuffer[i][0].id *= sign(nodeIDBuffer[i].id);
			}

//			if (i == 2 && sinkNodeIDBuffer[0][0].id == sinkNodeIDBuffer[1][0].id && sourceNodeIndex != id2Index(sinkNodeIDBuffer[0][0].id)
//				&& calcNodeCoverage(node[id2Index(sinkNodeIDBuffer[0][0].id)]) <= HOMO_COVERAGE_THRESHOLD) {
			if (i == 2 && sinkNodeIDBuffer[0][0].id == sinkNodeIDBuffer[1][0].id && sourceNodeIndex != id2Index(sinkNodeIDBuffer[0][0].id)) {
				bubbleNodePairID.emplace_back(nodeIDBuffer[0].id, nodeIDBuffer[1].id);
			}
		}
    }
}

void PairedDBG::setOppositeContigIDGappedConflicting(const long numThread)
{
	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	if (this->oppositeBubbleContigID.empty())
		this->oppositeBubbleContigID.resize(this->numContig, 0);

	vector<pair<long, long> > bubbleNodePairID;
	getGappedConflictingNodePairID(bubbleNodePairID, this->heteroCoverage);

	for (auto itr = bubbleNodePairID.begin(); itr != bubbleNodePairID.end(); ++itr) {
		GraphNode &node1 = this->node[id2Index(itr->first)];
		GraphNode &node2 = this->node[id2Index(itr->second)];

		if (calcNodeCoverage(node1) > COVERAGE_THRESHOLD || calcNodeCoverage(node2) > COVERAGE_THRESHOLD)
//		if (calcNodeCoverage(node1) > COVERAGE_THRESHOLD && calcNodeCoverage(node2) > COVERAGE_THRESHOLD)
			continue;

		long contigID2 = sign(itr->second) * node2.contig[0].id;
		for (auto contigItr = node1.contig.begin(); contigItr != node1.contig.end(); ++contigItr) {
			long contigID1 = sign(itr->first) * contigItr->id;
			if (this->oppositeBubbleContigID[id2Index(contigID1)] == 0) {
				this->oppositeBubbleContigID[id2Index(contigID1)] = sign(contigID1) * contigID2;
			}
		}

		long contigID1 = sign(itr->first) * node1.contig[0].id;
		for (auto contigItr = node2.contig.begin(); contigItr != node2.contig.end(); ++contigItr) {
			long contigID2 = sign(itr->second) * contigItr->id;
			if (this->oppositeBubbleContigID[id2Index(contigID2)] == 0) {
				this->oppositeBubbleContigID[id2Index(contigID2)] = sign(contigID2) * contigID1;
			}
		}
	}
}

void PairedDBG::getGappedConflictingNodePairID(vector<pair<long, long> > &bubbleNodePairID, double heteroNodeCoverage)
{
	const double MAX_HETERO_COVERAGE = HETERO_COVERAGE_THRESHOLD_FACTOR * heteroNodeCoverage;
	const double MAX_HOMO_COVERAGE = 3.0 * heteroNodeCoverage;
	vector<NodeIDWithGap> nodeIDBuffer;

	bubbleNodePairID.clear();
    for (long sourceNodeIndex = 0; sourceNodeIndex < numNode; ++sourceNodeIndex) {
		if (calcNodeCoverage(this->node[sourceNodeIndex]) > MAX_HOMO_COVERAGE)
			continue;

		for (char sourceDirection = -1; sourceDirection <= 1; sourceDirection += 2) {
			getUniqueConflictingNode(sourceNodeIndex, sourceDirection, nodeIDBuffer);
			if (nodeIDBuffer.size() == 2 && calcNodeCoverage(this->node[id2Index(nodeIDBuffer[0].id)]) <= MAX_HETERO_COVERAGE && calcNodeCoverage(this->node[id2Index(nodeIDBuffer[1].id)]) <= MAX_HETERO_COVERAGE) {
				bubbleNodePairID.emplace_back(nodeIDBuffer[0].id, nodeIDBuffer[1].id);
			}
		}
    }
}

void PairedDBG::getOverlappedNode(const long sourceNodeIndex, const char targetDirection, vector<long> &nodeIDBuffer)
{
	nodeIDBuffer.clear();
	GraphNode &node = this->node[sourceNodeIndex];
	for (long edgeIndex = 0; edgeIndex < node.numEdge; ++edgeIndex) {
		GraphEdge &edge = node.edge[edgeIndex];
		if (edge.direction == targetDirection && ((!(this->mode & NON_DBG_OVERLAP_MODE) && (edge.state & DBG_OVERLAP)) || ((this->mode & NON_DBG_OVERLAP_MODE) && edge.length < 0)))
			nodeIDBuffer.push_back(edge.end);
	}
}

long PairedDBG::getNumEdgeOneDirection(const GraphNode &targetNode, const char targetDirection)
{
	long n = 0;
	for (auto itr = targetNode.edge.begin(); itr != targetNode.edge.end(); ++itr) {
		if (itr->direction == targetDirection)
			++n;
	}
	return n;
}

void PairedDBG::markHeteroNode(const double maxHeteroCoverageFactor)
{
	const double maxHeteroCoverage = maxHeteroCoverageFactor * this->heteroCoverage;

	for (auto it = this->node.begin(); it != this->node.end(); ++it) {
		if (calcNodeCoverage(*it) <= maxHeteroCoverage)
			it->state |= DBG_HETERO;
	}
}

void PairedDBG::markBubbleHeteroNode(const vector<long> &candidateNodeIndex, const double maxHeteroCoverageFactor)
{
//	const double maxHeteroCoverage = maxHeteroCoverageFactor * this->heteroCoverage;

	for (auto it = candidateNodeIndex.begin(); it != candidateNodeIndex.end(); ++it) {
//		if (calcNodeCoverage(this->node[*it]) <= maxHeteroCoverage)
			this->node[*it].state |= DBG_HETERO;
	}
}

void PairedDBG::calculateHeteroCoverage(const vector<long> &bubbleNodeIndex)
{
	const long MIN_NUM_BUBBLE = 10000;
	const double TRUNCATION_FACTOR = 2.0; 

	vector<char> bubbleNodeFlag(this->node.size(), 0);

	long numBubble = 0;
	for (auto itr = bubbleNodeIndex.begin(); itr != bubbleNodeIndex.end(); ++itr) {
		bubbleNodeFlag[*itr] = 1;
		++numBubble;
	}

	vector<std::pair<int, int> > buffer(2);
	for (unsigned long i = 0; i < this->node.size(); ++i) {
		if (this->node[i].length <= this->contigMaxK)
			continue;
		
		if (bubbleNodeFlag[i])
			buffer.push_back(std::make_pair(static_cast<int>(calcNodeCoverage(this->node[i]) + 0.5), static_cast<int>(this->node[i].length)));
		else if (numBubble < MIN_NUM_BUBBLE)
			buffer.push_back(std::make_pair(static_cast<int>(calcNodeCoverage(this->node[i])/2.0 + 0.5), static_cast<int>(this->node[i].length)));
	}

	if (!buffer.empty()) {
		long sum = 0;
		long totalLength = 0;
		double weightedMean  = 0.0;
		for (auto itr = buffer.begin(); itr != buffer.end(); ++itr) {
			sum +=  itr->first * itr->second;
			totalLength += itr->second;
		}
		weightedMean  = std::round(static_cast<double>(sum) / totalLength);

		sum = totalLength = 0;
		for (auto itr = buffer.begin(); itr != buffer.end(); ++itr) {
			if (weightedMean / TRUNCATION_FACTOR  <= itr->first && itr->first <= weightedMean * TRUNCATION_FACTOR) {
				sum +=  itr->first * itr->second;
				totalLength += itr->second;
			}
		}
		if (totalLength > 0)
			this->heteroCoverage = std::round(static_cast<double>(sum) / totalLength);
		else
			this->heteroCoverage = weightedMean;
	}
	else {
		this->heteroCoverage = 1.0;
	}
	this->averageCoverage = 2.0 * this->heteroCoverage;

/*
	double weightedMedian = 1.0;
	if (!buffer.empty()) {
		std::sort(buffer.begin(), buffer.end(), platanus::PairFirstLess());

		long totalLength = 0;
		for (auto itr = buffer.begin(); itr != buffer.end(); ++itr)
			totalLength += itr->second;

		long sumLength = 0;
		for (auto itr = buffer.begin(); itr != buffer.end(); ++itr) {
			sumLength += itr->second;
			if (sumLength >= totalLength / 2) {
				weightedMedian = itr->first;
				break;
			}
		}
	}
	this->heteroCoverage = weightedMedian;
	this->averageCoverage = 2.0 * this->heteroCoverage;
*/

	cerr << "ESTIMATED_HETERO_COVERAGE = "<< this->heteroCoverage << endl;
}

void PairedDBG::calculateHeteroAndAverageCoverageUnphase()
{
	const double TRUNCATION_FACTOR = 2.0; 

	vector<std::pair<int, int> > buffer(2);
	unsigned long i;
	for (i = 0; i < this->node.size() - this->numInputBubbleContig; ++i) {
		if (this->node[i].length > this->contigMaxK)
			buffer.push_back(std::make_pair(static_cast<int>(calcNodeCoverage(this->node[i]) + 0.5), static_cast<int>(this->node[i].length)));
//			buffer.push_back(std::make_pair(static_cast<int>(calcNodeCoverage(this->node[i])/2.0 + 0.5), static_cast<int>(this->node[i].length)));
	}
	for (; i < this->node.size(); ++i) {
		if (this->node[i].length > this->contigMaxK)
//			buffer.push_back(std::make_pair(static_cast<int>(calcNodeCoverage(this->node[i]) + 0.5), static_cast<int>(this->node[i].length)));
			buffer.push_back(std::make_pair(static_cast<int>(calcNodeCoverage(this->node[i])*2.0 + 0.5), static_cast<int>(this->node[i].length)));
	}

	if (!buffer.empty()) {
		long sum = 0;
		long totalLength = 0;
		double weightedMean  = 0.0;
		for (auto itr = buffer.begin(); itr != buffer.end(); ++itr) {
			sum +=  itr->first * itr->second;
			totalLength += itr->second;
		}
		weightedMean  = std::round(static_cast<double>(sum) / totalLength);

		sum = totalLength = 0;
		for (auto itr = buffer.begin(); itr != buffer.end(); ++itr) {
			if (weightedMean / TRUNCATION_FACTOR  <= itr->first && itr->first <= weightedMean * TRUNCATION_FACTOR) {
				sum +=  itr->first * itr->second;
				totalLength += itr->second;
			}
		}
		if (totalLength > 0)
			this->averageCoverage = std::round(static_cast<double>(sum) / totalLength);
		else
			this->averageCoverage = weightedMean;
	}
	else {
		this->averageCoverage = 1.0;
	}
	this->heteroCoverage = this->heteroCoverage / 2.0;

	cerr << "ESTIMATED_AERAGE_COVERAGE = "<< this->averageCoverage << endl;
}

void PairedDBG::extractDBGBubbleInformation()
{
	const double MAX_HETERO_BUBBLE_COVERAGE_FACTOR = 2.0;
	vector<long> bubbleNodeIndex;

	getOverlappedBubbleNodeIndex(bubbleNodeIndex);
	if (this->heteroCoverage <= 0.0) {
		calculateHeteroCoverage(bubbleNodeIndex);
		this->averageCoverage = 2.0 * this->heteroCoverage;
	}
	markBubbleHeteroNode(bubbleNodeIndex, MAX_HETERO_BUBBLE_COVERAGE_FACTOR);
}

long PairedDBG::crushSimpleDBGBubble()
{
    const double coverageThreshold = heteroCoverage * 3.0;

    vector<long> ids;
	vector<long> nodeIDBuffer;
	vector<vector<long> > sinkNodeIDBuffer(2);

	long numCrush = 0;
    for (long sourceNodeIndex = 0; sourceNodeIndex < numNode; ++sourceNodeIndex) {
		for (char sourceDirection = -1; sourceDirection <= 1; sourceDirection += 2) {
			getOverlappedNode(sourceNodeIndex, sourceDirection, nodeIDBuffer);
			if (nodeIDBuffer.size() != 2)
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				if (sourceDirection == 1) {
					if (nodeIDBuffer[i] >= 0)
						getOverlappedNode(nodeIDBuffer[i] - 1, 1, sinkNodeIDBuffer[i]);
					else
						getOverlappedNode(-(nodeIDBuffer[i]) - 1, -1, sinkNodeIDBuffer[i]);
				}
				else {
					if (nodeIDBuffer[i] >= 0)
						getOverlappedNode(nodeIDBuffer[i] - 1, -1, sinkNodeIDBuffer[i]);
					else
						getOverlappedNode(-(nodeIDBuffer[i]) - 1, 1, sinkNodeIDBuffer[i]);
				}

				if (sinkNodeIDBuffer[i].size() != 1)
					break;

				sinkNodeIDBuffer[i][0] *= sign(nodeIDBuffer[i]);
			}

			long index1 = id2Index(nodeIDBuffer[0]);
			long index2 = id2Index(nodeIDBuffer[1]);
			GraphNode *node1 = &(node[index1]);
			GraphNode *node2 = &(node[index2]);
			double coverage1 = this->calcNodeCoverage(*node1);
			double coverage2 = this->calcNodeCoverage(*node2);

			if (i != 2 ||
				sinkNodeIDBuffer[0][0] != sinkNodeIDBuffer[1][0] ||
				(node1->state & SC_DEL) ||
				(node2->state & SC_DEL) ||
				coverage1 + coverage2 > coverageThreshold) {

				continue;
			}

			getOverlappedNode(id2Index(sinkNodeIDBuffer[0][0]), -(sourceDirection * sign(sinkNodeIDBuffer[0][0])), nodeIDBuffer);
			if (nodeIDBuffer.size() != 2)
				continue;

			if (coverage1 > coverage2) {
				std::swap(node1, node2);
				std::swap(index1, index2);
			}

			node1->state |= SC_DEL;
			for (long n = 0; n < node1->numEdge; ++n) {
				ids.push_back(index1 + 1);
				ids.push_back(node1->edge[n].end);
			}

			writeNodeSeq(*node1, bubbleFP);
			writeNodeSeq(*node2, bubbleOpositeFP);

			++numCrush;
		}
    }

    this->deleteEdges(ids);

	return numCrush;
}

void PairedDBG::deleteNonOverlapHomoEdge()
{
    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing non-hetero-informative edges..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
			if (node1.edge[edgeID].state & DBG_OVERLAP)
				continue;

            GraphNode &node2 = this->node[id2Index(node1.edge[edgeID].end)];
			if ((node1.state & DBG_HETERO) && (node2.state & DBG_HETERO))
				continue;

			ids.push_back(nodeID + 1);
			ids.push_back(node1.edge[edgeID].end);
			++numDelete;
        }
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

long PairedDBG::joinUnambiguousNodePair(const long numThread)
{
	vector<vector<GraphPath> > pathBuffer(numThread);

    # pragma omp parallel for schedule(static, 1)
	for (long t = 0; t < numThread; ++t) {
		vector<long> nodeIDBuffer;
		for (long initialNodeIndex = t; initialNodeIndex < numNode; initialNodeIndex += numThread) {
			if (node[initialNodeIndex].state & SC_DEL)
				continue;

			for (long i = 0; i < 2; ++i) {
				char initialDirection = 2*i - 1;
				getOverlappedNode(initialNodeIndex, initialDirection, nodeIDBuffer);
				if (nodeIDBuffer.size() != 1)
					continue;

				long altNodeID = nodeIDBuffer[0];

				getOverlappedNode(id2Index(altNodeID), -(initialDirection * sign(altNodeID)), nodeIDBuffer);
				if (nodeIDBuffer.size() != 1)
					continue;

				pathBuffer[t].emplace_back(initialNodeIndex, 2);
				if (i == 0) {
					pathBuffer[t].back().nodeID[0] = altNodeID;
					pathBuffer[t].back().nodeID[1] = initialNodeIndex + 1;
				}
				else {
					pathBuffer[t].back().nodeID[0] = initialNodeIndex + 1;
					pathBuffer[t].back().nodeID[1] = altNodeID;
				}

				break;
			}
		}
	}

	vector<GraphPath> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), GraphPathSelfIDLess());

	long numJoin = remakeGraphAccordingToPath(mergedPathBuffer);

	return numJoin;
}

long PairedDBG::joinUnambiguousNodePath(const long numThread)
{
	vector<vector<GraphPath> > pathBuffer(numThread);

    # pragma omp parallel for schedule(static, 1)
	for (long t = 0; t < numThread; ++t) {
		vector<long> nodeIDBuffer;
		vector<char> visitFlag(this->numNode, 0);
		for (long initialNodeIndex = t; initialNodeIndex < numNode; initialNodeIndex += numThread) {
			if (node[initialNodeIndex].state & SC_DEL)
				continue;

			visitFlag[initialNodeIndex] = 1;
			std::array<vector<long>, 2> bothSideNodeID;
			for (long i = 0; i < 2; ++i) {
				long currentNodeID = initialNodeIndex + 1;
				char initialDirection = 2*i - 1;
				while (1) {
					getOverlappedNode(id2Index(currentNodeID), initialDirection * sign(currentNodeID), nodeIDBuffer);
					if (nodeIDBuffer.size() != 1)
						break;

					long nextNodeID = nodeIDBuffer[0] * sign(currentNodeID);
					getOverlappedNode(id2Index(nextNodeID), -(initialDirection * sign(nextNodeID)), nodeIDBuffer);
					if (nodeIDBuffer.size() != 1)
						break;

					if (visitFlag[id2Index(nextNodeID)])
						break;

					visitFlag[id2Index(nextNodeID)] = 1;
					bothSideNodeID[i].push_back(nextNodeID);
					currentNodeID = nextNodeID;
				}
			}
			pathBuffer[t].emplace_back(initialNodeIndex, 0);
			GraphPath &path = pathBuffer[t].back();
			for (unsigned long i = 0; i < bothSideNodeID[0].size(); ++i)
				path.nodeID.push_back(bothSideNodeID[0][bothSideNodeID[0].size() - i - 1]);
			path.nodeID.push_back(initialNodeIndex + 1);
			for (unsigned long i = 0; i < bothSideNodeID[1].size(); ++i)
				path.nodeID.push_back(bothSideNodeID[1][i]);
		}
	}

	vector<GraphPath> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), GraphPathSelfIDLess());

	long numJoin = remakeGraphAccordingToPath(mergedPathBuffer);

	return numJoin;
}

long PairedDBG::remakeGraphAccordingToPath(vector<GraphPath> &pathBuffer)
{
	long numNewNode = 0;
	long newContigPoolSize = 0;
	long numJoin = 0;

    FILE *newContigFP = platanus::makeTemporaryFile();

	for (auto pathItr = pathBuffer.begin(); pathItr != pathBuffer.end(); ++pathItr) {
		if ((node[id2Index(pathItr->nodeID.front())].state & SC_INC) || (node[id2Index(pathItr->nodeID.back())].state & SC_INC))
			continue;

		numJoin += pathItr->nodeID.size();
		newContigPoolSize += writeAndMarkOverlappedNodes(pathItr->nodeID, newContigFP);
		++numNewNode;
	}

	this->writeSingletonNode(numNewNode, newContigPoolSize, newContigFP);
    this->remake(numNewNode, newContigPoolSize, newContigFP);
    fclose(newContigFP);

	return numJoin;
}

long PairedDBG::remakeGraphAccordingToPathPair(vector<GraphPath> &pathBuffer)
{
	long numNewNode = 0;
	long newContigPoolSize = 0;
	long numJoin = 0;
	vector<long> nodeIDForWriting;
	vector<long> nodeIDForWriting2;

    FILE *newContigFP = platanus::makeTemporaryFile();

	for (auto pathItr = pathBuffer.begin(); pathItr != pathBuffer.end(); pathItr += 2) {
		nodeIDForWriting.clear();
		auto nodeIDItr = pathItr->nodeID.begin();
		for (; nodeIDItr != pathItr->nodeID.end(); ++nodeIDItr) {
			if (node[id2Index(*nodeIDItr)].state & SC_INC)
				break;
			else
				nodeIDForWriting.push_back(*nodeIDItr);
		}
		if (nodeIDItr != pathItr->nodeID.end())
			continue;

		nodeIDForWriting2.clear();
		auto pathItr2 = pathItr + 1;
		nodeIDItr = pathItr2->nodeID.begin();
		for (; nodeIDItr != pathItr2->nodeID.end(); ++nodeIDItr) {
			if (node[id2Index(*nodeIDItr)].state & SC_INC)
				break;
			else
				nodeIDForWriting2.push_back(*nodeIDItr);
		}
		if (nodeIDItr != pathItr2->nodeID.end())
			continue;

		numJoin += (pathItr->nodeID.size() + pathItr2->nodeID.size());
		newContigPoolSize += writeAndMarkOverlappedNodes(nodeIDForWriting, newContigFP);
		newContigPoolSize += writeAndMarkOverlappedNodes(nodeIDForWriting2, newContigFP);
		numNewNode += 2;
	}

	this->writeSingletonNode(numNewNode, newContigPoolSize, newContigFP);
    this->remake(numNewNode, newContigPoolSize, newContigFP);
    fclose(newContigFP);

	return numJoin;
}

bool PairedDBG::isRedundantCrossComponent(const std::array<std::array<long, 2>, 2 > &externalNodeID, const long centerNodeID)
{
	std::array<long, 5> buffer;

	buffer[0] = abs(externalNodeID[0][0]);
	buffer[1] = abs(externalNodeID[0][1]);
	buffer[2] = abs(externalNodeID[1][0]);
	buffer[3] = abs(externalNodeID[1][1]);
	buffer[4] = abs(centerNodeID);
	
	for (long i = 0; i < 4; ++i) {
		for (long j = i + 1; j < 5; ++j) {
			if (buffer[i] == buffer[j])
				return true;
		}
	}

	return false;
}

bool PairedDBG::isRedundantGappedCrossComponent(const std::array<std::array<NodeIDWithGap, 2>, 2 > &externalNodeID, const long centerNodeID)
{
	std::array<long, 5> buffer;

	buffer[0] = abs(externalNodeID[0][0].id);
	buffer[1] = abs(externalNodeID[0][1].id);
	buffer[2] = abs(externalNodeID[1][0].id);
	buffer[3] = abs(externalNodeID[1][1].id);
	buffer[4] = abs(centerNodeID);
	
	for (long i = 0; i < 4; ++i) {
		for (long j = i + 1; j < 5; ++j) {
			if (buffer[i] == buffer[j])
				return true;
		}
	}

	return false;
}

long PairedDBG::solveSimpleCrossStructure(const double linkRateThreshold, const CROSS_RESOLUTION_MODE resolutionMode, const long numThread)
{
//	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;
//	const double COVERAGE_THRESHOLD = 2 * HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	vector<vector<GraphPath> > pathBuffer(numThread);

    # pragma omp parallel for schedule(static, 1)
	for (long t = 0; t < numThread; ++t) {
		vector<long> nodeIDBuffer;
		vector<long> nodeIDBufferForCheck;
		std::array<std::array<long, 2>, 2 > externalNodeID;

		for (long centerNodeIndex = t; centerNodeIndex < numNode; centerNodeIndex += numThread) {
//			if (2 * COVERAGE_THRESHOLD < calcNodeCoverage(this->node[centerNodeIndex]) || this->node[centerNodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE))
//				continue;

			long i = 0;
			for (; i < 2; ++i) {
				char initialDirection = 2*i - 1;
				getOverlappedNode(centerNodeIndex, initialDirection, nodeIDBuffer);
				if (nodeIDBuffer.size() != 2)
					break;

				long j = 0;
				for (; j < 2; ++j) {
					externalNodeID[i][j] = nodeIDBuffer[j];
					getOverlappedNode(id2Index(externalNodeID[i][j]), sign(externalNodeID[i][j]), nodeIDBufferForCheck);
					if (nodeIDBufferForCheck.size() >= 2)
						break;
				}
				if (j < 2)
					break;
			}
			if (i < 2)
				continue;

			if (isRedundantCrossComponent(externalNodeID, centerNodeIndex + 1))
				continue;

//			if (COVERAGE_THRESHOLD < std::min(std::min(calcNodeCoverage(this->node[id2Index(externalNodeID[0][0])]), calcNodeCoverage(this->node[id2Index(externalNodeID[0][1])])), std::min(calcNodeCoverage(this->node[id2Index(externalNodeID[1][0])]), calcNodeCoverage(this->node[id2Index(externalNodeID[1][1])]))))
//				continue;

			std::array<long, 2> sumLinkForHaplotype;
			sumLinkForHaplotype.fill(0);
			if (resolutionMode == TAG) {
				if (node[centerNodeIndex].length > maxFragmentLengthOfTag)
					continue;

				for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
					for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
						sumLinkForHaplotype[leftIndex == rightIndex] += this->getCommonTagBetweenNodePair(externalNodeID[0][leftIndex], externalNodeID[1][rightIndex]);
					}
				}
			}
			if (resolutionMode == SCORE) {
				for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
					for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
						sumLinkForHaplotype[leftIndex == rightIndex] += this->getScoreFromIDPair(externalNodeID[0][leftIndex], externalNodeID[1][rightIndex]);
					}
				}
			}
			else {
				for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
					for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
						sumLinkForHaplotype[leftIndex == rightIndex] += this->getNumLinkFromIDPair(externalNodeID[0][leftIndex], externalNodeID[1][rightIndex]);
					}
				}
			}

			if ((resolutionMode == LINK || resolutionMode == TAG)  && std::max(sumLinkForHaplotype[0], sumLinkForHaplotype[1]) < this->minLink)
				continue;

			char crossFlag;
			if (linkRateThreshold * sumLinkForHaplotype[0] >= sumLinkForHaplotype[1])
				crossFlag = 1;
			else if (linkRateThreshold * sumLinkForHaplotype[1] >= sumLinkForHaplotype[0])
				crossFlag = 0;
			else
				continue;

			for (long j = 0; j < 2; ++j) {
				pathBuffer[t].emplace_back(centerNodeIndex, 3);
				pathBuffer[t].back().nodeID[0] = externalNodeID[0][j];
				pathBuffer[t].back().nodeID[1] = centerNodeIndex + 1;
				pathBuffer[t].back().nodeID[2] = externalNodeID[1][(j + crossFlag) % 2];
			}
		}
	}

	vector<GraphPath> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), GraphPathSelfIDLess());

	long numSolvedCross = remakeGraphAccordingToPathPair(mergedPathBuffer);

	return numSolvedCross;
}

long PairedDBG::solveSimpleGappedCrossStructure(const double linkRateThreshold, const CROSS_RESOLUTION_MODE resolutionMode, const long numThread)
{
//	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;
//	const double COVERAGE_THRESHOLD = 2 * HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	vector<vector<GraphPathGapped> > pathBuffer(numThread);

    # pragma omp parallel for schedule(static, 1)
	for (long t = 0; t < numThread; ++t) {
		vector<NodeIDWithGap> nodeIDBuffer;
		vector<NodeIDWithGap> nodeIDBufferForCheck;
		std::array<std::array<NodeIDWithGap, 2>, 2 > externalNodeID;

		for (long centerNodeIndex = t; centerNodeIndex < numNode; centerNodeIndex += numThread) {
//			if (2 * COVERAGE_THRESHOLD < calcNodeCoverage(this->node[centerNodeIndex]) || this->node[centerNodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE))
//				continue;

			long i = 0;
			for (; i < 2; ++i) {
				char initialDirection = 2*i - 1;
				getUniqueConflictingNode(centerNodeIndex, initialDirection, nodeIDBuffer);
				if (nodeIDBuffer.size() != 2) {
					break;
				}

				long j = 0;
				for (; j < 2; ++j) {
					externalNodeID[i][j] = nodeIDBuffer[j];
					getLinkedNode(id2Index(externalNodeID[i][j].id), sign(externalNodeID[i][j].id), nodeIDBufferForCheck);
					if (nodeIDBufferForCheck.size() >= 3)
						break;
				}
				if (j < 2)
					break;
			}
			if (i < 2)
				continue;

			if (isRedundantGappedCrossComponent(externalNodeID, centerNodeIndex + 1))
				continue;

//			if (COVERAGE_THRESHOLD < std::min(std::min(calcNodeCoverage(this->node[id2Index(externalNodeID[0][0].id)]), calcNodeCoverage(this->node[id2Index(externalNodeID[0][1].id)])), std::min(calcNodeCoverage(this->node[id2Index(externalNodeID[1][0].id)]), calcNodeCoverage(this->node[id2Index(externalNodeID[1][1].id)]))))
//				continue;

			std::array<long, 2> sumLinkForHaplotype;
			sumLinkForHaplotype.fill(0);
			if (resolutionMode == TAG) {
				if (node[centerNodeIndex].length > maxFragmentLengthOfTag)
					continue;

				for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
					for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
						sumLinkForHaplotype[leftIndex == rightIndex] += this->getCommonTagBetweenNodePair(externalNodeID[0][leftIndex].id, externalNodeID[1][rightIndex].id);
					}
				}
			}
			if (resolutionMode == SCORE) {
				for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
					for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
						sumLinkForHaplotype[leftIndex == rightIndex] += this->getScoreFromIDPair(externalNodeID[0][leftIndex].id, externalNodeID[1][rightIndex].id);
					}
				}
			}
			else {
				for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
					for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
						sumLinkForHaplotype[leftIndex == rightIndex] += this->getNumLinkFromIDPair(externalNodeID[0][leftIndex].id, externalNodeID[1][rightIndex].id);
					}
				}
			}

			if ((resolutionMode == LINK || resolutionMode == TAG) && std::max(sumLinkForHaplotype[0], sumLinkForHaplotype[1]) < this->minLink)
				continue;

			char crossFlag;
			if (linkRateThreshold * sumLinkForHaplotype[0] >= sumLinkForHaplotype[1])
				crossFlag = 1;
			else if (linkRateThreshold * sumLinkForHaplotype[1] >= sumLinkForHaplotype[0])
				crossFlag = 0;
			else
				continue;
			
			for (long j = 0; j < 2; ++j) {
				pathBuffer[t].emplace_back(centerNodeIndex, 3);
				pathBuffer[t].back().node[0].id = externalNodeID[0][j].id;
				pathBuffer[t].back().node[0].gap = 0;
				pathBuffer[t].back().node[1].id = centerNodeIndex + 1;
				pathBuffer[t].back().node[1].gap = externalNodeID[0][j].gap;
				pathBuffer[t].back().node[2] = externalNodeID[1][(j + crossFlag) % 2];
			}
		}
	}

	vector<GraphPathGapped> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), GraphPathGappedSelfIDLess());

	long numSolvedCross = remakeGraphAccordingToGappedPathPair(mergedPathBuffer);

	return numSolvedCross;
}

long PairedDBG::writeAndMarkOverlappedNodes(const vector<long> &nodeID, FILE *storeFP)
{
	long numContig = 0;

	for (auto it = nodeID.begin(); it != nodeID.end(); ++it)
		numContig += node[id2Index(*it)].numContig;

	fwrite(&numContig, sizeof(long), 1, storeFP);

	long start = 0;
	for (auto it = nodeID.begin(); it != nodeID.end(); ++it) {
		long index;
		if (*it > 0) {
			index = *it - 1;
			node[index].state |= SC_INC;
			for (long i = 0; i < node[index].numContig; ++i) {
				ScaffoldPart scaffoldPartForWriting(node[index].contig[i].id, start + node[index].contig[i].start, start + node[index].contig[i].end);
				fwrite(&scaffoldPartForWriting, sizeof(ScaffoldPart), 1, storeFP);
			}
		}
		else {
			index = -(*it) - 1;
			node[index].state |= SC_INC;
			for (long i = node[index].numContig - 1; i >= 0; --i) {
				ScaffoldPart scaffoldPartForWriting(-(node[index].contig[i].id), start + node[index].length - node[index].contig[i].end, start + node[index].length - node[index].contig[i].start);
				fwrite(&scaffoldPartForWriting, sizeof(ScaffoldPart), 1, storeFP);
			}
		}

		start += node[index].length - this->minOverlap;
	}

	return numContig;
}

void PairedDBG::cutAndPrintSeq(const long minSeqLength, const unsigned long long readLength, const string &outFilename, const string &componentFilename)
{
    std::ofstream out(outFilename);
//    std::ofstream com(componentFilename);
    vector<long> maxLeftOverlap(numContig, 0);
    vector<long> maxRightOverlap(numContig, 0);

	vector<string> nodeName;
	long numOutputNode = 0;
    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node = this->node[nodeID];
		vector<char> nodeSeq;
		this->node2seq(node, nodeSeq);

		if (nodeSeq.size() < static_cast<unsigned long>(minSeqLength))
			continue;

		++numOutputNode;

		std::ostringstream oss;
		oss << ">seq" << numOutputNode << "_"
			<< "len" << nodeSeq.size() << "_"
			<< "cov" << calcNodeCoverage(node) << "_"
			<< "read" << readLength << "_"
			<< "maxK" << contigMaxK << endl;

		out << oss.str();

		unsigned long i;
		for (i = 0; i < nodeSeq.size(); ++i) {
			out << platanus::Bin2Char(nodeSeq[i]);
            if ((i + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                out.put('\n');
		}
        if (i % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
            out.put('\n');
    }

    out.close();
//    com.close();

    this->minOverlap = minOverlap;
}

long PairedDBG::solveUniquePathBetweenLinkedNodePair(const long numThread)
{
	vector<vector<GraphPath> > pathBuffer(numThread);
    # pragma omp parallel for schedule(static, 1)
	for (long t = 0; t < numThread; ++t) {
		for (long startNodeIndex = t; startNodeIndex < numNode; startNodeIndex += numThread) {
			GraphNode &startNode = this->node[startNodeIndex];
			for (long edgeIndex = 0; edgeIndex < startNode.numEdge; ++edgeIndex) {
				GraphEdge &edge = startNode.edge[edgeIndex];
				if (edge.state & DBG_OVERLAP)
					continue;

				searchUniqueOverlapPathGuidedByEdge(startNodeIndex + 1, edge, pathBuffer[t]); 
			}
		}
	}

	vector<GraphPath> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), GraphPathSelfIDLess());

	long numSolvedPath = remakeGraphAccordingToPath(mergedPathBuffer);
	return numSolvedPath;
}

void PairedDBG::searchUniqueOverlapPathGuidedByEdge(long startNodeID, const GraphEdge &guideEdge, vector<GraphPath> &pathBuffer)
{
	std::unordered_map<long, NodeInfoForPathSearch> nodeInfo;

	long endNodeID = guideEdge.end * sign(startNodeID);
	std::queue<long> nodeIDQueue;

	nodeIDQueue.push(startNodeID);
	nodeInfo[startNodeID] = NodeInfoForPathSearch(1, -(node[id2Index(startNodeID)].length), 0);
	while (nodeIDQueue.size() > 0) {
		long currentNodeID = nodeIDQueue.front();
		nodeIDQueue.pop();

		if (currentNodeID == endNodeID)
			continue;

		GraphNode &currentNode = node[id2Index(currentNodeID)];
		for (long edgeIndex = 0; edgeIndex < currentNode.numEdge; ++edgeIndex) {
			GraphEdge &edge = currentNode.edge[edgeIndex];
			if (!(edge.state & DBG_OVERLAP) || sign(currentNodeID) * edge.direction != guideEdge.direction)
				continue;

			long nextNodeID = edge.end * sign(currentNodeID);
			auto nodeInfoItr = nodeInfo.find(nextNodeID);
			long nextDistance = nodeInfo[currentNodeID].distance + currentNode.length + edge.length;
			if (nextDistance > guideEdge.length + this->tolerence)
				continue;

			if (nodeInfoItr == nodeInfo.end()) {
				nodeInfo[edge.end * sign(currentNodeID)] = NodeInfoForPathSearch(1, nextDistance, currentNodeID);
				nodeIDQueue.push(nextNodeID);
				continue;
			}

/*
			++(nodeInfoItr->second.numVisit);
			if (nextDistance < nodeInfoItr->second.distance) {
				nodeInfoItr->second.distance =  nextDistance;
				nodeInfoItr->second.preNodeID = currentNodeID;
				nodeIDQueue.push(nextNodeID);
			}
*/
		}
	}

	if (nodeInfo.size() <= 1)
		return;

	auto nodeInfoItr = nodeInfo.find(endNodeID);
	if (nodeInfoItr == nodeInfo.end())
		return;

	GraphPath path;
	while (nodeInfoItr->first != startNodeID) {
		 if (nodeInfoItr->second.numVisit != 1)
		 	return;
		path.nodeID.push_back(nodeInfoItr->first);
		nodeInfoItr = nodeInfo.find(nodeInfoItr->second.preNodeID);
	}
	path.nodeID.push_back(nodeInfoItr->first);
	if (guideEdge.direction > 0)
		std::reverse(path.nodeID.begin(), path.nodeID.end());

	pathBuffer.push_back(path);
	pathBuffer.back().sumLink = guideEdge.numLink;
	pathBuffer.back().selfID = startNodeID;
}

double PairedDBG::calcNodeCoverage(const GraphNode &node)
{
    unsigned long long sum = 0;
    unsigned long long num = 0;
    for (long i = 0; i < node.numContig; ++i) {
		unsigned long long contigIndex = id2Index(node.contig[i].id);
		if (contigPositionInScaffold[contigIndex].id != 0) {
			num += contig[contigIndex].length;
			sum += coverage[contigIndex] * contig[contigIndex].length;
		}
    }

    if (num == 0)
        return 0.0;
    else
        return static_cast<double>(sum) / num;

}

unsigned long long PairedDBG::crushHeteroBubble(const double averageCoverage)
{
    long numCrush = 0;
    const double homoCoverageThreshold = 2.0 * HETERO_COVERAGE_THRESHOLD_FACTOR * heteroCoverage;
    const double heteroCoverageThreshold = HETERO_COVERAGE_THRESHOLD_FACTOR * heteroCoverage;
    vector<long> ids;
    vector<GraphLayout> work, layout1, layout2;
    this->detectRepeat(averageCoverage);

    if (bubbleThreshold == 0.0) return 0;

	if (bubbleFP == NULL) 
		bubbleFP = platanus::makeTemporaryFile();

    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {
                const GraphEdge &edge1 = node[nodeID].edge[edgeID1];
                const GraphEdge &edge2 = node[nodeID].edge[edgeID2];

                work.clear();
                layout1.clear();
                layout2.clear();

                if (edge1.direction * edge2.direction < 0) continue;
                GraphNode *node1 = &(node[id2Index(edge1.end)]);
                if ((node1->state & SC_DEL) || edge1.length + node1->length <= edge2.length) continue;
                GraphNode *node2 = &(node[id2Index(edge2.end)]);
                if ((node2->state & SC_DEL) || edge2.length + node2->length <= edge1.length) continue;
                if (node1->isHomo && node2->isHomo) continue;
                long edgeEnd1, edgeEnd2;
                if (edge1.isForward()) {
                    edgeEnd1 = edge1.end;
                    edgeEnd2 = edge2.end;
                } else {
                    edgeEnd1 = edge2.end;
                    edgeEnd2 = edge1.end;
                }

                if (abs(edge1.length + node1->length - edge2.length) <= tolerence + this->getScaffoldOverlap(edgeEnd1, edgeEnd2)
                    || abs(edge2.length + node2->length - edge1.length) <= tolerence + this->getScaffoldOverlap(edgeEnd2, edgeEnd1)) continue;

                this->layoutNodes(node1, layout1, work);
                this->layoutNodes(node2, layout2, work);
                long rightEdge = std::min(layout1.size(), layout2.size());
                long leftEdge = rightEdge;
                long layoutID = 0;

                for (; layoutID < leftEdge; ++layoutID) {
                    if (layout1[layoutID].id != layout2[layoutID].id)
                        break;
                }
                if (layoutID == 0 || layoutID == leftEdge) continue;
                leftEdge = layoutID - 1;

//                if (calcNodeCoverage(node[abs(layout1[leftEdge].id) - 1]) >= homoCoverageThreshold) continue;
                for (layoutID = 1; layoutID <= rightEdge; ++layoutID) {
                    if (layout1[layout1.size() - layoutID].id != layout2[layout2.size() - layoutID].id)
                        break;
                }
                if (layoutID == 1) continue;
                rightEdge = layoutID - 1;

				if (abs(layout1[leftEdge].id) == abs(layout1[layout1.size() - rightEdge].id)) continue;

                if (calcNodeCoverage(node[id2Index(layout1[layout1.size() - rightEdge].id)]) > homoCoverageThreshold) continue;
                double coverage1 = this->layoutAverageCoverage(layout1, leftEdge + 1, layout1.size() - rightEdge - leftEdge - 1);
                double coverage2 = this->layoutAverageCoverage(layout2, leftEdge + 1, layout2.size() - rightEdge - leftEdge - 1);
                const vector<GraphLayout> &layoutRef = coverage1 < coverage2 ? layout1 : layout2;
                if (rightEdge + leftEdge + 1 >= static_cast<long>(layoutRef.size()) || coverage1 > heteroCoverageThreshold || coverage2 > heteroCoverageThreshold) continue;
                long numNodeInBubble = layoutRef.size() - rightEdge - leftEdge - 1;
                if (numNodeInBubble != 1 || coverage1 > heteroCoverageThreshold || coverage2 > heteroCoverageThreshold) continue;

                layoutID = leftEdge + 1;

				node1 = &(node[id2Index(layoutRef[layoutID].id)]);
				for (long n = 0; n < node1->numEdge; ++n) {
					ids.push_back(static_cast<long>(id2Index(layoutRef[layoutID].id)) + 1);
					ids.push_back(node1->edge[n].end);
				}
				for (long n = 0; n < node1->numContig; ++n)
					contigPositionInScaffold[id2Index(node1->contig[n].id)].id = 0;
				node1->state |= SC_DEL;

				if (coverage1 >= coverage2)
					std::swap(layout1, layout2);

				writeNodeSeq(node[id2Index(layout1[layoutID].id)], bubbleFP);
				writeNodeSeq(node[id2Index(layout2[layoutID].id)], bubbleOpositeFP);

                ++numCrush;
            }
        }
    }

    this->deleteEdges(ids);

    for (long i = 0; i < numNode; ++i)
        node[i].state &= ~SC_REP;

    cerr << "NUM_REMOVED_BUBBLES=" << numCrush << " (COVERAGE_THRESHOLD)" << endl;
    return numCrush;


}

long PairedDBG::deleteHeteroEdge(void)
{
    long numDelete = 0;
    const double homoCoverageThreshold = 2.0 * HETERO_COVERAGE_THRESHOLD_FACTOR * heteroCoverage;
    const double heteroCoverageThreshold = HETERO_COVERAGE_THRESHOLD_FACTOR * heteroCoverage;
    vector<long> ids;
    if (bubbleThreshold == 0.0) return 0;

    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {
                const GraphEdge &edge1 = node[nodeID].edge[edgeID1];
                const GraphEdge &edge2 = node[nodeID].edge[edgeID2];
                GraphNode *node1 = &(node[id2Index(edge1.end)]);
                GraphNode *node2 = &(node[id2Index(edge2.end)]);
                if (!((edge1.state & DBG_OVERLAP) && (edge2.state & DBG_OVERLAP) && edge1.direction * edge2.direction > 0) && !this->checkDeleteEdge(edge1, edge2, *node1, *node2)) continue;

                if (calcNodeCoverage(node[nodeID]) > homoCoverageThreshold) continue;

				long numExternalEdge1 = getNumEdgeOneDirection(*node1, edge1.direction * sign(edge1.end));
				long numExternalEdge2 = getNumEdgeOneDirection(*node2, edge2.direction * sign(edge2.end));
				if (numExternalEdge1 > 0 && numExternalEdge2 > 0)
					continue;

                double coverage1 = this->calcNodeCoverage(*node1);
                double coverage2 = this->calcNodeCoverage(*node2);
                long id = std::abs(edge1.end);
                if ((numExternalEdge1 > 0 && numExternalEdge2 == 0) || node1->length > node2->length) {
                    node1 = node2;
                    coverage1 = coverage2;
                    id = std::abs(edge2.end);
                }
                if (std::min(coverage1, coverage2) > heteroCoverageThreshold) continue;

                ++numDelete;
                node1->state |= SC_DEL;
                for (long edgeID = 0; edgeID < node1->numEdge; ++edgeID) {
                    ids.push_back(id);
                    ids.push_back(node1->edge[edgeID].end);
                }
                for (long contigID = 0; contigID < node1->numContig; ++contigID) {
                    contigPositionInScaffold[abs(node1->contig[contigID].id) - 1].id = 0;
                }
            }
        }
    }

    this->deleteEdges(ids);
    cerr << "NUM_DELETED_HETERO_EDGES = " << numDelete << endl;

    return numDelete;
}

void PairedDBG::loadResultSeq(const long minSeqLength, const unsigned long long readLength, const bool trimFlag, const bool clusterFlag)
{
	const long MIN_GAP_LENGTH = 10;
	const long MIN_OVERLAP_TO_JOIN = 32;

	resultSeq.clear();
	resultSeq.resize(node.size());

    vector<long> maxLeftOverlap(numContig, 0);
    vector<long> maxRightOverlap(numContig, 0);

    long defaultMinOverlap = this->minOverlap;
    this->minOverlap = MIN_OVERLAP_TO_JOIN;

    long maxNumContig = 0;
    long scaffoldLength = 0;

	if (clusterFlag)
		this->updateNodeClusterID();

	if (trimFlag) {
		if (!clusterFlag) {
			for (unsigned i = 0; i < TABLE_DIVID; ++i) {
				auto overlapIterator = overlapTable[i].begin();
				auto overlapEnd = overlapTable[i].end();
				for (; overlapIterator != overlapEnd; ++overlapIterator) {
					if (overlapIterator->second.id1 > 0) {
						if (overlapIterator->second.length > maxRightOverlap[overlapIterator->second.id1 - 1])
							maxRightOverlap[overlapIterator->second.id1 - 1] = overlapIterator->second.length;
					} else if (overlapIterator->second.length > maxLeftOverlap[-(overlapIterator->second.id1) - 1])
						maxLeftOverlap[-(overlapIterator->second.id1) - 1] = overlapIterator->second.length;

					if (overlapIterator->second.id2 > 0) {
						if (overlapIterator->second.length > maxLeftOverlap[overlapIterator->second.id2 - 1])
							maxLeftOverlap[overlapIterator->second.id2 - 1] = overlapIterator->second.length;
					} else if (overlapIterator->second.length > maxRightOverlap[-(overlapIterator->second.id2) - 1])
						maxRightOverlap[-(overlapIterator->second.id2) - 1] = overlapIterator->second.length;
				}
			}
		}
		else {
			for (unsigned i = 0; i < TABLE_DIVID; ++i) {
				auto overlapIterator = overlapTable[i].begin();
				auto overlapEnd = overlapTable[i].end();
				for (; overlapIterator != overlapEnd; ++overlapIterator) {
					long nodeID1 = contigPositionInScaffold[id2Index(overlapIterator->second.id1)];
					long nodeID2 = contigPositionInScaffold[id2Index(overlapIterator->second.id2)];

					long clusterID1 = 0;
					long clusterID2 = 0;
					if (nodeID1 != 0)
						clusterID1 = node[id2Index(nodeID1)].clusterID;
					if (nodeID2 != 0)
						clusterID2 = node[id2Index(nodeID2)].clusterID;

					if (clusterID1 == 0 || clusterID1 == clusterID2) {
						if (overlapIterator->second.id1 > 0) {
							if (overlapIterator->second.length > maxRightOverlap[overlapIterator->second.id1 - 1])
								maxRightOverlap[overlapIterator->second.id1 - 1] = overlapIterator->second.length;
						} else if (overlapIterator->second.length > maxLeftOverlap[-(overlapIterator->second.id1) - 1])
							maxLeftOverlap[-(overlapIterator->second.id1) - 1] = overlapIterator->second.length;
					}

					if (clusterID2 == 0 || clusterID2 == clusterID1) {
						if (overlapIterator->second.id2 > 0) {
							if (overlapIterator->second.length > maxLeftOverlap[overlapIterator->second.id2 - 1])
								maxLeftOverlap[overlapIterator->second.id2 - 1] = overlapIterator->second.length;
						} else if (overlapIterator->second.length > maxRightOverlap[-(overlapIterator->second.id2) - 1])
							maxRightOverlap[-(overlapIterator->second.id2) - 1] = overlapIterator->second.length;
					}
				}
			}
		}
	}

    for (long i = 0; i < numNode; ++i) {
        if (node[i].numContig > maxNumContig)
            maxNumContig = node[i].numContig;
    }
    vector<long> leftCut(maxNumContig, 0);
    vector<long> rightCut(maxNumContig, 0);
    vector<long> gap(maxNumContig, 0);

	long numOutputNode = 0;

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (node[nodeIndex].state & SC_DEL)
			continue;

        long contigID = 0;
        for (; contigID < node[nodeIndex].numContig; ++contigID) {
            if (contigPositionInScaffold[id2Index(node[nodeIndex].contig[contigID].id)].id != 0)
                break;
        }
        if (contigID == node[nodeIndex].numContig) continue;

        scaffoldLength = 0;

		if (trimFlag) {
			if (node[nodeIndex].contig[0].id > 0)
				leftCut[0] = maxLeftOverlap[node[nodeIndex].contig[0].id - 1] / 2;
			else 
				leftCut[0] = maxRightOverlap[-(node[nodeIndex].contig[0].id) - 1] / 2;
		}
		else {
			leftCut[0] = 0;
		}

        for (contigID = 0; contigID < node[nodeIndex].numContig - 1; ++contigID) {
            if (node[nodeIndex].contig[contigID].end > node[nodeIndex].contig[contigID + 1].start) {
				if (this->getOverlap(node[nodeIndex].contig[contigID].id, node[nodeIndex].contig[contigID + 1].id) > minOverlap) {
					rightCut[contigID] = 0;
					leftCut[contigID + 1] = this->getOverlap(node[nodeIndex].contig[contigID].id, node[nodeIndex].contig[contigID + 1].id);
					gap[contigID] = 0;
				}
				else {
					rightCut[contigID] = 0;
					leftCut[contigID + 1] = 0;
					gap[contigID] = MIN_GAP_LENGTH;
				}
            }
            else if (node[nodeIndex].contig[contigID + 1].start - node[nodeIndex].contig[contigID].end <= this->tolerence) {
				if (this->getOverlap(node[nodeIndex].contig[contigID].id, node[nodeIndex].contig[contigID + 1].id) > minOverlap) {
					rightCut[contigID] = 0;
					leftCut[contigID + 1] = this->getOverlap(node[nodeIndex].contig[contigID].id, node[nodeIndex].contig[contigID + 1].id);
					gap[contigID] = 0;
				}
				else {
					rightCut[contigID] = leftCut[contigID + 1] = 0;
					gap[contigID] = node[nodeIndex].contig[contigID + 1].start - node[nodeIndex].contig[contigID].end;
				}
			}
			else {
                rightCut[contigID] = leftCut[contigID + 1] = 0;
                gap[contigID] = node[nodeIndex].contig[contigID + 1].start - node[nodeIndex].contig[contigID].end;
            }

            scaffoldLength += (node[nodeIndex].contig[contigID].end - node[nodeIndex].contig[contigID].start) - leftCut[contigID] - rightCut[contigID] + gap[contigID];
        }

		if (trimFlag) {
			if (node[nodeIndex].contig[contigID].id > 0)
				rightCut[contigID] = maxRightOverlap[node[nodeIndex].contig[contigID].id - 1] / 2;
			else 
				rightCut[contigID] = maxLeftOverlap[-(node[nodeIndex].contig[contigID].id) - 1] / 2;
		}
		else {
			rightCut[contigID] = 0;
		}

        scaffoldLength += (node[nodeIndex].contig[contigID].end - node[nodeIndex].contig[contigID].start) - leftCut[contigID] - rightCut[contigID];
        gap[contigID] = 0;

        if (scaffoldLength < minSeqLength) continue;

		resultSeq[nodeIndex].seq.clear();
        for (contigID = 0; contigID < node[nodeIndex].numContig; ++contigID) {
            if (node[nodeIndex].contig[contigID].id > 0) {
                for (long seqID = leftCut[contigID]; seqID < contig[node[nodeIndex].contig[contigID].id - 1].length - rightCut[contigID]; ++seqID)
                    resultSeq[nodeIndex].seq.push_back(contig[node[nodeIndex].contig[contigID].id - 1].base[seqID]);
            }
			else {
                for (long seqID = contig[-(node[nodeIndex].contig[contigID].id) - 1].length - leftCut[contigID] - 1; seqID >= (long)rightCut[contigID]; --seqID) {
                    if (contig[-(node[nodeIndex].contig[contigID].id) - 1].base[seqID] != 4)
                        resultSeq[nodeIndex].seq.push_back(0x3 ^ contig[-(node[nodeIndex].contig[contigID].id) - 1].base[seqID]);
                    else
                        resultSeq[nodeIndex].seq.push_back(4);
                }
            }

            for (long k = 0; k < gap[contigID]; ++k)
				resultSeq[nodeIndex].seq.push_back(4);
        }

		++numOutputNode;

		std::ostringstream oss;
		oss << "scaffold" << numOutputNode << "_"
			<< "len" << resultSeq[nodeIndex].seq.size() << "_"
			<< "cov" << calcNodeCoverage(node[nodeIndex]) << "_"
			<< "read" << readLength << "_"
			<< "maxK" << contigMaxK;

		resultSeq[nodeIndex].name = oss.str();


		long currentPosition = 0;
		resultSeq[nodeIndex].component.clear();
        for (contigID = 0; contigID < node[nodeIndex].numContig; ++contigID) {

			long startPosition = currentPosition;
			long endPosition = startPosition + contig[id2Index(node[nodeIndex].contig[contigID].id)].length - leftCut[contigID] - rightCut[contigID];

			if (contigID > 0 && gap[contigID - 1] == 0)
				startPosition -= leftCut[contigID];

			if (contigID < node[nodeIndex].numContig - 1 && gap[contigID] == 0)
				endPosition += rightCut[contigID];

			std::ostringstream oss;
			oss	<< startPosition << '\t'
				<< endPosition << '\t'
				<< contigName[id2Index(node[nodeIndex].contig[contigID].id)] << '\t'
				<< 0 << '\t';

			if (node[nodeIndex].contig[contigID].id > 0)
				oss << '+';
			else
				oss << '-';

			resultSeq[nodeIndex].component += oss.str();

			currentPosition += contig[id2Index(node[nodeIndex].contig[contigID].id)].length - leftCut[contigID] - rightCut[contigID] + gap[contigID];
		}
    }

    this->minOverlap = defaultMinOverlap;
}

void PairedDBG::outputResultSeqWithBubble(const string filePrefix, const string &primaryBubbleSuffix, const string &secondaryBubbleSuffix, const string &nonBubbleHeteroSuffix, const string &nonBubbleOtherSuffix, const string &pairSuffix, const long minSeqLength, const long readLength)
{
    cerr << "writing scaffold and bubble files..." << endl;

    string outFilename;

    outFilename = filePrefix;
    outFilename += primaryBubbleSuffix;
    std::ofstream primaryBubbleOut(outFilename);

    outFilename = filePrefix;
    outFilename += secondaryBubbleSuffix;
    std::ofstream secondaryBubbleOut(outFilename);

    outFilename = filePrefix;
    outFilename += nonBubbleHeteroSuffix;
    std::ofstream nonBubbleHeteroOut(outFilename);

    outFilename = filePrefix;
    outFilename += nonBubbleOtherSuffix;
    std::ofstream nonBubbleOtherOut(outFilename);

    outFilename = filePrefix;
    outFilename += pairSuffix;
    std::ofstream pairOut(outFilename);

	vector<long> scaffoldIndexOrder(0);
	vector<long> bubbleIndexOrder(0);
	vector<char> pairFlag(numNode, false);
	vector<char> unprintFlag(numNode, false);

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (resultSeq[nodeIndex].seq.empty() || (node[nodeIndex].state & SC_DEL)) {
			unprintFlag[nodeIndex] = true;
			continue;
		}

		long oppositeNodeID = node[nodeIndex].oppositeBubbleNodeID;
		if (oppositeNodeID != 0) {
			if (node[id2Index(oppositeNodeID)].oppositeBubbleNodeID != 0 && id2Index(node[id2Index(oppositeNodeID)].oppositeBubbleNodeID) == nodeIndex)
				pairFlag[nodeIndex] = true;
		}
	}

	long numSeq = 0;

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (!(pairFlag[nodeIndex]) || (node[nodeIndex].state & DBG_SECONDARY_BUBBLE)) {
			continue;
		}

		long altNodeIndex = id2Index(node[nodeIndex].oppositeBubbleNodeID);
		if (resultSeq[nodeIndex].redundantFlag && resultSeq[altNodeIndex].redundantFlag)
			continue;

		if (node[nodeIndex].oppositeBubbleNodeID < 0)
			reverseComplement(resultSeq[altNodeIndex].seq);

        ++numSeq;

		primaryBubbleOut << ">primary_bubble" << numSeq << "_len" << resultSeq[nodeIndex].seq.size() << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeIndex]) + 0.5) << "_read" << readLength << "_maxK" << contigMaxK << endl;
		printBinSeqConstLineLength(resultSeq[nodeIndex].seq, primaryBubbleOut);

		secondaryBubbleOut << ">secondary_bubble" << numSeq << "_len" << resultSeq[altNodeIndex].seq.size() << "_cov" << static_cast<long>(calcNodeCoverage(node[altNodeIndex]) + 0.5) << "_read" << readLength << "_maxK" << contigMaxK << endl;
		printBinSeqConstLineLength(resultSeq[altNodeIndex].seq, secondaryBubbleOut);

		pairOut << "primary_bubble" << numSeq << "_len" << resultSeq[nodeIndex].seq.size() << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeIndex]) + 0.5) << "_read" << readLength << "_maxK" << contigMaxK << '\t';
		pairOut << "secondary_bubble" << numSeq << "_len" << resultSeq[altNodeIndex].seq.size() << "_cov" << static_cast<long>(calcNodeCoverage(node[altNodeIndex]) + 0.5) << "_read" << readLength << "_maxK" << contigMaxK << endl;
	}

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (unprintFlag[nodeIndex] || pairFlag[nodeIndex] || resultSeq[nodeIndex].redundantFlag)
			continue;

        ++numSeq;

		if ((node[nodeIndex].state & DBG_SECONDARY_BUBBLE)) {
			nonBubbleHeteroOut << ">non_bubble_hetero" << numSeq << "_len" << resultSeq[nodeIndex].seq.size() << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeIndex]) + 0.5) << "_read" << readLength << "_maxK" << contigMaxK << endl;
			printBinSeqConstLineLength(resultSeq[nodeIndex].seq, nonBubbleHeteroOut);
		}
		else {
			nonBubbleOtherOut << ">non_bubble_other" << numSeq << "_len" << resultSeq[nodeIndex].seq.size() << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeIndex]) + 0.5) << "_read" << readLength << "_maxK" << contigMaxK << endl;
			printBinSeqConstLineLength(resultSeq[nodeIndex].seq, nonBubbleOtherOut);
		}
	}

	primaryBubbleOut.close();
	secondaryBubbleOut.close();
	nonBubbleHeteroOut.close();
	nonBubbleOtherOut.close();
	pairOut.close();
}

void PairedDBG::outputResultSeqComponent(const string filePrefix, const string &fileSuffix)
{
    string outFilename;

    outFilename = filePrefix;
    outFilename += fileSuffix;
    std::ofstream out(outFilename);

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (!(resultSeq[nodeIndex].seq.empty() || (node[nodeIndex].state & SC_DEL) || resultSeq[nodeIndex].redundantFlag))
			out << resultSeq[nodeIndex].name << '\t' << resultSeq[nodeIndex].component << '\n';
	}

	out.close();
}

void PairedDBG::solveSimpleCrossStructureIterative(const CROSS_RESOLUTION_MODE resolutionMode, const bool clusterFlag, const long numThread)
{
	const unsigned long currentMinLink = this->minLink;

	long total = 0;
	long num;
	long iteration = 0;
	cerr << endl << "solving simple cross-structurs in a de Bruijn graph.." << endl;
	do {
		setMinLink(1);
		makeGraph(numThread);

		if (clusterFlag)
			deleteChimericEdge();
		else
			deleteInterOTUEdge();

		setMinLink(currentMinLink);
		if (resolutionMode == TAG)
			countMappedTagForEachScaffold(numThread);

		if (resolutionMode == SCORE)
			num = solveSimpleCrossStructure(CROSS_SCORE_RATE_THRESHOLD, resolutionMode, numThread);
		else
			num = solveSimpleCrossStructure(CROSS_LINK_RATE_THRESHOLD, resolutionMode, numThread);

		total += num;
		++iteration;
		cerr << "NUM_SOLVED_CROSSES = " << num << endl;
	} while (num > 0 && iteration < MAX_ITERATION_OF_CROSS_SOLUTION);
	cerr << "TOTAL_NUM_NODES_IN_SOLVED_CROSSES=" << total << endl << endl;
}

void PairedDBG::solveSimpleCrossStructureAllLibrariesIterative(const CROSS_RESOLUTION_MODE resolutionMode, const bool clusterFlag, const long numThread)
{
	const unsigned long currentMinLink = this->minLink;

	long total = 0;
	long num;
	long iteration = 0;
	cerr << endl << "solving simple cross-structurs in a de Bruijn graph.." << endl;
	do {
		setMinLink(1);
		makeGraphAllLibraries(numThread);

		if (clusterFlag)
			deleteChimericEdge();
		else
			deleteInterOTUEdge();

		setMinLink(currentMinLink);
		if (resolutionMode == TAG)
			countMappedTagForEachScaffold(numThread);
		num = solveSimpleCrossStructure(CROSS_LINK_RATE_THRESHOLD, resolutionMode, numThread);

		total += num;
		++iteration;
		cerr << "NUM_SOLVED_CROSSES = " << num << endl;
	} while (num > 0 && iteration < MAX_ITERATION_OF_CROSS_SOLUTION);
	cerr << "TOTAL_NUM_NODES_IN_SOLVED_CROSSES=" << total << endl << endl;
}

void PairedDBG::solveSimpleGappedCrossStructureAllLibrariesIterative(const CROSS_RESOLUTION_MODE resolutionMode, bool clusterFlag, const long numThread)
{
	const unsigned long currentMinLink = this->minLink;

	long total = 0;
	long num;
	long iteration = 0;
	cerr << endl << "solving simple cross-structurs in a scaffold graph using all libraries..." << endl;
	do {
		setMinLink(1);
		makeGraphAllLibraries(numThread);
		if (!((*allLibraryMT).empty()))
			deleteEdgeFromShortNode(2 * this->tolerenceFactor * (*allLibraryMT)[targetLibraryIndex][0].getSDInsSize());
		deleteLongEdge((*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize());

		if (clusterFlag)
			deleteChimericEdge();
		else
			deleteInterOTUEdge();

		setMinLink(currentMinLink);
		if (resolutionMode == TAG)
			countMappedTagForEachScaffold(numThread);

		if (resolutionMode == SCORE)
			num = solveSimpleGappedCrossStructure(CROSS_SCORE_RATE_THRESHOLD, resolutionMode, numThread);
		else
			num = solveSimpleGappedCrossStructure(CROSS_LINK_RATE_THRESHOLD, resolutionMode, numThread);

		total += num;
		++iteration;
		cerr << "NUM_SOLVED_CROSSES = " << num << endl;
	} while (num > 0 && iteration < MAX_ITERATION_OF_CROSS_SOLUTION);
	cerr << "TOTAL_NUM_NODES_IN_SOLVED_CROSSES=" << total << endl << endl;
}

void PairedDBG::solveSimpleGappedCrossStructureIterative(const CROSS_RESOLUTION_MODE resolutionMode, bool clusterFlag, const long numThread)
{
	long total = 0;
	long num;
	long iteration = 0;
	const unsigned long currentMinLink = this->minLink;
	cerr << endl << "solving simple cross-structurs in a scaffold graph..." << endl;

	do {
		setMinLink(1);
		makeGraph(numThread);
		setOppositeBubbleNodeIDForEachNode(numThread);
		setMinLink(currentMinLink);

		if (clusterFlag)
			deleteChimericEdge();
		else
			deleteInterOTUEdge();

		if (resolutionMode == TAG)
			countMappedTagForEachScaffold(numThread);

		if (resolutionMode == SCORE)
			num = solveSimpleGappedCrossStructure(CROSS_SCORE_RATE_THRESHOLD, resolutionMode, numThread);
		else
			num = solveSimpleGappedCrossStructure(CROSS_LINK_RATE_THRESHOLD, resolutionMode, numThread);

		total += num;
		++iteration;
		cerr << "NUM_SOLVED_CROSSES = " << num << endl;
	} while (num > 0 && iteration < MAX_ITERATION_OF_CROSS_SOLUTION);

	cerr << "TOTAL_NUM_NODES_IN_SOLVED_CROSSES=" << total << endl << endl;
}

void PairedDBG::joinUnambiguousNodePairIterative(const long numThread)
{
//	const double MAX_HETERO_COVERAGE_FACTOR = 1.5;

	unsigned currentMode = getMode();
	setMode(OVERLAP_MODE);

	long total = 0;
	long num;

	cerr << endl << "joining unambiguous pair of nodes in a de Bruijn graph.." << endl;
	do {
		makeGraph(numThread);
		extractDBGBubbleInformation();
//		markHeteroNode(MAX_HETERO_COVERAGE_FACTOR);
//		deleteNonOverlapHomoEdge();
		num = joinUnambiguousNodePair(numThread);
		total += num;
		cerr << "NUM_JOINED = " << num << endl;
	} while (num > 0);

	cerr << "TOTAL_NUM_NODES_IN_JOINED_PAIRS =" << total << endl << endl;
	setMode(currentMode);
}

void PairedDBG::joinUnambiguousNodePairGappedIterativeAllLibraries(const long numThread)
{
	long total = 0;
	long num;
	long iteration = 0;
	const unsigned long currentMinLink = this->minLink;
	cerr << endl << "joining unambiguous pair of nodes in a scaffold graph using all libraries.." << endl;

	do {
		setMinLink(1);
		makeGraphAllLibraries(numThread);
		setMinLink(currentMinLink);
//		if (!(libraryMT.empty()))
//			pairedDBG.deleteEdgeFromShortNode(2 * this->tolerenceFactor * (*allLibraryMT)[targetLibraryIndex][0].getSDInsSize());
//		pairedDBG.deleteLongEdge((*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize());

		num = joinUnambiguousNodePairGapped(numThread);
		total += num;
		++iteration;
		cerr << "NUM_JOINED = " << num << endl;
	} while (num > 0);

	cerr << "TOTAL_NUM_NODES_IN_JOINED_PAIRS =" << total << endl << endl;
}

void PairedDBG::solveUniquePathBetweenLinkedNodePairIterative(const long numThread)
{
	const double MAX_HETERO_COVERAGE_FACTOR = HETERO_COVERAGE_THRESHOLD_FACTOR;
	const double MAX_ITERATION = 2;

	long total = 0;
	long num;
	long numIteration = 0;
	cerr << endl << "solving unique paths guided by linked path in a de Bruijn graph.." << endl;
	do {
		makeGraph(numThread);
		extractDBGBubbleInformation();
		markHeteroNode(MAX_HETERO_COVERAGE_FACTOR);
		deleteNonOverlapHomoEdge();
		num = solveUniquePathBetweenLinkedNodePair(numThread);
		total += num;
		cerr << "NUM_JOINED = " << num << endl;

		++numIteration;
	} while (num > 0 && numIteration < MAX_ITERATION);
	cerr << "TOTAL_NUM_NODES_IN_SOLVED_PATHS =" << total << endl << endl;
}

void PairedDBG::solveUniquePathBetweenLinkedNodePairAllLibrariesIterative(const long numThread)
{
	const double MAX_HETERO_COVERAGE_FACTOR = HETERO_COVERAGE_THRESHOLD_FACTOR;
	const double MAX_ITERATION = 2;

	long total = 0;
	long num;
	long numIteration = 0;
	cerr << endl << "solving unique paths guided by linked path in a de Bruijn graph using all libraries.." << endl;
	do {
		makeGraphAllLibraries(numThread);
		extractDBGBubbleInformation();
		markHeteroNode(MAX_HETERO_COVERAGE_FACTOR);
		deleteNonOverlapHomoEdge();
		num = solveUniquePathBetweenLinkedNodePair(numThread);
		total += num;
		cerr << "NUM_JOINED = " << num << endl;

		++numIteration;
	} while (num > 0 && numIteration < MAX_ITERATION);
	cerr << "TOTAL_NUM_NODES_IN_SOLVED_PATHS =" << total << endl << endl;
}

void PairedDBG::setOppositeBubbleContigIDOverlapped(const long numThread)
{
	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	if (this->oppositeBubbleContigID.empty())
		this->oppositeBubbleContigID.resize(this->numContig, 0);

	vector<pair<long, long> > bubbleNodePairID;
	getOverlappedBubbleNodePairID(bubbleNodePairID);

	for (auto itr = bubbleNodePairID.begin(); itr != bubbleNodePairID.end(); ++itr) {
		GraphNode &node1 = this->node[id2Index(itr->first)];
		GraphNode &node2 = this->node[id2Index(itr->second)];

		if (calcNodeCoverage(node1) > COVERAGE_THRESHOLD || calcNodeCoverage(node2) > COVERAGE_THRESHOLD)
//		if (calcNodeCoverage(node1) > COVERAGE_THRESHOLD && calcNodeCoverage(node2) > COVERAGE_THRESHOLD)
			continue;

		long contigID2 = sign(itr->second) * node2.contig[0].id;
		for (auto contigItr = node2.contig.begin(); contigItr != node2.contig.end(); ++contigItr) {
			if (this->contigPositionInScaffold[id2Index(contigItr->id)].id != 0) {
				contigID2 = sign(itr->second) * contigItr->id;
				break;
			}
		}
		for (auto contigItr = node1.contig.begin(); contigItr != node1.contig.end(); ++contigItr) {
			long contigID1 = sign(itr->first) * contigItr->id;
			if (this->contigPositionInScaffold[id2Index(contigID1)] != 0 && this->oppositeBubbleContigID[id2Index(contigID1)] == 0) {
				this->oppositeBubbleContigID[id2Index(contigID1)] = sign(contigID1) * contigID2;
			}
		}

		long contigID1 = sign(itr->first) * node1.contig[0].id;
		for (auto contigItr = node1.contig.begin(); contigItr != node1.contig.end(); ++contigItr) {
			if (this->contigPositionInScaffold[id2Index(contigItr->id)].id != 0) {
				contigID1 = sign(itr->first) * contigItr->id;
				break;
			}
		}
		for (auto contigItr = node2.contig.begin(); contigItr != node2.contig.end(); ++contigItr) {
			long contigID2 = sign(itr->second) * contigItr->id;
			if (this->contigPositionInScaffold[id2Index(contigID2)] != 0 && this->oppositeBubbleContigID[id2Index(contigID2)] == 0) {
				this->oppositeBubbleContigID[id2Index(contigID2)] = sign(contigID2) * contigID1;
			}
		}
	}
}

void PairedDBG::setOppositeBubbleContigIDGapped(const long numThread)
{
	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	if (this->oppositeBubbleContigID.empty())
		this->oppositeBubbleContigID.resize(this->numContig, 0);

	vector<pair<long, long> > bubbleNodePairID;
	getGappedBubbleNodePairID(bubbleNodePairID);

	for (auto itr = bubbleNodePairID.begin(); itr != bubbleNodePairID.end(); ++itr) {
		GraphNode &node1 = this->node[id2Index(itr->first)];
		GraphNode &node2 = this->node[id2Index(itr->second)];

		if (calcNodeCoverage(node1) > COVERAGE_THRESHOLD || calcNodeCoverage(node2) > COVERAGE_THRESHOLD)
//		if (calcNodeCoverage(node1) > COVERAGE_THRESHOLD && calcNodeCoverage(node2) > COVERAGE_THRESHOLD)
			continue;

		long contigID2 = sign(itr->second) * node2.contig[0].id;
		for (auto contigItr = node2.contig.begin(); contigItr != node2.contig.end(); ++contigItr) {
			if (this->contigPositionInScaffold[id2Index(contigItr->id)].id != 0) {
				contigID2 = sign(itr->second) * contigItr->id;
				break;
			}
		}
		for (auto contigItr = node1.contig.begin(); contigItr != node1.contig.end(); ++contigItr) {
			long contigID1 = sign(itr->first) * contigItr->id;
			if (this->contigPositionInScaffold[id2Index(contigID1)] != 0 && this->oppositeBubbleContigID[id2Index(contigID1)] == 0) {
				this->oppositeBubbleContigID[id2Index(contigID1)] = sign(contigID1) * contigID2;
			}
		}

		long contigID1 = sign(itr->first) * node1.contig[0].id;
		for (auto contigItr = node1.contig.begin(); contigItr != node1.contig.end(); ++contigItr) {
			if (this->contigPositionInScaffold[id2Index(contigItr->id)].id != 0) {
				contigID1 = sign(itr->first) * contigItr->id;
				break;
			}
		}
		for (auto contigItr = node2.contig.begin(); contigItr != node2.contig.end(); ++contigItr) {
			long contigID2 = sign(itr->second) * contigItr->id;
			if (this->contigPositionInScaffold[id2Index(contigID2)] != 0 && this->oppositeBubbleContigID[id2Index(contigID2)] == 0) {
				this->oppositeBubbleContigID[id2Index(contigID2)] = sign(contigID2) * contigID1;
			}
		}
	}
}

long PairedDBG::divideNodeUsingBubbleContigPair(const long numThread)
{
//	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

    cerr << "dividing scaffolds to adjust bubble-breakpoints..." << endl;

    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();
	long totalNumDivision = 0;

	vector<long> numContigUsed(this->numContig, 0);

    omp_set_num_threads(numThread);
//    pragma omp parallel for schedule(dynamic) reduction(+: numNewNode, newContigPoolSize, totalNumDivision)
	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;
		vector<long> oppositeBubbleNodeID;
		setOppositeBubbleNodeID(oppositeBubbleNodeID, contigRef);

		pair<long, long> ends(0, oppositeBubbleNodeID.size());
		vector<pair<long, long> > endsStack(1, ends);

		while (endsStack.size() > 0) {
			ends = endsStack.back();
			endsStack.pop_back();

			pair<long, long> newEnds = fillMajorityIDRun(oppositeBubbleNodeID, this->node[nodeIndex], ends, 1.0);
			if (newEnds != ends) {
				endsStack.emplace_back(ends.first, newEnds.first);
				endsStack.emplace_back(newEnds.second, ends.second);
			}
		}

		vector<char> breakpointFlag(oppositeBubbleNodeID.size() + 1, 0);
		breakpointFlag.front() = 1;
		breakpointFlag.back() = 1;

		for (unsigned long i = 1; i < oppositeBubbleNodeID.size(); ++i) {
			if (oppositeBubbleNodeID[i - 1] != oppositeBubbleNodeID[i]) {
				breakpointFlag[i] = 1;
				++totalNumDivision;
			}
		}

		long i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;

			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
	//			pragma omp critical
	//			{
					fwrite(&k, sizeof(long), 1, scaffoldFP);
					for (k = j; k < i; ++k) {
						fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
						++numContigUsed[id2Index(contigRef[k].id)];
					}
	//			}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);

	return totalNumDivision;
}

long PairedDBG::divideNodeUsingBubbleContigPairStrandAware(const long numThread)
{
    cerr << "dividing scaffolds to adjust bubble-breakpoints distiguishing strands ..." << endl;

    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();
	long totalNumDivision = 0;

	vector<long> numContigUsed(this->numContig, 0);

    omp_set_num_threads(numThread);
//    pragma omp parallel for schedule(dynamic) reduction(+: numNewNode, newContigPoolSize, totalNumDivision)
	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;
		vector<long> oppositeBubbleNodeID;
		setOppositeBubbleNodeID(oppositeBubbleNodeID, contigRef);

		pair<long, long> ends(0, oppositeBubbleNodeID.size());
		vector<pair<long, long> > endsStack(1, ends);

		fillMajorityIDRun(oppositeBubbleNodeID, this->node[nodeIndex], ends, 1.0);
		flipOppositeBubbleNodeID(oppositeBubbleNodeID, contigRef);

		while (endsStack.size() > 0) {
			ends = endsStack.back();
			endsStack.pop_back();

			pair<long, long> newEnds = fillMajorityIDRun(oppositeBubbleNodeID, this->node[nodeIndex], ends, 1.0);
			if (newEnds != ends) {
				endsStack.emplace_back(ends.first, newEnds.first);
				endsStack.emplace_back(newEnds.second, ends.second);
			}
		}

		vector<char> breakpointFlag(oppositeBubbleNodeID.size() + 1, 0);
		breakpointFlag.front() = 1;
		breakpointFlag.back() = 1;

		for (unsigned long i = 0; i < oppositeBubbleNodeID.size(); ++i) {
			if (oppositeBubbleNodeID[i] == oppositeBubbleNodeID.back()) {
				breakpointFlag[i] = 1;
				if (i != 0)
					++totalNumDivision;
				break;
			}
		}

		long i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;

			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
	//			pragma omp critical
	//			{
					fwrite(&k, sizeof(long), 1, scaffoldFP);
					for (k = j; k < i; ++k) {
						fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
						++numContigUsed[id2Index(contigRef[k].id)];
					}
	//			}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);

	return totalNumDivision;
}

void PairedDBG::setOppositeBubbleNodeID(vector<long> &nodeIDVector, const vector<ScaffoldPart> &partVector)
{
	nodeIDVector.assign(partVector.size(), 0);
	for (auto partItr = partVector.begin(); partItr != partVector.end(); ++partItr) {
		if (this->oppositeBubbleContigID[id2Index(partItr->id)] == 0 || abs(this->contigPositionInScaffold[id2Index(partItr->id)].id) == abs(this->contigPositionInScaffold[id2Index(this->oppositeBubbleContigID[id2Index(partItr->id)])].id))
			nodeIDVector[partItr - partVector.begin()] = 0;
		else
			nodeIDVector[partItr - partVector.begin()] = abs(this->contigPositionInScaffold[id2Index(this->oppositeBubbleContigID[id2Index(partItr->id)])].id);
	}
}

void PairedDBG::setOppositeBubbleNodeIDStrandAware(vector<long> &nodeIDVector, const vector<ScaffoldPart> &partVector)
{
	nodeIDVector.assign(partVector.size(), 0);
	for (auto partItr = partVector.begin(); partItr != partVector.end(); ++partItr) {
		if (this->oppositeBubbleContigID[id2Index(partItr->id)] == 0 || abs(this->contigPositionInScaffold[id2Index(partItr->id)].id) == abs(this->contigPositionInScaffold[id2Index(this->oppositeBubbleContigID[id2Index(partItr->id)])].id))
			nodeIDVector[partItr - partVector.begin()] = 0;
		else
			nodeIDVector[partItr - partVector.begin()] = this->contigPositionInScaffold[id2Index(this->oppositeBubbleContigID[id2Index(partItr->id)])].id;
	}
}

void PairedDBG::flipOppositeBubbleNodeID(vector<long> &nodeIDVector, const vector<ScaffoldPart> &partVector)
{
	for (auto partItr = partVector.begin(); partItr != partVector.end(); ++partItr) {
		long oppositeNodeID;
		if (this->oppositeBubbleContigID[id2Index(partItr->id)] != 0)
			oppositeNodeID = this->contigPositionInScaffold[id2Index(this->oppositeBubbleContigID[id2Index(partItr->id)])].id;
		else
			oppositeNodeID = 0;

		if (oppositeNodeID == 0 || abs(oppositeNodeID) == abs(nodeIDVector[partItr - partVector.begin()]))
			nodeIDVector[partItr - partVector.begin()] = sign(partItr->id) * oppositeNodeID;
	}
}

double PairedDBG::calcNodeCoveragePartial(const GraphNode &node, const long start, const long end)
{
    unsigned long long sum = 0;
    unsigned long long num = 0;
    for (long i = start; i < end; ++i) {
		unsigned long long contigIndex = id2Index(node.contig[i].id);
		if (contigPositionInScaffold[contigIndex].id != 0) {
			num += contig[contigIndex].length;
			sum += coverage[contigIndex] * contig[contigIndex].length;
		}
    }

    if (num == 0)
        return 0.0;
    else
        return static_cast<double>(sum) / num;
}


long PairedDBG::maxLengthContigID(vector<long> &IDs, const long start, const long end)
{
	std::unordered_map<long, long> countMap;

	for (long i = start; i < end; ++i) {
		if (IDs[i] == 0)
			continue;

		auto countItr = countMap.find(IDs[i]);
		if (countItr != countMap.end())
			countItr->second += this->node[id2Index(IDs[i])].length;
		else
			countMap[IDs[i]] = this->node[id2Index(IDs[i])].length;
	}

	long maxCount = 0;
	long maxID = 0;
	for (auto countItr = countMap.begin(); countItr != countMap.end(); ++countItr) {
		if (maxCount < countItr->second) {
			maxID = countItr->first;
			maxCount = countItr->second;
		}
	}

	return maxID;
}

pair<long, long> PairedDBG::fillMajorityIDRun(vector<long> &IDs, const GraphNode &targetNode , const pair<long, long> &ends, const double scoreFactor)
{
	long maxID = maxLengthContigID(IDs, ends.first, ends.second);

	pair<long, long> newEnds(ends);
	if (maxID == 0)
		return newEnds;

	for (long i = ends.first; i < ends.second; ++i) {
		if (IDs[i] == maxID) {
			newEnds.first = i;
			break;
		}
	}
	for (long i = ends.second - 1; i >= ends.first; --i) {
		if (IDs[i] == maxID) {
			newEnds.second = i + 1;
			break;
		}
	}

	long score = 0;
	long maxScore = 0;
	long maxScoreI = newEnds.first;
	for (long i = newEnds.first; i < newEnds.second; ++i) {
		if (score > maxScore) {
			score = maxScore;
			maxScoreI = i;
		}
			
		if (IDs[i] != 0) {
			if (IDs[i] == maxID)
				score -= this->node[id2Index(IDs[i])].length;
			else
				score += this->node[id2Index(IDs[i])].length * scoreFactor;
		}
	}
	newEnds.first = maxScoreI;

	score = 0;
	maxScore = 0;
	maxScoreI = newEnds.second - 1;
	for (long i = newEnds.second - 1; i > newEnds.first; --i) {
		if (score > maxScore) {
			score = maxScore;
			maxScoreI = i;
		}
			
		if (IDs[i] != 0) {
			if (IDs[i] == maxID)
				score -= this->node[id2Index(IDs[i])].length;
			else
				score += this->node[id2Index(IDs[i])].length * scoreFactor;
		}
	}
	newEnds.second = maxScoreI + 1;

/*
	long maxIDCount = 0;
	for (long i = newEnds.first; i < newEnds.second; ++i) {
		if (IDs[i] == maxID)
			++maxIDCount;
	}
	if (maxIDCount * MIN_BUBBLE_COUNT_FACTOR >= (targetNode.contig[newEnds.second - 1].end - targetNode.contig[newEnds.first].start))
		std::fill(IDs.begin() + newEnds.first, IDs.begin() + newEnds.second, maxID);
*/

/*
	long regionLength = targetNode.contig[newEnds.second - 1].end - targetNode.contig[newEnds.first].start;
	double coverageThreshold =  this->heteroCoverage * std::max(1.25, HETERO_COVERAGE_THRESHOLD_FACTOR - (0.25 * 0.00001 * regionLength));
	if (calcNodeCoveragePartial(targetNode, newEnds.first, newEnds.second) <= coverageThreshold)
*/
		std::fill(IDs.begin() + newEnds.first, IDs.begin() + newEnds.second, maxID);

	return newEnds;
}

void PairedDBG::setOppositeBubbleNodeIDForEachNode(const long numThread)
{
	for (long nodeIndex = 0; nodeIndex < static_cast<long>(this->node.size()); ++nodeIndex) {
		this->node[nodeIndex].oppositeBubbleNodeID = 0;
		this->node[nodeIndex].state &= ~(DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE);
	}

	if (this->oppositeBubbleContigID.empty())
		return;

//	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

    cerr << "pairing bubble nodes..." << endl;

    omp_set_num_threads(numThread);
	# pragma omp parallel for schedule(static, 1)
	for (long nodeIndex = 0; nodeIndex < static_cast<long>(this->node.size()); ++nodeIndex) {
		vector<long> oppositeBubbleNodeID;
		setOppositeBubbleNodeID(oppositeBubbleNodeID, this->node[nodeIndex].contig);

		long oppositeNodeID = maxLengthContigID(oppositeBubbleNodeID, 0, oppositeBubbleNodeID.size());
//		if (oppositeNodeID == 0 || calcNodeCoverage(this->node[nodeIndex]) > COVERAGE_THRESHOLD || calcNodeCoverage(this->node[id2Index(oppositeNodeID)]) > COVERAGE_THRESHOLD)
		if (oppositeNodeID == 0 ||
			calcNodeCoverage(this->node[nodeIndex]) > this->heteroCoverage * std::max(1.25, HETERO_COVERAGE_THRESHOLD_FACTOR - (0.25 * 0.00001 * this->node[nodeIndex].length)) ||
			calcNodeCoverage(this->node[id2Index(oppositeNodeID)]) > this->heteroCoverage * std::max(1.25, HETERO_COVERAGE_THRESHOLD_FACTOR - (0.25 * 0.00001 * this->node[id2Index(oppositeNodeID)].length)))
			continue;

		if (id2Index(oppositeNodeID) == nodeIndex)
			continue;

		this->node[nodeIndex].oppositeBubbleNodeID = oppositeNodeID;
	}
}

void PairedDBG::setOppositeBubbleNodeIDAndStateForEachNode()
{
	for (long nodeIndex = 0; nodeIndex < static_cast<long>(this->node.size()); ++nodeIndex) {
		this->node[nodeIndex].oppositeBubbleNodeID = 0;
		this->node[nodeIndex].state &= ~(DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE);
	}

	for (unsigned long i = 0; i < this->contigState.size(); ++i)
		contigState[i] &= ~(DBG_CONTIG_PRIMARY_BUBBLE | DBG_CONTIG_SECONDARY_BUBBLE);

	if (this->oppositeBubbleContigID.empty())
		return;

//	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

    cerr << "pairing bubble nodes..." << endl;

	for (long nodeIndex = 0; nodeIndex < static_cast<long>(this->node.size()); ++nodeIndex) {
		vector<long> oppositeBubbleNodeID;
		setOppositeBubbleNodeID(oppositeBubbleNodeID, this->node[nodeIndex].contig);

		long oppositeNodeID = maxLengthContigID(oppositeBubbleNodeID, 0, oppositeBubbleNodeID.size());
//		if (oppositeNodeID == 0 || calcNodeCoverage(this->node[nodeIndex]) > COVERAGE_THRESHOLD || calcNodeCoverage(this->node[id2Index(oppositeNodeID)]) > COVERAGE_THRESHOLD)
		if (oppositeNodeID == 0 ||
			calcNodeCoverage(this->node[nodeIndex]) > this->heteroCoverage * std::max(1.25, HETERO_COVERAGE_THRESHOLD_FACTOR - (0.25 * 0.00001 * this->node[nodeIndex].length)) ||
			calcNodeCoverage(this->node[id2Index(oppositeNodeID)]) > this->heteroCoverage * std::max(1.25, HETERO_COVERAGE_THRESHOLD_FACTOR - (0.25 * 0.00001 * this->node[id2Index(oppositeNodeID)].length)))
			continue;

		long oppositeNodeIndex = id2Index(oppositeNodeID);
		if (oppositeNodeIndex == nodeIndex)
			continue;

		this->node[nodeIndex].oppositeBubbleNodeID = oppositeNodeID;

		long numEdgeDirection = getNumEdgeDirectionOfNode(this->node[nodeIndex]);
		long oppositeNumEdgeDirection = getNumEdgeDirectionOfNode(this->node[oppositeNodeIndex]);
		if (numEdgeDirection > oppositeNumEdgeDirection) {
			this->node[oppositeNodeIndex].state |= DBG_SECONDARY_BUBBLE;
		}
		else if (numEdgeDirection < oppositeNumEdgeDirection) {
			this->node[nodeIndex].state |= DBG_SECONDARY_BUBBLE;
		}
		else {
			long nonGapLength = getNonGapContigLengthOfNode(this->node[nodeIndex]);
			long oppositeNonGapLength = getNonGapContigLengthOfNode(this->node[oppositeNodeIndex]);
			if (nonGapLength > oppositeNonGapLength) {
				this->node[oppositeNodeIndex].state |= DBG_SECONDARY_BUBBLE;
			}
			else if (nonGapLength < oppositeNonGapLength) {
				this->node[nodeIndex].state |= DBG_SECONDARY_BUBBLE;
			}
			else {
				double nodeCoverage = calcNodeCoverage(this->node[nodeIndex]);
				double oppositeNodeCoverage = calcNodeCoverage(this->node[oppositeNodeIndex]);
				if (nodeCoverage > oppositeNodeCoverage) {
					this->node[oppositeNodeIndex].state |= DBG_SECONDARY_BUBBLE;
				}
				else if (nodeCoverage < oppositeNodeCoverage) {
					this->node[nodeIndex].state |= DBG_SECONDARY_BUBBLE;
				}
				else {
					if (nodeIndex < oppositeNodeIndex) {
						this->node[oppositeNodeIndex].state |= DBG_SECONDARY_BUBBLE;
					}
					else {
						this->node[nodeIndex].state |= DBG_SECONDARY_BUBBLE;
					}
				}
			}
		}

		if (this->node[nodeIndex].state & DBG_SECONDARY_BUBBLE)
			this->node[oppositeNodeIndex].state |= DBG_PRIMARY_BUBBLE;
		else
			this->node[nodeIndex].state |= DBG_PRIMARY_BUBBLE;
	}


/*
	vector<char> pairFlag(numNode, false);
    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		long oppositeNodeID = node[nodeIndex].oppositeBubbleNodeID;
		if (oppositeNodeID != 0) {
			if (node[id2Index(oppositeNodeID)].oppositeBubbleNodeID != 0 && id2Index(node[id2Index(oppositeNodeID)].oppositeBubbleNodeID) == nodeIndex)
				pairFlag[nodeIndex] = true;
		}
	}
    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (!(pairFlag[nodeIndex]) && (node[nodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)))
			node[nodeIndex].state &= ~(DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE);
	}
	pairFlag.clear();
*/


	for (unsigned long i = 0; i < this->contigState.size(); ++i) {
		if (contigPositionInScaffold[i].id == 0)
			continue;
		else if (node[id2Index(contigPositionInScaffold[i].id)].state & DBG_PRIMARY_BUBBLE)
			contigState[i] |= DBG_CONTIG_PRIMARY_BUBBLE;
		else if (node[id2Index(contigPositionInScaffold[i].id)].state & DBG_SECONDARY_BUBBLE)
			contigState[i] |= DBG_CONTIG_SECONDARY_BUBBLE;
	}
}

long PairedDBG::getNonGapContigLengthOfNode(const GraphNode &targetNode)
{
	if (targetNode.contig.empty())
		return 0;

	long gapLength = 0;
	for (unsigned i = 0; i < targetNode.contig.size() - 1; ++i)
		gapLength += targetNode.contig[i + 1].start - targetNode.contig[i].end;
	
	return targetNode.contig.back().end - gapLength;
}

long PairedDBG::getNumEdgeDirectionOfNode(const GraphNode &targetNode)
{
	long numLeft = 0;
	long numRight = 0;

	for (auto itr = targetNode.edge.begin(); itr != targetNode.edge.end(); ++itr) {
		if (itr->direction > 0)
			numLeft = 1;
		else
			numRight = 1;
	}

	return (numLeft + numRight);
}

long PairedDBG::deleteDifferentBubbleEdge(const long numThread)
{
    vector<long> ids;
    long numDelete = 0;

    omp_set_num_threads(numThread);

    #pragma omp parallel for schedule(dynamic) reduction(+: numDelete)
    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		const GraphNode &sourceNode = node[nodeIndex];
        for (long edgeID1 = 0; edgeID1 < sourceNode.numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeIndex].numEdge; ++edgeID2) {

                const GraphEdge &edge1 = (sourceNode.edge[edgeID1]);
                const GraphEdge &edge2 = (sourceNode.edge[edgeID2]);
                const GraphNode &node1 = (node[id2Index(edge1.end)]);
                const GraphNode &node2 = (node[id2Index(edge2.end)]);
                if (!this->checkDeleteEdge(edge1, edge2, node1, node2) || sourceNode.oppositeBubbleNodeID == 0) continue;

                if (sourceNode.oppositeBubbleNodeID != sign(edge1.end)*node1.oppositeBubbleNodeID && sourceNode.oppositeBubbleNodeID == sign(edge2.end)*node2.oppositeBubbleNodeID) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeIndex + 1);
                        ids.push_back(edge1.end);
                    }
                    ++numDelete;
                } else if (sourceNode.oppositeBubbleNodeID != sign(edge2.end)*node2.oppositeBubbleNodeID && sourceNode.oppositeBubbleNodeID == sign(edge2.end)*node1.oppositeBubbleNodeID) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeIndex + 1);
                        ids.push_back(edge2.end);
                    }
                    ++numDelete;
                }
            }
        }
    }
    this->deleteEdges(ids);
    return numDelete;
}

void PairedDBG::deleteDifferentBubbleEdgeIterative(const long numThread)
{
    cerr << "removing edges between nodes with different contig-bubble assignments..." << endl << endl;;

	this->setOppositeBubbleNodeIDForEachNode(numThread);

    long totalDelete = 0;
    long numDelete;
    do {
        numDelete = this->deleteDifferentBubbleEdge(numThread);
        totalDelete += numDelete;
        cerr << "NUM_REMOVED_EDGES =" << numDelete << endl;
    } while (numDelete > 0);

    cerr << "TOTAL_REMOVED_EDGES =" << totalDelete << endl << endl;;
}

void PairedDBG::deleteThinEdgeCostantKeepingOverlap(const long linkThreshold)
{
    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing thin edges (NUM_LINK < " <<  linkThreshold << ") ..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
            if (!(node1.edge[edgeID].state & DBG_OVERLAP) && node1.edge[edgeID].numLink < linkThreshold) {
				ids.push_back(nodeID + 1);
				ids.push_back(node1.edge[edgeID].end);
                ++numDelete;
            }
        }
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

long PairedDBG::deleteConflictingBubbleEdge(const long numThread)
{
	const double CROSS_LINK_RATE_THRESHOLD = 0.25;

    omp_set_num_threads(numThread);
	this->setOppositeBubbleNodeIDForEachNode(numThread);

    vector<long> ids;
    long numDelete = 0;

    cerr << "removing conflicting edges between bubbles..." << endl;

    #pragma omp parallel for schedule(dynamic) reduction(+: numDelete)
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
		if (node[nodeID].oppositeBubbleNodeID == 0)
			continue;

        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {

                const GraphEdge &edge1 = (node[nodeID].edge[edgeID1]);
                const GraphEdge &edge2 = (node[nodeID].edge[edgeID2]);
                const GraphNode &node1 = (node[abs(edge1.end) - 1]);
                const GraphNode &node2 = (node[abs(edge2.end) - 1]);
                if (!this->checkDeleteEdge(edge1, edge2, node1, node2))
					continue;

				if (node1.oppositeBubbleNodeID != 0 && node1.oppositeBubbleNodeID != sign(edge1.end) * edge2.end)
					continue;

                if (edge1.numLink < CROSS_LINK_RATE_THRESHOLD * edge2.numLink) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeID + 1);
                        ids.push_back(edge1.end);
                    }
                    ++numDelete;
                } else if (edge2.numLink < CROSS_LINK_RATE_THRESHOLD * edge1.numLink) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeID + 1);
                        ids.push_back(edge2.end);
                    }
                    ++numDelete;
                }
            }
        }
    }
    cerr << "TOTAL_NUM_DELETE=" << numDelete << endl << endl;
    this->deleteEdges(ids);
    return numDelete;
}

long PairedDBG::deleteSecondaryBubbleNodeAndEdge(const long numThread)
{
    omp_set_num_threads(numThread);
	this->setOppositeBubbleNodeIDAndStateForEachNode();

    cerr << "removing secondary bubbles from scaffold graph..." << endl;

    vector<long> ids;
    unsigned long long numDelete = 0;

    for (long nodeIndex = 0 ; nodeIndex < numNode; ++nodeIndex) {
		if (!(this->node[nodeIndex].state & DBG_SECONDARY_BUBBLE))
			continue;

        GraphNode &bubbleNode = this->node[nodeIndex];

		++numDelete;
		bubbleNode.state |= SC_DEL;
        for (long edgeIndex = 0; edgeIndex < bubbleNode.numEdge; ++edgeIndex) {
			ids.push_back(nodeIndex + 1);
			ids.push_back(bubbleNode.edge[edgeIndex].end);
        }
    }
    std::cerr << "TOTAL_NUM_DELETED_NODES =" << numDelete << std::endl;
    this->deleteEdges(ids);

	return numDelete;
}

long PairedDBG::deleteShortAndLowCoverageBranch(const long lengthThreshold, const double coverageThreshold, const long numThread)
{
    omp_set_num_threads(numThread);

    vector<long> ids;
    unsigned long long numDelete = 0;

    for (long nodeIndex = 0 ; nodeIndex < numNode; ++nodeIndex) {
        GraphNode &targetNode = this->node[nodeIndex];
		if (targetNode.length >= lengthThreshold || calcNodeCoverage(targetNode) >= coverageThreshold || getNumEdgeDirectionOfNode(this->node[nodeIndex]) >= 2)
			continue;

		++numDelete;
		targetNode.state |= SC_DEL;
        for (long edgeIndex = 0; edgeIndex < targetNode.numEdge; ++edgeIndex) {
			ids.push_back(nodeIndex + 1);
			ids.push_back(targetNode.edge[edgeIndex].end);
        }
    }
    this->deleteEdges(ids);

	return numDelete;
}

void PairedDBG::deleteShortAndLowCoverageBranchIterative(const long lengthThreshold, const double coverageThreshold, const long numThread)
{
	long total = 0;
	long num;
    cerr << "removing short and low-coverage branches..." << endl;
	do {
		joinUnambiguousNodePairIterative(numThread);

		makeGraph(numThread);
		num = deleteShortAndLowCoverageBranch(lengthThreshold, coverageThreshold, numThread);

		total += num;
		cerr << "NUM_SOLVED_CROSSES = " << num << endl;
	} while (num > 0);

	cerr << "TOTAL_NUM_REMOVED_NODES =" << total << endl << endl;
}

void PairedDBG::setBubbleJunctionContigIDOverlapped()
{
	if (this->joinedBubbleContigID.empty()) {
		this->joinedBubbleContigID.resize(this->numContig);
		for (auto itr = this->joinedBubbleContigID.begin(); itr != joinedBubbleContigID.end(); ++itr)
			(*itr)[0] = (*itr)[1] = 0;
	}

	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	if (this->contigState.empty())
		this->contigState.resize(this->numContig, 0);
	
	for (auto itr = this->contigState.begin(); itr != this->contigState.end(); ++itr)
		*itr &= ~DBG_CONTIG_BUBBLE_JUNCTION;

	vector<pair<long, long> > bubbleNodePairID;
	getOverlappedBubbleNodePairID(bubbleNodePairID);

	vector<char> bubbleFlag(this->node.size(), false);
	for (auto itr = bubbleNodePairID.begin(); itr != bubbleNodePairID.end(); ++itr) {
		if (calcNodeCoverage(this->node[id2Index(itr->first)]) <= COVERAGE_THRESHOLD && calcNodeCoverage(this->node[id2Index(itr->second)]) <= COVERAGE_THRESHOLD) {
//		if (calcNodeCoverage(this->node[id2Index(itr->first)]) <= COVERAGE_THRESHOLD || calcNodeCoverage(this->node[id2Index(itr->second)]) <= COVERAGE_THRESHOLD) {
			bubbleFlag[id2Index(itr->first)] = true;
			bubbleFlag[id2Index(itr->second)] = true;
		}
	}

	vector<long> nodeIDBuffer;

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		for (char direction = -1; direction <= 1; direction += 2) {
			getOverlappedNode(nodeIndex, direction, nodeIDBuffer);
			if (nodeIDBuffer.size() != 2)
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				if (!bubbleFlag[id2Index(nodeIDBuffer[i])])
					break;
			}
			if (i != 2)
				continue;

			if (direction > 0) {
				this->contigState[id2Index(this->node[nodeIndex].contig.back().id)] |= DBG_CONTIG_BUBBLE_JUNCTION;
				if (this->node[nodeIndex].contig.back().id > 0) {
					if (nodeIDBuffer[0] > 0)
						this->joinedBubbleContigID[id2Index(this->node[nodeIndex].contig.back().id)][1] = this->node[id2Index(nodeIDBuffer[0])].contig.front().id;
					else
						this->joinedBubbleContigID[id2Index(this->node[nodeIndex].contig.back().id)][1] = -(this->node[id2Index(nodeIDBuffer[0])].contig.back().id);
				}
				else {
					if (nodeIDBuffer[0] > 0)
						this->joinedBubbleContigID[id2Index(this->node[nodeIndex].contig.back().id)][0] = this->node[id2Index(nodeIDBuffer[0])].contig.front().id;
					else
						this->joinedBubbleContigID[id2Index(this->node[nodeIndex].contig.back().id)][0] = -(this->node[id2Index(nodeIDBuffer[0])].contig.back().id);
				}
			}
			else {
				this->contigState[id2Index(this->node[nodeIndex].contig.front().id)] |= DBG_CONTIG_BUBBLE_JUNCTION;
				if (this->node[nodeIndex].contig.front().id > 0) {
					if (nodeIDBuffer[0] > 0)
						this->joinedBubbleContigID[id2Index(this->node[nodeIndex].contig.front().id)][0] = this->node[id2Index(nodeIDBuffer[0])].contig.back().id;
					else
						this->joinedBubbleContigID[id2Index(this->node[nodeIndex].contig.front().id)][0] = -(this->node[id2Index(nodeIDBuffer[0])].contig.front().id);
				}
				else {
					if (nodeIDBuffer[0] > 0)
						this->joinedBubbleContigID[id2Index(this->node[nodeIndex].contig.front().id)][1] = this->node[id2Index(nodeIDBuffer[0])].contig.back().id;
					else
						this->joinedBubbleContigID[id2Index(this->node[nodeIndex].contig.front().id)][1] = -(this->node[id2Index(nodeIDBuffer[0])].contig.front().id);
				}
			}
		}
	}
}

void PairedDBG::markJunctionContigJoinedToBubble()
{
	vector<char> bubbleEdgeFlag(this->numContig, 0);

	for (auto itr = this->contigState.begin(); itr != this->contigState.end(); ++itr)
		*itr &= ~DBG_CONTIG_BUBBLE_JUNCTION;

	for (auto nodeIt = node.begin(); nodeIt != node.end(); ++nodeIt) {
		if (nodeIt->state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)) {
			bubbleEdgeFlag[id2Index(nodeIt->contig.front().id)] = 1;
			bubbleEdgeFlag[id2Index(nodeIt->contig.back().id)] = 1;
		}
	}

	for (auto nodeIt = node.begin(); nodeIt != node.end(); ++nodeIt) {
		for (auto contigIt = nodeIt->contig.begin(); contigIt != nodeIt->contig.end(); ++contigIt) {
			for (long i = 0; i < 2; ++i) {
				if (joinedBubbleContigID[id2Index(contigIt->id)][i] != 0 && bubbleEdgeFlag[id2Index(joinedBubbleContigID[id2Index(contigIt->id)][i])])
					this->contigState[id2Index(contigIt->id)] |= DBG_CONTIG_BUBBLE_JUNCTION;
			}
		}
	}
}


long PairedDBG::divideBubbleJunctionNode(const bool gapDivideFlag)
{
    cerr << "dividing scaffolds at bubble-junctions..." << endl;

	this->setOppositeBubbleNodeIDAndStateForEachNode();

	this->markJunctionContigJoinedToBubble();

	const long MIN_OVERLAP_TO_JOIN = 20;
    long defaultMinOverlap = this->minOverlap;
	if (gapDivideFlag) 
		this->minOverlap = MIN_OVERLAP_TO_JOIN;
	else
		this->minOverlap = this->contigMaxK - 1;

	long numDivision = 0;
    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		long i;
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.front() = 1;
		breakpointFlag.back() = 1;

		if (gapDivideFlag) {
			if (this->node[nodeIndex].oppositeBubbleNodeID == 0) {
				for (i = 1; i < static_cast<long>(contigRef.size()); ++i) {
					if ((((this->contigState[id2Index(contigRef[i - 1].id)] | this->contigState[id2Index(contigRef[i].id)]) & DBG_CONTIG_BUBBLE_JUNCTION) && this->getOverlap(contigRef[i - 1].id, contigRef[i].id) < this->contigMaxK) ||
						(contigRef[i].start > contigRef[i - 1].end || this->getOverlap(contigRef[i - 1].id, contigRef[i].id) < minOverlap)) {

						breakpointFlag[i] = 1;
						++numDivision;
					}
				}
			}
		}
		else {
			for (i = 1; i < static_cast<long>(contigRef.size()); ++i) {
				if (((this->contigState[id2Index(contigRef[i - 1].id)] | this->contigState[id2Index(contigRef[i].id)]) & DBG_CONTIG_BUBBLE_JUNCTION) 
					&& this->getOverlap(contigRef[i - 1].id, contigRef[i].id) < this->contigMaxK - 1) {

					breakpointFlag[i] = 1;
					++numDivision;
				}
			}
		}

		i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;


			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->minOverlap = defaultMinOverlap;
    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);

	return numDivision;
}

long PairedDBG::divideBubbleContigInNonHeteroNode()
{
    cerr << "dividing scaffolds at bubble-contig in non-hetero scaffolds ..." << endl;

	this->setOppositeBubbleNodeIDAndStateForEachNode();

	long numDivision = 0;
    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		long i;
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.front() = 1;
		breakpointFlag.back() = 1;

		if (this->node[nodeIndex].oppositeBubbleNodeID == 0 && contigRef.size() > 1) {
			for (i = 0; i < static_cast<long>(contigRef.size()); ++i) {
				if (this->oppositeBubbleContigID[id2Index(contigRef[i].id)] != 0) {
					breakpointFlag[i] = 1;
					breakpointFlag[i + 1] = 1;
					++numDivision;
				}
			}
		}

		i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;


			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);

	return numDivision;
}

void PairedDBG::divideGappedNode(const long minGapSize)
{
    cerr << "dividing scaffolds at gaps..." << endl;

	const long MIN_OVERLAP_TO_JOIN = 20;
    long defaultMinOverlap = this->minOverlap;
    this->minOverlap = MIN_OVERLAP_TO_JOIN;

    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		if (this->node[nodeIndex].state & SC_DEL)
			continue;

		long i;
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.back() = 1;

		for (i = 1; i < static_cast<long>(contigRef.size()); ++i) {
//			if (contigRef[i].start > contigRef[i - 1].end || this->getOverlap(contigRef[i - 1].id, contigRef[i].id) < minOverlap)
			if (contigRef[i].start - contigRef[i - 1].end > minGapSize)
				breakpointFlag[i] = 1;
		}

		i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;

			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->minOverlap = defaultMinOverlap;
    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::trimSparseEnd()
{
    cerr << "trimming sparse ends of scaffolds..." << endl;

    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		if (this->node[nodeIndex].state & SC_DEL)
			continue;

		long i;
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.back() = 1;


		if (contigRef.size() > 1) {
			if (contigRef[1].start - contigRef[0].end > contigRef[0].end - contigRef[0].start)
				breakpointFlag[1] = 1;

			if (contigRef[contigRef.size() - 1].start - contigRef[contigRef.size() - 2].end > contigRef[contigRef.size() - 1].end - contigRef[contigRef.size() - 1].start)
				breakpointFlag[contigRef.size() - 2] = 1;
		}


		i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;

			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::trimRepeatEnd()
{
    cerr << "trimming repeat-ends of scaffolds..." << endl;

	vector<long> currentNumContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		if (this->node[nodeIndex].state & SC_DEL)
			continue;

		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		for(long i = 0; i < contigRef.size(); ++i)
			++(currentNumContigUsed[id2Index(contigRef[i].id)]);
	}


    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<long> newNumContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		if (this->node[nodeIndex].state & SC_DEL)
			continue;

		long i;
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.back() = 1;

		for (long i = 0; i < contigRef.size(); ++i) {
			if (currentNumContigUsed[id2Index(contigRef[i].id)] < 2) {
				breakpointFlag[i] = 1;
				break;
			}
		}

		for (long i = contigRef.size() - 1; i > 0; --i) {
			if (currentNumContigUsed[id2Index(contigRef[i].id)] < 2) {
				breakpointFlag[i + 1] = 1;
				break;
			}
		}

		i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;

			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (newNumContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++newNumContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::copyAllNodes(PairedDBG &targetGraph)
{
    targetGraph.numNode = this->numNode;
    targetGraph.numContig = this->numContig;
    targetGraph.node= this->node;
    targetGraph.contigPositionInScaffold = this->contigPositionInScaffold;
	targetGraph.oppositeBubbleContigID = this->oppositeBubbleContigID;
	targetGraph.contigState = this->contigState;
	targetGraph.joinedBubbleContigID = this->joinedBubbleContigID;
}

void PairedDBG::divideUsingGuideGraph(PairedDBG &guideGraph, const bool bubbleDivideFlag)
{
	this->setOppositeBubbleNodeIDAndStateForEachNode();

	vector<long> numContigUsed(this->numContig, 0);

	vector<vector<char> > breakpointFlag(this->node.size());
	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		breakpointFlag[nodeIndex].assign(this->node[nodeIndex].contig.size() + 1, 0);
		breakpointFlag[nodeIndex].front() = 1;
		breakpointFlag[nodeIndex].back() = 1;
	}

	for (unsigned long nodeIndex = 0; nodeIndex < guideGraph.node.size(); ++nodeIndex) {
		vector<ScaffoldPart> &contigRef = guideGraph.node[nodeIndex].contig;
		vector<platanus::Position> altNodePosition;
		this->setCorrespondingNodePosition(altNodePosition, contigRef);

		vector<long> smoothedAltNodeID(altNodePosition.size(), 0);
		for (unsigned long i = 0; i < altNodePosition.size(); ++i)
			smoothedAltNodeID[i] = altNodePosition[i].id;

		smoothNodeIDVector(smoothedAltNodeID, this->node[nodeIndex], 0.0);

		long preID = 0;
		for (unsigned long i = 0; i < altNodePosition.size(); ++i) {
			if (preID != altNodePosition[i].id && smoothedAltNodeID[i] != 0 && altNodePosition[i].id != 0 && smoothedAltNodeID[i] != altNodePosition[i].id) {
				long contigIndex1 = id2Index(contigRef[i].id);
				long contigIndex2 = id2Index(this->node[id2Index(altNodePosition[i].id)].contig[altNodePosition[i].offset].id);

				if (bubbleDivideFlag || !((this->contigState[contigIndex1] | this->contigState[contigIndex2]) & (DBG_CONTIG_PRIMARY_BUBBLE | DBG_CONTIG_SECONDARY_BUBBLE))) {
					if (this->contigPositionInScaffold[contigIndex1].id != 0)
						breakpointFlag[id2Index(this->contigPositionInScaffold[contigIndex1].id)][this->contigPositionInScaffold[contigIndex1].offset] = 1;

					if (altNodePosition[i].id > 0)
						breakpointFlag[id2Index(altNodePosition[i].id)][altNodePosition[i].offset] = 1;
					else
						breakpointFlag[id2Index(altNodePosition[i].id)][altNodePosition[i].offset + 1] = 1;
				}
			}

			preID = altNodePosition[i].id;
		}
	}

    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();


	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;
		long i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[nodeIndex][i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;

			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::setCorrespondingNodePosition(vector<platanus::Position> &positionVector, const vector<ScaffoldPart> &partVector)
{
	positionVector.resize(partVector.size());
	for (auto partItr = partVector.begin(); partItr != partVector.end(); ++partItr) {
		if (this->contigPositionInScaffold[id2Index(partItr->id)].id != 0) {
			positionVector[partItr - partVector.begin()].id = sign(partItr->id) * this->contigPositionInScaffold[id2Index(partItr->id)].id;
			positionVector[partItr - partVector.begin()].offset = this->contigPositionInScaffold[id2Index(partItr->id)].offset;
		}
		else if (this->contigPositionInScaffold[id2Index(this->oppositeBubbleContigID[id2Index(partItr->id)])].id != 0) {
			positionVector[partItr - partVector.begin()].id = sign(partItr->id) * this->contigPositionInScaffold[id2Index(this->oppositeBubbleContigID[id2Index(partItr->id)])].id;
			positionVector[partItr - partVector.begin()].offset = this->contigPositionInScaffold[id2Index(this->oppositeBubbleContigID[id2Index(partItr->id)])].offset;
		}
	}
}

void PairedDBG::deleteEdgesUsingGuideGraph(PairedDBG &guideGraph)
{
	vector<unsigned char> deleteNodeFlag(this->numNode, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < guideGraph.node.size(); ++nodeIndex) {
		vector<ScaffoldPart> &contigRef = guideGraph.node[nodeIndex].contig;
		vector<platanus::Position> altNodePosition;
		this->setCorrespondingNodePosition(altNodePosition, contigRef);

		vector<long> smoothedAltNodeID(altNodePosition.size(), 0);
		for (unsigned long i = 0; i < altNodePosition.size(); ++i)
			smoothedAltNodeID[i] = altNodePosition[i].id;
		smoothNodeIDVector(smoothedAltNodeID, this->node[nodeIndex], 0.0);

		for (unsigned long i = 0; i < altNodePosition.size(); ++i) {
			if (smoothedAltNodeID[i] != 0 && altNodePosition[i].id != 0 && smoothedAltNodeID[i] != altNodePosition[i].id) {
				deleteNodeFlag[id2Index(altNodePosition[i].id)] = 1;
			}
		}
	}


    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing edges from nodes unaligned to the guide consensus graph ..." << std::endl;

    for (long nodeIndex = 0 ; nodeIndex < numNode; ++nodeIndex) {
		if (!(deleteNodeFlag[nodeIndex]))
			continue;

        GraphNode &node1 = this->node[nodeIndex];
        for (long edgeIndex = 0; edgeIndex < node1.numEdge; ++edgeIndex) {
			ids.push_back(nodeIndex + 1);
			ids.push_back(node1.edge[edgeIndex].end);
			++numDelete;
        }
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

void PairedDBG::smoothNodeIDVector(vector<long> &nodeIDVector, const GraphNode &targetNode, const double scoreFactor)
{
	pair<long, long> ends(0, nodeIDVector.size());
	vector<pair<long, long> > endsStack(1, ends);

	while (endsStack.size() > 0) {
		ends = endsStack.back();
		endsStack.pop_back();

		pair<long, long> newEnds = fillMajorityIDRun(nodeIDVector, targetNode, ends, scoreFactor);
		if (newEnds != ends) {
			endsStack.emplace_back(ends.first, newEnds.first);
			endsStack.emplace_back(newEnds.second, ends.second);
		}
	}
}


void PairedDBG::extendUsingGuideGraph(PairedDBG &guideGraph)
{
	vector<vector<char> > breakpointFlag(this->node.size());
	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		breakpointFlag[nodeIndex].assign(this->node[nodeIndex].contig.size() + 1, 0);
		breakpointFlag[nodeIndex].front() = 1;
		breakpointFlag[nodeIndex].back() = 1;
	}

	for (unsigned long nodeIndex = 0; nodeIndex < guideGraph.node.size(); ++nodeIndex) {
		vector<ScaffoldPart> &contigRef = guideGraph.node[nodeIndex].contig;
		vector<platanus::Position> altNodePosition;
		this->setCorrespondingNodePosition(altNodePosition, contigRef);

		for (unsigned long i = 0; i < altNodePosition.size(); ++i) {
			char c = 'N';
			if (this->node[id2Index(altNodePosition[i].id)].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE))
				c = 'B';
			std::cout << altNodePosition[i].id << c << ' ';
		}
		std::cout << endl;


		vector<long> smoothedAltNodeID(altNodePosition.size(), 0);
		for (unsigned long i = 0; i < altNodePosition.size(); ++i)
			smoothedAltNodeID[i] = altNodePosition[i].id;
		smoothNodeIDVector(smoothedAltNodeID, this->node[nodeIndex], 0.0);


		for (unsigned long i = 0; i < altNodePosition.size(); ++i) {
			char c = 'N';
			if (this->node[id2Index(smoothedAltNodeID[i])].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE))
				c = 'B';
/*
			if (altNodePosition[i].offset == 0)
				c = 'L';
			else if (altNodePosition[i].offset == this->node[id2Index(altNodePosition[i].id)].numContig - 1)
				c = 'R';
*/
			std::cout << smoothedAltNodeID[i] << c << ' ';
		}
		std::cout << endl;


		vector<std::pair<unsigned long, unsigned long> > boundary;
		unsigned long i = 0;
		while (i < altNodePosition.size()) {
			unsigned long j = i + 1;
			while (j < altNodePosition.size() && smoothedAltNodeID[i] == smoothedAltNodeID[j])
				++j;
			boundary.push_back(std::make_pair(i, j));
			i = j;
		}

		for (auto it = boundary.begin(); it != boundary.end(); ++it) {
			long alignedNodeID = altNodePosition[it->first].id;

/*
			long contigIndex1 = id2Index(this->node[id2Index(alignedNodeID)].contig[altNodePosition[it->first].offset].id);
			long contigIndex2 = id2Index(this->node[id2Index(altNodePosition[i].id)].contig[altNodePosition[i].offset].id);

			if (!((this->contigState[contigIndex1] | this->contigState[contigIndex2]) & (DBG_CONTIG_PRIMARY_BUBBLE | DBG_CONTIG_SECONDARY_BUBBLE))) {
				if (this->contigPositionInScaffold[contigIndex1].id != 0)
					breakpointFlag[id2Index(this->contigPositionInScaffold[contigIndex1].id)][this->contigPositionInScaffold[contigIndex1].offset] = 1;

				if (altNodePosition[i].id > 0)
					breakpointFlag[id2Index(altNodePosition[i].id)][altNodePosition[i].offset] = 1;
				else
					breakpointFlag[id2Index(altNodePosition[i].id)][altNodePosition[i].offset + 1] = 1;
			}
*/

			std::cout << altNodePosition[it->first].id << ' ';
			if (this->node[id2Index(alignedNodeID)].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE))
				std::cout << 'B';
			if (alignedNodeID > 0) {
				if (altNodePosition[it->first].offset == 0 && altNodePosition[it->second - 1].offset == this->node[id2Index(alignedNodeID)].numContig - 1)
					std::cout << "E(" << altNodePosition[it->first].offset << ", " << altNodePosition[it->second - 1].offset << ") n: " <<  this->node[id2Index(alignedNodeID)].numContig;
				else
					std::cout << "(" << altNodePosition[it->first].offset << ", " << altNodePosition[it->second - 1].offset << ") n: " <<  this->node[id2Index(alignedNodeID)].numContig;
			}
			else if (alignedNodeID < 0) {
				if (altNodePosition[it->first].offset == this->node[id2Index(alignedNodeID)].numContig - 1 && altNodePosition[it->second - 1].offset == 0 )
					std::cout << "E(" << altNodePosition[it->first].offset << ", " << altNodePosition[it->second - 1].offset << ") n: " <<  this->node[id2Index(alignedNodeID)].numContig;
				else
					std::cout << "(" << altNodePosition[it->first].offset << ", " << altNodePosition[it->second - 1].offset << ") n: " <<  this->node[id2Index(alignedNodeID)].numContig;
			}
			else {
				std::cout << "0(" << altNodePosition[it->first].offset << ", " << altNodePosition[it->second - 1].offset << ")";
			}
			std::cout << endl;
		}
	}
}

long PairedDBG::getQuartileLengthOfBubble(const unsigned long quartileNumber)
{
	vector<long> lengthBuffer;

    for (long nodeIndex = 0; nodeIndex < this->numNode; ++nodeIndex) {
		if (this->node[nodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE))
			lengthBuffer.push_back(this->node[nodeIndex].length);
	}

    sort(lengthBuffer.begin(), lengthBuffer.end());

	for (unsigned i = 0; i < lengthBuffer.size(); ++i) {
		if (i/lengthBuffer.size() >= quartileNumber/4)
			return lengthBuffer[i];
	}

	return 1;
}

void PairedDBG::detectRepeat(const double averageCoverage)
{
//    const double coverageThreshold = averageCoverage * 1.75;

    // # pragma omp for schedule(dynamic)
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
//        if (node[nodeID].numContig == 1 && calcNodeCoverage(node[nodeID]) > coverageThreshold) {
//            node[nodeID].state |= SC_REP;
//            continue;
//        }
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {
                long edgeEnd1, edgeEnd2;
                const GraphEdge &edge1 = (node[nodeID].edge[edgeID1]);
                const GraphEdge &edge2 = (node[nodeID].edge[edgeID2]);

                if (edge1.direction * edge2.direction < 0) continue;
                GraphNode &node1 = (node[abs(edge1.end) - 1]);
                if (edge1.length + node1.length <= edge2.length) continue;
                GraphNode &node2 = (node[abs(edge2.end) - 1]);
                if (edge2.length + node2.length <= edge1.length) continue;

                if (edge1.isForward()) {
                    edgeEnd1 = edge1.end;
                    edgeEnd2 = edge2.end;
                } else {
                    edgeEnd1 = edge2.end;
                    edgeEnd2 = edge1.end;
                }
                if ((abs(edge1.length + node1.length - edge2.length) <= tolerence + this->getScaffoldOverlap(edgeEnd1, edgeEnd2)
                    || abs(edge2.length + node2.length - edge1.length) <= tolerence + this->getScaffoldOverlap(edgeEnd2, edgeEnd1))
					&& abs(edge1.end) != abs(node2.oppositeBubbleNodeID)) continue;

//std::cout << "detect " << nodeID + 1 << ", " << edge1.end << ", " << edge2.end << std::endl;
                node[nodeID].state |= SC_REP;
                edgeID1 = node[nodeID].numEdge;
                break;
            }
        }
    }
}

void PairedDBG::deleteRepeatEdge()
{
//    const double coverageThreshold = averageCoverage * 1.75;

    vector<long> ids;

    cerr << "deleting edges from repeat contigs..." << endl;

    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        if (node[nodeID].numContig == 1) continue;
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {

                const GraphEdge &edge1 = (node[nodeID].edge[edgeID1]);
                const GraphEdge &edge2 = (node[nodeID].edge[edgeID2]);
                const GraphNode &node1 = (node[abs(edge1.end) - 1]);
                const GraphNode &node2 = (node[abs(edge2.end) - 1]);

				bool collisionFlag = true; 
                if (!this->checkDeleteEdge(edge1, edge2, node1, node2) && abs(edge1.end) != abs(node2.oppositeBubbleNodeID))
					collisionFlag = false;

                for (long m = 0; m < node[nodeID].numContig; ++m) {
//                    if (coverage[id2Index(node[nodeID].contig[m].id)] <= coverageThreshold && !collisionFlag)
                    if (!collisionFlag)
						continue;

					long numLinkToSubtract = std::min(edge1.breakdown[m], edge2.breakdown[m]);

                    for (long n = 0; n < node[nodeID].numEdge; ++n) {
                        node[nodeID].edge[n].numLink -= numLinkToSubtract;
                        node[nodeID].edge[n].breakdown[m] -= numLinkToSubtract;
                    }
                    contigPositionInScaffold[abs(node[nodeID].contig[m].id) - 1].id = 0;
                }
            }
        }
    }

    for (long i = 0; i < numNode; ++i) {
        for (long j = 0; j < node[i].numEdge; ++j) {
            if (node[i].edge[j].numLink < minLink) {
                ids.push_back(i + 1);
                ids.push_back(node[i].edge[j].end);
            }
        }
    }
    this->deleteEdges(ids);

}

void PairedDBG::deleteEdgeFromDifferentPreviousParent(const long numThread)
{
	if (this->contigPreviousParentNodeID.size() != this->contig.size())
		return;

    vector<long> ids;

    cerr << "deleting edges from previously divided points ..." << endl;

    omp_set_num_threads(numThread);
    #pragma omp parallel for schedule(dynamic)
    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		const GraphNode &node1 = node[nodeIndex];
        for (long edgeIndex = 0; edgeIndex < node[nodeIndex].numEdge - 1; ++edgeIndex) {
			const GraphNode &node2 = node[id2Index(node[nodeIndex].edge[edgeIndex].end)];
			for (long contigIndex1 = 0; contigIndex1 < node1.numContig; ++contigIndex1) {
				if (contigPreviousParentNodeID[id2Index(node1.contig[contigIndex1].id)] == 0)
					break;

				for (long contigIndex2 = 0; contigIndex2 < node2.numContig; ++contigIndex2) {
					if (contigPreviousParentNodeID[id2Index(node2.contig[contigIndex2].id)] == 0)
						break;

					if (contigPreviousParentNodeID[id2Index(node1.contig[contigIndex1].id)] == contigPreviousParentNodeID[id2Index(node2.contig[contigIndex2].id)]) {
                        node[nodeIndex].edge[edgeIndex].numLink -= node1.edge[edgeIndex].breakdown[contigIndex1];
						node[nodeIndex].edge[edgeIndex].breakdown[contigIndex1] = 0;
					}
				}
			}
		}
	}

    long numDelete = 0;
    for (long i = 0; i < numNode; ++i) {
        for (long j = 0; j < node[i].numEdge; ++j) {
            if (node[i].edge[j].numLink < minLink) {
                ids.push_back(i + 1);
                ids.push_back(node[i].edge[j].end);
				++numDelete;
            }
        }
    }
    this->deleteEdges(ids);

	cerr << "NUM_REMOVED_EDGES =" << numDelete / 2 << endl;
}

void PairedDBG::getUniqueConflictingNode(const long sourceNodeIndex, const char targetDirection, vector<NodeIDWithGap> &nodeIDBuffer)
{
	nodeIDBuffer.clear();
	GraphNode &sourceNode = this->node[sourceNodeIndex];

	for (long edgeIndex1 = 0; edgeIndex1 < sourceNode.numEdge - 1; ++edgeIndex1) {
		GraphEdge &edge1 = sourceNode.edge[edgeIndex1];
		if (edge1.direction != targetDirection || edge1.numLink < minLink)
			continue;

		for (long edgeIndex2 = edgeIndex1 + 1; edgeIndex2 < sourceNode.numEdge; ++edgeIndex2) {
			GraphEdge &edge2 = sourceNode.edge[edgeIndex2];
			if (edge2.direction != targetDirection || edge2.numLink < minLink)
				continue;

			if (this->checkDeleteEdge(edge1, edge2, this->node[id2Index(edge1.end)], this->node[id2Index(edge2.end)]) || abs(edge1.end) == abs(this->node[id2Index(edge2.end)].oppositeBubbleNodeID)) {
				if (!nodeIDBuffer.empty()) {
					nodeIDBuffer.clear();
					return;
				}
				nodeIDBuffer.emplace_back(edge1.end, edge1.length);
				nodeIDBuffer.emplace_back(edge2.end, edge2.length);
			}
		}
	}
}

void PairedDBG::divideNestedBubbleNode(const long numThread)
{
    cerr << "dividing nested bubble scaffolds..." << endl;

	setOppositeBubbleNodeIDForEachNode(numThread);

	const long MIN_OVERLAP_TO_JOIN = 20;
    long defaultMinOverlap = this->minOverlap;
    this->minOverlap = MIN_OVERLAP_TO_JOIN;

	vector<char> nestedBubbleFlag(this->node.size(), false);

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (node[nodeIndex].state & SC_DEL)
			continue;

		long oppositeNodeID = node[nodeIndex].oppositeBubbleNodeID;
		if (oppositeNodeID != 0) {
			long oppositeNodeIndex = id2Index(oppositeNodeID);
			if (node[oppositeNodeIndex].oppositeBubbleNodeID != 0 && id2Index(node[oppositeNodeIndex].oppositeBubbleNodeID) != nodeIndex)
				nestedBubbleFlag[nodeIndex] = nestedBubbleFlag[oppositeNodeIndex] = true;
		}
	}


    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (node[nodeIndex].state & SC_DEL)
			continue;

		long oppositeNodeID = node[nodeIndex].oppositeBubbleNodeID;
		if (oppositeNodeID != 0 && nestedBubbleFlag[id2Index(oppositeNodeID)])
			nestedBubbleFlag[nodeIndex] =  true;
	}


    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		long i;
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.back() = 1;

		if (nestedBubbleFlag[nodeIndex]) {
			for (i = 1; i < static_cast<long>(contigRef.size()); ++i) {
				if ((this->contigState[id2Index(contigRef[i].id)] & DBG_CONTIG_BUBBLE_JUNCTION) || contigRef[i].start > contigRef[i - 1].end || this->getOverlap(contigRef[i - 1].id, contigRef[i].id) < minOverlap)
					breakpointFlag[i] = 1;
			}
		}

		i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;


			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->minOverlap = defaultMinOverlap;
    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::deleteLongEdge(const long maxEdgeLength)
{
    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing long-gap edges (GAP_LENGTH > " <<  maxEdgeLength << ") ..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
            if (node1.edge[edgeID].length > maxEdgeLength) {
				ids.push_back(nodeID + 1);
				ids.push_back(node1.edge[edgeID].end);
                ++numDelete;
            }
        }
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

long PairedDBG::deleteErroneousEdgeNumTagRate(const bool overlapFlag, const long numThread)
{
    vector<long> ids;
    long numDelete = 0;
	const double RATE_THRESHOLD = 0.25;

    omp_set_num_threads(numThread);

    # pragma omp parallel for schedule(dynamic) reduction(+: numDelete)
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {

                const GraphEdge &edge1 = (node[nodeID].edge[edgeID1]);
                const GraphEdge &edge2 = (node[nodeID].edge[edgeID2]);
                const GraphNode &node1 = (node[abs(edge1.end) - 1]);
                const GraphNode &node2 = (node[abs(edge2.end) - 1]);

				if (overlapFlag) {
					if (!(edge1.state & edge2.state & DBG_OVERLAP))
						continue;
				}
				else {
					if (!this->checkDeleteEdge(edge1, edge2, node1, node2))
						continue;
				}

				long numTag1 = this->getCommonTagBetweenNodePair(nodeID + 1, edge1.end);
				long numTag2 = this->getCommonTagBetweenNodePair(nodeID + 1, edge2.end);

                if (numTag1 < RATE_THRESHOLD * numTag2) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeID + 1);
                        ids.push_back(edge1.end);
                    }
                    ++numDelete;
                } else if (numTag2 < RATE_THRESHOLD * numTag1) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeID + 1);
                        ids.push_back(edge2.end);
                    }
                    ++numDelete;
                }
            }
        }
    }
    this->deleteEdges(ids);

    return numDelete;
}

void PairedDBG::deleteErroneousEdgeNumTagRateIterative(const bool overlapFlag, const long numThread)
{
	if (this->tagLibraryMT == NULL)
		return;

	countMappedTagForEachScaffold(numThread);

    cerr << "removing erroneous edges using linked-reads (barcodes)..." << endl << endl;;
    long totalDelete = 0;
    long numDelete;
    do {
        numDelete = this->deleteErroneousEdgeNumTagRate(overlapFlag, numThread);
        totalDelete += numDelete;
        cerr << "NUM_REMOVED_EDGES =" << numDelete << endl;
    } while (numDelete > 0);

    cerr << "TOTAL_NUM_REMOVED_EDGES_BY_BARCODE =" << totalDelete << endl << endl;
}

void PairedDBG::deleteEdgeFromShortNodeKeepingBubble(const long lengthThreshold)
{
    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing edges from short nodes (LENGTH < " <<  lengthThreshold << ") ..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
		if (node1.length < lengthThreshold && node1.oppositeBubbleNodeID == 0) {
			for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
				ids.push_back(nodeID + 1);
				ids.push_back(node1.edge[edgeID].end);
				++numDelete;
			}
		}
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

long PairedDBG::divideInconsistentBubbleEnd()
{
    cerr << "dividing indonsistent bubble-ends..." << endl;

	this->setOppositeBubbleNodeIDAndStateForEachNode();

    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<vector<char> > breakpointFlag(this->node.size());
	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		breakpointFlag[nodeIndex].assign(this->node[nodeIndex].contig.size() + 1, 0);
		breakpointFlag[nodeIndex].front() = 1;
		breakpointFlag[nodeIndex].back() = 1;
	}

	vector<char> leftEndFlag(this->numContig, 0);
	vector<char> rightEndFlag(this->numContig, 0);

	long totalDivision = -1;
	long numDivision = 1;
	while (numDivision > 0) {
		totalDivision += numDivision;
		numDivision = 0;
		for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
			if (!(node[nodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)))
				continue;

			vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;
			for (unsigned long i = 0; i < contigRef.size(); ++i) {
				if (breakpointFlag[nodeIndex][i]) {
					if (contigRef[i].id > 0)
						leftEndFlag[id2Index(contigRef[i].id)] = 1;
					else
						rightEndFlag[id2Index(contigRef[i].id)] = 1;
				}

				if (breakpointFlag[nodeIndex][i + 1]) {
					if (contigRef[i].id > 0)
						rightEndFlag[id2Index(contigRef[i].id)] = 1;
					else
						leftEndFlag[id2Index(contigRef[i].id)] = 1;
				}
			}
		}

		for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
			if (!(node[nodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)))
				continue;

			vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;
			for (unsigned long i = 0; i < contigRef.size(); ++i) {
				long oppositeID = sign(contigRef[i].id) *  this->oppositeBubbleContigID[id2Index(contigRef[i].id)];
//				if (oppositeID == 0 || contigPositionInScaffold[id2Index(oppositeID)].id != abs(node[nodeIndex].oppositeBubbleNodeID))
				if (oppositeID == 0 || (abs(contigPositionInScaffold[id2Index(oppositeID)].id) != (nodeIndex + 1) && abs(contigPositionInScaffold[id2Index(oppositeID)].id) != abs(node[nodeIndex].oppositeBubbleNodeID)))
					continue;

				if (oppositeID > 0) {
					if (leftEndFlag[id2Index(oppositeID)] && !breakpointFlag[nodeIndex][i]) {
						breakpointFlag[nodeIndex][i] = 1;
						++numDivision;
					}
					if (rightEndFlag[id2Index(oppositeID)] && !breakpointFlag[nodeIndex][i + 1]) {
						breakpointFlag[nodeIndex][i + 1] = 1;
						++numDivision;
					}
				}
				else {
					if (leftEndFlag[id2Index(oppositeID)] && !breakpointFlag[nodeIndex][i + 1]) {
						breakpointFlag[nodeIndex][i + 1] = 1;
						++numDivision;
					}
					if (rightEndFlag[id2Index(oppositeID)] && !breakpointFlag[nodeIndex][i]) {
						breakpointFlag[nodeIndex][i] = 1;
						++numDivision;
					}
				}
			}
		}
	}


	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		long i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[nodeIndex][i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;


			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);

	return totalDivision;
}

void PairedDBG::adjustOppositeBubbleNodeIDDirection()
{
	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
//		if (!(node[nodeIndex].state & DBG_PRIMARY_BUBBLE))
		if (!(node[nodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)))
			continue;

		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;
		for (long j = 0; j < 2; ++j) {
			long contigID = j == 0 ? contigRef.front().id : contigRef.back().id;
			long oppositeContigID = this->oppositeBubbleContigID[id2Index(contigID)];
			if (oppositeContigID != 0) {
				long oppositeNodeID = this->contigPositionInScaffold[id2Index(oppositeContigID)].id;
				node[nodeIndex].oppositeBubbleNodeID = sign(contigID) * sign(oppositeContigID) * sign(oppositeNodeID) * abs(node[nodeIndex].oppositeBubbleNodeID);
				break;
			}
		}
	}
}

void PairedDBG::deleteEdgeFromSecondaryBubble()
{
	vector<char> secondaryBubbleFlag(numNode, 0);
    for (long nodeIndex = 0 ; nodeIndex < numNode; ++nodeIndex) {
        GraphNode &node1 = this->node[nodeIndex];

		for (long contigIndex = 0; contigIndex < node1.numContig; ++contigIndex) {
			if (contigState[id2Index(node1.contig[contigIndex].id)] & DBG_CONTIG_SECONDARY_BUBBLE) {
				secondaryBubbleFlag[nodeIndex] = 1;
				break;
			}
		}
	}
				
    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing edges from bubble nodes ..." << std::endl;

    for (long nodeIndex = 0 ; nodeIndex < numNode; ++nodeIndex) {
        GraphNode &node1 = this->node[nodeIndex];

		for (long edgeIndex = 0; edgeIndex < node1.numEdge; ++edgeIndex) {
//			if (node1.oppositeBubbleNodeID == 0 && node[id2Index(node1.edge[edgeIndex].end)].oppositeBubbleNodeID == 0)
			if (!(secondaryBubbleFlag[nodeIndex]) && !(secondaryBubbleFlag[id2Index(node1.edge[edgeIndex].end)]))
				continue;

			ids.push_back(nodeIndex + 1);
			ids.push_back(node1.edge[edgeIndex].end);
			++numDelete;
		}
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

void PairedDBG::divideNodeBasedOnBubblesIterative(const bool strandFlag, const long numThread)
{
	unsigned currentMode = getMode();
	setMode(0);

	long total = 0;
	long num;

	cerr << endl << "dividing nodes based on bubbles ..." << endl;
	do {
		num = divideNodeUsingBubbleContigPair(numThread);
//		num += divideBubbleJunctionNode(false);
		num += divideInconsistentBubbleEnd();
//		num += divideBubbleContigInNonHeteroNode();
		if (strandFlag) 
			num += divideNodeUsingBubbleContigPairStrandAware(numThread);

		total += num;
		cerr << "NUM_DIVISION = " << num << endl;
	} while (num > 0);

	cerr << "TOTAL_NUM_DIVISIONS =" << total << endl << endl;
	setMode(currentMode);
}

void PairedDBG::divideNodeByPrimaryBubbleBoundary(const PairedDBG &bubbleGraph)
{
    cerr << "dividing scaffolds by bubble-boundaries ..." << endl;

	std::vector<std::vector<long> > bubbleContigID;
    for (long nodeIndex = 0; nodeIndex < bubbleGraph.numNode; ++nodeIndex) {
		if (!(bubbleGraph.node[nodeIndex].state & DBG_PRIMARY_BUBBLE))
			continue;

		bubbleContigID.resize(bubbleContigID.size() + 1);
		bubbleContigID.back().resize(bubbleGraph.node[nodeIndex].contig.size());
		for (unsigned i = 0; i < bubbleGraph.node[nodeIndex].contig.size(); ++i)
			bubbleContigID.back()[i] = bubbleGraph.node[nodeIndex].contig[i].id;
	}

	vector<vector<long> > bubbleHeadIndex(this->numContig);
	for (unsigned long i = 0; i < bubbleContigID.size(); ++i)
		bubbleHeadIndex[id2Index(bubbleContigID[i][0])].push_back(i);


    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();
	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.front() = 1;
		breakpointFlag.back() = 1;

		for (unsigned long i = 0; i < contigRef.size(); ++i) {
			long contigIndex = id2Index(contigRef[i].id);
			if (bubbleHeadIndex[contigIndex].empty())
				continue;

			for (unsigned long j = 0; j < bubbleHeadIndex[contigIndex].size(); ++j) {
				std::vector<long> &bubbleContigIDRef = bubbleContigID[bubbleHeadIndex[contigIndex][j]];
				unsigned long k;
				if (contigRef[i].id == bubbleContigIDRef[0]) {
					for (k = 1; k < bubbleContigIDRef.size(); ++k) {
						if (contigRef[i + k].id != bubbleContigIDRef[k])
							break;
					}
					if (k == bubbleContigIDRef.size()) {
						breakpointFlag[i] = 1;
						breakpointFlag[i + k] = 1;
					}	
				}
				else {
					for (k = 1; k < bubbleContigIDRef.size(); ++k) {
						if (contigRef[i - k].id != -(bubbleContigIDRef[k]))
							break;
					}
					if (k == bubbleContigIDRef.size()) {
						breakpointFlag[i + 1] = 1;
						breakpointFlag[i - k + 1] = 1;
					}	
				}
			}
		}

		long i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;

			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::remakeGraphRecoveringSecondaryBubble(PairedDBG &bubbleGraph)
{
    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

    for (long i = 0; i < numNode; ++i) {
        if (node[i].state & SC_DEL)
			continue;

		fwrite(&(node[i].numContig), sizeof(long), 1, scaffoldFP);
		for (long j = 0; j < node[i].numContig; ++j)
			fwrite(&(node[i].contig[j]), sizeof(ScaffoldPart), 1, scaffoldFP);
		newContigPoolSize += node[i].numContig;
		++numNewNode;
    }

    for (long i = 0; i < bubbleGraph.numNode; ++i) {
        if ((bubbleGraph.node[i].state & SC_DEL) || !(bubbleGraph.node[i].state & DBG_SECONDARY_BUBBLE))
			continue;

		fwrite(&(bubbleGraph.node[i].numContig), sizeof(long), 1, scaffoldFP);
		for (long j = 0; j < bubbleGraph.node[i].numContig; ++j)
			fwrite(&(bubbleGraph.node[i].contig[j]), sizeof(ScaffoldPart), 1, scaffoldFP);
		newContigPoolSize += bubbleGraph.node[i].numContig;
		++numNewNode;
    }

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::divideErroneousLink(const vector<vector<unsigned> >& numErroneousPair, const vector<vector<unsigned> >& numSpanningPair, const vector<vector<double> >& sumExpectedLink, std::unordered_map<std::pair<int, int>, bool, platanus::PairHash, platanus::PairEqual> &errorLink, const long minLink, const DIVISION_MODE divisionMode, const long numThread, const long maxGapSize)
{
	const long MIN_OVERLAP_TO_JOIN = 32;

    long defaultMinOverlap = this->minOverlap;
    this->minOverlap = MIN_OVERLAP_TO_JOIN;

	if (this->contigPreviousParentNodeID.size() != this->contig.size())
		this->contigPreviousParentNodeID.resize(this->contig.size(), 0);

    unsigned long numDivided = 0;
    long numNewNode = 0;
    long newContigPoolSize = 0;
    bool isSplit = false;
    FILE *scaffoldFP = platanus::makeTemporaryFile();
    vector<char> breakpoint;
    for (long i = 0; i < numNode; ++i) {
        if (node[i].numContig == 1) {
            fwrite(&(node[i].numContig), sizeof(long), 1, scaffoldFP);
            fwrite(&(node[i].contig[0]), sizeof(ScaffoldPart), 1, scaffoldFP);
            ++numNewNode;
            ++newContigPoolSize;
            continue;
        }

        isSplit = false;
        breakpoint.resize(node[i].numContig + 1, 0);
        std::fill(breakpoint.begin(), breakpoint.end(), 0);
        breakpoint.back() = 1;
        for (unsigned long j = 0; j < numSpanningPair[i].size(); ++j) {
			if (divisionMode == SWITCH) {
//				if (numErroneousPair[i][j] <= numSpanningPair[i][j])
				if (numErroneousPair[i][j] <= numSpanningPair[i][j] || numErroneousPair[i][j] < minLink)
					continue;
			}
			else if (divisionMode == GAP) {
				if (node[i].contig[j].end - node[i].contig[j + 1].start < maxGapSize || numErroneousPair[i][j] < minLink)
					continue;
			}
			else {
//				long overlapLen = this->getOverlap(node[i].contig[j].id, node[i].contig[j + 1].id);
//				if (!(!(overlapLen >= minOverlap && node[i].contig[j].end - node[i].contig[j + 1].start <= 0) && numErroneousPair[i][j] >= minLink &&  numSpanningPair[i][j] < minLink && sumExpectedLink[i][j] > 1.0 && (double)numSpanningPair[i][j] < sumExpectedLink[i][j] * CHECK_USING_LONGER_LIB_TH))
				if (numErroneousPair[i][j] <= numSpanningPair[i][j])
					continue;
			}

			++numDivided;
			breakpoint[j + 1] = 1;
			isSplit = true;
        }
        if (isSplit) {
            long j = 0;
            while (j < node[i].numContig) {
                long start = node[i].contig[j].start;
                long k = j;
                while (breakpoint[j + 1] == 0) {
                    node[i].contig[j].start -= start;
                    node[i].contig[j].end -= start;
                    ++j;
                }
                node[i].contig[j].start -= start;
                node[i].contig[j].end -= start;
                ++j;
                long tmp = j - k;
                fwrite(&tmp, sizeof(long), 1, scaffoldFP);
                for (tmp = k; tmp < j; ++tmp)
                    fwrite(&(node[i].contig[tmp]), sizeof(ScaffoldPart), 1, scaffoldFP);
                ++numNewNode;
                newContigPoolSize += tmp;

				for (auto itr = node[i].contig.begin(); itr != node[i].contig.end(); ++itr) {
					if (this->contigPositionInScaffold[id2Index(itr->id)].id != 0)
						this->contigPreviousParentNodeID[id2Index(itr->id)] = i + 1;
					else
						this->contigPreviousParentNodeID[id2Index(itr->id)] = 0;
				}
            }
        }
        else {
            fwrite(&(node[i].numContig), sizeof(long), 1, scaffoldFP);
            auto contigIterator = node[i].contig.begin();
            auto contigEnd = node[i].contig.end();
            for (; contigIterator != contigEnd; ++contigIterator)
                fwrite(&(*contigIterator), sizeof(ScaffoldPart), 1, scaffoldFP);
            ++numNewNode;
            newContigPoolSize += node[i].numContig;
        }
    }

    std::cerr << "NUM_DIVIDED_SWITCH_ERROR_CANDIDATES = " << numDivided << endl;

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);

    this->minOverlap = defaultMinOverlap;
}

void PairedDBG::countPairsSpanningGap(vector<vector<unsigned> >& numSpanningPair, const DIVISION_MODE divisionMode, long numThread)
{
    const long averageInsSize = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();

    omp_set_num_threads(numThread);

    vector<vector<vector<unsigned> > > numSpanningPairThread(numThread);
	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		numSpanningPairThread[threadID].resize(numSpanningPair.size());
		for (unsigned i = 0; i < numSpanningPair.size(); ++i)
			numSpanningPairThread[threadID][i].resize(numSpanningPair[i].size(), 0);
	}

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		platanus::Position forwardResult;
		platanus::Position reverseResult;
		unsigned numPairLink;
		unsigned numReadLink;

		rewind((*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
		while (fread(&numPairLink, sizeof(unsigned), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);

			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);

				unsigned long forwardIndex = id2Index(forwardResult.id);
				unsigned long reverseIndex = id2Index(reverseResult.id);
				if (divisionMode == SWITCH && (this->oppositeBubbleContigID[forwardIndex] == 0 || this->oppositeBubbleContigID[reverseIndex] == 0))
					continue;

				if (contigPositionInScaffold[forwardIndex].id == 0)
					continue;
				forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[forwardIndex].id : -(contigPositionInScaffold[forwardIndex].id);
				forwardResult.offset = contigPositionInScaffold[forwardIndex].id > 0 ? forwardResult.offset : contig[forwardIndex].length - forwardResult.offset - 1;
				forwardResult.offset += node[abs(forwardResult.id) - 1].contig[contigPositionInScaffold[forwardIndex].offset].start;

				if (contigPositionInScaffold[reverseIndex].id == 0)
					continue;
				reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[reverseIndex].id : -(contigPositionInScaffold[reverseIndex].id);
				reverseResult.offset = contigPositionInScaffold[reverseIndex].id > 0 ? reverseResult.offset : contig[reverseIndex].length - reverseResult.offset - 1;
				reverseResult.offset += node[abs(reverseResult.id) - 1].contig[contigPositionInScaffold[reverseIndex].offset].start;

				if (forwardResult.id != -reverseResult.id || abs(forwardResult.offset - reverseResult.offset) - averageInsSize > tolerence)
					continue;

				unsigned long leftMostContigNum = std::min(contigPositionInScaffold[forwardIndex].offset, contigPositionInScaffold[reverseIndex].offset);
				unsigned long rightMostContigNum = std::max(contigPositionInScaffold[forwardIndex].offset, contigPositionInScaffold[reverseIndex].offset);

				for (unsigned long i = leftMostContigNum; i < rightMostContigNum; ++i) {
					std::pair<int, int> redundancyCheckKey = std::make_pair(abs(forwardResult.id) - 1, i);
					if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
						redundancyCheckSet.insert(redundancyCheckKey);
						++(numSpanningPairThread[threadID][abs(forwardResult.id) - 1][i]);
					}
				}
			}

			for (unsigned pairIndex = 0; pairIndex < numReadLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
			}
        }
    }

    for (long threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned i = 0; i < numNode; ++i) {
			for (unsigned j = 0; j < numSpanningPair[i].size(); ++j) {
				numSpanningPair[i][j] += numSpanningPairThread[threadID][i][j];
			}
		}
	}
}

void PairedDBG::countLongReadSpanningGap(vector<vector<unsigned> >& numSpanningPair, const DIVISION_MODE divisionMode, long numThread)
{
    omp_set_num_threads(numThread);

    vector<vector<vector<unsigned> > > numSpanningPairThread(numThread);
	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		numSpanningPairThread[threadID].resize(numSpanningPair.size());
		for (unsigned i = 0; i < numSpanningPair.size(); ++i)
			numSpanningPairThread[threadID][i].resize(numSpanningPair[i].size(), 0);
	}

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		unsigned numResult;
		long score;
		vector<platanus::Position> positionBuffer;

        rewind((*longReadLibraryMT)[threadID].mappedReadFP);
		while (fread(&numResult, sizeof(unsigned), 1, (*longReadLibraryMT)[threadID].mappedReadFP)) {
			positionBuffer.resize(numResult);
			for (unsigned long i = 0; i < numResult; ++i) {
				fread(&(positionBuffer[i]), sizeof(platanus::Position), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
				fread(&score, sizeof(long), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
			}

			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (auto positionItr1 = positionBuffer.begin(); positionItr1 != positionBuffer.end() - 1; ++positionItr1) {
				for (auto positionItr2 = positionItr1 + 1; positionItr2 != positionBuffer.end(); ++positionItr2) {
					unsigned long forwardIndex = id2Index(positionItr1->id);
					unsigned long reverseIndex = id2Index(positionItr2->id);
					if (divisionMode == SWITCH && (this->oppositeBubbleContigID[forwardIndex] == 0 || this->oppositeBubbleContigID[reverseIndex] == 0))
						continue;

					if (contigPositionInScaffold[forwardIndex].id == 0) continue;
					positionItr1->id = positionItr1->id > 0 ? contigPositionInScaffold[forwardIndex].id : -(contigPositionInScaffold[forwardIndex].id);
					positionItr1->offset = contigPositionInScaffold[forwardIndex].id > 0 ? positionItr1->offset : contig[forwardIndex].length - positionItr1->offset - 1;
					positionItr1->offset += node[abs(positionItr1->id) - 1].contig[contigPositionInScaffold[forwardIndex].offset].start;

					if (contigPositionInScaffold[reverseIndex].id == 0) continue;
					positionItr2->id = positionItr2->id > 0 ? contigPositionInScaffold[reverseIndex].id : -(contigPositionInScaffold[reverseIndex].id);
					positionItr2->offset = contigPositionInScaffold[reverseIndex].id > 0 ? positionItr2->offset : contig[reverseIndex].length - positionItr2->offset - 1;
					positionItr2->offset += node[abs(positionItr2->id) - 1].contig[contigPositionInScaffold[reverseIndex].offset].start;

					if (positionItr1->id != positionItr2->id) continue;

					unsigned long leftMostContigNum = std::min(contigPositionInScaffold[forwardIndex].offset, contigPositionInScaffold[reverseIndex].offset);
					unsigned long rightMostContigNum = std::max(contigPositionInScaffold[forwardIndex].offset, contigPositionInScaffold[reverseIndex].offset);

					for (unsigned long i = leftMostContigNum; i < rightMostContigNum; ++i) {
						std::pair<int, int> redundancyCheckKey = std::make_pair(abs(positionItr1->id) - 1, i);
						if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
							redundancyCheckSet.insert(redundancyCheckKey);
							++(numSpanningPairThread[threadID][abs(positionItr1->id) - 1][i]);
						}
					}
				}
			}
		}
    }

    for (long threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned i = 0; i < numNode; ++i) {
			for (unsigned j = 0; j < numSpanningPair[i].size(); ++j) {
				numSpanningPair[i][j] += numSpanningPairThread[threadID][i][j];
			}
		}
	}
}

void PairedDBG::countLinksInsideContigs(vector<vector<unsigned> >& numPair, const long numThread)
{
	const long averageInsSize = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();

    omp_set_num_threads(numThread);

    vector<vector<vector<unsigned> > > numPairThread(numThread);
	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		numPairThread[threadID].resize(numPair.size());
		for (unsigned i = 0; i < numPair.size(); ++i)
			numPairThread[threadID][i].resize(numPair[i].size(), 0);
	}

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		platanus::Position contigResultF;
		platanus::Position contigResultR;
		platanus::Position scaffoldResultF;
		platanus::Position scaffoldResultR;
		long scaffoldOverhangF;
		long scaffoldOverhangR;
		long contigOverhangF;
		long contigOverhangR;
		unsigned long contigIndexF;
		unsigned long contigIndexR;
		unsigned numPairLink;
		unsigned numReadLink;

		FILE *mappedFP = (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP;
        rewind(mappedFP);
		while (fread(&numPairLink, sizeof(unsigned), 1, mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, mappedFP);

			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&contigResultF, sizeof(platanus::Position), 1, mappedFP);
				fread(&contigResultR, sizeof(platanus::Position), 1, mappedFP);

				contigIndexF = abs(contigResultF.id) - 1;
				if (contigPositionInScaffold[contigIndexF].id == 0) continue;
				scaffoldResultF.id = contigResultF.id > 0 ? contigPositionInScaffold[contigIndexF].id : -(contigPositionInScaffold[contigIndexF].id);
				scaffoldResultF.offset = contigPositionInScaffold[contigIndexF].id > 0 ? contigResultF.offset : contig[contigIndexF].length - contigResultF.offset - 1;
				scaffoldResultF.offset += node[abs(scaffoldResultF.id) - 1].contig[contigPositionInScaffold[contigIndexF].offset].start;
				scaffoldOverhangF = scaffoldResultF.id > 0  ? node[scaffoldResultF.id - 1].length - scaffoldResultF.offset : scaffoldResultF.offset;

				contigIndexR = abs(contigResultR.id) - 1;
				if (contigPositionInScaffold[contigIndexR].id == 0) continue;
				scaffoldResultR.id = contigResultR.id > 0 ? contigPositionInScaffold[contigIndexR].id : -(contigPositionInScaffold[contigIndexR].id);
				scaffoldResultR.offset = contigPositionInScaffold[contigIndexR].id > 0 ? contigResultR.offset : contig[contigIndexR].length - contigResultR.offset - 1;
				scaffoldResultR.offset += node[abs(scaffoldResultR.id) - 1].contig[contigPositionInScaffold[contigIndexR].offset].start;
				scaffoldOverhangR = scaffoldResultR.id > 0  ? node[scaffoldResultR.id - 1].length - scaffoldResultR.offset : scaffoldResultR.offset;

				if (scaffoldResultF.id == -scaffoldResultR.id || scaffoldOverhangF + scaffoldOverhangR <= averageInsSize + tolerence)
					continue;

				contigOverhangF = contigResultF.id > 0  ? contig[contigIndexF].length - contigResultF.offset : contigResultF.offset;
				if (contigOverhangF > averageInsSize + tolerence) {
					if (scaffoldResultF.id > 0) {
						for (long i = contigPositionInScaffold[contigIndexF].offset; i < node[scaffoldResultF.id - 1].numContig - 1; ++i) {
							if (node[scaffoldResultF.id - 1].contig[i].end - scaffoldResultF.offset <= averageInsSize) {
								std::pair<int, int> redundancyCheckKey = std::make_pair(scaffoldResultF.id - 1, i);
								if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
									redundancyCheckSet.insert(redundancyCheckKey);
									++(numPairThread[threadID][scaffoldResultF.id - 1][i]);
								}
							}
							else {
								break;
							}
						}
					}
					else {
						for (long i = contigPositionInScaffold[contigIndexF].offset; i > 0; --i) {
							if (scaffoldResultF.offset - node[-(scaffoldResultF.id) - 1].contig[i].start <= averageInsSize) {
								std::pair<int, int> redundancyCheckKey = std::make_pair(-(scaffoldResultF.id) - 1, i - 1);
								if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
									redundancyCheckSet.insert(redundancyCheckKey);
									++(numPairThread[threadID][-(scaffoldResultF.id) - 1][i - 1]);
								}
							}
							else {
								break;
							}
						}
					}
				}

				contigOverhangR = contigResultR.id > 0  ? contig[contigIndexR].length - contigResultR.offset : contigResultR.offset;
				if (contigOverhangR > averageInsSize + tolerence) {
					if (scaffoldResultR.id > 0) {
						for (long i = contigPositionInScaffold[contigIndexR].offset; i < node[scaffoldResultR.id - 1].numContig - 1; ++i) {
							if (node[scaffoldResultR.id - 1].contig[i].end - scaffoldResultR.offset <= averageInsSize) {
								std::pair<int, int> redundancyCheckKey = std::make_pair(scaffoldResultR.id - 1, i);
								if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
									redundancyCheckSet.insert(redundancyCheckKey);
									++(numPairThread[threadID][scaffoldResultR.id - 1][i]);
								}
							}
							else {
								break;
							}
						}
					}
					else {
						for (long i = contigPositionInScaffold[contigIndexR].offset; i > 0; --i) {
							if (scaffoldResultR.offset - node[-(scaffoldResultR.id) - 1].contig[i].start <= averageInsSize) {
								std::pair<int, int> redundancyCheckKey = std::make_pair(-(scaffoldResultR.id) - 1, i - 1);
								if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
									redundancyCheckSet.insert(redundancyCheckKey);
									++(numPairThread[threadID][-(scaffoldResultR.id) - 1][i - 1]);
								}
							}
							else {
								break;
							}
						}
					}
				}
			}

			for (unsigned pairIndex = 0; pairIndex < numReadLink; ++pairIndex) {
				fread(&contigResultF, sizeof(platanus::Position), 1, mappedFP);
				fread(&contigResultR, sizeof(platanus::Position), 1, mappedFP);
			}
		}
    }

    for (long threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned i = 0; i < numNode; ++i) {
			for (unsigned j = 0; j < numPair[i].size(); ++j) {
				numPair[i][j] += numPairThread[threadID][i][j];
			}
		}
	}
}

void PairedDBG::countLongReadLinksInsideContigs(vector<vector<unsigned> >& numPair, const long numThread)
{
	const long averageInsSize = (*longReadLibraryMT)[0].getAverageInsSize();

    omp_set_num_threads(numThread);

    vector<vector<vector<unsigned> > > numPairThread(numThread);
	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		numPairThread[threadID].resize(numPair.size());
		for (unsigned i = 0; i < numPair.size(); ++i)
			numPairThread[threadID][i].resize(numPair[i].size(), 0);
	}

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		platanus::Position scaffoldResultF;
		platanus::Position scaffoldResultR;
		long scaffoldOverhangF;
		long scaffoldOverhangR;
		long contigOverhangF;
		long contigOverhangR;
		unsigned long contigIndexF;
		unsigned long contigIndexR;

		FILE *mappedFP = (*longReadLibraryMT)[threadID].mappedReadFP;
		unsigned numResult;
		long score;
		vector<platanus::Position> positionBuffer;

        rewind(mappedFP);
		while (fread(&numResult, sizeof(unsigned), 1, mappedFP)) {
			positionBuffer.resize(numResult);
			for (unsigned long i = 0; i < numResult; ++i) {
				fread(&(positionBuffer[i]), sizeof(platanus::Position), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
				fread(&score, sizeof(long), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
			}

			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (auto positionItr1 = positionBuffer.begin(); positionItr1 != positionBuffer.end() - 1; ++positionItr1) {
				for (auto positionItr2 = positionItr1 + 1; positionItr2 != positionBuffer.end(); ++positionItr2) {
					contigIndexF = abs(positionItr1->id) - 1;
					if (contigPositionInScaffold[contigIndexF].id == 0) continue;
					scaffoldResultF.id = positionItr1->id > 0 ? contigPositionInScaffold[contigIndexF].id : -(contigPositionInScaffold[contigIndexF].id);
					scaffoldResultF.offset = contigPositionInScaffold[contigIndexF].id > 0 ? positionItr1->offset : contig[contigIndexF].length - positionItr1->offset - 1;
					scaffoldResultF.offset += node[abs(scaffoldResultF.id) - 1].contig[contigPositionInScaffold[contigIndexF].offset].start;
					scaffoldOverhangF = scaffoldResultF.id > 0  ? node[scaffoldResultF.id - 1].length - scaffoldResultF.offset : scaffoldResultF.offset;

					contigIndexR = abs(positionItr2->id) - 1;
					if (contigPositionInScaffold[contigIndexR].id == 0) continue;
					scaffoldResultR.id = positionItr2->id > 0 ? contigPositionInScaffold[contigIndexR].id : -(contigPositionInScaffold[contigIndexR].id);
					scaffoldResultR.offset = contigPositionInScaffold[contigIndexR].id > 0 ? positionItr2->offset : contig[contigIndexR].length - positionItr2->offset - 1;
					scaffoldResultR.offset += node[abs(scaffoldResultR.id) - 1].contig[contigPositionInScaffold[contigIndexR].offset].start;
					scaffoldOverhangR = scaffoldResultR.id > 0  ? node[scaffoldResultR.id - 1].length - scaffoldResultR.offset : scaffoldResultR.offset;

					if (scaffoldResultF.id == scaffoldResultR.id || scaffoldOverhangF + scaffoldOverhangR <= averageInsSize + tolerence) continue;

					contigOverhangF = positionItr1->id > 0  ? contig[contigIndexF].length - positionItr1->offset : positionItr1->offset;
					if (contigOverhangF > averageInsSize + tolerence) {
						if (scaffoldResultF.id > 0) {
							for (long i = contigPositionInScaffold[contigIndexF].offset; i < node[scaffoldResultF.id - 1].numContig - 1; ++i) {
								if (node[scaffoldResultF.id - 1].contig[i].end - scaffoldResultF.offset <= averageInsSize) {
									std::pair<int, int> redundancyCheckKey = std::make_pair(scaffoldResultF.id - 1, i);
									if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
										redundancyCheckSet.insert(redundancyCheckKey);
										++(numPairThread[threadID][scaffoldResultF.id - 1][i]);
									}
								}
								else {
									break;
								}
							}
						}
						else {
							for (long i = contigPositionInScaffold[contigIndexF].offset; i > 0; --i) {
								if (scaffoldResultF.offset - node[-(scaffoldResultF.id) - 1].contig[i].start <= averageInsSize) {
									std::pair<int, int> redundancyCheckKey = std::make_pair(-(scaffoldResultF.id) - 1, i - 1);
									if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
										redundancyCheckSet.insert(redundancyCheckKey);
										++(numPairThread[threadID][-(scaffoldResultF.id) - 1][i - 1]);
									}
								}
								else {
									break;
								}
							}
						}
					}

					contigOverhangR = positionItr2->id > 0  ? contig[contigIndexR].length - positionItr2->offset : positionItr2->offset;
					if (contigOverhangR > averageInsSize + tolerence) {
						if (scaffoldResultR.id > 0) {
							for (long i = contigPositionInScaffold[contigIndexR].offset; i < node[scaffoldResultR.id - 1].numContig - 1; ++i) {
								if (node[scaffoldResultR.id - 1].contig[i].end - scaffoldResultR.offset <= averageInsSize) {
									std::pair<int, int> redundancyCheckKey = std::make_pair(scaffoldResultR.id - 1, i - 1);
									if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
										redundancyCheckSet.insert(redundancyCheckKey);
										++(numPairThread[threadID][scaffoldResultR.id - 1][i]);
									}
								}
								else {
									break;
								}
							}
						}
						else {
							for (long i = contigPositionInScaffold[contigIndexR].offset; i > 0; --i) {
								if (scaffoldResultR.offset - node[-(scaffoldResultR.id) - 1].contig[i].start <= averageInsSize) {
									std::pair<int, int> redundancyCheckKey = std::make_pair(-(scaffoldResultR.id) - 1, i - 1);
									if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
										++(numPairThread[threadID][-(scaffoldResultR.id) - 1][i - 1]);
									}
								}
								else {
									break;
								}
							}
						}
					}
				}
			}
        }
    }

    for (long threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned i = 0; i < numNode; ++i) {
			for (unsigned j = 0; j < numPair[i].size(); ++j) {
				numPair[i][j] += numPairThread[threadID][i][j];
			}
		}
	}
}

void PairedDBG::countSwitchErrorLinks(vector<vector<unsigned> >& numPair, const long numThread)
{
    const long averageInsSize = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();
	setOppositeBubbleNodeIDAndStateForEachNode();

    omp_set_num_threads(numThread);

    vector<vector<vector<unsigned> > > numPairThread(numThread);
	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		numPairThread[threadID].resize(numPair.size());
		for (unsigned i = 0; i < numPair.size(); ++i)
			numPairThread[threadID][i].resize(numPair[i].size(), 0);
	}

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		std::array<std::array<platanus::Position, 2>, 2> contigResult;
		std::array<std::array<platanus::Position, 2>, 2> scaffoldResult;
		std::array<std::array<unsigned long, 2>, 2> contigIndex;
		std::array<std::array<unsigned long, 2>, 2> scaffoldIndex;
		unsigned numPairLink;
		unsigned numReadLink;
		unsigned i, j;

		FILE *mappedFP;
		mappedFP = (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP;

        rewind(mappedFP);
		while (fread(&numPairLink, sizeof(unsigned), 1, mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, mappedFP);

			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&(contigResult[0][0]), sizeof(platanus::Position), 1, mappedFP);
				fread(&(contigResult[0][1]), sizeof(platanus::Position), 1, mappedFP);

				if (abs(contigResult[0][0].id) == abs(this->oppositeBubbleContigID[id2Index(contigResult[0][1].id)]))
					continue;

				for (i = 0; i < 2; ++i) {
					contigIndex[0][i] = abs(contigResult[0][i].id) - 1;
					if (contigPositionInScaffold[contigIndex[0][i]].id == 0 || oppositeBubbleContigID[contigIndex[0][i]] == 0)
						break;

					contigResult[1][i].id = sign(contigResult[0][i].id) * oppositeBubbleContigID[contigIndex[0][i]];
					contigResult[1][i].offset = contigResult[0][i].offset;
					contigIndex[1][i] = abs(contigResult[1][i].id) - 1;
					if (contigPositionInScaffold[contigIndex[1][i]].id == 0)
						break;
				}
				if (i != 2)
					continue;

				for (i = 0; i < 2; ++i) {
					for (j = 0; j < 2; ++j) {
						scaffoldResult[j][i].id = contigResult[j][i].id > 0 ? contigPositionInScaffold[contigIndex[j][i]].id : -(contigPositionInScaffold[contigIndex[j][i]].id);
						scaffoldResult[j][i].offset = contigPositionInScaffold[contigIndex[j][i]].id > 0 ? contigResult[j][i].offset : contig[contigIndex[j][i]].length - contigResult[j][i].offset - 1;
						scaffoldResult[j][i].offset += node[abs(scaffoldResult[j][i].id) - 1].contig[contigPositionInScaffold[contigIndex[j][i]].offset].start;

						scaffoldIndex[j][i] = abs(scaffoldResult[j][i].id) - 1;
					}
				}
				if (scaffoldIndex[0][0] != scaffoldIndex[1][1] || scaffoldIndex[1][0] != scaffoldIndex[0][1])
					continue;

				for (j = 0; j < 2; ++j) {
					if (scaffoldIndex[j][0] == scaffoldIndex[j][1] || (int)scaffoldIndex[j][0] != abs(node[id2Index(scaffoldResult[j][1].id)].oppositeBubbleNodeID) - 1)
						break;
				}
				if (j != 2)
					continue;


				for (j = 0; j < 2; ++j) {
					long insertLength = abs(scaffoldResult[j][0].offset - scaffoldResult[j^1][1].offset);
					if (insertLength - averageInsSize > tolerence)
						continue;

					unsigned long leftMostContigNum = std::min(contigPositionInScaffold[contigIndex[j][0]].offset, contigPositionInScaffold[contigIndex[j^1][1]].offset);
					unsigned long rightMostContigNum = std::max(contigPositionInScaffold[contigIndex[j][0]].offset, contigPositionInScaffold[contigIndex[j^1][1]].offset);

					for (unsigned long k = leftMostContigNum; k < rightMostContigNum; ++k) {
						std::pair<int, int> redundancyCheckKey = std::make_pair(scaffoldIndex[j][0], k);
						if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
							redundancyCheckSet.insert(redundancyCheckKey);
							++(numPairThread[threadID][scaffoldIndex[j][0]][k]);
						}
					}
				}
			}

			for (unsigned pairIndex = 0; pairIndex < numReadLink; ++pairIndex) {
				fread(&(contigResult[0][0]), sizeof(platanus::Position), 1, mappedFP);
				fread(&(contigResult[0][1]), sizeof(platanus::Position), 1, mappedFP);
			}
		}
    }

    for (long threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned i = 0; i < numNode; ++i) {
			for (unsigned j = 0; j < numPair[i].size(); ++j) {
				numPair[i][j] += numPairThread[threadID][i][j];
			}
		}
	}
}

void PairedDBG::countLongReadSwitchErrorLinks(vector<vector<unsigned> >& numPair, const long numThread)
{
	setOppositeBubbleNodeIDAndStateForEachNode();

    omp_set_num_threads(numThread);

    vector<vector<vector<unsigned> > > numPairThread(numThread);
	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		numPairThread[threadID].resize(numPair.size());
		for (unsigned i = 0; i < numPair.size(); ++i)
			numPairThread[threadID][i].resize(numPair[i].size(), 0);
	}

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		std::array<std::array<platanus::Position, 2>, 2> contigResult;
		std::array<std::array<platanus::Position, 2>, 2> scaffoldResult;
		std::array<std::array<unsigned long, 2>, 2> contigIndex;
		std::array<std::array<unsigned long, 2>, 2> scaffoldIndex;
		unsigned i, j;

		long score;
		unsigned numResult;
		FILE *mappedFP;
		vector<platanus::Position> positionBuffer;
		mappedFP = (*longReadLibraryMT)[threadID].mappedReadFP;

        rewind(mappedFP);
		while (fread(&numResult, sizeof(unsigned), 1, mappedFP)) {
			positionBuffer.resize(numResult);
			for (i = 0; i < numResult; ++i) {
				fread(&(positionBuffer[i]), sizeof(platanus::Position), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
				fread(&score, sizeof(long), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
			}

			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (auto positionItr1 = positionBuffer.begin(); positionItr1 != positionBuffer.end() - 1; ++positionItr1) {
				for (auto positionItr2 = positionItr1 + 1; positionItr2 != positionBuffer.end(); ++positionItr2) {
					contigResult[0][0] = *positionItr1;
					contigResult[0][1] = *positionItr2;

					if (abs(contigResult[0][0].id) == abs(this->oppositeBubbleContigID[id2Index(contigResult[0][1].id)]))
						continue;

					for (i = 0; i < 2; ++i) {
						contigIndex[0][i] = abs(contigResult[0][i].id) - 1;
						if (contigPositionInScaffold[contigIndex[0][i]].id == 0 || oppositeBubbleContigID[contigIndex[0][i]] == 0)
							break;

						contigResult[1][i].id = sign(contigResult[0][i].id) * oppositeBubbleContigID[contigIndex[0][i]];
						contigResult[1][i].offset = contigResult[0][i].offset;
						contigIndex[1][i] = abs(contigResult[1][i].id) - 1;
						if (contigPositionInScaffold[contigIndex[1][i]].id == 0)
							break;
					}
					if (i != 2)
						continue;

					for (i = 0; i < 2; ++i) {
						for (j = 0; j < 2; ++j) {
							scaffoldResult[j][i].id = contigResult[j][i].id > 0 ? contigPositionInScaffold[contigIndex[j][i]].id : -(contigPositionInScaffold[contigIndex[j][i]].id);
							scaffoldResult[j][i].offset = contigPositionInScaffold[contigIndex[j][i]].id > 0 ? contigResult[j][i].offset : contig[contigIndex[j][i]].length - contigResult[j][i].offset - 1;
							scaffoldResult[j][i].offset += node[abs(scaffoldResult[j][i].id) - 1].contig[contigPositionInScaffold[contigIndex[j][i]].offset].start;

							scaffoldIndex[j][i] = abs(scaffoldResult[j][i].id) - 1;
						}
					}
					if (scaffoldIndex[0][0] != scaffoldIndex[1][1] || scaffoldIndex[1][0] != scaffoldIndex[0][1])
						continue;

					for (j = 0; j < 2; ++j) {
						if (scaffoldIndex[j][0] == scaffoldIndex[j][1] || (int)scaffoldIndex[j][0] != abs(node[id2Index(scaffoldResult[j][1].id)].oppositeBubbleNodeID) - 1)
							break;
					}
					if (j != 2)
						continue;


					for (j = 0; j < 2; ++j) {
						unsigned long leftMostContigNum = std::min(contigPositionInScaffold[contigIndex[j][0]].offset, contigPositionInScaffold[contigIndex[j^1][1]].offset);
						unsigned long rightMostContigNum = std::max(contigPositionInScaffold[contigIndex[j][0]].offset, contigPositionInScaffold[contigIndex[j^1][1]].offset);

						for (unsigned long k = leftMostContigNum; k < rightMostContigNum; ++k) {
							std::pair<int, int> redundancyCheckKey = std::make_pair(scaffoldIndex[j][0], k);
							if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
								redundancyCheckSet.insert(redundancyCheckKey);
								++(numPairThread[threadID][scaffoldIndex[j][0]][k]);
							}
						}
					}
				}
			}
        }
    }

    for (long threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned i = 0; i < numNode; ++i) {
			for (unsigned j = 0; j < numPair[i].size(); ++j) {
				numPair[i][j] += numPairThread[threadID][i][j];
			}
		}
	}
}

void PairedDBG::divideErroneousNode(const long minLink, const DIVISION_MODE divisionMode, const long numThread, const long maxGapSize)
{
	long currentLibraryIndex = this->targetLibraryIndex;

	cerr << "dividing erroneous scaffolds..." << endl;

    vector<vector<unsigned> > numSpanningPair(numNode);
    vector<vector<unsigned> > numErroneousPair(numNode);
    for (unsigned i = 0; i < numNode; ++i) {
        numSpanningPair[i].resize(node[i].numContig - 1, 0);
        numErroneousPair[i].resize(node[i].numContig - 1, 0);
	}

    std::unordered_map<std::pair<int, int>, bool, platanus::PairHash, platanus::PairEqual> errorLink;
    vector<vector<double> > sumExpectedLink(numNode);
    omp_set_num_threads(numThread);

	if (allLibraryMT != NULL) {
		for (unsigned libraryID = 0; libraryID < (*allLibraryMT).size(); ++libraryID) {
			this->setTargetLibraryIndex(libraryID);
			this->setTolerence(minTolerenceFactor * (*allLibraryMT)[libraryID][0].getSDInsSize());

			this->countPairsSpanningGap(numSpanningPair, divisionMode, numThread);
			if (divisionMode == SWITCH)
				this->countSwitchErrorLinks(numErroneousPair, numThread);
			else
				this->countLinksInsideContigs(numErroneousPair, numThread);

			if (divisionMode != SWITCH) {
				for (unsigned i = 0; i < numNode; ++i) {
					sumExpectedLink[i].resize(node[i].numContig - 1, 0.0);
					for (unsigned j = 0; j < node[i].numContig - 1; ++j) {
						sumExpectedLink[i][j] += this->calcExpectedLink((*allLibraryMT)[libraryID][0].getAverageCoverage(), node[i].contig[j].end, node[i].length - node[i].contig[j + 1].start, node[i].contig[j + 1].start - node[i].contig[j].end);
					}
				}
			}
		}
	}

	if (longReadLibraryMT != NULL) {
		this->setTolerence((*longReadLibraryMT)[0].getAverageInsSize());

		this->countLongReadSpanningGap(numSpanningPair, divisionMode, numThread);
		if (divisionMode == SWITCH)
			this->countLongReadSwitchErrorLinks(numErroneousPair, numThread);
		else
			this->countLongReadLinksInsideContigs(numErroneousPair, numThread);

		if (divisionMode != SWITCH) {
			for (unsigned i = 0; i < numNode; ++i) {
				sumExpectedLink[i].resize(node[i].numContig - 1, 0.0);
				for (unsigned j = 0; j < node[i].numContig - 1; ++j) {
					sumExpectedLink[i][j] += this->calcExpectedLink((*longReadLibraryMT)[0].getAverageCoverage(), node[i].contig[j].end, node[i].length - node[i].contig[j + 1].start, node[i].contig[j + 1].start - node[i].contig[j].end);
				}
			}
		}
	}


    cerr << "dividing low coverage links..." << endl;
	this->divideErroneousLink(numErroneousPair, numSpanningPair, sumExpectedLink, errorLink, minLink, divisionMode, numThread, maxGapSize);

	this->setTargetLibraryIndex(currentLibraryIndex);
}

void PairedDBG::makeScaffold(void)
{

    long numNewContig = 0;
    long numNewNode = 0;
    long newContigPoolSize = 0;
	long currentClusterID;
    vector<GraphLayout> include, candidate;
    include.reserve(10000);
    candidate.reserve(10000);

    cerr << "scaffolding..." << endl;

    FILE *scaffoldFP = platanus::makeTemporaryFile();

    for (long i = 0; i < numNode; ++i) {
        if (node[i].state & (SC_INC | SC_REP | SC_DEL)) continue;
        include.clear();
        candidate.clear();
        include.push_back(GraphLayout());
        include[0].start =  include[0].distance = 0;
        include[0].end = node[i].length;
        include[0].id = i + 1;
        numNewContig = node[i].numContig;
        node[i].state |= SC_INC;
		currentClusterID = node[i].clusterID;

        for (long j = 0; j < node[i].numEdge; ++j) {
            long tmpNodeID = abs(node[i].edge[j].end) - 1;
            if ((node[tmpNodeID].state & SC_INC) && !(node[tmpNodeID].state & SC_REP)) continue;

			if (currentClusterID != node[tmpNodeID].clusterID && currentClusterID * node[tmpNodeID].clusterID != 0) continue;

            candidate.push_back(GraphLayout());
            unsigned long long candidateEnd = candidate.size() - 1;
            if (node[i].edge[j].direction > 0) {
                candidate[candidateEnd].start = include[0].end + node[i].edge[j].length;
                candidate[candidateEnd].end = candidate[candidateEnd].start + node[tmpNodeID].length;
            } else {
                candidate[candidateEnd].end = -(node[i].edge[j].length);
                candidate[candidateEnd].start = candidate[candidateEnd].end - node[tmpNodeID].length;
            }
            candidate[candidateEnd].id = node[i].edge[j].end;
            candidate[candidateEnd].distance = 1;
            candidate[candidateEnd].numLink = node[i].edge[j].numLink;
        }

        while (candidate.size() > 0) {
            long minDistance = candidate[0].distance;
            long minNumLink = candidate[0].start;
            long minCandidateID = 0;
            for (unsigned j = 1; j < candidate.size(); ++j) {
                if (candidate[j].distance < minDistance || (candidate[j].distance == minDistance && abs(candidate[j].start) < minNumLink)) {
                    minDistance = candidate[j].distance;
                    minNumLink = abs(candidate[j].start);
                    minCandidateID = j;
                }
            }

            long tmpNodeID = abs(candidate[minCandidateID].id) - 1;
            if (((node[tmpNodeID].state & SC_INC) && !((node[tmpNodeID].state & SC_INC) & SC_REP)) || (currentClusterID != node[tmpNodeID].clusterID && currentClusterID * node[tmpNodeID].clusterID != 0)) {
                auto candidateIterator = candidate.begin();
                candidateIterator = candidate.erase(candidateIterator + minCandidateID);
                continue;
            }
            unsigned j = 0;
            for (; j < include.size(); ++j) {
                if (candidate[minCandidateID].end <= include[j].start
                 || candidate[minCandidateID].start >= include[j].end
                 || abs(candidate[minCandidateID].start - include[j].end) <= tolerence + this->getScaffoldOverlap(include[j].id, candidate[minCandidateID].id)
                 || abs(candidate[minCandidateID].end - include[j].start) <= tolerence + this->getScaffoldOverlap(candidate[minCandidateID].id, include[j].id))
                    continue;
                break;
            }
            if (j == include.size()) {
                include.push_back(candidate[minCandidateID]);

                GraphNode &newNode = node[abs(include[include.size() - 1].id) - 1];

				if (currentClusterID == 0)
					currentClusterID = newNode.clusterID;

                if (~(newNode.state) & SC_REP) {
                    for (long k = 0; k < newNode.numEdge; ++k) {
                        long tmpNodeID = abs(newNode.edge[k].end) - 1;
                        if ((node[tmpNodeID].state & SC_INC) && !(node[tmpNodeID].state & SC_REP)) continue;

						if (currentClusterID != node[tmpNodeID].clusterID && currentClusterID * node[tmpNodeID].clusterID != 0) continue;

                        GraphLayout tmpLayout;
                        if (include[include.size() - 1].id * newNode.edge[k].direction > 0) {
                            tmpLayout.start = include[include.size() - 1].end + newNode.edge[k].length;
                            tmpLayout.end = tmpLayout.start + node[tmpNodeID].length;
                        }
                        else {
                            tmpLayout.end = include[include.size() - 1].start - newNode.edge[k].length;
                            tmpLayout.start = tmpLayout.end - node[tmpNodeID].length;
                        }
                        tmpLayout.id = include[include.size() - 1].id > 0 ? newNode.edge[k].end : -(newNode.edge[k].end);
                        tmpLayout.distance = include[include.size() - 1].distance + 1;
                        tmpLayout.numLink = newNode.edge[k].numLink;
                        candidate.push_back(tmpLayout);
                    }
                }

                numNewContig += newNode.numContig;
                if (!(newNode.state & SC_REP))
                    newNode.state |= SC_INC;
            }


            auto candidateIterator = candidate.begin();
            candidateIterator = candidate.erase(candidateIterator + minCandidateID);
        }

        sort(include.begin(), include.end());
        long includeSizeMover = include.size();
        long j = 0;
        for (; node[abs(include[j].id) - 1].state & SC_REP; ++j)
            numNewContig -= node[abs(include[j].id) - 1].numContig;
        for (; node[abs(include[includeSizeMover - 1].id) - 1].state & SC_REP; --includeSizeMover)
            numNewContig -= node[abs(include[includeSizeMover - 1].id) - 1].numContig;
        fwrite(&numNewContig, sizeof(long), 1, scaffoldFP);
        long minStart = include[j].start;
        for (; j < includeSizeMover; ++j) {
            long tmpNodeID = abs(include[j].id) - 1;

            node[tmpNodeID].state |= SC_INC;

            include[j].start -= minStart;
            include[j].end -= minStart;

            if (include[j].start != 0) {
                long overlapLength = this->getScaffoldOverlap(include[j-1].id, include[j].id);
                if (overlapLength > 0 && overlapLength + include[j].start - include[j-1].end <= tolerence) {
                    overlapLength = include[j-1].end - include[j].start - overlapLength;
                    for (long k = j; k < includeSizeMover; ++k) {
                        include[k].end += overlapLength;
                        include[k].start += overlapLength;
                    }
                }
            }

            if (include[j].id > 0) {
                for (long k = 0; k < node[tmpNodeID].numContig; ++k) {
                    ScaffoldPart tmpScaffoldPart(node[tmpNodeID].contig[k].id, include[j].start + node[tmpNodeID].contig[k].start, include[j].start + node[tmpNodeID].contig[k].end);
                    fwrite(&tmpScaffoldPart, sizeof(ScaffoldPart), 1, scaffoldFP);
                }
            }
            else {
                for (long k = node[tmpNodeID].numContig - 1; k >= 0; --k) {
                    ScaffoldPart tmpScaffoldPart(-(node[tmpNodeID].contig[k].id), include[j].start + node[tmpNodeID].length - node[tmpNodeID].contig[k].end, include[j].start + node[tmpNodeID].length - node[tmpNodeID].contig[k].start);
                    fwrite(&tmpScaffoldPart, sizeof(ScaffoldPart), 1, scaffoldFP);
                }
            }
        }

        newContigPoolSize += numContig;
        ++numNewNode;
    }
    include.clear();
    candidate.clear();

    for (long i = 0; i < numNode; ++i) {
        if (!(node[i].state & SC_REP)) continue;
        long numContig = node[i].numContig;
        if (node[i].state & SC_INC) {
            for (long j = 0; j < numContig; ++j)
                contigPositionInScaffold[abs(node[i].contig[j].id) - 1].id = 0;
        } else {
            fwrite(&(node[i].numContig), sizeof(long), 1, scaffoldFP);
            for (long j = 0; j < numContig; ++j)
                fwrite(&(node[i].contig[j]), sizeof(ScaffoldPart), 1, scaffoldFP);
            newContigPoolSize += node[i].numContig;
            ++numNewNode;
        }
    }
    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::makeScaffoldCombine(void)
{
    long numNewContig = 0;
    long numNewNode = 0;
    long newContigPoolSize = 0;
	long currentClusterID;
    vector<GraphLayout> include, candidate;
    include.reserve(10000);
    candidate.reserve(10000);

    cerr << "scaffolding..." << endl;

    FILE *scaffoldFP = platanus::makeTemporaryFile();

    for (long i = 0; i < numNode; ++i) {
        if (node[i].state & (SC_INC | SC_REP | SC_DEL)) continue;
        include.clear();
        candidate.clear();
        include.push_back(GraphLayout());
        include[0].start =  include[0].distance = 0;
        include[0].end = node[i].length;
        include[0].id = i + 1;
        numNewContig = node[i].numContig;
        node[i].state |= SC_INC;
		currentClusterID = node[i].clusterID;

        for (long j = 0; j < node[i].numEdge; ++j) {
            long tmpNodeID = abs(node[i].edge[j].end) - 1;
            if ((node[tmpNodeID].state & SC_INC) && !(node[tmpNodeID].state & SC_REP)) continue;

			if (currentClusterID != node[tmpNodeID].clusterID && currentClusterID * node[tmpNodeID].clusterID != 0) continue;

            candidate.push_back(GraphLayout());
            unsigned long long candidateEnd = candidate.size() - 1;
            if (node[i].edge[j].direction > 0) {
                candidate[candidateEnd].start = include[0].end + node[i].edge[j].length;
                candidate[candidateEnd].end = candidate[candidateEnd].start + node[tmpNodeID].length;
            } else {
                candidate[candidateEnd].end = -(node[i].edge[j].length);
                candidate[candidateEnd].start = candidate[candidateEnd].end - node[tmpNodeID].length;
            }
            candidate[candidateEnd].id = node[i].edge[j].end;
            candidate[candidateEnd].distance = 1;
            candidate[candidateEnd].numLink = node[i].edge[j].numLink;
            candidate[candidateEnd].score = node[i].edge[j].score;
        }


        while (candidate.size() > 0) {
            long maxScore = candidate[0].score;
            long maxCandidateID = 0;

			bool repFlag = false;
			if (node[id2Index(candidate[0].id)].state & SC_REP)
				repFlag = true;
			
            for (unsigned j = 1; j < candidate.size(); ++j) {
				if (repFlag && !(node[id2Index(candidate[j].id)].state & SC_REP)) {
                    maxScore = candidate[j].score;
                    maxCandidateID = j;
					repFlag = false;
				}
                else if (candidate[j].score > maxScore) {
					maxScore = candidate[j].score;
                    maxCandidateID = j;
                }
            }


            long tmpNodeID = abs(candidate[maxCandidateID].id) - 1;

            if ((node[tmpNodeID].state & SC_INC) && (!(node[tmpNodeID].state & SC_REP) || (node[tmpNodeID].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)))) {
                auto candidateIterator = candidate.begin();
                candidateIterator = candidate.erase(candidateIterator + maxCandidateID);
                continue;
            }
            unsigned j = 0;
            for (; j < include.size(); ++j) {
				long overlapTolerence = tolerence;
//				long overlapTolerence = std::min(tolerence, std::min(node[tmpNodeID].length, include[j].end - include[j].start)/2);
				if (node[tmpNodeID].state & SC_REP)
					overlapTolerence = 0;

				if (candidate[maxCandidateID].end <= include[j].start
				 || candidate[maxCandidateID].start >= include[j].end
				 || abs(candidate[maxCandidateID].start - include[j].end) <= overlapTolerence + this->getScaffoldOverlap(include[j].id, candidate[maxCandidateID].id)
				 || abs(candidate[maxCandidateID].end - include[j].start) <= overlapTolerence + this->getScaffoldOverlap(candidate[maxCandidateID].id, include[j].id))
					continue;

                break;
            }
            if (j == include.size()) {
                include.push_back(candidate[maxCandidateID]);

                GraphNode &newNode = node[abs(include[include.size() - 1].id) - 1];
                if (!(newNode.state & SC_REP)) {
                    for (long k = 0; k < newNode.numEdge; ++k) {
                        long tmpNodeID = abs(newNode.edge[k].end) - 1;
                        if ((node[tmpNodeID].state & SC_INC) && !(node[tmpNodeID].state & SC_REP)) continue;

                        GraphLayout tmpLayout;
                        if (include[include.size() - 1].id * newNode.edge[k].direction > 0) {
                            tmpLayout.start = include[include.size() - 1].end + newNode.edge[k].length;
                            tmpLayout.end = tmpLayout.start + node[tmpNodeID].length;
                        }
                        else {
                            tmpLayout.end = include[include.size() - 1].start - newNode.edge[k].length;
                            tmpLayout.start = tmpLayout.end - node[tmpNodeID].length;
                        }
                        tmpLayout.id = include[include.size() - 1].id > 0 ? newNode.edge[k].end : -(newNode.edge[k].end);
                        tmpLayout.distance = include[include.size() - 1].distance + 1;
                        tmpLayout.numLink = newNode.edge[k].numLink;
                        tmpLayout.score = newNode.edge[k].score;
                        candidate.push_back(tmpLayout);
                    }
                }

                numNewContig += newNode.numContig;
                if (!(newNode.state & SC_REP))
                    newNode.state |= SC_INC;
            }

            auto candidateIterator = candidate.begin();
            candidateIterator = candidate.erase(candidateIterator + maxCandidateID);
        }

        sort(include.begin(), include.end());
        long includeSizeMover = include.size();
        long j = 0;
        for (; node[abs(include[j].id) - 1].state & SC_REP; ++j)
            numNewContig -= node[abs(include[j].id) - 1].numContig;
        for (; node[abs(include[includeSizeMover - 1].id) - 1].state & SC_REP; --includeSizeMover)
            numNewContig -= node[abs(include[includeSizeMover - 1].id) - 1].numContig;
        fwrite(&numNewContig, sizeof(long), 1, scaffoldFP);
        long minStart = include[j].start;
        for (; j < includeSizeMover; ++j) {
            long tmpNodeID = abs(include[j].id) - 1;

            node[tmpNodeID].state |= SC_INC;

            include[j].start -= minStart;
            include[j].end -= minStart;

            if (include[j].start != 0) {
                long overlapLength = this->getScaffoldOverlap(include[j-1].id, include[j].id);
                if (overlapLength > 0 && overlapLength + include[j].start - include[j-1].end <= tolerence) {
                    overlapLength = include[j-1].end - include[j].start - overlapLength;
                    for (long k = j; k < includeSizeMover; ++k) {
                        include[k].end += overlapLength;
                        include[k].start += overlapLength;
                    }
                }
            }

            if (include[j].id > 0) {
                for (long k = 0; k < node[tmpNodeID].numContig; ++k) {
                    ScaffoldPart tmpScaffoldPart(node[tmpNodeID].contig[k].id, include[j].start + node[tmpNodeID].contig[k].start, include[j].start + node[tmpNodeID].contig[k].end);
                    fwrite(&tmpScaffoldPart, sizeof(ScaffoldPart), 1, scaffoldFP);
                }
            }
            else {
                for (long k = node[tmpNodeID].numContig - 1; k >= 0; --k) {
                    ScaffoldPart tmpScaffoldPart(-(node[tmpNodeID].contig[k].id), include[j].start + node[tmpNodeID].length - node[tmpNodeID].contig[k].end, include[j].start + node[tmpNodeID].length - node[tmpNodeID].contig[k].start);
                    fwrite(&tmpScaffoldPart, sizeof(ScaffoldPart), 1, scaffoldFP);
                }
            }
        }

        newContigPoolSize += numContig;
        ++numNewNode;
    }
    include.clear();
    candidate.clear();

    for (long i = 0; i < numNode; ++i) {
        if (!(node[i].state & SC_REP)) continue;
        long numContig = node[i].numContig;
        if (node[i].state & SC_INC) {
            for (long j = 0; j < numContig; ++j)
                contigPositionInScaffold[abs(node[i].contig[j].id) - 1].id = 0;
        } else {
            fwrite(&(node[i].numContig), sizeof(long), 1, scaffoldFP);
            for (long j = 0; j < numContig; ++j)
                fwrite(&(node[i].contig[j]), sizeof(ScaffoldPart), 1, scaffoldFP);
            newContigPoolSize += node[i].numContig;
            ++numNewNode;
        }
    }
    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::clearContigPreviousParentNodeID()
{
	this->contigPreviousParentNodeID.clear();
}

void PairedDBG::setOppositeBubbleContigIDByEndMatch()
{
	if (this->oppositeBubbleContigID.empty())
		this->oppositeBubbleContigID.resize(this->numContig, 0);

	std::unordered_map<string, std::array<long, 3> > endSeqMap;

	for (unsigned long i = this->node.size() - this->numInputBubbleContig; i < this->node.size(); ++i) {
		if (this->node[i].length >= this->contigMaxK) {
			string endSeq = contig[i].base.substr(0, contigMaxK - 1) + contig[i].base.substr(this->contig[i].base.size() - (contigMaxK - 1), contigMaxK - 1);
			auto mapItr = endSeqMap.find(endSeq);

			if (mapItr == endSeqMap.end()) {
				std::array<long, 3> endSeqInfo;
				endSeqInfo[0] = i + 1;
				endSeqInfo[1] = 0;
				endSeqInfo[2] = 1;
				endSeqMap[endSeq] = endSeqInfo;
			}
			else {
				if ((mapItr->second)[2] == 1)
					(mapItr->second)[1] = i + 1;
				++(mapItr->second)[2];
			}
		}
	}

	for (auto mapItr = endSeqMap.begin(); mapItr != endSeqMap.end(); ++mapItr) {
		if ((mapItr->second)[2] == 2) {
			this->oppositeBubbleContigID[(mapItr->second)[0] - 1] = (mapItr->second)[1];
			this->oppositeBubbleContigID[(mapItr->second)[1] - 1] = (mapItr->second)[0];
		}
	}
}

long PairedDBG::getScoreFromIDPair(long leftNodeID, long rightNodeID)
{
	GraphNode &leftNode = this->node[abs(leftNodeID) - 1];
	char directionFromLeft = sign(leftNodeID);
	char rightStrand = directionFromLeft * sign(rightNodeID);

	for (long i = 0; i < leftNode.numEdge; ++i) {
		if (abs(leftNode.edge[i].end) == abs(rightNodeID) &&
			leftNode.edge[i].direction == directionFromLeft &&
			sign(leftNode.edge[i].end) == rightStrand) {

			return leftNode.edge[i].score;
		}
	}

	return 0;
}

void PairedDBG::insertSizeDistribution(vector<SeqLib>& library, vector<long>& distribution, const long numThread)
{
    long insertLength;
    platanus::Position forwardResult;
    platanus::Position reverseResult;

    for (long threadID = 0; threadID < numThread; ++threadID) {
		unsigned numPairLink;
		unsigned numReadLink;

        rewind(library[threadID].mappedFP);
        while (fread(&numPairLink, sizeof(unsigned), 1, library[threadID].mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, library[threadID].mappedFP);

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, library[threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, library[threadID].mappedFP);

				if (pairIndex > 0)
					continue;

				unsigned long forwardIndex = abs(forwardResult.id) - 1;
				if (contigPositionInScaffold[forwardIndex].id == 0) continue;
				forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[forwardIndex].id : -(contigPositionInScaffold[forwardIndex].id);
				forwardResult.offset = contigPositionInScaffold[forwardIndex].id > 0 ? forwardResult.offset : contig[forwardIndex].length - forwardResult.offset - 1;
				forwardResult.offset += node[abs(forwardResult.id) - 1].contig[contigPositionInScaffold[forwardIndex].offset].start;

				unsigned long reverseIndex = abs(reverseResult.id) - 1;
				if (contigPositionInScaffold[reverseIndex].id == 0) continue;
				reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[reverseIndex].id : -(contigPositionInScaffold[reverseIndex].id);
				reverseResult.offset = contigPositionInScaffold[reverseIndex].id > 0 ? reverseResult.offset : contig[reverseIndex].length - reverseResult.offset - 1;
				reverseResult.offset += node[abs(reverseResult.id) - 1].contig[contigPositionInScaffold[reverseIndex].offset].start;

				if (forwardResult.id != -reverseResult.id) continue;

				if (forwardResult.id > 0 && forwardResult.offset < reverseResult.offset)
					insertLength = static_cast<long>(reverseResult.offset - forwardResult.offset + 1);
				else if (reverseResult.id > 0 && reverseResult.offset < forwardResult.offset)
					insertLength = static_cast<long>(forwardResult.offset - reverseResult.offset + 1);
				else
					continue;

				if (insertLength >= static_cast<long>(distribution.size())) {
					distribution.resize(insertLength + 1, static_cast<long>(0));
				}
				++distribution[insertLength];
			}

			for (unsigned pairIndex = 0; pairIndex < numReadLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, library[threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, library[threadID].mappedFP);
			}
        }
    }
}

void PairedDBG::updateInsertLengthFP(vector<SeqLib>& lib, const long numThread)
{
    platanus::Position forwardResult;
    platanus::Position reverseResult;

	fseek(lib[0].insertLengthFP, 0, SEEK_END);

    for (long threadID = 0; threadID < numThread; ++threadID) {
		unsigned numPairLink;
		unsigned numReadLink;

        rewind(lib[threadID].mappedFP);
        while (fread(&numPairLink, sizeof(unsigned), 1, lib[threadID].mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, lib[threadID].mappedFP);

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, lib[threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, lib[threadID].mappedFP);

				if (pairIndex > 0)
					continue;

				unsigned long forwardIndex = abs(forwardResult.id) - 1;
				if (contigPositionInScaffold[forwardIndex].id == 0) continue;

				unsigned long reverseIndex = abs(reverseResult.id) - 1;
				if (contigPositionInScaffold[reverseIndex].id == 0) continue;

				if (forwardIndex == reverseIndex) continue;

				forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[forwardIndex].id : -(contigPositionInScaffold[forwardIndex].id);
				forwardResult.offset = contigPositionInScaffold[forwardIndex].id > 0 ? forwardResult.offset : contig[forwardIndex].length - forwardResult.offset - 1;
				forwardResult.offset += node[abs(forwardResult.id) - 1].contig[contigPositionInScaffold[forwardIndex].offset].start;

				reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[reverseIndex].id : -(contigPositionInScaffold[reverseIndex].id);
				reverseResult.offset = contigPositionInScaffold[reverseIndex].id > 0 ? reverseResult.offset : contig[reverseIndex].length - reverseResult.offset - 1;
				reverseResult.offset += node[abs(reverseResult.id) - 1].contig[contigPositionInScaffold[reverseIndex].offset].start;

				if (forwardResult.id != -reverseResult.id) continue;

				long insertLength = abs(forwardResult.offset - reverseResult.offset) ;
				fwrite(&insertLength, sizeof(long), 1, lib[0].insertLengthFP);
			}

			for (unsigned pairIndex = 0; pairIndex < numReadLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, lib[threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, lib[threadID].mappedFP);
			}
        }
    }
}

long PairedDBG::solveSimpleCrossStructureCoverage(const double minCoverage, const long minLength, const double sameCoverageThreshold, const double diffCoverageThreshold, const long numThread)
{
	vector<vector<GraphPath> > pathBuffer(numThread);

    # pragma omp parallel for schedule(static, 1)
	for (long t = 0; t < numThread; ++t) {
		vector<long> nodeIDBuffer;
		vector<long> nodeIDBufferForCheck;
		std::array<std::array<long, 2>, 2 > externalNodeID;

		for (long centerNodeIndex = t; centerNodeIndex < numNode; centerNodeIndex += numThread) {
			if (calcNodeCoverage(this->node[centerNodeIndex]) < minCoverage || this->node[centerNodeIndex].length < minLength)
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				char initialDirection = 2*i - 1;
				getOverlappedNode(centerNodeIndex, initialDirection, nodeIDBuffer);
				if (nodeIDBuffer.size() != 2)
					break;

				long j = 0;
				for (; j < 2; ++j) {
					externalNodeID[i][j] = nodeIDBuffer[j];
					if (calcNodeCoverage(this->node[id2Index(externalNodeID[i][j])]) < minCoverage || this->node[id2Index(externalNodeID[i][j])].length < minLength)
						break;

					getOverlappedNode(id2Index(externalNodeID[i][j]), sign(externalNodeID[i][j]), nodeIDBufferForCheck);
					if (nodeIDBufferForCheck.size() >= 2)
						break;
				}
				if (j < 2)
					break;
			}
			if (i < 2)
				continue;

			std::array<char, 2> inconsistentHaplotypeFlag;
			inconsistentHaplotypeFlag.fill(0);
			for (long j = 0; j < 2; ++j) {
				GraphNode &leftNode = this->node[id2Index(externalNodeID[0][0])];
				GraphNode &rightNode = this->node[id2Index(externalNodeID[1][j])];
				GraphNode &leftAltNode = this->node[id2Index(externalNodeID[0][1])];
				GraphNode &rightAltNode = this->node[id2Index(externalNodeID[1][abs(j - 1)])];

				if (!(judgeCloseness(calcNodeCoverage(leftNode), calcNodeCoverage(rightNode), sameCoverageThreshold))  ||
					!(judgeCloseness(calcNodeCoverage(leftAltNode), calcNodeCoverage(rightAltNode), sameCoverageThreshold)) ||
					judgeCloseness(calcNodeCoverage(leftNode), calcNodeCoverage(leftAltNode), diffCoverageThreshold) ||
					judgeCloseness(calcNodeCoverage(leftNode), calcNodeCoverage(rightAltNode), diffCoverageThreshold) ||
					judgeCloseness(calcNodeCoverage(leftAltNode), calcNodeCoverage(rightNode), diffCoverageThreshold) ||
					judgeCloseness(calcNodeCoverage(rightNode), calcNodeCoverage(rightAltNode), diffCoverageThreshold)) {

					inconsistentHaplotypeFlag[j] = 1;
				}
			}

			if (inconsistentHaplotypeFlag[0] == inconsistentHaplotypeFlag[1])
				continue;

			char crossFlag;
			if (inconsistentHaplotypeFlag[1])
				crossFlag = 0;
			else if (inconsistentHaplotypeFlag[0])
				crossFlag = 1;
			else
				continue;
			
			for (long j = 0; j < 2; ++j) {
				pathBuffer[t].emplace_back(centerNodeIndex, 3);
				pathBuffer[t].back().nodeID[0] = externalNodeID[0][j];
				pathBuffer[t].back().nodeID[1] = centerNodeIndex + 1;
				pathBuffer[t].back().nodeID[2] = externalNodeID[1][(j + crossFlag) % 2];
			}
		}
	}

	vector<GraphPath> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), GraphPathSelfIDLess());

	long numSolvedCross = remakeGraphAccordingToPathPair(mergedPathBuffer);

	return numSolvedCross;
}

void PairedDBG::solveSimpleCrossStructureCoverageIterative(const bool clusterFlag, const long numThread)
{
	const double MIN_COVERAGE = 20.0;
//	const long MIN_LENGTH = this->contigMaxK * 4 + 1;
	const long MIN_LENGTH = 0;
	const double SAME_COVERAGE_THRESHOLD = 1.25;
	const double DIFF_COVERAGE_THRESHOLD = 1.50;

	const unsigned long currentMinLink = this->minLink;

	unsigned currentMode = getMode();
	this->setMode(OVERLAP_MODE);

	long total = 0;
	long num;
	long iteration = 0;
	cerr << endl << "solving simple cross-structurs in a de Bruijn graph according to coverage depth ..." << endl;
	do {
		setMinLink(1);
		makeGraph(numThread);

		if (clusterFlag)
			deleteChimericEdge();
		else
			deleteInterOTUEdge();

		setMinLink(currentMinLink);
		num = solveSimpleCrossStructureCoverage(MIN_COVERAGE, MIN_LENGTH, SAME_COVERAGE_THRESHOLD, DIFF_COVERAGE_THRESHOLD, numThread);

		total += num;
		++iteration;
		cerr << "NUM_SOLVED_CROSSES = " << num << endl;
	} while (num > 0 && iteration < MAX_ITERATION_OF_CROSS_SOLUTION);
	cerr << "TOTAL_NUM_NODES_IN_SOLVED_CROSSES=" << total << endl << endl;

	setMode(currentMode);
}



void PairedDBG::deleteChimericEdge()
{
    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing edges bridging different clusters..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
            GraphNode &node2 = this->node[abs(node1.edge[edgeID].end) - 1];
			if (node1.clusterID != node2.clusterID && node1.clusterID * node2.clusterID != 0) {
				ids.push_back(nodeID + 1);
				ids.push_back(node[nodeID].edge[edgeID].end);
				++numDelete;
			}
        }
    }

    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}


void PairedDBG::setMgcEmpDistr()
{
	Mgc (4).calcDistributionFromContig(contig, MIN_LENGTH_FOR_MGC, mgcEmpDistr);
}


void PairedDBG::setCoverageEmpDistr(const string &kmerOccurrenceFile)
{
    unsigned long long kmerLength = platanus::getKmerLengthFromBinary(kmerOccurrenceFile);

    if (kmerLength <= 32) {
        Counter<Kmer31> counter(kmerLength);
        counter.readOccurrenceTableBinary(kmerOccurrenceFile);
		counter.getCoverageDiffRateDistributionFromContig(contig, covEmpDistr, MIN_LENGTH_FOR_COV_ESTIMATE);
    } else if (kmerLength <= 64) {
        Counter<KmerN<Binstr63> > counter(kmerLength);
        counter.readOccurrenceTableBinary(kmerOccurrenceFile);
		counter.getCoverageDiffRateDistributionFromContig(contig, covEmpDistr, MIN_LENGTH_FOR_COV_ESTIMATE);
    } else if (kmerLength <= 96) {
        Counter<KmerN<Binstr95> > counter(kmerLength);
        counter.readOccurrenceTableBinary(kmerOccurrenceFile);
		counter.getCoverageDiffRateDistributionFromContig(contig, covEmpDistr, MIN_LENGTH_FOR_COV_ESTIMATE);
    } else if (kmerLength <= 128) {
        Counter<KmerN<Binstr127> > counter(kmerLength);
        counter.readOccurrenceTableBinary(kmerOccurrenceFile);
		counter.getCoverageDiffRateDistributionFromContig(contig, covEmpDistr, MIN_LENGTH_FOR_COV_ESTIMATE);
    } else if (kmerLength <= 160) {
        Counter<KmerN<Binstr159> > counter(kmerLength);
        counter.readOccurrenceTableBinary(kmerOccurrenceFile);
		counter.getCoverageDiffRateDistributionFromContig(contig, covEmpDistr, MIN_LENGTH_FOR_COV_ESTIMATE);
    } else {
        Counter<KmerN<binstr_t> > counter(kmerLength);
        counter.readOccurrenceTableBinary(kmerOccurrenceFile);
		counter.getCoverageDiffRateDistributionFromContig(contig, covEmpDistr, MIN_LENGTH_FOR_COV_ESTIMATE);
    }
}


bool PairedDBG::setClusterID(const platanus::Contig &contig, const string ClusterTsvFilename)
{
	std::ifstream ifs(ClusterTsvFilename);
	if (!ifs) {
		cerr << "Warning: " << ClusterTsvFilename << " cannot be open!" << endl;
		return false;
	}

	clusterID.assign(numContig, 0);

	string oneLine;
	string seqName;
	long ID;
	std::ostringstream outStream;
    while (ifs && getline(ifs, oneLine)) {
		std::istringstream ss(oneLine);
		ss >> seqName >> ID;

		auto it = contig.nameIndex.find(seqName);
		if (it != contig.nameIndex.end()) {
			clusterID[it->second] = ID;
		}
	}
	ifs.close();

	updateNodeClusterID();

	return true;
}


void PairedDBG::setClusterIDfromSeqName(const platanus::Contig &contig)
{
	static const string CLUSTER_PREFIX = "cluster";

	for (long i = 0; i < this->numContig; ++i) {
		long ID = 0;
		size_t pos = contig.name[i].find(CLUSTER_PREFIX);
		if (pos != std::string::npos || pos + CLUSTER_PREFIX.size() + 1 < contig.name[i].size())
			ID = atol(contig.name[i].substr(pos + CLUSTER_PREFIX.size()).c_str());

		clusterID[i] = ID;
	}
}


void PairedDBG::deleteInterOTUEdge()
{
	Mgc mgc(4);
    vector<long> ids;
	const unsigned minLengthForMgc = MIN_LENGTH_FOR_MGC;
    unsigned long long numDelete = 0;

    std::cerr << "removing inter-species edge candidates..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
            GraphNode &node2 = this->node[abs(node1.edge[edgeID].end) - 1];

			double cov1 = this->calcNodeCoverage(node1);
			double cov2 = this->calcNodeCoverage(node2);
			double covRatio = (cov1 > cov2 ? cov1/cov2 : cov2/cov1);

			vector<char> scaffoldSeq1,  scaffoldSeq2;
			this->node2seq(node1, scaffoldSeq1);
			this->node2seq(node2, scaffoldSeq2);

			vector<char> seq1, seq2;
			this->node2seq(node1, seq1);
			this->node2seq(node2, seq2);

			double p;

			if (seq1.size() >= minLengthForMgc && seq2.size() >= minLengthForMgc) {
				mgc.clear();
				mgc.addSeq(seq1);
				mgc.addSeq(seq2);
				mgc.kmerStat();
				mgc.setMatrix(false);
				double kmerFreqDistance =  mgc.getDistance(0, 1);

				p = covEmpDistr.probGreater(covRatio) * mgcEmpDistr.probGreater(kmerFreqDistance);
			}
			else {
				p = covEmpDistr.probGreater(covRatio) * 0.5;
			}

			if (p < condProbThreshold) {
				ids.push_back(nodeID + 1);
				ids.push_back(node[nodeID].edge[edgeID].end);
				++numDelete;
			}
        }
    }

    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}


void PairedDBG::splitResultSeqForClusters(vector<vector<string > > &seq,  vector<vector<string> > &name)
{
	this->updateNodeClusterID();

	long numCluster = 0;
    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if ((node[nodeIndex].state & SC_DEL) || resultSeq[nodeIndex].redundantFlag)
			continue;

		if (numCluster < node[nodeIndex].clusterID + 1)
			numCluster = node[nodeIndex].clusterID + 1;
	}
	seq.resize(numCluster);
	name.resize(numCluster);

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if ((node[nodeIndex].state & SC_DEL) || resultSeq[nodeIndex].redundantFlag)
			continue;

		seq[node[nodeIndex].clusterID].push_back(resultSeq[nodeIndex].seq);

		std::ostringstream oss;
		oss << resultSeq[nodeIndex].name << "_cluster" << node[nodeIndex].clusterID;
		resultSeq[nodeIndex].name = oss.str();
		name[node[nodeIndex].clusterID].push_back(resultSeq[nodeIndex].name);
    }
}


void PairedDBG::storeGraphLinkFromTagReadPair(vector<GraphLinkWithFlag> &graphLinkPool, const bool scaffoldFlag, long numThread)
{
	const long MAX_INTERVAL_LENGTH = 10000;
	const long MIN_NUM_TO_CALCULATE_MEAN_LENGTH = 20;
	const long DEFAULT_MEAN_LENGTH = 50000;

	if (!(this->mode & TAG_SCAFFOLD_MODE) || tagLibraryMT == NULL)
		return;

	cerr << "calculating links based on bardodes information ..." << endl;

	vector<vector<MappedTagInfo> > mappedTagInfoBuffer(numThread);

    omp_set_num_threads(numThread);
	# pragma omp parallel for schedule(static, 1)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		platanus::Position forwardResult;
		platanus::Position reverseResult;
		int tagID;

		rewind((*tagLibraryMT)[threadID].mappedFP);
		while (fread(&forwardResult, sizeof(platanus::Position), 1, (*tagLibraryMT)[threadID].mappedFP)) {
			fread(&reverseResult, sizeof(platanus::Position), 1, (*tagLibraryMT)[threadID].mappedFP);
			fread(&tagID, sizeof(int), 1, (*tagLibraryMT)[threadID].mappedFP);

			long i = id2Index(forwardResult.id);
			if (contigPositionInScaffold[i].id == 0) continue;
			forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[i].id : -(contigPositionInScaffold[i].id);
			forwardResult.offset = contigPositionInScaffold[i].id > 0 ? forwardResult.offset : contig[i].length - forwardResult.offset - 1;
			forwardResult.offset += node[id2Index(forwardResult.id)].contig[contigPositionInScaffold[i].offset].start;

			long j = id2Index(reverseResult.id);
			if (contigPositionInScaffold[j].id == 0) continue;
			reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[j].id : -(contigPositionInScaffold[j].id);
			reverseResult.offset = contigPositionInScaffold[j].id > 0 ? reverseResult.offset : contig[j].length - reverseResult.offset - 1;
			reverseResult.offset += node[id2Index(reverseResult.id)].contig[contigPositionInScaffold[j].offset].start;

			if (forwardResult.id == 0 && reverseResult.id == 0)
				continue;

			if (forwardResult.id == -(reverseResult.id)) {
				# pragma omp critical (push)
				{
					mappedTagInfoBuffer[tagID % numThread].emplace_back(tagID, (forwardResult.offset + reverseResult.offset)/2, abs(forwardResult.id));
				}
			}
			else {
				if (forwardResult.id == reverseResult.id)
					continue;

				if (forwardResult.id != 0) {
					# pragma omp critical (push)
					{
						mappedTagInfoBuffer[tagID % numThread].emplace_back(tagID, forwardResult.offset, abs(forwardResult.id));
					}
				}

				if (reverseResult.id != 0) {
					# pragma omp critical (push)
					{
						mappedTagInfoBuffer[tagID % numThread].emplace_back(tagID, reverseResult.offset, abs(reverseResult.id));
					}
				}
			}

		}
	}

	vector<vector<MappedTagInfo> > moleculePositionBuffer(numThread);
	double weightedSumForMean = 0.0;
	double weightedNumForMean = 0.0;
	long numForMean = 0;

	# pragma omp parallel for schedule(static, 1) reduction(+: weightedSumForMean, weightedNumForMean, numForMean)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		sort(mappedTagInfoBuffer[threadID].begin(), mappedTagInfoBuffer[threadID].end());
		mappedTagInfoBuffer[threadID].emplace_back(0, 0, 0);

		unsigned long beginIndex = 0;
		for (unsigned long i = 0; i < mappedTagInfoBuffer[threadID].size() - 1; ++i) {
			if (mappedTagInfoBuffer[threadID][i + 1].tagID != mappedTagInfoBuffer[threadID][i].tagID ||
				mappedTagInfoBuffer[threadID][i + 1].position.id != mappedTagInfoBuffer[threadID][i].position.id ||
				mappedTagInfoBuffer[threadID][i + 1].position.offset - mappedTagInfoBuffer[threadID][i].position.offset > MAX_INTERVAL_LENGTH) {

				if (i + 1 - beginIndex >= this->minTagReadToLoad) {
					platanus::Position &leftPosition = mappedTagInfoBuffer[threadID][beginIndex].position;
					platanus::Position &rightPosition = mappedTagInfoBuffer[threadID][i].position;

					moleculePositionBuffer[threadID].emplace_back(mappedTagInfoBuffer[threadID][i].tagID, leftPosition.offset, leftPosition.id);
					moleculePositionBuffer[threadID].emplace_back(mappedTagInfoBuffer[threadID][i].tagID, rightPosition.offset, rightPosition.id);

					if (leftPosition.offset > MAX_INTERVAL_LENGTH && node[id2Index(rightPosition.id)].length - rightPosition.offset > MAX_INTERVAL_LENGTH) {
						weightedSumForMean += ((double)rightPosition.offset - leftPosition.offset) * (rightPosition.offset - leftPosition.offset);
						weightedNumForMean += (double)rightPosition.offset - leftPosition.offset;
						++numForMean;
					}
				}

				beginIndex = i + 1;
			}
		}
	}


	double weightedMeanMoleculeLength = DEFAULT_MEAN_LENGTH;

	if (numForMean >= MIN_NUM_TO_CALCULATE_MEAN_LENGTH)
		weightedMeanMoleculeLength = weightedSumForMean / weightedNumForMean;

	if (numForMean >= MIN_NUM_TO_CALCULATE_MEAN_LENGTH || (*tagLibraryMT)[0].getAverageInsSize() == 0)
		(*tagLibraryMT)[0].setAverageInsSize(weightedMeanMoleculeLength);
	else
		weightedMeanMoleculeLength = (*tagLibraryMT)[0].getAverageInsSize();

	cerr << "LENGTH_WEIGHTED_MEAN_BARCODED_MOLECULE_LENGTH =  " << weightedMeanMoleculeLength << endl;

	if (!scaffoldFlag)
		return;


	# pragma omp parallel for schedule(static, 1)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		moleculePositionBuffer[threadID].emplace_back(0, 0, 0);
		vector<platanus::Position> endPositionBuffer;

		for (unsigned long i = 0; i < moleculePositionBuffer[threadID].size() - 1; i += 2) {

			if (i > 0 && (moleculePositionBuffer[threadID][i].tagID != moleculePositionBuffer[threadID][i - 1].tagID)) {
				for (long j = 0; j < static_cast<long>(endPositionBuffer.size()) - 1; ++j) {
					for (long k = j + 1; k < endPositionBuffer.size(); ++k) {
						if (abs(endPositionBuffer[j].id) == abs(endPositionBuffer[k].id))
							continue;

						platanus::Position &forwardResult = endPositionBuffer[j];
						platanus::Position &reverseResult = endPositionBuffer[k];
						GraphLinkWithFlag graphLink;

						if (abs(forwardResult.id) < abs(reverseResult.id)) {
							graphLink.id1 = forwardResult.id;
							graphLink.offset1 = forwardResult.id > 0 ? node[id2Index(forwardResult.id)].numContig - 1 : 0;
							graphLink.id2 = -(reverseResult.id);
							graphLink.offset2 = reverseResult.id > 0 ? node[id2Index(reverseResult.id)].numContig - 1 : 0;
						}
						else {
							graphLink.id1 = reverseResult.id;
							graphLink.offset1 =  reverseResult.id > 0 ? node[id2Index(reverseResult.id)].numContig - 1 : 0;
							graphLink.id2 = -(forwardResult.id);
							graphLink.offset2 =  forwardResult.id > 0 ? node[id2Index(forwardResult.id)].numContig - 1 : 0;
						}

						graphLink.gap = weightedMeanMoleculeLength;
						long gapMinus;
						if (forwardResult.id > 0) {
							gapMinus = node[forwardResult.id - 1].length - forwardResult.offset;
						} else {
							gapMinus = forwardResult.offset + 1;
						}
						graphLink.gap -= gapMinus;

						if (reverseResult.id > 0) {
							if (abs(forwardResult.id) == abs(reverseResult.id)) continue;
							graphLink.gap -= node[reverseResult.id-1].length - reverseResult.offset;
						} else {
							if (abs(forwardResult.id) == abs(reverseResult.id)) continue;
							graphLink.gap -= reverseResult.offset + 1;
						}

						# pragma omp critical (push)
						{
							graphLinkPool.push_back(graphLink);
						}
					}
				}

				endPositionBuffer.clear();
			}

			const platanus::Position &leftPosition = moleculePositionBuffer[threadID][i].position;
			const platanus::Position &rightPosition = moleculePositionBuffer[threadID][i + 1].position;

			if (rightPosition.offset - leftPosition.offset < this->minTagIslandLengthToScaffold)
				continue;

			if (rightPosition.offset < node[id2Index(leftPosition.id)].length - leftPosition.offset) {
//				if (rightPosition.offset < weightedMeanMoleculeLength && node[id2Index(rightPosition.id)].length - rightPosition.offset > MAX_INTERVAL_LENGTH) {
				if (rightPosition.offset < weightedMeanMoleculeLength) {
					endPositionBuffer.push_back(rightPosition);
					endPositionBuffer.back().id *= -1;
//std::cout << "island: " << "(" << leftPosition.id << ", " << leftPosition.offset  << ") (" << rightPosition.id << ", " << rightPosition.offset  << ") length " << rightPosition.offset - leftPosition.offset << std::endl;
				}
			}
			else {
//				if (node[id2Index(leftPosition.id)].length - leftPosition.offset < weightedMeanMoleculeLength && leftPosition.offset > MAX_INTERVAL_LENGTH) {
				if (node[id2Index(leftPosition.id)].length - leftPosition.offset < weightedMeanMoleculeLength) {
					endPositionBuffer.push_back(leftPosition);
//std::cout << "island: " << "(" << leftPosition.id << ", " << leftPosition.offset  << ") (" << rightPosition.id << ", " << rightPosition.offset  << ") length " << rightPosition.offset - leftPosition.offset << std::endl;
				}
			}
		}
	}
}


void PairedDBG::estimateTaggedMoleculeLength(const long numThread)
{
	if (tagLibraryMT == NULL)
		return;

    vector<GraphLinkWithFlag> graphLinkPool;
	storeGraphLinkFromTagReadPair(graphLinkPool, false, numThread);
}


void PairedDBG::dumpAllEdges(const string &outputFilename)
{
//    std::ofstream out(outputFilename.c_str());

    for (long nodeIndex = 0 ; nodeIndex < numNode; ++nodeIndex) {
        GraphNode &nodeRef = this->node[nodeIndex];
        for (long edgeIndex = 0; edgeIndex < nodeRef.numEdge; ++edgeIndex) {
			GraphEdge &edgeRef = nodeRef.edge[edgeIndex];

			if (edgeRef.direction > 0) {
				if (nodeRef.state & SC_REP)
					std::cout << "R ";
				std::cout << nodeIndex + 1;
				std::cout << "(";

				for (long contigIndex = 0; contigIndex < nodeRef.numContig; ++contigIndex) {
					if (nodeRef.contig[contigIndex].id < 0)
						std::cout << "-";
					std::cout << contigName[id2Index(nodeRef.contig[contigIndex].id)] << ", ";
				}
				std::cout << ") " << edgeRef.length << ", " << edgeRef.numLink << " ";

				GraphNode &endNodeRef = this->node[id2Index(edgeRef.end)];
				if (endNodeRef.state & SC_REP)
					std::cout << "R ";
				std::cout << edgeRef.end;
				std::cout << "(";

				if (edgeRef.end > 0) {
					for (long contigIndex = 0; contigIndex < endNodeRef.numContig; ++contigIndex) {
						if (endNodeRef.contig[contigIndex].id < 0)
							std::cout << "-";
						std::cout << contigName[id2Index(endNodeRef.contig[contigIndex].id)] << ", ";
					}
				}
				else {
					for (long contigIndex = 0; contigIndex < endNodeRef.numContig; ++contigIndex) {
						if (endNodeRef.contig[endNodeRef.numContig - contigIndex - 1].id > 0)
							std::cout << "-";
						std::cout << contigName[id2Index(endNodeRef.contig[endNodeRef.numContig - contigIndex - 1].id)] << ", ";
					}
				}
			}
			else {
				GraphNode &endNodeRef = this->node[id2Index(edgeRef.end)];
				if (endNodeRef.state & SC_REP)
					std::cout << "R ";
				std::cout << edgeRef.end;
				std::cout << "(";

				if (edgeRef.end > 0) {
					for (long contigIndex = 0; contigIndex < endNodeRef.numContig; ++contigIndex) {
						if (endNodeRef.contig[contigIndex].id < 0)
							std::cout << "-";
						std::cout << contigName[id2Index(endNodeRef.contig[contigIndex].id)] << ", ";
					}
				}
				else {
					for (long contigIndex = 0; contigIndex < endNodeRef.numContig; ++contigIndex) {
						if (endNodeRef.contig[endNodeRef.numContig - contigIndex - 1].id > 0)
							std::cout << "-";
						std::cout << contigName[id2Index(endNodeRef.contig[endNodeRef.numContig - contigIndex - 1].id)] << ", ";
					}
				}
				std::cout << ") " << edgeRef.length << ", " << edgeRef.numLink << " ";
				if (nodeRef.state & SC_REP)
					std::cout << "R ";
				std::cout << nodeIndex + 1;
				std::cout << "(";

				for (long contigIndex = 0; contigIndex < nodeRef.numContig; ++contigIndex) {
					if (nodeRef.contig[contigIndex].id < 0)
						std::cout << "-";
					std::cout << contigName[id2Index(nodeRef.contig[contigIndex].id)] << ", ";
				}
			}
			std::cout << std::endl;
		}
	}

//    out.close();
}


long PairedDBG::deleteConflictingEdgeToSameNode(const long numThread)
{
	const double DIFF_RATE_THRESHOLD = 2.0;
    vector<long> ids;
    long numDelete = 0;

    omp_set_num_threads(numThread);

    #pragma omp parallel for schedule(dynamic) reduction(+: numDelete)
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {
                const GraphEdge &edge1 = (node[nodeID].edge[edgeID1]);
                const GraphEdge &edge2 = (node[nodeID].edge[edgeID2]);

				if (abs(edge1.end) == abs(edge2.end)) {
					if (judgeCloseness(edge1.numLink, edge2.numLink, DIFF_RATE_THRESHOLD)) {
						#pragma omp critical (delete_edge)
						{
							ids.push_back(nodeID + 1);
							ids.push_back(edge1.end);
							ids.push_back(nodeID + 1);
							ids.push_back(edge2.end);
						}
						numDelete += 2;
					}
					else if (edge1.numLink < edge2.numLink) {
						#pragma omp critical (delete_edge)
						{
							ids.push_back(nodeID + 1);
							ids.push_back(edge1.end);
						}
						++numDelete;
					}
					else if (edge2.numLink < edge1.numLink) {
						#pragma omp critical (delete_edge)
						{
							ids.push_back(nodeID + 1);
							ids.push_back(edge2.end);
						}
						++numDelete;
					}
                }
            }
        }
    }
    this->deleteEdges(ids);
    return numDelete;
}


void PairedDBG::markRedundantResultSeq(const long numThread, const bool clusterID)
{
	unsigned long prefixLength = this->contigMaxK;

	for (unsigned long seqIndex = 0; seqIndex < resultSeq.size(); ++seqIndex) {
		if (resultSeq[seqIndex].seq.empty() || (node[seqIndex].state & SC_DEL))
			continue;

		if (prefixLength > resultSeq[seqIndex].seq.size())
			prefixLength = resultSeq[seqIndex].seq.size();
	}

	std::unordered_map<string, vector<long> > prefixToSeqIndex;

	for (unsigned long seqIndex = 0; seqIndex < resultSeq.size(); ++seqIndex) {
		auto mapItr = prefixToSeqIndex.find(resultSeq[seqIndex].seq.substr(0, prefixLength));
		if (mapItr == prefixToSeqIndex.end())
			prefixToSeqIndex[resultSeq[seqIndex].seq.substr(0, prefixLength)] = std::vector<long>(1, seqIndex);
		else
		 	mapItr->second.push_back(seqIndex);
	}
	
    omp_set_num_threads(numThread);
	vector<vector<char> > redundantFlag(numThread);

	# pragma omp parallel for schedule(static, 1)
	for (unsigned threadID = 0; threadID < numThread; ++threadID) {
		redundantFlag[threadID].assign(resultSeq.size(), false);

		for (unsigned long seqIndex = threadID; seqIndex < resultSeq.size(); seqIndex += numThread) {
			if (resultSeq[seqIndex].seq.empty() || (node[seqIndex].state & SC_DEL))
				continue;

			for (unsigned strandIndex = 0; strandIndex < 2; ++strandIndex) {
				string targetSeq;
				if (strandIndex == 0) {
					targetSeq = resultSeq[seqIndex].seq;
				}
				else {
					targetSeq.resize(resultSeq[seqIndex].seq.size());
					for (unsigned long baseIndex = 0; baseIndex < resultSeq[seqIndex].seq.size(); ++baseIndex) {
						if (resultSeq[seqIndex].seq[baseIndex] < 4)
							targetSeq[targetSeq.size() - baseIndex - 1] = 3^(resultSeq[seqIndex].seq[baseIndex]);
						else
							targetSeq[targetSeq.size() - baseIndex - 1] = 4;
					}
				}

				for (unsigned long baseIndex = 0; baseIndex < targetSeq.size() - prefixLength + 1; ++baseIndex) {
					auto mapItr = prefixToSeqIndex.find(targetSeq.substr(baseIndex, prefixLength));
					if (mapItr == prefixToSeqIndex.end())
						continue;

					for (auto vecItr = mapItr->second.begin(); vecItr != mapItr->second.end(); ++vecItr) {
						if (*vecItr != seqIndex &&
							!(node[*vecItr].state & SC_DEL) &&
							(targetSeq.size() > resultSeq[*vecItr].seq.size() || (targetSeq.size() == resultSeq[*vecItr].seq.size() && seqIndex < *vecItr)) &&
							targetSeq.size() - baseIndex >= resultSeq[*vecItr].seq.size() &&
							std::search(resultSeq[*vecItr].seq.begin()+prefixLength, resultSeq[*vecItr].seq.end(), targetSeq.begin()+baseIndex+prefixLength, targetSeq.begin()+baseIndex+resultSeq[*vecItr].seq.size()) == resultSeq[*vecItr].seq.begin()+prefixLength) {

							if (clusterID && (node[*vecItr].clusterID != 0 && node[*vecItr].clusterID != node[seqIndex].clusterID))
								continue;

							redundantFlag[threadID][*vecItr] = true;
						}
					}
				}
			}
		}
	}

	for (unsigned threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned long seqIndex = 0; seqIndex < resultSeq.size(); ++seqIndex)
			if (resultSeq[seqIndex].redundantFlag == false && redundantFlag[threadID][seqIndex]) 
				resultSeq[seqIndex].redundantFlag = true;
	}
}


void PairedDBG::divideChimericClusterNode()
{
    cerr << "dividing chimeric-cluster scaffolds..." << endl;

    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		if (this->node[nodeIndex].state & SC_DEL)
			continue;

		long i;
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.back() = 1;

		long preClusterID = clusterID[id2Index(contigRef[0].id)];
		long preEnd = 1;
		for (i = 1; i < static_cast<long>(contigRef.size()); ++i) {
			if(clusterID[id2Index(contigRef[i].id)] != 0) {
				if (clusterID[id2Index(contigRef[i].id)] != preClusterID) {
					breakpointFlag[preEnd] = 1;
					breakpointFlag[i] = 1;
				}
				preClusterID = clusterID[id2Index(contigRef[i].id)];
				preEnd = i + 1;
			}
		}

		i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;

			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}
