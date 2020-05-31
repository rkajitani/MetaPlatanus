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

#ifndef PAIRED_DBG_H
#define PAIRED_DBG_H
#define STATIC_VERSION

#include "seqlib.h"
#include "mapper.h"
#include "scaffoldGraph.h"
#include "counter.h"
#include "mgc.h"
#include <unordered_map>
#include <vector>
#include <array>


class PairedDBG : public ScaffoldGraph
{
public:
    enum CROSS_RESOLUTION_MODE {LINK, TAG, SCORE};
    enum DIVISION_MODE {SWITCH, GAP, MIS};

protected:
	
    long contigMaxK;
	long minTagIslandLengthToScaffold;
	double heteroCoverage;
	double tolerenceFactor;
	double cutoffLengthFactor;
	std::vector<long> contigPreviousParentNodeID;
	std::vector<long> oppositeBubbleContigID;
	std::vector<std::array<long, 2> > joinedBubbleContigID;
	unsigned mode;
	unsigned long numInputBubbleContig;
	double condProbThreshold;
	platanus::empiricalDistribution mgcEmpDistr;
	platanus::empiricalDistribution covEmpDistr;

	struct MapInfoForGraph
	{
		platanus::Position position;	
		long contigIndex;
		long score;
		
		struct PositionIDLessScoreGreater
		{
			bool operator()(const MapInfoForGraph &a, const MapInfoForGraph &b)
			{
				return (a.position.id != b.position.id) ? (a.position.id < b.position.id) : ((a.score != b.score) ? (a.score > b.score) : (a.contigIndex < b.contigIndex));
			}
		};

		struct PositionIDLess
		{
			bool operator()(const MapInfoForGraph &a, const MapInfoForGraph &b)
			{
				return (a.position.id < b.position.id);
			}
		};

//		bool operator <(const MapInfoForGraph &a) const { return (position < a.position); }
		MapInfoForGraph(): position(0, 0), contigIndex(0), score(0) {}
		MapInfoForGraph(const platanus::Position p, const long c, const long s): position(p), contigIndex(c), score(s) {}
		~MapInfoForGraph() {}
	};

    struct GraphLinkWithFlag : public GraphLink
    {
		bool overlapFlag;
		long score;

        GraphLinkWithFlag(): GraphLink(), overlapFlag(false), score(1) {}
        ~GraphLinkWithFlag() = default;
	};

    struct GraphLinkWithFlagPoolIndex : public GraphLinkPoolIndex
    {
		bool overlapFlag;
		long gap;

        GraphLinkWithFlagPoolIndex(): GraphLinkPoolIndex(), overlapFlag(false), gap(0) {}
        GraphLinkWithFlagPoolIndex(unsigned long idx): GraphLinkPoolIndex(idx), overlapFlag(false), gap(0) {}
        ~GraphLinkWithFlagPoolIndex() = default;
    };

    struct GraphLinkWithFlagPoolIndexGreater
    {
        bool operator() (const GraphLinkWithFlagPoolIndex& link1, const GraphLinkWithFlagPoolIndex& link2) const
        { return link1.numLink > link2.numLink; }
    };

    struct GraphPath
    {
		long selfID;
		std::vector<long> nodeID;
		unsigned sumLink;

        GraphPath(): selfID(), nodeID(), sumLink(0) {}
        GraphPath(long ID, unsigned long size): selfID(ID), nodeID(size) {}
        GraphPath(long ID, unsigned long size, unsigned sum): selfID(ID), nodeID(size), sumLink(sum) {}
	};	

    struct GraphPathSelfIDLess
    {
        bool operator() (const GraphPath& path1, const GraphPath& path2) const
        { return path1.selfID < path2.selfID; }
    };

    struct NodeInfoForPathSearch
    {
		unsigned numVisit;
		long distance;
		long preNodeID;

        NodeInfoForPathSearch(): numVisit(0), distance(0) {}
        NodeInfoForPathSearch(unsigned n, long d, long p):  numVisit(n), distance(d), preNodeID(p) {}
	};	

    struct MappedTagInfo 
    {
		int tagID;
		platanus::Position position;

        MappedTagInfo(): tagID(0), position() {}
        MappedTagInfo(const int id, const int offset, const int seqID): tagID(id), position(offset, seqID) {}
        ~MappedTagInfo() = default;

        bool operator<(const MappedTagInfo &a) const
        {
            if (tagID != a.tagID)
                return (tagID < a.tagID);
			else
                return (position < a.position);
        }
    };


	static const unsigned MIN_LENGTH_FOR_MGC;
	static const unsigned MIN_LENGTH_FOR_COV_ESTIMATE;
	static const double EMP_DISTR_BIN_SIZE;

	// for node-state
    static const unsigned DBG_HETERO;
	static const unsigned DBG_PRIMARY_BUBBLE;
    static const unsigned DBG_SECONDARY_BUBBLE;

	// for edge-state
    static const char DBG_OVERLAP;

	// for contig-state
    static const char DBG_CONTIG_BUBBLE_JUNCTION;
	static const char DBG_CONTIG_PRIMARY_BUBBLE;
	static const char DBG_CONTIG_SECONDARY_BUBBLE;


	void storeGraphLinkFromOverlap(std::vector<GraphLinkWithFlag> &graphLinkPool);
	void storeGraphLinkFromMappedPair(std::vector<GraphLinkWithFlag> &graphLinkPool, long numThread);
	void storeGraphLinkFromMappedLongRead(std::vector<GraphLinkWithFlag> &graphLinkPool, long numThread);
    virtual void calcLink(const long libraryIndex, const long linkThreshold, const long numThread);
    virtual double calcNodeCoverage(const GraphNode &node);
	void calcLinkAndWriteGraphLinkWithFlagFile(const std::vector<GraphLinkWithFlag>& graphLinkPool, const GraphLinkWithFlagPoolIndex& index, const long libraryIndex, const bool multiThreadFlag);
	void getOverlappedBubbleNodeIndex(std::vector<long> &bubbleNodeIndex);
	void getOverlappedBubbleNodePairID(std::vector<std::pair<long, long> > &bubbleNodePairID);
	void getGappedBubbleNodePairID(std::vector<std::pair<long, long> > &bubbleNodePairID);
	void getGappedConflictingNodePairID(std::vector<std::pair<long, long> > &bubbleNodePairID, double heteroNodeCoverage);
	void markBubbleHeteroNode(const std::vector<long> &candidateNodeIndex, const double maxHeteroCoverageFactor);
	void calculateHeteroCoverage(const std::vector<long> &bubbleNodeIndex);
	long writeAndMarkOverlappedNodes(const std::vector<long> &nodeID, FILE *storeFP);
	long remakeGraphAccordingToPath(std::vector<GraphPath> &pathBufferForEachThread);
	long remakeGraphAccordingToPathPair(std::vector<GraphPath> &pathBufferForEachThread);
	void setOppositeBubbleNodeID(std::vector<long> &nodeIDVector, const std::vector<ScaffoldPart> &partVector);
	void flipOppositeBubbleNodeID(std::vector<long> &nodeIDVector, const std::vector<ScaffoldPart> &partVector);
	void setOppositeBubbleNodeIDStrandAware(std::vector<long> &nodeIDVector, const std::vector<ScaffoldPart> &partVector);
	long maxLengthContigID(std::vector<long> &IDs, const long start, const long end);
	std::pair<long, long> fillMajorityIDRun(std::vector<long> &IDs, const GraphNode &targetNode, const std::pair<long, long> &ends, const double scoreFactor);
	long getNonGapContigLengthOfNode(const GraphNode &targetNode);
	long getNumEdgeDirectionOfNode(const GraphNode &targetNode);
	long deleteDifferentBubbleEdge(const long numThread);
	double calcNodeCoveragePartial(const GraphNode &node, const long start, const long end);
	void setCorrespondingNodePosition(std::vector<platanus::Position> &positionVector, const std::vector<ScaffoldPart> &partVector);
	void smoothNodeIDVector(std::vector<long> &nodeIDVector, const GraphNode &targetNode, const double scoreFactor);
	long deleteErroneousEdgeNumTagRate(const bool overlapFlag, const long numThread);
    void divideErroneousLink(const std::vector<std::vector<unsigned> >& numErroneousPair, const std::vector<std::vector<unsigned> >& numSpanningPair, const std::vector<std::vector<double> >& sumExpectedLink, std::unordered_map<std::pair<int, int>, bool, platanus::PairHash, platanus::PairEqual> &errorLink, const long minLink, const DIVISION_MODE divisionMode, const long numThread, const long maxGapSize=0);
    void countPairsSpanningGap(std::vector<std::vector<unsigned> >& numSpanningPair, const DIVISION_MODE divisionMode, const long numThread);
	void countLongReadSpanningGap(std::vector<std::vector<unsigned> >& numSpanningPair, const DIVISION_MODE divisionMode, long numThread);
    void countLinksInsideContigs(std::vector<std::vector<unsigned> >& numPair, const long numThread);
	void countLongReadLinksInsideContigs(std::vector<std::vector<unsigned> >& numPair, const long numThread);
	void countSwitchErrorLinks(std::vector<std::vector<unsigned> >& numPair, const long numThread);
	void countLongReadSwitchErrorLinks(std::vector<std::vector<unsigned> >& numPair, const long numThread);
	void markJunctionContigJoinedToBubble();
	long getScoreFromIDPair(long leftNodeID, long rightNodeID);
	void storeGraphLinkFromTagReadPair(std::vector<GraphLinkWithFlag> &graphLinkPool, const bool scaffoldFlag, long numThread);


    inline bool judgeCloseness(const double val1, const double val2, const double rateThreshold) const
    {
        return (val1 * rateThreshold >= val2) && (val2 * rateThreshold >= val1);
    }


public:

    PairedDBG();
    void setContigMaxK(const long len) { contigMaxK = len; }
    void setMinTagIslandLengthToScaffold(const long num) { minTagIslandLengthToScaffold = num; }
	void setTolerenceFactor(const double factor) { tolerenceFactor = factor; }
	void setCutoffLengthFactor(const double factor) { cutoffLengthFactor = factor; }
	void setMode(const unsigned bits) { mode = bits; }
	void unsetMode(unsigned bits) { mode &= ~bits; }
	void setHeteroCoverage(const double cov) { heteroCoverage = cov; }
	void setNumInputBubbleContig(const unsigned long num) { numInputBubbleContig = num; }
	void setUondProbThreshold(const double prob) { condProbThreshold = prob; }

	unsigned getMode() { return mode; }

    virtual void makeGraph(const long numThread);
	void makeGraphAllLibraries(const long numThread);
	void resetGraph();
	void getOverlappedNode(const long sourceNodeIndex, const char targetDirection, std::vector<long> &nodeIDBuffer);
	long getNumEdgeOneDirection(const GraphNode &targetNode, const char targetDirection);
	void extractDBGBubbleInformation();
	void calculateHeteroAndAverageCoverageUnphase();
	void markHeteroNode(const double maxHeteroCoverageFactor);
	void deleteNonOverlapHomoEdge();
	long joinUnambiguousNodePair(const long numThread);
	long joinUnambiguousNodePath(const long numThread);
	bool isRedundantCrossComponent(const std::array<std::array<long, 2>, 2 > &externalNodeID, const long centerNodeID);
	bool isRedundantGappedCrossComponent(const std::array<std::array<NodeIDWithGap, 2>, 2 > &externalNodeID, const long centerNodeID);
	long solveSimpleCrossStructure(const double linkRateThreshold, const CROSS_RESOLUTION_MODE resolutionMode, const long numThread);
	long solveUniquePathBetweenLinkedNodePair(const long numThread);
	void solveUniquePathBetweenLinkedNodePairAllLibrariesIterative(const long numThread);
	void searchUniqueOverlapPathGuidedByEdge(const long startNodeID, const GraphEdge &guideEdge, std::vector<GraphPath> &pathBuffer);
    virtual void cutAndPrintSeq(long minSeqLength, const unsigned long long readLength, const std::string &outFilename, const std::string &componentFilename);
	long crushSimpleDBGBubble();
    virtual unsigned long long crushHeteroBubble(const double averageCoverage);
    virtual long deleteHeteroEdge(void);
	virtual void loadResultSeq(const long minSeqLength, const unsigned long long readLength, const bool trimFlag, const bool clusterFlag);
	void markRedundantResultSeq(const long numThread, const bool clusterFlag=false);
	void solveSimpleCrossStructureIterative(const CROSS_RESOLUTION_MODE resolutionMode, const bool clusterFlag, const long numThread);
	void solveSimpleCrossStructureAllLibrariesIterative(const CROSS_RESOLUTION_MODE resolutionMode, const bool clusterFlag, const long numThread);
	void joinUnambiguousNodePairIterative(const long numThread);
	void joinUnambiguousNodePairGappedIterativeAllLibraries(const long numThread);
	long solveSimpleGappedCrossStructure(const double linkRateThreshold, const CROSS_RESOLUTION_MODE resolutionMode, const long numThread);
	void solveSimpleGappedCrossStructureIterative(const CROSS_RESOLUTION_MODE resolutionMode, const bool clusterFlag, const long numThread);
	void solveSimpleGappedCrossStructureAllLibrariesIterative(const CROSS_RESOLUTION_MODE resolutionMode, const bool clusterFlag, const long numThread);
	void solveUniquePathBetweenLinkedNodePairIterative(const long numThread);
	void etOppositeContigIDGappedConflicting(const long numThread);
	void setOppositeBubbleContigIDOverlapped(const long numThread);
	void setOppositeBubbleContigIDGapped(const long numThread);
	void setOppositeContigIDGappedConflicting(const long numThread);
	long divideNodeUsingBubbleContigPair(const long numThread);
	long divideNodeUsingBubbleContigPairStrandAware(const long numThread);
	void setOppositeBubbleNodeIDForEachNode(const long numThread);
	void setOppositeBubbleNodeIDAndStateForEachNode();
	void outputResultSeqWithBubble(const std::string filePrefix, const std::string &primaryBubbleSuffix, const std::string &secondaryBubbleSuffix, const std::string &nonBubbleHeteroSuffix, const std::string &nonBubbleOtherSuffix, const std::string &pairSuffix, const long minSeqLength, const long readLength);
	void outputResultSeqComponent(const std::string filePrefix, const std::string &fileSuffix);
	void deleteDifferentBubbleEdgeIterative(const long numThread);
	void deleteThinEdgeCostantKeepingOverlap(const long linkThreshold);
	long deleteConflictingBubbleEdge(const long numThread);
	long deleteSecondaryBubbleNodeAndEdge(const long numThread);
	long deleteShortAndLowCoverageBranch(const long lengthThreshold, const double coverageThreshold, const long numThread);
	void deleteShortAndLowCoverageBranchIterative(const long lengthThreshold, const double coverageThreshold, const long numThread);
	void setBubbleJunctionContigIDOverlapped();
	long divideBubbleJunctionNode(const bool gapDivideFlag);
	long divideBubbleContigInNonHeteroNode();
	void divideGappedNode(const long minGapSize);
	void copyAllNodes(PairedDBG &targetGraph);
	void divideUsingGuideGraph(PairedDBG &guideGraph, const bool bubbleDivideFlag);
	void deleteEdgesUsingGuideGraph(PairedDBG &guideGraph);
	void extendUsingGuideGraph(PairedDBG &guideGraph);
	long getQuartileLengthOfBubble(const unsigned long quartileNumber);
    virtual void detectRepeat(const double averageCoverage);
    virtual void deleteRepeatEdge(void);
	void getUniqueConflictingNode(const long sourceNodeIndex, const char targetDirection, std::vector<NodeIDWithGap> &nodeIDBuffer);
	void divideNestedBubbleNode(const long numThread);
	void deleteLongEdge(const long maxEdgeLength);
	void deleteErroneousEdgeNumTagRateIterative(const bool overlapFlag, const long numThread);
	void deleteEdgeFromShortNodeKeepingBubble(const long lengthThreshold);
	void trimSparseEnd();
	void trimRepeatEnd();
	long divideInconsistentBubbleEnd();
	void adjustOppositeBubbleNodeIDDirection();
	void saveBubble();
	void deleteEdgeFromSecondaryBubble();
	void divideNodeBasedOnBubblesIterative(const bool strandFlag, const long numThread);
	void divideNodeByBubbleBoundary(const long numThread);
	void remakeGraphRecoveringSecondaryBubble(PairedDBG &bubbleGraph);
	void divideNodeByPrimaryBubbleBoundary(const PairedDBG &bubbleGraph);
	void divideErroneousNode(const long minLink, const DIVISION_MODE divisionMode, const long numThread, const long maxGapSize=0);
    virtual void makeScaffold(void);
    void makeScaffoldCombine(void);
	void deleteEdgeFromDifferentPreviousParent(const long numThread);
	void clearContigPreviousParentNodeID();
	void setOppositeBubbleContigIDByEndMatch();
    virtual void insertSizeDistribution(std::vector<SeqLib>& library, std::vector<long>& distribution, const long numThread);
	virtual void updateInsertLengthFP(std::vector<SeqLib>& lib, const long numThread);
	long solveSimpleCrossStructureCoverage(const double minCoverage, const long minLength, const double sameCoverageThreshold, const double diffCoverageThreshold, const long numThread);
	void solveSimpleCrossStructureCoverageIterative(const bool clusterFlag, const long numThread);
	void deleteChimericEdge();
	void setMgcEmpDistr();
	void setCoverageEmpDistr(const std::string &kmerOccurrenceFile);
	bool setClusterID(const platanus::Contig &contig, const std::string ClusterTsvFilename);
	void setClusterIDfromSeqName(const platanus::Contig &contig);
	void deleteInterOTUEdge();
	void splitResultSeqForClusters(std::vector<std::vector<std::string > > &seq,  std::vector<std::vector<std::string> > &name);
	void estimateTaggedMoleculeLength(const long numThread);
	virtual void dumpAllEdges(const std::string &outputFilename);
	long deleteConflictingEdgeToSameNode(const long numThread);

    long getTolerence(void) const { return tolerence; }
    long getNumNode(void) const { return numNode; }
	double getHeteroCoverage() { return heteroCoverage; } 
    unsigned long long getSeqTotalLength(void) const
    {
        unsigned long long total = 0;
        for (auto itr = contig.begin(), end = contig.end(); itr != end; ++itr) {
            total += itr->length;
        }
        return total;
    }



	static const unsigned OVERLAP_MODE;
    static const unsigned PAIRED_END_LINK_MODE; 
    static const unsigned LENGTH_CUTOFF_MODE;
    static const unsigned NON_DBG_OVERLAP_MODE;
	static const unsigned LONG_READ_LINK_MODE;
	static const unsigned BUBBLE_AWARE_MODE;
	static const unsigned SECONDARY_BUBBLE_REMOVAL_MODE;
	static const unsigned PREVIOUS_DIVISION_AWARE_MODE;
	static const unsigned CLUSTERING_MODE;
	static const unsigned INTRA_OTU_MODE;
	static const unsigned TAG_SCAFFOLD_MODE;

	static const long MAX_ITERATION_OF_CROSS_SOLUTION;

    static const double NON_PAIRED_END_TOLERENCE_FACTOR;
    static const double HETERO_COVERAGE_THRESHOLD_FACTOR;
	static const double CROSS_LINK_RATE_THRESHOLD;
	static const double CROSS_SCORE_RATE_THRESHOLD;
	static const double MIN_BUBBLE_COUNT_FACTOR;
};


#endif
