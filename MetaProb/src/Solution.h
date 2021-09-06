#ifndef	__SOLUTION_H__
#define __SOLUTION_H__ 

#include "Utilities.h"
#include "SequencePE.h"
#include "Group.h"
#include "HandlerInput/DNAFile.h"
#include "HashTable.h"
#include "SeqGraph.h"
#include "Cluster.h"
#include "Clustering/Kmeans.h"

#include "Clustering/ClusterUtility.h"
#include "Clustering/Gmeans.h"
#include "HandlerInput/Parameter.h"
#include "Normalization.h"

class clsSolut
{
public:
	typedef shared_ptr<clsSolut> Ptr;

	//input variables
	Parameter::Ptr param;
	clsDNAFile::Ptr mFile;
	lMer_type feature_size; //number of feature

	//temporary variables
	clsSeqGraph::Ptr cSeqGrph;
	vector<clsGroup::Ptr> vGroup;//List of group of sequence
	vector<clsCluster::Ptr> vClusts;//List of clusters

	vector<string> info;
	
public:
	bool InitData(Parameter::Ptr param);
	bool DoBinning(Normalization::Norm norm_type);//start to bin

	void PrintResult();
	void SaveInfo(string elapsedTimeTot, const string& dir_output);
	void SaveBinning(const string& dir_output);

	enum ResetType{Group_and_Cluster, Cluster};
	void resetComputation(ResetType reset_type);

private:
	bool FirstPhase(); //Grouping read
	bool Grouping();

	bool IsThereAdjBetweenGroup(); //there is adjacent groups
	void UnionGrpIfPossible(); //Union the groups if possible (use paired end information)
	void SetStateSeed();

	bool SecondPhase(Normalization::Norm norm_type);//The Second phase, binning of contigs
	void MergeGroup();//Merge groups
	void AssignsSequenceToCluster(); //Assign the sequence to a cluster looking the group cluster assignement

	void SaveClusters(const string& dir_output, MapIDFile_Header& map_idxFile_header);
	void SaveGroups(const string& dir_output, MapIDFile_Header& map_idxFile_header);
};
#endif
