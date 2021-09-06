#ifndef	__CLUSTER_H__
#define __CLUSTER_H__ 

#include"Utilities.h"
#include"SequencePE.h"
#include"Group.h"

class clsCluster
{
private:
	id_cluster_type ID;//ID of the cluster
	vector<clsGroup::Ptr> vGroup;

	vector<double> vMean;//Frequency vector

public:
	typedef shared_ptr<clsCluster> Ptr;

	void InitID(id_cluster_type ID);
	void InitMean(lMer_type feature_size);
	void ClearMean();
	void AddItem(clsGroup::Ptr c_Grp);//Add a group to this cluster
	double ComputeDistance(clsGroup::Ptr g);//Compute the distance between a group g and the mean
	bool ReComputeMean();//reCompute mean of this cluster, return false if the mean is not changed

	void RemoveAllItem();//Remove all items

	void AssignSequenceToCluster();

	id_cluster_type GetID();
	vector<double>& GetMean();
	id_seq_type GetSize();
	bool isEmpty();
	const vector<clsGroup::Ptr>& getGroup() const;
};
#endif
