#include "Cluster.h"

void clsCluster::InitID(id_cluster_type ID)
{
	this->ID=ID;
	this->ClearMean();
}

void clsCluster::InitMean(lMer_type feature_size) {
	this->vMean = move(vector<double>(feature_size, 0));
}

void clsCluster::ClearMean(){
	this->vMean.clear();
}

//Addition a group's pointer to this cluster
void clsCluster::AddItem(clsGroup::Ptr c_Grp)
{
	this->vGroup.push_back(c_Grp);
}

//Compute the distance between a group g and this cluster's mean
double clsCluster::ComputeDistance(clsGroup::Ptr g)
{
	vector<double> dist;
	for(lMer_type i = 0; i < this->vMean.size(); i++)
		dist.push_back(abs(this->vMean[i] - g->vFreq[i]));

	//Passaggio per evitare overflow
	double y = 0;
	for(lMer_type i = 0; i < this->vMean.size(); i++)
		if(y < dist[i])
			y = dist[i];

	//Calcolo con divisione per y
	double kq=0;
	for(lMer_type i=0;i<this->vMean.size();i++)
		kq += pow(dist[i]/y,2);

	//Rimoltiplico y
	kq = y * sqrt(kq);

	return kq;
}

//ReCompute mean of this cluster, return false if the mean is not changed
bool clsCluster::ReComputeMean()
{
	bool change = false;//not change
	id_grp_type iGSize = this->vGroup.size();

//	//Threaded recomputed mean
//	vector<future<double>> fut;
//	for(lMer_type j = 0; j < this->vMean.size(); j++)
//	{
//		fut.push_back(async(launch::async, [iGSize,j,this]
//													   {
//														double dVal=0.0;
//														for(id_grp_type i = 0; i < iGSize; i++)
//															dVal += this->vGroup[i]->vFreq[j];
//														dVal = dVal/(double)iGSize;
//														return dVal;
//													   }));
//	}
//
//	for(lMer_type j = 0; j < this->vMean.size(); j++)
//	{
//		double dVal=fut[j].get();
//		if(this->vMean[j] != dVal)
//		{
//			this->vMean[j] = dVal;
//			change = true;//change
//		}
//	}

	//Not threaded recomputed mean
	for(lMer_type j = 0; j < this->vMean.size(); j++)
	{
		double dVal=0.0;
		for(id_grp_type i = 0; i < iGSize; i++)
			dVal += this->vGroup[i]->vFreq[j];
		dVal = dVal/(double)iGSize;
		if(this->vMean[j] != dVal)
		{
			this->vMean[j] = dVal;
			change = true;//change
		}
	}

//	//Not threaded Density recomputed mean
//	id_seq_type tot = 0;
//	for(id_grp_type i = 0; i < this->vGroup.size(); i++)
//		tot += this->vGroup[i]->GetSize();
//
//	bool flag = false;//not change
//	for(lMer_type j = 0;j < this->vMean.size(); j++)
//	{
//		double dVal=0.0;
//		for(id_grp_type i = 0; i < this->vGroup.size(); i++)
//			dVal+= this->vGroup[i]->vFreq[j] *  this->vGroup[i]->GetSize();
//		dVal=dVal/(double)tot;
//		if(this->vMean[j] != dVal)
//		{
//			this->vMean[j]=dVal;
//			flag=true;//change
//		}
//	}

	return change;
}

//Remove all items
void clsCluster::RemoveAllItem()
{
	this->vGroup.clear();
}

void clsCluster::AssignSequenceToCluster() {
	for(id_grp_type iGr = 0; iGr < this->vGroup.size();iGr++)
		this->vGroup[iGr]->SetIDCluster(this->ID);
}

id_cluster_type clsCluster::GetID() {
	return this->ID;
}

vector<double>& clsCluster::GetMean() {
	return this->vMean;
}

id_seq_type clsCluster::GetSize() {
	id_seq_type size = 0;
	for(id_grp_type iGr = 0; iGr < this->vGroup.size();iGr++)
		size += this->vGroup[iGr]->GetSize();
	return size;
}

bool clsCluster::isEmpty()
{
	return this->GetSize() == 0;
}

const vector<clsGroup::Ptr>& clsCluster::getGroup() const {
	return vGroup;
}
