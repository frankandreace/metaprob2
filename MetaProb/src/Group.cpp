
#include "Group.h"
void clsGroup::InitData(id_grp_type ID)
{
	this->SetID(ID);
	this->idxGrp=0;
	this->idxSeGrp=0;

	this->choose = false;
	this->SetIDCluster(numeric_limits<id_cluster_type>::max()); //valore senza senso (riconoscibile)
}

//Addition a read's pointer to this cluster, to seed (2) or not (1)
void clsGroup::AddItem(clsSequencePE::Ptr c_SeqPE, SeedState seed_state)
{
	if(seed_state == Seed)
	{
		c_SeqPE->setIdGrp(this->ID);
		c_SeqPE->setSeed(Seed);
		this->vSeedSeqPE.push_back(c_SeqPE);
	}
	else if (seed_state == No_Seed)
	{
		c_SeqPE->setIdGrp(this->ID);
		c_SeqPE->setSeed(No_Seed);
		this->vSeqPE.push_back(c_SeqPE);
	}
}

//Add a group of reads to this group
void clsGroup::AddItemG(clsGroup::Ptr c_Grp, size_seq i_thres)
{
	//check whether each read will be assign into this group's seed or not
	for(id_seq_type i=0;i<c_Grp->vSeedSeqPE.size();i++)
	{
		bool flag=false;
		for(id_seq_type j=0;j<this->vSeedSeqPE.size();j++)
			if(c_Grp->vSeedSeqPE[i]->IsAdjacent(this->vSeedSeqPE[j], i_thres))
			{
				flag=true;
				break;
			}
		if(flag == false)//read in the group's seed
			this->AddItem(c_Grp->vSeedSeqPE[i], Seed);
		else //read not in the group's seed
			this->AddItem(c_Grp->vSeedSeqPE[i], No_Seed);
	}

	for(id_seq_type i=0;i < c_Grp->vSeqPE.size();i++)
	{
		bool flag=false;
		for(id_seq_type j=0;j<this->vSeedSeqPE.size();j++)
			if(c_Grp->vSeqPE[i]->IsAdjacent(this->vSeedSeqPE[j], i_thres))
			{
				flag=true;
				break;
			}
		if(flag == false)//read in the group's seed
			this->AddItem(c_Grp->vSeqPE[i], Seed);
		else //read not in the group's seed
			this->AddItem(c_Grp->vSeqPE[i], No_Seed);
	}
}

void clsGroup::AddAdjacenceGrp(id_grp_type ID_grp)
{
	auto ID_it = this->adj.find(ID_grp);
	if(ID_it == this->adj.end()) //Non trovato
		this->adj.insert(MapAdjGrp::value_type(ID_grp, 1));
	else
		ID_it->second++;
}

//check whether a read can be add to this group (in seed (2), or not in seed (1)) or not (0)
SeedState clsGroup::IsCandidateR(clsSequencePE::Ptr c_SeqPE, size_seq i_thres)
{
	//Prima controllo se è adiacente a nodi Seed (se si non può essere seed)
	for(id_seq_type i=0;i<this->vSeedSeqPE.size();i++)
		if(this->vSeedSeqPE[i]->IsAdjacent(c_SeqPE, i_thres))
			return No_Seed;//not in seed

	//Controllo ora adiacenza a nodi non seed e se è adiacente posso inserirlo nei nodi seed
	//perchè il controllo fatto prima esclude sia adiacente a nodi Seed
	for(id_seq_type i=0;i<this->vSeqPE.size();i++)
		if(this->vSeqPE[i]->IsAdjacent(c_SeqPE, i_thres))//We still do not concern the number of adjacent reads
			return Seed;//in seed

	//Non è adiacente a sequenze nel gruppo
	return No_Other_Read;//can not be in this group
}

//check whether a group can be add to this group or not
bool clsGroup::IsCandidateG(clsGroup::Ptr c_Grp, size_seq i_thres)
{
	for(id_seq_type i=0;i<c_Grp->vSeedSeqPE.size();i++)
		if(this->IsCandidateR(c_Grp->vSeedSeqPE[i], i_thres) != 0)//We still only check that there is ONLY ONE a common read between two groups
			return true;

	for(id_seq_type i=0;i<c_Grp->vSeqPE.size();i++)
		if(this->IsCandidateR(c_Grp->vSeqPE[i], i_thres) != 0)//We still only check that there is ONLY ONE a common read between two groups
			return true;
	return false;
}

//Select the next read which will be used to find other reads
SeedState clsGroup::SelectNextR(id_seq_type& iCRead)
{
	if(this->idxGrp < this->vSeqPE.size())
	{
		iCRead=this->vSeqPE[idxGrp]->getSeqId();
		++this->idxGrp;
		return No_Seed;//not in seed
	}
	if(this->idxSeGrp < this->vSeedSeqPE.size())
	{
		iCRead=this->vSeedSeqPE[idxSeGrp]->getSeqId();
		++this->idxSeGrp;
		return Seed;//in seed
	}
	return No_Other_Read;//there is no any reads
}

size_seq_tot clsGroup::get_Seeds_L_mer_count(size_seq l_mer_freq)
{
	size_seq_tot L_mer_count = 0;
	for(id_seq_type i=0;i<this->vSeedSeqPE.size();i++)
		L_mer_count += (size_seq_tot)this->vSeedSeqPE[i]->get_L_mer_count(l_mer_freq);
	return L_mer_count;
}

size_seq_tot clsGroup::get_Seeds_Size()
{
	size_seq_tot seed_size = 0;
	for(id_seq_type i=0;i<this->vSeedSeqPE.size();i++)
		seed_size += (size_seq_tot)this->vSeedSeqPE[i]->get_size();
	return seed_size;
}

size_seq_tot clsGroup::get_NoSeeds_Size()
{
	size_seq_tot noSeed_size = 0;
	for(id_seq_type i=0;i<this->vSeqPE.size();i++)
		noSeed_size += (size_seq_tot)this->vSeqPE[i]->get_size();
	return noSeed_size;
}

void clsGroup::SetIDCluster(id_cluster_type ID_cluster)
{
	this->IDClus = ID_cluster;
	//Assigns the no seed in a group to a cluster
	for(id_seq_type iRed = 0; iRed < this->vSeqPE.size(); iRed++)
		this->vSeqPE[iRed]->setIdCluster(ID_cluster);

	//Assigns the seed in a group to a cluster
	for(id_seq_type iRed = 0; iRed < this->vSeedSeqPE.size(); iRed++)
		this->vSeedSeqPE[iRed]->setIdCluster(ID_cluster);
}

id_grp_type clsGroup::GetIDCluster() const {
	return this->IDClus;
}

void clsGroup::SetID(id_grp_type ID)
{
	this->ID = ID;
	for(id_grp_type i = 0; i < this->vSeqPE.size(); i++)
		this->vSeqPE[i]->setIdGrp(ID);
	for(id_grp_type i = 0; i < this->vSeedSeqPE.size(); i++)
		this->vSeedSeqPE[i]->setIdGrp(ID);
}

id_grp_type clsGroup::GetID() const {
	return this->ID;
}

id_seq_type clsGroup::GetSize() const {
	id_seq_type tot = 0;
	for(id_seq_type i = 0; i < this->vSeqPE.size(); i++)
	{
		if(this->vSeqPE[i]->isPairedEnd())
			tot += 2;
		else
			tot++;
	}
	for(id_seq_type i = 0; i < this->vSeedSeqPE.size(); i++)
	{
		if(this->vSeedSeqPE[i]->isPairedEnd())
			tot += 2;
		else
			tot++;
	}
	return tot;
}

bool clsGroup::isEmpty() const {
	return this->GetSize() == 0;
}

bool clsGroup::isChoose() const {
	return choose;
}

void clsGroup::setChoose(bool choose) {
	this->choose = choose;
}
