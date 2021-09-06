#include "SequencePE.h"

clsSequencePE::clsSequencePE(PairFiles& File, id_seq_type iSeqID) : clsVirSeq(iSeqID), File(File)
{
	this->reset(Info_reset);
	this->reset(Seed_Reset);
	this->reset(Group_reset);
	this->reset(Cluster_reset);
	this->reset(Graph_reset);
}

clsSequencePE::~clsSequencePE(){}

void clsSequencePE::AddSequence(string& seq, int flagEnd)
{
	if(flagEnd == 1)
		this->seqPE.first = true;
	else if(flagEnd == 2)
		this->seqPE.second = true;
	if(flagEnd == 1 || flagEnd == 2)
		this->Count_Nucleotide(seq);
}

bool clsSequencePE::isPairedEnd()
{
	return this->seqPE.first && this->seqPE.second;
}

id_cluster_type clsSequencePE::getIdCluster() const {
	return idCluster;
}

void clsSequencePE::setIdCluster(id_cluster_type idCluster) {
	this->idCluster = idCluster;
}

id_grp_type clsSequencePE::getIdGrp() const {
	return idGrp;
}

void clsSequencePE::setIdGrp(id_grp_type idGrp) {
	this->idGrp = idGrp;
}

const CountNucleotide<size_seq>& clsSequencePE::getInfo() const {
	return info;
}

SeedState clsSequencePE::getSeed() const {
	return seed;
}

void clsSequencePE::setSeed(SeedState seed) {
	this->seed = seed;
}

//Count dna nucleotide in a read
void clsSequencePE::Count_Nucleotide(string& seq)
{
	for(size_seq i = 0; i < seq.size(); ++i)
	{
		switch(seq[i])
		{
		case 'A':
			++this->info.a;
			break;
		case 'C':
			++this->info.c;
			break;
		case 'G':
			++this->info.g;
			break;
		case 'T':
			++this->info.t;
			break;
		}
	}
}

//Add a adjacent vertice (read) to this read
inline void clsSequencePE::AddAdjVertice(clsSequencePE::Ptr c_SeqPE)
{
	this->AddAdjVertice(c_SeqPE, 1);
}

//Add a adjacent vertice (read) to this read
inline void clsSequencePE::AddAdjVertice(clsSequencePE::Ptr c_SeqPE, size_seq_tot count)
{
//	#pragma omp critical(updateAdj)
	{
		MapAdj::mapped_type& num = (*this->vAdjSeqPE)[c_SeqPE->getSeqId()];
		num += count;
	}
}

//Check whether c_SeqPE is adjacent with this read
inline bool clsSequencePE::IsAdjacent(clsSequencePE::Ptr c_SeqPE, size_seq i_thres)
{
	auto ID = this->vAdjSeqPE->find(c_SeqPE->getSeqId());
	if(ID == this->vAdjSeqPE->end()) //Non trovato
		return false;
	else
	{
		if(ID->second >= i_thres)//The number of shared l-mers is at least equal to i_thres
			return true;
		else
			return false;
	}
}

size_seq clsSequencePE::get_L_mer_count(size_seq l_mer_freq)
{
	return this->info.get_L_mer_count(this->isPairedEnd(), l_mer_freq);
}

size_seq clsSequencePE::get_size()
{
	return this->info.get_size();
}

void clsSequencePE::reset(ResetType type)
{
	switch(type)
	{
	case(Info_reset):
		this->info = move(CountNucleotide<size_seq>());
		break;
	case(Seed_Reset):
		this->seed = No_State_Seed;
		break;
	case(Group_reset):
		this->idGrp = numeric_limits<id_cluster_type>::max(); //valore senza senso (riconoscibile)
		break;
	case(Cluster_reset):
		this->idCluster = numeric_limits<id_cluster_type>::max(); //valore senza senso (riconoscibile)
		break;
	case(Graph_reset):
		this->vAdjSeqPE.reset(new MapAdj);
		break;
	}
}

const PairFiles& clsSequencePE::getFile() const {
	return File;
}
