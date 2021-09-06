#ifndef	__GROUP_H__
#define __GROUP_H__ 

#include "Utilities.h"
#include "SequencePE.h"

class clsGroup
{
private:
	id_grp_type ID;//ID of the group

	id_seq_type idxSeGrp;//Index order in seed
	id_seq_type idxGrp;//Index order in list of reads which are not in seed

	id_cluster_type IDClus; //ID of Cluster assigned

	bool choose;

public:
	typedef shared_ptr<clsGroup> Ptr;
	typedef vector<long double> VectorFeature;
	typedef unordered_map<id_grp_type, id_grp_type> MapAdjGrp; //<gruppo adiacente, quante volte lo è stato>

	vector<clsSequencePE::Ptr>vSeedSeqPE;//List of clsVirSeq's pointer belong to this group
	vector<clsSequencePE::Ptr>vSeqPE;//List of clsVirSeq's pointer belong to this group

	MapAdjGrp adj;
	
	VectorFeature vCountMer;//Count Mer vector
	VectorFeature vFreq;//Frequency vector

public:
	void InitData(id_grp_type ID);
	void AddItem(clsSequencePE::Ptr c_SeqPE, SeedState seed_state);//Addition a read's pointer to this group, to seed (2) or not (1)
	void AddItemG(clsGroup::Ptr c_Grp, size_seq i_thres);//Add a group of reads to this group
	void AddAdjacenceGrp(id_grp_type ID_grp); //Add Adjacent group

	SeedState IsCandidateR(clsSequencePE::Ptr c_SeqPE, size_seq i_thres); //check whether a read can be add to this group (in seed (2), or not in seed (1)) or not (0), number of shared l-mer i_thres
	bool IsCandidateG(clsGroup::Ptr c_Grp, size_seq i_thres);//check whether a group can be add to this group or not
	SeedState SelectNextR(id_seq_type& iCRead);//Select the next read which will be used to find other reads
	
	size_seq_tot get_Seeds_L_mer_count(size_seq l_mer_freq);
	size_seq_tot get_Seeds_Size();
	size_seq_tot get_NoSeeds_Size();

	id_seq_type GetSize() const;
	bool isEmpty() const;

	//Getter and Setter
	void SetID(id_grp_type  ID);
	id_grp_type GetID() const;

	void SetIDCluster(id_cluster_type  ID_cluster);
	id_grp_type GetIDCluster() const;

	bool isChoose() const;
	void setChoose(bool choose);
};
#endif
