#ifndef	__SEQUENCEPE_H__
#define __SEQUENCEPE_H__ 

#include "HandlerInput/Input.h"
#include "Utilities.h"
#include "VirSeq.h"
#include "CountNucleotide.h"

//Is also the Node of the SeqGraph implicitly
class clsSequencePE:public clsVirSeq
{
private:
//	typedef pair<string, bool> Sequence_isAdd;
	typedef bool Sequence_isAdd;
	typedef pair<Sequence_isAdd,Sequence_isAdd> SequencePresence;

	//Variable for sequence details
	PairFiles& File; //File
	SeedState seed; //this read is a seed (TRUE) or not (FALSE) of a group
	id_grp_type idGrp; //ID grp assign
	id_cluster_type idCluster; //ID cluster assign
	CountNucleotide<size_seq> info; //Nucleotide num

	SequencePresence seqPE; //Contains if sequences is present in file
public:
	typedef shared_ptr<clsSequencePE> Ptr;
	typedef unordered_map<id_seq_type, size_seq_tot> MapAdj;
	typedef unique_ptr<MapAdj> MapAdj_Ptr;

	//Variable for graph details
	MapAdj_Ptr vAdjSeqPE;//List of adjacent vertices' pointer

	clsSequencePE(PairFiles& File, id_seq_type iSeqID);//Initialize the sequence
	virtual ~clsSequencePE();
	void AddSequence(string& seq, int flagEnd);
	bool isPairedEnd();

	//Method to SeqGraph
	virtual void AddAdjVertice(clsSequencePE::Ptr c_SeqPE);//Add a adjacent vertice (read) to this read
	virtual void AddAdjVertice(clsSequencePE::Ptr c_SeqPE, size_seq_tot count);//Add a adjacent vertice (read) to this read
	virtual bool IsAdjacent(clsSequencePE::Ptr c_SeqPE, size_seq i_thres);//Check whether c_SeqPE is adjacent (TRUE) or not (FALSE) with this read

	//Method to Norm calc
	size_seq get_L_mer_count(size_seq l_mer_freq);
	size_seq get_size();

	enum ResetType{Seed_Reset, Group_reset, Cluster_reset, Info_reset, Graph_reset};
	void reset(ResetType type);

	//Getter and Setter variable
	id_cluster_type getIdCluster() const;
	void setIdCluster(id_cluster_type idCluster);

	id_grp_type getIdGrp() const;
	void setIdGrp(id_grp_type idGrp);

	const CountNucleotide<size_seq>& getInfo() const;
	const PairFiles& getFile() const;

	SeedState getSeed() const;
	void setSeed(SeedState seed);

private:
	void Count_Nucleotide(string &seq);
};
#endif
