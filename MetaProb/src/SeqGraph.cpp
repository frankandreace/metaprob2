#include"SeqGraph.h"

clsSeqGraph::clsSeqGraph(clsDNAFile::Ptr mFile) {
	this->mFile = mFile;
	this->qMerFilter = 0;
}

bool clsSeqGraph::FilteredQmer(MapAdjSeqID& map)
{
	id_seq_type num_different_read = map.size();
	if(num_different_read <= 1 )//|| num_different_read >= 10) //Filtro solo da 1, quindi Q-mer che individuano un singolo read
	{
		++this->qMerFilter;
		return true;
	}
	else
		return false;
}

//OLD, memoria occupata eccessiva
//Create graph from hashtable of q-mers, eliminate q-mer Ptr (remain q-mer key)
//void clsSeqGraph::CreateGraph()
//{
//	cout << endl << "Building Graph... " << flush;
//	auto key_kmer = this->mFile->getHash()->arrMyHash.begin();
//	while(key_kmer != this->mFile->getHash()->arrMyHash.end())
//	{
////		cout << "\r" << "Building Graph... Analyzing q-mer: " << key + 1 << flush;
//		MapAdjSeqID vcrIDSeq;
//		GetMapSeqID(vcrIDSeq, key_kmer->second);
//		key_kmer = this->mFile->getHash()->arrMyHash.erase(key_kmer);
//		if(!this->FilteredQmer(vcrIDSeq))
//		{
//			auto vcrIDSeq_it_end = vcrIDSeq.end();
//			auto vcrIDSeq_it_curr = vcrIDSeq.begin();
//			while(vcrIDSeq_it_curr != vcrIDSeq_it_end)
//			{
//				auto vcrIDSeq_it_succ = next(vcrIDSeq_it_curr);
//				while(vcrIDSeq_it_succ != vcrIDSeq_it_end)
//				{
//					//Evita archi su se stessi con poco senso
//					if(vcrIDSeq_it_curr->first != vcrIDSeq_it_succ->first)
//					{
//						size_seq_tot count = (size_seq_tot)vcrIDSeq_it_curr->second * (size_seq_tot)vcrIDSeq_it_succ->second;
//						this->mFile->GetAt(vcrIDSeq_it_curr->first)->AddAdjVertice(this->mFile->GetAt(vcrIDSeq_it_succ->first), count);
//						this->mFile->GetAt(vcrIDSeq_it_succ->first)->AddAdjVertice(this->mFile->GetAt(vcrIDSeq_it_curr->first), count);
//					}
//					++vcrIDSeq_it_succ;
//				}
//				++vcrIDSeq_it_curr;
//			}
//		}
//	}
//
////	//Elimino adj che portano a read paired-end
////	//Control if is other End of paired End, otherwise jump. Paired-end not overlap
////	auto seqPE_it = this->mFile->begin();
////	while(seqPE_it != this->mFile->end())
////	{
////		id_seq_type ID_other_end = this->mFile->GetIDotherEnd((*seqPE_it)->getIdFile(), (*seqPE_it)->getSeqId());
////		(*seqPE_it)->vAdjSeqPE->erase(ID_other_end);
////		++seqPE_it;
////	}
//
//	cout << "Complete" << flush;
//}

//Create graph from hashtable of q-mers
void clsSeqGraph::CreateGraph()
{
	this->mFile->CreateGraph();
}

void clsSeqGraph::ResetAdj() {
	auto ID_it = this->mFile->begin();
	while(ID_it != this->mFile->end())
	{
		(*ID_it)->reset(clsSequencePE::Graph_reset);
		++ID_it;
	}
}

hash_type clsSeqGraph::getQMerFilter() const {
	return qMerFilter;
}

const clsDNAFile::Ptr& clsSeqGraph::getFile() const {
	return mFile;
}
