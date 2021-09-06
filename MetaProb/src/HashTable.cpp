#include"HashTable.h"

//Initialize some optional variables and the hash table
bool clsHashTable::InitOptionData(size_seq Q_SIZE, hash_type expected_element)
{
	this->Q_SIZE=Q_SIZE;
	this->arrMyHash.reserve(expected_element);
	return true;
}

size_seq clsHashTable::getQSize() const {
	return this->Q_SIZE;
}

//Extract l-mer from a string which belongs to sequence having id iSeqInd
void clsHashTable::ExtractFromString(id_seq_type iSeqInd, string& s_Str,
		RepeatBloomFilter::Ptr repeatQmer) {
	vector<bool> repeat;
	if(repeatQmer)
		repeatQmer->areRepeats(s_Str, repeat);

	vector<HashCorrect> vHash;
	GetHashes(s_Str, this->Q_SIZE, vHash, CharToInt);

	//Insert q-mer in hash and save ID read where is present
	#pragma omp parallel for ordered
	for(int i = 0; i < vHash.size(); ++i)
	{
		if(vHash[i].second)
		{
			#pragma omp ordered
			if(repeatQmer && repeat[i])
			{
				MapHash::mapped_type& node = this->arrMyHash[vHash[i].first];
				node.push_back(iSeqInd);
			}
			else if(!repeatQmer)
			{
				MapHash::mapped_type& node = this->arrMyHash[vHash[i].first];
				node.push_back(iSeqInd);
			}
		}
	}
}

void GetMapSeqID(MapAdjSeqID& map, KmNode& kmnode) {
	map.reserve(kmnode.size());
	auto id_it = kmnode.begin();
	while(id_it != kmnode.end())
	{
		++map[*id_it];
		++id_it;
	}
}
