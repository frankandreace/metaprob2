#ifndef	__HASHTABLE_H__
#define __HASHTABLE_H__ 

#include "Utilities.h"
#include "BloomFilter/RepeatBloomFilter.h"

typedef vector<id_seq_type> KmNode; //Kmer Node used in HashTable

class clsHashTable
{
private:
	size_seq Q_SIZE;

public:
	typedef shared_ptr<clsHashTable> Ptr;
	typedef unordered_map<hash_type, KmNode> MapHash;

	MapHash arrMyHash;

	bool InitOptionData(size_seq Q_SIZE, hash_type expected_element);//Initialize some optional variables and the hash table
	void ExtractFromString(id_seq_type iSeqInd, string& s_Str, RepeatBloomFilter::Ptr repeatQmer);//Extract l-mer from a string which belongs to sequence having id iSeqInd

	size_seq getQSize() const;
};

void GetMapSeqID(MapAdjSeqID& map, KmNode& kmnode); //add in Seq in vector

#endif
