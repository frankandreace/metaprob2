#ifndef SRC_USATO_GENOMICBLOOMFILTER_H_
#define SRC_USATO_GENOMICBLOOMFILTER_H_

#include <boost/dynamic_bitset.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <string>
#include <vector>
#include <memory>

using namespace std;

//template is the type used for hash value
class GenomicBloomFilter {
public:
	typedef shared_ptr<GenomicBloomFilter> Ptr;

	typedef uint64_t TYPEHASH;
	typedef boost::dynamic_bitset<uint64_t, allocator<uint64_t>> DynamicBitSet;

	GenomicBloomFilter(TYPEHASH num_element, double prob_false_positive, size_t Q);
	virtual ~GenomicBloomFilter();

	virtual void contains(string& read, vector<bool>& posCont);
	virtual void insert(string& read);
protected:
	typedef TYPEHASH (GenomicBloomFilter::*HashFunction)(char&); //puntatori a funzione
	DynamicBitSet bitset;
	size_t Q;
	vector<HashFunction> function; //Vettore di puntatori a funzione

	virtual bool contains(vector<size_t>& index, DynamicBitSet& bit);
	virtual void insert(vector<size_t>& index, DynamicBitSet& bit);
	virtual void GetPositionWithHash(vector<TYPEHASH>& sevenHash, vector<size_t>& pos, DynamicBitSet& bit);
	virtual void GetHashes(string& read, vector<vector<TYPEHASH>>& sevenHashes);
	virtual void sevenHash(string& read, size_t begin, vector<TYPEHASH>& sevenHash);
//	virtual TYPEHASH f_hash(char& ch, int f_type);

	virtual TYPEHASH f_hash1(char& ch);
	virtual TYPEHASH f_hash2(char& ch);
	virtual TYPEHASH f_hash3(char& ch);
	virtual TYPEHASH f_hash4(char& ch);
	virtual TYPEHASH f_hash5(char& ch);
	virtual TYPEHASH f_hash6(char& ch);
	virtual TYPEHASH f_hash7(char& ch);
};

#endif /* SRC_USATO_GENOMICBLOOMFILTER_H_ */
