/*
 * GenomicBloomFilter.cpp
 *
 *  Created on: 17/nov/2015
 *      Author: samuele
 */

#include "GenomicBloomFilter.h"

GenomicBloomFilter::GenomicBloomFilter(TYPEHASH num_element, double prob_false_positive, size_t Q) {
	double num_hash = 7;
	const TYPEHASH size = -num_hash*num_element/log(1.0 - pow(prob_false_positive, 1.0/num_hash));
	DynamicBitSet(size).swap(this->bitset);
	this->Q = Q;
	this->function.push_back(&GenomicBloomFilter::f_hash1);
	this->function.push_back(&GenomicBloomFilter::f_hash2);
	this->function.push_back(&GenomicBloomFilter::f_hash3);
	this->function.push_back(&GenomicBloomFilter::f_hash4);
	this->function.push_back(&GenomicBloomFilter::f_hash5);
	this->function.push_back(&GenomicBloomFilter::f_hash6);
	this->function.push_back(&GenomicBloomFilter::f_hash7);
}

GenomicBloomFilter::~GenomicBloomFilter() {
}

//vector<bool>& posCont after call contains if the Q-mer starting at position index is present in bloom
inline void GenomicBloomFilter::contains(string& read, vector<bool>& posCont) {
	vector<vector<TYPEHASH>> vSevenHashes;
	this->GetHashes(read, vSevenHashes);
	posCont = vector<bool>(vSevenHashes.size(), false);
	#pragma omp parallel for
	for(size_t pos = 0; pos < vSevenHashes.size(); ++pos)
	{
		vector<size_t> index;
		this->GetPositionWithHash(vSevenHashes[pos], index, this->bitset);
		posCont[pos] = this->contains(index, this->bitset);
	}
}

inline void GenomicBloomFilter::insert(string& read) {
	vector<vector<TYPEHASH>> vSevenHashes;
	this->GetHashes(read, vSevenHashes);
	#pragma omp parallel for
	for(size_t pos = 0; pos < vSevenHashes.size(); ++pos)
	{
		vector<size_t> index;
		this->GetPositionWithHash(vSevenHashes[pos], index, this->bitset);
		this->insert(index, this->bitset);
	}
}

//sevenHash contains the seven hash value calcolated with sevenHash function
//vector pos after call contains the position in bitset rappresented the hash
inline void GenomicBloomFilter::GetPositionWithHash(vector<TYPEHASH>& sevenHash,
		vector<size_t>& pos, DynamicBitSet& bit) {
	pos = vector<size_t>(sevenHash.size(), 0);
	#pragma omp parallel for
	for(size_t i = 0; i < sevenHash.size(); ++i)
		pos[i] = sevenHash[i]%bit.size();
}

//Compute
inline void GenomicBloomFilter::GetHashes(string& read, vector<vector<TYPEHASH>>& sevenHashes)
{
	if(read.size() >= this->Q)
	{
		#pragma omp parallel for ordered
		for(size_t pos=0; pos <= (size_t)(read.size() - this->Q); ++pos)
		{
			#pragma omp ordered
			{
				vector<GenomicBloomFilter::TYPEHASH> vSevenHash;
				if(pos == 0)
					sevenHash(read, pos, vSevenHash);
				else
				{
					vSevenHash = sevenHashes[pos-1];//Copy constructor call
	//				TYPEHASH base = 2;
					#pragma omp parallel for firstprivate(pos)
					for(size_t func = 0; func < this->function.size(); ++func)
					{
						vSevenHash[func] -= (this->*this->function[func])(read[pos-1]);
						vSevenHash[func] >>= 1; //Divido per 2, spostando a dx i bit
						vSevenHash[func] |= (this->*this->function[func])(read[pos + this->Q - 1]) <<(this->Q - 1);//(TYPEHASH)boost::multiprecision::pow((boost::multiprecision::uint1024_t)base, this->Q - 1);
					}
				}
				sevenHashes.push_back(vSevenHash);
			}
		}
	}
}

//Compute seven hash function for one Q-mer starting from begin pos in a read
inline void GenomicBloomFilter::sevenHash(string& read, size_t begin,
		vector<TYPEHASH>& sevenHash)
{
	   for(int func = 0; func < 7; ++func)
		   sevenHash.push_back(0);

//	   TYPEHASH base = 2;
	   #pragma omp parallel for
	   for(size_t pos = begin; pos<begin + this->Q; ++pos)
	   {
		   #pragma omp parallel for firstprivate(pos)
		   for(size_t func = 0; func < 7; ++func)
		   {
			   #pragma omp atomic
			   sevenHash[func] |= (this->*this->function[func])(read[pos]) << (pos - begin);//(TYPEHASH)boost::multiprecision::pow((boost::multiprecision::uint1024_t)base, pos - begin);
		   }
	   }
}

inline bool GenomicBloomFilter::contains(vector<size_t>& index,
		DynamicBitSet& bit) {
	bool ret = true;
	#pragma omp parallel for
	for(size_t i = 0; i < index.size(); ++i)
		#pragma omp atomic
		ret &= bit.test(index[i]);
	return ret;
}

inline void GenomicBloomFilter::insert(vector<size_t>& index,
		DynamicBitSet& bit) {
	#pragma omp parallel for
	for(size_t i = 0; i < index.size(); ++i)
		bit[index[i]] = 1;
}

//inline GenomicBloomFilter::TYPEHASH GenomicBloomFilter::f_hash(char& ch, int f_type)
//{
//	   switch(f_type)
//	   {
//	   case 0:
//		   if(ch == 'A' || ch == 'C')
//			   return 0;
//		   else
//			   return 1;
//		   break;
//	   case 1:
//		   if(ch == 'A' || ch == 'G')
//			   return 0;
//		   else
//			   return 1;
//		   break;
//	   case 2:
//		   if(ch == 'A' || ch == 'T')
//			   return 0;
//		   else
//			   return 1;
//		   break;
//	   case 3:
//		   if(ch == 'A' )
//			   return 0;
//		   else
//			   return 1;
//		   break;
//	   case 4:
//		   if(ch == 'C')
//			   return 0;
//		   else
//			   return 1;
//		   break;
//	   case 5:
//		   if(ch == 'G')
//			   return 0;
//		   else
//			   return 1;
//		   break;
//	   case 6:
//		   if(ch == 'T')
//			   return 0;
//		   else
//			   return 1;
//		   break;
//	   default:
//		   return 0;
//		   break;
//	   }
//}

inline GenomicBloomFilter::TYPEHASH GenomicBloomFilter::f_hash1(char& ch) {
	   if(ch == 'A' || ch == 'C')
		   return 0;
	   else
		   return 1;
}

inline GenomicBloomFilter::TYPEHASH GenomicBloomFilter::f_hash2(char& ch) {
	   if(ch == 'A' || ch == 'G')
		   return 0;
	   else
		   return 1;
}

inline GenomicBloomFilter::TYPEHASH GenomicBloomFilter::f_hash3(char& ch) {
	   if(ch == 'A' || ch == 'T')
		   return 0;
	   else
		   return 1;
}

inline GenomicBloomFilter::TYPEHASH GenomicBloomFilter::f_hash4(char& ch) {
	   if(ch == 'A' )
		   return 0;
	   else
		   return 1;
}

inline GenomicBloomFilter::TYPEHASH GenomicBloomFilter::f_hash5(char& ch) {
	   if(ch == 'C' )
		   return 0;
	   else
		   return 1;
}

inline GenomicBloomFilter::TYPEHASH GenomicBloomFilter::f_hash6(char& ch) {
	   if(ch == 'G' )
		   return 0;
	   else
		   return 1;
}

inline GenomicBloomFilter::TYPEHASH GenomicBloomFilter::f_hash7(char& ch) {
	   if(ch == 'T' )
		   return 0;
	   else
		   return 1;
}
