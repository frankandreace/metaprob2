/*
 * RepeatBloomFilter.cpp
 *
 *  Created on: 18/nov/2015
 *      Author: samuele
 */

#include "RepeatBloomFilter.h"

RepeatBloomFilter::RepeatBloomFilter(TYPEHASH num_element, double prob_false_positive,
		size_t Q, double perc_repeat_element) : GenomicBloomFilter(num_element, prob_false_positive, Q) {
	const TYPEHASH size_repeat = -(num_element*perc_repeat_element)*log(prob_false_positive)/pow(log(2), 2);
	DynamicBitSet(size_repeat).swap(this->repeat);
}

RepeatBloomFilter::~RepeatBloomFilter() {
}

void RepeatBloomFilter::areRepeats(string& read, vector<bool>& posCont) {
	vector<vector<TYPEHASH>> vSevenHashes;
	GenomicBloomFilter::GetHashes(read, vSevenHashes);
	posCont = vector<bool>(vSevenHashes.size(), false);
	#pragma omp parallel for
	for(size_t pos = 0; pos < vSevenHashes.size(); ++pos)
	{
		vector<size_t> index;
		GenomicBloomFilter::GetPositionWithHash(vSevenHashes[pos], index, this->repeat);
		posCont[pos] = GenomicBloomFilter::contains(index, this->repeat);
	}
}

void RepeatBloomFilter::insertForRepeatSearch(string& read) {
	vector<vector<TYPEHASH>> vSevenHashes;
	GenomicBloomFilter::GetHashes(read, vSevenHashes);
	#pragma omp parallel for ordered
	for(size_t pos = 0; pos < vSevenHashes.size(); ++pos)
	{
		vector<size_t> index_bitset;
		GenomicBloomFilter::GetPositionWithHash(vSevenHashes[pos], index_bitset, this->bitset);
		#pragma omp ordered
		{
			bool already_contains = GenomicBloomFilter::contains(index_bitset, this->bitset);
			if(already_contains)
			{
				vector<size_t> index_repeat;
				GenomicBloomFilter::GetPositionWithHash(vSevenHashes[pos], index_repeat, this->repeat);
				GenomicBloomFilter::insert(index_repeat, this->repeat);
			}
			else
				GenomicBloomFilter::insert(index_bitset, this->bitset);
		}
	}
}
