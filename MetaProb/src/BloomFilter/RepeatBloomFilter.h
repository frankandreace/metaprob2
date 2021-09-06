#ifndef SRC_REPEATBLOOMFILTER_H_
#define SRC_REPEATBLOOMFILTER_H_

#include "GenomicBloomFilter.h"

class RepeatBloomFilter: public GenomicBloomFilter {
public:
	typedef shared_ptr<RepeatBloomFilter> Ptr;

	DynamicBitSet repeat;

	RepeatBloomFilter(TYPEHASH num_element, double prob_false_positive, size_t Q, double perc_repeat_element);
	virtual ~RepeatBloomFilter();

	void areRepeats(string& read, vector<bool>& posCont);
	void insertForRepeatSearch(string& read);
};

#endif /* SRC_REPEATBLOOMFILTER_H_ */
