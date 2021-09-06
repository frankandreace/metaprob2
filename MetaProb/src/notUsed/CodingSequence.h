/*
 * CodingSequence.h
 *
 *  Created on: 27/nov/2015
 *      Author: samuele
 */

#ifndef SRC_CODINGSEQUENCE_H_
#define SRC_CODINGSEQUENCE_H_

#include "../Utilities.h"
#include <boost/dynamic_bitset.hpp>

using namespace std;

class CodingSequence {
public:
	typedef CodingSequence* Ptr;
	typedef boost::dynamic_bitset<> VectorCode;
	typedef vector<hash_type> VectorHash;
	typedef unordered_set<hash_type> SetHash;

	CodingSequence(); //How many characters that you want in hash > 0
	virtual ~CodingSequence();
	static bool IsOkToBeCoded(string& seq);
	bool Encode(string& seq);
	void Decode(string& seq);

	void GetAllHashes(VectorHash& hashes, size_t HashSize);

private:
	VectorCode code;
	size_t lastHashSize;
	VectorHash hashes;
};

#endif /* SRC_CODINGSEQUENCE_H_ */
