#ifndef	__VIRSEQ_H__
#define __VIRSEQ_H__ 

#include "Utilities.h"

using namespace std;

class clsVirSeq
{
public:
	virtual id_seq_type getSeqId() const;
	clsVirSeq(id_seq_type iSeqID);
	virtual ~clsVirSeq();

private:
	id_seq_type iSeqID;

};
#endif
