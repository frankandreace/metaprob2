#include"VirSeq.h"

id_seq_type clsVirSeq::getSeqId() const {
	return iSeqID;
}

clsVirSeq::clsVirSeq(id_seq_type iSeqID) {
	this->iSeqID = iSeqID;
}

clsVirSeq::~clsVirSeq(){}
