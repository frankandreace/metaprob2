#ifndef	__SEQGRAPH_H__
#define __SEQGRAPH_H__ 

#include"Utilities.h"
#include"HandlerInput/DNAFile.h"
#include"HashTable.h"

class clsSeqGraph
{
public:
	typedef shared_ptr<clsSeqGraph> Ptr;

	clsSeqGraph(clsDNAFile::Ptr mFile);
	bool FilteredQmer(MapAdjSeqID& map); //Filtraggio q-mer troppo condivisi tra diversi read (non rappresentano sovrapposizioni alle estremità dei read)
	void CreateGraph();//Create graph from hashtable of l-mers
	void ResetAdj(); //Remove Adj
	const clsDNAFile::Ptr& getFile() const;
	hash_type getQMerFilter() const;

private:
	clsDNAFile::Ptr mFile;
	hash_type qMerFilter;
};
#endif
