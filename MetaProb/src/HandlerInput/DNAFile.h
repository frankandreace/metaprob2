#ifndef	__DNAFILE_H__
#define __DNAFILE_H__ 

#include "../Utilities.h"
#include "../SequencePE.h"
#include "../Group.h"
#include "../HashTable.h"
#include "../BloomFilter/RepeatBloomFilter.h"

class clsDNAFile
{
public:
	typedef shared_ptr<clsDNAFile> Ptr;

	typedef vector<clsSequencePE::Ptr> Vector_Sequence;
	typedef boost::bimaps::bimap<id_seq_type, id_seq_type> Map_IDFir_IDSec;

	clsDNAFile(vector<PairFiles> filename, size_seq Q_SIZE, TypeGraph graph, size_seq m_thres);
	bool GetFileInfo(); //Get info from file, necessary to call after creation
//	bool CreateBloomFilter(); //Create BloomFilter
	bool ExtractQMerFromFile(); //Read fasta file containing paired-end sequences
	bool CreateGraph(); //Create graph thanks to qmers and thresh
	bool CountSeedLmerFromFile(vector<clsGroup::Ptr>& vGroup, LMerVectorCompress::Ptr compr_vLmer);//Read fasta file of reads in groups' seed

	void GetHeaderFromFile(MapIDFile_Header& map_idxFile_Header); //Save header in strings

	Vector_Sequence::iterator begin();
	Vector_Sequence::iterator end();

	Map_IDFir_IDSec::iterator beginMapID();
	Map_IDFir_IDSec::iterator endMapID();
//	id_seq_type GetIDotherEnd(id_seq_type IDAppSeqCurr);

	const id_seq_type GetSize() const;
	const clsSequencePE::Ptr GetAt(id_seq_type ID) const;
	const clsHashTable::Ptr& getHash() const;

	void resetComputation(clsSequencePE::ResetType type);
	void resetHash();

	const id_seq_type getSeqTot() const;
	vector<PairFiles>& getFilename();

private:
	vector<PairFiles> filename;
	size_seq Q_SIZE;
	TypeGraph graph;
	size_seq m_thres;//Number of shared l-mers between reads

	id_seq_type seq_tot; //N. sequence correct readed
	size_seq_tot qMer_Max; //N. of max q-mer present in file

	Vector_Sequence vSeqPE; //ID 0 non esiste, si parte da 1
	Map_IDFir_IDSec map_IDF_IDS; //Mappa ID First end e Second end

	RepeatBloomFilter::Ptr repeatQmer; //q-mer repeat filter
	clsHashTable::Ptr mHash; //Hash delle sequenze nel file

	bool ScanSingleEndFileToCallFunction(PairFiles& fileScan, int (clsDNAFile::*functocall)(PairFiles&, id_seq_type&, string&, string&));
	bool ScanPairedEndFileToCallFunction(PairFiles& fileScan, int (clsDNAFile::*functocall)(PairFiles&, id_seq_type& , string&, string&));

	int storeSequencesInfo(PairFiles& fileScan, id_seq_type& idxFile, string& seq1, string& seq2); //Get n. sequence and max qMer number
//	int extractQMerFromFileForBloomConstruction(PairFiles& fileScan, id_seq_type& idxFile, string& seq1, string& seq2); //Create BloomFilter
	int extractQMerFromSequence(PairFiles& fileScan, id_seq_type& idxFile, string& seq1, string& seq2);
	int createGraph(PairFiles& fileScan, id_seq_type& idxFile, string& seq1, string& seq2);

	bool CountSeedLmerFromSingleEndFile(PairFiles& fileScan, vector<clsGroup::Ptr>& vGroup, LMerVectorCompress::Ptr compr_vLmer);//Read fna file of reads in groups' seed
	bool CountSeedLmerFromPairedEndFile(PairFiles& fileScan, vector<clsGroup::Ptr>& vGroup, LMerVectorCompress::Ptr compr_vLmer);//Read fastq file of reads in groups' seed
	int countLMerIfSeed(PairFiles& fileScan, id_seq_type& idxFile, string& seq1, string& seq2, vector<clsGroup::Ptr>& vGroup, LMerVectorCompress::Ptr compr_vLmer);
	void count_l_mer(string& seq, id_grp_type id_grp, vector<clsGroup::Ptr>& vGroup, LMerVectorCompress::Ptr compr_vLmer); //count l-mer of end string in vGroup[id_group]

	void parserHeader(string& header);
	int storeHeader(PairFiles& fileScan, id_seq_type& idxFile, string& seq1, string& seq2,
			string& header1, string& header2, MapIDFile_Header& map_idxFile_Header);
	void GetHeaderFromSingleEndFile(PairFiles& fileScan, MapIDFile_Header& map_idxFile_Header);
	void GetHeaderFromPairedEndFile(PairFiles& fileScan, MapIDFile_Header& map_idxFile_Header);
};
#endif
