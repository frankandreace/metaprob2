#include "DNAFile.h"

clsDNAFile::clsDNAFile(vector<PairFiles> filename, size_seq Q_SIZE, TypeGraph graph, size_seq m_thres):
filename(filename), Q_SIZE(Q_SIZE), graph(graph), m_thres(m_thres)
{
	this->seq_tot = 0;
	this->qMer_Max = 0;
}

bool clsDNAFile::GetFileInfo() {
	cout << endl << "Loading Sequences... " << flush;

	//Initialize Map
	Map_IDFir_IDSec().swap(this->map_IDF_IDS);

	bool ret = false;
	for(size_t i = 0; i < this->filename.size(); ++i)
	{
		if(!this->filename[i].isCorrect())
			continue;
		if(this->filename[i].getType() == SingleEnd)
			ret = this->ScanSingleEndFileToCallFunction(this->filename[i], &clsDNAFile::storeSequencesInfo);
		else if(this->filename[i].getType() == PairedEnd)
			ret = this->ScanPairedEndFileToCallFunction(this->filename[i], &clsDNAFile::storeSequencesInfo);
	}
	cout << "Complete" << flush;
	cout << endl << "Loaded sequences: " + to_string(this->getSeqTot()) << flush;
	return ret;
}

////Create BloomFilter in from file
//bool clsDNAFile::CreateBloomFilter() {
//	cout << endl << "Create BloomFilter... " << flush;
//
//	this->repeatQmer = make_shared<RepeatBloomFilter>(this->qMer_Max, 0.001, this->Q_SIZE, 0.6);
//
//	bool ret = false;
//	for(size_t i = 0; i < this->filename.size(); ++i)
//	{
//		if(!this->filename[i].isCorrect())
//			continue;
//		if(this->filename[i].getType() == SingleEnd)
//			ret = this->ScanSingleEndFileToCallFunction(this->filename[i], &clsDNAFile::extractQMerFromFileForBloomConstruction);
//		else if(this->filename[i].getType() == PairedEnd)
//			ret = this->ScanPairedEndFileToCallFunction(this->filename[i], &clsDNAFile::extractQMerFromFileForBloomConstruction);
//	}
//	cout << "Complete" << flush;
//	return ret;
//}

bool clsDNAFile::ExtractQMerFromFile()
{
	//Initialize MapHash
	this->mHash = make_shared<clsHashTable>();
	this->mHash->InitOptionData(this->Q_SIZE, this->qMer_Max*0.8);

	cout << endl << "Extract q-mer from Sequences... " << flush;

	bool ret = false;
	for(size_t i = 0; i < this->filename.size(); ++i)
	{
		if(!this->filename[i].isCorrect())
			continue;
		if(this->filename[i].getType() == SingleEnd)
			ret = this->ScanSingleEndFileToCallFunction(this->filename[i], &clsDNAFile::extractQMerFromSequence);
		else if(this->filename[i].getType() == PairedEnd)
			ret = this->ScanPairedEndFileToCallFunction(this->filename[i], &clsDNAFile::extractQMerFromSequence);
		else
			ret = false;
	}
	cout << "Complete" << flush;
	return ret;
}

bool clsDNAFile::CreateGraph()
{
	cout << endl << "Building Graph..." << flush;

	bool ret = false;
	#pragma omp parallel for
	for(size_t i = 0; i < this->filename.size(); ++i)
	{
		if(!this->filename[i].isCorrect())
			continue;
		if(this->filename[i].getType() == SingleEnd)
			ret = this->ScanSingleEndFileToCallFunction(this->filename[i], &clsDNAFile::createGraph);
		else if(this->filename[i].getType() == PairedEnd)
			ret = this->ScanPairedEndFileToCallFunction(this->filename[i], &clsDNAFile::createGraph);
		else
			ret = false;
	}
	cout << "Complete" << flush;
	return ret;
}

bool clsDNAFile::CountSeedLmerFromFile(vector<clsGroup::Ptr>& vGroup, LMerVectorCompress::Ptr compr_vLmer) {
	cout << endl << "Count " << compr_vLmer->getL() << "-mer from file..." << flush;

	bool ret = false;
	for(size_t i = 0; i < this->filename.size(); ++i)
	{
		if(!this->filename[i].isCorrect())
			continue;
		if(this->filename[i].getType() == SingleEnd)
		ret = this->CountSeedLmerFromSingleEndFile(this->filename[i], vGroup, compr_vLmer);
		else if(this->filename[i].getType() == PairedEnd)
		ret = this->CountSeedLmerFromPairedEndFile(this->filename[i], vGroup, compr_vLmer);
		else
			ret = false;
	}
		
	cout << " Complete" << flush;
	return ret;
}

clsDNAFile::Vector_Sequence::iterator clsDNAFile::begin() {
	return this->vSeqPE.begin();
}

clsDNAFile::Vector_Sequence::iterator clsDNAFile::end() {
	return this->vSeqPE.end();
}


clsDNAFile::Map_IDFir_IDSec::iterator clsDNAFile::beginMapID() {
	return this->map_IDF_IDS.begin();
}

clsDNAFile::Map_IDFir_IDSec::iterator clsDNAFile::endMapID() {
	return this->map_IDF_IDS.end();
}

//id_seq_type clsDNAFile::GetIDotherEnd(id_seq_type IDAppSeqCurr) {
//	auto it_find_fir = this->map_IDF_IDS.left.find(IDAppSeqCurr);
//	if(it_find_fir != this->map_IDF_IDS.left.end()) //Trovato a sx
//		return it_find_fir->second;
//	auto it_find_sec = this->map_IDF_IDS.right.find(IDAppSeqCurr);
//	if (it_find_sec != this->map_IDF_IDS.right.end()) //Trovato a dx
//		return it_find_sec->first;
//	return 0;
//}

const id_seq_type clsDNAFile::GetSize() const{
	return this->vSeqPE.size();
}

const id_seq_type clsDNAFile::getSeqTot() const {
	return seq_tot;
}

const clsSequencePE::Ptr clsDNAFile::GetAt(id_seq_type ID) const {
	 try
	 {
		 clsSequencePE::Ptr seqSaved = this->vSeqPE.at(ID-1);
		 return seqSaved;      // vector::at throws an out-of-range
		 }
	 catch (const std::out_of_range& oor) {
		 return NULL;
	 }
}

const clsHashTable::Ptr& clsDNAFile::getHash() const {
	return mHash;
}

vector<PairFiles>& clsDNAFile::getFilename() {
	return filename;
}

void clsDNAFile::resetComputation(clsSequencePE::ResetType type)
{
	auto it_ID = this->begin();
	while(it_ID != this->end())
	{
		(*it_ID)->reset(type);
		++it_ID;
	}
}

void clsDNAFile::resetHash() {
	this->mHash.reset();
	this->mHash = make_shared<clsHashTable>();
}

//General scan of FnaFile
inline bool clsDNAFile::ScanSingleEndFileToCallFunction(PairFiles& fileScan, int (clsDNAFile::*functocall)(PairFiles&, id_seq_type&, string&, string&))
{
	string empty = "";
	string line; //store the read data
	ifstream file(fileScan.first.getPath());
	if(file.is_open())
	{
		string seq = "";
		id_seq_type idxFile = 0;
		bool isLineQuality = false;
		while(getline(file,line))
		{
			if(!line.empty())
			{
				if(line[0] == fileScan.first.getReadDelimiter())
				{
					//due casi:
					//o delimitatore non fa parte di qualità se isLineQuality == false, quindi line contiene la descrizione del read
					//o delimitatore fa parte di qualità se isLineQuality == true e line contiene la riga qualità
					if(!isLineQuality && seq != "")
						(this->*functocall)(fileScan, idxFile, seq, empty);//Store previus data
					//altrimenti line contiene la sequenza qualità che inizia per il delimitatore, quindi salta
					isLineQuality = false;
				}
				else
				{
					//caso line non inizia per il delimitatore tre casi:
					//o line contiene il + quindi prossima line sarà qualità, si imposta isLineQuality == true
					//o line contiene la line di qualità quindi isLineQuality == true
					//o line contiene la line di sequenza quindi isLineQuality == false
					if(fileScan.first.getFileType() == Fastq && line[0] == '+')
						isLineQuality = true; //prossima line sarà linea qualità
					else
					{
						if(!isLineQuality) //Non è line qualità quindi ha la sequenza da memorizzare
							seq.append(line);
						isLineQuality = false;
					}
				}
			}
		}
		(this->*functocall)(fileScan, idxFile, seq, empty);//Store last data
		#pragma omp taskwait
		file.close();
		return true;
	}
	else
	{
		cout<<"Cannot open input Fasta file "<< fileScan.first.getPath();
		return false;
	}
}

//General scan of FastqFile
inline bool clsDNAFile::ScanPairedEndFileToCallFunction(PairFiles& fileScan, int (clsDNAFile::*functocall)(PairFiles&, id_seq_type&, string&, string&))
{
	string line1; //store the read data1
	ifstream file1(fileScan.first.getPath());
	string line2; //store the read data2
	ifstream file2(fileScan.second.getPath());
	if(file1.is_open() && file2.is_open())
	{
		string seq1 = "";
		string seq2 = "";
		id_seq_type idxFile = 0;
		bool isLineQuality = false;
		while(getline(file1,line1) && getline(file2,line2))
		{
			if(!line1.empty() && !line2.empty())
			{
				if(line1[0] == fileScan.first.getReadDelimiter() && line2[0] == fileScan.second.getReadDelimiter())
				{
					//due casi:
					//o delimitatore non fa parte di qualità se isLineQuality == false, quindi line contiene la descrizione del read
					//o delimitatore fa parte di qualità se isLineQuality == true e line contiene la riga qualità
					if(!isLineQuality && seq1 != "" && seq2 != "")
						(this->*functocall)(fileScan, idxFile, seq1, seq2);//Store previus data
					//altrimenti line contiene la sequenza qualità che inizia per il delimitatore, quindi salta
					isLineQuality = false;
				}
				else
				{
					//caso line non inizia per il delimitatore tre casi:
					//o line contiene il + quindi prossima line sarà qualità, si imposta isLineQuality == true
					//o line contiene la line di qualità quindi isLineQuality == true
					//o line contiene la line di sequenza quindi isLineQuality == false
					if(fileScan.first.getFileType() == Fastq && line1[0] == '+' && fileScan.second.getFileType() == Fastq && line2[0] == '+')
						isLineQuality = true; //prossima line sarà linea qualità
					else
					{
						if(!isLineQuality) //Non è line qualità quindi ha la sequenza da memorizzare
						{
							seq1.append(line1);
							seq2.append(line2);
						}
						isLineQuality = false;
					}
				}
			}
		}
		(this->*functocall)(fileScan, idxFile, seq1, seq2);//Store previus data
		#pragma omp taskwait
		file1.close();
		file2.close();
		return true;
	}
	else
	{
		if(!file1.is_open())
			cout<<"Cannot open input Fasta file "<< fileScan.first.getPath();
		if(!file2.is_open())
			cout<<"Cannot open input Fasta file "<< fileScan.second.getPath();
		return false;
	}
}

//Read fasta file and insert sequence in vector, do not store reads in memory,
//only the presence and the IDFile and ID in Application
inline int clsDNAFile::storeSequencesInfo(PairFiles& fileScan, id_seq_type& idxFile, string& seq1, string& seq2)
{
	++idxFile; //Incrementa ad ogni chiamata, lo 0 non esiste, parte da 1 così facendo
	clsSequencePE::Ptr mySeqPE1;
	clsSequencePE::Ptr mySeqPE2;

	mySeqPE1 = make_shared<clsSequencePE>(fileScan, this->vSeqPE.size()+1);//ID APP. L'elemento 0 non deve esistere
	mySeqPE1->AddSequence(seq1, 1);
	this->vSeqPE.push_back(mySeqPE1);

	++this->seq_tot;
	this->qMer_Max += seq1.length() - this->Q_SIZE + 1;
	fileScan.first.AddReadMap(idxFile, mySeqPE1->getSeqId()); //Mappa idxFile in idxApp

	if(fileScan.getType() == PairedEnd)
	{
		if(this->graph == Single || this->graph == SingleUnion)//Salva singoli
		{
			//Aggiungo anche seq2 in un altra sequenza
			mySeqPE2 = make_shared<clsSequencePE>(fileScan, this->vSeqPE.size()+1);//ID APP. L'elemento 0 non deve esistere
			mySeqPE2->AddSequence(seq2, 2);
			this->vSeqPE.push_back(mySeqPE2);

			++this->seq_tot;
			this->qMer_Max += seq2.length() - this->Q_SIZE + 1;

			if(this->graph == SingleUnion) //Salva link tra due paired per unione
				this->map_IDF_IDS.insert(Map_IDFir_IDSec::value_type(mySeqPE1->getSeqId(), mySeqPE2->getSeqId()));
		}
		else if (this->graph == Paired) //Salva a coppie
		{
			//Aggiungo anche seq2 stessa sequenza
			mySeqPE2 = mySeqPE1;
			mySeqPE2->AddSequence(seq2, 2);

			++this->seq_tot;
			this->qMer_Max += seq2.length() - this->Q_SIZE + 1;
		}
		fileScan.second.AddReadMap(idxFile, mySeqPE2->getSeqId());//Mappa idxFile in idxApp
	}
	//reset for new sequence
	seq1 = "";
	seq2 = "";
	return 0;
}

////Read fasta file and extract list of q-mers to create BloomFilter
//inline int clsDNAFile::extractQMerFromFileForBloomConstruction(PairFiles& fileScan, id_seq_type& idxFile, string& seq1, string& seq2)
//{
//	++idxFile;
//	if(seq1.size() >= this->Q_SIZE)
//		this->repeatQmer->insertForRepeatSearch(seq1);
//	if(seq2.size() >= this->Q_SIZE)
//		this->repeatQmer->insertForRepeatSearch(seq2);
//	//reset for new sequence
//	seq1 = "";
//	seq2 = "";
//	return 0;
//}

//Read fasta file and extract list of q-mers, do not store reads in memory
//L'elemento 0 non deve esistere, ID applicazione assegnato con this->vSeqPE.size()+1 per evitare di avere 0
inline int clsDNAFile::extractQMerFromSequence(PairFiles& fileScan, id_seq_type& idxFile, string& seq1, string& seq2)
{
	++idxFile;
	id_seq_type ID1 = fileScan.first.getIdxApp(idxFile);
	id_seq_type ID2 = fileScan.second.getIdxApp(idxFile);

	//Estrai q-mers
	if(ID1 != 0)
		this->mHash->ExtractFromString(ID1, seq1, this->repeatQmer);
	if(ID2 != 0)
		this->mHash->ExtractFromString(ID2, seq2, this->repeatQmer);

//		//Stampa
//		cout << "\r" << "Sequence: " << this->vSeqPE->size() << flush;

	//reset for new sequence
	seq1 = "";
	seq2 = "";
	return 0;
}

inline int clsDNAFile::createGraph(PairFiles& fileScan, id_seq_type& idxFile, string& _seq1, string& _seq2)
{
	++idxFile;
	id_seq_type ID1 = fileScan.first.getIdxApp(idxFile);
	id_seq_type ID2 = fileScan.second.getIdxApp(idxFile);

	//Stampa
	//cout << "\r" << "Building graph... Analyzed sequences: " << idxFile << " " << flush;

	string seq1 = _seq1;
	string seq2 = _seq2;
	auto hash = [this](string& str, unordered_set<hash_type>& set_hash){
		vector<HashCorrect> vHash;
		GetHashes(str, this->Q_SIZE, vHash, CharToInt);

		//Insert q-mer in hash and save ID read where is present
		for(size_seq i = 0; i < vHash.size(); ++i)
			if(vHash[i].second)
				set_hash.insert(vHash[i].first);
	};
	auto adj = [this](id_seq_type ID, unordered_set<hash_type>& set_hash) {
		clsSequencePE::Ptr seq = this->GetAt(ID);
		clsSequencePE::MapAdj tmp_mapAdj;
		auto it_qmer_seq = set_hash.begin();
		while(it_qmer_seq != set_hash.end())
		{
			KmNode& kmnode = this->mHash->arrMyHash[*it_qmer_seq];
			MapAdjSeqID adj_for_qmer;
			GetMapSeqID(adj_for_qmer, kmnode);
			if(adj_for_qmer.size() > 1)
			{
				size_seq& size_qmer = adj_for_qmer[ID];
				auto it_adj = adj_for_qmer.begin();
				while(it_adj != adj_for_qmer.end())
				{
					if(ID != it_adj->first)
					{
						size_seq_tot count = (size_seq_tot)size_qmer * (size_seq_tot)it_adj->second;
						tmp_mapAdj[it_adj->first] += count;
					}
					++it_adj;
				}
			}
			++it_qmer_seq;
		}
		auto it_final_adj = tmp_mapAdj.begin();
		while(it_final_adj != tmp_mapAdj.end())
		{
			if(it_final_adj->second >= this->m_thres)
			{
				clsSequencePE::Ptr seq_adj = this->GetAt(it_final_adj->first);
				seq->AddAdjVertice(seq_adj, it_final_adj->second);
			}
			++it_final_adj;
		}
	};
	#pragma omp task firstprivate (ID1, ID2, seq1, seq2)
	{
		if(ID1 != 0 && ID1 == ID2)
		{
			unordered_set<hash_type> set_hash;
			hash(seq1, set_hash);
			hash(seq2, set_hash);
			adj(ID1, set_hash);
		}
		else
		{
			//Estrai q-mers
			if(ID1 != 0)
			{
				unordered_set<hash_type> set_hash;
				hash(seq1, set_hash);
				adj(ID1, set_hash);
			}
			if(ID2 != 0)
			{
				unordered_set<hash_type> set_hash;
				hash(seq2, set_hash);
				adj(ID2, set_hash);
			}
		}
	}

	//reset for new sequence
	_seq1 = "";
	_seq2 = "";
	return 0;
}

//Read fasta file of reads in groups' seed
//(funziona sia con Paired End che con Single End)
bool clsDNAFile::CountSeedLmerFromSingleEndFile(PairFiles& fileScan, vector<clsGroup::Ptr>& vGroup, LMerVectorCompress::Ptr compr_vLmer)
{
	string line; //store the read data
	ifstream file(fileScan.first.getPath());
	if(file.is_open())
	{
		string empty = "";
		string seq = "";
		id_seq_type idxFile = 0;
		bool isLineQuality = false;
		while(getline(file,line))
		{
			if(!line.empty())
			{
				if(line[0] == fileScan.first.getReadDelimiter())
				{
					//due casi:
					//o delimitatore non fa parte di qualità se isLineQuality == false, quindi line contiene la descrizione del read
					//o delimitatore fa parte di qualità se isLineQuality == true e line contiene la riga qualità
					if(!isLineQuality && seq != "")
						this->countLMerIfSeed(fileScan, idxFile, seq, empty, vGroup, compr_vLmer);//Store previus data
					//altrimenti line contiene la sequenza qualità che inizia per il delimitatore, quindi salta
					isLineQuality = false;
				}
				else
				{
					//caso line non inizia per il delimitatore tre casi:
					//o line contiene il + quindi prossima line sarà qualità, si imposta isLineQuality == true
					//o line contiene la line di qualità quindi isLineQuality == true
					//o line contiene la line di sequenza quindi isLineQuality == false
					if(fileScan.first.getFileType() == Fastq && line[0] == '+')
						isLineQuality = true; //prossima line sarà linea qualità
					else
					{
						if(!isLineQuality) //Non è line qualità quindi ha la sequenza da memorizzare
							seq.append(line);
						isLineQuality = false;
					}
				}
			}
		}
		this->countLMerIfSeed(fileScan, idxFile, seq, empty, vGroup, compr_vLmer);//Store last data
		file.close();
		return true;
	}
	else
	{
		cout<<"Cannot open input Fasta file "<< fileScan.first.getPath();
		return false;
	}
}

//Read fasta file of reads in groups' seed
//(funziona sia con Paired End che con Single End)
bool clsDNAFile::CountSeedLmerFromPairedEndFile(PairFiles& fileScan, vector<clsGroup::Ptr>& vGroup, LMerVectorCompress::Ptr compr_vLmer)
{
	string line1; //store the read data1
	ifstream file1(fileScan.first.getPath());
	string line2; //store the read data2
	ifstream file2(fileScan.second.getPath());
	if(file1.is_open() && file2.is_open())
	{
		string seq1 = "";
		string seq2 = "";
		id_seq_type idxFile = 0;
		bool isLineQuality = false;
		while(getline(file1,line1) && getline(file2,line2))
		{
			if(!line1.empty() && !line2.empty())
			{
				if(line1[0] == fileScan.first.getReadDelimiter() && line2[0] == fileScan.second.getReadDelimiter())
				{
					//due casi:
					//o delimitatore non fa parte di qualità se isLineQuality == false, quindi line contiene la descrizione del read
					//o delimitatore fa parte di qualità se isLineQuality == true e line contiene la riga qualità
					if(!isLineQuality && seq1 != "" && seq2 != "")
						this->countLMerIfSeed(fileScan, idxFile, seq1, seq2, vGroup, compr_vLmer);//Store previus data
					//altrimenti line contiene la sequenza qualità che inizia per il delimitatore, quindi salta
					isLineQuality = false;
				}
				else
				{
					//caso line non inizia per il delimitatore tre casi:
					//o line contiene il + quindi prossima line sarà qualità, si imposta isLineQuality == true
					//o line contiene la line di qualità quindi isLineQuality == true
					//o line contiene la line di sequenza quindi isLineQuality == false
					if(fileScan.first.getFileType() == Fastq && line1[0] == '+' && fileScan.second.getFileType() == Fastq && line2[0] == '+')
						isLineQuality = true; //prossima line sarà linea qualità
					else
					{
						if(!isLineQuality) //Non è line qualità quindi ha la sequenza da memorizzare
						{
							seq1.append(line1);
							seq2.append(line2);
						}
						isLineQuality = false;
					}
				}
			}
		}
		this->countLMerIfSeed(fileScan, idxFile, seq1, seq2, vGroup, compr_vLmer);//Store previus data
		file1.close();
		file2.close();
		return true;
	}
	else
	{
		if(!file1.is_open())
			cout<<"Cannot open input Fasta file "<< fileScan.first.getPath();
		if(!file2.is_open())
			cout<<"Cannot open input Fasta file "<< fileScan.second.getPath();
		return false;
	}
}

inline int clsDNAFile::countLMerIfSeed(PairFiles& fileScan, id_seq_type& idxFile, string& seq1, string& seq2,
		vector<clsGroup::Ptr>& vGroup, LMerVectorCompress::Ptr compr_vLmer)
{
	++idxFile;
	id_seq_type ID1 = fileScan.first.getIdxApp(idxFile);
	id_seq_type ID2 = fileScan.second.getIdxApp(idxFile);
	clsSequencePE::Ptr seqS1 = this->GetAt(ID1);
	clsSequencePE::Ptr seqS2 = this->GetAt(ID2);
	if(seqS1 != NULL && seqS1->getSeed() == Seed)
	{
		id_grp_type id_grp_1 = seqS1->getIdGrp();
		this->count_l_mer(seq1, id_grp_1, vGroup, compr_vLmer);
	}
	if(seqS2 != NULL && seqS2->getSeed() == Seed)
	{
		id_grp_type id_grp_2 = seqS2->getIdGrp();
		this->count_l_mer(seq2, id_grp_2, vGroup, compr_vLmer);
	}
	//reset for new sequence
	seq1 = "";
	seq2 = "";
	return 0;
}

inline void clsDNAFile::count_l_mer(string& seq, id_grp_type id_grp, vector<clsGroup::Ptr>& vGroup,
		LMerVectorCompress::Ptr compr_vLmer)
{
	size_seq l_mer_size = compr_vLmer->getL();
	if(seq.length() >= l_mer_size)
	{
		vector<HashCorrect> hashes;
		GetHashes(seq, l_mer_size, hashes, CharToInt);
		for(size_seq i = 0; i < hashes.size(); ++i)
			if(hashes[i].second)
				++vGroup[id_grp]->vCountMer[compr_vLmer->GetIndexWithHash(hashes[i].first)];
	}
}

void clsDNAFile::GetHeaderFromFile(MapIDFile_Header& map_idxFile_Header)
{
	for(size_t i = 0; i < this->filename.size(); ++i)
	{
		if(!this->filename[i].isCorrect())
			continue;
		if(this->filename[i].getType() == SingleEnd)
			return this->GetHeaderFromSingleEndFile(this->filename[i], map_idxFile_Header);
		else if(this->filename[i].getType() == PairedEnd)
			return this->GetHeaderFromPairedEndFile(this->filename[i], map_idxFile_Header);
	}
}

inline void clsDNAFile::parserHeader(string& header)
{
	vector<string> delimiter = {"\t", "\n", " "};
	size_t min_pos = numeric_limits<size_t>::max();
	for(size_t i = 0; i < delimiter.size(); ++i)
	{
		size_t tmp_pos = header.find(delimiter[i]);
		if(min_pos > tmp_pos)
			min_pos = tmp_pos;
	}
	header = header.substr(0,min_pos);
}

inline int clsDNAFile::storeHeader(PairFiles& fileScan, id_seq_type& idxFile, string& seq1, string& seq2,
		string& header1, string& header2, MapIDFile_Header& map_idxFile_Header)
{
	++idxFile;
	MapIDFile_Header::mapped_type& pairheader = map_idxFile_Header[idxFile];

	//Parse header
	this->parserHeader(header1);
	this->parserHeader(header2);

	pairheader.first = header1;
	pairheader.second = header2;

	//reset for new sequence
	seq1 = "";
	seq2 = "";
	header1 = "";
	header2 = "";
	return 0;
}

void clsDNAFile::GetHeaderFromSingleEndFile(PairFiles& fileScan, MapIDFile_Header& map_idxFile_Header) {
	string line; //store the read data
	ifstream file(fileScan.first.getPath());
	if(file.is_open())
	{
		string empty = "";
		string seq = "";
		string header = "";
		id_seq_type idxFile = 0;
		bool isLineQuality = false;
		while(getline(file,line))
		{
			if(!line.empty())
			{
				if(line[0] == fileScan.first.getReadDelimiter())
				{
					//due casi:
					//o delimitatore non fa parte di qualità se isLineQuality == false, quindi line contiene la descrizione del read
					//o delimitatore fa parte di qualità se isLineQuality == true e line contiene la riga qualità
					if(!isLineQuality)
					{
						if(seq != "")
							this->storeHeader(fileScan, idxFile, seq, empty, header, empty, map_idxFile_Header);//Store previus data
						header = line; //Store next header
					}
					//altrimenti line contiene la sequenza qualità che inizia per il delimitatore, quindi salta
					isLineQuality = false;
				}
				else
				{
					//caso line non inizia per il delimitatore tre casi:
					//o line contiene il + quindi prossima line sarà qualità, si imposta isLineQuality == true
					//o line contiene la line di qualità quindi isLineQuality == true
					//o line contiene la line di sequenza quindi isLineQuality == false
					if(fileScan.first.getFileType() == Fastq && line[0] == '+')
						isLineQuality = true; //prossima line sarà linea qualità
					else
					{
						if(!isLineQuality) //Non è line qualità quindi ha la sequenza da memorizzare
							seq.append(line);
						isLineQuality = false;
					}
				}
			}
		}
		this->storeHeader(fileScan, idxFile, seq, empty, header, empty, map_idxFile_Header);//Store last data
		file.close();
	}
	else
	{
		cout<<"Cannot open input Fasta file "<< fileScan.first.getPath();
	}
}

void clsDNAFile::GetHeaderFromPairedEndFile(PairFiles& fileScan, MapIDFile_Header& map_idxFile_Header) {
	string line1; //store the read data1
	ifstream file1(fileScan.first.getPath());
	string line2; //store the read data2
	ifstream file2(fileScan.second.getPath());
	if(file1.is_open() && file2.is_open())
	{
		string seq1 = "";
		string seq2 = "";
		string header1 = "";
		string header2 = "";
		id_seq_type idxFile = 0;
		bool isLineQuality = false;
		while(getline(file1,line1) && getline(file2,line2))
		{
			if(!line1.empty() && !line2.empty())
			{
				if(line1[0] == fileScan.first.getReadDelimiter() && line2[0] == fileScan.second.getReadDelimiter())
				{
					//due casi:
					//o delimitatore non fa parte di qualità se isLineQuality == false, quindi line contiene la descrizione del read
					//o delimitatore fa parte di qualità se isLineQuality == true e line contiene la riga qualità
					if(!isLineQuality)
					{
						if(seq1 != "" && seq2 != "")
							this->storeHeader(fileScan, idxFile, seq1, seq2, header1, header2, map_idxFile_Header);//Store previus data
						header1 = line1;
						header2 = line2;
					}
					//altrimenti line contiene la sequenza qualità che inizia per il delimitatore, quindi salta
					isLineQuality = false;
				}
				else
				{
					//caso line non inizia per il delimitatore tre casi:
					//o line contiene il + quindi prossima line sarà qualità, si imposta isLineQuality == true
					//o line contiene la line di qualità quindi isLineQuality == true
					//o line contiene la line di sequenza quindi isLineQuality == false
					if(fileScan.first.getFileType() == Fastq && line1[0] == '+' && fileScan.second.getFileType() == Fastq && line2[0] == '+')
						isLineQuality = true; //prossima line sarà linea qualità
					else
					{
						if(!isLineQuality) //Non è line qualità quindi ha la sequenza da memorizzare
						{
							seq1.append(line1);
							seq2.append(line2);
						}
						isLineQuality = false;
					}
				}
			}
		}
		this->storeHeader(fileScan, idxFile, seq1, seq2, header1, header2, map_idxFile_Header);//Store previus data
		file1.close();
		file2.close();
	}
	else
	{
		if(!file1.is_open())
			cout<<"Cannot open input Fasta file "<< fileScan.first.getPath();
		if(!file2.is_open())
			cout<<"Cannot open input Fasta file "<< fileScan.second.getPath();
	}
}
