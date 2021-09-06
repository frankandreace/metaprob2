#ifndef	__UTILITIES_CPP__
#define __UTILITIES_CPP__

#include"Utilities.h"

LMerVectorCompress::LMerVectorCompress(size_seq L) {
	this->L = L;
	vector<string> allLmers = GetAllKmers(L);
	for(lMer_type i = 0; i < allLmers.size(); i++)
	{
		Lmer::Ptr l_mer = make_shared<Lmer>();
		l_mer->l_mer = allLmers[i];
		this->vLmer.push_back(l_mer);
	}

	for(lMer_type i = 0; i < this->vLmer.size(); ++i)
	{
		//Sicuramente hashLmer assegnato a indice i
		HashCorrect hashLmer;
		GetHash(this->vLmer[i]->l_mer, 0, L, hashLmer, CharToInt);
		this->mapHash[hashLmer.first] = i;

		//Vediamo dove inserire il reverse, creo hash del reverse
		string revLm = "";
		reverse_copy(this->vLmer[i]->l_mer.begin(), this->vLmer[i]->l_mer.end(), revLm.begin());
		HashCorrect hashReverseLmer;
		GetHash(revLm, 0, L, hashReverseLmer, CharToIntComplement);

		for(lMer_type j = i + 1; j < this->vLmer.size(); ++j)
		{
			HashCorrect hashLmerNext;
			GetHash(this->vLmer[j]->l_mer, 0, L, hashLmerNext, CharToInt);
			if(hashLmerNext.first == hashReverseLmer.first)//if(revLm == this->vLmer[j]->l_mer) //Se il reverse è presente lo cancello e metto tutto su indice hashLmer
			{
				++this->vLmer[i]->count;
				this->vLmer.erase(this->vLmer.begin() + j);
				this->mapHash[hashReverseLmer.first] = i;
				break;
			}
		}
	}
}


const Lmer::Ptr& LMerVectorCompress::GetWithIndex(lMer_type index) {
	return this->vLmer[index];
}

const Lmer::Ptr& LMerVectorCompress::GetWithHash(lMer_type hash) {
	return this->vLmer[this->mapHash.at(hash)];
}

lMer_type LMerVectorCompress::GetIndexWithHash(lMer_type hash) {
	return this->mapHash.at(hash);
}

size_seq LMerVectorCompress::getL() const {
	return L;
}

const LMerVectorCompress::Map_HashLMer_IndexVector& LMerVectorCompress::getMapHash() const {
	return mapHash;
}

const vector<Lmer::Ptr>& LMerVectorCompress::getLmer() const {
	return vLmer;
}

void GetKmer(hash_type index, size_seq K, string& Kmer) {
	if(index >= ((hash_type)1 << K*2))//(hash_type)boost::multiprecision::pow(boost::multiprecision::uint1024_t(4), K)) //Superato indice massimo possibile
		return;
	vector<string> nucleotide = {"A","C","G","T"};
	vector<int> indexNucleotide;
	for(size_seq i = 0; i < K; ++i)
	{
		indexNucleotide.push_back(index%4);
		index /= 4;
	}
	for(size_seq i = K; i > 0; --i)
		Kmer.append(nucleotide[indexNucleotide[i-1]]);
}

vector<string> GetAllKmers(size_seq K) {
	hash_type max = (hash_type)1 << K*2;//(hash_type)boost::multiprecision::pow(boost::multiprecision::uint1024_t(4), K);
	vector<string> allKmer;
	for(hash_type index = 0; index < max; ++index)
	{
		string LMer = "";
		GetKmer(index, K, LMer);
		allKmer.push_back(LMer);
	}
	return allKmer;
}

void createDirAndSubDir(string path) {
	string dir = "";
	string delimiter = "/";
	size_t pos = path.find(delimiter);
	while(pos != path.npos)
	{
		dir += path.substr(0, pos + delimiter.length());
		path.erase(0, pos + delimiter.length());
		mkdir(dir.c_str(), S_IRWXU);
		pos = path.find(delimiter);
	}
}

int parseLineForMemory(char* line){
    int i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i-3] = '\0';
    i = atoi(line);
    return i;
}


int getVirtualMemoryUsed(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLineForMemory(line);
            break;
        }
    }
    fclose(file);
    return result;
}

int getPeakVirtualMemoryUsed() { //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmPeak:", 7) == 0){
            result = parseLineForMemory(line);
            break;
        }
    }
    fclose(file);
    return result;
}

int getPhysicalMemoryUsed(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLineForMemory(line);
            break;
        }
    }
    fclose(file);
    return result;
}

#endif
