#ifndef	__UTILITIES_H__
#define __UTILITIES_H__ 

#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <forward_list>
#include <chrono>
#include <cstdint>
#include <memory>
#include <unordered_set>
#include <map>
#include <unordered_map>
//#include <boost/multiprecision/cpp_int.hpp>
#include <cmath>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <boost/bimap/bimap.hpp> //Bidirectional map

using namespace std;

using namespace std;
using namespace std::chrono;

//Type definition for all project
typedef uint16_t size_seq; //max 2^16=65.536
typedef uint64_t size_seq_tot; //max 2^64=1,844674407×10¹⁹

typedef uint32_t id_seq_type;
typedef uint16_t id_specie_type;
typedef uint32_t id_grp_type;
typedef uint16_t id_cluster_type;

//typedef uint128_t hash_type; //max 2^128=4^64=3,402823669×10³⁸
typedef uint64_t hash_type; //max 2^64=1,844674407×10¹⁹
typedef pair<hash_type, bool> HashCorrect;
typedef uint16_t lMer_type;

//For sequences
enum SeedState {No_State_Seed, No_Seed, Seed, No_Other_Read};
typedef unordered_map<id_seq_type, size_seq> MapAdjSeqID;
enum TypeGraph {Paired = 0, Single, SingleUnion};

typedef pair<string, string> SequenceHeader;
typedef map<id_seq_type, SequenceHeader> MapIDFile_Header;

//////////////////////////////////////////////////////////////////////////////////////////////////////
//Necessari per compilazione codice Assessment.h e Assessment.cpp
typedef map<id_specie_type, id_seq_type> Map_Specie__NumRead;

typedef map<id_grp_type, id_seq_type> Map_Grp__Size;
typedef map<id_grp_type, Map_Specie__NumRead> Map_Grp__Map_Specie__IDSeq;
typedef map<id_grp_type, Map_Specie__NumRead::value_type> Map_Grp__Max_Pair_Specie__IDSeq;

typedef map<id_cluster_type, id_seq_type> Map_Cluster__Size;
typedef map<id_cluster_type, Map_Specie__NumRead> Map_Cluster__Map_Specie__IDSeq;
typedef map<id_cluster_type, Map_Specie__NumRead::value_type> Map_Cluster__Max_Pair_Specie__IDSeq;

typedef map<id_specie_type, unordered_set<id_seq_type>> Map_Specie__SetIdSeq;
//////////////////////////////////////////////////////////////////////////////////////////////////////

struct Lmer
{
	typedef shared_ptr<Lmer> Ptr;
	string l_mer;
	size_seq count = 1; //1 = non ha in sé il complemento inverso, 2 = ha il complemento inverso
	double prob = 0;
};

struct LMerVectorCompress
{
	typedef shared_ptr<LMerVectorCompress> Ptr;
	typedef map<lMer_type, lMer_type> Map_HashLMer_IndexVector;
	LMerVectorCompress(size_seq L);
	const Lmer::Ptr& GetWithIndex(lMer_type index);
	const Lmer::Ptr& GetWithHash(lMer_type hash);
	lMer_type GetIndexWithHash(lMer_type hash);
	size_seq getL() const;
	const Map_HashLMer_IndexVector& getMapHash() const;
	const vector<Lmer::Ptr>& getLmer() const;

private:
	size_seq L;
	Map_HashLMer_IndexVector mapHash;
	vector<Lmer::Ptr> vLmer;
};

inline hash_type CharToInt(char ch)
{
	if(ch == 'A')
		return 0;
	if(ch == 'C')
		return 1;
	if(ch == 'G')
		return 2;
	if(ch == 'T')
		return 3;
	return 4; //ERROR CODE
}

inline hash_type CharToIntComplement(char ch)
{
	if(ch == 'A')
		return 3;
	if(ch == 'C')
		return 2;
	if(ch == 'G')
		return 1;
	if(ch == 'T')
		return 0;
	return 4; //ERROR CODE
}

inline void GetHash(const string& s_Str, size_seq startQmer, size_seq length,
		HashCorrect& hash, hash_type (*fConvertion)(char))
{
	hash.first = 0;
	hash.second = true;
	#pragma omp parallel for ordered
	for(size_seq i = startQmer; i < startQmer + length; ++i)
	{
		hash_type ch = (*fConvertion)(s_Str[i]);
		#pragma omp ordered
		if(hash.second)
		{
			if(ch == 4) //Errore conversione
				hash.second = false;
			if(hash.second)
				hash.first |= ch << ((i - startQmer) * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
		}
	}
	if(!hash.second)
		hash.first =  0;
}

inline void GetHashes(const string& s_Str, size_seq length, vector<HashCorrect>& vHash, hash_type (*fConvertion)(char))
{
	if(s_Str.size() >= length)
	{
		size_t n_hashes = s_Str.size() - length + 1;
		vHash.resize(n_hashes, HashCorrect(0, true)); //Crea vettore
		vector<size_seq> err(n_hashes, 0); //Vettore che conta errori su qmer
		#pragma omp parallel for ordered
		for(size_t pos=0; pos < s_Str.size(); ++pos)//Computa vettore che mi indica quanti errori ci sono nei qmer, quindi poi posso decidere se calcolare o no.
		{
			bool newErr = (*fConvertion)(s_Str[pos]) == 4;
			#pragma omp ordered
			if(pos < length)
			{
				if(newErr)
					++err[0];
			}
			else
			{
				size_t actual = pos-length+1;
				size_t prev = pos-length;
				bool exitErr = (*fConvertion)(s_Str[prev]) == 4;
				err[actual] = err[prev];
				if(exitErr)
					--err[actual];
				if(newErr)
					++err[actual];
			}
		}
		#pragma omp parallel for
		for(int pos=0; pos < vHash.size(); ++pos)
		{
			if(err[pos] > 0)
			{
				vHash[pos].first = 0;
				vHash[pos].second = false;
			}
		}
		#pragma omp parallel for ordered
		for(int pos=0; pos < vHash.size(); ++pos)
		{
			if(vHash[pos].second) //Se devo computare
			{
				#pragma omp ordered
				if(pos == 0 || !vHash[pos-1].second)
					GetHash(s_Str, pos, length, vHash[pos], fConvertion);
				else
				{
					vHash[pos].first = vHash[pos - 1].first;
					vHash[pos].first -= (*fConvertion)(s_Str[pos - 1]); //sottrai primo elemento che viene eliminato
					vHash[pos].first >>= 2; //dividi per 4, sposta 2 bit
					vHash[pos].first |= ((*fConvertion)(s_Str[pos + length - 1]) << ((length - 1) * 2));	//aggiungi ultimo elemento OR possibile perchè prima ho
																											//diviso per 4 e la posizione dove scrivo ha sicuramente 0
				}
			}
		}
	}
}

void GetKmer(hash_type index, size_seq K, string& Kmer);
vector<string> GetAllKmers(size_seq K);

void createDirAndSubDir(string path);

int parseLineForMemory(char* line);
int getVirtualMemoryUsed();
int getPeakVirtualMemoryUsed();
int getPhysicalMemoryUsed();

#endif
