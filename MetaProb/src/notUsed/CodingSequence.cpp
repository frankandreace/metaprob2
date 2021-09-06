#include "../notUsed/CodingSequence.h"

CodingSequence::CodingSequence() {
	this->lastHashSize = 0;
}

CodingSequence::~CodingSequence() {
	// TODO Auto-generated destructor stub
}

bool CodingSequence::IsOkToBeCoded(string& seq) {
	//Non diverso da ACGT solo queste parole permesse
	bool ret = true;
	bool exact = true;
	for(size_t i = 0; i < seq.size(); ++i)
	{
		switch(seq[i])
		{
		case 'A':
			break;
		case 'C':
			break;
		case 'G':
			break;
		case 'T':
			break;
		default:
			exact = false;
			break;
		}
		if(!exact)
			break;
	}
	ret = seq != "" && exact;
	return ret;
}

bool CodingSequence::Encode(string& seq) {
	if(!this->IsOkToBeCoded(seq))
		return false;
	this->code.clear();
	this->code.resize(seq.length()*2);
	#pragma omp parallel for ordered //Anche in disordine, ma funziona solo con ordered o critical section. Preferisco Ordered
	for(size_t i = 0; i < seq.size(); ++i)
	{
		#pragma omp ordered
		if(seq[i] == 'A')//00 salvato prima bit meno significativo
		{
			this->code.set(i*2, 0);
			this->code.set(i*2+1, 0);
		}
		else if(seq[i] == 'C') //01 salvato prima bit meno significativo
		{
			this->code.set(i*2, 1);
			this->code.set(i*2+1, 0);
		}
		else if(seq[i] == 'G') //10 salvato prima bit meno significativo
		{
			this->code.set(i*2, 0);
			this->code.set(i*2+1, 1);
		}
		else if(seq[i] == 'T') //11 salvato prima bit meno significativo
		{
			this->code.set(i*2, 1);
			this->code.set(i*2+1, 1);
		}
	}
	return true;
}

void CodingSequence::Decode(string& seq)
{
	seq.clear();
	seq.shrink_to_fit();
	seq.resize(this->code.size()*2);
	#pragma omp parallel for
	for(size_t i = 0; i < this->code.size(); ++i)
	{
		if(i%2 == 0) //ogni 2
		{
			switch((char)this->code[i+1])
			{
				case 0:
					switch((char)this->code[i])
					{
					case 0:
						seq[i/2] = 'A';
						break;
					case 1:
						seq[i/2] = 'C';
						break;
					}
					break;
				case 1:
					switch((char)this->code[i])
					{
					case 0:
						seq[i/2] = 'G';
						break;
					case 1:
						seq[i/2] = 'T';
						break;
					}
					break;
			}
		}
	}
}

void CodingSequence::GetAllHashes(VectorHash& hashes, size_t HashSize) {
	if(this->lastHashSize != HashSize)
	{
		this->lastHashSize = HashSize;

		//Dimensiono il vettore in modo da occupare la memoria corretta
		this->hashes.clear();
		this->hashes.shrink_to_fit();
		size_t size_hashes = this->code.size()/2 - this->lastHashSize + 1;
		size_hashes = size_hashes > 0 ? size_hashes : 1;
		this->hashes.resize(size_hashes , 0);

		//Primo hash da computare a parte
		#pragma omp parallel for ordered
		for(size_t i = 0; i < HashSize * 2; ++i)
			if(i < this->code.size())
				#pragma omp ordered
				this->hashes[0] |= ((hash_type)this->code[i].operator bool() << i);

		//Computazione degli altri hash
		hash_type last_hash_compute = this->hashes[0];

		//Se si verifica questa condizione non ci son altri hash da computare
		if(HashSize * 2 >= this->code.size())
			return;

		#pragma omp parallel for ordered
		for(size_t i = (HashSize * 2); i < this->code.size(); ++i)
		{
			#pragma omp ordered
			if(i%2 == 0) //Nuova lettera
			{
				last_hash_compute >>= 2; //dividi per 4 quello appena salvato
				switch((char)this->code[i+1].operator bool())
				{
					case 0:
						switch((char)this->code[i].operator bool())
						{
						case 0: //A=00
							last_hash_compute |= ((hash_type)0 << ((HashSize - 1) * 2));
							break;
						case 1: //C=01
							last_hash_compute |= ((hash_type)1 << ((HashSize - 1) * 2));
							break;
						}
						break;
					case 1:
						switch((char)this->code[i].operator bool())
						{
						case 0: //G=10
							last_hash_compute |= ((hash_type)2 << ((HashSize - 1) * 2));
							break;
						case 1: //T=11
							last_hash_compute |= ((hash_type)3 << ((HashSize - 1) * 2));
							break;
						}
						break;
				}
				this->hashes[i/2 - this->lastHashSize + 1] = last_hash_compute;
			}
		}
	}
	hashes = this->hashes;
}
