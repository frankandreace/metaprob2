#include "Normalization.h"

const LMerVectorCompress::Ptr& Normalization::getLmerCompress() const {
	return compr_vLmer;
}

//Generate l-mer vector
//TODO: creare classe che fa il vettore e salva anche
//le posizioni nel vettore in base ad un hash
void Normalization::GenerateLmVector()
{
	this->compr_vLmer = make_shared<LMerVectorCompress>(this->l_mer_size);
}

Normalization::Normalization(clsDNAFile::Ptr file, vector<clsGroup::Ptr>& vGroup, Norm norm_type, size_seq l_mer_size)
{
	this->vGroup = vGroup; //Copy constructor of the vector
	this->norm_type = norm_type;
	this->l_mer_size = l_mer_size;
	this->GenerateLmVector();

	//Count l-mer in groups from file
	this->CountLmerGroups(file);

	//Compute probability
	if(this->norm_type == NORM_D2star_All_Read_Prob_Lmer || this->norm_type == NORM_D2star_All_Read_Prob_Lmer_Euclidian)
	{
		clsGroup::VectorFeature vFreq;
		for(id_grp_type i = 0; i <this->vGroup.size(); i++)
		{
			if(i == 0)
				vFreq = this->vGroup[i]->vFreq;
			else
			{
				for(lMer_type j = 0; j < this->vGroup[i]->vFreq.size(); j++)
					vFreq[j] += this->vGroup[i]->vFreq[j];
			}
		}
		long double tot = 0;
		for(lMer_type i = 0; i < vFreq.size(); i++)
			tot += vFreq[i];
		for(lMer_type i = 0; i < vFreq.size(); i++)
			this->compr_vLmer->GetWithIndex(i)->prob = (double)((long double)vFreq[i]/tot);
	}
	CountNucleotide<size_seq_tot> count;
	if(this->norm_type == NORM_D2star_All_Read_Prob_Bernoulli || this->norm_type == NORM_D2star_All_Read_Prob_Bernoulli_Euclidian)
	{
		for(id_grp_type i = 0; i < this->vGroup.size(); i++)
		{
			for(id_seq_type j = 0; j < this->vGroup[i]->vSeedSeqPE.size(); j++)
			{
				const CountNucleotide<size_seq>& info = this->vGroup[i]->vSeedSeqPE[j]->getInfo();
				count += CountNucleotide<size_seq_tot>(info);
			}
			for(id_seq_type j = 0; j < this->vGroup[i]->vSeqPE.size(); j++)
			{
				const CountNucleotide<size_seq>& info = this->vGroup[i]->vSeqPE[j]->getInfo();
				count += CountNucleotide<size_seq_tot>(info);
			}
		}
	}
	if(this->norm_type == NORM_D2star_All_Seed_Read_Prob_Bernoulli || this->norm_type == NORM_D2star_All_Seed_Read_Prob_Bernoulli_Euclidian)
	{
		for(id_grp_type i = 0; i < this->vGroup.size(); i++)
			for(id_seq_type j = 0; j < this->vGroup[i]->vSeedSeqPE.size(); j++)
			{
				const CountNucleotide<size_seq>& info = this->vGroup[i]->vSeedSeqPE[j]->getInfo();
				count += CountNucleotide<size_seq_tot>(info);
			}
	}
	if(this->norm_type == NORM_D2star_All_Read_Prob_Bernoulli || this->norm_type == NORM_D2star_All_Seed_Read_Prob_Bernoulli
		|| this->norm_type == NORM_D2star_All_Read_Prob_Bernoulli_Euclidian || this->norm_type == NORM_D2star_All_Seed_Read_Prob_Bernoulli_Euclidian)
	{
		for(lMer_type i = 0; i < this->compr_vLmer->getLmer().size(); i++)
			this->compr_vLmer->GetWithIndex(i)->prob = count.p_Lmer(this->compr_vLmer->GetWithIndex(i)->l_mer);
	}
}

void Normalization::Compute()
{
	cout << endl << flush;
	cout << "Compute feature of group... " << flush;

	//Compute freq for all group
	for(id_grp_type i=0;i<this->vGroup.size();i++)
	{
//		//Stampa
//		cout << "\r" << "Compute feature of group: " << i + 1 << flush;
		this->ComputeFreq(this->vGroup[i]);
	}

	cout << " Complete" << flush;
}

void Normalization::ComputeFreq(clsGroup::Ptr grp)
{
	//l-mer is already counted in grp->freq (and memorized in grp->vCountMer), compute freq with norm type
	//Compute different prob
	if(this->norm_type == NORM_D2star_Group_Prob_Lmer || this->norm_type == NORM_D2star_Group_Prob_Lmer_Euclidian)
	{
		double tot = 0;
		for(lMer_type i = 0; i < this->compr_vLmer->getLmer().size(); i++)
			tot += grp->vFreq[i];

		if(tot == 0)
			tot = 1; //E' 0 se son tutti zero, così ho il vettore a zero

		for(lMer_type i = 0; i < this->compr_vLmer->getLmer().size(); i++)
			this->compr_vLmer->GetWithIndex(i)->prob = grp->vFreq[i]/tot;
	}
	if(this->norm_type == NORM_D2star_Group_Prob_Bernulli_Lmer || this->norm_type == NORM_D2star_Group_Prob_Bernulli_Lmer_Euclidian)
	{
		CountNucleotide<size_seq_tot> g_count;
		for(id_seq_type i = 0; i < grp->vSeedSeqPE.size(); i++)
		{
			const CountNucleotide<size_seq>& info = grp->vSeedSeqPE[i]->getInfo();
			g_count += CountNucleotide<size_seq_tot>(info);
		}
		for(id_seq_type i = 0; i < grp->vSeqPE.size(); i++)
		{
			const CountNucleotide<size_seq>& info = grp->vSeqPE[i]->getInfo();
			g_count += CountNucleotide<size_seq_tot>(info);
		}
		for(lMer_type i = 0; i < this->compr_vLmer->getLmer().size(); i++)
			this->compr_vLmer->GetWithIndex(i)->prob = g_count.p_Lmer(this->compr_vLmer->GetWithIndex(i)->l_mer);
	}

	//Compute the standardization
	if(this->norm_type == NORM_BiMeta || this->norm_type == NORM_BiMeta_Euclidian)
	{
		for(lMer_type i = 0; i < this->compr_vLmer->getLmer().size(); i++)
			grp->vFreq[i] = grp->vFreq[i]/(double)grp->get_Seeds_L_mer_count(this->l_mer_size);
	}
	if(	this->norm_type == NORM_D2star_Group_Prob_Lmer || this->norm_type == NORM_D2star_Group_Prob_Lmer_Euclidian
		|| this->norm_type == NORM_D2star_Group_Prob_Bernulli_Lmer || this->norm_type == NORM_D2star_Group_Prob_Bernulli_Lmer_Euclidian
		|| this->norm_type == NORM_D2star_All_Read_Prob_Lmer || this->norm_type == NORM_D2star_All_Read_Prob_Lmer_Euclidian
		|| this->norm_type == NORM_D2star_All_Read_Prob_Bernoulli || this->norm_type == NORM_D2star_All_Read_Prob_Bernoulli_Euclidian
		|| this->norm_type == NORM_D2star_All_Seed_Read_Prob_Bernoulli || this->norm_type == NORM_D2star_All_Seed_Read_Prob_Bernoulli_Euclidian)
	{
		for(lMer_type i = 0; i < this->compr_vLmer->getLmer().size(); i++)
			if(this->compr_vLmer->GetWithIndex(i)->prob != 0)
				grp->vFreq[i] = (grp->vFreq[i] - (double)grp->get_Seeds_L_mer_count(this->l_mer_size)*this->compr_vLmer->GetWithIndex(i)->prob)/sqrt((double)grp->get_Seeds_L_mer_count(this->l_mer_size)*(1-this->compr_vLmer->GetWithIndex(i)->prob)*this->compr_vLmer->GetWithIndex(i)->prob);
	}

	//Compute the normalization
	if(	this->norm_type == NORM_BiMeta_Euclidian
		|| this->norm_type == NORM_D2star_Group_Prob_Lmer_Euclidian
		|| this->norm_type == NORM_D2star_Group_Prob_Bernulli_Lmer_Euclidian
		|| this->norm_type == NORM_D2star_All_Read_Prob_Lmer_Euclidian
		|| this->norm_type == NORM_D2star_All_Read_Prob_Bernoulli_Euclidian
		|| this->norm_type == NORM_D2star_All_Seed_Read_Prob_Bernoulli_Euclidian)
	{
		//Passaggio per evitare overflow
		double y = 0;
		for(lMer_type i = 0; i < this->compr_vLmer->getLmer().size(); i++)
			if(y < abs(grp->vFreq[i]))
				y = abs(grp->vFreq[i]);
		//Calcolo con divisione per y
		double tot = 0;
		for(lMer_type i = 0; i < this->compr_vLmer->getLmer().size(); i++)
			tot += this->compr_vLmer->GetWithIndex(i)->count * pow(grp->vFreq[i]/(y*this->compr_vLmer->GetWithIndex(i)->count), 2); //TODO:count perché si tiene conto anche del reverse
		//Rimoltiplico y
		tot = y*sqrt(tot);

		if(tot == 0)
			tot = 1; //E' 0 se son tutti zero, così ho il vettore finale a zero senza errori

		for(lMer_type i = 0; i  < this->compr_vLmer->getLmer().size(); i++)
			grp->vFreq[i] = grp->vFreq[i]/tot;
	}
}

//Count l-mer in groups
void Normalization::CountLmerGroups(clsDNAFile::Ptr file)
{
	//Verifica se già contati
	//(contati se il vettore conteggio è almeno grande quanto il vettore vLmer per ogni gruppo)
	bool count_exe = true;
	for(id_grp_type i = 0; i < this->vGroup.size(); i++)
		count_exe = count_exe && (this->vGroup[i]->vCountMer.size() == this->compr_vLmer->getLmer().size());

	if(!count_exe) //Non contati
	{
		for(id_grp_type i = 0; i < this->vGroup.size(); i++)
		{
			//Azzera valori anche se già a zero
			this->vGroup[i]->vCountMer.clear();
			this->vGroup[i]->vFreq.clear();

			//Inizializza conteggio
			this->vGroup[i]->vCountMer.resize(this->compr_vLmer->getLmer().size(), 0);
		}
		file->CountSeedLmerFromFile(this->vGroup, this->compr_vLmer);
	}

	//Salva in vFreq il contatore, così il conteggio non viene toccato e può essere riutilizzato
	for(id_grp_type i = 0; i < this->vGroup.size(); i++)
		this->vGroup[i]->vFreq = this->vGroup[i]->vCountMer; //Copy costructor
}

string Normalization::enum_string(unsigned type)
{
	string returnState = "";
	switch(type)
	{
		case (NORM_BiMeta):
			returnState = "BiMeta";
		break;
		case (NORM_D2star_Group_Prob_Lmer):
			returnState = "D2star_Group_Prob_Lmer";
		break;
		case (NORM_D2star_Group_Prob_Bernulli_Lmer):
			returnState = "D2star_Group_Prob_Bernulli_Lmer";
		break;
		case (NORM_D2star_All_Read_Prob_Lmer):
			returnState = "D2star_All_Read_Prob_Lmer";
		break;
		case (NORM_D2star_All_Read_Prob_Bernoulli):
			returnState = "D2star_All_Read_Prob_Bernoulli";
		break;
		case (NORM_D2star_All_Seed_Read_Prob_Bernoulli):
			returnState = "D2star_All_Seed_Read_Prob_Bernoulli";
		break;
		case (NORM_BiMeta_Euclidian):
			returnState = "BiMeta_Euclidian";
		break;
		case(NORM_D2star_Group_Prob_Lmer_Euclidian):
			returnState = "D2star_Group_Prob_Kmer_Euclidian";
		break;
		case (NORM_D2star_Group_Prob_Bernulli_Lmer_Euclidian):
			returnState = "D2star_Group_Prob_Bernulli_Lmer_Euclidian";
		break;
		case (NORM_D2star_All_Read_Prob_Lmer_Euclidian):
			returnState = "D2star_All_Read_Prob_Lmer_Euclidian";
		break;
		case(NORM_D2star_All_Read_Prob_Bernoulli_Euclidian):
			returnState = "D2star_All_Read_Prob_Bernoulli_Euclidian";
		break;
		case(NORM_D2star_All_Seed_Read_Prob_Bernoulli_Euclidian):
			returnState = "D2star_All_Seed_Read_Prob_Bernoulli_Euclidian";
		break;
	}
	return returnState;
}

Normalization::Norm Normalization::int_to_enum(unsigned type)
	{
		Norm returnState = Last;
		switch(type)
		{
			case (NORM_BiMeta):
				returnState = NORM_BiMeta;
			break;
			case (NORM_D2star_Group_Prob_Lmer):
				returnState = NORM_D2star_Group_Prob_Lmer;
			break;
			case (NORM_D2star_Group_Prob_Bernulli_Lmer):
				returnState = NORM_D2star_Group_Prob_Bernulli_Lmer;
			break;
			case (NORM_D2star_All_Read_Prob_Lmer):
				returnState = NORM_D2star_All_Read_Prob_Lmer;
			break;
			case (NORM_D2star_All_Read_Prob_Bernoulli):
				returnState = NORM_D2star_All_Read_Prob_Bernoulli;
			break;
			case (NORM_D2star_All_Seed_Read_Prob_Bernoulli):
				returnState = NORM_D2star_All_Seed_Read_Prob_Bernoulli;
			break;
			case (NORM_BiMeta_Euclidian):
				returnState = NORM_BiMeta_Euclidian;
			break;
			case (NORM_D2star_Group_Prob_Lmer_Euclidian):
				returnState = NORM_D2star_Group_Prob_Lmer_Euclidian;
			break;
			case (NORM_D2star_Group_Prob_Bernulli_Lmer_Euclidian):
				returnState = NORM_D2star_Group_Prob_Bernulli_Lmer_Euclidian;
			break;
			case (NORM_D2star_All_Read_Prob_Lmer_Euclidian):
				returnState = NORM_D2star_All_Read_Prob_Lmer_Euclidian;
			break;
			case (NORM_D2star_All_Read_Prob_Bernoulli_Euclidian):
				returnState = NORM_D2star_All_Read_Prob_Bernoulli_Euclidian;
			break;
			case (NORM_D2star_All_Seed_Read_Prob_Bernoulli_Euclidian):
				returnState = NORM_D2star_All_Seed_Read_Prob_Bernoulli_Euclidian;
			break;
		}
		return returnState;
	}

