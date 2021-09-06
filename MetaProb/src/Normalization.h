#ifndef SRC_NORMALIZATION_H_
#define SRC_NORMALIZATION_H_

#include "Utilities.h"
#include "Group.h"
#include "HandlerInput/DNAFile.h"

class Normalization {
public:
	typedef shared_ptr<Normalization> Ptr;
	enum Norm {
		NORM_D2star_All_Read_Prob_Lmer_Euclidian = 1, //NORM_D2star_All_Read_Prob_Kmer + norma euclidea
		NORM_D2star_Group_Prob_Bernulli_Lmer_Euclidian, //NORM_D2star_Group_Prob_Bernulli_Kmer + norma euclidea
		NORM_BiMeta, //Norma del paper
		NORM_D2star_Group_Prob_Lmer, //Probabilità del K-mer dato il vettore del conteggio presenza K-mer in read seed
		NORM_D2star_All_Read_Prob_Lmer,
		NORM_D2star_Group_Prob_Bernulli_Lmer, //Probabilità del K-mer tenendo conto di quante volte compaiono i nucleotidi in un gruppo
		NORM_D2star_All_Read_Prob_Bernoulli, //Probabilità del K-mer tenendo conto di quante volte compaiono i nucleotidi in tutti i read e calcolo probabilità K-mer con modello Bernulli
		NORM_D2star_All_Seed_Read_Prob_Bernoulli, //Probabilità del K-mer tenendo conto di quante volte compaiono i nucleotidi nei read seed e calcolo probabilità K-mer con modello Bernulli
		NORM_BiMeta_Euclidian, //NORM_Size_Seed + norma euclidea
		NORM_D2star_Group_Prob_Lmer_Euclidian, //NORM_D2star_Group_Prob_Kmer + norma euclidea
		NORM_D2star_All_Read_Prob_Bernoulli_Euclidian, //NORM_D2star_Prob_Bernoulli_All_Read + norma euclidea
		NORM_D2star_All_Seed_Read_Prob_Bernoulli_Euclidian, //NORM_D2star_Prob_Bernoulli_All_Seed_Read + norma euclidea
		Last //Isn't a norm type, is only a indicator
	}; //Norm Type

	static string enum_string(unsigned type);
	static Norm int_to_enum(unsigned type);

	Normalization(clsDNAFile::Ptr file, vector<clsGroup::Ptr>& vGroup, Norm norm_type, size_seq l_mer_size);
	void Compute();

	const LMerVectorCompress::Ptr& getLmerCompress() const;

private:
	Norm norm_type;
	size_seq l_mer_size;
	LMerVectorCompress::Ptr compr_vLmer; //Contains Vector of l-mers
	vector<clsGroup::Ptr> vGroup; //Gruppi da computare

	void CountLmerGroups(clsDNAFile::Ptr file);//Compute L-mer in groups
	void ComputeFreq(clsGroup::Ptr grp);//Compute and Normalize vFreq

	void GenerateLmVector();//Generate l-mer vector
};

#endif /* SRC_NORMALIZATION_H_ */
