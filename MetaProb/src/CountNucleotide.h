#ifndef SRC_COUNTNUCLEOTIDE_H_
#define SRC_COUNTNUCLEOTIDE_H_

#include "Utilities.h"

template<typename T>
class CountNucleotide
{
public:
	typedef shared_ptr<CountNucleotide<T>> Ptr;
	T a, c, g, t;
	CountNucleotide();
	CountNucleotide(T a, T c, T g, T t);
	template<typename Z> CountNucleotide(const CountNucleotide<Z> count);
	CountNucleotide<T>& operator +=(const CountNucleotide<T>& b);
	T get_L_mer_count(bool is_PE, size_seq l_mer_freq);
	T get_size();
	double get_probability(char ch);
	double get_reserse_probability(char ch);
	double p_Lmer(string& Lmer);//Compute prop L-mer
};

#include "CountNucleotide.cpp"

#endif /* SRC_COUNTNUCLEOTIDE_H_ */
