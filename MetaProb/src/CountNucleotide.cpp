#ifndef SRC_COUNTNUCLEOTIDE_CPP_
#define SRC_COUNTNUCLEOTIDE_CPP_

#include "CountNucleotide.h"

template<typename T>
CountNucleotide<T>::CountNucleotide()
{
	this->a = 0;
	this->c = 0;
	this->g = 0;
	this->t = 0;
}

template<typename T>
CountNucleotide<T>::CountNucleotide(T a, T c, T g, T t)
{
	this->a = a;
	this->c = c;
	this->g = g;
	this->t = t;
}

template<typename T>
template<typename Z>
CountNucleotide<T>::CountNucleotide(const CountNucleotide<Z> count)
{
	this->a = (T)count.a;
	this->c = (T)count.c;
	this->g = (T)count.g;
	this->t = (T)count.t;
	}

template<typename T>
CountNucleotide<T>& CountNucleotide<T>::operator +=(const CountNucleotide<T>& b)
{
	this->a += b.a;
	this->c += b.c;
	this->g += b.g;
	this->t += b.t;
	return *this;
}

template<typename T>
T CountNucleotide<T>::get_L_mer_count(bool is_PE, size_seq l_mer_freq)
{
	if(is_PE)
		return this->get_size() + 2 * (1 - l_mer_freq);
	else
		return this->get_size() + 1 - l_mer_freq ;
}

template<typename T>
T CountNucleotide<T>::get_size()
{
	return this->a + this->c + this->g + this->t;
}

template<typename T>
double CountNucleotide<T>::get_probability(char ch)
{
	double prob = 0;
	double tot = this->get_size();
	switch(ch)
	{
	case 'A':
		prob = (double)this->a/tot;
		break;
	case 'C':
		prob = (double)this->c/tot;
		break;
	case 'G':
		prob = (double)this->g/tot;
		break;
	case 'T':
		prob = (double)this->t/tot;
		break;
	}
	return prob;
}

template<typename T>
double CountNucleotide<T>::get_reserse_probability(char ch)
{
	double prob = 0;
	double tot = this->get_size();
	switch(ch)
	{
	case 'A':
		prob = (double)this->t/tot;
		break;
	case 'C':
		prob = (double)this->g/tot;
		break;
	case 'G':
		prob = (double)this->c/tot;
		break;
	case 'T':
		prob = (double)this->a/tot;
		break;
	}
	return prob;
}

template<typename T>
double CountNucleotide<T>::p_Lmer(string& Lmer)
{
	double p_Lmer = 1;
	for(lMer_type j = 0; j < Lmer.size(); j++)
		p_Lmer *= this->get_probability(Lmer[j]);

	//Non serve calcolare il reverse, basta invertire solo le basi, il verso di lettura non cambia la probabilità
	double p_rev_Lmer = 1;
	for(lMer_type j = 0; j < Lmer.size(); j++)
		p_rev_Lmer *= this->get_reserse_probability(Lmer[j]);

	return p_Lmer + p_rev_Lmer;
}

#endif /* SRC_COUNTNUCLEOTIDE_CPP_ */
