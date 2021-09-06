/*
 * Parameter.h
 *
 *  Created on: 05/feb/2016
 *      Author: samuele
 */

#ifndef SRC_PARAMETER_H_
#define SRC_PARAMETER_H_

#include "Input.h"
#include "../Utilities.h"
#include "../Normalization.h"

class Parameter
{
public:
	typedef shared_ptr<Parameter> Ptr;

	vector<PairFiles> input_files;
	id_cluster_type nClus;//number of clusters
	size_seq q;//length of q-mers of overlapping
	size_seq m_thres;//Number of shared l-mers between reads
	size_seq_tot SEEDSIZE;//The maximum number of reads in each seed
	size_seq l_mer_freq; //size L-mer for count frequency
	size_t iteration; //Kmeans iteration
	size_t time_max; //Kmeans time max
	Normalization::Norm norm_type;
	TypeGraph graph;
	bool estimateK;

	Parameter():input_files(), nClus(0),
				q(30), m_thres(5), SEEDSIZE(20), l_mer_freq(4),
				iteration(100), time_max(3600),
				norm_type(Normalization::NORM_D2star_All_Read_Prob_Lmer_Euclidian), graph(Paired), estimateK(true){}

	Parameter(	id_cluster_type nClus,
				size_seq q, size_seq m_thres, size_seq_tot SEEDSIZE, size_seq l_mer_freq,
				size_t iteration, size_t time_max,
				Normalization::Norm norm_type, TypeGraph graph, bool estimateK);

	Parameter(const Parameter::Ptr param);

	PairFiles& insertSingleEndFile(string path); //return last insert
	PairFiles& insertPairedEndFiles(string path_end1, string path_end2); //return last insert

	void printFiles();
	void printParameter();
	vector<string> GetParametersString();
	string GetGraphTypeString();
};

#endif /* SRC_PARAMETER_H_ */
