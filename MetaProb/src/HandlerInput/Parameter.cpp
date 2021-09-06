/*
 * Parameter.cpp
 *
 *  Created on: 05/feb/2016
 *      Author: samuele
 */

#include "Parameter.h"

Parameter::Parameter(id_cluster_type nClus,
					size_seq q, size_seq m_thres, size_seq_tot SEEDSIZE,
					size_seq l_mer_freq, size_t iteration, size_t time_max,
					Normalization::Norm norm_type, TypeGraph graph, bool estimateK) {
	this->nClus = nClus;
	this->q = q;
	this->m_thres = m_thres;
	this->SEEDSIZE = SEEDSIZE;
	this->l_mer_freq = l_mer_freq;
	this->iteration = iteration;
	this->time_max = time_max;
	this->norm_type = norm_type;
	this->graph = graph;
	this->estimateK = estimateK;
}

Parameter::Parameter(const Parameter::Ptr param)
{
	this->input_files = param->input_files;
	this->nClus = param->nClus;
	this->q = param->q;
	this->m_thres = param->m_thres;
	this->SEEDSIZE = param->SEEDSIZE;
	this->l_mer_freq = param->l_mer_freq;
	this->iteration = param->iteration;
	this->time_max = param->time_max;
	this->norm_type = param->norm_type;
	this->graph = param->graph;
	this->estimateK = param->estimateK;
}

PairFiles& Parameter::insertSingleEndFile(string path) {
	this->input_files.push_back(move(PairFiles(path,"")));
	return this->input_files.back();
}

PairFiles& Parameter::insertPairedEndFiles(string path_end1, string path_end2) {
	this->input_files.push_back(move(PairFiles(path_end1, path_end2)));
	return this->input_files.back();
}

string Parameter::GetGraphTypeString() {
	string graphType = "";
	switch(this->graph)
	{
	case Paired:
		graphType = "Paired";
		break;
	case Single:
		graphType = "Single";
		break;
	case SingleUnion:
		graphType = "Single with union Paired";
		break;
	}
	return graphType;
}

void Parameter::printFiles() {
	for(size_t i = 0; i < this->input_files.size(); ++i)
	{
		if(this->input_files[i].first.isCorrect())
			cout << endl << "Files: " << get<0>(this->input_files[i]).getPath() << flush;
		if(this->input_files[i].second.isCorrect())
			cout << endl << "Files: " << get<1>(this->input_files[i]).getPath() << flush;
	}
}

void Parameter::printParameter()
{
	cout << endl << "N. Cluster: " << this->nClus << flush;
	cout << endl << "Parameter q: " << this->q << flush;
	cout << endl << "Parameter m: " << this->m_thres << flush;
	cout << endl << "Parameter SeedSize: " << this->SEEDSIZE << flush;
	cout << endl << "Parameter lmerfreq: " << this->l_mer_freq << flush;
	cout << endl << "Parameter Kmeans iteration max: " << this->iteration << flush;
	cout << endl << "Norm: " << Normalization::enum_string(this->norm_type) << flush;
	cout << endl << "Graph Type: " << this->GetGraphTypeString() << flush;
}

vector<string> Parameter::GetParametersString() {
	vector<string> parameterStrings;
	for(size_t i = 0; i < this->input_files.size(); ++i)
	{
		if(get<0>(this->input_files[i]).isCorrect())
			parameterStrings.push_back("File: " + get<0>(this->input_files[i]).getPath());
		if(get<1>(this->input_files[i]).isCorrect())
			parameterStrings.push_back("File: " + get<1>(this->input_files[i]).getPath());
	}
	parameterStrings.push_back("Norm: " + Normalization::enum_string(this->norm_type));
	parameterStrings.push_back("N. Cluster: " + to_string(this->nClus));
	parameterStrings.push_back("Parameter q: " + to_string(this->q));
	parameterStrings.push_back("Parameter m: " + to_string(this->m_thres));
	parameterStrings.push_back("Parameter SeedSize: " + to_string(this->SEEDSIZE));
	parameterStrings.push_back("Graph Type: " + this->GetGraphTypeString());
	return parameterStrings;
}
