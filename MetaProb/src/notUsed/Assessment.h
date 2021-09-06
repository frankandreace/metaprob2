/*
 * Assessment.h
 *
 *  Created on: 08/apr/2015
 *      Author: Samuele Girotto
 */

#ifndef SRC_ASSESSMENT_H_
#define SRC_ASSESSMENT_H_

#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <map>
#include <set>
#include "../Utilities.h"
#include "../Normalization.h"
#include "../HandlerInput/DNAFile.h"

struct ValueAssessmentCluster
{
	typedef shared_ptr<ValueAssessmentCluster> Ptr;
	double precision;
	double recall;
	double f_measure;
};

struct ValueAssessmentGroup
{
	typedef shared_ptr<ValueAssessmentGroup> Ptr;
	double precision_max_achievable;
	double group_of_one_specie;
	double read_in_group_of_one_specie;
};

struct ResultGroups
{
	typedef shared_ptr<ResultGroups> Ptr;

	Map_Grp__Size group_size;
	Map_Grp__Map_Specie__IDSeq group_species_size;
	Map_Grp__Max_Pair_Specie__IDSeq group_max_species_size;
	ValueAssessmentGroup::Ptr assesGroups;
	ResultGroups();
	void ResetMapGroups();
};

struct ResultClusters
{
	typedef shared_ptr<ResultClusters> Ptr;

	Map_Cluster__Size cluster_size;
	Map_Cluster__Map_Specie__IDSeq cluster_species_size;
	Map_Cluster__Max_Pair_Specie__IDSeq cluster_max_species_size;
	ValueAssessmentCluster::Ptr assesClusters;
	ResultClusters();
	void ResetMapClusters();
};

struct Result : public ResultGroups, public ResultClusters
{
	typedef shared_ptr<Result> Ptr;

	Normalization::Norm norm_type; //Normalizzazione computata
	id_seq_type num_tot_read;
	id_seq_type num_assigned_read;
	Map_Specie__SetIdSeq map_specie_setID;
//	clsDNAFile::Map_Specie__IDSpecie_Ptr map_specie_idSpecie;
	Result(Normalization::Norm norm_type);
	void ResetMap();
};

class Assessment {
public:
	typedef shared_ptr<Assessment> Ptr;

	void Compute(Result::Ptr result, id_seq_type numTotRead); //Solution Assessment

private:
	//Result groups
	void compute_precision_max_achievable(Result::Ptr result);
	void compute_group_of_one_specie(Result::Ptr result);
	void compute_read_in_group_of_one_specie(Result::Ptr result);
	//Result clusters
	void compute_precision_clusters(Result::Ptr result);
	void compute_recall_clusters(Result::Ptr result);
	void compute_F_measure_clusters(Result::Ptr result);
};

#endif /* SRC_ASSESSMENT_H_ */
