/*
 * Assessment.cpp
 *
 *  Created on: 08/apr/2015
 *      Author: Samuele Girotto
 */

#include "../notUsed/Assessment.h"

void ResultGroups::ResetMapGroups() {
	this->group_size = Map_Grp__Size();
	this->group_species_size = Map_Grp__Map_Specie__IDSeq();
	this->group_max_species_size = Map_Grp__Max_Pair_Specie__IDSeq();
}

void ResultClusters::ResetMapClusters() {
	this->cluster_size = Map_Cluster__Size();
	this->cluster_species_size = Map_Cluster__Map_Specie__IDSeq();
	this->cluster_max_species_size = Map_Cluster__Max_Pair_Specie__IDSeq();
}

void Result::ResetMap() {
	this->ResetMapGroups();
	this->ResetMapClusters();
	this->map_specie_setID = Map_Specie__SetIdSeq();
}

ResultGroups::ResultGroups() {
	this->assesGroups = make_shared<ValueAssessmentGroup>();
}

ResultClusters::ResultClusters() {
	this->assesClusters = make_shared<ValueAssessmentCluster>();
}

Result::Result(Normalization::Norm norm_type) : ResultGroups(), ResultClusters() {
	this->norm_type = norm_type;
}

void Assessment::Compute(Result::Ptr result, id_seq_type numTotRead)
{
	cout << endl << "Compute Assessment... " << flush;

	//Compute assigned read (read assigned to a cluster)
	result->num_tot_read = numTotRead;
	result->num_assigned_read = 0;
	auto ID = result->cluster_size.begin();
	while(ID != result->cluster_size.end())
	{
		result->num_assigned_read += ID->second;
		++ID;
	}

	//Compute solution precision, recall, f-measure groups
	this->compute_precision_max_achievable(result);
	this->compute_group_of_one_specie(result);
	this->compute_read_in_group_of_one_specie(result);

	//Compute solution precision, recall, f-measure solution
	this->compute_precision_clusters(result);
	this->compute_recall_clusters(result);
	this->compute_F_measure_clusters(result);

	cout << "Complete" << flush;
}

//-------------------------------------------------GROUPS----------------------------------------------------------//

void Assessment::compute_precision_max_achievable(Result::Ptr result)
{
	//Eseguo sommatorie su gruppi e specie massima su quel gruppo
	double max_tot_grp = 0;
	auto grp = result->group_max_species_size.begin();
	while(grp != result->group_max_species_size.end())
	{
		max_tot_grp +=  grp->second.second; //Funziona anche se c'è elemento vuoto pair<"",0>
		++grp;
	}
	result->assesGroups->precision_max_achievable = max_tot_grp/(double)result->num_assigned_read;

//	cout << "Complete" << flush;
}


void Assessment::compute_group_of_one_specie(Result::Ptr result)
{
//	cout << endl << "Compute group of one specie... " << flush;

	//Eseguo sommatorie su gruppi e specie unica su quel gruppo
	double group_of_one_specie = 0;
	auto grp = result->group_species_size.begin();
	while(grp != result->group_species_size.end())
	{
		if(grp->second.size() == 1)
			++group_of_one_specie;
		++grp;
	}
	result->assesGroups->group_of_one_specie = group_of_one_specie/(double)result->group_species_size.size();
}

void Assessment::compute_read_in_group_of_one_specie(Result::Ptr result)
{
	//Eseguo sommatorie su gruppi e specie unica su quel gruppo
	double group_of_one_specie = 0;
	auto grp = result->group_species_size.begin();
	while(grp != result->group_species_size.end())
	{
		if(grp->second.size() == 1)
			group_of_one_specie += grp->second.begin()->second;
		++grp;
	}
	result->assesGroups->read_in_group_of_one_specie = group_of_one_specie/(double)result->num_assigned_read;
}

//-------------------------------------------------CLUSTERS----------------------------------------------------------//

void Assessment::compute_precision_clusters(Result::Ptr result)
{
	//Eseguo sommatorie su cluster e specie massima su quel cluster
	double max_tot_cluster = 0;
	auto cluster = result->cluster_max_species_size.begin();
	while(cluster != result->cluster_max_species_size.end())
	{
		max_tot_cluster +=  cluster->second.second; //Funziona anche se c'è elemento vuoto pair<"",0>
		++cluster;
	}
	result->assesClusters->precision = max_tot_cluster/(double)result->num_assigned_read;
}

void Assessment::compute_recall_clusters(Result::Ptr result)
{
	//Salvo specie e il cluster dove la specie è massima e quanti read contiene di quella specie quel cluster
	map<id_specie_type, pair<id_cluster_type, id_seq_type>> specie_max_cluster;
	auto cluster = result->cluster_species_size.begin();
	while(cluster != result->cluster_species_size.end())
	{
		auto species_cluster_i = cluster->second.begin();
		while(species_cluster_i != cluster->second.end())
		{
			auto species_max = specie_max_cluster.find(species_cluster_i->first);
			if(species_max == specie_max_cluster.end()) //Se non trovato
				specie_max_cluster.insert(make_pair(species_cluster_i->first, make_pair(cluster->first, species_cluster_i->second)));
			else
			{
				if(species_max->second.second < species_cluster_i->second)
					species_max->second = make_pair(cluster->first, species_cluster_i->second);
			}
			++species_cluster_i;
		}
		++cluster;
	}

	//Eseguo sommatorie sui read delle specie nel cluster dove si trovano in maggioranza
	double max_tot_species = 0;
	auto specie_max_cluster_it = specie_max_cluster.begin();
	while(specie_max_cluster_it != specie_max_cluster.end())
	{
		max_tot_species += specie_max_cluster_it->second.second;
		++specie_max_cluster_it;
	}

	result->assesClusters->recall = max_tot_species/(double)result->num_tot_read;
}

void Assessment::compute_F_measure_clusters(Result::Ptr result) {
	result->assesClusters->f_measure = 2 * result->assesClusters->precision * result->assesClusters->recall/(result->assesClusters->precision + result->assesClusters->recall);
}
