/*
 * Gmeans.cpp
 *
 *  Created on: 29/feb/2016
 *      Author: samuele
 */

#include "Gmeans.h"

Gmeans::Gmeans(MatrixXld& observation):kmeans(make_shared<Kmeans>(observation)), observation(observation) {
	this->estimateK = 0;
	this->data_to_cluster.resize(this->observation.cols());
}

Gmeans::~Gmeans() {
	// TODO Auto-generated destructor stub
}

Gmeans::Cluster Gmeans::getEstimateK() const {
	return estimateK;
}

Gmeans::VectorCluster& Gmeans::compute(	Cluster start_k, uint64_t iter_split_max,
										uint64_t iter_kmeans, uint64_t time_in_second, TypeVector eps_ch_thres, TypeVector p_value)
{
	//Genera centroidi e computa prima iterazione kmeans
	Kmeans::Vector_VectorXld centroids;
	Kmeans::GenerateFirstRandomCentroidsFromData(this->observation, centroids, start_k);
	Kmeans::VectorIndex vector_index(this->observation.cols()); //Indici rimanenti, quando computo kmeans mapClusterToIndex
																//tiene un riferimento a questo vettore, cioè se
																//mapClusterToIndex[cls][i] da un indice questo è rispetto
																//a questo vettore
	#pragma omp parallel for //Inizializza con i rimantenti tutti i punti di observation
	for(size_t i = 0; i < (size_t)this->observation.cols(); ++i)
		vector_index[i] = i;

	Kmeans::MatrixXld observation_iter; //observation di vector_index
	Kmeans::mapClusterIndex mapClusterToIndex; //clusterizzazione di vector_index
	uint64_t iter_split = 0;
	bool there_is_splitted = true;
	while(there_is_splitted && iter_split < iter_split_max)
	{
		ClusterUtility::GetSubMatrix(this->observation, vector_index, observation_iter);
		kmeans = make_shared<Kmeans>(observation_iter);
		kmeans->compute(iter_kmeans, time_in_second, eps_ch_thres, centroids);
		kmeans->GetSplitIndex(mapClusterToIndex);

		//Splitta se possibile
		vector<ClusterSplitter::Result> v_result_split(centroids.size());
		#pragma omp parallel for
		for(size_t cls = 0; cls < centroids.size(); ++cls)
		{
			ClusterSplitter cls_stack(observation_iter, mapClusterToIndex[cls]);
			cls_stack.isClusterSplitted(iter_kmeans, time_in_second, eps_ch_thres, p_value, v_result_split[cls]);
		}

		there_is_splitted = false;
		bool any_no_splitted = false;
		for(size_t i = 0; i < v_result_split.size(); i++)
		{
			there_is_splitted |= v_result_split[i].isSplitted;
			any_no_splitted |= !v_result_split[i].isSplitted;
		}
//		auto comp = [](ClusterSplitter::Result r1, ClusterSplitter::Result r2){return r1.p_value < r2.p_value;};
//		auto it_min = min_element(v_result_split.begin(), v_result_split.end(), comp);
//		auto it_max = max_element(v_result_split.begin(), v_result_split.end(), comp);
//		size_t pos_min = distance(v_result_split.begin(), it_min);
//		size_t pos_max = distance(v_result_split.begin(), it_max);
//		there_is_splitted &= v_result_split[pos_min].isSplitted;
//		bool any_no_splitted = !v_result_split[pos_max].isSplitted;

		size_t splitted_max = 1; //Max is centroids.size()
		size_t splitted = 0;
		Kmeans::Vector_VectorXld new_centroids;
		Kmeans::VectorIndex new_vector_index;
		for(size_t cls = 0; cls < centroids.size(); ++cls)
		{
			if(v_result_split[cls].isSplitted)
			{
				if(splitted < splitted_max && !any_no_splitted)
				{
					++splitted;
					for(size_t i = 0; i < v_result_split[cls].centroidSplit.size(); ++i)
						new_centroids.push_back(v_result_split[cls].centroidSplit[i]);
				}
				else
					new_centroids.push_back(centroids[cls]);
				for(size_t i = 0; i < mapClusterToIndex[cls].size(); ++i)
					new_vector_index.push_back(vector_index[mapClusterToIndex[cls][i]]);
			}
			else//Non si verifica mai ma lo metto perché basta togliere if(there_is_no_splitted) per splittare e inserire nuovi cluster nella stessa iterazione
			{
				this->centroid_saved.push_back(centroids[cls]);
				#pragma omp parallel for
				for(size_t i = 0; i < mapClusterToIndex[cls].size(); ++i)
					this->data_to_cluster[vector_index[mapClusterToIndex[cls][i]]] = this->centroid_saved.size()-1;
			}
		}
		vector_index = new_vector_index;
		centroids = new_centroids;

		++iter_split;
		cout << "\r" << "Estimate K and clustering... Cluster found: " << this->centroid_saved.size() << flush;
	}
	this->estimateK = this->centroid_saved.size();
//	this->data_to_cluster = kmeans.getDataToCluster();
	return this->data_to_cluster;
}

ClusterSplitter::ClusterSplitter(MatrixXld& observation, VectorIndex& indexCluster)
{
	ClusterUtility::GetSubMatrix(observation, indexCluster, this->observation);
	ClusterUtility::Mean(this->observation, this->centroid);
}

//Accepted split if distance.norm() > this->AD_thres
void ClusterSplitter::isClusterSplitted(uint64_t iter, uint64_t second_in_time, TypeVector eps_ch_thres,
		TypeVector p_value_thresh, Result& result) {

	if(this->observation.cols() < 25)
	{
		result.isSplitted = false;
		return;
	}

	this->GetSplitCentroid(result.centroidSplit);
	Kmeans kmeans(this->observation);
	kmeans.compute(iter, second_in_time, eps_ch_thres, result.centroidSplit);

	mapClusterIndex mapIdx;
	kmeans.GetSplitIndex(mapIdx);
	bool separateIdx = mapIdx.size() > 1;

	VectorXld distance = result.centroidSplit[0] - result.centroidSplit[1];
	MatrixXld proiection = distance.transpose() * this->observation/distance.squaredNorm();

	//Compute mean and covariance (1x1 only variance)
	VectorXld mean;
	MatrixXld covariance;
	ClusterUtility::Mean(proiection, mean);
	ClusterUtility::CovarianceMatrix(proiection, covariance);

	//standardize
	VectorXld proiection_ev = proiection.row(0);
	proiection_ev -= VectorXld::Constant(proiection_ev.size(), mean(0));
	TypeVector std_deviation = sqrt((TypeVector)covariance(0,0));
	proiection_ev /= std_deviation;

	//Sort
	vector<TypeVector> v_proiection(proiection_ev.data(), proiection_ev.data() + proiection_ev.rows() * proiection_ev.cols());
	sort(v_proiection.begin(), v_proiection.end());

//	//KSONE
//	TypeVector d;
//	ksone(v_proiection, normalCDF, d, result.p_value);

	//KSTWO
	//Empirical_cdf
	vector<TypeVector> v_proiection_Empirical_cdf(v_proiection.size());
	#pragma omp parallel for
	for(size_t i = 0; i < v_proiection_Empirical_cdf.size(); ++i)
		v_proiection_Empirical_cdf[i] = (TypeVector)i/(TypeVector)v_proiection.size();

	//Real_cdf
	boost::math::normal_distribution<TypeVector> normal;
	vector<TypeVector> v_proiection_Real_cdf(v_proiection.size());
	#pragma omp parallel for
	for(size_t i = 0; i < v_proiection_Real_cdf.size(); ++i)
		v_proiection_Real_cdf[i] =  boost::math::cdf(normal, v_proiection[i]);

	TypeVector d;
	kstwo(v_proiection_Empirical_cdf, v_proiection_Real_cdf, d, result.p_value);

	result.isSplitted = result.p_value < p_value_thresh && separateIdx;
}

////Accepted split if distance.norm() > this->AD_thres
//bool ClusterStack::isClusterSplitted(uint64_t iter, uint64_t second_in_time, TypeVector p_value,
//		Vector_VectorXld& centroidSplit) {
//
//	if(this->observation.cols() < 8)
//		return false;
//
//	this->GetSplitCentroid(centroidSplit);
//	Kmeans kmeans(this->observation);
//	kmeans.compute(iter, second_in_time, centroidSplit);
//
//	VectorXld distance = centroidSplit[0] - centroidSplit[1];
//	MatrixXld proiection = distance.transpose() * this->observation/distance.squaredNorm();
//
//	//Compute mean and covariance (1x1 only variance)
//	VectorXld mean;
//	MatrixXld covariance;
//	ClusterUtility::Mean(proiection, mean);
//	ClusterUtility::CovarianceMatrix(proiection, covariance);
//
//	//Normalize
//	VectorXld proiection_ev = proiection.row(0);
//	proiection_ev -= VectorXld::Constant(proiection_ev.size(), mean(0));
//	TypeVector std_deviation = sqrt((TypeVector)covariance(0,0));
//	proiection_ev /= std_deviation;
//
//	//Sort
//	vector<TypeVector> v_proiection(proiection_ev.data(), proiection_ev.data() + proiection_ev.rows() * proiection_ev.cols());
//	sort(v_proiection.begin(), v_proiection.end());
//	Data data(v_proiection.data(), v_proiection.size());
//	proiection_ev = data;
//
//	boost::math::normal_distribution<TypeVector> normal;
//	TypeVector AD_2 = 0;
//	for(size_t i = 1; i <= (size_t)proiection_ev.size(); ++i)
//	{
//		TypeVector i_1 = proiection_ev(i-1);
//		TypeVector n_i = proiection_ev(proiection_ev.size() - i);
//		TypeVector cdf_i_1 = boost::math::cdf(normal, i_1);
//		TypeVector cdf_n_i = boost::math::cdf(normal, n_i);
//		AD_2 += (2*i - 1)*(log(cdf_i_1) + log(1 - cdf_n_i));
//	}
//	AD_2 /= -(proiection_ev.size());
//	AD_2 -= proiection_ev.size();
//	TypeVector AD_2_star = AD_2 * (1 + 0.75/proiection_ev.size() + 2.25/pow(proiection_ev.size(), 2));
//	this->p_value = this->convertAD2_star_To_PValue(AD_2_star);
//	return this->p_value < p_value;
//}

//ClusterSplitter::TypeVector ClusterSplitter::convertAD2_star_To_PValue(TypeVector AD_2_star)
//{
////	If AD*=>0.6, then p = exp(1.2937 - 5.709(AD*)+ 0.0186(AD*)2
////	If 0.34 < AD* < .6, then p = exp(0.9177 - 4.279(AD*) - 1.38(AD*)2
////	If 0.2 < AD* < 0.34, then p = 1 - exp(-8.318 + 42.796(AD*)- 59.938(AD*)2)
////	If AD* <= 0.2, then p = 1 - exp(-13.436 + 101.14(AD*)- 223.73(AD*)2)
//	if(AD_2_star >= 0.6)
//		return exp(1.2937 - 5.709*AD_2_star + 0.0186*pow(AD_2_star,2));
//	else if(0.34 <= AD_2_star && AD_2_star  < 0.6)
//		return exp(0.9177 - 4.279*AD_2_star - 1.38*pow(AD_2_star,2));
//	else if(0.2 < AD_2_star && AD_2_star < 0.34)
//		return 1 - exp(-8.318 + 42.796*AD_2_star - 59.938*pow(AD_2_star,2));
//	else //if(AD_2_star <= 0.2)
//		return 1 - exp(-13.436 + 101.14*AD_2_star- 223.73*pow(AD_2_star,2));
//}

void ClusterSplitter::GetSplitCentroid(
		Vector_VectorXld& centroidSplit) {
	MatrixXld covariance;
	ClusterUtility::CovarianceMatrix(this->observation, covariance);

	//Trova autovettore e autovalori
	SelfAdjointEigenSolver_ eigensolver(covariance);
	if (eigensolver.info() != Success)
	{
		cout << endl << "Error to compute eigenvectors and eigenvalues." << flush;
		abort();
	}
	//Find max eigenvalue
	VectorXld eigenvalues = eigensolver.eigenvalues();
	size_t indexMax = 0;
	TypeVector eigenvaluesMax = eigenvalues.maxCoeff(&indexMax);
	size_t indexMin = 0;
	TypeVector eigenvaluesMin = eigenvalues.minCoeff(&indexMin);

	size_t indexToUse = 0;
	TypeVector eigenvaluesToUse = 0;
	if(abs(eigenvaluesMin) < abs(eigenvaluesMax))
	{
		indexToUse = indexMax;
		eigenvaluesToUse = eigenvaluesMax;
	}
	else
	{
		indexToUse = indexMin;
		eigenvaluesToUse = eigenvaluesMin;
	}

	VectorXld m = eigensolver.eigenvectors().col(indexToUse) * sqrt(2*eigenvaluesToUse/M_PI);
	centroidSplit.push_back(this->centroid + m);
	centroidSplit.push_back(this->centroid - m);
}
