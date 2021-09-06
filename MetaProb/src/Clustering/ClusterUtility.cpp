/*
 * ClusterUtility.cpp
 *
 *  Created on: 29/feb/2016
 *      Author: samuele
 */

#include "ClusterUtility.h"

ClusterUtility::ClusterUtility() {
	this->dim_feature = 0;
}

ClusterUtility::~ClusterUtility() {
	// TODO Auto-generated destructor stub
}

void ClusterUtility::addData(vector<long double>* data) {
	this->v_data.push_back(data);
	if(this->dim_feature < (*data).size())
		this->dim_feature = (*data).size();
}

void ClusterUtility::Convert(MatrixXld& observation) {
	observation.setZero(this->dim_feature, this->v_data.size());
	#pragma omp parallel for
	for(size_t col = 0; col < this->v_data.size(); ++col)
		observation.col(col) = Data(this->v_data[col]->data(), this->v_data[col]->size()) ;
}

void ClusterUtility::Convert(Vector_VectorXld& observation)
{
	observation.resize(this->v_data.size());
	#pragma omp parallel for
	for(size_t col = 0; col < this->v_data.size(); ++col)
		observation[col] = Data(this->v_data[col]->data(), this->v_data[col]->size()) ;
}

void ClusterUtility::GetSubMatrix(MatrixXld& observation, vector<size_t>& index,
		MatrixXld& sub_observation) {
	size_t num = index.size() <= (size_t)observation.cols() ? index.size() : observation.cols();
	sub_observation.setZero(observation.rows(), num);
	#pragma omp parallel for
	for(size_t i = 0; i < num; ++i)
		sub_observation.col(i) = observation.col(index[i]);
}

void ClusterUtility::GetSubMatrix(Vector_VectorXld& observation, vector<size_t>& index,
		MatrixXld& sub_observation) {
	size_t feature = 0;
	for(size_t i = 0; i < observation.size(); ++i)
		if((size_t)observation[i].size() > feature)
			feature = observation[i].size();
	size_t num = index.size() <= (size_t)observation.size() ? index.size() : observation.size();

	sub_observation.setZero(feature, num);
	#pragma omp parallel for
	for(size_t i = 0; i < num; ++i)
		sub_observation.col(i) = observation[index[i]];
}

void ClusterUtility::GetSubMatrix(MatrixXld& observation, vector<size_t>& index, Vector_VectorXld& sub_observation)
{
	size_t num = index.size() <= (size_t)observation.cols() ? index.size() : observation.cols();
	sub_observation.resize(num);
	for(size_t i = 0; i < num; ++i)
		sub_observation[i] = observation.col(index[i]);
}

void ClusterUtility::GetSubMatrix(Vector_VectorXld& observation, vector<size_t>& index, Vector_VectorXld& sub_observation)
{
	sub_observation.resize(index.size());
	for(size_t i = 0; i < index.size(); ++i)
		sub_observation[i] = observation[index[i]];
}

void ClusterUtility::Convert(MatrixXld& from_observation, Vector_VectorXld& to_observation)
{
	if(from_observation.size() == 0)
		return;

	to_observation.resize(from_observation.cols());
	#pragma omp parallel for
	for(size_t i = 0; i < (size_t)from_observation.cols(); ++i)
		to_observation[i] = from_observation.col(i);
}

void ClusterUtility::Convert(Vector_VectorXld& from_observation, MatrixXld& to_observation)
{
	if(from_observation.size() == 0)
		return;
	size_t feature = from_observation[0].size();

	to_observation.setZero(feature, from_observation.size());
	#pragma omp parallel for
	for(size_t i = 0; i < from_observation.size(); ++i)
		to_observation.col(i) = from_observation[i];
}

void ClusterUtility::Mean(MatrixXld& observation, VectorIndex& v_index_col, VectorXld& mean)
{
	if(observation.size() == 0 || v_index_col.empty())
		return;
	mean.setZero(observation.rows());
	for(size_t i = 0; i < v_index_col.size(); ++i)
		mean += observation.col(v_index_col[i]);
	mean /= v_index_col.size();
}

void ClusterUtility::Mean(Vector_VectorXld& observation, VectorIndex& v_index_col, VectorXld& mean)
{
	if(observation.size() == 0 || v_index_col.empty())
		return;
	mean.setZero(observation[0].size());
	for(size_t i = 0; i < v_index_col.size(); ++i)
		mean += observation[v_index_col[i]];
	mean /= v_index_col.size();
}

void ClusterUtility::Mean(MatrixXld& observation, VectorXld& mean)
{
	mean = observation.rowwise().mean();
}

void ClusterUtility::Mean(Vector_VectorXld& observation, VectorXld& mean)
{
	if(observation.size() == 0)
		return;
	mean.setZero(observation[0].size());
	for(size_t i = 0; i < observation.size(); ++i)
		mean += observation[i];
	mean /= observation.size();
}

void ClusterUtility::CovarianceMatrix(MatrixXld& observation, MatrixXld& covariance)
{
	VectorXld mean;
	ClusterUtility::Mean(observation, mean);

	VectorXld one;
	one.setOnes(observation.cols());

	MatrixXld valueXOne = mean * one.transpose();
	MatrixXld observation_sub_valueXOne = observation - valueXOne;
	covariance = observation_sub_valueXOne * observation_sub_valueXOne.transpose();
	covariance /= (observation.cols() - 1);
}

void ClusterUtility::CovarianceMatrix(Vector_VectorXld& observation, MatrixXld& covariance)
{
	MatrixXld observation_mtrx;
	ClusterUtility::Convert(observation, observation_mtrx);
	CovarianceMatrix(observation_mtrx, covariance);
}
