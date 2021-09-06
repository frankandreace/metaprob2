/*
 * Gmeans.h
 *
 *  Created on: 29/feb/2016
 *      Author: samuele
 */

#ifndef SRC_CLUSTERING_GMEANS_H_
#define SRC_CLUSTERING_GMEANS_H_

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stack>
#include <boost/math/distributions/normal.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include "KStest.h"
#include "ClusterUtility.h"
#include "Kmeans.h"

using namespace std;
using namespace Eigen;

struct ClusterSplitter
{
public:
	typedef ClusterUtility::TypeVector TypeVector;
	typedef ClusterUtility::Dimension Dimension;
	typedef ClusterUtility::Cluster Cluster;
	typedef ClusterUtility::VectorXld VectorXld;
	typedef ClusterUtility::Vector_VectorXld Vector_VectorXld;
	typedef ClusterUtility::Data Data;
	typedef ClusterUtility::VectorData VectorData;
	typedef ClusterUtility::MatrixXld MatrixXld;
	typedef ClusterUtility::VectorCluster VectorCluster;
	typedef ClusterUtility::VectorIndex VectorIndex;
	typedef ClusterUtility::mapClusterIndex mapClusterIndex;
	typedef SelfAdjointEigenSolver<MatrixXld> SelfAdjointEigenSolver_;

	struct Result
	{
	public:
		TypeVector p_value = 1;
		Vector_VectorXld centroidSplit;
		bool isSplitted = false;
	};

	ClusterSplitter(MatrixXld& observation, VectorIndex& indexCluster); //Call with all observation and index that you want to considerer

private:
	MatrixXld observation;
	VectorXld centroid;
	void GetSplitCentroid(Vector_VectorXld& centroidSplit);
//	TypeVector convertAD2_star_To_PValue(TypeVector AD_2_star);

public:
	void isClusterSplitted(uint64_t iter, uint64_t time_in_second, TypeVector eps_ch_thres, TypeVector p_value_thresh, Result& result); //Return if split execute correct
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class Gmeans {
public:
	typedef ClusterUtility::TypeVector TypeVector;
	typedef ClusterUtility::Dimension Dimension;
	typedef ClusterUtility::Cluster Cluster;
	typedef ClusterUtility::VectorXld VectorXld;
	typedef ClusterUtility::Vector_VectorXld Vector_VectorXld;
	typedef ClusterUtility::Data Data;
	typedef ClusterUtility::VectorData VectorData;
	typedef ClusterUtility::MatrixXld MatrixXld;
	typedef ClusterUtility::VectorCluster VectorCluster;
	typedef ClusterUtility::VectorIndex VectorIndex;
	typedef ClusterUtility::mapClusterIndex mapClusterIndex;

	Gmeans(MatrixXld& observation);
	virtual ~Gmeans();
	Kmeans::Ptr kmeans;

	VectorCluster& compute(Cluster start_k, uint64_t iter_split_max, uint64_t iter_kmeans, uint64_t time_in_second, TypeVector eps_ch_thres, TypeVector p_value);
	Cluster getEstimateK() const;

private:
	MatrixXld& observation; //Row is feature size, cols is vector observed
	Cluster estimateK;
	Vector_VectorXld centroid_saved;
	VectorCluster data_to_cluster;
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

#endif /* SRC_CLUSTERING_GMEANS_H_ */
