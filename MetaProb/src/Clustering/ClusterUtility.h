/*
 * ClusterUtility.h
 *
 *  Created on: 29/feb/2016
 *      Author: samuele
 */

#ifndef SRC_CLUSTERING_CLUSTERUTILITY_H_
#define SRC_CLUSTERING_CLUSTERUTILITY_H_

#include <unordered_map>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
using namespace std;
using namespace Eigen;

class ClusterUtility {
public:
	typedef long double TypeVector;
	typedef uint16_t Cluster;
	typedef uint16_t Dimension;
	typedef Matrix< TypeVector, Dynamic, 1 > VectorXld;
	typedef vector<VectorXld, aligned_allocator<VectorXld>> Vector_VectorXld;
	typedef Map<VectorXld> Data;
	typedef vector<Data, aligned_allocator<Data>> VectorData;
	typedef Matrix< TypeVector, Dynamic, Dynamic > MatrixXld;
	typedef vector<Cluster, aligned_allocator<Cluster>> VectorCluster;
	typedef vector<size_t> VectorIndex;
	typedef unordered_map<Cluster, VectorIndex> mapClusterIndex;

	ClusterUtility();
	virtual ~ClusterUtility();

	void addData(vector<long double>* data);
	void Convert(MatrixXld& observation);
	void Convert(Vector_VectorXld& observation);

	static void GetSubMatrix(MatrixXld& observation, vector<size_t>& index, MatrixXld& sub_observation);
	static void GetSubMatrix(Vector_VectorXld& observation, vector<size_t>& index,MatrixXld& sub_observation);
	static void GetSubMatrix(MatrixXld& observation, vector<size_t>& index, Vector_VectorXld& sub_observation);
	static void GetSubMatrix(Vector_VectorXld& observation, vector<size_t>& index, Vector_VectorXld& sub_observation);
	static void Convert(MatrixXld& from_observation, Vector_VectorXld& to_observation);
	static void Convert(Vector_VectorXld& from_observation, MatrixXld& to_observation);
	static void Mean(MatrixXld& observation, VectorIndex& v_index_col, VectorXld& mean);
	static void Mean(Vector_VectorXld& observation, VectorIndex& v_index_col, VectorXld& mean);
	static void Mean(MatrixXld& observation, VectorXld& mean);
	static void Mean(Vector_VectorXld& observation, VectorXld& mean);
	static void CovarianceMatrix(MatrixXld& observation, MatrixXld& covariance);
	static void CovarianceMatrix(Vector_VectorXld& observation, MatrixXld& covariance);
private:
	vector<vector<long double>*> v_data;
	size_t dim_feature;

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

#endif /* SRC_CLUSTERING_CLUSTERUTILITY_H_ */
