#ifndef SRC_CLUSTERING_KMEANS_H_
#define SRC_CLUSTERING_KMEANS_H_

#include <memory>
#include <vector>
#include <iostream>
#include <chrono>
#include <unordered_map>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include "ClusterUtility.h"

using namespace std;
using namespace std::chrono;
using namespace Eigen;

class Kmeans {
public:
	typedef shared_ptr<Kmeans> Ptr;
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

	Kmeans(MatrixXld& observation);
	virtual ~Kmeans();

	static void GenerateFirstRandomCentroidsFromData(MatrixXld& observation, Vector_VectorXld& centroid, size_t K);
	VectorCluster& compute(uint64_t iter, uint64_t time_in_second, TypeVector eps_ch_thres, Vector_VectorXld& centroid);

	//Get. After computation help to identify the cause of interrupt's iterations
	uint32_t getNumStep() const;
	bool isChangedCls() const;
	double getTotTime() const;
	void GetSplitIndex(mapClusterIndex& map_Observation);
	const VectorCluster& getDataToCluster() const;

private:
	MatrixXld& observation;
	Dimension feature;

	VectorCluster data_to_cluster;
	uint32_t numStep;
	double tot_time;
	TypeVector eps_ch;
	bool changed_cls;

	void computeCentroid(Vector_VectorXld& centroid);
	bool isChangeCentroid(Vector_VectorXld& centroid, Vector_VectorXld& previus_centroid, TypeVector eps_ch_thres, VectorXld& min_distanceBetweenCentroid);
	void compudeDistanceBetweenCentroid(Vector_VectorXld& centroid, MatrixXld& matrix);
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

#endif /* SRC_CLUSTERING_KMEANS_H_ */
