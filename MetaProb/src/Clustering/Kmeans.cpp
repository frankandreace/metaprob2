#include "Kmeans.h"
#include <random>
using namespace std;

Kmeans::Kmeans(MatrixXld& observation):observation(observation) {
	if(observation.size() == 0)
	{
		cout << endl << "Zero data to compute K-means." << flush;
		abort();
	}
	this->feature = observation.rows();

	//Reset
	this->numStep = 0;
	this->changed_cls = true;
	this->tot_time = 0;
	this->eps_ch = 0;
}

Kmeans::~Kmeans() {
}

Kmeans::VectorCluster& Kmeans::compute(uint64_t iter, uint64_t time_in_second, TypeVector eps_ch_thres, Vector_VectorXld& centroid)
{
	if(centroid.size() == 0)
	{
		cout << endl << "Give some centroid for K-means computation." << flush;
		abort();
	}
	for(size_t i = 0; i < centroid.size(); ++i)
		if(centroid[i].size() != this->feature)
		{
			cout << endl << "Centroid of different size compared with data for K-means." << flush;
			abort();
		}

	//Reset previus computation
	this->data_to_cluster.clear();
	this->data_to_cluster.resize(this->observation.cols(), numeric_limits<Cluster>::max());
	Vector_VectorXld previus_centroid(centroid.size(), Vector_VectorXld::value_type::Constant(this->feature, 0));

	this->numStep = 0;
	this->changed_cls = true;
	this->tot_time = 0;
	while(this->changed_cls == true && this->numStep < iter && this->tot_time < time_in_second)
	{
		////////////////////////////////////////
		auto start = high_resolution_clock::now();

		++this->numStep;
//		cout << "\r" << "Clustering step: " << this->numStep << " Time: " << this->tot_time << flush;
		////////////////////////////////////////

		//Copia previus assign for fast computation
		VectorCluster data_to_cluster_previus = this->data_to_cluster;

		//Reset cluster assign
		this->data_to_cluster.clear();
		this->data_to_cluster.resize(this->observation.cols(), numeric_limits<Cluster>::max());

		//Distanza minima tra centroidi
		MatrixXld distanceBetwCentroid;
		this->compudeDistanceBetweenCentroid(centroid, distanceBetwCentroid);
		VectorXld min = distanceBetwCentroid.rowwise().minCoeff();

		//Limite circonferenza intorno a centroidi. Prendo la metà in modo da non aver insiemi sovrapposti
		min /= 2;

		//Inserisco i punti che stanno all'interno della circonferenza nei rispettivi cluster
		//Inizializzo distanze a -1
		MatrixXld distance_centroid_to_data = MatrixXld::Constant(centroid.size(), this->observation.cols(), -1);

		//Computo prima la distanza con cluster assegnato precedentemente, se non era assegnato salto passaggio
		#pragma omp parallel for
		for(size_t i = 0; i < data_to_cluster_previus.size(); ++i)
		{
			if(data_to_cluster_previus[i] != numeric_limits<Cluster>::max()) //Se era già assegnato
			{
				//Computo distanza con cluster precedente
				VectorXld v_dist = centroid[data_to_cluster_previus[i]] - this->observation.col(i);
				long double dist = v_dist.norm();
				distance_centroid_to_data(data_to_cluster_previus[i],i) = dist; //Aggiorno distanza
				//Se distanza minore del min/2 inserisco in quel cluster
				if(distance_centroid_to_data(data_to_cluster_previus[i],i) < min(data_to_cluster_previus[i]))
					this->data_to_cluster[i] = data_to_cluster_previus[i];
			}
		}

		//Cluster sequenziali per evitare conti inutili se riesco ad assegnare prima,
		//salto già assegnati grazie a passo precedente
		for(Cluster cls = 0; cls < centroid.size(); ++cls)
		{
			//Posso eseguire in parallelo perché assegno cluster in locazioni diverse (indice i)
			#pragma omp parallel for
			for(size_t i = 0; i < this->data_to_cluster.size(); ++i)
			{
				if(this->data_to_cluster[i] == numeric_limits<Cluster>::max()) //Se non già assegnato
				{
					VectorXld v_dist = centroid[cls] - this->observation.col(i);
					long double dist = v_dist.norm();
					distance_centroid_to_data(cls,i) = dist; //Aggiorno distanza
					//Se distanza minore del min/2 inserisco in quel cluster
					if(distance_centroid_to_data(cls,i) < min(cls))
						this->data_to_cluster[i] = cls;
				}
			}
		}

		//Inserisco punti fuori dalla circonferenza e non assegnati prima
		#pragma omp parallel for
		for(size_t i = 0; i < (size_t)this->observation.cols(); ++i)
		{
			if(this->data_to_cluster[i] == numeric_limits<Cluster>::max()) //Se non già assegnato
			{
				//Cerca minimo tra cluster
				double minDataValue = numeric_limits<double>::max();
				Cluster minDataRow = 0;
				for(Cluster cls = 0; cls < centroid.size(); cls++)
				{
					if(distance_centroid_to_data(cls,i) < 0) //Non computata, computa, altrimenti salta
					{
						VectorXld v_dist = centroid[cls] - this->observation.col(i);
						long double dist = v_dist.norm();
						distance_centroid_to_data(cls,i) = dist;
					}
					if(distance_centroid_to_data(cls,i) < minDataValue)
					{
						minDataValue = distance_centroid_to_data(cls,i);
						minDataRow = cls;
					}
				}
				this->data_to_cluster[i] = minDataRow;
			}
		}

		previus_centroid = centroid; //Store previus centroid
		this->computeCentroid(centroid);
		this->changed_cls = this->isChangeCentroid(centroid, previus_centroid, eps_ch_thres, min);

		////////////////////////////////////////
		auto end = high_resolution_clock::now();
		float elapsedSeconds = duration_cast<duration<float>>(end-start).count();

		this->tot_time += elapsedSeconds;
		////////////////////////////////////////
	}
	return this->data_to_cluster;
}

//Return if change
//Call GenerateFirstRandomCentroid before compute
void Kmeans::computeCentroid(Vector_VectorXld& centroid) {
	//Salvo precedenti
	size_t K = centroid.size();
	centroid.clear();
	centroid.resize(K, Vector_VectorXld::value_type::Constant(this->feature, 0));

	mapClusterIndex map_cluster_index;
	this->GetSplitIndex(map_cluster_index);

	#pragma omp parallel for
	for(size_t cls = 0; cls < centroid.size(); ++cls)
		ClusterUtility::Mean(this->observation, map_cluster_index[cls], centroid[cls]);
}

bool Kmeans::isChangeCentroid(Vector_VectorXld& centroid, Vector_VectorXld& previus_centroid, TypeVector eps_ch_thres, VectorXld& min_distanceBetweenCentroid)
{
	TypeVector min = min_distanceBetweenCentroid.minCoeff()*2; //Perché lo passo diviso per 2
	VectorXld norm(centroid.size());
	#pragma omp parallel for
	for(size_t i = 0; i < centroid.size(); i++)
		norm(i) = (centroid[i] - previus_centroid[i]).norm();
	this->eps_ch = norm.sum();
	TypeVector thresh = eps_ch_thres*min;
	return this->eps_ch > thresh;
}

bool Kmeans::isChangedCls() const {
	return changed_cls;
}

double Kmeans::getTotTime() const {
	return tot_time;
}

uint32_t Kmeans::getNumStep() const {
	return numStep;
}

void Kmeans::GetSplitIndex(mapClusterIndex& map_Observation) {
	map_Observation.clear();
	for(size_t i = 0; i < this->data_to_cluster.size(); ++i)
		map_Observation[this->data_to_cluster[i]].push_back(i);
}


const Kmeans::VectorCluster& Kmeans::getDataToCluster() const {
	return data_to_cluster;
}


void Kmeans::compudeDistanceBetweenCentroid(Vector_VectorXld& centroid, MatrixXld& matrix) {
	matrix.setZero(centroid.size(), centroid.size());
	#pragma omp parallel for
	for(size_t i = 0; i < centroid.size(); i++)
		#pragma omp parallel for firstprivate(i)
		for(size_t j = i; j < centroid.size(); j++)
		{
			long double distance;
			if(i == j)
				distance = numeric_limits<TypeVector>::max();
			else
			{
				VectorXld v_dist = centroid[i] - centroid[j];
				distance = v_dist.norm();
			}
			matrix(i,j) = distance;
			matrix(j,i) = distance;
		}
}

//Call this before compute
void Kmeans::GenerateFirstRandomCentroidsFromData(MatrixXld& observation, Vector_VectorXld& centroid, size_t K)
{
	centroid.clear();
	centroid.resize(K);

	//Initialize random generator
	vector<Cluster> v_size_data(centroid.size(), 1); //Initialize vector weight
	default_random_engine generator(time(0));
	discrete_distribution<size_t> rnd_cluster(v_size_data.begin(), v_size_data.end());

	//Init centroid
	mapClusterIndex map_cluster_index;
	for(size_t i = 0; i < (size_t)observation.cols(); ++i)
		map_cluster_index[rnd_cluster(generator)].push_back(i);

	#pragma omp parallel for
	for(size_t cls = 0; cls < centroid.size(); ++cls)
		ClusterUtility::Mean(observation, map_cluster_index[cls], centroid[cls]);
}
