#include "Solution.h"

#include <chrono>

#include<sys/stat.h>
#include "HandlerInput/Parameter.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[])
{
	Eigen::initParallel();
	//Initialize parameter DEFAULT
	bool grp_mode = false;
	string dir_output = "output/";
	Parameter::Ptr param = make_shared<Parameter>();

	//Default parameter
	param->nClus = 0;
	param->q = 30;
	param->m_thres = 5;
	param->SEEDSIZE = 9000;
	param->l_mer_freq = 4;
	param->iteration = 100;
	param->time_max = 60*60; //in second
	param->graph = Paired;
	param->norm_type = Normalization::NORM_D2star_All_Read_Prob_Lmer_Euclidian;
	param->estimateK = true;

	//----------------------GET INPUT----------------------//

	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i], "-si") == 0)
		{
			i++;
			PairFiles& file = param->insertSingleEndFile(argv[i]);
			if(!file.isCorrect())
			{
				cerr<<endl<<"Please enter an input filename single-end: -si <AbsPathFile>"<<flush;
				return 0;
			}
		}
		else if(strcmp(argv[i], "-pi") == 0)
		{
			i++;
			PairFiles& file = param->insertPairedEndFiles(argv[i], argv[i+1]);
			if(!file.isCorrect())
			{
				cerr<<endl<<"Please enter an input filenames paired-end: -pi <AbsPathFile1> <AbsPathFile2>"<<flush;
				return 0;
			}
		}
		else if(strcmp(argv[i], "-dirOutput") == 0)
		{
			i++;
			dir_output.assign(argv[i]);
			if(dir_output == "")
			{
				cerr<<endl<<"Please enter an output directory if you specify -dirOutput"<<flush;
				return 0;
			}
		}
		else if(strcmp(argv[i], "-graphType") == 0)
		{
			switch(atoll(argv[i]))
			{
			case 0:
				param->graph = Paired;
				break;
			case 1:
				param->graph = Single;
				break;
			case 2:
				param->graph = SingleUnion;
				break;
			}
		}
		else if(strcmp(argv[i], "-numSp") == 0)
		{
			i++;
			long long int var = atoll(argv[i]);
			if(var <= 0)
			{
				cerr<<endl<<"Please enter -numSp > 0"<<flush;
				return 0;
			}
			else
			{
				param->nClus = var;
				param->estimateK = false;
			}
		}
		else if(strcmp(argv[i], "-q") == 0)
		{
			i++;
			long long int var = atoll(argv[i]);
			if(var <= 0)
			{
				cerr<<endl<<"Please enter -q > 0"<<flush;
				return 0;
			}
			else
				param->q = var;
		}
		else if(strcmp(argv[i], "-m") == 0)
		{
			i++;
			long long int var = atoll(argv[i]);
			if(var <= 0)
			{
				cerr<<endl<<"Please enter -m > 0"<<flush;
				return 0;
			}
			else
				param->m_thres = var;
		}
		else if(strcmp(argv[i], "-ssize") == 0)
		{
			i++;
			long long int var = atoll(argv[i]);
			if(var <= 0)
			{
				cerr<<endl<<"Please enter -ssize > 0"<<flush;
				return 0;
			}
			else
				param->SEEDSIZE = var;
		}
		else if(strcmp(argv[i], "-lmerFreq") == 0)
		{
			i++;
			long long int var = atoll(argv[i]);
			if(var <= 0)
			{
				cerr<<endl<<"Please enter -lmerFreq > 0"<<flush;
				return 0;
			}
			else
				param->l_mer_freq = var;
		}
		else if(strcmp(argv[i], "-iterMaxKmeans") == 0)
		{
			i++;
			long long int var = atoll(argv[i]);
			if(var <= 0)
			{
				cerr<<endl<<"Please enter -iterMaxKmeans > 0"<<flush;
				return 0;
			}
			else
				param->iteration = var;
		}
		else if(strcmp(argv[i], "-timeMaxKmeans") == 0)
		{
			i++;
			long long int var = atoll(argv[i]);
			if(var <= 0)
			{
				cerr<<endl<<"Please enter -timeMaxKmeans > 0"<<flush;
				return 0;
			}
			else
				param->time_max = var;
		}
		else if(strcmp(argv[i], "-feature") == 0)
		{
			i++;
			long long int var = atoll(argv[i]);
			if(var <= 0 || var >= Normalization::Last)
			{
				cerr<<endl<<"-feature != [1, ..., 12]. Applied default parameter. "<< Normalization::enum_string(Normalization::NORM_D2star_All_Read_Prob_Lmer_Euclidian) <<flush;
				param->norm_type = Normalization::NORM_D2star_All_Read_Prob_Lmer_Euclidian;
			}
			else
				param->norm_type = Normalization::int_to_enum(var);
		}
		else if(strcmp(argv[i], "-mg") == 0)
		{
			grp_mode = true;
			cout << endl << "Only group mode activated" << flush;
		}
		else if(strcmp(argv[i], "-eK") == 0)
		{
			param->estimateK = true;
		}
	}
	//Creo cartella output se non presente
	createDirAndSubDir(dir_output);

	//Print parameter input
	cout << endl << "Directory output: " << dir_output << flush;
	param->printFiles();
	param->printParameter();

	//--------------------------------------------------------------------//

	auto start = high_resolution_clock::now();
	//Class Solution init
	clsSolut::Ptr mSol = make_shared<clsSolut>();
	bool init_correct = mSol->InitData(param);

	if(init_correct && !grp_mode)
	{
		if(mSol->DoBinning(param->norm_type) == true)
			cout<<endl<<"Binning successfully!"<<flush;
		else
		{
			cout<<endl<<"-----------------------------------------"<<flush;
			cout<<endl<<"Can not bin data!"<<flush;
			return 0;
		}
	}
	auto end = high_resolution_clock::now();
	float elapsedSeconds = duration_cast<duration<float>>(end-start).count();

	//--------------------------------------------------------------------//
	//Stampa risultati a video
	mSol->PrintResult();
	mSol->SaveInfo(to_string(elapsedSeconds), dir_output);
	mSol->SaveBinning(dir_output);

	return 0;
}
