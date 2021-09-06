#include"Solution.h"

//Init data
bool clsSolut::InitData(Parameter::Ptr param)
{
	//-----------------------------------------------------//
	this->param = param;

	//-----------------------------------------------------//
	//File DNA da non eliminare, da resettare se cambia file output solo perchè cambia Norma o Run
	this->mFile = make_shared<clsDNAFile>(this->param->input_files, this->param->q, this->param->graph, this->param->m_thres);

	////////////////////////////////////////////////////////////////////////////////////////////////////
	//Informazioni su file e numero di read corrette nel file da memorizzare

	auto start = high_resolution_clock::now();
	if(!this->mFile->GetFileInfo())
		return false; //Lettura errata
	auto end = high_resolution_clock::now();
	float elapsedSeconds = duration_cast<duration<float>>(end-start).count();
	info.push_back("Read info File (seconds): " + to_string(elapsedSeconds));
	info.push_back("Loaded sequences: " + to_string(this->mFile->getSeqTot()));
	if(this->mFile->getSeqTot() == 0)
		return false;

	////////////////////////////////////////////////////////////////////////////////////////////////////
	//Crea BloomFilter per scoprire i q-mer ripetuti e poter inserire solo quelli

//	start = high_resolution_clock::now();
//	if(!this->mFile->CreateBloomFilter())
//		return false; //Lettura errata
//	end = high_resolution_clock::now();
//	elapsedSeconds = duration_cast<duration<float>>(end-start).count();
//	info.push_back("Creation BloomFilter (seconds): " + to_string(elapsedSeconds));

	////////////////////////////////////////////////////////////////////////////////////////////////////
	//Estraggo i q-mer ripetuti, o tutti se il BloomFilter non è applicato

	start = high_resolution_clock::now();
	if(!this->mFile->ExtractQMerFromFile())
		return false; //Lettura errata
	end = high_resolution_clock::now();
	elapsedSeconds = duration_cast<duration<float>>(end-start).count();
	info.push_back("Creation q-mer graph (seconds): " + to_string(elapsedSeconds));
	info.push_back("hash q-mer before filter: " + to_string(this->mFile->getHash()->arrMyHash.size()));
	//-----------------------------------------------------//

	//Initialize graph (creation is after)
	this->cSeqGrph = make_shared<clsSeqGraph>(this->mFile);

	//Build graph of reads, and remove the hash table
	start = high_resolution_clock::now();
	this->cSeqGrph->CreateGraph();
	end = high_resolution_clock::now();
	elapsedSeconds = duration_cast<duration<float>>(end-start).count();
	info.push_back("Creation adj graph (seconds): " +  to_string(elapsedSeconds));
	info.push_back("hash q-mer filtered: " + to_string(this->cSeqGrph->getQMerFilter()));

	info.push_back("Peak Virtual Memory used (MB): " + to_string(getPeakVirtualMemoryUsed()/1024));
	info.push_back("Virtual Memory used (MB): " + to_string(getVirtualMemoryUsed()/1024));

	//eliminated hash no longer necessary
	this->mFile->resetHash();

	bool firstPhase = this->FirstPhase();
	return firstPhase;
}

bool clsSolut::FirstPhase() {
	//First phase
	cout<<endl<<"-----------------------------------------"<<flush;
	cout<<endl<<"PHASE I: CREATE GROUPS"<<flush;

	auto start = high_resolution_clock::now();
	bool createGrp = this->Grouping();
	auto end = high_resolution_clock::now();
	float elapsedSeconds = duration_cast<duration<float>>(end-start).count();
	info.push_back("Creation groups (seconds): " +  to_string(elapsedSeconds));

	if(!createGrp)
		return false;

	//Use paired end information to union the groups
	if(this->param->graph == SingleUnion)
	{
		start = high_resolution_clock::now();
		this->UnionGrpIfPossible();
		end = high_resolution_clock::now();
		elapsedSeconds = duration_cast<duration<float>>(end-start).count();
		info.push_back("Union groups (seconds): " +  to_string(elapsedSeconds));
	}

	this->SetStateSeed();
	info.push_back("Groups created: " + to_string(this->vGroup.size()));
	return true;
}

//start to bin
bool clsSolut::DoBinning(Normalization::Norm norm_type)
{
	if(this->vClusts.empty())
	{
		//Third phase
		cout<<endl<<"-----------------------------------------"<<flush;
		cout<<endl<<"PHASE II: CREATE CLUSTERS"<<flush;

		auto start = high_resolution_clock::now();
		bool clusterExe = this->SecondPhase(norm_type);
		auto end = high_resolution_clock::now();
		float elapsedSeconds = duration_cast<duration<float>>(end-start).count();
		info.push_back("Clustering (seconds): " +  to_string(elapsedSeconds));

		if(!clusterExe)
			return false;
	}
	cout<<endl<<"-----------------------------------------"<<flush;
	return true;
}

//The second phase, grouping reads
bool clsSolut::Grouping()
{
	cout << endl << flush;
	cout << "Create group... " << flush;

	//GROUPING (ID from 0)
	//----------Step 1: build groups--------------//
	auto nextR = this->mFile->begin();
	while(nextR !=  this->mFile->end())
	{
		clsGroup::Ptr cGrp = make_shared<clsGroup>();
		cGrp->InitData(this->vGroup.size()); //Attribuisce ID

		//Frontiera del gruppo, analizzo sempre nodi seed se inseriti
		typedef pair<size_seq, id_seq_type> Weight_Recurrence__Border; //Salva il peso dell'arco e quante volte viene visto analizzando la frontiera
		typedef unordered_map<id_seq_type, Weight_Recurrence__Border> MapAdjSeqID__Weight_Recurrence__Border;
		typedef unordered_map<SeedState, MapAdjSeqID__Weight_Recurrence__Border, hash<int>> Border;
		Border borderOfGroups;

		//Inserisco il primo read nel Bordo tra i nodi seed
		borderOfGroups[Seed][(*nextR)->getSeqId()].first =  1; //Peso ininfluente, è il primo
		borderOfGroups[Seed][(*nextR)->getSeqId()].second = 1; //Trovato una volta

		while(cGrp->get_Seeds_Size() < this->param->SEEDSIZE) //Taglio gruppo a dato numero di nucleotidi
		{
			//Elimino i read che sono nelle frontiere ma son già stati assegnati e prendo i massimi
			//Frontiera Seed
			auto BorderSeed_it = borderOfGroups[Seed].begin();
			auto BorderMaxSeed = BorderSeed_it;
			while(BorderSeed_it != borderOfGroups[Seed].end())
			{
				clsSequencePE::Ptr seq = this->mFile->GetAt(BorderSeed_it->first);
				if(seq->getSeed() != No_State_Seed) //Quindi è o seed o no_seed. Già assegnato. Da eliminare
				{
					if(BorderMaxSeed == BorderSeed_it)
					{
						BorderSeed_it = borderOfGroups[Seed].erase(BorderSeed_it);
						BorderMaxSeed = BorderSeed_it;
					}
					else
						BorderSeed_it = borderOfGroups[Seed].erase(BorderSeed_it);
				}
				else
				{
					if(BorderMaxSeed->second.second < BorderSeed_it->second.second)
						BorderMaxSeed = BorderSeed_it;
					++BorderSeed_it;
				}
			}
			//Frontiera NoSeed
			auto BorderNoSeed_it = borderOfGroups[No_Seed].begin();
			auto BorderMaxNoSeed = BorderNoSeed_it;
			while(BorderNoSeed_it != borderOfGroups[No_Seed].end())
			{
				clsSequencePE::Ptr seq = this->mFile->GetAt(BorderNoSeed_it->first);
				if(seq->getSeed() != No_State_Seed) //Quindi è o seed o no_seed. Già assegnato. Da eliminare
				{
					if(BorderMaxNoSeed == BorderNoSeed_it)
					{
						BorderNoSeed_it = borderOfGroups[No_Seed].erase(BorderNoSeed_it);
						BorderMaxNoSeed = BorderNoSeed_it;
					}
					else
						BorderNoSeed_it = borderOfGroups[No_Seed].erase(BorderNoSeed_it);
				}
				else
				{
					if(BorderMaxNoSeed->second.second < BorderNoSeed_it->second.second)
						BorderMaxNoSeed = BorderNoSeed_it;
					++BorderNoSeed_it;
				}
			}

			///////////////////////////////////////////////////////////////////////////////
			//Seleziono prossimo da inserire prendendo il massimo tra i due insiemi
			MapAdjSeqID__Weight_Recurrence__Border::iterator max;
			SeedState max_state = No_Other_Read;
			if(BorderMaxSeed != borderOfGroups[Seed].end() && BorderMaxNoSeed != borderOfGroups[No_Seed].end())
			{
				if(BorderMaxSeed->second.second < BorderMaxNoSeed->second.second)
				{
					max = BorderMaxNoSeed;
					max_state = No_Seed;
				}
				else
				{
					max = BorderMaxSeed;
					max_state = Seed;
				}
			}
			else if (BorderMaxSeed == borderOfGroups[Seed].end() && BorderMaxNoSeed != borderOfGroups[No_Seed].end())
			{
				max = BorderMaxNoSeed;
				max_state = No_Seed;
			}
			else if (BorderMaxSeed != borderOfGroups[Seed].end() && BorderMaxNoSeed == borderOfGroups[No_Seed].end())
			{
				max = BorderMaxSeed;
				max_state = Seed;
			}

			if (max_state == No_Other_Read)
				break;

			id_seq_type invalidID = 0; //ID = 0 so che non esiste
			id_seq_type iCRead = invalidID;
			SeedState flagiCRead = No_State_Seed;
			if(max_state == No_Seed)
			{
				flagiCRead = max_state;
				iCRead = max->first;
			}
			else if(max_state == Seed)
			{
				//Verifica non sia presente anche su No_Seed, se presente anche li allora è un No_Seed
				auto NoSeed_it = borderOfGroups[No_Seed].find(max->first);
				if(NoSeed_it == borderOfGroups[No_Seed].end()) //Non trovato su frontiera NoSeed, quindi è Seed
					flagiCRead = max_state;
				else //Trovato anche su frontiera no seed. E' un non seed e devo eliminarlo anche da quell'insieme
					flagiCRead = No_Seed;
				iCRead = max->first;
			}

			///////////////////////////////////////////////////////////////////////////////

			//Aggiungi a gruppo come flagiCRead(Seed o Non Seed dipende da quale seleziono)
			clsSequencePE::Ptr seq = this->mFile->GetAt(iCRead);
			cGrp->AddItem(seq, flagiCRead);

			//Prendi adiacenti nodo appena inserito nel gruppo. Arco su sè stesso eliminato in GetAdiacence
			clsSequencePE::MapAdj& adj = *seq->vAdjSeqPE;

			//Incrementa frontiera con gli adiacenti del read appena inserito nel gruppo
			auto adj_analize_it = adj.begin();
			while(adj_analize_it != adj.end())
			{
				//Se lato supera soglia. Arco su sè stesso eliminato in GetAdiacence
				if(adj_analize_it->second >= this->param->m_thres)
				{
					clsSequencePE::Ptr seq = this->mFile->GetAt(adj_analize_it->first);
					//Se non ha alcun stato di seed o non seed (vicini non analizzati)
					if(seq->getSeed() == No_State_Seed)
					{
						auto border_noSeed_it = borderOfGroups[No_Seed].find(adj_analize_it->first);
						bool isInBorderNoSeed = border_noSeed_it != borderOfGroups[No_Seed].end();

						auto border_Seed_it = borderOfGroups[Seed].find(adj_analize_it->first);
						bool isInBorderSeed = border_Seed_it != borderOfGroups[Seed].end();

						//Aggiorna pesi e read condivisi con bordo
						if(isInBorderNoSeed) //Evita un accesso che è già stato fatto
						{
							border_noSeed_it->second.first += adj_analize_it->second;
							++border_noSeed_it->second.second;
						}
						if(isInBorderSeed)//For evaluate most adj  read
						{
							border_Seed_it->second.first += adj_analize_it->second;
							++border_Seed_it->second.second;
						}
						if(!isInBorderSeed || !isInBorderNoSeed)
						{
							//Vedere dove inserire
							SeedState whereAddBorder = No_Seed;
							//Se precedente non è seed e non è sul bordo come NoSeed allora
							//questo può essere Seed
							if(!isInBorderNoSeed && flagiCRead == No_Seed)
								whereAddBorder = Seed;

							if(!isInBorderSeed && !isInBorderNoSeed)
							{
								Weight_Recurrence__Border& weight_recurrence = borderOfGroups[whereAddBorder][adj_analize_it->first];
								weight_recurrence.first += adj_analize_it->second;
								++weight_recurrence.second;
							}
							else
							{
								//previously update
								if(isInBorderSeed && whereAddBorder == No_Seed)
								{
									Weight_Recurrence__Border& weight_recurrence = borderOfGroups[whereAddBorder][adj_analize_it->first];
									weight_recurrence.first = border_Seed_it->second.first;
									weight_recurrence.second = border_Seed_it->second.second;
								}
							}
						}
					}
				}
				++adj_analize_it;
			}
		}
		this->vGroup.push_back(cGrp);

//		//Stampa
//		cout << "\r" << "Create group: " << this->vGroup.size() << flush;

		//find the next read which will be assigned into a new group at first
		++nextR; //parti da nextR+1
		while(nextR != this->mFile->end())
		{
			if((*nextR)->getSeed() == No_State_Seed)
				break;
			++nextR;
		}
	}

	cout << " Complete. Create: " << this->vGroup.size() << " groups" << flush;

	//Ho isChoose = true su tutti i read
	//Bisogna impostare i Seed (Funzione SetStateSeed)

	return true;
}

bool clsSolut::IsThereAdjBetweenGroup()
{
	bool findAdj = false;
	//Creo grafo di adiacenze tra gruppi e poi li unisco
	auto ID_pair = this->mFile->beginMapID();
	while(ID_pair != this->mFile->endMapID())
	{
		//Se tutti e due diversi da 0 ho informazione paired end
		if(ID_pair->left != 0 && ID_pair->right != 0) //ID 0 non esiste e valore default
		{
			id_grp_type id_first = this->mFile->GetAt(ID_pair->left)->getIdGrp();
			id_grp_type id_second = this->mFile->GetAt(ID_pair->right)->getIdGrp();

			//Sono su gruppi diversi e l'unione non supera soglia
			if(id_first != id_second && (this->vGroup[id_first]->GetSize() + this->vGroup[id_second]->GetSize()) < this->param->SEEDSIZE)
			{
				findAdj=true;
				//Aggiungo adiacenze
				this->vGroup[id_first]->AddAdjacenceGrp(id_second);
				this->vGroup[id_second]->AddAdjacenceGrp(id_first);
			}
		}
//		//verifica e stampa Single end se presenti in Paired End
//		else if(this->param->flagRead == ReadType::Paired_End && (ID_pair->second.first == 0 || ID_pair->second.second == 0))
//			cout << endl << ID_pair->first << " " << ID_pair->second.first << " " << ID_pair->second.second << flush;
		++ID_pair;
	}
	return findAdj;
}

void clsSolut::UnionGrpIfPossible()
{
	bool change = true;
	bool possibleUnion = this->IsThereAdjBetweenGroup();
	if(possibleUnion)
	{
		cout << endl << flush;
		cout << "Union groups... " << flush;

		while(possibleUnion && change)
		{
			cout << endl << flush;
			change = false; //Si presuppone che non cambi
			vector<clsGroup::Ptr> new_grps;
			for(id_grp_type grp = 0; grp < this->vGroup.size(); grp++)
			{
				if(this->vGroup[grp]->isChoose() == false) //Gruppo non ancora analizzato
				{
					//Crea nuovo gruppo
					//I gruppi nuovi avranno choose = false
					clsGroup::Ptr newGrp = make_shared<clsGroup>();
					newGrp->InitData(new_grps.size());

					//Aggiungi primo elemento
					newGrp->AddItemG(this->vGroup[grp], this->param->m_thres);
					this->vGroup[grp]->setChoose(true);

					//Prendi i suoi vicini come margine esterno e scorrili
					clsGroup::MapAdjGrp border = this->vGroup[grp]->adj; //Costruttore di copia
					auto adj_border_it = border.begin();
					while(adj_border_it != border.end())
					{
						bool possibleToAdd = this->vGroup[adj_border_it->first]->isChoose() == false;
						bool threshSize = newGrp->get_Seeds_Size() + this->vGroup[adj_border_it->first]->get_Seeds_Size() < this->param->SEEDSIZE;
						//Se è stato scelto già in un altro gruppo salta e cancella
						if(possibleToAdd && threshSize)
						{
							//Il gruppo si è allargato
							change = true;
							//Aggiungi a gruppo e setta già inserito in qualche gruppo
							newGrp->AddItemG(this->vGroup[adj_border_it->first], this->param->m_thres);
							this->vGroup[adj_border_it->first]->setChoose(true);

							//Scorro suoi vicini e aggiungo a bordo quelli con choose = false e aggiorno pesi bordo nuovo
							auto adj_current_border_grp_it = this->vGroup[adj_border_it->first]->adj.begin();
							while(adj_current_border_grp_it != this->vGroup[adj_border_it->first]->adj.end())
							{
								if(this->vGroup[adj_current_border_grp_it->first]->isChoose() == false)
								{
									auto find_in_border = border.find(adj_current_border_grp_it->first);
									if(find_in_border == border.end()) //Non trovato
										border.insert(clsGroup::MapAdjGrp::value_type(*adj_current_border_grp_it));
									else
										find_in_border->second += adj_current_border_grp_it->second;
								}
								++adj_current_border_grp_it;
							}
						}

						//Elimina elemento analizzato
						border.erase(adj_border_it);
						//Inizia con il primo elemento rimanente
						adj_border_it = border.begin();
					}

					//Aggiungi nuovo gruppo
					new_grps.push_back(newGrp);

//					//Stampa
//					cout << "\r" << "Creation group union: " << new_grps.size() << flush;
				}
			}
			//riassegno correttamente ID gruppi
			this->vGroup = new_grps; //Copy costructor
			possibleUnion = this->IsThereAdjBetweenGroup();
		}

		cout << " Complete. After union: " << this->vGroup.size() << " groups" << flush;
	}
}

void clsSolut::SetStateSeed()
{
	//Imposta read seed
	for(id_grp_type i=0;i<this->vGroup.size();i++)
		for(id_seq_type j=0;j<this->vGroup[i]->vSeedSeqPE.size();j++)
			this->vGroup[i]->vSeedSeqPE[j]->setSeed(Seed);

	//Imposta read non seed
	for(id_grp_type i=0;i<this->vGroup.size();i++)
		for(id_seq_type j=0;j<this->vGroup[i]->vSeqPE.size();j++)
			this->vGroup[i]->vSeqPE[j]->setSeed(No_Seed);
}

//The third phase, binning of contigs
bool clsSolut::SecondPhase(Normalization::Norm norm_type)
{
	//Count l-mers frequency from file, then Compute l-mers frequency feature
	Normalization::Ptr norm = make_shared<Normalization>(this->mFile, this->vGroup, norm_type, this->param->l_mer_freq);
	norm->Compute();

	//Get feature vector size
	this->feature_size = norm->getLmerCompress()->getLmer().size();

	//Merging Group into cluster by k-means
	this->MergeGroup();

	//Assign all read to a cluster
	this->AssignsSequenceToCluster();
	return true;
}

//SAMU EDIT
//Merge groups, with fixed number of cluster (species) and Assigning groups into the clusters
void clsSolut::MergeGroup()
{
	if(this->param->nClus == 0 && !this->param->estimateK)
	{
		cout << endl << "Not enough cluster. Set -numSp or -eK." << flush;
		abort();
	}
	//Create ClusterUtility to convert data
	ClusterUtility::MatrixXld observation;
	ClusterUtility cls_utility;
	for(id_grp_type grp = 0; grp < this->vGroup.size(); ++grp)
		cls_utility.addData(&this->vGroup[grp]->vFreq);
	cls_utility.Convert(observation);

	//Extra parameter
	uint64_t max_iter_split = observation.cols() - 1; //al massimo cluster come numero di dati
	Gmeans::TypeVector p_value_thresh = 0.05;
	Kmeans::TypeVector eps_ch_thres = 0.0001; //Spostamento minimo per cui ritengo i centroidi cambiati = eps_ch_thres * min_distance_between_centroid

	ClusterUtility::VectorCluster data_to_cluster;
	if(this->param->estimateK)
	{
		//Estimate K
		size_t start_K = (this->param->nClus != 0) ? this->param->nClus : 1;
		cout << endl << "Estimate K and clustering... " << flush;
		Gmeans gmeans(observation);
		data_to_cluster = gmeans.compute(start_K, max_iter_split, this->param->iteration, this->param->time_max, eps_ch_thres, p_value_thresh);//gmeans.compute(start_K, this->param->iteration, this->param->time_max, A_2_star_thresh);
		this->param->nClus = gmeans.getEstimateK();
		cout << "\r"<< "Complete. K is: " << gmeans.getEstimateK() << flush;

		info.push_back("N. Cluster: " + to_string(this->param->nClus));
		info.push_back("N. iter clustering: " + to_string(gmeans.kmeans->getNumStep()));
		string conv = gmeans.kmeans->isChangedCls() ? "No" : "Si";
		info.push_back("Converged: " + conv);
	}
	else
	{
		cout << endl << "Clustering... " << flush;
		Kmeans::Vector_VectorXld centroid;
		Kmeans::GenerateFirstRandomCentroidsFromData(observation, centroid, this->param->nClus);
		Kmeans kmeans(observation);
		data_to_cluster = kmeans.compute(this->param->iteration, this->param->time_max, eps_ch_thres, centroid);
		cout << "Complete" << flush;

		//Save output kmeans
		info.push_back("N. Cluster: " + to_string(this->param->nClus));
		info.push_back("N. iter clustering: " + to_string(kmeans.getNumStep()));
		string conv = kmeans.isChangedCls() ? "No" : "Si";
		info.push_back("Converged: " + conv);
	}

	for(id_cluster_type i = 0; i < this->param->nClus; ++i)
	{
		clsCluster::Ptr newClus = make_shared<clsCluster>();
		newClus->InitID(i);
		this->vClusts.push_back(newClus);
	}
	for(id_grp_type grp = 0; grp < this->vGroup.size(); ++grp)
		this->vClusts[data_to_cluster[grp]]->AddItem(this->vGroup[grp]);
}

//Assign ID of sequence to a cluster in mFile
void clsSolut::AssignsSequenceToCluster()
{
	//assigns ID sequence to a cluster
	for(id_cluster_type iClus = 0; iClus < this->param->nClus; iClus++)
		this->vClusts[iClus]->AssignSequenceToCluster();
}
 
void clsSolut::PrintResult()
{
	cout << endl << "-----------------------------------------" << flush;
	cout << endl << "Number of groups: " << this->vGroup.size();
	cout << endl << "-----------------------------------------" << flush;
	for(id_cluster_type iClus = 0; iClus < this->vClusts.size(); iClus++)
		cout << endl << "Cluster "<< iClus << ": " << this->vClusts[iClus]->GetSize() <<" reads";
	cout << endl << "-----------------------------------------" << endl << flush;
}

void clsSolut::SaveInfo(string elapsedTimeTot, const string& dir_output)
{
	vector<string> param = this->param->GetParametersString();
	param.push_back("");
	info.insert(info.begin(), param.begin(), param.end()); //Inserisci all'inizio
	//Poi tutte le altre info computate in ordine
	info.push_back("Total Computation (second): " + elapsedTimeTot);

	createDirAndSubDir(dir_output);

	//Stampa su file
	string info_file = dir_output + "binning.info";
	ofstream result(info_file, ofstream::out);
	for(size_t i = 0; i < this->info.size(); i++)
		result << info[i] << endl << flush;
	result << endl << flush;

	//Dettagli sui cluster, numero read
	if(this->vClusts.size() > 0)
	{
		result << "Cluster details (id_cls,num_seq):" << endl << flush;
		for(id_cluster_type iClus = 0; iClus < this->vClusts.size(); iClus++)
			result << iClus << "," << this->vClusts[iClus]->GetSize() << endl << flush;
	}
}

void clsSolut::SaveBinning(const string& dir_output)
{
	createDirAndSubDir(dir_output);

	//Prendi header read
	MapIDFile_Header map_idxFile_header;
	this->mFile->GetHeaderFromFile(map_idxFile_header);
	this->SaveGroups(dir_output, map_idxFile_header);
	if(this->vClusts.size() > 0)
		this->SaveClusters(dir_output, map_idxFile_header);
}

void clsSolut::SaveClusters(const string& dir_output, MapIDFile_Header& map_idxFile_header)
{
	vector<PairFiles>& filename = this->mFile->getFilename();
	for(size_t idxFile = 0; idxFile < filename.size(); ++idxFile)
	{
		if(!filename[idxFile].isCorrect())
			continue;
		if(filename[idxFile].getType() == SingleEnd || filename[idxFile].getType() == PairedEnd)
		{
			string cluster_out = dir_output + filename[idxFile].first.getFilename() + ".clusters.csv";
			ofstream result1(cluster_out, ofstream::out);
			auto it_ID = map_idxFile_header.begin();
			while(it_ID != map_idxFile_header.end())
			{
				clsSequencePE::Ptr seq1 = this->mFile->GetAt(filename[idxFile].first.getIdxApp(it_ID->first));
				result1 << it_ID->second.first << "," << seq1->getIdCluster() << endl << flush;
				++it_ID;
			}
		}

		if(filename[idxFile].getType() == PairedEnd)
		{
			string cluster_out = dir_output + filename[idxFile].second.getFilename() + ".clusters.csv";
			ofstream result2(cluster_out, ofstream::out);
			auto it_ID = map_idxFile_header.begin();
			while(it_ID != map_idxFile_header.end())
			{
				clsSequencePE::Ptr seq2 = this->mFile->GetAt(filename[idxFile].second.getIdxApp(it_ID->first));
				result2 << it_ID->second.second << "," << seq2->getIdCluster() << endl << flush;
				++it_ID;
			}
		}
	}
}
void clsSolut::SaveGroups(const string& dir_output, MapIDFile_Header& map_idxFile_header)
{
	vector<PairFiles>& filename = this->mFile->getFilename();
	for(size_t idxFile = 0; idxFile < filename.size(); ++idxFile)
	{
		if(!filename[idxFile].isCorrect())
			continue;
		if(filename[idxFile].getType() == SingleEnd || filename[idxFile].getType() == PairedEnd)
		{
			string gpr_out = dir_output + filename[idxFile].first.getFilename() + ".groups.csv";
			ofstream result1(gpr_out, ofstream::out);
			auto it_ID = map_idxFile_header.begin();
			while(it_ID != map_idxFile_header.end())
			{
				clsSequencePE::Ptr seq1 = this->mFile->GetAt(filename[idxFile].first.getIdxApp(it_ID->first));
				result1 << it_ID->second.first << "," << seq1->getIdGrp() << endl << flush;
				++it_ID;
			}
		}

		if(filename[idxFile].getType() == PairedEnd)
		{
			string gpr_out = dir_output + filename[idxFile].second.getFilename() + ".groups.csv";
			ofstream result2(gpr_out, ofstream::out);
			auto it_ID = map_idxFile_header.begin();
			while(it_ID != map_idxFile_header.end())
			{
				clsSequencePE::Ptr seq2 = this->mFile->GetAt(filename[idxFile].second.getIdxApp(it_ID->first));
				result2 << it_ID->second.second << "," << seq2->getIdGrp() << endl << flush;
				++it_ID;
			}
		}
	}
}

//Reset e libera memoria
void clsSolut::resetComputation(ResetType reset_type)
{
	switch(reset_type)
	{
	case(Group_and_Cluster):
		//Size (quindi oggetto info) si lascia inalterato perché serve a MergeGroup_Fix()
		//Il grafo serve ancora per costruire i gruppi
		this->mFile->resetComputation(clsSequencePE::Group_reset);
		this->mFile->resetComputation(clsSequencePE::Seed_Reset);
		this->mFile->resetComputation(clsSequencePE::Cluster_reset);
		this->vGroup.clear();
		this->vClusts.clear();
		break;
	case(Cluster):
		//Size (quindi oggetto info) si lascia inalterato perché serve a MergeGroup_Fix()
		//Choose non resettato perchè solo perdita di tempo
		//vettore gruppi necessario per computare il merge
		//isSeed necessario
		this->mFile->resetComputation(clsSequencePE::Cluster_reset); //Deve essere ricomputato
		this->mFile->resetComputation(clsSequencePE::Graph_reset); //Non più necessario (Necessario per costruire gruppi)
		this->vClusts.clear();
		break;
	}
}
