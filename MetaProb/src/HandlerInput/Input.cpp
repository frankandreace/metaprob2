/*
 * InputClass.cpp
 *
 *  Created on: 10/apr/2015
 *      Author: Samuele Girotto
 */

#include "Input.h"

PairFiles::PairFiles(string path_end1, string path_end2)
{
	bool init1 = false;
	bool init2 = false;
	if(!path_end1.empty())
	{
		init1 = get<0>(*this).init(path_end1);
		if(init1)
		{
			this->pair_type = SingleEnd;
			if(!path_end2.empty())
			{
				init2 = get<1>(*this).init(path_end2);
				if(init2)
					this->pair_type = PairedEnd;
			}
		}
		else
		{
			if(!path_end2.empty())
			{
				init2 = get<0>(*this).init(path_end2);
				if(init2)
					this->pair_type = SingleEnd;
			}
		}
		if(init1 || init2)
			this->identify = get<0>(*this).getFilenameWithoutExt();
	}
}

bool PairFiles::isCorrect() {
	if(this->pair_type == SingleEnd)
		return get<0>(*this).isCorrect();
	else if(this->pair_type == PairedEnd)
		return get<0>(*this).isCorrect() && get<1>(*this).isCorrect();
	else
		return false;
}

TypePair PairFiles::getType() const {
	return pair_type;
}

SingleEndFile::SingleEndFile()
{
	this->reset();
}

bool SingleEndFile::init(string path)
{
	this->reset();
	string line; //store the read data
	ifstream file(path);
	if(file.is_open())
	{
		//Prendi prima riga e guarda primo elemento per capire la tipologia di file
		getline(file,line);
		if(!line.empty())
		{
			this->path = path;
			this->parsePath();
			this->open_correct = true;
			if(line[0] == '>')
			{
				this->file_type = FileType::Fasta;
				this->read_delimiter = '>';
			}
			if(line[0] == '@')
			{
				this->file_type = FileType::Fastq;
				this->read_delimiter = '@';
			}
			file.close();
			return true;
		}
		else
		{
			file.close();
			return false;
		}
	}
	else
	{
		if(!path.empty())
			cerr << endl << "Fail to open: " << path << flush;
		return false;
	}
}

bool SingleEndFile::isCorrect() const {
	return this->open_correct && this->file_type != UnknownFile;
}

void SingleEndFile::AddReadMap(id_seq_type idxFile, id_seq_type idxApp) {
	this->map_idxF_idxA[idxFile] = idxApp; //Accedo direttamente perché devo solo inserire
}

const string& SingleEndFile::getPath() const {
	return path;
}

const vector<string>& SingleEndFile::getPathParse() const {
	return path_parse;
}

const string SingleEndFile::getDirectory() const
{
	string dir = "";
	for(size_t i = 0; i < this->path_parse.size() - 1; ++i)
		dir += this->path_parse[i] + "/";
	return dir;
}

const string& SingleEndFile::getFilename() const {
	return this->path_parse.back();
}

const string SingleEndFile::getFilenameWithoutExt() const {
	string filename_without_ext = this->path_parse.back();
	string delimiter = ".";
	size_t pos = filename_without_ext.rfind(delimiter);
	filename_without_ext.erase(pos, filename_without_ext.length());
	return filename_without_ext;
}

const string SingleEndFile::getExt() const {
	string filename_without_ext = this->path_parse.back();
	string delimiter = ".";
	size_t pos = filename_without_ext.rfind(delimiter);
	return filename_without_ext.substr(pos, filename_without_ext.length());
}

FileType SingleEndFile::getFileType() const {
	return file_type;
}

const string& PairFiles::getIdentify() const {
	return identify;
}

const char& SingleEndFile::getReadDelimiter() const {
	return read_delimiter;
}

id_seq_type SingleEndFile::getIdxApp(id_seq_type idxFile) {
	//Uso find perchè se non inserito non è presente e non voglio occupar memoria per niente
	auto it = this->map_idxF_idxA.find(idxFile);
	if(it != this->map_idxF_idxA.end()) //trovato
		return it->second;
	else
		return 0;
}

void SingleEndFile::reset()
{
	this->path = "";
	vector<string>().swap(this->path_parse);
	this->file_type = UnknownFile;
	this->read_delimiter = '>';
	this->open_correct = false;
	Map_IdxFile_IdxApp().swap(this->map_idxF_idxA);
}

SingleEndFile::Map_IdxFile_IdxApp::iterator SingleEndFile::beginMapID() {
	return this->map_idxF_idxA.begin();
}

SingleEndFile::Map_IdxFile_IdxApp::iterator SingleEndFile::endMapID() {
	return this->map_idxF_idxA.end();
}

void SingleEndFile::parsePath()
{
	string tmp = path;
	string delimiter = "/";
	size_t pos = tmp.find(delimiter);
	while(pos != tmp.npos)
	{
		this->path_parse.push_back(tmp.substr(0, pos));
		tmp.erase(0, pos + delimiter.length());
		pos = tmp.find(delimiter);
	}
	this->path_parse.push_back(tmp);
}
