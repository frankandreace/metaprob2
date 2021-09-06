#ifndef SRC_INPUTCLASS_H_
#define SRC_INPUTCLASS_H_

#include "../Utilities.h"

enum FileType {Fasta, Fastq, UnknownFile};
enum TypePair {SingleEnd, PairedEnd, UnknownPairState};

struct SingleEndFile
{
	typedef shared_ptr<SingleEndFile> Ptr;
	typedef unordered_map<id_seq_type, id_seq_type> Map_IdxFile_IdxApp;
	SingleEndFile();
	bool init(string path);
	bool isCorrect() const;
	void AddReadMap(id_seq_type idxFile, id_seq_type idxApp); //idFile is ordered read sequence number in file
	//Get
	const string& getPath() const;
	const vector<string>& getPathParse() const;
	const string getDirectory() const;
	const string& getFilename() const;
	const string getFilenameWithoutExt() const;
	const string getExt() const;
	FileType getFileType() const;
	const char& getReadDelimiter() const;
	id_seq_type getIdxApp(id_seq_type idxFile);
	void reset();
	Map_IdxFile_IdxApp::iterator beginMapID();
	Map_IdxFile_IdxApp::iterator endMapID();

private:
	string path;
	vector<string> path_parse;
	FileType file_type;
	char read_delimiter;
	bool open_correct;
	Map_IdxFile_IdxApp map_idxF_idxA;

	void parsePath();
};

struct PairFiles : public pair<SingleEndFile, SingleEndFile>
{
	PairFiles(string path_end1, string path_end2);
	bool isCorrect(); //return if the file or files opened isCorrect
	TypePair getType() const;
	const string& getIdentify() const;

private:
	TypePair pair_type;
	string identify;
};

#endif /* SRC_INPUTCLASS_H_ */
