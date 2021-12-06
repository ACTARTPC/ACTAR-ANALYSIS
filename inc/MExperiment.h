#ifndef MEXPERIMENT_H
#define MEXPERIMENT_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string.h>
#include <TObject.h>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifndef WIN32
#include <unistd.h>
#endif
#ifdef WIN32
#define stat _stat
#endif


using namespace std;

#include <MEvent.h>
#include <MFile.h>
#include <MUnpacker.h>
#include <MVisu.h>
#include <Parameters.h>
#include <limits>

class MExperiment
{
	public:
	MExperiment(int, char**);
	~MExperiment();
	
	int TreatRuns();
	void NewRunAction(long);
	bool TreatEvent(MFMCommonFrame *);
	void OpenNewTree(long);
	void CloseTree();
	void EndOfExperiment();
  
	void SetPartialConversionList();

	int Experiment_Number;
	TString Experiment_Name;
	TString Base_Path;
	TString Tree_Base_Path;
	TString SingleFileName;

	MEvent* Event;
	MFile* File;
	MUnpacker* Unpacker;
	MVisu* Visu;
	TTree* Tree;
	TFile* TreeFile;
	
	int runf;
	int runl;	

	bool VisuOpt;
	bool TreeOpt;
	bool FullTreeOpt;
	bool SingleFile;
  bool isOnline;
	int IsFastPeak;
	
	bool ImplantationExtract;
	
	int Run_number_limit;
	int Cur_run_number;

	long int SplitRuns;

	bool isRawMode;
	int fLunRaw;
	
	bool PartialConversion;
	char* PartialConversionBaseFile;
	TString PartialConversionFile;
	std::vector<long int> PartialConversionList;
	
	bool isTSmode;
	
	int NFragments;
	unsigned long int* last_EN;
	unsigned long int* last_TS;

	//for the Time measurement
	struct timeval  tv1, tv2, tv3, tv4;
	double tdiff1, tdiff2, tdiff3, tdiff4;

	
};



#endif
