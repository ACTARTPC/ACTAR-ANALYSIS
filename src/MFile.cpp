///////////////////////////////////////
//                                   //
// T. Roger - GANIL 2015             //
// email: roger@ganil.fr             //
//                                   //
// MFile class file:                 //
//     - file opening & closing      //
//     - frame buffering             //
//                                   //
///////////////////////////////////////

#include <MFile.h>
#include <iostream>

using namespace std;

MFile::MFile(TString path, int runf, int runl)
{
#ifdef CC_IN2P3
	LSbase="ils ";
#endif
#ifndef CC_IN2P3
	LSbase="ls -1 -v ";
#endif

	runpath=path;
	for(int i=runf;i<=runl;i++)
	{
		char runnumber[256];
		sprintf(runnumber,"run_%04d",i);
		TString command=LSbase + runpath + " |grep " + runnumber;
		TString rawLS = gSystem->GetFromPipe(command.Data());
		TObjArray* orderedLS = rawLS.Tokenize("\n");
		//orderedLS->Sort();
		for(int j=0;j<orderedLS->GetEntries();j++)
 		{
			if(strstr(((TObjString*)(*orderedLS)[j])->GetName(),runnumber)==((TObjString*)(*orderedLS)[j])->GetName())
			{
				List.push_back(path);
				List.back()+=((TObjString*)(*orderedLS)[j])->GetName();
			}
		}
		delete orderedLS;
	}
	
	
	ACTIONpath=runpath;
	
	EOFreached=false;
	isSingleFile=false;
	
	int minsizeheader = MFM_BLOB_HEADER_SIZE; // =8	
	vectorsize = minsizeheader;
	vector = (char*) (malloc(vectorsize));
	pvector = &vector;

	nFrames=0;
}


MFile::MFile(TString run_name)
{

	List.push_back(run_name);
	
	EOFreached=false;
	isSingleFile=true;
	
	int minsizeheader = MFM_BLOB_HEADER_SIZE; // =8	
	vectorsize = minsizeheader;
	vector = (char*) (malloc(vectorsize));
	pvector = &vector;
}


MFile::~MFile()
{
	free(vector);
	List.clear();
}


long int MFile::OpenNext(TString run_file)
{
	fLun = open(run_file.Data(),O_RDONLY);

	if(fLun<=0 && run_file.Length()>1)
	{
		cout << "Error open file: " << run_file.Data() << endl;
		return(-1);
	}
	
	else if(run_file.Length()<=1)
	{
		return(0);
	}
	
	else if(run_file.Length()>1)
	{
		char *pch;
  		pch =strstr ((char *)(run_file.Data()),"/run_");
		
		if(!isSingleFile) run_number = strtol(pch+5,NULL,10);
		else run_number = 0;
		
		if(run_number!=prev_run_number)
		{
			if(prev_run_number>=0) cout << "------------------------------" << endl;
			isNewRun=true;
			
			char runnumber[256];
			sprintf(runnumber,"run_%04ld",run_number);
			
			prev_run_number=run_number;
			
			ACTIONfileName = LSCommand(run_number, ACTIONpath, "ACTION");			
		}
		
		else isNewRun=false;
		
		cout << "Treating run " << run_number << " : " << run_file.Data() << endl;

		EOFreached=false;
		return(run_number);
	}
	else return(0);
}



TString MFile::LSCommand(int run, TString basepath, TString filetype)
{
	char runnumber[256];
	sprintf(runnumber,"run_%04d",run);
	
	TString command;
	TString rawLS;
	TObjArray* orderedLS;
	TString fileName;
	
	command=LSbase + basepath;
	rawLS = gSystem->GetFromPipe(command.Data());
	orderedLS = rawLS.Tokenize("\n");
	orderedLS->Sort();
	for(int j=0;j<orderedLS->GetEntries();j++)
 	{
		if(((TObjString*)(*orderedLS)[j])->GetString().Contains(filetype) && ((TObjString*)(*orderedLS)[j])->GetString().Contains(runnumber))
		{
			cout << "Found " << filetype << " file: " << ((TObjString*)(*orderedLS)[j])->GetName() << endl;
			fileName = basepath + ((TObjString*)(*orderedLS)[j])->GetName();
		}
	}
	delete orderedLS;
	return(fileName);
}


int MFile::GetNextBuffer(MFMCommonFrame * frame)
{
	int framesize = frame->ReadInFile(&fLun, pvector, &vectorsize);
	nFrames++;
	
	if(framesize<=0) EOFreached=true;
	if(EOFreached) cout << "Analyzed " << nFrames << " Frames" << endl;
	return(framesize);
}

int MFile::GetPrevBuffer(MFMCommonFrame * frame, int prev_framesize)
{
	lseek(fLun,-prev_framesize,SEEK_CUR);
	int framesize = frame->ReadInFile(&fLun, pvector, &vectorsize);
		
	if(framesize<=0) EOFreached=true;
	return(framesize);
}
