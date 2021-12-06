#include <MExperiment.h>
#include <iostream>

using namespace std;


MExperiment::MExperiment(int argc,char **argv)
{
	if(atoi(argv[1])) runf=atoi(argv[1]);
	else
	{
		cout << "usage: >Analyse.exe first_run [last_run] [-options]" << endl;
		exit(1);
	}
	if(atoi(argv[2])) runl=atoi(argv[2]);
	else runl=runf;

	Base_Path="/adata/eactar/e791/Full/";
	
	VisuOpt=false;
	TreeOpt=false;
	FullTreeOpt=false;
	SingleFile=false;
	IsFastPeak=1;
	isTSmode=false;
	isRawMode=false;
	PartialConversion=false;
	ImplantationExtract=false;
  isOnline=false;
	SplitRuns=-1;
	Cur_run_number=0;
	Run_number_limit=std::numeric_limits<int>::max();
	
	Event=new MEvent();
	Unpacker=new MUnpacker();

	Event->hasSpecificTreatment=false;
	
	int opt;// = getopt(argc, argv, "P:V:T:F");
	while ((opt = getopt(argc, argv, "VT:FP:f:S:B:p:bL:M:CIRO"))!=-1)
	{
		switch(opt)
		{ 
			case 'V':                                        // visu mode
				VisuOpt=true;
				if(Visu==NULL) Visu=new MVisu(argc,argv);
				// Visu->SetListToDraw(optarg);
				break;
				
			case 'P':                                        // data path change
				if(optarg[0]=='-')
				{
					cout << "P option requires path name" << endl;
					exit(1);
				}
				if(optarg[0]=='+')
				{
					for(int cc=1;cc<sizeof(optarg)-1;cc++)
						Base_Path+=optarg[cc];
					char* pch=strstr(optarg,"CoBo");
					if(pch!=NULL) SplitRuns=strtol(pch+4,NULL,10);
				}
				else Base_Path=optarg;
				break;
				
			case 'T':                                        // tree writing mode
				if(!TreeOpt) TreeOpt=true;
				Tree=NULL;
				Tree_Base_Path=optarg;
				break;
				
			case 'F':                                        // full tree mode: all samples
				if(TreeOpt) FullTreeOpt=true;
				break;
			
			case 'f':                                        // single file treatment (for debug mainly)
				SingleFileName=optarg;
				SingleFile=true;
				break;
				
			case 'S':                                        // "S" fitted charge treatment -Much slower or "S+" for integral above threshold -Faster
				if(optarg[0]=='+')
				{
					IsFastPeak=2;
					cout << "Integral charge treatment method ON" << endl;
				}
				else
				{
					IsFastPeak=0;
					cout << "Fitted charge treatment option ON" << endl;
				}
				break;

			case 'B':                                        // specific algo (zero suppress) file
				Event->hasSpecificTreatment=true;
				Event->SpecificTreatmentFile=optarg;
				cout << "Specific treatment ON " << endl;
				Event->SetSpecificTreatment();
				break;

			case 'O':                                        // specific algo (zero suppress) file
        if(runf!=runl)
        {
          cout << "ERROR: ONLINE option must be used for a single run number" << endl;
          cout << "---> ONLINE option disactivated!!!" << endl;
        }
        else isOnline=true;
				break;

			case 'b':                                        // specific algo (zero suppress) file
				Event->hasDeconv=true;
				Event->DeconvCondLabel=30;
				Event->DeconvCondValue=450;
				Event->hasResponseFunction=true;
				Event->GetAndFillRespFunc();
				cout << "Specific deconvolution ON for label=" << Event->DeconvCondLabel << " with value" << Event->DeconvCondValue << endl;
				break;

			case 'p':                                        // partial conversion (event number) - works well for 1 run & list not too long...
				PartialConversion=true;
				Event->GetAndFillRespFunc();
				PartialConversionBaseFile=optarg;
				cout << "Partial conversion ON" << endl;
				//SetPartialConversionList();
				break;
			
			case 'L':
				Run_number_limit=atoi(optarg);
				cout << "Limiting analysis to " << Run_number_limit << " run files" << endl;
				break;
				
			case 'M':                                        // partial conversion (event number) - works well for 1 run & list not too long...
				NFragments=atoi(optarg);
				isTSmode=true;
				last_EN=(unsigned long int*)malloc(NFragments*sizeof(unsigned long int));
				last_TS=(unsigned long int*)malloc(NFragments*sizeof(unsigned long int));
				cout << "TimeStamp mode ON: event number relies on this program..." << endl;
				break;
			
			case 'C':                                        //Calib mode: forces deconvolution in case of unmerged data
				Event->isCalmode=true;
				cout << "Calib mode ON" << endl;
				break;
			
			case 'I':
 				ImplantationExtract=true;
				Event->ImplantCondLabel1=28;
				Event->ImplantCondLabel2=30;
				Event->ImplantCondValueLow1=2700;
				Event->ImplantCondValueHigh1=4000;
				Event->ImplantCondValueLow2=5950;
				Event->ImplantCondValueHigh2=6400;
				break;
				
			case 'R':                                        //Calib mode: forces deconvolution in case of unmerged data
				if(!PartialConversion && !ImplantationExtract)
				{
					cout << "Raw data mode should be used with partial conversion mode: -p event.list" << endl;
					break;
				}
				isRawMode=true;
				cout << "Writing Raw data in a file" << endl;
				break;
			
			
			default:
				exit(1);
		}
	}
	
	
	cout << "Analyzing run(s): " << runf << " " << runl << " in path " << Base_Path << endl;
	
	
	if(!SingleFile) File=new MFile(Base_Path,runf,runl);
	else File=new MFile(SingleFileName);
}


MExperiment::~MExperiment()
{
	delete Unpacker;
//	delete Event;
	delete File;
}


void MExperiment::OpenNewTree(long run_number)
{
	TString TreeName;
	TreeName.Form("Tree_Run_%04ld",run_number);

	TString TreeNameFile;
	if(FullTreeOpt) cout << "Full ";
	TreeNameFile.Form("%s%s.root",Tree_Base_Path.Data(),TreeName.Data());
	cout << "TTree " << TreeNameFile.Data() << " will be written..." << endl;
	TreeFile=new TFile(TreeNameFile.Data(),"recreate");

	Tree=new TTree("ACTAR_TTree","1st level Tree",0);
	Tree->Branch("data","MEventReduced",&Event->ReducedEvent,16000,0);
	Tree->Branch("TimeStamp",&Event->ReducedEvent.timestamp,"TimeStamp/L");
}

void MExperiment::CloseTree()
{
	if(Tree!=NULL)
	{
		Tree->Write();
		TreeFile->Close();
		printf("\nWriting tree in file of size = %.1f Mo\n",TreeFile->GetSize()/1024./1024.);
		if(Tree!=NULL)
		{
// 			delete Tree;
			Tree=NULL;
		}
		if(TreeFile!=NULL)
		{
			delete TreeFile;
			TreeFile=NULL;
		}
	}
}

int MExperiment::TreatRuns()
{
	cout << "Treating " << File->List.size() << " Files" << endl;
  cout << "Options Are:" << endl;
  cout << "TreeOpt: " << TreeOpt << endl;
  cout << "FullTreeOpt: " << FullTreeOpt << endl;
  cout << "SingleFile: " << SingleFile << endl;
  cout << "IsFastPeak: " << IsFastPeak << endl;
  cout << "isRawMode: " << isRawMode << endl;
  cout << "PartialConversion: " << PartialConversion << endl;
  cout << "ImplantationExtract: " << ImplantationExtract << endl;
  cout << "isOnline: " << isOnline << endl;
	
  if(!isOnline)
  {
    for(TString run_file :File->List) if(strcmp(run_file.Data(),"/adata/eactar/e791/Full/run_0276.dat.19-05-21_18h01m51s.27"))
	  {
		  long run_number=File->OpenNext(run_file);
		  Cur_run_number++;
		  if(Cur_run_number>Run_number_limit) return(0);
    
		  if(File->isNewRun)
        NewRunAction(run_number);
		
		  long int total_size=0;
		  MFMCommonFrame * frame = new MFMCommonFrame();

		  do
		  {
			  bool isOK=false;
			  long int fsize=0;
			  fsize =  File->GetNextBuffer(frame);
			  total_size+= fsize;
			  isOK=TreatEvent(frame);

			  if(total_size/(1024*1024)%100==0)
			  {
 				  cout << "total read size = " << total_size/(1024*1024) << " Mo      Good data: " << Event->stat_good_frame << " Corrupted data stat: " << Event->stat_bad_frame << " Recovered data stat: " << Event->stat_recovered_frame << "     Deconv frame: " << Event->stat_deconv_frame*100./Event->stat_all_frame << "%                            \r";
  			  cout.flush();
			  }
		  }
		  while(!File->EOFreached);
    
		  delete frame;
	  }
  }
  
  
  else
  {

  	TString LSbase;
#ifdef CC_IN2P3
	  LSbase="ils ";
#endif
#ifndef CC_IN2P3
	  LSbase="ls -1 -v ";
#endif
    
    NewRunAction(runf);
    
    bool end=false;
    bool hasAnalyzed=false;
    struct stat result;
    
    int cur_list_pos=0;
    
    char runnumber[256];
    sprintf(runnumber,"run_%04d",runf);
    TObjArray* orderedLS_new;
    TString command=LSbase + Base_Path + "|grep -v -e \"ACTIONS\" -e \"configure\" -e \"describe\" |grep " + runnumber;
    while(!end)
    {
      TString rawLS = gSystem->GetFromPipe(command.Data());
      TObjArray* orderedLS = rawLS.Tokenize("\n");

      if(orderedLS->GetEntries()==0) usleep(ONLINE_WAIT_NEWRUN*1.E6);
      else
      {
        gettimeofday(&tv1, NULL);
        time_t mod_time0;
        time_t mod_time;
        auto Filename =strstr(((TObjString*)(*orderedLS)[cur_list_pos])->GetName(),runnumber);
        TString filename=Base_Path+Filename;
//        cout << "filename = " << filename << " at pos " << cur_list_pos << " / " << orderedLS->GetEntries() << endl;
        do
        {
          if(stat(filename.Data(), &result)==0) mod_time0 = result.st_mtime;
          do
          {
            usleep(ONLINE_WAIT_NEWRUN*1.E6/4.);
            if(stat(filename.Data(), &result)==0) mod_time = result.st_mtime;
            cout << "File " << filename.Data() << " currently opened (" << result.st_size/1000./1000. << " Mo): Wait " << ONLINE_WAIT_NEWRUN/40. << " seconds                 \r";
            cout.flush();
            gettimeofday(&tv2, NULL);
            //cout << " mod_time= " << mod_time << "   mod_time0= " << mod_time0 << "    "  << "times are: " << tv2.tv_sec << " " << tv1.tv_sec << endl;
          }
          while(mod_time!=mod_time0 && (tv2.tv_sec-tv1.tv_sec)<ONLINE_WAIT_NEWRUN);
          //cout << "OUT: mod_time= " << mod_time << "   mod_time0= " << mod_time0 << endl;
        }
        while(mod_time!=mod_time0);
        
        cout << endl;
        cout << "Starting Run Analysis..." << endl;
        long run_number=File->OpenNext(filename.Data());
        long int total_size=0;
        MFMCommonFrame * frame = new MFMCommonFrame();
        do
        {
          bool isOK=false;
          long int fsize=0;
          fsize =  File->GetNextBuffer(frame);
          total_size+= fsize;
          isOK=TreatEvent(frame);
          if(total_size/(1024*1024)%100==0)
          {
            cout << "total read size = " << total_size/(1024*1024) << " Mo      Good data: " << Event->stat_good_frame << " Corrupted data stat: " << Event->stat_bad_frame << " Recovered data stat: " << Event->stat_recovered_frame << "     Deconv frame: " << Event->stat_deconv_frame*100./Event->stat_all_frame << "%                            \r";
            cout.flush();
          }
        }
        while(!File->EOFreached);
        cur_list_pos++;
        delete frame;
      }
    
      rawLS = gSystem->GetFromPipe(command.Data());
      TObjArray* orderedLS_new = rawLS.Tokenize("\n");
      if(orderedLS_new->GetEntries()==orderedLS->GetEntries() && orderedLS->GetEntries()!=0 && cur_list_pos==orderedLS->GetEntries()) end=true; 
    }
	}
  
  
	return(1);
}


bool MExperiment::TreatEvent(MFMCommonFrame * frame)
{
	Event->CoboAsad.clear();
	if(Event->isCalmode) Event->isDeconvOK=true;
	else Event->isDeconvOK=false;
	Unpacker->Unpack(frame,Event);

	if(isTSmode)
	{
		int global_fragment_number;
		if(Event->CoboAsad[0].cobo_number!=31 && Event->CoboAsad[0].global_asad_number<NFragments)  global_fragment_number=Event->CoboAsad[0].global_asad_number;
		else if(Event->CoboAsad[0].cobo_number==31 && NFragments%4==1) global_fragment_number=NFragments-1;
		else global_fragment_number=-1;
		if(global_fragment_number>=0)
		{
			if(last_TS[global_fragment_number]<Event->TS)
			{
				last_TS[global_fragment_number]=Event->TS;
				Event->EN=last_EN[global_fragment_number];
				last_EN[global_fragment_number]++;
			}
			else cout << "WARNING: Event Number not well defined in TimeStamp mode: cur_TS = " << Event->TS << " and last_TS = " << last_TS[global_fragment_number] << " For CoBo " << Event->CoboAsad[0].cobo_number << " and AsAd " << Event->CoboAsad[0].asad_number << endl;
		}
	}
	
	if(!PartialConversion && !ImplantationExtract)
	{//Event->isDeconvOK=true;
		if(Event->hasDeconv && Event->isDeconvOK)
		{
			Event->TreatFullSignal(FullTreeOpt);
			Event->isDeconvOK=false;
			Event->stat_deconv_frame++;
			Event->stat_all_frame++;
		}
		else
		{
			Event->TreatBaseline(FullTreeOpt,IsFastPeak);
			Event->stat_all_frame++;
		}
		if(VisuOpt) Visu->Draw(Event,100);
		if(TreeOpt) Tree->Fill();
	}
	else if(ImplantationExtract)
	{
		Event->stat_all_frame++;
		if(Event->ImplantationExtract)
		{
			int checksize = write(fLunRaw, frame->GetPointHeader(), frame->GetFrameSize());
			if(checksize!=frame->GetFrameSize())
			{
				cout << "Problems while writing raw file" << endl;
				exit(1);
			}
			Event->ImplantationExtract=false;
			Event->stat_deconv_frame++;
		}		
	}
	else
	{
		bool isInList=false;
		for(long int NEV : PartialConversionList)
		{
			if(Event->EN==NEV)
			{
				isInList=true;
				break;
			}
		}
		if(isInList)
		{	
			if(isRawMode)
			{
				int checksize = write(fLunRaw, frame->GetPointHeader(), frame->GetFrameSize());
				if(checksize!=frame->GetFrameSize())
				{
					cout << "Problems while writing raw file" << endl;
					exit(1);
				}
			}
			else
			{
				if(Event->hasResponseFunction)
				{
					Event->TreatFullSignal(FullTreeOpt);
					Event->stat_deconv_frame++;
					Event->stat_all_frame++;
				}
				else
				{
					Event->TreatBaseline(FullTreeOpt,IsFastPeak);
					Event->stat_all_frame++;
				}
				if(TreeOpt) Tree->Fill();
			}
		}
	}
	Event->ReducedEvent.CoboAsad.clear();


	return(true);
}


void MExperiment::EndOfExperiment()
{
	if(TreeOpt)
	{
		Tree->Write();
		TreeFile->Close();
		printf("\nWriting tree in file of size = %.1f Mo\n",TreeFile->GetSize()/1024./1024.);
	}
	if(VisuOpt)
	{
		Visu->End();
	}
  
  cout << "Frame counts:" << endl;
  for(int i=0;i<NB_COBO;i++) cout << i << "\t";
  cout << endl;
  for(int i=0;i<NB_COBO;i++) cout << Event->INFO.FrameCounter[i] << "\t";
  cout << endl;
}


void MExperiment::SetPartialConversionList()
{
	PartialConversionList.clear();
	FILE* flist=fopen(PartialConversionFile.Data(),"r");
	if(flist==NULL)
	{
		cout << "Partial conversion file " << PartialConversionFile.Data() << " does not exist... STOPPING" << endl;
		exit(0);
	}
	long int tmp_EN;
	while(fscanf(flist,"%ld",&tmp_EN)!=EOF)
		PartialConversionList.push_back(tmp_EN);

	cout << "found " << PartialConversionList.size() << " events to convert from: " << PartialConversionFile << endl;
	
	fclose(flist);
}



void MExperiment::NewRunAction(long int run_number)
{
	if((PartialConversion && isRawMode) || ImplantationExtract)
	{
		PartialConversionFile.Form("%s%04ld.list",PartialConversionBaseFile,run_number);
		SetPartialConversionList();
		if(fLunRaw>0)
		{
			close(fLunRaw);
			fLunRaw=-1;
		}
		if(fLunRaw<=0)
		{
			fLunRaw=open(Form("raw/run_%04ld-part.dat",run_number),(O_RDWR | O_CREAT | O_TRUNC),0644);
			cout << "Opening raw reduced output file: " << Form("raw/run_%04ld-part.dat",run_number) << " with logic unit: " << fLunRaw << endl;
		}
	}
	if(TreeOpt) CloseTree();
	if(TreeOpt) OpenNewTree(run_number);
	Event->GetAndFillCalibCoeffs();
	Event->run_number=run_number;
}
