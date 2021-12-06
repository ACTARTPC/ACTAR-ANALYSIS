#include <MEvent.h>

MEvent::MEvent()
{
	hasConfigFile=false;
	hasCalBaseline=false;
	isBaselineCorr=true;
	isCalmode=false;
	ImplantationExtract=false;
	
	implantCounter=0;
	
	
	Gsig=new TGraph(NB_SAMPLES);
	ffitSig=new TF1("ffitSig",SignalFitFunction,0,511,2);
	stat_good_frame = stat_bad_frame = stat_deconv_frame = stat_all_frame = 0;
	CoboAsad.resize(65);
	CoboAsad.shrink_to_fit();
  
  INFO.Nevents=0;
  for(int i=0;i<NB_COBO;i++) INFO.FrameCounter[i]=0;
	
	f_out=fopen("check.dat","w");
}


MEvent::~MEvent()
{
	delete Gsig;
	delete ffitSig;
	fclose(f_out);
}



void MEvent::GetAndFillCalibCoeffs()
{	
	int co, as, ag, ch;
	float val[2];
	FILE* fcal;
	TString calfile={"cal/cal_gen_.dat"};

	cout << "No config file found -> try to find cal/cal_gen.dat file..." << endl;
	fcal=fopen(calfile.Data(),"r");
	if(fcal==NULL)
	{
		cout << "Generic calibration file not found. Working with raw data" << endl;
		for(int co=0;co<NB_COBO;co++) for(int as=0;as<NB_ASAD;as++) for(int ag=0;ag<NB_AGET;ag++) for(int ch=0;ch<NB_CHANNEL;ch++)
		{
			CalibCoefs[co][as][ag][ch][0]=1.;
			CalibCoefs[co][as][ag][ch][1]=0.;
		}
	}
	else
	{
 	cout << "Using generic calibration coefficients" << endl;
		while(fscanf(fcal,"%d %d %d %d %f %f",&co,&as,&ag,&ch,&val[0],&val[1])!=EOF)
		{
			if(co<NB_COBO && as<NB_ASAD && ch<NB_CHANNEL)
				for(int i=0;i<2;i++)
					CalibCoefs[co][as][ag][ch][i]=val[i];
			else cout << "calibration file corrupted: cobo/asad/aget/channel exceeds max value of 15/3/3/67" << endl;
		}
		fclose(fcal);
	}
	
	fcal=fopen("cal/BL_.dat","r");
	if(fcal==NULL)
	{
		cout << "Baseline calibration file not found -> No Baseline subtraction." << endl;
		for(int co=0;co<NB_COBO;co++) for(int as=0;as<NB_ASAD;as++) for(int ag=0;ag<NB_AGET;ag++) for(int ch=0;ch<NB_CHANNEL;ch++) for(int bu=0;bu<NB_SAMPLES;bu++)
			BaselineCal[co][as][ag][ch][bu]=0;
	}
	else
	{
		cout << "Baseline subtraction ON from cal/BL.dat" << endl;
		hasCalBaseline=true;
		int co, as, ag, ch;
		while(fscanf(fcal,"%d %d %d %d",&co,&as,&ag,&ch)!=EOF)
		{
			if(co<NB_COBO && as<NB_ASAD && ag<NB_AGET && ch<NB_CHANNEL)
			{
				for(int bu=0;bu<NB_SAMPLES;bu++)
					fscanf(fcal,"%f",&BaselineCal[co][as][ag][ch][bu]);
			}
			else cout << "Wrong global channel number for baseline calibration: " << co << " " << as << " " << ag << " " << ch << endl;
		}
		fclose(fcal);
	}
}



void MEvent::InitVXIParameters(char* ActionFName,std::vector<TString> ParNames,std::vector<int> Numbers)
{
	TString PName;
	parametersVXI.FillFromActionFile(ActionFName);
	auto ind = parametersVXI.GetLabel("Hello");
	
	for(int i=0;i<16384;i++)
		labelVXI[i]=-1;
	
	for(int it=0;it<ParNames.size();it++)
	{
		if(Numbers[it])
		{
			for(int noc=0;noc<Numbers[it];noc++)
			{
       				PName.Form("%s%d",ParNames[it].Data(),noc);
				
				
				if(parametersVXI.GetLabel(std::string(PName.Data()))>=0 && parametersVXI.GetLabel(std::string(PName.Data()))<16384)
          			labelVXI[parametersVXI.GetLabel(std::string(PName.Data()))]=it*1000+noc;
			}
		}
		else
    		{
			if(parametersVXI.GetLabel(std::string(ParNames[it].Data()))>=0 && parametersVXI.GetLabel(std::string(ParNames[it].Data()))<16384)
				labelVXI[parametersVXI.GetLabel(std::string(ParNames[it].Data()))]=it*1000;
		}
	}
}


void MEvent::SetSpecificTreatment()
{
	int co, as, ag, ch;
	FILE* fcal;
	for(int co=0;co<NB_COBO;co++) for(int as=0;as<NB_ASAD;as++) for(int ag=0;ag<NB_AGET;ag++) for(int ch=0;ch<NB_CHANNEL;ch++)
		SpecificTreatment[co][as][ag][ch]=0;
	fcal=fopen(SpecificTreatmentFile,"r");
	int algonum;
	cout << "Filling specific treatment arguments from " << SpecificTreatmentFile << endl;
	if(fcal==NULL)
	{
		hasSpecificTreatment=false;
		cout << "Warning: specific treatment file : " << SpecificTreatmentFile << " does not exist" << endl;
	}
	else while(fscanf(fcal,"%d %d %d %d %d",&co,&as,&ag,&ch, &algonum)!=EOF)
		SpecificTreatment[co][as][ag][ch]=algonum;
	if(fcal!=NULL)
		fclose(fcal);
}


void MEvent::GetAndFillRespFunc()
{
	const int NCO=16;

	const char* AsAdList[NCO*NB_ASAD]={"0011114D","0011114C","001110FE","001110FF","001110F1","001110F8","0011111D","00111120","00111C96","00111C9C","001110E1","0011111A","00111119","001110DB","001110DF","001110DC","001110E0","001110DE","001110EA","00110D1E","00110D10","00110D0C","00110D04","00110D04","00110D07","00110D0F","00110D0B","00110D08","00110D00","00110D03","00110D02","00110D12","00110D01","00110D06","00110D0A","001110F4","001110FC","001110ED","001110EE","001110FA","001110F9","001110F2","001110F3","001110EB","001110F7","001110F6","0011111F","00111113","00111115","001110E9","00111121","001110E3","001110DA","00111116","00110D0D","00110D09","001110DD","001110E4","00110C97","00110C98","00110C99","00110C9A","0011011B","001100E2"};
	
	TString Path="/home/roger/Desktop/ANALYSIS/E690/CALIB_ASAD/";
	TString Supl="/";
// 	TString Path="/media/roger/ACTAR_USB_2/E690_CALIB/";
// 	TString Supl="/";
	
	for(int co=0;co<NCO;co++) for(int as=0;as<NB_ASAD;as++)
	{
		TString RootFileName=Path+AsAdList[co*NB_ASAD+as]+Supl+"ResponseFFT.root";
		TFile* fFunc=new TFile(RootFileName.Data(),"read");
		for(int ag=0;ag<NB_AGET;ag++)
		{
			for(int ch=0;ch<NB_CHANNEL;ch++)
			{
				TString TH1Name;
				TH1Name.Form("%s_G%d_Ch%02d",AsAdList[co*NB_ASAD+as],ag,ch);

				TH1D* hRespF=(TH1D*)(fFunc->Get(TH1Name.Data()));
				for(int Bu=0;Bu<NB_SAMPLES;Bu++)
					ResponseFunction[co][as][ag][ch][Bu]=hRespF->GetBinContent(Bu+1);
			}
		}
		
		fFunc->Close();
		delete fFunc;
	}
}


