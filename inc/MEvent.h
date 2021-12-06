#ifndef MEVENT_H
#define MEVENT_H

#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string.h>
#include <cstdlib>
#include <complex>
#include <TObject.h>
#include <TString.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TFile.h>


#include <MCoboAsad.h>
#include <MEventReduced.h>
#include <Parameters.h>
#include <DataParameters.h>

#include "GFFTDianelsonLanczos.hh"

using namespace std;


struct RunInfo{
  int Nevents;
  int FrameCounter[NB_COBO];
};

class MEvent
{
	public:
	int run_number;
	unsigned long int EN;
	unsigned long int TS;
	std::vector<MCoboAsad> CoboAsad;
  
  RunInfo INFO;

	
	bool hasConfigFile;
	float CalibCoefs[NB_COBO][NB_ASAD][NB_AGET][NB_CHANNEL][2];
	bool hasCalBaseline;
	bool isBaselineCorr;
	float BaselineCal[NB_COBO][NB_ASAD][NB_AGET][NB_CHANNEL][NB_SAMPLES];
	bool hasSpecificTreatment;
	char* SpecificTreatmentFile;
	int SpecificTreatment[NB_COBO][NB_ASAD][NB_AGET][NB_CHANNEL];
	
	bool hasResponseFunction;
	double ResponseFunction[NB_COBO][NB_ASAD][NB_AGET][NB_CHANNEL][NB_SAMPLES];
	bool hasDeconv;
	bool isDeconvOK;
	bool isCalmode;
	int DeconvCondLabel;
	int DeconvCondValue;
  
	bool ImplantationExtract;
	unsigned long implantCounter;
	int ImplantCondLabel1;
	int ImplantCondLabel2;
	int ImplantCondValueLow1;
	int ImplantCondValueHigh1;
	int ImplantCondValueLow2;
	int ImplantCondValueHigh2;
	
	unsigned long stat_good_frame;
	unsigned long stat_recovered_frame;
	unsigned long stat_bad_frame;	
	unsigned long stat_deconv_frame;
	unsigned long stat_all_frame;

	
	MEvent();
	~MEvent();
	void GetAndFillCalibCoeffs();
	void SetSpecificTreatment();
	void GetAndFillRespFunc();
	void InitVXIParameters(char*,std::vector<TString>,std::vector<int>);
	int TreatBaseline(bool,int);
	int TreatFullSignal(bool);
   
	
	MEventReduced ReducedEvent;
	short labelVXI[16384];
  
  	FILE* f_out;
	
	protected:
	DataParameters parametersVXI;
// 	short labelVXI[16384];
	
	TF1* ffitSig;
	TGraph* Gsig;
	float Charge[NB_SAMPLES];
  

};

double SignalFitFunction(double*, double*);


#endif
