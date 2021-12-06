#include <TFile.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TSpectrum.h>

#include <iostream>

using namespace std;

int comp (const void *, const void *);


int gainMatch()
{
	const int NB_COBO=16;
	const int NB_ASAD=4;
	const int NB_AGET=4;
	const int NB_CHANNEL=68;

	TFile* f=new TFile("pulser_run97.root","read");
	
	TH2F* PadSummary=(TH2F*)(f->Get("PadSummary"));
	
	TCanvas* C=new TCanvas("C","C",600,600);
	
	TF1* f2=new TF1("f2","pol1",0,16385);
	
	PadSummary->Draw("colz");
	
	int co, as, ag, ch;
	int nPeaks;
	bool isOk=false;
	
	TH1D* h=new TH1D("h","h",PadSummary->GetNbinsY(),PadSummary->GetYaxis()->GetXmin(),PadSummary->GetYaxis()->GetXmax());
  
  TH2F* PadSummary_cal=new TH2F("PadSummary_cal","PadSummary_cal",PadSummary->GetNbinsX(),PadSummary->GetXaxis()->GetXmin(),PadSummary->GetXaxis()->GetXmax(),PadSummary->GetNbinsY(),PadSummary->GetYaxis()->GetXmin(),PadSummary->GetYaxis()->GetXmax());
	
	do
	{
		cout << "choose reference channel (co as ag ch): " << endl;
		cin >> co >> as >> ag >> ch;
		C->cd();
		PadSummary->ProjectionY("h",co*NB_ASAD*NB_AGET*NB_CHANNEL+as*NB_AGET*NB_CHANNEL+ag*NB_CHANNEL+ch+1,co*NB_ASAD*NB_AGET*NB_CHANNEL+as*NB_AGET*NB_CHANNEL+ag*NB_CHANNEL+ch+1)->Draw();
		C->Update();
		cout << "Is reference channel Ok? (0/1)" << endl;
		cin >> isOk;
	}
	while(!isOk);
	cout << "Number of peaks?" << endl;
	cin >> nPeaks;
	
	int threshold=20;
	int nint=50;
	
	TSpectrum* S=new TSpectrum(100,1);
	S->Search(h,5,"goff",0.3);
	cout << "found " << S->GetNPeaks() << " peaks ";
	if(S->GetNPeaks()!=nPeaks) cout << "instead of " << nPeaks << endl;
	else cout << endl;
	
	const int npeak=nPeaks;
	float ref[npeak];
	float* REF=(float*)S->GetPositionX();
	for(int j=0;j<npeak;j++)
		ref[j]=REF[j];
	qsort(ref,npeak,sizeof(float),comp);
	
	float P[NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL][npeak];
	
	FILE* fout=fopen("calib_coefs.dat","w");
	
	for(int cobo=0;cobo<NB_COBO;cobo++) for(int asad=0;asad<NB_ASAD;asad++) for(int aget=0;aget<NB_AGET;aget++) for(int channel=0;channel<NB_CHANNEL;channel++)
	{
		int globalI=cobo*NB_ASAD*NB_AGET*NB_CHANNEL+asad*NB_AGET*NB_CHANNEL+aget*NB_CHANNEL+channel;
		PadSummary->ProjectionY("h",globalI+1,globalI+1);
		S->Search(h,5,"goff",0.3);
		float* RREF=(float*)S->GetPositionX();
		if(S->GetNPeaks()!=npeak)
		{
			cout << "abnormal number of peaks: " << S->GetNPeaks() << " for global channel " << globalI << endl;
			fprintf(fout,"%d %d %d %d\t\t%.4f %.1f\n",cobo,asad,aget,channel,0.,0.);
		}
		else
		{
			for(int j=0;j<npeak;j++)
				P[globalI][j]=RREF[j];
			qsort(P[globalI],npeak,sizeof(float),comp);
			
			TGraph* G1=new TGraph(S->GetNPeaks(),P[globalI],ref);
			G1->Fit(f2,"RQ","",-100,20000);
			fprintf(fout,"%d %d %d %d\t\t%.4f %.1f\n",cobo,asad,aget,channel,f2->GetParameter(1),f2->GetParameter(0));
		}
	}
	
	fclose(fout);
	
	TH2F* hcoef_a=new TH2F("hcoef_a","hcoef_a",128,0,128,128,0,128);
  TH2F* hcoef_b=new TH2F("hcoef_b","hcoef_b",128,0,128,128,0,128);
  TH2F* init=new TH2F("init","init",128,0,128,128,0,128);
  
 	float CAL[NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL][2];
 	fout=fopen("calib_coefs.dat","r");
  
  FILE* fTable=fopen("../dat/LT_GANIL_NewCF.dat","r");
  int TABLE[6][NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL];
	for(int i=0;i<NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL;i++)
  {
		fscanf(fTable,"%d %d %d %d %d %d",&TABLE[0][i],&TABLE[1][i],&TABLE[2][i],&TABLE[3][i],&TABLE[4][i],&TABLE[5][i]);
  }
  fclose(fTable);

 	for(int i=0;i<NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL;i++)
 	{
 		float a,b;
 		
    fscanf(fout,"%d %d %d %d %f %f",&co,&as,&ag,&ch,&a,&b);
    int where=co*NB_ASAD*NB_AGET*NB_CHANNEL + as*NB_AGET*NB_CHANNEL + ag*NB_CHANNEL + ch;
    
 		for(int j=0;j<PadSummary->GetNbinsY();j++)
    {
 			float Y=PadSummary->GetYaxis()->GetBinCenter(j+1);
      PadSummary_cal->Fill(i,a*Y+b,PadSummary->GetBinContent(i+1,j+1));
      
      if(TABLE[4][where]!=-1 && TABLE[5][where]!=-1)
      {
        hcoef_a->Fill(TABLE[4][where],TABLE[5][where],a);
 	      hcoef_b->Fill(TABLE[4][where],TABLE[5][where],b);
        if(Y>2000) init->Fill(TABLE[4][where],TABLE[5][where],PadSummary->GetBinContent(i+1,j+1));
      }
    }		
 	}
  
  TCanvas* CC=new TCanvas("CC","CC",700,700);
  CC->Divide(2,2);
  CC->cd(1);
  init->Draw("colz");
  CC->cd(2);
	PadSummary_cal->Draw("colz");
  CC->cd(3);
  hcoef_a->Draw("colz");
  CC->cd(4);
  hcoef_b->Draw("colz");	
 	fclose(fout);
	
	return(0);
}


int comp (const void * a, const void * b)
{
	return ( *(float*)a - *(float*)b );
}
