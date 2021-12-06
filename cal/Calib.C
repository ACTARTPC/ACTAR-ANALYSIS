#include <TFile.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TSpectrum.h>
#include <TCanvas.h>
#include <iostream>
#include <cstdlib>

#define NBCOBO 16
#define NBASAD 4
#define NBAGET 4
#define NBCHANNEL 68

using namespace std;

int comp (const void *, const void *);

void Calib()
{
	TFile* f=new TFile("../root_e796/PadSummary.root","READ");
	
	FILE* fcal=fopen("allign_run_8-12_NoDeconv_2nd.dat","w");
	TH2F* Padsummary=new TH2F(*(TH2F*)(f->Get("hPADSummary")));
//	TH2I* Padsummary_cal=new TH2I(*(TH2I*)(f->Get("PadSummary")));

	TCanvas* C=new TCanvas("C","C",600,600);
	
	
	int cor=2;
	int asr=0;
	int agr=0;
	int chr=0;
	
	TH1D* refspectrum;
	
	char isRefOK='N';
	float* REF;
	int NpeaksRef;
	
	TSpectrum* S=new TSpectrum(200,1);
	
	do
	{
		do
		{
			cout << "choose reference channel: cobo asad aget channel: ";
			cout.flush();
			cin >> cor >> asr >> agr >> chr; 
		
			refspectrum=Padsummary->ProjectionY("",cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + agr*NBCHANNEL + chr +1,cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + agr*NBCHANNEL + chr +1);
		
			refspectrum->Draw();
			C->Update();
		
			cout << "Is reference channel OK?"; cout.flush();
			cin >> isRefOK;
		}
		while(isRefOK=='n' || isRefOK=='N');
		
		cout << "Number of peaks to find: "; cout.flush();
		cin >> NpeaksRef;
	
	
		S->Search(refspectrum,0.01,"goff",0.1);
		if(S->GetNPeaks()!=NpeaksRef) cout << "found " << S->GetNPeaks() << " peaks in reference spectrum instead of " << NpeaksRef << endl;
		refspectrum->Draw(); C->Update();
	}		
	while(S->GetNPeaks()!=NpeaksRef);

  int order=1;
  cout << "Calib Order (1/2): ";
  cout.flush();
  cin >> order; 

	const int npeak=NpeaksRef;
	float ref[npeak];
	for(int k=0;k<npeak;k++) ref[k]=0;
	REF=(float*)S->GetPositionX();
	for(int i=0;i<npeak;i++)
		ref[i]=REF[i];		
	
	qsort(ref,npeak,sizeof(float),comp);
	
	delete S;
	
	float pos[npeak];
	TH1D* projspectrum;

	TF1* f1=new TF1("f1","gaus",0,4095);
	TF1* f2;
  if(order==2) f2=new TF1("f2","pol2",0,4095);
  else f2=new TF1("f2","pol1",0,4095);
  
  
	for(int co=0;co<NBCOBO;co++)for(int as=0;as<NBASAD;as++)for(int ag=0;ag<NBAGET;ag++)for(int ch=0;ch<NBCHANNEL;ch++) if(ch!=11 && ch!=22 && ch!=45 && ch!=56)
	{
		int nch=co*NBASAD*NBAGET*NBCHANNEL + as*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch;
		
		S=new TSpectrum(200,0.5);
		projspectrum=Padsummary->ProjectionY("",nch+1,nch+1);
		projspectrum->Rebin(2);
		
		Float_t *POS=new Float_t[400];
		S->Search(projspectrum,0.01,"",0.2);			
		POS=S->GetPositionX();
		for(int k=0;k<npeak;k++) pos[k]=0;

		if (POS)
		{
			if(S->GetNPeaks()!=npeak)
			{
				cout << "abnormal number of peaks: " << S->GetNPeaks() << " cobo/asad/aget/channel: " << co << " " << as << " " << ag << " " << ch << endl;
				projspectrum->Draw();
				if(S->GetNPeaks()>1)
				{
					C->Update();
					C->WaitPrimitive();
				}
				if(order==1) fprintf(fcal,"%d\t%d\t%d\t%d\t\t%.5e\t%.5e\n",co,as,ag,ch,0.0,0.0);
        if(order==2) fprintf(fcal,"%d\t%d\t%d\t%d\t\t%.5e\t%.5e\t%.5e\n",co,as,ag,ch,0.0,0.0,0.0);
			}
				
			else if(S->GetNPeaks()==npeak)
			{
				for(int k=0;k<min(npeak,S->GetNPeaks());k++)
					pos[k]=POS[k];
				qsort(pos,S->GetNPeaks(),sizeof(float),comp);
				
				projspectrum->Fit(f1,"RQ","",pos[0]-10,pos[0]+10);
				TGraph* G1=new TGraph(S->GetNPeaks(),pos,ref);
				f2->SetParameter(0,0);
				f2->SetParameter(1,1);
				f2->SetParameter(2,0);
				G1->Fit(f2,"RQ","",pos[0]-1,pos[npeak-1]+1);
				if(order==1) fprintf(fcal,"%d\t%d\t%d\t%d\t\t%.8e\t%.8e\n",co,as,ag,ch,f2->GetParameter(1),f2->GetParameter(0));
        if(order==2) fprintf(fcal,"%d\t%d\t%d\t%d\t\t%.8e\t%.8e\t%.8e\n",co,as,ag,ch,f2->GetParameter(2),f2->GetParameter(1),f2->GetParameter(0));
// 				G1->Draw("AP*");
// 				C->Update();
// 				C->WaitPrimitive();
				delete G1;
			}
		}
// 			
		else
		{
			cout << "NO PEAKS FOUND in cobo/asad/aget/channel: " << co << " " << as << " " << ag << " " << ch << endl;
			if(order==1) fprintf(fcal,"%d\t%d\t%d\t%d\t\t%.5e\t%.5e\n",co,as,ag,ch,0.0,0.0);
      if(order==2) fprintf(fcal,"%d\t%d\t%d\t%d\t\t%.5e\t%.5e\t%.5e\n",co,as,ag,ch,0.0,0.0,0.0);
			projspectrum->Draw();
			C->Update();
			C->WaitPrimitive();
		}		
// 		fprintf(cal[i],"\n");
// 		fflush(cal[i]);		
		
		
		delete S;
	}
	else 
  {
    if(order==1) fprintf(fcal,"%d\t%d\t%d\t%d\t\t%.5e\t%.5e\n",co,as,ag,ch,0.0,0.0);
    if(order==2) fprintf(fcal,"%d\t%d\t%d\t%d\t\t%.5e\t%.5e\t%.5e\n",co,as,ag,ch,0.0,0.0,0.0);
	}
	
	
	C->Update();
	
// 	int ind=0;
// 	do
// 	{
// 		S->Search(h[i][(int)(nch/2)+(int)pow((int)(ind/2),ind%2)-1],5,"goff",0.05);
// 		cout << "found " << S->GetNPeaks() << " peaks in detector " << dec <<(int)(nch/2)+(int)pow((int)(ind/2),ind%2)-1 << endl;
// 		h[i][(int)(nch/2)+(int)pow((int)(ind/2),ind%2)-1]->Draw();
// 		CC->WaitPrimitive();
// 		ind++;
// 	}
// 	while(S->GetNPeaks()!=npeak && ind<nch);
// 	
// 	for(int ag=0;ag<4;ag++)
// 		for(int ch=0;ch<68;ch++)
// 		{
// 			if(ag<2 && ch!=11 && ch!=22 && ch!=45 && ch!=56)
// 			{
// 				Padsummary->ProjectionY("",cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1,cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1)->Draw();
// 				Padsummary->ProjectionY("",cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1,cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1)->GetXaxis()->SetRangeUser(700,1200);
// 				float M=Padsummary->ProjectionY("",cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1,cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1)->GetMean(1);
// 				
// 				ffit->SetParLimits(1,M-120,M+120);
// 				ffit->SetParLimits(4,M-120,M+120);
// 				ffit->SetParLimits(7,M-120,M+120);
// 				
// 				for(int i=0;i<3;i++)
// 					for(int j=0;j<3;j++)
// 					{
// 						if(Padsummary->ProjectionY("",cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1,cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1)->GetMaximum()>100)
// 						{
// 							if(j==0) ffit->SetParameter(i*3+j,par1[i][j]);
// 							else  ffit->SetParameter(i*3+j,par1[i][j]*M/Mref1);
// 						}
// 						else
// 						{
// 							if(j==0) ffit->SetParameter(i*3+j,par2[i][j]);
// 							else  ffit->SetParameter(i*3+j,par2[i][j]*M/Mref2);
// 						}
// 					}
// 
// 				Padsummary->ProjectionY("",cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1,cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1)->Fit("ffit","RQ","",700,1200);
// 				
// 				for(int i=0;i<3;i++)
// 					G->SetPoint(i,ffit->GetParameter(i*3+1),par1[i][1]);
// 				G->Draw("AP*");
// 				G->Fit("fpol1","RQ","",700,1200);
// 				fprintf(fcal,"%d\t%d\t%d\t%d\t%f\t%f\n",cor,asr,ag,ch,fpol1->GetParameter(0),fpol1->GetParameter(1));
// // 				C->WaitPrimitive();
// 			}
// 			else fprintf(fcal,"%d\t%d\t%d\t%d\t%.2f\t%f\n",cor,asr,ag,ch,0.,1.);
// 		}
// 	
		
// 	for(int co=0;co<NBCOBO;co++)
// 		for(int as=0;as<NBASAD;as++)
// 			for(int ag=0;ag<NBAGET;ag++)
// 				for(int ch=0;ch<NBCHANNEL;ch++)
// 				{
// 					if(co<2) fprintf(fcal,"%d\t%d\t%d\t%d\t0.000\t%.3f\n",co,as,ag,ch,ref/Padsummary->ProjectionY("",co*NBASAD*NBAGET*NBCHANNEL + as*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1,co*NBASAD*NBAGET*NBCHANNEL + as*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1)->GetMaximumBin()*1.);
// 					else  fprintf(fcal,"%d\t%d\t%d\t%d\t0.000\t1.000\n",co,as,ag,ch);
// 				}
// 	
	fclose(fcal);
	
// 	Padsummary->Draw("colz");
}


int comp (const void * a, const void * b)
{
	return ( *(float*)a - *(float*)b );
}
