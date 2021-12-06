#define NBCOBO 4
#define NBASAD 4
#define NBAGET 4
#define NBCHANNEL 68

void AllignPadSummary()
{
	TFile* f=new TFile("PadSummary_run0331.root","READ");
	
	FILE* fcal=fopen("allign_0331.dat","w");
	TH2I* Padsummary=new TH2I(*(TH2I*)(f->Get("PadSummary")));
//	TH2I* Padsummary_cal=new TH2I(*(TH2I*)(f->Get("PadSummary")));

	TCanvas* C=new TCanvas("C","C",600,600);
	
	int cor=2;
	int asr=0;
	int agr=0;
	int chr=0;
	
	float Mref1=978.998;
	float par1[3][3]={{1.58631e+03,9.23820e+02,1.00594e+01},{1.64006e+03,9.86654e+02,9.96197e+00},{1.43729e+03,1.04602e+03,9.20739e+00}};
	
	float Mref2=902.4;
	float par2[3][3]={{2.47150e+0,8.50409e+02,2.54180e+01},{2.81073e+01,9.18360e+02,1.63662e+01},{3.07996e+01,9.76165e+02,1.33941e+01}};

	float ref=Padsummary->ProjectionY("",cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + agr*NBCHANNEL + chr +1,cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + agr*NBCHANNEL + chr +1)->GetMaximumBin()*1.;
	
	TGraph* G=new TGraph(3);
	
	TF1* ffit=new TF1("ffit","gaus(0)+gaus(3)+gaus(6)",700,1200);
	TF1* fpol1=new TF1("fpol1","pol1",700,1200);
	ffit->SetParLimits(2,5,25);
	ffit->SetParLimits(5,5,25);
	ffit->SetParLimits(8,5,25);
	
	for(int ag=0;ag<4;ag++)
		for(int ch=0;ch<68;ch++)
		{
			if(ag<2 && ch!=11 && ch!=22 && ch!=45 && ch!=56)
			{
				Padsummary->ProjectionY("",cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1,cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1)->Draw();
				Padsummary->ProjectionY("",cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1,cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1)->GetXaxis()->SetRangeUser(700,1200);
				float M=Padsummary->ProjectionY("",cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1,cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1)->GetMean(1);
				
				ffit->SetParLimits(1,M-120,M+120);
				ffit->SetParLimits(4,M-120,M+120);
				ffit->SetParLimits(7,M-120,M+120);
				
				for(int i=0;i<3;i++)
					for(int j=0;j<3;j++)
					{
						if(Padsummary->ProjectionY("",cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1,cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1)->GetMaximum()>100)
						{
							if(j==0) ffit->SetParameter(i*3+j,par1[i][j]);
							else  ffit->SetParameter(i*3+j,par1[i][j]*M/Mref1);
						}
						else
						{
							if(j==0) ffit->SetParameter(i*3+j,par2[i][j]);
							else  ffit->SetParameter(i*3+j,par2[i][j]*M/Mref2);
						}
					}

				Padsummary->ProjectionY("",cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1,cor*NBASAD*NBAGET*NBCHANNEL + asr*NBAGET*NBCHANNEL + ag*NBCHANNEL + ch +1)->Fit("ffit","RQ","",700,1200);
				
				for(int i=0;i<3;i++)
					G->SetPoint(i,ffit->GetParameter(i*3+1),par1[i][1]);
				G->Draw("AP*");
				G->Fit("fpol1","RQ","",700,1200);
				fprintf(fcal,"%d\t%d\t%d\t%d\t%f\t%f\n",cor,asr,ag,ch,fpol1->GetParameter(0),fpol1->GetParameter(1));
// 				C->WaitPrimitive();
			}
			else fprintf(fcal,"%d\t%d\t%d\t%d\t%.2f\t%f\n",cor,asr,ag,ch,0.,1.);
		}
	
		
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
