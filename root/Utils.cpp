///////////////////////////////////////
//                                   //
// T. Roger - GANIL 2015             //
// email: roger@ganil.fr             //
//                                   //
// Utils misc. functions file:       //
//     - Automatic pad calibration   //
//     - Pad treatment functions     //
//     - ....                        //
//                                   //
///////////////////////////////////////

#include <TCanvas.h>
#include <TGraph.h>
#include <TSpectrum.h>
#include <math.h>

#include <stdio.h>
#include <iostream>
#include <Parameters.h>

#include <Utils.h>


#define R2D 57.2957795
#define D2R 0.01745329
#define PI  3.14159265

using namespace std;



void Calibrate(TH2I* hPSummary)
{
	TCanvas* Ccal=new TCanvas("Ccal","Ccal",600,600);
	TH2I* padsummary_cal=new TH2I("padsummary_cal","padsummary_cal",NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL,0,NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL, hPSummary->GetNbinsY(),-16,4096);
	TH1D* hPad;
	FILE* fcal=fopen("cal/allign_charge_new.dat","w");
	FILE* fdat=fopen("cal/peaks_new.dat","w");	
	int ncobo, nasad, naget, nchannel;
	int npeaks;
	char OK;
	int cobo0, cobo1;
	
	float min_peak_pct=0.3;
	
	do
	{
		cout << "which cobos need to be calibrated? (first last): ";
		cin >> cobo0 >> cobo1; 
		cout << "choose a reference (cobo asad aget channel): ";
		cin >> ncobo >> nasad >> naget >> nchannel;
		cout << "reference is now cobo " << ncobo << " asad " << nasad << " aget " << naget << "  channel " << nchannel << endl;
		hPad=hPSummary->ProjectionY("",ncobo*NB_ASAD*NB_AGET*NB_CHANNEL + nasad*NB_AGET*NB_CHANNEL + naget*NB_CHANNEL + nchannel +1,ncobo*NB_ASAD*NB_AGET*NB_CHANNEL + nasad*NB_AGET*NB_CHANNEL + naget*NB_CHANNEL + nchannel +1);
		hPad->Draw();
		Ccal->Update();
		cout << "Is reference channel OK (y/n)?" << endl;
		cin >> OK;
	}
	while(OK!='y');
	cout << "how many peaks should be used? ";
	cin >> npeaks;
	
	float CAL[NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL][2];
	
	TSpectrum* S;
	TF1* fitgaus=new TF1("fitgaus","gaus",10,4090);
	TF1* fitpol1=new TF1("fitpol1","pol1",10,4090);
	TGraph* Gallign=new TGraph(npeaks);


	float *REF=0;
	S=new TSpectrum(100,2);
	bool isOK=false;
	do
	{
		S->Search(hPad,10,"nobackground noMarkov",min_peak_pct);
		REF=S->GetPositionX();
		Ccal->Update();
		for (int i=0;i<npeaks;i++)
			for (int j=0;j<npeaks-1;j++)
				if (REF[j]>REF[j+1])
				{
					float var=REF[j];
					REF[j]=REF[j+1];
					REF[j+1]=var;
				}
		if(npeaks==S->GetNPeaks()) isOK=true;
		else
		{
			cout << "Could not find " << npeaks << " peaks with the current setting \nPlease enter new search limit" << endl;
			cin >> min_peak_pct;
		}
	}
	while(!isOK);
// 	Ccal->WaitPrimitive();
	
	for(int i=0;i<npeaks;i++)
	{
		hPad->Fit(fitgaus,"RQ","",REF[i]-10,REF[i]+10);
// 		if(i==0) hPad->Fit(fitgaus,"RQ","",0,REF[0]+(REF[1]-REF[0])/2.);
// 		else if(i==npeaks-1) hPad->Fit(fitgaus,"RQ","",REF[npeaks-1]-(REF[npeaks-1]-REF[npeaks-2])/2.,4095);
// 		else hPad->Fit(fitgaus,"RQ","",REF[i]-(REF[i]-REF[i-1])/2.,REF[i]+(REF[i+1]-REF[i])/2.);
		REF[i]=fitgaus->GetParameter(1);
		cout << "Unordered Reference peak #" << i+1 << ": mean = " << fitgaus->GetParameter(1) << "  sigma = " << fitgaus->GetParameter(2) << endl;
	}
	for (int i=0;i<npeaks;i++)
		for (int j=0;j<npeaks-1;j++)
			if (REF[j]>REF[j+1])
			{
				float var=REF[j];
				REF[j]=REF[j+1];
				REF[j+1]=var;
			}
	for(int i=0;i<npeaks;i++) cout << "Reference peak #" << i+1 << ": at " << REF[i] << endl;
	Ccal->Update();
	

	TSpectrum* s=new TSpectrum(100,2);
	for(int cobo=cobo0;cobo<=cobo1;cobo++)
		for(int asad=0;asad<NB_ASAD ;asad++)
	for(int aget=0;aget<NB_AGET;aget++)
		for(int channel=0;channel<NB_CHANNEL ;channel++)
		{
			float *pos=0;
			hPad=hPSummary->ProjectionY("",cobo*NB_ASAD*NB_AGET*NB_CHANNEL + asad*NB_AGET*NB_CHANNEL + aget*NB_CHANNEL + channel +1,cobo*NB_ASAD*NB_AGET*NB_CHANNEL + asad*NB_AGET*NB_CHANNEL + aget*NB_CHANNEL + channel +1);
			s->Search(hPad,10,"nobackground noMarkov",min_peak_pct);
			pos=s->GetPositionX();
			
// 			hPad->Draw();
// 			Ccal->Update();
// 			Ccal->WaitPrimitive();
			if(s->GetNPeaks()==npeaks)
			{
				for(int loop=0;loop<2;loop++)
					for (int i=0;i<npeaks;i++)
						for (int j=0;j<npeaks-1;j++)
							if (pos[j]>pos[j+1])
							{
								float var=pos[j];
								pos[j]=pos[j+1];
								pos[j+1]=var;
							}
// 				for(int i=0;i<npeaks;i++) cout << "Found peak #" << i+1 << ": at " << pos[i] << endl;
							
				for(int i=0;i<npeaks;i++)
				{
					hPad->Fit(fitgaus,"RQ","",pos[i]-10,pos[i]+10);
// 					if(i==0) hPad->Fit(fitgaus,"RQ","",0,pos[0]+(pos[1]-pos[0])/2.);
// 					else if(i==npeaks-1) hPad->Fit(fitgaus,"RQ","",pos[npeaks-1]-(pos[npeaks-1]-pos[npeaks-2])/2.,4095);
// 					else hPad->Fit(fitgaus,"RQ","",pos[i]-(pos[i]-pos[i-1])/2.,pos[i]+(pos[i+1]-pos[i])/2.);
					Gallign->SetPoint(i,fitgaus->GetParameter(1),REF[i]);
					fprintf(fdat,"%d %d %d %d %d\t%f\n",cobo,asad,aget,channel,i+1,fitgaus->GetParameter(1));
				}
				Gallign->Fit(fitpol1,"RQ","",0,4095);
// 				Gallign->Draw("AP*");
// 				Ccal->Update();
// 				Ccal->WaitPrimitive();
				fprintf(fcal,"%d\t%d\t%d\t%d\t%.2f\t%.6f\n",cobo,asad,aget,channel,fitpol1->GetParameter(0),fitpol1->GetParameter(1));

			}
			else
			{
				if(channel!=11 && channel !=22 && channel!=45 && channel!=56)
				{
					cout << "cobo " << cobo << " asad " << asad << " aget " << aget << "  channel " << channel << " cannot be calibrated: found " << s->GetNPeaks() << " peaks instead of " << npeaks << endl;
					hPad->Draw();
					Ccal->Update();
					Ccal->WaitPrimitive();
				}
				fprintf(fcal,"%d\t%d\t%d\t%d\t%.2f\t%.6f\n",cobo,asad,aget,channel,0.,1.);

			}
			CAL[cobo*NB_ASAD*NB_AGET*NB_CHANNEL + asad*NB_AGET*NB_CHANNEL + aget*NB_CHANNEL + channel][0]=fitpol1->GetParameter(0);
			CAL[cobo*NB_ASAD*NB_AGET*NB_CHANNEL + asad*NB_AGET*NB_CHANNEL + aget*NB_CHANNEL + channel][1]=fitpol1->GetParameter(1);
		}
	
	fclose(fcal);
	fclose(fdat);
	for(int i=1;i<=NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL;i++)
		for(int j=1;j<=padsummary_cal->GetNbinsY();j++)
			padsummary_cal->Fill(i-1,CAL[i-1][0]+CAL[i-1][1]*hPSummary->GetYaxis()->GetBinCenter(j),hPSummary->GetBinContent(i,j));
	padsummary_cal->Draw("colz");
	Ccal->Update();
	Ccal->WaitPrimitive();
}


void CleanPad(float PAD[NPADY][NPADX],float TIME[NPADY][NPADX])
{
	short NEIGHBOUR[NPADY][NPADX];
	for(int r=0;r<NPADY;r++)
		for(int c=0;c<NPADX;c++)
		{
			NEIGHBOUR[r][c]=0;
			if(r>0)
			{
				if(c>0 && PAD[r-1][c-1]) NEIGHBOUR[r][c]++;
				if(c<NPADX-1 && PAD[r-1][c+1]) NEIGHBOUR[r][c]++;
				if(PAD[r-1][c]) NEIGHBOUR[r][c]++;
			}
			if(r<NPADY-1)
			{
				if(c>0 && PAD[r+1][c-1]) NEIGHBOUR[r][c]++;
				if(c<NPADX-1 && PAD[r+1][c+1]) NEIGHBOUR[r][c]++;
				if(PAD[r+1][c]) NEIGHBOUR[r][c]++;
			}
			if(c>0 && PAD[r][c-1]) NEIGHBOUR[r][c]++; 
			if(c<NPADX-1 && PAD[r][c+1]) NEIGHBOUR[r][c]++; 
			
// 			if(c>=40 && c<=47 && r<=10) NEIGHBOUR[r][c]=3;
// 			if(c<30 && r<12) NEIGHBOUR[r][c]=0;
		}
	
	for(int r=0;r<NPADY;r++)
		for(int c=0;c<NPADX;c++)
		{
			if(r>0 && c>0 && r<NPADY-1 && c<NPADX-1 && NEIGHBOUR[r][c]<2) PAD[r][c]=TIME[r][c]=0;
			if((r==0 || c==0 || r==NPADY-1 || c==NPADX-1) && NEIGHBOUR[r][c]<1) PAD[r][c]=TIME[r][c]=0;
		}
}

bool GetVertexFromRMS(TH2F* visu_charge,int& i_vertex,int& i_end)
{
	float BLM=0;
	float BLRMS=0;
	float RMSdev=2;
	int cpt=0;
	i_vertex=NPADX;
	i_end=0;
	for(int c=NPADX-2;c>=NPADX-9;c--)
	{
		BLRMS+=visu_charge->ProjectionY("",c,c)->GetRMS();
		BLM+=visu_charge->ProjectionY("",c,c)->GetMean();
	}
	BLRMS/=7.;
	BLM/=7.;
	
	if(BLM<15 || BLM>20) return(false);
	
	for(int c=NPADX-10;c>=0;c--)
	{
		if(visu_charge->ProjectionY("",c,c)->GetRMS()>RMSdev*BLRMS)
		{
			cpt++;
			i_vertex=c;
		}
		else
		{
			cpt=0;
			i_vertex=NPADX;
		}
		if(cpt>2) break;
	}
	cpt=0;
	for(int c=i_vertex+2;c>=0;c--)
	{
		if(visu_charge->ProjectionY("",c,c)->GetRMS()<RMSdev*BLRMS)
		{
			cpt++;
			i_end=c;
		}
		else
		{
			cpt=0;
			i_end=0;
		}
		if(cpt>2) break;
	}
	if(cpt>2) return(true);
	else return(false);
}


bool GetVerticalConsistensy(MTrack* T1,MTrack* T2)
{
	double a1=(T1->Yh-T1->Ym)/(T1->Xh-T1->Xm);
	double a2=(T2->Yh-T2->Ym)/(T2->Xh-T2->Xm);
	double b1=T1->Ym-a1*T1->Xm;
	double b2=T2->Ym-a2*T2->Xm;
	
	double x=(b2-b1)/(a1-a2);
	double y=a1*x+b1;
	
	double alpha1=(y-T1->Yh)/(T1->Yh-T1->Ym);
	double alpha2=(y-T2->Yh)/(T2->Yh-T2->Ym);
	
	double Dz=alpha2*(T2->Zh-T2->Zm)+T2->Zh-(alpha1*(T1->Zh-T1->Zm)+T1->Zh);
	
	if(fabs(Dz)<5) return(true);
	else return(false);
}


float AverageEloss(MTrack* T, float PAD[NPADY][NPADX])
{
	float Sum=0;
	float a=(T->Yh-T->Ym)/(T->Xh-T->Xm);
	float b=T->Ym-a*T->Xm;
	for(int r=T->zy_s;r<=NPADY;r++)
		for(int c=0;c<NPADX;c++)
			if((a*(c+1.)+b-(r+1.))/sqrt(a*a+1.) < 5.)
				Sum+=PAD[r][c];
	
	float length;
	if(T->Ym-(T->Yh-T->Ym)*T->Xm/(T->Xh-T->Xm) > NPADY)
	{
		length=sqrt(pow(T->Xm+(T->Xh-T->Xm)*(T->zy_s-T->Ym)/(T->Yh-T->Ym) - (T->Xm+(T->Xh-T->Xm)*(NPADY-T->Ym)/(T->Yh-T->Ym)),2)+pow(NPADY-T->zy_s,2));
		T->L2DXY->SetX1(T->Xm+(T->Xh-T->Xm)*(T->zy_s-T->Ym)/(T->Yh-T->Ym));
		T->L2DXY->SetX2(T->Xm+(T->Xh-T->Xm)*(NPADY-T->Ym)/(T->Yh-T->Ym));
		T->L2DXY->SetY1(T->zy_s);
		T->L2DXY->SetY2(NPADY);
	}
	else
	{
		length=sqrt(pow(T->Xm+(T->Xh-T->Xm)*(T->zy_s-T->Ym)/(T->Yh-T->Ym),2)+pow(T->Ym-(T->Yh-T->Ym)*T->Xm/(T->Xh-T->Xm)- T->zy_s,2));
		T->L2DXY->SetX1(T->Xm+(T->Xh-T->Xm)*(T->zy_s-T->Ym)/(T->Yh-T->Ym));
		T->L2DXY->SetX2(0);
		T->L2DXY->SetY1(T->zy_s);
		T->L2DXY->SetY2(T->Ym-(T->Yh-T->Ym)*T->Xm/(T->Xh-T->Xm));
	}
	
	length=sqrt(pow(length,2)+pow(T->Zm+(T->Zh-T->Zm)*(T->L2DXY->GetY1()-T->Ym)/(T->Yh-T->Ym) - T->Zm+(T->Zh-T->Zm)*(T->L2DXY->GetY2()-T->Ym)/(T->Yh-T->Ym) ,2));
// 	cout << length << endl;
	
	return(Sum/length);
}


void FitMat(float PAD[NPADY][NPADX], int Rmin, int Rmax, int Cmin, int Cmax, float threshold, float &a, float &b)
{
	int Row, Col;
	float A, B, C, UEV, Q, X, Xg, Y, Yg;
	A=B=C=UEV=Q=X=Y=Xg=Yg=0.;
	
	for (Row=Rmin;Row<=Rmax;Row++)
		for (Col=Cmin;Col<=Cmax;Col++)
			if(PAD[Row][Col]>threshold)
			{
				X= Col*2.+1.;
				Y= Row*2.+1.;
				Q+=PAD[Row][Col];
				Xg+=X*PAD[Row][Col];
				Yg+=Y*PAD[Row][Col];
			}
	Xg/=Q;
	Yg/=Q;
	
	for (Row=Rmin;Row<=Rmax;Row++)
		for(Col=Cmin;Col<=Cmax;Col++)
			if(PAD[Row][Col]>threshold)
			{
				X= Col*2.+1.;
				Y= Row*2.+1.;
				A+=PAD[Row][Col]*(X-Xg)*(X-Xg);
				B+=PAD[Row][Col]*(X-Xg)*(Y-Yg);
				C+=PAD[Row][Col]*(Y-Yg)*(Y-Yg);
			}
	UEV=0.5*(A+C+sqrt((A+C)*(A+C)-4*(A*C-B*B)));
	a=B/(UEV-C);
	b=Yg-a*Xg;
}


float FitMat3D(float PAD[NPADY][NPADX], float TIME[NPADY][NPADX], float threshold, MTrack* T)
{
	// adapted from: http://fr.scribd.com/doc/31477970/Regressions-et-trajectoires-3D
	
	int Rmin=T->zx_s;
	int Rmax=T->zx_e;
	int Cmin=T->zy_s;
	int Cmax=T->zy_e;
	int R, C;
	double Q,X,Y,Z;
	double Xm,Ym,Zm;
	double Sxx,Sxy,Syy,Sxz,Szz,Syz;
	double theta;
	double K11,K22,K12,K10,K01,K00;
	double c0,c1,c2;
	double p,q,r,dm2;
	double rho,phi;
	double a,b;
	
	Q=Xm=Ym=Zm=0.;
	Sxx=Syy=Szz=Sxy=Sxz=Syz=0.;
	
	for (R=Rmin;R<=Rmax;R++)
		for (C=Cmin;C<=Cmax;C++)
			if(PAD[R][C]>threshold && TIME[R][C])
			{
				X= C+0.5;//*2.+1.;
				Y= R+0.5;//*2.+1.;
				Z= TIME[R][C];//*VDRIFT;
				Q+=PAD[R][C]/10.;
				Xm+=X*PAD[R][C]/10.;
				Ym+=Y*PAD[R][C]/10.;
				Zm+=Z*PAD[R][C]/10.;
				Sxx+=X*X*PAD[R][C]/10.;
				Syy+=Y*Y*PAD[R][C]/10.;
				Szz+=Z*Z*PAD[R][C]/10.;
				Sxy+=X*Y*PAD[R][C]/10.;
				Sxz+=X*Z*PAD[R][C]/10.;
				Syz+=Y*Z*PAD[R][C]/10.;
			}
	Xm/=Q;
	Ym/=Q;
	Zm/=Q;
	Sxx/=Q;
	Syy/=Q;
	Szz/=Q;
	Sxy/=Q;
	Sxz/=Q;
	Syz/=Q;
	Sxx-=(Xm*Xm);
	Syy-=(Ym*Ym);
	Szz-=(Zm*Zm);
	Sxy-=(Xm*Ym);
	Sxz-=(Xm*Zm);
	Syz-=(Ym*Zm);
	
	theta=0.5*atan((2.*Sxy)/(Sxx-Syy));
	
	K11=(Syy+Szz)*pow(cos(theta),2)+(Sxx+Szz)*pow(sin(theta),2)-2.*Sxy*cos(theta)*sin(theta);
	K22=(Syy+Szz)*pow(sin(theta),2)+(Sxx+Szz)*pow(cos(theta),2)+2.*Sxy*cos(theta)*sin(theta);
	K12=-Sxy*(pow(cos(theta),2)-pow(sin(theta),2))+(Sxx-Syy)*cos(theta)*sin(theta);
	K10=Sxz*cos(theta)+Syz*sin(theta);
	K01=-Sxz*sin(theta)+Syz*cos(theta);
	K00=Sxx+Syy;
	
	c2=-K00-K11-K22;
	c1=K00*K11+K00*K22+K11*K22-K01*K01-K10*K10;
	c0=K01*K01*K11+K10*K10*K22-K00*K11*K22;
	
	
	p=c1-pow(c2,2)/3.;
	q=2.*pow(c2,3)/27.-c1*c2/3.+c0;
	r=pow(q/2.,2)+pow(p,3)/27.;
	
	
	if(r>0) dm2=-c2/3.+pow(-q/2.+sqrt(r),1./3.)+pow(-q/2.-sqrt(r),1./3.);
	if(r<0)
	{
		rho=sqrt(-pow(p,3)/27.);
		phi=acos(-q/(2.*rho));
		dm2=min(-c2/3.+2.*pow(rho,1./3.)*cos(phi/3.),min(-c2/3.+2.*pow(rho,1./3.)*cos((phi+2.*M_PI)/3.),-c2/3.+2.*pow(rho,1./3.)*cos((phi+4.*M_PI)/3.)));
	}
	a=-K10*cos(theta)/(K11-dm2)+K01*sin(theta)/(K22-dm2);
	b=-K10*sin(theta)/(K11-dm2)-K01*cos(theta)/(K22-dm2);

	T->Xm=Xm;
	T->Ym=Ym;
	T->Zm=Zm;
	T->Xh=((1.+b*b)*Xm-a*b*Ym+a*Zm)/(1.+a*a+b*b);
	T->Yh=((1.+a*a)*Ym-a*b*Xm+b*Zm)/(1.+a*a+b*b);
	T->Zh=((a*a+b*b)*Zm+a*Xm+b*Ym)/(1.+a*a+b*b);
	
	T->L2DXY->SetX1(T->Xm);
	T->L2DXY->SetY1(T->Ym);
	T->L2DXY->SetX2(T->Xh);
	T->L2DXY->SetY2(T->Yh);
	
	T->L2DXZ->SetX1(T->Xm);
	T->L2DXZ->SetY1(T->Zm);
	T->L2DXZ->SetX2(T->Xh);
	T->L2DXZ->SetY2(T->Zh);
	
	T->L2DYZ->SetX1(T->Ym);
	T->L2DYZ->SetY1(T->Zm);
	T->L2DYZ->SetX2(T->Yh);
	T->L2DYZ->SetY2(T->Zh);
	
	return(dm2/Q);
}


float DEdX(MTrack* T, float PAD[NPADY][NPADX], float dmax, TH1D* profile)
{
	float thetaX_Y=atan((T->Yh-T->Ym)/(T->Xh-T->Xm));
	float thetaXY_Z=atan((T->Zh-T->Zm)/sqrt(pow(T->Yh-T->Ym,2)+pow(T->Xh-T->Xm,2)));
	float a=(T->Yh-T->Ym)/(T->Xh-T->Xm);
	float b=T->Yh-T->Xh*(T->Yh-T->Ym)/(T->Xh-T->Xm);
	float Qtot=0;

	for (int R=T->zy_s;R<=T->zy_e;R++)
		for (int C=T->zx_s;C<=T->zx_e;C++)
			if(sqrt((a*(C+0.5) + b - (R+0.5))/(a*a+1.))<dmax)
			{
				Qtot+=PAD[R][C];
				profile->Fill((cos(thetaX_Y)*(C+0.5) + sin(thetaX_Y)*(R+0.5))/sin(thetaXY_Z),PAD[R][C]);
			}
	return(Qtot);
}


void MakeConfigFileThr(char* f_base_name, char* f_name, float BL[NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL],float THR)
{
	FILE* f_in=fopen(f_base_name,"r");
	FILE* f_out=fopen(f_name,"w");
	
	char line[2048];
	if(f_in!=NULL && f_out!=NULL)
	{
		for(int i=0;i<2;i++)
		{
			fgets(line,2048,f_in);
			fputs(line,f_out);
		}
		do
		{
			fgets(line,2048,f_in);
			if(strcmp(line,"</Setup>")!=1) fputs(line,f_out);
		}
		while(!feof(f_in));
		
		char instance[256];
		
		fprintf(f_out,"\t<Node id=\"CoBo\">\n");
		for(int i1=0;i1<2;i1++)
		{
			sprintf(instance,"<Instance id=\"Crate00_Slot0%d\">",i1);
			fprintf(f_out,"\t\t%s\n",instance);
			for(int i2=0;i2<NB_ASAD;i2++)
			{
				sprintf(instance,"<AsAd id=\"%d\">",i2);
				fprintf(f_out,"\t\t\t%s\n",instance);
				for(int i3=0;i3<NB_AGET;i3++)
				{
					sprintf(instance,"<Aget id=\"%d\">",i3);
					fprintf(f_out,"\t\t\t\t%s\n",instance);
					for(int i4=0;i4<NB_CHANNEL;i4++)
					{
						sprintf(instance,"<channel id=\"%d\">",i4);
						fprintf(f_out,"\t\t\t\t\t%s\n",instance);
						sprintf(instance,"<zeroSuppressionThreshold>%d</zeroSuppressionThreshold>",(int)(BL[i1*NB_ASAD*NB_AGET*NB_CHANNEL + i2*NB_AGET*NB_CHANNEL + i3*NB_CHANNEL + i4]+THR));
						fprintf(f_out,"\t\t\t\t\t\t%s\n",instance);
						fprintf(f_out,"\t\t\t\t\t</channel>\n");
					}
					fprintf(f_out,"\t\t\t\t</Aget>\n");
				}
				fprintf(f_out,"\t\t\t</AsAd>\n");
			}
			fprintf(f_out,"\t\t</Instance>\n");
		}
		fprintf(f_out,"\t</Node>\n");
		fprintf(f_out,"</Setup>\n");
		fclose(f_in);
		fclose(f_out);
	}
	else cout << "file " << f_base_name << " does not exists" << endl;
}

