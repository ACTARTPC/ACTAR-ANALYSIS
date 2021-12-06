#include <MEvent.h>


///////////////////////////////////////////////////////////
//
// TreatBaseline method:
// Calibration coefficients availables with CalibCoefs[cobo][asad][aget][channel][2]
// Filling ReducedEvent (Tree branch class) with reduced data
//
// T. Roger 30/07/2017
// email: roger@ganil.fr
//
//////////////////////////////////////////////////////////





int MEvent::TreatBaseline(bool FullTreeOpt, int IsFastPeak)
{//cout << NB_SAMPLES << endl;
	int   BLCALC=  30;
	float THR_PAD= 50.;
	float THR=     50.;
	int   MIN_INT= 50;
	ReducedEvent.event=EN;
	ReducedEvent.timestamp=TS;
	for(MCoboAsad CoAs : CoboAsad)
	{
		int co = CoAs.cobo_number;
		int as = CoAs.asad_number;

		for(unsigned short Ch : CoAs.hit_pattern)
		{
			int ag = (int)(Ch)/NB_CHANNEL;
			int ch = (int)(Ch)%NB_CHANNEL;

			ReducedData DataReduced;
			DataReduced.globalchannelid=ch+(ag<<7)+(as<<9)+(co<<11);
//if(EN==10) cout << co << " " << as << " " << ag << " " << ch << " " << CoAs.Channel[Ch].Sample_Number.size() << endl;
			if(CoAs.Channel[Ch].Sample_Number.size()>0)//==NB_SAMPLES)
			{
				int NMAX=0;
				int P_CHARGE=-1;
				int P_TIME=-1;
				int NSUP=0;
				
				float Baseline=0;
				float BL1=0;
				float BL2=0;
				
				
				if(co!=31) CoAs.Channel[Ch].Raw_Sample[0]=(CoAs.Channel[Ch].Raw_Sample[1]+CoAs.Channel[Ch].Raw_Sample[2])/2.;
								
////////////////////////////////////////////////////////////////////////////
//       Baseline calculation if !VXI and !DSSD && 512 SAMPLES        //////				
////////////////////////////////////////////////////////////////////////////

				if(!FullTreeOpt && co!=31)
				{
					if(CoAs.Channel[Ch].Sample_Number.size()==NB_SAMPLES)
					{
						for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
						{ 
						
							if(Bu<BLCALC) BL1+=(CoAs.Channel[Ch].Raw_Sample[Bu]-BaselineCal[co][as][ag][ch][Bu])/BLCALC;
							if(Bu>NB_SAMPLES-BLCALC-1) BL2+=(CoAs.Channel[Ch].Raw_Sample[Bu]-BaselineCal[co][as][ag][ch][Bu])/BLCALC;
						}
						Baseline=min(BL1,BL2); //Baseline=BL1;
					}
				}
				
////////////////////////////////////////////////////
//      Signal fit if 512 samples             //////				
////////////////////////////////////////////////////
				if(!FullTreeOpt && co!=31 && CoAs.Channel[Ch].Sample_Number.size()==NB_SAMPLES && !IsFastPeak)
				{	
					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
					{
						if(CoAs.Channel[Ch].Raw_Sample[Bu]<4095)
						{
							/*if(hasCalBaseline)*/ Charge[Bu]=(CoAs.Channel[Ch].Raw_Sample[Bu]-BaselineCal[co][as][ag][ch][Bu]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
							/*else Charge[Bu]=(CoAs.Channel[Ch].Raw_Sample[Bu]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];*/
						}
						else
						{
							Charge[Bu]=4096;
							DataReduced.hasSaturation=true;
						}
						Gsig->SetPoint(Bu+1,Bu,Charge[Bu]);
					}
					
					int Bu0=0;
					while(Bu0<NB_SAMPLES-1)
					{
						float ChargeMax=0;
						if(DataReduced.hasSaturation) THR=1000.;
            
						for(unsigned short Bu=Bu0;Bu<NB_SAMPLES;Bu++)
						{	
							Bu0=Bu;
							if(Charge[Bu]>ChargeMax) ChargeMax=Charge[Bu];
// cout << DataReduced.globalchannelid << " " << Bu << " " << Charge[Bu] << " " << ChargeMax << "    " << THR<< endl;
							if(Charge[Bu]<ChargeMax && ChargeMax>THR && ChargeMax<4096)
							{
								ffitSig->SetParameter(0,ChargeMax);
								ffitSig->SetParLimits(0,(ChargeMax)*0.8,(ChargeMax)*1.2);
								ffitSig->SetParameter(1,Bu-143.);
								ffitSig->SetParLimits(1,min((Bu-143.)*0.8,(Bu-143.)*1.2),max((Bu-143.)*0.8,(Bu-143.)*1.2));

								Gsig->Fit(ffitSig,"RQ","",Bu-15,Bu+10);
								DataReduced.peakheight.push_back(ffitSig->GetParameter(0));
								DataReduced.peaktime.push_back(ffitSig->GetParameter(1)+143.);
								//cout << "FIT " << DataReduced.globalchannelid << endl;
					
								for(int bu=0;bu<NB_SAMPLES;bu++)
								{
									Charge[bu]-=ffitSig->Eval(bu);
									Gsig->SetPoint(bu+1,bu,Charge[bu]);
								}
								Bu=NB_SAMPLES;; 
							}
							else if(Charge[Bu]<ChargeMax && ChargeMax>THR && ChargeMax==4096)
							{
								DataReduced.peakheight.push_back(4096);
								DataReduced.peaktime.push_back(Bu0);
								for(int bu=0;bu<NB_SAMPLES;bu++)
								{
									if(bu<Bu0 || Charge[bu]>=THR) Charge[bu]=0;
									else bu=NB_SAMPLES;
									Gsig->SetPoint(bu+1,bu,Charge[bu]);
								}
								Bu=NB_SAMPLES;
							}
						}
					}
				}

///////////////////////////////////////////////
//    Peak find if  FastPeak option ON   //////				
///////////////////////////////////////////////
				else if(!FullTreeOpt && co!=31 && IsFastPeak==1 && (CoAs.Channel[Ch].Sample_Number.size()==NB_SAMPLES))
				{
					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
					{	
						float charge;
						/*if(hasCalBaseline)*/ 
						  charge=(CoAs.Channel[Ch].Raw_Sample[Bu]-BaselineCal[co][as][ag][ch][Bu]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
						/*else
						  charge=(CoAs.Channel[Ch].Raw_Sample[Bu]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];*/
           
//             				cout << Bu << " " << CoAs.Channel[Ch].Raw_Sample[Bu] << " " << Baseline << " " << charge << endl;
						if(Bu>1 && Bu<NB_SAMPLES-1) // pad plane: take only the maxs
						{
							float prev1_charge=(CoAs.Channel[Ch].Raw_Sample[Bu-1]-BaselineCal[co][as][ag][ch][Bu-1]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
							float prev2_charge=(CoAs.Channel[Ch].Raw_Sample[Bu-2]-BaselineCal[co][as][ag][ch][Bu-2]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
							if(charge>P_CHARGE && charge>THR_PAD && charge > prev1_charge && charge > prev2_charge)
								if(!NMAX || (NMAX && charge > (CoAs.Channel[Ch].Raw_Sample[Bu-MIN_INT]-BaselineCal[co][as][ag][ch][Bu-MIN_INT]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1]))
								{
									P_CHARGE=charge;
									P_TIME=Bu;
									NSUP=0;
								}
				
							if(charge<P_CHARGE && P_CHARGE>0) NSUP++;
							if(NSUP==MIN_INT)
							{
								NMAX++;
	
								DataReduced.peakheight.push_back(P_CHARGE);
								DataReduced.peaktime.push_back(P_TIME);
								P_CHARGE=P_TIME=-1;
								NSUP=0;
							}
						}
					}
				}
        
////////////////////////////////////////////////////////////////
//    Peak find if  FastPeak option 2 : Integral method   //////				
////////////////////////////////////////////////////////////////
				else if(!FullTreeOpt && co!=31 && IsFastPeak==2 && (CoAs.Channel[Ch].Sample_Number.size()==NB_SAMPLES))
				{
					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
					{	
						float charge= (CoAs.Channel[Ch].Raw_Sample[Bu]-BaselineCal[co][as][ag][ch][Bu]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
					
            				if(charge>THR_PAD)
						{
							DataReduced.peakheight.push_back(charge);
							DataReduced.peaktime.push_back(Bu);
						}		
					}
				}
        
/////////////////////////////////////////////////
//    Peak find if zero-suppress enabled   //////				
/////////////////////////////////////////////////
        
 				else if(!FullTreeOpt && co!=31 && (CoAs.Channel[Ch].Sample_Number.size()<NB_SAMPLES))
				{
          			if(!hasCalBaseline)
          			{
           				Baseline=0;
						for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
           					BaselineCal[co][as][ag][ch][Bu]=0;
          			}
          
          			for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
					{	
						float charge=(CoAs.Channel[Ch].Raw_Sample[Bu]-BaselineCal[co][as][ag][ch][Bu]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
            						//cout << CoAs.Channel[Ch].Sample_Number.size() << " " << Bu << " " << charge << endl;
						if(Bu>1 && Bu<NB_SAMPLES-2) // pad plane: take only the maxs
						{
							float prev1_charge=(CoAs.Channel[Ch].Raw_Sample[Bu-1]-BaselineCal[co][as][ag][ch][Bu-1]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
							float prev2_charge=(CoAs.Channel[Ch].Raw_Sample[Bu-2]-BaselineCal[co][as][ag][ch][Bu-2]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
							float next1_charge=(CoAs.Channel[Ch].Raw_Sample[Bu+1]-BaselineCal[co][as][ag][ch][Bu+1]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
              float next2_charge=(CoAs.Channel[Ch].Raw_Sample[Bu+2]-BaselineCal[co][as][ag][ch][Bu+2]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
		    					
							if(charge>THR_PAD && charge > prev1_charge && charge > prev2_charge)
								if(charge> next1_charge && charge> next2_charge)
								{	
									DataReduced.peakheight.push_back(charge);
									DataReduced.peaktime.push_back(Bu);
								}
						}
					}
				}



////////////////////////////////////////////////////
//   Full tree option ON: take all samples    //////				
////////////////////////////////////////////////////				
				else if(FullTreeOpt && co!=31)
				{
					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
					{
						DataReduced.peakheight.push_back(CoAs.Channel[Ch].Raw_Sample[Bu]);//-BaselineCal[co][as][ag][ch][Bu]-Baseline);
						DataReduced.peaktime.push_back(Bu);
					}
				}

////////////////////////////////////////////////////
//         VXI parameters treatment           //////				
////////////////////////////////////////////////////				
				else if(co==31)
				{
					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
					{
						DataReduced.peakheight.push_back((CoAs.Channel[Ch].Raw_Sample[Bu]<<16)>>16);
						DataReduced.peaktime.push_back(CoAs.Channel[Ch].Raw_Sample[Bu]>>16);
					}
				}

			  
        
       // cout <<"event#" << EN << "  - ReducedEvent.CoboAsad.size()= " <<ReducedEvent.CoboAsad.size() << "   - DataReduced.peakheight.size()= " << DataReduced.peakheight.size() << "  at add = 0x" << &DataReduced << "    from cobo/asad= " << DataReduced.globalchannelid << endl;
        
				if((co<31 && ch!=11 && ch!=22 && ch!=45 && ch!=56 && !FullTreeOpt) || co>=31 || FullTreeOpt) ReducedEvent.CoboAsad.push_back(DataReduced);
		  }
	  }
  }
	
	return(0);
}



double SignalFitFunction(double* X, double* P)
{
	double x=X[0];
	float Yfit[512]={0.001477,0.001407,0.001472,0.001406,0.001521,0.001585,0.001579,0.001543,0.001497,0.001656,0.001586,0.001603,0.001573,0.001663,0.001569,0.001571,0.001679,0.001521,0.001693,0.001533,0.001633,0.001683,0.001646,0.001371,0.001482,0.001224,0.001208,0.001317,0.001186,0.001113,0.001077,0.001090,0.001106,0.000968,0.000767,0.000883,0.000868,0.000867,0.000768,0.000720,0.000745,0.000720,0.000748,0.000629,0.000732,0.000419,0.000518,0.000564,0.000590,0.000505,0.000280,0.000134,0.000339,0.000343,0.000401,0.000262,0.000247,0.000162,0.000203,0.000071,0.000090,-0.000029,0.000114,-0.000091,0.000017,-0.000092,-0.000058,-0.000147,-0.000046,-0.000603,-0.000265,-0.000254,-0.000267,-0.000376,-0.000351,-0.000522,-0.000249,-0.000488,-0.000383,-0.000534,-0.000482,-0.000501,-0.000419,-0.000465,-0.000673,-0.000580,-0.000294,-0.000411,-0.000313,-0.000334,-0.000153,-0.000407,-0.000525,-0.000339,-0.000196,-0.000201,-0.000147,-0.000188,-0.000094,-0.000222,0.000025,-0.000126,-0.000348,-0.000083,0.000005,-0.000088,0.000210,0.000028,0.000115,0.000299,0.000431,0.000691,0.001814,0.003569,0.007420,0.013836,0.023587,0.038147,0.057499,0.083505,0.115805,0.154736,0.199902,0.250644,0.304893,0.365852,0.424410,0.487061,0.543518,0.603403,0.657186,0.705947,0.756410,0.799159,0.839419,0.875242,0.906395,0.932616,0.955065,0.972800,0.985855,0.994410,1.000000,0.999844,0.997589,0.991925,0.983941,0.972873,0.960011,0.944966,0.928146,0.909500,0.889462,0.868718,0.845567,0.822120,0.797981,0.773560,0.748810,0.723077,0.697922,0.672224,0.646955,0.621514,0.596243,0.571732,0.547277,0.523626,0.499621,0.476996,0.454503,0.432575,0.411516,0.390676,0.370709,0.351772,0.333649,0.315428,0.298565,0.281598,0.266198,0.251358,0.236785,0.223210,0.209615,0.197283,0.185161,0.173997,0.163452,0.152918,0.143277,0.134474,0.125833,0.117888,0.110211,0.102613,0.096264,0.089538,0.083693,0.078248,0.072701,0.067839,0.063431,0.058883,0.054914,0.051056,0.048007,0.044819,0.041488,0.039177,0.036267,0.034116,0.031906,0.029858,0.028215,0.026227,0.024996,0.023302,0.022159,0.020705,0.019665,0.018584,0.017789,0.016865,0.016229,0.015578,0.014788,0.014531,0.013880,0.013419,0.012990,0.012576,0.012538,0.012295,0.012016,0.011665,0.011760,0.011323,0.011271,0.011372,0.011360,0.011191,0.011198,0.011053,0.010999,0.011240,0.011147,0.010901,0.011077,0.010920,0.011017,0.011101,0.011349,0.011269,0.011168,0.011272,0.011418,0.011547,0.011308,0.011507,0.011395,0.011479,0.011515,0.011554,0.011457,0.011437,0.011827,0.011557,0.011815,0.011715,0.011980,0.011743,0.011773,0.011864,0.011694,0.011894,0.011901,0.011663,0.011763,0.011562,0.011847,0.011760,0.011881,0.011684,0.011981,0.011649,0.011777,0.011875,0.012075,0.011832,0.012022,0.011868,0.011939,0.012058,0.012160,0.011784,0.012068,0.011670,0.011801,0.011686,0.011923,0.011848,0.011680,0.011807,0.011662,0.011343,0.011752,0.011668,0.011395,0.011596,0.011243,0.011571,0.011534,0.011133,0.011246,0.011231,0.011307,0.011097,0.011118,0.010872,0.010858,0.010576,0.010894,0.010956,0.011072,0.010840,0.010547,0.010621,0.010739,0.010664,0.010425,0.010883,0.010759,0.010367,0.010559,0.010422,0.010733,0.010220,0.010357,0.010115,0.010052,0.009857,0.010014,0.009539,0.009942,0.010062,0.009977,0.009559,0.009881,0.009749,0.009409,0.009362,0.009561,0.009121,0.009269,0.009170,0.009015,0.008981,0.009099,0.008960,0.009012,0.008963,0.008561,0.008897,0.008740,0.008539,0.008639,0.008431,0.008278,0.008643,0.008237,0.008659,0.008440,0.008289,0.008365,0.008052,0.008523,0.007672,0.007867,0.007656,0.007698,0.007626,0.007804,0.007766,0.007804,0.007858,0.007673,0.007704,0.007570,0.007412,0.007576,0.007290,0.007439,0.007335,0.007525,0.007535,0.007506,0.007048,0.007195,0.006924,0.006601,0.007259,0.007054,0.006885,0.006872,0.006690,0.006336,0.006523,0.006614,0.006497,0.006028,0.005774,0.005791,0.005830,0.005806,0.005910,0.006127,0.006189,0.006185,0.006054,0.006415,0.005597,0.005597,0.005569,0.005812,0.005930,0.005327,0.005781,0.005595,0.005421,0.005719,0.005479,0.005755,0.005588,0.005103,0.005505,0.005521,0.005945,0.005668,0.005594,0.005316,0.005541,0.005625,0.005659,0.005255,0.005213,0.005406,0.005595,0.004987,0.005339,0.005514,0.005260,0.005391,0.004863,0.005632,0.005164,0.005059,0.004963,0.004957,0.004722,0.004952,0.004613,0.004680,0.004283,0.004968,0.004715,0.004304,0.004739,0.004877,0.004508,0.004926,0.004506,0.004550,0.004841,0.004734,0.004667,0.004773,0.004799,0.004708,0.004259,0.004284,0.004268,0.004417,0.004236,0.004087,0.004361,0.004765,0.003824,0.004488,0.004377,0.004173,0.004125,0.004605,0.004373,0.004186,0.004471,0.003989,0.004202,0.004669,0.004678,0.004143,0.004567,0.004292,0.004060,0.004115,0.004107,0.004345,0.004014,0.004388,0.004051,0.003493,0.003712};
	
	double res;
	int ind=(int)((x-P[1]));
	if(ind<1) res=P[0]*Yfit[0];
	else if(ind>511) res=P[0]*Yfit[511];
	else res=P[0]*(Yfit[ind-1] + (Yfit[ind]-Yfit[ind-1])*((x-P[1])-ind) ); // INTERPOLATION A FAIRE!!!
	return(res);
}

