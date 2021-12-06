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



std::vector<int> readalwayslist={1703,1707,1888,1890,1898,1902,2142,2146,2150,2155,2165,2167,2262,2264,2431,2433,2267,2269,2707,2710,3230,3232,3234,3236};



int MEvent::TreatFullSignal(bool FullTreeOpt)
{
	int   BLCALC=         30;
	float THR_PAD_DECONV= 0.5;
	float THR=            50.;
	int   MIN_INT=        50;
	
	TH2F* hBaseline=new TH2F("hBaseline","hBaseline",NB_SAMPLES,0,NB_SAMPLES,4096,0,4096);
	
	ReducedEvent.event=EN;
	ReducedEvent.timestamp=TS;
	

// /**/	fprintf(f_out,"%d ",EN);
// /**/	int nch=0;
	
	for(MCoboAsad CoAs : CoboAsad)
	{
		int co = CoAs.cobo_number;
		int as = CoAs.asad_number;

		if(co!=16 && co!=31) for(unsigned short Ch : CoAs.hit_pattern)
		{
			int ag = (int)(Ch)/NB_CHANNEL;
			int ch = (int)(Ch)%NB_CHANNEL;

			int globalchannelid=ch+(ag<<7)+(as<<9)+(co<<11);
			
// /**/			nch++;
						
			if(CoAs.Channel[Ch].Sample_Number.size()>0)//==NB_SAMPLES)
			{
				int NMAX=0;
				int P_CHARGE=-1;
				int P_TIME=-1;
				int NSUP=0;
				
				float Baseline=0;
				float BL1=0;
				float BL2=0;
				
				
////////////////////////////////////////////////////////////////////////////
//       Baseline calculation if !VXI and !DSSD && 512 SAMPLES        //////				
////////////////////////////////////////////////////////////////////////////

				for(auto readalwayspad:readalwayslist)
				{
					if(globalchannelid==readalwayspad)
					{
						CoAs.Channel[Ch].Raw_Sample[0]=(CoAs.Channel[Ch].Raw_Sample[1]+CoAs.Channel[Ch].Raw_Sample[2])/2.;
						for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
						{ 
							if(Bu<BLCALC) BL1+=(CoAs.Channel[Ch].Raw_Sample[Bu]/BLCALC);
							if(Bu>NB_SAMPLES-BLCALC-2) BL2+=(CoAs.Channel[Ch].Raw_Sample[Bu]/BLCALC);
						}
						Baseline=min(BL1,BL2); //Baseline=BL1;
						
						for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
							hBaseline->Fill(Bu,CoAs.Channel[Ch].Raw_Sample[Bu]-Baseline);
					}
				}					
			}
		}
	}
	

	float BASELINE[NB_SAMPLES];
	for(int i=0;i<NB_SAMPLES;i++)
	{
		if(isCalmode)BASELINE[i]=0;//hBaseline->ProjectionY("",i+1,i+1)->GetMean();
		else BASELINE[i]=hBaseline->ProjectionY("",i+1,i+1)->GetMean();
	
	}

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
		
			if(co!=16 && co!=31)
			{
				double SIGNAL[NB_SAMPLES];
			
				float Baseline=0;
				float BL1=0;
				float BL2=0;
				CoAs.Channel[Ch].Raw_Sample[0]=(CoAs.Channel[Ch].Raw_Sample[1]+CoAs.Channel[Ch].Raw_Sample[2])/2.;
								
				for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
				{ 
					if(Bu<BLCALC) BL1+=(CoAs.Channel[Ch].Raw_Sample[Bu])/BLCALC;
					if(Bu>NB_SAMPLES-BLCALC-2) BL2+=(CoAs.Channel[Ch].Raw_Sample[Bu])/BLCALC;
				}
				Baseline=min(BL1,BL2); //Baseline=BL1;
			
			
				for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
				{
					SIGNAL[Bu]= CoAs.Channel[Ch].Raw_Sample[Bu]-BASELINE[Bu]-Baseline-BASELINE[Bu];
				}
				
			
				double SIGNALFFT[NB_SAMPLES];
				double SIGNALDECONVFFT[NB_SAMPLES];
				GFFTReal(SIGNAL,SIGNALFFT,NB_SAMPLES,false);
				
				int fc=60; // cut frequency of the low-pass filter
				double Qf=1/2.; // quality factor
				for(unsigned short Bu=0;Bu<NB_SAMPLES;Bu+=2)
				{
					if(Bu<2)
					{
						SIGNALDECONVFFT[Bu]=SIGNALFFT[Bu]/ResponseFunction[co][as][ag][ch][Bu];
						SIGNALDECONVFFT[Bu+1]=SIGNALFFT[Bu+1]/ResponseFunction[co][as][ag][ch][Bu+1];
					}
					else
					{
						std::complex <double> sig(SIGNALFFT[Bu],SIGNALFFT[Bu+1]);
						std::complex <double> resp(ResponseFunction[co][as][ag][ch][Bu],ResponseFunction[co][as][ag][ch][Bu+1]);
						std::complex <double> deconv=sig/resp;
						
						std::complex<double> filter2(1-pow(Bu/2./fc,2.),Bu/2./fc/Qf);
						std::complex<double> final=deconv/pow(filter2,8); // 16th order filter...
			
						SIGNALDECONVFFT[Bu]=final.real();
						SIGNALDECONVFFT[Bu+1]=final.imag();
					}
				}
			
				GFFTReal(SIGNALDECONVFFT,SIGNAL,NB_SAMPLES,true);
		
				double ChargeMax=0;
				unsigned short TChargeMax=0;
				
				for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
				{	
           			if((((SIGNAL[Bu]>THR_PAD_DECONV && THR_PAD_DECONV>0) || (THR_PAD_DECONV<=0) ) && !isCalmode ) || FullTreeOpt)
					{
						DataReduced.peakheight.push_back(SIGNAL[Bu]);
						DataReduced.peaktime.push_back(Bu);
						//if(co==1 && as==0 && ag==0 && ch==0) cout << SIGNAL[Bu] << " " ;
					}
					else if(isCalmode && !FullTreeOpt && SIGNAL[Bu]>ChargeMax)
					{
						ChargeMax=SIGNAL[Bu];
						TChargeMax=Bu;
					}	
				}
				if(isCalmode && !FullTreeOpt)
				{
					DataReduced.peakheight.push_back(ChargeMax);
					DataReduced.peaktime.push_back(TChargeMax);
				}
				//if(co==1 && as==0 && ag==0 && ch==0) cout << endl;
			}


			else if(co==31)
			{
				for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
				{
					DataReduced.peakheight.push_back((CoAs.Channel[Ch].Raw_Sample[Bu]<<16)>>16);
					DataReduced.peaktime.push_back(CoAs.Channel[Ch].Raw_Sample[Bu]>>16);
				}
			}

			else if(co==16)
			{          
				TrapezoidAnalyser analyser(313,35,25);
				if((as==0 && ag==3) || (as==3 && ag==0)) analyser.SetDecayTime(125);
				analyser.Analyse(CoAs.Channel[Ch].Raw_Sample,NB_SAMPLES);
				double signal = analyser.GetSignal();
				DataReduced.peakheight.push_back((float)(abs(signal)));
				DataReduced.peaktime.push_back(1);
			}			
		
			if((co<16 && ch!=11 && ch!=22 && ch!=45 && ch!=56 && !FullTreeOpt) || co>=16 || FullTreeOpt) ReducedEvent.CoboAsad.push_back(DataReduced);
		}
	}


				
// ////////////////////////////////////////////////////
// //      Signal fit if 512 samples             //////				
// ////////////////////////////////////////////////////
// 				if(!FullTreeOpt && co!=31 && co!=16 && CoAs.Channel[Ch].Sample_Number.size()==NB_SAMPLES && !IsFastPeak)
// 				{	
// 					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
// 					{
// 						if(CoAs.Channel[Ch].Raw_Sample[Bu]<4095)
// 						{
// 							/*if(hasCalBaseline)*/ Charge[Bu]=(CoAs.Channel[Ch].Raw_Sample[Bu]-BaselineCal[co][as][ag][ch][Bu]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
// 							/*else Charge[Bu]=(CoAs.Channel[Ch].Raw_Sample[Bu]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];*/
// 						}
// 						else
// 						{
// 							Charge[Bu]=4096;
// 							DataReduced.hasSaturation=true;
// 						}
// 						Gsig->SetPoint(Bu+1,Bu,Charge[Bu]);
// 					}
// 					
// 					int Bu0=0;
// 					while(Bu0<NB_SAMPLES-1)
// 					{
// 						float ChargeMax=0;
// 						if(DataReduced.hasSaturation) THR=1000.;
//             
// 						for(unsigned short Bu=Bu0;Bu<NB_SAMPLES;Bu++)
// 						{	
// 							Bu0=Bu;
// 							if(Charge[Bu]>ChargeMax) ChargeMax=Charge[Bu];
// // cout << DataReduced.globalchannelid << " " << Bu << " " << Charge[Bu] << " " << ChargeMax << "    " << THR<< endl;
// 							if(Charge[Bu]<ChargeMax && ChargeMax>THR && ChargeMax<4096)
// 							{
// 								ffitSig->SetParameter(0,ChargeMax);
// 								ffitSig->SetParLimits(0,(ChargeMax)*0.8,(ChargeMax)*1.2);
// 								ffitSig->SetParameter(1,Bu-143.);
// 								ffitSig->SetParLimits(1,min((Bu-143.)*0.8,(Bu-143.)*1.2),max((Bu-143.)*0.8,(Bu-143.)*1.2));
// 
// 								Gsig->Fit(ffitSig,"RQ","",Bu-15,Bu+10);
// 								DataReduced.peakheight.push_back(ffitSig->GetParameter(0));
// 								DataReduced.peaktime.push_back(ffitSig->GetParameter(1)+143.);
// 								//cout << "FIT " << DataReduced.globalchannelid << endl;
// 					
// 								for(int bu=0;bu<NB_SAMPLES;bu++)
// 								{
// 									Charge[bu]-=ffitSig->Eval(bu);
// 									Gsig->SetPoint(bu+1,bu,Charge[bu]);
// 								}
// 								Bu=NB_SAMPLES;; 
// 							}
// 							else if(Charge[Bu]<ChargeMax && ChargeMax>THR && ChargeMax==4096)
// 							{
// 								DataReduced.peakheight.push_back(4096);
// 								DataReduced.peaktime.push_back(Bu0);
// 								for(int bu=0;bu<NB_SAMPLES;bu++)
// 								{
// 									if(bu<Bu0 || Charge[bu]>=THR) Charge[bu]=0;
// 									else bu=NB_SAMPLES;
// 									Gsig->SetPoint(bu+1,bu,Charge[bu]);
// 								}
// 								Bu=NB_SAMPLES;
// 							}
// 						}
// 					}
// 				}
// 
// ///////////////////////////////////////////////
// //    Peak find if  FastPeak option ON   //////				
// ///////////////////////////////////////////////
// 				else if(!FullTreeOpt && co!=31 && co!=16 && IsFastPeak && (CoAs.Channel[Ch].Sample_Number.size()==NB_SAMPLES))
// 				{
// 					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
// 					{	
// 						float charge;
// 						/*if(hasCalBaseline)*/ 
// 						  charge=(CoAs.Channel[Ch].Raw_Sample[Bu]-BaselineCal[co][as][ag][ch][Bu]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
// 						/*else
// 						  charge=(CoAs.Channel[Ch].Raw_Sample[Bu]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];*/
//            
//             						
// 						if(Bu>1 && Bu<NB_SAMPLES-1) // pad plane: take only the maxs
// 						{
// 							float prev1_charge=(CoAs.Channel[Ch].Raw_Sample[Bu-1]-BaselineCal[co][as][ag][ch][Bu-1]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
// 							float prev2_charge=(CoAs.Channel[Ch].Raw_Sample[Bu-2]-BaselineCal[co][as][ag][ch][Bu-2]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
// 							if(charge>P_CHARGE && charge>THR_PAD && charge > prev1_charge && charge > prev2_charge)
// 								if(!NMAX || (NMAX && charge > (CoAs.Channel[Ch].Raw_Sample[Bu-MIN_INT]-BaselineCal[co][as][ag][ch][Bu-MIN_INT]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1]))
// 								{
// 									P_CHARGE=charge;
// 									P_TIME=Bu;
// 									NSUP=0;
// 								}
// 				
// 							if(charge<P_CHARGE && P_CHARGE>0) NSUP++;
// 							if(NSUP==MIN_INT)
// 							{
// 								NMAX++;
// 	
// 								DataReduced.peakheight.push_back(P_CHARGE);
// 								DataReduced.peaktime.push_back(P_TIME);
// 								P_CHARGE=P_TIME=-1;
// 								NSUP=0;
// 							}
// 						}
// 					}
// 				}
//         
//         
// /////////////////////////////////////////////////
// //    Peak find if zero-suppress enabled   //////				
// /////////////////////////////////////////////////
//         
//  				else if(!FullTreeOpt && co!=31 && co!=16 && (CoAs.Channel[Ch].Sample_Number.size()!=NB_SAMPLES))
// 				{
//           			if(!hasCalBaseline)
//           			{
//            				Baseline=0;
// 						for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
//            					BaselineCal[co][as][ag][ch][Bu]=0;
//           			}
//           
//           			for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
// 					{	
// 						float charge=(CoAs.Channel[Ch].Raw_Sample[Bu]-BaselineCal[co][as][ag][ch][Bu]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
//             						
// 						if(Bu>1 && Bu<NB_SAMPLES-2) // pad plane: take only the maxs
// 						{
// 							float prev1_charge=(CoAs.Channel[Ch].Raw_Sample[Bu-1]-BaselineCal[co][as][ag][ch][Bu-1]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
// 							float prev2_charge=(CoAs.Channel[Ch].Raw_Sample[Bu-2]-BaselineCal[co][as][ag][ch][Bu-2]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
// 							float next1_charge=(CoAs.Channel[Ch].Raw_Sample[Bu+1]-BaselineCal[co][as][ag][ch][Bu+1]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
//               float next2_charge=(CoAs.Channel[Ch].Raw_Sample[Bu+2]-BaselineCal[co][as][ag][ch][Bu+2]-Baseline)*CalibCoefs[co][as][ag][ch][0]+CalibCoefs[co][as][ag][ch][1];
// 		    					
// 							if(charge>THR_PAD && charge > prev1_charge && charge > prev2_charge)
// 								if(charge> next1_charge && charge> next2_charge)
// 								{	
// 									DataReduced.peakheight.push_back(charge);
// 									DataReduced.peaktime.push_back(Bu);
// 								}
// 						}
// 					}
// 				}
// 
// 
// 
// ////////////////////////////////////////////////////
// //   Full tree option ON: take all samples    //////				
// ////////////////////////////////////////////////////				
// 				else if(FullTreeOpt && co!=31 && co!=16)
// 				{
// 					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
// 					{
// 						DataReduced.peakheight.push_back(CoAs.Channel[Ch].Raw_Sample[Bu]);//-BaselineCal[co][as][ag][ch][Bu]-Baseline);
// 						DataReduced.peaktime.push_back(Bu);
// 					}
// 				}
// 
// ////////////////////////////////////////////////////
// //         VXI parameters treatment           //////				
// ////////////////////////////////////////////////////				
// 				else if(co==31)
// 				{
// 					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
// 					{
// 						DataReduced.peakheight.push_back((CoAs.Channel[Ch].Raw_Sample[Bu]<<16)>>16);
// 						DataReduced.peaktime.push_back(CoAs.Channel[Ch].Raw_Sample[Bu]>>16);
// 					}
// 				}
// 
// //////////////////////////////////////////
// //         DSSD treatment           //////				
// //////////////////////////////////////////				
// 				else if(co==16)
// 				{
//          // vector<int> asads = {2,2,0,2,0,2,0,3,1,0};
//          // vector<int> agets = {1,0,2,3,0,2,3,0,0,1};
//          // e780 setting: DSSD 0,1,2,3,4,5,8,9 with DecayTime=313 - DSSD 6,7 with DecayTime=125
//           
// 					TrapezoidAnalyser analyser(313,35,25);
// 					if((as==0 && ag==3) || (as==3 && ag==0)) analyser.SetDecayTime(125);
// 					analyser.Analyse(CoAs.Channel[Ch].Raw_Sample,NB_SAMPLES);
// 					double signal = analyser.GetSignal();
// 					DataReduced.peakheight.push_back((float)(abs(signal)));
// 					DataReduced.peaktime.push_back(1);
// 				}			
// 			  
//         
//         
//         
// 				ReducedEvent.CoboAsad.push_back(DataReduced);
// 			}
// 		}
// 	}
	
	
	delete hBaseline;
	
	return(0);
}

