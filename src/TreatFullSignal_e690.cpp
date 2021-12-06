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




std::vector<int> readalwayslist_b= {1703,1707,1888,1890,1898,1902,2262,2264,2431,2433,2267,2269};
std::vector<int> readalwayslist_t= {2142,2146,2150,2155,2165,2167,2707,2710,3230,3232,3234,3236};
std::vector<int> readalwayslist= {1703,1707,1888,1890,1898,1902,2262,2264,2431,2433,2267,2269 , 2142,2146,2150,2155,2165,2167,2707,2710,3230,3232,3234,3236};


int MEvent::TreatFullSignal(bool FullTreeOpt)
{
	int   BLCALC=         30;
	
	float THR_PAD_DECONV= -0.5;
	int fc=60; // cut frequency of the low-pass filter
	double Qf=1/2.; // quality factor
	
	float THR=            50.;
	int   MIN_INT=        50;
	
	TH2F* hBaseline=new TH2F("hBaseline","hBaseline",NB_SAMPLES,0,NB_SAMPLES,4096,0,4096);
	
	ReducedEvent.event=EN;
	ReducedEvent.timestamp=TS;
	
// /**/	fprintf(f_out,"%d ",EN);
// /**/	int nch=0;
// 	cout << endl;

	int cptRAL=0;
	
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
				
				float tmp_max=0;
				bool isRAL_bOK=true;
				bool isRAL_tOK=true;
				for(int ind=0;ind<readalwayslist.size();ind++)
				{
					if( co*NB_ASAD*NB_AGET*NB_CHANNEL + as*NB_AGET*NB_CHANNEL + ag*NB_CHANNEL+ch == readalwayslist[ind])
					{
						cptRAL++;
						CoAs.Channel[Ch].Raw_Sample[0]=(CoAs.Channel[Ch].Raw_Sample[1]+CoAs.Channel[Ch].Raw_Sample[2])/2.;
						for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
						{ 
							if(Bu<BLCALC) BL1+=(CoAs.Channel[Ch].Raw_Sample[Bu]/BLCALC);
							if(Bu>NB_SAMPLES-BLCALC-2) BL2+=(CoAs.Channel[Ch].Raw_Sample[Bu]/BLCALC);
						}
						Baseline=min(BL1,BL2);

						for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
							if(CoAs.Channel[Ch].Raw_Sample[Bu]-Baseline>tmp_max) tmp_max=CoAs.Channel[Ch].Raw_Sample[Bu]-Baseline;
						if(tmp_max>THR)
						{
							if(ind<readalwayslist.size()/2.) isRAL_bOK=false;
							if(ind>=readalwayslist.size()/2.) isRAL_tOK=false;
						} 
					}
				}	
				
				if(isRAL_bOK) for(auto readalwayspad:readalwayslist_b) if(co*NB_ASAD*NB_AGET*NB_CHANNEL + as*NB_AGET*NB_CHANNEL + ag*NB_CHANNEL+ch == readalwayspad)
				{
					Baseline=BL1=BL2=0;
					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
					{ 
						if(Bu<BLCALC) BL1+=(CoAs.Channel[Ch].Raw_Sample[Bu]/BLCALC);
						if(Bu>NB_SAMPLES-BLCALC-2) BL2+=(CoAs.Channel[Ch].Raw_Sample[Bu]/BLCALC);
					}
					Baseline=min(BL1,BL2);
					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
					{
						hBaseline->Fill(Bu,CoAs.Channel[Ch].Raw_Sample[Bu]-Baseline);
					}
				}
				if(isRAL_tOK) for(auto readalwayspad:readalwayslist_t) if(co*NB_ASAD*NB_AGET*NB_CHANNEL + as*NB_AGET*NB_CHANNEL + ag*NB_CHANNEL+ch == readalwayspad)
				{
					Baseline=BL1=BL2=0;
					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
					{ 
						if(Bu<BLCALC) BL1+=(CoAs.Channel[Ch].Raw_Sample[Bu]/BLCALC);
						if(Bu>NB_SAMPLES-BLCALC-2) BL2+=(CoAs.Channel[Ch].Raw_Sample[Bu]/BLCALC);
					}
					Baseline=min(BL1,BL2);
					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
						hBaseline->Fill(Bu,CoAs.Channel[Ch].Raw_Sample[Bu]-Baseline);
				}
			}
		}
	}
	
	if(cptRAL!=readalwayslist.size())
	{
		cout << "ReadAlwaysList is not correct: " << cptRAL << "/" << readalwayslist.size() << " channels in event - skipping to NO BASELINE CORRECTION MODE" << endl;
		isBaselineCorr=false;
	}	
	
	float BASELINE[NB_SAMPLES];
	for(int i=0;i<NB_SAMPLES;i++)
	{
		if(isCalmode || !isBaselineCorr)BASELINE[i]=0;
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
// if(co==9 && as==1 && ag==2 && ch==64 && EN==52400)
// {
// 	cout << endl; 
// 	for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number) cout << CoAs.Channel[Ch].Raw_Sample[Bu] << " " ;
// 	cout << endl; 
// }
			
			if(co!=31)
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
					if(CoAs.Channel[Ch].Raw_Sample[Bu]>=4095) DataReduced.hasSaturation=true;
					SIGNAL[Bu]= CoAs.Channel[Ch].Raw_Sample[Bu]-Baseline-BASELINE[Bu];
				}
				
				
			
				double SIGNALFFT[NB_SAMPLES];
				double SIGNALDECONVFFT[NB_SAMPLES];
				GFFTReal(SIGNAL,SIGNALFFT,NB_SAMPLES,false);
				
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
				int max_sample_block=0;
				
				unsigned short bucket_start[(int)(NB_SAMPLES/2)]={0};
				unsigned short bucket_length[(int)(NB_SAMPLES/2)]={0};
				int cur_sample_block=0;
	
////////////////// NORMAL MODE //////////////////////////			
				if(!FullTreeOpt && ch!=11 && ch!=22 && ch!=45 && ch!=56)
				{
					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
					{	
						if(SIGNAL[Bu]>THR_PAD_DECONV && THR_PAD_DECONV>0)
						{
							if(bucket_start[cur_sample_block]==0)
								bucket_start[cur_sample_block]=Bu;
							bucket_length[cur_sample_block]++;
						}
						
						if(SIGNAL[Bu]>ChargeMax)
						{
							ChargeMax=SIGNAL[Bu];
							TChargeMax=Bu;
							max_sample_block=cur_sample_block;
						}
						if((SIGNAL[Bu]<=THR_PAD_DECONV && THR_PAD_DECONV>0) && Bu>0)
							if(SIGNAL[Bu-1]>THR_PAD_DECONV) cur_sample_block++;						
					}
				}
				if(ch!=11 && ch!=22 && ch!=45 && ch!=56)for(unsigned short Bu=bucket_start[max_sample_block];Bu<bucket_start[max_sample_block]+bucket_length[max_sample_block];Bu++)
				{
					DataReduced.peakheight.push_back(SIGNAL[Bu]);
					DataReduced.peaktime.push_back(Bu);
				}
				
////////////////// CALIB MODE ///////////////////////////				
				if(isCalmode && !FullTreeOpt)
				{
					DataReduced.peakheight.push_back(ChargeMax);
					DataReduced.peaktime.push_back(TChargeMax);
				}


////////////////// FULL READOUT MODE ////////////////////								
				if(((THR_PAD_DECONV<=0) && !isCalmode ) || FullTreeOpt)
					for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
					{
						DataReduced.peakheight.push_back(SIGNAL[Bu]);
						DataReduced.peaktime.push_back(Bu);
					}
			}


			else if(co==31)
			{
				for(unsigned short Bu : CoAs.Channel[Ch].Sample_Number)
				{
					DataReduced.peakheight.push_back((CoAs.Channel[Ch].Raw_Sample[Bu]<<16)>>16);
					DataReduced.peaktime.push_back(CoAs.Channel[Ch].Raw_Sample[Bu]>>16);
				}
			}

			if((co<31 && ch!=11 && ch!=22 && ch!=45 && ch!=56 && !FullTreeOpt) || co>=31 || FullTreeOpt) ReducedEvent.CoboAsad.push_back(DataReduced);
		}
	}


				

	
	
	delete hBaseline;
	
	return(0);
}

