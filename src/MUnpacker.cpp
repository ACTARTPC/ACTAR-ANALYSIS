#include <iostream>
#include <MUnpacker.h>
#include <MEvent.h>

using namespace std;

MUnpacker::MUnpacker()
{
	insideframe=new MFMCommonFrame();
	mergeframe=new MFMMergeFrame();
	coboframe=new MFMCoboFrame();
	mutantframe=new MFMMutantFrame();
	ebyedatframe=new MFMEbyedatFrame();
}

MUnpacker::~MUnpacker()
{
	delete insideframe;
	delete mergeframe;
	delete coboframe;
	delete mutantframe;
	delete ebyedatframe;
}


long int MUnpacker::Unpack(MFMCommonFrame* frame, MEvent* Event)
{
	frame->SetAttributs();
	int type=frame->GetFrameType();
	int headersize = frame->GetHeaderSize();
//  	cout << "frame type: " << type << " and size " << frame->GetFrameSizeAttribut() << endl;
		
	if((type == MFM_MERGE_EN_FRAME_TYPE) || (type == MFM_MERGE_TS_FRAME_TYPE))
	{
		mergeframe->SetAttributs(frame->GetPointHeader());
		
		int nbinsideframe=mergeframe->GetNbItems();
		Event->EN = mergeframe->GetEventNumber();
		Event->TS = mergeframe->GetTimeStamp();
		mergeframe->ResetReadInMem();
// 		cout << "Merge frame EN: " << mergeframe->GetEventNumber() << " with " << nbinsideframe << " inside frames and size " << mergeframe->GetFrameSize() << endl;

		for(int i1=0;i1<nbinsideframe;i1++)
		{
			mergeframe->ReadInFrame(insideframe);
// 			cout << "--- " << i1 << "  ";
			Unpack(insideframe,Event);
		}
		//if(Event->CoboAsad.size()==nbinsideframe) return(Event->EN);
		//else return(-1);
    		return(nbinsideframe);
	}
	
	else if ((type==MFM_COBO_FRAME_TYPE) || (type==MFM_COBOF_FRAME_TYPE))
	{	
		coboframe->SetAttributs(frame->GetPointHeader());
		
		nbitems=coboframe->GetNbItems();
		
		MCoboAsad TheCoboAsad;

		TheCoboAsad.cobo_number=coboframe->CoboGetCoboIdx();
		TheCoboAsad.asad_number=coboframe->CoboGetAsaIdx();
		TheCoboAsad.global_asad_number=TheCoboAsad.cobo_number*COBO_NB_ASAD+TheCoboAsad.asad_number;
		TheCoboAsad.EN=Event->EN=coboframe->GetEventNumber();
		TheCoboAsad.TS=Event->TS=coboframe->GetTimeStamp();
    
    if(coboframe->CoboGetAsaIdx()==0) Event->INFO.FrameCounter[coboframe->CoboGetCoboIdx()]++;

		for(int ag=0;ag<4;ag++)
		{
 			for(int j=0;j<4;j++)
 				if(((int)(*(coboframe->CoboGetHitPat(ag)))>>j)&1)
 					TheCoboAsad.hit_pattern.push_back(ag*NB_CHANNEL + 3-j);
 			for(int i=1;i<9;i++)
 				for(int j=0;j<8;j++)
 					if(((int)(*(coboframe->CoboGetHitPat(ag)+i))>>j)&1)
 						TheCoboAsad.hit_pattern.push_back(ag*NB_CHANNEL + (i-1)*8+11-j);
			
			TheCoboAsad.last_cell[ag]=((*((uint16_t*)(coboframe->CoboGetLastCell(ag))))&511);
			TheCoboAsad.multiplicity[ag]=*((uint16_t*)(coboframe->CoboGetMultip(ag)));
// 			cout << TheCoboAsad.EN << " -  ";
// 			for(int ih=0;ih<2;ih++) for(int jh=0;jh<8;jh++) cout <<(((int)(*(coboframe->CoboGetLastCell(ag)+ih))>>jh)&1) << " " ;
// 			cout << "     " << TheCoboAsad.last_cell[ag] << endl;
		}	
		
		
		short iChan=0;
		short iBuck=0;
		short iAget=0;
		
		
		for(int i2=0;i2<nbitems;i2++)
		{
			coboframe->CoboGetParameters(i2,&sample, &buckidx,&chanidx,&agetidx);

			if(agetidx>3 || chanidx>67 || buckidx>511)
			{
				if(!CheckAndRepairCorruptedFrame(frame,TheCoboAsad)) Event->stat_bad_frame++;
				else Event->stat_recovered_frame++;
				break;			
			}
			
			else if(type==MFM_COBO_FRAME_TYPE)
			{
				Event->stat_good_frame++;
				TheCoboAsad.Channel[agetidx*NB_CHANNEL+chanidx].Raw_Sample[buckidx]=sample;
				TheCoboAsad.Channel[agetidx*NB_CHANNEL+chanidx].Sample_Number.push_back(buckidx);
				//if(!buckidx) TheCoboAsad.hit_pattern.push_back(agetidx*NB_CHANNEL+chanidx);
			}

			else if(type==MFM_COBOF_FRAME_TYPE)
			{							
				Event->stat_good_frame++;
				TheCoboAsad.Channel[agetidx*NB_CHANNEL+iChan].Raw_Sample[iBuck]=sample;
				TheCoboAsad.Channel[agetidx*NB_CHANNEL+iChan].Sample_Number.push_back(iBuck);
				if(!iBuck) TheCoboAsad.hit_pattern.push_back(agetidx*NB_CHANNEL+iChan);
				iChan++;
				if(i2%2==1)
				{
					iAget++;
					iChan-=2;
				}
				if(iAget>=NB_AGET)
				{
					iAget=0;
					iChan+=2;
				}
				if(iChan>=NB_CHANNEL)
				{
					iBuck++;
					iChan=0;
				}
			}
		}
		
		Event->CoboAsad.push_back(TheCoboAsad);		

		return(TheCoboAsad.EN);
	}
	
	else if(type == MFM_MUTANT_FRAME_TYPE)
	{
		mutantframe->SetAttributs(frame->GetPointHeader());
//  		event = mutantframe->GetEventNumber();
// 		timestamp = mutantframe->GetTimeStamp();
		return(mutantframe->GetEventNumber());
		
	}
	
	else if(type == MFM_EBY_EN_TS_FRAME_TYPE)
	{
		ebyedatframe->SetAttributs(frame->GetPointHeader());
		
		MCoboAsad TheCoboAsad;
		
		TheCoboAsad.hit_pattern.push_back(0);
		
		TheCoboAsad.cobo_number=31;
		TheCoboAsad.asad_number=TheCoboAsad.global_asad_number=0;
		TheCoboAsad.global_asad_number=TheCoboAsad.cobo_number*COBO_NB_ASAD+TheCoboAsad.asad_number;
		TheCoboAsad.EN=Event->EN=ebyedatframe->GetEventNumber();
		TheCoboAsad.TS=Event->TS=ebyedatframe->GetTimeStamp();
		Event->INFO.FrameCounter[16]++;

		uint16_t label, value;
		
		int isImplantOK=0;
		
		for(int i=0;i<ebyedatframe->GetNbItems();i++)
		{
			ebyedatframe->EbyedatGetParameters(i,&label,&value);

			TheCoboAsad.Channel[(int)(i/NB_SAMPLES)].Raw_Sample[i%NB_SAMPLES]=value+(label<<16);
			TheCoboAsad.Channel[(int)(i/NB_SAMPLES)].Sample_Number.push_back(i%NB_SAMPLES); 
			if(Event->hasResponseFunction && ((label==Event->DeconvCondLabel && value<Event->DeconvCondValue) || Event->DeconvCondValue<0))
				Event->isDeconvOK=true;
			
			if(label==Event->ImplantCondLabel1 && value>Event->ImplantCondValueLow1 && value<Event->ImplantCondValueHigh1) isImplantOK+=1;
			if(label==Event->ImplantCondLabel2 && value>Event->ImplantCondValueLow2 && value<Event->ImplantCondValueHigh2) isImplantOK+=2;
		}
		
		if(isImplantOK==3)
		{
			if(Event->implantCounter%10==0) Event->ImplantationExtract=true;
			Event->implantCounter++;
		}
		
		Event->CoboAsad.push_back(TheCoboAsad);
		
		return(ebyedatframe->GetEventNumber());
	}

	
	else if (type != 0xff00 && type != 0xff11 && type != 0x7)
	{
		printf("\nUnknown frame of type %d (0x%x) and size %d (0x%x)\n",type,type,frame->GetFrameSize(),frame->GetFrameSize());
		frame->HeaderDisplay(NULL);
		frame->DumpRaw(128,0);
		cout << endl;
		return(-666);
	}
	
	else return(0);
}


bool MUnpacker::CheckAndRepairCorruptedFrame(MFMCommonFrame* frame, MCoboAsad& TheCoboAsad)
{
	std::vector <short> aget_channel_present={};
	coboframe->CoboGetParameters(0,&sample, &buckidx,&chanidx,&agetidx);
	TheCoboAsad.Channel[agetidx*NB_CHANNEL+chanidx].Raw_Sample[buckidx]=sample;
	TheCoboAsad.Channel[agetidx*NB_CHANNEL+chanidx].Sample_Number.push_back(buckidx);
	aget_channel_present.push_back(agetidx*NB_CHANNEL+chanidx);
	unsigned short first_bucket=buckidx;
	int index=0;
	for(int i2=1;i2<nbitems;i2++)
	{
		coboframe->CoboGetParameters(i2,&sample, &buckidx,&chanidx,&agetidx);
		index++;
		if(buckidx==first_bucket) aget_channel_present.push_back(agetidx*NB_CHANNEL+chanidx);
		if(chanidx>67) chanidx=aget_channel_present[index%aget_channel_present.size()]%NB_CHANNEL;
		else if(agetidx>3) return(false);
		else if(buckidx>511) return(false);
		TheCoboAsad.Channel[agetidx*NB_CHANNEL+chanidx].Raw_Sample[buckidx]=sample;
		TheCoboAsad.Channel[agetidx*NB_CHANNEL+chanidx].Sample_Number.push_back(buckidx);
	}
	return(true);
}

