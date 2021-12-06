#include <iostream>
#include <math.h>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3S.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <MTreeStructureOld.h>
#include <Parameters.h>

#include <Utils.h>
#include <MTrack.h>
#include <MEventReduced.h>
#include <MEvent.h>
#include <DataParameters.h>


int TreeReadOld(int run)
{
	
	char path[512];
// 	cin >> path;

// 	TFile treeFile(Form("/data/get_testX/test_get/acquisition/run/Root_6Li_OK/Run_%d.root",atoi(argv[1])));
	TFile treeFile(Form("/run/media/roger/ACTARTPC3/12C/recorrelated/Run_%d.root",run));
	TTree* tree = (TTree*)treeFile.Get("Tree");
		
// 	TApplication *theApp = new TApplication("ROOT example", &argc, argv, NULL, 0);
	
	MEventOld* S=new MEventOld();
	TBranch* event=tree->GetBranch("data");
	event->SetAddress(&S);
	
	char line[256];
	int TABLE[7][NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL];
	FILE* fTable=fopen("../../actar_analysis_IPNO/dat/LookupTableSi.dat","r");
	fgets(line,256,fTable);
	for(int i=0;i<NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL;i++)
		fscanf(fTable,"%d %d %d %d %d %d %d",&TABLE[0][i],&TABLE[1][i],&TABLE[2][i],&TABLE[3][i],&TABLE[4][i],&TABLE[5][i],&TABLE[6][i]);
	fclose(fTable);
	

	int globalasadnumber;
	int cobo, asad, aget, channel;
	int row, col;


	TCanvas* C=new TCanvas("C","C",1200,600);
	C->Divide(3,2);
	
	TH3S* visu_voxel=new TH3S("visu_voxel","visu_voxel",64,0,64,32,0,32,512,0,512);
	TH2S* visu_channel=new TH2S("visu_channel","visu_channel",NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL,0,NB_COBO*NB_ASAD*NB_AGET*NB_CHANNEL,1000,0,4096);
	
	TGraph* G01=new TGraph(10000);
     TGraph* G02=new TGraph(10000);
	TGraph* G03=new TGraph(10000);
	TGraph* G12=new TGraph(10000);
	TGraph* G13=new TGraph(10000);
	TGraph* G23=new TGraph(10000);
	long int TS[NB_COBO];
	
	std::cout << "\nRE-ORDERED FRAGMENTS READ FROM TREE with " << tree->GetEntries() <<  " entries" << endl;

// 
	int prev_event=-1;
	double QTOTSUM=0;
	int counter=0;
	
	const int nRow=32;
	const int nCol=64;
	float PAD[nRow][nCol];
	float TIME[nRow][nCol];
// 	MTrack* T1=new MTrack();
	TH1F* hThetaXY=new TH1F("hThetaXY","hThetaXY",600,-90,90);
	TH1F* hThetaXZ=new TH1F("hThetaXZ","hThetaXZ",600,-90,90);
	TH1F* hThetaYZ=new TH1F("hThetaYZ","hThetaYZ",600,-1,95);
	TH2F* visu_charge=new TH2F("visu_charge","visu_charge",64,0,64,32,0,32);


	MEventReduced* ReducedEvent=new MEventReduced();
	TString TreeNameFile;
// 	TreeNameFile.Form("/run/media/roger/ACTARTPC4/Tree_Run_%04d_Merged.root",run_number);
	TreeNameFile.Form("/run/media/roger/ACTARTPC3/12C/NewFormat/Tree_Run_%04d.root",run);
	TFile* TreeFile=new TFile(TreeNameFile.Data(),"recreate");		
	TTree* MergedTree=new TTree("ACTAR_TTree","2nd level Tree",0);
	MergedTree->Branch("data","MEventReduced",&ReducedEvent,16000,0);
// 	TBranch* event=tree->GetBranch("data");
// 	event->SetAddress(&S);
	
	for(int e = 0/*tree->GetEntries()-10000*/; e<tree->GetEntries(); ++e)
	{
		tree->GetEntry(e);
		
		
		if(prev_event!=-1 && prev_event>S->event_number) cout << "problem in event merging! with " << prev_event << " " << S->event_number << " ----> " << prev_event-S->event_number << endl;
		
		ReducedEvent->CoboAsad.clear();
		
		for(int f=0;f<S->Fragments.size();f++)
		{
			
			cobo=(int)(S->Fragments[f].globalasadnumber/NB_ASAD);
			asad=(int)(S->Fragments[f].globalasadnumber%NB_ASAD);
			TS[cobo]=S->Fragments[f].time_stamp;			
			
			for(int c=0;c<S->Fragments[f].Channels.size();c++)
			{
				aget=(int)(S->Fragments[f].Channels[c].globalchannelid/NB_CHANNEL);
				channel=(int)(S->Fragments[f].Channels[c].globalchannelid%NB_CHANNEL);
				row=TABLE[5][cobo*NB_ASAD*NB_AGET*NB_CHANNEL + asad*NB_AGET*NB_CHANNEL + aget*NB_CHANNEL + channel];
				col=TABLE[6][cobo*NB_ASAD*NB_AGET*NB_CHANNEL + asad*NB_AGET*NB_CHANNEL + aget*NB_CHANNEL + channel];
				short Qmax=0;
				short Tmax=0;
				
				ReducedEvent->event=S->event_number;
				ReducedEvent->timestamp=S->Fragments[f].time_stamp;
				
				ReducedData DataReduced;
				DataReduced.globalchannelid=channel+(aget<<7)+(asad<<9)+(cobo<<11);

				for(int p=0;p<S->Fragments[f].Channels[c].peaknumber.size();p++)
				{
					if(Qmax<S->Fragments[f].Channels[c].peakheight[p])
					{
						Qmax=S->Fragments[f].Channels[c].peakheight[p];
						Tmax=S->Fragments[f].Channels[c].peaktime[p];//*VDRIFT;
					}
					DataReduced.peakheight.push_back(S->Fragments[f].Channels[c].peakheight[p]);
					DataReduced.peaktime.push_back(S->Fragments[f].Channels[c].peaktime[p]);
				}
				
				ReducedEvent->CoboAsad.push_back(DataReduced);				
				
// 				PAD[row][col]=Qmax;
// 				TIME[row][col]=Tmax;
// 				visu_charge->Fill(col,row,Qmax);
// 				
// 				QTOTSUM+=Qmax;
// 				visu_channel->Fill(cobo*NB_ASAD*NB_AGET*NB_CHANNEL + asad*NB_AGET*NB_CHANNEL + aget*NB_CHANNEL + channel,Qmax);
			}
		}
		
		MergedTree->Fill();
		
		
// 		T1->zx_e=NPADX-1;
// 		T1->zx_s=1;
// 		T1->zy_s=0;
// 		T1->zy_e=13;//NPADY-1;
	
// 		FitMat3D(PAD,TIME,T1->zy_s,T1->zy_e,T1->zx_s,T1->zx_e,8.,T1);
// 		hThetaXY->Fill(R2D*atan((T1->Yh-T1->Ym)/(T1->Xh-T1->Xm)));
// 		hThetaXZ->Fill(R2D*atan((T1->Zh-T1->Zm)/(T1->Xh-T1->Xm)));


// 		visu_charge->Draw("colz");
// 		C->WaitPrimitive();
		
// 		for(int r=0;r<nRow;r++)
// 			for(int c=0;c<nCol;c++)
// 				PAD[r][c]=TIME[r][c]=0;
// 		visu_charge->Reset();
		
		
// 		G01->SetPoint(counter,counter,TS[0]-TS[1]);
// 		G02->SetPoint(counter,counter,TS[0]-TS[2]);
// 		G03->SetPoint(counter,counter,TS[0]-TS[3]);
// 		G12->SetPoint(counter,counter,TS[1]-TS[2]);
// 		G13->SetPoint(counter,counter,TS[1]-TS[3]);
// 		G23->SetPoint(counter,counter,TS[2]-TS[3]);
// 		
		counter++;
		prev_event=S->event_number;
				
	}
	
	MergedTree->Write();
	TreeFile->Close();

	
	C->cd(1);
	G01->Draw("AP");
	C->cd(2);
	G02->Draw("AP");
	C->cd(3);
	G03->Draw("AP");
	C->cd(4);
	G12->Draw("AP");
	C->cd(5);
	G13->Draw("AP");
	C->cd(6);
	G23->Draw("AP");
	
// 	C->SaveAs(Form("figs/CheckRun%d.png",atoi(argv[1])));
	
	cout << "QTOTSUM is " << QTOTSUM << endl;
	
	
	TCanvas* CC=new TCanvas("CC","CC",800,800);
	CC->Divide(2,2);
	CC->cd(1);
	hThetaXY->Draw();
	CC->cd(2);
	hThetaXZ->Draw();

// 	visu_channel->Draw("colz");
	
	// Run interactive interface
// 	theApp->Run();
	return(0);
}

