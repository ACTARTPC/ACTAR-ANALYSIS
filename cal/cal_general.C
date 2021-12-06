void cal_general()
{

	TFile* f=new TFile("Calib_CATS1Y.root","READ");
	
	FILE* fcal=fopen("Calib_CATS1Y.dat","w");
	TH2I* Padsummary=new TH2I(*(TH2I*)(f->Get("hCATS1Y_raw")));

	for(int i=0;i<28;i++)
	{
		
		for(int p=0;p<5;p++)
	}
}
