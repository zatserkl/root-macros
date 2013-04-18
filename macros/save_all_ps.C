void save_all_ps(const char psfile[] = "xxx.ps")
{
// Int_t ps_type = 111;   // portrait  ps
Int_t ps_type = 112;      // landscape ps
// Int_t ps_type = 113;   // eps
TPostScript* postscript = new TPostScript(psfile, ps_type);
postscript->Clear();
TSeqCollection* clist =  gROOT->GetListOfCanvases();
TCanvas* canvas;
for(int i=0;i<clist->GetSize();i++){canvas = (TCanvas*) clist->At(i); postscript->NewPage(); canvas->Modified(1); canvas->Update();}
postscript->Close();
cout<< "Output file: " << psfile <<endl;
}
