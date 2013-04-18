// zatserkl@fnal.gov

/*
  NB: inverse quotes ``!
  g++ `$ROOTSYS/bin/root-config --cflags --glibs` -o sp3 sp3.C
*/
#include <TRint.h>          // the only include you need

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>
#include <sstream>

using std::cout;     using std::endl;

//-- global pointer to the tree
TTree* t;

//
// namespace Tree
//
namespace Tree
{
   Int_t evt;
   Int_t t1;
   Int_t t2;
   Int_t t3;
   Int_t a1;
   Int_t a2;
   Int_t a3;
   
   void clear()
   {
      evt      = 0;
      t1       = 0;
      t2       = 0;
      t3       = 0;
      a1       = 0;
      a2       = 0;
      a3       = 0;
   }
   
   void book(TTree* tree) {
      tree->Branch("evt",         &evt,        "evt/I");
      tree->Branch("t1",          &t1,         "t1/I");
      tree->Branch("t2",          &t2,         "t2/I");
      tree->Branch("t3",          &t3,         "t3/I");
      tree->Branch("a1",          &a1,         "a1/I");
      tree->Branch("a2",          &a2,         "a2/I");
      tree->Branch("a3",          &a3,         "a3/I");
   }
   void connect(TTree* tree)                                // need for event-by-event analysis
   {   
      // connects tree buffers with variables to use for event-by-event analysis
      // tree->SetBranchAddress("evt",           &evt);
      tree->SetBranchAddress("t1",            &t1);
      tree->SetBranchAddress("t2",            &t2);
      tree->SetBranchAddress("t3",            &t3);
      tree->SetBranchAddress("a1",            &a1);
      tree->SetBranchAddress("a2",            &a2);
      tree->SetBranchAddress("a3",            &a3);
   }
}  // namespace Tree

TTree* sp3(const char* ifname="data.out", const char* tree_name="t")
{
   TTree* tree_old = (TTree*) gDirectory->Get(tree_name);
   if (tree_old) {
      delete tree_old;
   }

   TTree* tree = new TTree(tree_name,"3 channels tree");
   Tree::book(tree);

   //-- assign to global pointer
   t = tree;

   std::ifstream ifile(ifname);
   if (!ifile) {
      cout<< "***ERROR addrun: file not found: " << ifname <<endl;
      return 0;
   }

   Tree::clear();
   --Tree::evt;
   while (ifile
         >> Tree::t1
         >> Tree::t2
         >> Tree::t3
         >> Tree::a1
         >> Tree::a2
         >> Tree::a3
         )
   {
      ++Tree::evt;
      tree->Fill();
   }
   ifile.close();

   TCanvas* can = new TCanvas;
   // can->Divide(2,1);
   // // can->cd(1); tree->Draw("t1", "t1>0&&t1<6000"); htemp()->Fit("gaus");
   // can->cd(1); tree->Draw("t1", "t1>0&&t1<6000");

   // // can->cd(2); tree->Draw("a1", "a1>1"); htemp()->Fit("gaus");
   // can->cd(2); tree->Draw("a1", "a1>1");
   // // Double_t gmean1 = htemp()->GetFunction("gaus")->GetParameter(1);
   // // Double_t gsigma1 = htemp()->GetFunction("gaus")->GetParameter(2);
   // // Double_t Npe1 = ((gmean1 - ped1)/gsigma1)*((gmean1 - ped1)/gsigma1);
   // // //cout<< "--> Fit parameters: gmean0 = " << gmean0 << " gsigma0 = " << gsigma0 <<endl;

   can->Divide(3,2);
   can->cd(1); tree->Draw("t1", "t1>0&&t1<6000");
   can->cd(2); tree->Draw("t2", "t2>0&&t2<6000");
   can->cd(3); tree->Draw("t3", "t3>0&&t3<6000");
   can->cd(4); tree->Draw("a1", "a1>1");
   can->cd(5); tree->Draw("a2", "a2>1");
   can->cd(6); tree->Draw("a3", "a3>1");
   can->cd();

   return tree;
}

// void ampl(Int_t plotmin=0, Int_t plotmax=0, Int_t fitmin=0, Int_t fitmax=0, Double_t ped1=50, Double_t ped2=37, Double_t ped3=3)
// {
//    // TCanvas* can = new TCanvas(cname.c_str(), cname.c_str(), 1000,500);
//    // TCanvas* can = new TCanvas(cname.c_str(), cname.c_str(), 1000,300);
//    TCanvas* can = new TCanvas;
//    can->Divide(3,1);
//    TCut cut;
//    if (plotmin > 0) cut += Form("a1>%d", plotmin);
//    if (plotmax > 0) cut += Form("a1<%d", plotmax);
//    t->Draw("a1", cut);
//    // t->Draw("a1", cut); htemp()->Fit("gaus");
//    // Double_t gmean1 = htemp()->GetFunction("gaus")->GetParameter(1);
//    // Double_t gsigma1 = htemp()->GetFunction("gaus")->GetParameter(2);
//    // Double_t Npe1 = ((gmean1 - ped1)/gsigma1)*((gmean1 - ped1)/gsigma1);
//    // cout<< "Npe1 = " << Npe1 <<endl;
//    // //cout<< "--> Fit parameters: gmean0 = " << gmean0 << " gsigma0 = " << gsigma0 <<endl;
// }
// 
// void time(Int_t plotmin=0, Int_t plotmax=6000, Int_t fitmin=0, Int_t fitmax=0)
// {
//    // TCanvas* can = new TCanvas(cname.c_str(), cname.c_str(), 1000,500);
//    // TCanvas* can = new TCanvas(cname.c_str(), cname.c_str(), 1000,300);
//    TCanvas* can = new TCanvas;
//    // can->Divide(3,1);
//    TCut cut;
//    if (plotmin > -1) cut += Form("t1>%d", plotmin);
//    if (plotmax > -1) cut += Form("t1<%d", plotmax);
//    t->Draw("t1", cut);
//    // t->Draw("t1", cut); htemp()->Fit("gaus");
//    // Double_t gmean1 = htemp()->GetFunction("gaus")->GetParameter(1);
//    // Double_t gsigma1 = htemp()->GetFunction("gaus")->GetParameter(2);
//    // cout<< "--> Fit parameters: gmean0 = " << gmean0 << " gsigma0 = " << gsigma0 <<endl;
// }

#if !defined(__CINT__) && !defined(__MAKECINT__)
int main(int argc, char *argv[])                   // w/o int will be warning ~ "no type"
{
  TRint* theApp = new TRint("Rint", &argc, argv, 0, 0, 1);  // do not show splash screen

  // code starts here

  //const char runlist_fname[] = "runtree.lst";
  //const char ofname[] = "runtree.root";

  //tree_write(runlist_fname, ofname);
  sp3();
  
  // TH1F* h = new TH1F("h", "h", 100, -3, 3);

  // h->FillRandom("gaus", 10000);
  // 
  // cout<< "h->GetMean() = " << h->GetMean() << " h->GetRMS() = " << h->GetRMS() <<endl;

  // TCanvas* can = new TCanvas();
  // 
  // h->Draw();

  // work in command line mode
  theApp->Run();
  delete theApp;
}
#endif
