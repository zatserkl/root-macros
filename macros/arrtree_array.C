// zatserkl@fnal.gov

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>

using std::cout;     using std::endl;

//
// namespace Tree
//
namespace Tree
{
   Int_t    v[19];
   
   void clear()
   {
      for (int i=0; i<19; ++i) v[i] = 0;
   }
   
   void book(TTree* tree) {
      tree->Branch("v",          v,         "v[19]/I");     // NB: not &v. About variable length see User Guide, "Adding a Branch to Hold a List of Variables"
   }

   void connect(TTree* tree)
   {   
      // connects tree buffers with variables for event-by-event analysis
      tree->SetBranchAddress("v",            v);            // NB: not &v
   }
}  // namespace Tree

TTree* arrtree_array(const char* ifname="run_158.erik", const char* ofname="")
{
   //bool debug = true;
   bool debug = false;
   if (debug) cout<< "debug is on" <<endl;  // to use identificator debug
   
   if (*ofname == '\0') ofname = Form("%s.root", ifname);

   // create tree
   TFile* ofile = TFile::Open(ofname,"recreate");
   TTree* tree = new TTree("t","tree example");
   Tree::book(tree);
   
   std::ifstream ifile(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return 0;
   }

   // fill tree
   Tree::clear();
   while (
         ifile 
         >> Tree::v[0] >> Tree::v[1] >> Tree::v[2] >> Tree::v[3] >> Tree::v[4]
         >> Tree::v[5] >> Tree::v[6] >> Tree::v[7] >> Tree::v[8] >> Tree::v[9]
         >> Tree::v[10] >> Tree::v[11] >> Tree::v[12] >> Tree::v[13] >> Tree::v[14]
         >> Tree::v[15] >> Tree::v[16] >> Tree::v[17] >> Tree::v[18]
         )
   {
      tree->Fill();
   }
   
   ofile = tree->GetCurrentFile();  // update current file pointer
   ofile->Write();
   cout<< "Wrote tree " << tree->GetName() << " in file " << ofile->GetName() <<endl;
   return tree;
}

// void tree_draw(const char* ifname="arrtree.root")
// {
//    // uses TTree::Draw
//    
//    //bool debug = true;
//    bool debug = false;
//    if (debug) cout<< "debug is on" <<endl;  // to use identificator debug
//     
//    // arrange tree
//    TFile* ifile = TFile::Open(ifname);
//    TTree* tree  = (TTree*) ifile->Get("t");
// 
//    TH1F* h = new TH1F("h", "h", 100, 0, 100);
//    if (debug) cout<< "Fill histogram " << h->GetName() <<endl;  // to use h
//    
//    new TCanvas;
//    tree->Draw("a>>h", "");
// }

void tree_read(const char* ifname="run_158.erik.root")
{
   // calls Tree::connect(tree) to connect tree buffers with variables for event-by-event analysis
   
   //bool debug = true;
   bool debug = false;
   if (debug) cout<< "debug is on" <<endl;  // to use identificator debug
    
   // arrange tree
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return;
   }
   TTree* tree  = (TTree*) ifile->Get("t");
   Tree::connect(tree);

   TH1F* h_time = new TH1F("TAC1_158", "", 100, 1700, 2100);
   h_time->SetFillStyle(3001);
   h_time->SetFillColor(4);
 
   // for (int ievent=0; ievent<tree->GetEntries(); ++ievent)  // loop over events
   // {
   //    tree->GetEntry(ievent);
   //    h->Fill(Tree::v[11]);
   // }

   new TCanvas;
   gPad->SetGrid(0,0);
   gStyle->SetOptStat(0);

   tree->Project(h_time->GetName(), "v[11]", "v[11]>1700&&v[11]<2100");
   h_time->Draw();

   TF1* ftime = new TF1("ftime", "gaus", 1900, 1975);
   h_time->Fit(ftime, "R");

   TH1F* h_ampl = new TH1F("ADC1_158", "", 100, 0, 1000);
   h_ampl->SetFillStyle(3001);
   h_ampl->SetFillColor(4);

   tree->Project(h_ampl->GetName(), "v[9]", "v[9]<1000");

   new TCanvas;
   gPad->SetGrid(0,0);
   h_ampl->Draw();
}

#ifdef __CINT__
void rightgaus();
void leftgaus();

void NIM_plots(const char* ifname="run_158.erik.root")
{
   // calls Tree::connect(tree) to connect tree buffers with variables for event-by-event analysis
   
   //bool debug = true;
   bool debug = false;
   if (debug) cout<< "debug is on" <<endl;  // to use identificator debug
    
   // arrange tree
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return;
   }
   TTree* tree  = (TTree*) ifile->Get("t");
   Tree::connect(tree);

   TH1F* h_time = new TH1F("TAC1_158", "", 100, 1700, 2100);
   h_time->SetFillStyle(3001);
   h_time->SetFillColor(4);
 
   // for (int ievent=0; ievent<tree->GetEntries(); ++ievent)  // loop over events
   // {
   //    tree->GetEntry(ievent);
   //    h->Fill(Tree::v[11]);
   // }

   new TCanvas;
   gPad->SetGrid(0,0);
   gStyle->SetOptStat(0);

   tree->Project(h_time->GetName(), "v[11]", "v[11]>1700&&v[11]<2100");
   h_time->Draw();

   TF1* ftime = new TF1("ftime", "gaus", 1900, 1975);
   h_time->Fit(ftime, "R");
   gPad->Modified();
   gPad->Update();
   // leftgaus();                 // 'No stat box found' w/o Modified() and Update()
   rightgaus();                 // 'No stat box found' w/o Modified() and Update()

   TH1F* h_ampl = new TH1F("ADC1_158", "", 100, 0, 1000);
   h_ampl->SetFillStyle(3001);
   h_ampl->SetFillColor(4);

   tree->Project(h_ampl->GetName(), "v[9]", "v[9]<1000");

   new TCanvas;
   gPad->SetGrid(0,0);
   h_ampl->Draw();
}
#endif

/*
NB: inverse quotes ``!
g++ `$ROOTSYS/bin/root-config --cflags --glibs` arrtree.C
*/
#if !defined(__CINT__)
#include <TRint.h>         // include TRint.h to build stand-along ROOT program
#include <TROOT.h>
#include <cstdarg>

int main(int argc, char *argv[])                      // w/o int ACLiC warning ~ "no type"
{
   TRint* theApp = new TRint("Rint", 0, 0, 0, 0, 1);  // do not show splash screen

   const char* ifname = "run_158.erik";
   if (argc > 1) ifname = argv[1];
   const char* ofname = "";
   if (argc > 2) ofname = argv[2];

   // arrtree(ifname); // NB: CINT will not see the return value

   // work in command line mode
   cout<< argv[0] << ": load macro arrtree.C" <<endl;
   gROOT->LoadMacro("arrtree.C");                        // add "+" to compile with ACLiC
   if (*ifname != '\0') gROOT->ProcessLine(Form("arrtree(\"%s\", \"%s\")", ifname, ofname));
   else                 gROOT->ProcessLine("arrtree()");

   gROOT->ProcessLine("tree_read()");

   theApp->Run();
   delete theApp;
}
#endif
