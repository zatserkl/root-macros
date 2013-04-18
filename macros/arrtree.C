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
#include <string>
#include <sstream>

using std::cout;     using std::endl;

//
// namespace Tree
//
namespace Tree
{
   const Int_t MAXDIM = 100;
   Float_t v[MAXDIM];
   
   void clear() {
      for (int i=0; i<MAXDIM; ++i) v[i] = 0;
   }
   
   void book(TTree* tree, Int_t nbra=1) {
      for (int i=0; i<nbra; ++i) {
         tree->Branch(Form("v%d",i), &v[i], Form("v%d/F",i));
      }
   }
   void book(TTree* tree, const char* bnames="v") {
      std::stringstream ss(bnames);
      std::string bname;
      for (int i=0; ss >> bname; ++i) {
         // cout<< "book branch " << bname <<endl;
         tree->Branch(Form("%s",bname.c_str()), &v[i], Form("%s/F",bname.c_str()));
      }
   }

   void connect(TTree* tree) {   
      // connects tree buffers with variables for event-by-event analysis
      for (int ibra=0; ibra<tree->GetListOfBranches()->GetEntries(); ++ibra) {
         tree->SetBranchAddress(tree->GetListOfBranches()->At(ibra)->GetName(), &v[ibra]);
      }
   }
}  // namespace Tree

TTree* arrtree(const char* ifname="run_158.erik", const char* ofname="")
{
   //bool debug = true;
   bool debug = false;
   if (debug) cout<< "debug is on" <<endl;  // to use identificator debug
   
   if (*ofname == '\0') ofname = Form("%s.root", ifname);

   // create tree
   TFile* ofile = TFile::Open(ofname,"recreate");
   TTree* tree = new TTree("t","tree example");

   // Int_t ncol = 19;
   // Tree::book(tree, ncol);

   // Tree::book(tree, "v0 v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18");
   Tree::book(tree, "runType runNum blkCnt tuplNum scaler1 scaler2 scaler3 ramp1 ramp2 adc1 adc2 tac1 tac2 tac3 adc3 adc4 adc5 hrClock lrClock");
   
   std::ifstream ifile(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return 0;
   }
 
   Tree::clear();
   for (std::string line; std::getline(ifile,line);)
   {
      std::istringstream ss(line);
      for (int ibranch=0; ibranch<tree->GetListOfBranches()->GetEntries(); ++ibranch) {
         ss >> Tree::v[ibranch];
      }
      if (ss.rdstate() & ss.failbit) continue;

      // successfully read ncol numbers from the line
      if (debug) {
         for (int ibranch=0; ibranch<tree->GetListOfBranches()->GetEntries(); ++ibranch) {
            cout<< Tree::v[ibranch] << " ";
         }
         cout<<endl;
      }

      tree->Fill();
   }
   ifile.close();
   cout<< "Read " << tree->GetEntries() << " data lines" <<endl;
   
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
   //    h_time->Fill(Tree::v[11]);
   // }

   new TCanvas;
   gPad->SetGrid(0,0);
   gStyle->SetOptStat(0);

   tree->Project(h_time->GetName(), "tac1", "tac1>1700&&tac1<2100");
   h_time->Draw();

   TF1* ftime = new TF1("ftime", "gaus", 1900, 1975);
   h_time->Fit(ftime, "R");

   TH1F* h_ampl = new TH1F("ADC1_158", "", 100, 0, 1000);
   h_ampl->SetFillStyle(3001);
   h_ampl->SetFillColor(4);

   tree->Project(h_ampl->GetName(), "adc1", "adc1<1000");

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
   //    h->Fill(Tree::tac1);
   // }

   new TCanvas;
   gPad->SetGrid(0,0);
   gStyle->SetOptStat(0);

   tree->Project(h_time->GetName(), "tac1", "tac1>1700&&tac1<2100");
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

   tree->Project(h_ampl->GetName(), "adc1", "adc1<1000");

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
