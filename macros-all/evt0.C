#include <TROOT.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TCut.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TTimer.h>

#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>

using std::cout;     using std::endl;

// utils.C stuff
TGraph* gtempClone(const char* name=0, const char* title=0);

TMultiGraph* evt(Int_t ch1=-1, Int_t ch2=-1, TCut cut="", Int_t evtNo=1)
{
   if (ch1==-1 || ch2==-1) {
      cout<< "Draws:    y vs x for tree with name \"t\" for channels ch1 and ch2" <<endl;
      cout<< "Usage:    evt(Int_t ch1=-1, Int_t ch2=-1, TCut cut=\"\", Int_t evtNo=1)" <<endl;
      cout<< "Example:  evt(1,2, \"x>0&&x<20\", 9)" <<endl;
      return 0;
   }

   TTree *t = (TTree*) gDirectory->Get("t");
   if (!t) {
      cout<< "No tree \"t\" found" <<endl;
      return 0;
   }

   Long64_t nselected = 0;
   nselected = t->Draw("y:x",cut+Form("ch==%d &&evt==%d",ch1,evtNo),"");
   if (nselected == 0) {
      cout<< "--> ch = " << ch1 << ": nselected == 0" <<endl;
      return 0;
   }
   TGraph* gr1 = gtempClone("gr1");
   gr1->SetMarkerStyle(7);
   gr1->SetMarkerColor(2);
   gr1->SetLineColor(2);

   nselected = t->Draw("y:x",cut+Form("ch==%d &&evt==%d",ch2,evtNo),"");
   if (nselected == 0) {
      cout<< "--> ch = " << ch2 << ": nselected == 0" <<endl;
      return 0;
   }
   TGraph* gr2 = gtempClone("gr2");
   gr2->SetMarkerStyle(7);
   gr2->SetMarkerColor(4);
   gr2->SetLineColor(4);

   if (gROOT->GetListOfCanvases()->Last()) gPad->Clear();
   TMultiGraph* mg = new TMultiGraph();
   mg->SetTitle(Form("evt %d, ch %d and %d %s", evtNo,ch1,ch2,cut.GetTitle()));
   mg->Add(gr1,"pl");
   mg->Add(gr2,"pl");
   mg->Draw("aw");

   gPad->Modified();
   gPad->Update();

   return mg;
}

void evtloop(Int_t ch1=-1, Int_t ch2=-1, TCut cut="", Int_t evtNo=1)
{
   if (ch1==-1 || ch2==-1) {
      cout<< "Draws:    (in loop) y vs x for tree with name \"t\" for channels ch1 and ch2" <<endl;
      cout<< "Usage:    evtloop(Int_t ch1=-1, Int_t ch2=-1, TCut cut=\"\", Int_t evtNo=0)" <<endl;
      cout<< "Example:  evtloop(1,2, \"x>0&&x<20\", 9)" <<endl;
      return;
   }

   TTimer timer("gSystem->ProcessEvents();",50,kFALSE);  // process mouse events every 50 ms

   Int_t ievtNo = evtNo;
   while (ievtNo >= 0)
   {
      TMultiGraph *mg = evt(ch1,ch2,cut,ievtNo);
      if (!mg) return;

      timer.Start();    // start processing of mouse events on TCanvas

      ++ievtNo;
      // cout<< "Number: show entry, <CR>: show next evt " << ievtNo << ", N: new TCanvas,  Q: stop, !command: execute command   ";
      cout<< "evt (<CR>: next evt " << ievtNo << ", N: new TCanvas,  Q: stop, !command: execute command) ";
      std::string line;
      std::getline(cin, line);

      timer.Stop();     // disable timer after <CR>

      if (line.size() == 0) continue;                                                           // <CR>
      if (line.find("!") != std::string::npos) {                                                // !command
         line.erase(line.find("!"), 1);
         gROOT->ProcessLine(line.c_str());
         --ievtNo;
         continue;
      }
      if (line.find("N") != std::string::npos || line.find("n") != std::string::npos) {         // N or n
         new TCanvas;
         --ievtNo;
         continue;
      }
      if (line.find("Q") != std::string::npos || line.find("q") != std::string::npos) break;    // Q or q
      std::stringstream ss(line);                                                               // number
      Int_t number = -1;
      if (ss >> number) ievtNo = number;
      else break;
   }
}
