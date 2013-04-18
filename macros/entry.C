#include <TROOT.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TCut.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TPad.h>

#include <iostream>

using std::cout;     using std::endl;

// utils.C stuff
TGraph* gtempClone(const char* name=0, const char* title=0);

TMultiGraph *entry(const char *y_vs_x="help", TCut case1="", TCut case2="", TCut common="", Int_t entryNo=0
      , Int_t marker=7, Int_t color1=2, Int_t color2=4)
{
   if (strcmp(y_vs_x,"help")==0) {
      cout<< "Draws:    y vs x for tree with name \"t\" for conditions case1 and case2" <<endl;
      cout<< "Usage:    entry(y_vs_x, TCut case1, TCut case2, TCut common, Int_t entryNo, Int_t marker=7, Int_t color1=2, Int_t color2=4)" <<endl;
      cout<< "Example:  entry(\"y:x\", \"ch==1\", \"ch==2\", \"x>0&&x<20\", 9)" <<endl;
      return 0;
   }

   TTree *t = (TTree*) gDirectory->Get("t");
   if (!t) {
      cout<< "No tree \"t\" found" <<endl;
      return 0;
   }

   Long64_t nselected = 0;
   nselected = t->Draw(y_vs_x, common+case1, "", 1, entryNo);
   if (nselected == 0) {
      cout<< "--> case1 " << case1.GetTitle() << ": nselected == 0" <<endl;
      return 0;
   }
   TGraph* gr1 = gtempClone("gr1");
   gr1->SetMarkerStyle(marker);
   gr1->SetMarkerColor(color1);
   gr1->SetLineColor(color1);

   nselected = t->Draw(y_vs_x, common+case2, "", 1, entryNo);
   if (nselected == 0) {
      cout<< "--> case2 " << case2.GetTitle() << ": nselected == 0" <<endl;
      return 0;
   }
   TGraph* gr2 = gtempClone("gr2");
   gr2->SetMarkerStyle(marker);
   gr2->SetMarkerColor(color2);
   gr2->SetLineColor(color2);

   if (gROOT->GetListOfCanvases()->Last()) gPad->Clear();
   TMultiGraph* mg = new TMultiGraph();
   mg->SetTitle(Form("entry %d for %s and %s   %s",entryNo,case1.GetTitle(),case2.GetTitle(),common.GetTitle()));
   mg->Add(gr1,"pl");
   mg->Add(gr2,"pl");
   mg->Draw("aw");

   gPad->Modified();
   gPad->Update();

   return mg;
}
