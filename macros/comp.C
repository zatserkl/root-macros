#include <TROOT.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TCut.h>
#include <iostream>
using std::cout;     using std::endl;

//
// functor implementation
//
class CompareBulkies {
public:
   Int_t sci1;
   Double_t thres;
   TCut cut;
   std::string wname;
   CompareBulkies(Int_t sci1_=0, Double_t thres_=100., TCut cut_="", std::string wname_="bulkies")
      :sci1(sci1_)
       , thres(thres_)
       , cut(cut_)
       , wname(wname_)
   {}
   void operator() (Int_t sci1_) {
      sci1 = sci1_;
      operator()();
   }
   void operator() () 
   {
      TTree* T = (TTree*) gDirectory->Get("T");
      if (!T) {
         cout<< "Could not find tree \"T\" in the current directory" <<endl;
         return;
      }

      Int_t ndivx = 1;
      Int_t ndivy = 2;

      TCanvas* can = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(wname.c_str());
      if (!can) {
         can = new TCanvas(wname.c_str(),wname.c_str(), 0,0, ndivx*gStyle->GetCanvasDefW(), ndivy*gStyle->GetCanvasDefH());
      }
      can->Clear();
      can->Divide(ndivx,ndivy);

      can->cd(1);
      gPad->SetLogy();
      T->Draw(Form("PhCh%d",sci1),cut+Form("PhCh%d>%f&&PhCh%d<%f",sci1,thres,sci1+1,thres));
      can->cd(2);
      gPad->SetLogy();
      T->Draw(Form("PhCh%d",sci1+1),cut+Form("PhCh%d>%f",sci1+1,thres));
      can->cd(0);

      const char* pave_name = "dirname";
      TPaveText* pave = (TPaveText*) gPad->GetListOfPrimitives()->FindObject(pave_name);
      if (pave) delete pave;
      pave = new TPaveText();
      gPad->GetListOfPrimitives()->Add(pave);

      pave->AddText(gDirectory->GetName());
      pave->SetTextSize(0.024);

      gPad->Modified();
      gPad->Update();
   }
};

//
// functor implementation
//
class CompareBulkiesThree {
public:
   Int_t sci1;
   Double_t thres;
   TCut cut;
   std::string wname;
   CompareBulkiesThree(Int_t sci1_=0, Double_t thres_=100., TCut cut_="", std::string wname_="bulkies")
      :sci1(sci1_)
       , thres(thres_)
       , cut(cut_)
       , wname(wname_)
   {}
   void operator() (Int_t sci1_) {
      sci1 = sci1_;
      operator()();
   }
   void operator() (const TCut& cut_) {
      cut = cut_;
      operator()();
   }
   void operator()() 
   {
      TTree* T = (TTree*) gDirectory->Get("T");
      if (!T) {
         cout<< "Could not find tree \"T\" in the current directory" <<endl;
         return;
      }

      Int_t ndivx = 1;
      Int_t ndivy = 3;

      TCanvas* can = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(wname.c_str());
      if (!can) {
         can = new TCanvas(wname.c_str(),wname.c_str(), 0,0, ndivx*gStyle->GetCanvasDefW(), ndivy*gStyle->GetCanvasDefH());
      }
      can->Clear();
      can->Divide(ndivx,ndivy);

      can->cd(1);
      //-- gPad->SetLogy();
      T->Draw(Form("PhCh%d",sci1),cut+Form("PhCh%d>%f&&PhCh%d<%f",sci1,thres,sci1+1,thres));
      can->cd(2);
      //-- gPad->SetLogy();
      T->Draw(Form("PhCh%d",sci1+1),cut+Form("PhCh%d>%f&&PhCh%d<%f",sci1+1,thres,sci1+2,thres));
      can->cd(3);
      //-- gPad->SetLogy();
      T->Draw(Form("PhCh%d",sci1+2),cut+Form("PhCh%d>%f",sci1+2,thres));
      can->cd(0);

      const char* pave_name = "dirname";
      TPaveText* pave = (TPaveText*) gPad->GetListOfPrimitives()->FindObject(pave_name);
      if (pave) delete pave;
      pave = new TPaveText();
      gPad->GetListOfPrimitives()->Add(pave);

      pave->AddText(gDirectory->GetName());
      pave->SetTextSize(0.025);

      gPad->Modified();
      gPad->Update();
   }
};

//
// instance of CompareBulkies
//
CompareBulkiesThree comp;
