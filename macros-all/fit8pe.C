#include <TROOT.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>

#include <iostream>

using std::cout;     using std::endl;

#include "peakchan.C"

TF1* fpulse(Double_t xmin, Double_t xmax
      , Double_t A, Double_t x0, Double_t tau1, Double_t tau2, Double_t T, Double_t sigma
      , const char* name="fpulse"
      );

void fit8pe(Int_t iRC=0, Double_t sigma=-1.)
{
   const char* ifname = "ketek50-100Hz-lf_100pF_330ohm-26.8V-1.dat.root.peaks.root-sel_0.0320-0.0360.root";
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }

   TTree* speaks = (TTree*) ifile->Get("speaks");
   if (!speaks) {
      cout<< "Could not find tree 'speaks'" <<endl;
      return;
   }

   Peak* peak = 0;
   speaks->SetBranchAddress("speak", &peak);

   speaks->GetEntry(1);
   Peak* peak1 = new Peak(*peak);
   new TCanvas;
   peak1->Plot(iRC, sigma);
   // TGraphErrors* g1 = peak1->Plot();
   // TF1* f1 = fpulse(65,81, 0.2, 70, 1e4, 2, 0, 0.8, "f1");
   // f1->FixParameter(2, 1e4);
   // f1->FixParameter(4, 0);
   // g1->Fit(f1, "R", "", 69,80);

   speaks->GetEntry(2);
   Peak* peak2 = new Peak(*peak);
   new TCanvas;
   peak2->Plot(iRC, sigma);
   // TGraphErrors* g2 = peak2->Plot();
   // TF1* f2 = fpulse(200,218, 0.2, 205, 1e4, 2, 0, 0.8, "f2");
   // f2->FixParameter(2, 1e4);
   // f2->FixParameter(4, 0);
   // g2->Fit(f2, "R", "", 203,215);

   speaks->GetEntry(3);
   Peak* peak3 = new Peak(*peak);
   new TCanvas;
   peak3->Plot(iRC, sigma);
   // TGraphErrors* g3 = peak3->Plot();
   // TF1* f3 = fpulse(100,120, 0.2, 107, 1e4, 2, 0, 0.8, "f3");
   // f3->FixParameter(2, 1e4);
   // f3->FixParameter(4, 0);
   // g3->Fit(f3, "R", "", 103,115);

   speaks->GetEntry(4);
   Peak* peak4 = new Peak(*peak);
   new TCanvas;
   peak4->Plot(iRC, sigma);
   // TGraphErrors* g4 = peak4->Plot();
   // TF1* f4 = fpulse(150,170, 0.2, 157, 1e4, 2, 0, 0.8, "f4");
   // f4->FixParameter(2, 1e4);
   // f4->FixParameter(4, 0);
   // g4->Fit(f4, "R", "", 155,165);
}
