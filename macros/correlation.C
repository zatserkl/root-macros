#include <TROOT.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TGraph.h>
#include <TMath.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>

#include <iostream>
#include <sstream>

using std::cout;     using std::endl;

Int_t nbi(Int_t n=0);
TH1* htemp(const char* name=0, const char* title=0);
void fitgaus(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void fitpol(Double_t xmin=0, Double_t xmax=0, Int_t power=1, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void left();
void right();

/*
root -l
.L correlation.C
corr(10)
*/

// to use with TTree: correlation(t->GetSelectedRows(),t->GetV1(),t->GetV2())
Double_t correlation(Int_t np, Double_t u[], Double_t v[])
{
   // sigma12 = sqrt(sigma1*sigma1 + sigma2*sigma2 - 2*rho*sigma1*sigma2

   // mean
   Double_t umean = 0;
   Double_t vmean = 0;
   for (int i=0; i<np; ++i) {
      umean += u[i];
      vmean += v[i];
   }
   umean /= np;
   vmean /= np;
   
   // variance as an average deviation from the mean (no averaging, really...)
   Double_t uVar = 0;
   Double_t vVar = 0;
   for (int i=0; i<np; ++i) {
      Double_t dmean;
      dmean = u[i] - umean;
      uVar += dmean*dmean;
      dmean = v[i] - vmean;
      vVar += dmean*dmean;
   }
   Double_t csum = 0;
   for (int i=0; i<np; ++i) csum += (u[i] - umean) * (v[i] - vmean);
   Double_t corr = csum / TMath::Sqrt(uVar*vVar);
   return corr;
}

TGraph* corr(Float_t gsig=10)
{
   TRandom3 rand;
   Double_t u[100];
   Double_t v[100];
   
   Int_t np = 100;
   
   for (int i=0; i<np; ++i) {
      u[i] = i;
      v[i] = rand.Gaus(i,gsig);
   }
   
   TGraph* gr = new TGraph(np,u,v);
   gr->SetName("gr");
   TCanvas* can = new TCanvas();
   gr->Draw("aw*");

   Float_t corr = correlation(np,u,v);
   cout<< "correlation coefficient = " << corr <<endl;
   
   std::stringstream ss;
   ss.str("");
   ss << "correlation coefficient = " << corr;
   gr->SetTitle(ss.str().c_str());
   can->Modified();
   
   return gr;
}

/*
root -l
.L correlation.C
sigmasum()
*/
TTree* sigmasum()
{
   // model experiment:
   // light source with Poisson(Gaussian) distributed number of photons
   // two photodetectors with quantum efficiencies eff1 and eff2

   Float_t nphot;
   Float_t npe1, npe2;

   TTree *tpe = new TTree("tpe", "Poisson model");
   tpe->SetMarkerStyle(7);
   tpe->SetMarkerColor(2);
   tpe->Branch("nphot", &nphot, "nphot/F");
   tpe->Branch("npe1", &npe1, "npe1/F");
   tpe->Branch("npe2", &npe2, "npe2/F");

   Double_t eff1 = 0.6;
   Double_t eff2 = 0.5;

   Double_t nphot_mean = 1000;

   Int_t nevents = 10000;
   TRandom3 rand;

   for (int ievent=0; ievent<nevents; ++ievent)
   {
      // generate the number of photons
      nphot = rand.Poisson(nphot_mean);

      // the number of photoelectrons in the detectors
      npe1 = rand.Poisson(eff1*nphot);
      npe2 = rand.Poisson(eff2*nphot);

      tpe->Fill();
   }

   // analysis

   Int_t nbins = nbi();
   nbi(40);

   new TCanvas;
   tpe->Draw("npe1", "");
   fitgaus(); right();
   Double_t sigma1 = htemp()->GetFunction("gaus")->GetParameter("Sigma");
   new TCanvas;
   tpe->Draw("npe2", "");
   fitgaus(); right();
   Double_t sigma2 = htemp()->GetFunction("gaus")->GetParameter("Sigma");
   new TCanvas;
   tpe->Draw("npe1-npe2", "");
   fitgaus(); right();
   Double_t sigma12 = htemp()->GetFunction("gaus")->GetParameter("Sigma");
   new TCanvas;
   tpe->Draw("npe2:npe1","");
   fitpol();
   Double_t rho = correlation(tpe->GetSelectedRows(),tpe->GetV1(),tpe->GetV2());

   Double_t sigma12_calc = TMath::Sqrt(sigma1*sigma1 + sigma2*sigma2 - 2*rho*sigma1*sigma2);

   cout<<endl;
   cout<< "sigma1 = " << sigma1 <<endl;
   cout<< "sigma2 = " << sigma2 <<endl;
   cout<< "sigma12 = " << sigma12 <<endl;
   cout<< "correlation coefficient = " << rho <<endl;
   cout<< "sigma12_calc = " << sigma12_calc <<endl;
   cout<<endl;

   nbi(nbins);
   return tpe;
}
