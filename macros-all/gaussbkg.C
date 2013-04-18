#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include <iostream>
#include <cmath>

using std::cout;        using std::endl;

Double_t fsig(Double_t* xx, Double_t* par)
{
   Double_t x = *xx;
   Double_t a = par[0];
   Double_t mean = par[1];
   Double_t sigma = par[2];

   Double_t fsig = 0;
   Double_t arg = std::pow( (x-mean)/2/sigma, 2);
   if (std::abs(arg) > 50) return fsig;

   fsig = a*exp(-arg);
   return fsig;
}

Double_t fbkg(Double_t* xx, Double_t* par)
{
   *xx = *xx;
   return par[0];
}

Double_t fsigbkg(Double_t* xx, Double_t* par)
{
   Double_t fsigbkg = fbkg(xx,par) + fsig(xx, &par[1]);
   return fsigbkg;
}

void gaussbkg()
{
   TH1F* h = new TH1F("h","h",100,0,100);

   TRandom3 rand;           // default parameter is UInt_t seed = 4357
   Double_t sig_mean = 60;
   Double_t sig_sigma = 5;

   Int_t Ntot = 10000;
   for (int i=0; i<Ntot; ++i) {
      Int_t igauss = rand.Gaus(sig_mean,sig_sigma) + 0.5;
      h->Fill(igauss);
   }

   // background level
   Double_t bkg_mean = 100;
   Double_t bkg_sigma = 10;

   for (int i=1; i<=h->GetNbinsX(); ++i) {
      Int_t ibkg = rand.Gaus(bkg_mean, bkg_sigma) + 0.5;
      h->SetBinContent(i, ibkg + h->GetBinContent(i));
   }

   new TCanvas;
   h->Draw();

   Int_t npar_bkg = 1;
   Int_t npar_sig = 3;
   TF1* ffit = new TF1("ffit", fsigbkg, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(), npar_bkg+npar_sig);

   Double_t par[10];
   par[0] = 1000;                   // some problem with par[0] = 100 (which is bkg_mean)???
   Int_t ipar_sig = npar_bkg;

   Double_t a = h->GetMaximum();
   Double_t mean = (h->GetXaxis()->GetXmin() + h->GetXaxis()->GetXmax())/2;
   Double_t sigma = 1;
   par[ipar_sig++] = a;
   par[ipar_sig++] = mean;
   par[ipar_sig++] = sigma;
   ffit->SetParameters(par);
   cout<< "ipar_sig = " << ipar_sig <<endl;
   
   h->Fit(ffit);
}
