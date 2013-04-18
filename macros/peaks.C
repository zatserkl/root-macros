// How to find peaks in histograms
//root.cern.ch/root/htmldoc/examples/peaks.C.html

// Example to illustrate the peak finder (class TSpectrum).
// This script generates a random number of gaussian peaks
// on top of a linear background.
// The position of the peaks is found via TSpectrum and injected
// as initial values of parameters to make a global fit
// To execute this example, do
//  root > .x peaks.C  (generate 10 peaks by default)
//  root > .x peaks.C++ (use the compiler)
//  root > .x peaks.C++(30) (generates 30 peaks)
//
// Author: Rene Brun

#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom.h>
#include <TSpectrum.h>
#include <TVirtualFitter.h>

#include <TFile.h>
#include <TMath.h>
#include <iostream>

using std::cout;           using std::endl;
   
Int_t peaks_NPEAKS = 30;
Double_t fpeaks(Double_t *x, Double_t *par) {
   Double_t result = par[0] + par[1]*x[0];
   for (Int_t p=0;p<peaks_NPEAKS;p++) {
      Double_t norm  = par[3*p+2];
      Double_t mean  = par[3*p+3];
      Double_t sigma = par[3*p+4];
      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}

TH1* peaks(TH1* h, Float_t thres=0.2, Int_t np=10)
{
   peaks_NPEAKS = np;
//    TH1 *h = new TH1("h","test",500,0,1000);
//    //generate n peaks at random
   Double_t par[3000];
//    par[0] = 0.8;
//    par[1] = -0.6/1000;
   Int_t p;
//    for (p=0;p<peaks_NPEAKS;p++) {
//       par[3*p+2] = 1;
//       par[3*p+3] = 10+gRandom->Rndm()*980;
//       par[3*p+4] = 3+2*gRandom->Rndm();
//    }
//    TF1 *f = new TF1("f",fpeaks,0,1000,2+3*peaks_NPEAKS);
//    f->SetNpx(1000);
//    f->SetParameters(par);
   TCanvas *c1 = new TCanvas("n1_c1","n1_c1",10,10,1000,900);
   c1->Divide(1,2);
   c1->cd(1);
//    h->FillRandom("f",200000);
   h->Draw();
   TH1 *h2 = (TH1*)h->Clone("h2");
   
   // Use TSpectrum to find the peak candidates
   TSpectrum *s = new TSpectrum(2*peaks_NPEAKS);
   //-- Int_t nfound = s->Search(h,1,"new");
   Int_t nfound = s->Search(h,1,"new", thres);
   printf("Found %d candidate peaks to fit\n",nfound);
   c1->Update();
   c1->cd(2);

   // estimate linear background
   TF1 *fline = new TF1("fline","pol1",h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
   h->Fit("fline","qn");
   
   // Loop on all found peaks. Eliminate peaks at the background level
   par[0] = fline->GetParameter(0);
   par[1] = fline->GetParameter(1);
   peaks_NPEAKS = 0;
   Float_t *xpeaks = s->GetPositionX();
   for (p=0;p<nfound;p++) {
      Float_t xp = xpeaks[p];
      Int_t bin = h->GetXaxis()->FindBin(xp);
      Float_t yp = h->GetBinContent(bin);
      //-------------- if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;
      if (yp - fline->Eval(xp) < thres) continue;
      par[3*peaks_NPEAKS+2] = yp;
      par[3*peaks_NPEAKS+3] = xp;
      // par[3*peaks_NPEAKS+4] = 3;
      par[3*peaks_NPEAKS+4] = 1;
      //par[3*peaks_NPEAKS+4] = 2;
      peaks_NPEAKS++;
   }
   printf("Found %d useful peaks to fit\n",peaks_NPEAKS);
   printf("Now fitting: Be patient\n");
   TF1 *fit = new TF1("fit",fpeaks,h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),2+3*peaks_NPEAKS);
   fit->SetLineWidth(1);
   fit->SetLineColor(2);
   TVirtualFitter::Fitter(h2,10+3*peaks_NPEAKS); // we may have more than the default 25 parameters
   fit->SetParameters(par);
   fit->SetNpx(1000);
   for (int i=2; i<fit->GetNpar(); i+=3) {
      fit->SetParName(i, "ampl");
      fit->SetParName(i+1, "mean");
      fit->SetParName(i+2, "#sigma");
   }
//    for (int i=2; i<fit->GetNpar(); i+=3) {
//       cout<< i << " name = " << fit->GetParName(i) <<endl;
//       cout<< i+1 << " name = " << fit->GetParName(i+1) <<endl;
//       cout<< i+2 << " name = " << fit->GetParName(i+2) <<endl;
//    }
   cout<< "h2->Fit" <<endl;
   h2->Fit("fit");
   
   h2->SetName(h->GetName());
   new TCanvas();
   h2->Draw();

   return h2;
}

TH1* peaks(const char* ifname="slice_2048_event_763.root", const char* hname="h436505368_signal_763", Float_t thres=0.2, Int_t np=10)
{
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return 0;
   }
   TH1* h = (TH1*) ifile->Get(hname);
   return peaks(h, thres, np);
}

void Add2File(TH1* h, const char* ofname="peaks.root", bool recreate=0)
{
   if (!h) {
      cout<< "h == 0" <<endl;
      return;
   }
   
   TFile* ofile;
   if (recreate) ofile = TFile::Open(ofname, "recreate");
   else          ofile = TFile::Open(ofname, "update");
   h->Clone()->Write();
   ofile->Close();
}

