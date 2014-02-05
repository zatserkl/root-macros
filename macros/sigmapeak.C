#include <TROOT.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMath.h>
#include <TF1.h>
#include <TPaveText.h>

#include <iostream>

using std::cout;     using std::endl;

// utils.C stuff
// TPaveText* titmax(Double_t x1NDC=0, Double_t x2NDC=1., Double_t y1NDC=0.920, Double_t y2NDC=0.999);

Double_t Nsigma(Double_t y[], Int_t i1, Int_t i2, Double_t bkg_sigma, Double_t bkg_mean=0)
{
   // calculate sum between ii1 and ii2 (use actual values)
   Double_t sum = 0;
   Int_t n = 0;
   for (int i=i1; i<=i2; ++i) {
      n++;
      sum += y[i] - bkg_mean;
   }
   if (n == 0) {
      cout<< "Nsigmas: n == 0" <<endl;
      return 0;
   }
   Double_t sum_ratio = sum / (bkg_sigma*TMath::Sqrt(n));
   // cout<< "sum_ratio = " << sum_ratio <<endl;
   return sum_ratio;
}

void smooth5a(Int_t np, Double_t y[], Double_t ys[]) {
   ys[0] = (3.*y[0] + 2.*y[0+1] + 1.*y[0+2]) / 6.;
   ys[1] = (2.*y[1-1] + 3.*y[1] + 2.*y[1+1] + 1.*y[1+2]) / 8.;
   for (int i=2; i<np-2; ++i) {
      // smooth on 5 points
      ys[i] = (
            -3.*y[i-2] + 
            12.*y[i-1] +
            17.*y[i] +
            12.*y[i+1] +
            -3.*y[i+2]
            ) / 35.;
   }
   ys[np-2] = (2.*y[np-2+1] + 3.*y[np-2] + 2.*y[np-2-1] + 1.*y[np-2-2]) / 8.;
   ys[np-1] = (3.*y[np-1] + 2.*y[np-1-1] + 1.*y[np-1-2]) / 6.;
}

void smooth5g(Int_t np, Double_t y[], Double_t ys[]) {
   ys[0] = (3.*y[0] + 2.*y[0+1] + 1.*y[0+2]) / 6.;
   ys[1] = (2.*y[1-1] + 3.*y[1] + 2.*y[1+1] + 1.*y[1+2]) / 8.;
   for (int i=2; i<np-2; ++i) {
      // smooth on 5 points
      ys[i] = (
            1.*y[i-2] + 
            3.*y[i-1] +
            4.*y[i] +
            3.*y[i+1] +
            1.*y[i+2]
            ) / 12.;
   }
   ys[np-2] = (2.*y[np-2+1] + 3.*y[np-2] + 2.*y[np-2-1] + 1.*y[np-2-2]) / 8.;
   ys[np-1] = (3.*y[np-1] + 2.*y[np-1-1] + 1.*y[np-1-2]) / 6.;
}

void smooth5(Int_t np, Double_t y[], Double_t ys[]) {
   ys[0] = (3.*y[0] + 2.*y[0+1] + 1.*y[0+2]) / 6.;
   ys[1] = (2.*y[1-1] + 3.*y[1] + 2.*y[1+1] + 1.*y[1+2]) / 8.;
   for (int i=2; i<np-2; ++i) {
      // smooth on 5 points
      ys[i] = (
            1.*y[i-2] + 
            2.*y[i-1] +
            3.*y[i] +
            2.*y[i+1] +
            1.*y[i+2]
            ) / 9.;
   }
   ys[np-2] = (2.*y[np-2+1] + 3.*y[np-2] + 2.*y[np-2-1] + 1.*y[np-2-2]) / 8.;
   ys[np-1] = (3.*y[np-1] + 2.*y[np-1-1] + 1.*y[np-1-2]) / 6.;
}

void smooth7(Int_t np, Double_t y[], Double_t ys[]) {
   ys[0] = (3.*y[0] + 2.*y[0+1] + 1.*y[0+2]) / 6.;
   ys[1] = (2.*y[1-1] + 3.*y[1] + 2.*y[1+1] + 1.*y[1+2]) / 8.;
   ys[2] = (1.*y[2-2] + 2.*y[2-1] + 3.*y[2] + 2.*y[2+1] + 1*y[2+2]) / 9.;
   for (int i=3; i<np-3; ++i) {
      // smooth on 7 points
      ys[i] = (
            1.*y[i-3] + 
            2.*y[i-2] + 
            3.*y[i-1] +
            4.*y[i] +
            3.*y[i+1] +
            2.*y[i+2] +
            1.*y[i+3]
            ) / 16.;
   }
   ys[np-3] = (1.*y[np-5] + 2.*y[np-4] + 3.*y[np-3] + 2.*y[np-2] + 1.*y[np-1]) / 9.;
   ys[np-2] = (2.*y[np-2+1] + 3.*y[np-2] + 2.*y[np-2-1] + 1.*y[np-2-2]) / 8.;
   ys[np-1] = (3.*y[np-1] + 2.*y[np-1-1] + 1.*y[np-1-2]) / 6.;
}

void smooth7a(Int_t np, Double_t y[], Double_t ys[]) {
   ys[0] = (3.*y[0] + 2.*y[0+1] + 1.*y[0+2]) / 6.;
   ys[1] = (2.*y[1-1] + 3.*y[1] + 2.*y[1+1] + 1.*y[1+2]) / 8.;
   ys[2] = (1.*y[2-2] + 2.*y[2-1] + 3.*y[2] + 2.*y[2+1] + 1*y[2+2]) / 9.;
   for (int i=3; i<np-3; ++i) {
      // smooth on 7 points
      ys[i] = (
            -2.*y[i-3] + 
            3.*y[i-2] + 
            6.*y[i-1] +
            7.*y[i] +
            6.*y[i+1] +
            3.*y[i+2] +
            -2.*y[i+3]
            ) / 21.;
   }
   ys[np-3] = (1.*y[np-5] + 2.*y[np-4] + 3.*y[np-3] + 2.*y[np-2] + 1.*y[np-1]) / 9.;
   ys[np-2] = (2.*y[np-2+1] + 3.*y[np-2] + 2.*y[np-2-1] + 1.*y[np-2-2]) / 8.;
   ys[np-1] = (3.*y[np-1] + 2.*y[np-1-1] + 1.*y[np-1-2]) / 6.;
}

/*
root -l 'sigmapeak.C+(0.025)'
// 
Processing sigmapeak.C+(0.025)...
background only: bkg_nsigma_0_99 = -0.352803
sig_nsigma_0_59  = -0.608621
sig_nsigma_0_99  = 2.1472
sig_nsigma_60_99 = 4.14042
sig_nsigma_70_90 = 5.57689
sig_nsigma_75_85 = 8.056
*/
void sigmapeak(Double_t signal_area=0.025)
{
   Double_t bkg_mean = 0;
   Double_t bkg_sigma = 0.001;

   Double_t x[100];
   Double_t y[100];
   Int_t np = 100;
   
   TRandom3 rand;
   
   for (int i=0; i<np; ++i) {
      x[i] = i;
      y[i] = rand.Gaus(bkg_mean,bkg_sigma);
   }
   
   TGraph* gr = new TGraph(np,x,y);
   gr->SetNameTitle("gr",Form("#sigma = %0.1f",bkg_sigma));
   gr->SetMarkerStyle(7);
   gr->SetMarkerColor(46);
   new TCanvas();
   gr->Draw("ap");

   Double_t bkg_nsigma_0_99 = Nsigma(gr->GetY(), 0,99, bkg_sigma);
   cout<< "background only: bkg_nsigma_0_99 = " << bkg_nsigma_0_99 <<endl;

   // add signal

   Double_t signal_mean = 80;
   Double_t signal_sigma = 3;

   for (int i=0; i<gr->GetN(); ++i) {
      Double_t xx = (gr->GetX()[i] - signal_mean)/signal_sigma;
      Double_t arg = 0.5*xx*xx;
      Double_t exp = arg < 50? TMath::Exp(-arg): 0;
      Double_t signal = signal_area/(TMath::Sqrt(TMath::TwoPi())*signal_sigma) * exp;
      gr->SetPoint(i, gr->GetX()[i], gr->GetY()[i] + signal);
   }

   gr->SetTitle(Form("#sigma_{bkg} = %0.3f signal: area = %0.3f mean = %0.0f sigma = %0.1f", bkg_sigma,signal_area,signal_mean,signal_sigma));

   gr->Draw("apl");
   gr->Fit("gaus", "R", "", signal_mean - 5*signal_sigma, signal_mean+5*signal_sigma);
   Double_t fit_area = 2.5 * gr->GetFunction("gaus")->GetParameter("Constant") * gr->GetFunction("gaus")->GetParameter("Sigma");
   cout<< "Area under fitted gaussian = " << fit_area <<endl;

   // titmax();
   gPad->Modified();    // to create box (NB: the pad was not drawn yet at this point!)
   gPad->Update();
   TPaveText* tit = (TPaveText*)gPad->GetPrimitive("title");
   tit->SetX1NDC(0.);
   tit->SetX2NDC(1.);
   tit->SetY1NDC(0.9);
   tit->SetY2NDC(1.);
   gPad->Modified();    // to update the pad
   gPad->Update();

   Double_t sig_nsigma_0_59 = Nsigma(gr->GetY(), 0,59, bkg_sigma);
   cout<< "sig_nsigma_0_59  = " << sig_nsigma_0_59 <<endl;

   Double_t sig_nsigma_0_99 = Nsigma(gr->GetY(), 0,99, bkg_sigma);
   cout<< "sig_nsigma_0_99  = " << sig_nsigma_0_99 <<endl;

   Double_t sig_nsigma_60_99 = Nsigma(gr->GetY(), 60,99, bkg_sigma);
   cout<< "sig_nsigma_60_99 = " << sig_nsigma_60_99 <<endl;

   Double_t sig_nsigma_70_90 = Nsigma(gr->GetY(), 70,90, bkg_sigma);
   cout<< "sig_nsigma_70_90 = " << sig_nsigma_70_90 <<endl;

   Double_t sig_nsigma_75_85 = Nsigma(gr->GetY(), 75,85, bkg_sigma);
   cout<< "sig_nsigma_75_85 = " << sig_nsigma_75_85 <<endl;

   Double_t ys5[100];
   smooth5(np, gr->GetY(), ys5);
   TGraph* gr5 = new TGraph(np, x, ys5);
   gr5->SetNameTitle("gr5","smoothed on 5 points");
   gr5->SetMarkerStyle(7);
   gr5->SetMarkerColor(2);
   gr5->SetLineColor(2);
   new TCanvas;
   gr5->Draw("apl");

   Double_t ys7[100];
   smooth7(np, gr->GetY(), ys7);
   TGraph* gr7 = new TGraph(np, x, ys7);
   gr7->SetNameTitle("gr7","smoothed on 7 points");
   gr7->SetMarkerStyle(7);
   gr7->SetMarkerColor(8);
   gr7->SetLineColor(8);
   new TCanvas;
   gr7->Draw("apl");

   Double_t ys7a[100];
   smooth7a(np, gr->GetY(), ys7a);
   TGraph* gr7a = new TGraph(np, x, ys7a);
   gr7a->SetNameTitle("gr7a","smoothed on 7a points");
   gr7a->SetMarkerStyle(7);
   gr7a->SetMarkerColor(3);
   gr7a->SetLineColor(3);
   new TCanvas;
   gr7a->Draw("apl");

   Double_t ys5g[100];
   smooth5g(np, gr->GetY(), ys5g);
   TGraph* gr5g = new TGraph(np, x, ys5g);
   gr5g->SetNameTitle("gr5g","smoothed on 5g points");
   gr5g->SetMarkerStyle(7);
   gr5g->SetMarkerColor(4);
   gr5g->SetLineColor(4);
   new TCanvas;
   gr5g->Draw("apl");

   Double_t ys5a[100];
   smooth5a(np, gr->GetY(), ys5a);
   TGraph* gr5a = new TGraph(np, x, ys5a);
   gr5a->SetNameTitle("gr5a","smoothed on 5a points");
   gr5a->SetMarkerStyle(7);
   gr5a->SetMarkerColor(6);
   gr5a->SetLineColor(6);
   new TCanvas;
   gr5a->Draw("apl");
}
