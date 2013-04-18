#include <TROOT.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>

#include <iostream>

using std::cout;     using std::endl;

// // utils.C stuff
// void fitpol(Double_t xmin=0, Double_t xmax=0, Int_t power=1, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);

/*
root -l nfitpoints.C+
*/
void nfitpoints()
{
   gStyle->SetOptFit();

   Double_t x3[1000], y3[1000], ex3[1000], ey3[1000];
   Int_t n3 = 160;
   Double_t x1[1000], y1[1000], ex1[1000], ey1[1000];
   Int_t n1 = 40;
   Double_t x2[1000], y2[1000], ex2[1000], ey2[1000];
   Int_t n2 = 10;

   TRandom3 random;
   Double_t x0, delta;
   Double_t dy = 10;
   Double_t sigma = 10;

   // parameters of line y = k*x + b
   Double_t k;
   Double_t b;
   // errors on parameters
   Double_t dk;
   Double_t db;
   // intersection with x-axis
   Double_t x;
   Double_t dx;

   TF1 *fitfunction = 0;

   delta = 100./n3;
   x0 = 0 - delta;
   for (int i=0; i<n3; ++i) {
      x0 += delta;
      x3[i] = x0;
      //-- y3[i] = random.Gaus(x3[i], 0.1*x3[i]);
      y3[i] = random.Gaus(x3[i], sigma);
      ex3[i] = 0;
      // ey3[i] = TMath::Sqrt(y3[i]);
      ey3[i] = dy;
   }

   TGraphErrors *gr3 = new TGraphErrors(n3, x3,y3, ex3, ey3);
   gr3->SetNameTitle("gr3", "gr3");
   gr3->SetMarkerStyle(20);
   gr3->SetMarkerColor(8);
   gr3->SetLineColor(8);

   new TCanvas;
   gr3->Draw("ap");
   // fitpol();
   gr3->Fit("pol1");

   fitfunction = gr3->GetFunction("pol1");
   b = fitfunction->GetParameter(0);
   db = fitfunction->GetParError(0);
   k = fitfunction->GetParameter(1);
   dk = fitfunction->GetParError(1);
   x = -b / k;
   dx = TMath::Abs(b/k) * TMath::Sqrt( (db/b)*(db/b) + (dk/k)*(dk/k) );
   gr3->SetTitle(Form("intersection with x-axis: %0.3f #pm %0.3f", x,dx));
   cout<< "--> intersection with x-axis: " << x << " +- " << dx <<endl<<endl;

   delta = 100./n1;
   x0 = 0 - delta;
   for (int i=0; i<n1; ++i) {
      x0 += delta;
      x1[i] = x0;
      // y1[i] = random.Gaus(x1[i], 0.1*x1[i]);
      y1[i] = random.Gaus(x1[i], sigma);
      ex1[i] = 0;
      // ey1[i] = TMath::Sqrt(y1[i]);
      ey1[i] = dy;
   }

   TGraphErrors *gr1 = new TGraphErrors(n1, x1,y1, ex1, ey1);
   gr1->SetNameTitle("gr1", "gr1");
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerColor(2);
   gr1->SetLineColor(2);

   new TCanvas;
   gr1->Draw("ap");
   // fitpol();
   gr1->Fit("pol1");

   fitfunction = gr1->GetFunction("pol1");
   b = fitfunction->GetParameter(0);
   db = fitfunction->GetParError(0);
   k = fitfunction->GetParameter(1);
   dk = fitfunction->GetParError(1);
   x = -b / k;
   dx = TMath::Abs(b/k) * TMath::Sqrt( (db/b)*(db/b) + (dk/k)*(dk/k) );
   gr1->SetTitle(Form("intersection with x-axis: %0.3f #pm %0.3f", x,dx));
   cout<< "--> intersection with x-axis: " << x << " +- " << dx <<endl<<endl;

   delta = 100./n2;
   x0 = 0 - delta;
   for (int i=0; i<n2; ++i) {
      x0 += delta;
      x2[i] = x0;
      //-- y2[i] = random.Gaus(x2[i], 0.1*x2[i]);
      y2[i] = random.Gaus(x2[i], sigma);
      ex2[i] = 0;
      // ey2[i] = TMath::Sqrt(y2[i]);
      ey2[i] = dy;
   }

   TGraphErrors *gr2 = new TGraphErrors(n2, x2,y2, ex2, ey2);
   gr2->SetNameTitle("gr2", "gr2");
   gr2->SetMarkerStyle(20);
   gr2->SetMarkerColor(4);
   gr2->SetLineColor(4);

   new TCanvas;
   gr2->Draw("ap");
   // fitpol();
   gr2->Fit("pol1");

   fitfunction = gr2->GetFunction("pol1");
   b = fitfunction->GetParameter(0);
   db = fitfunction->GetParError(0);
   k = fitfunction->GetParameter(1);
   dk = fitfunction->GetParError(1);
   x = -b / k;
   dx = TMath::Abs(b/k) * TMath::Sqrt( (db/b)*(db/b) + (dk/k)*(dk/k) );
   gr2->SetTitle(Form("intersection with x-axis: %0.3f #pm %0.3f", x,dx));
   cout<< "--> intersection with x-axis: " << x << " +- " << dx <<endl<<endl;
}
