// zatserkl@fnal.gov

// NB: the interpetor fails to execute this macro

/*
   --> to compile with MathMore
g++ `$ROOTSYS/bin/root-config --cflags --glibs` -lMathMore -o gauss gauss.C
*/

#include <TRint.h>          // the only include you need to call ROOT
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>
#include "Math/GSLIntegrator.h"
#include "Math/WrappedTF1.h"

#include <iostream>
#include <cstdarg>

using std::cout;      using std::endl;

Double_t fgauss(Double_t* x, Double_t* par)
{
   // direct implementation of the gaussian in ROOT style

   Double_t arg = x[0];
   Double_t ampl = par[0];
   Double_t mean = par[1];
   Double_t sigma = par[2];

   Double_t res = 0;
   Double_t exparg = ((arg-mean)*(arg-mean)) / (2.*sigma*sigma);
   if (TMath::Abs(exparg) < 50) res = ampl*exp(-exparg);
   return res;
}

Double_t fgaussGSL(Double_t x, void* param)
{
   // direct implementation of the gaussian in GSL style

   Double_t arg = x;
   Double_t* par = (Double_t*) param;

   Double_t ampl = par[0];
   Double_t mean = par[1];
   Double_t sigma = par[2];

   Double_t res = 0;
   Double_t exparg = ((arg-mean)*(arg-mean)) / (2.*sigma*sigma);
   if (TMath::Abs(exparg) < 50) res = ampl*exp(-exparg);
   return res;
}

void gauss()
{
   Double_t xlow = 0;
   Double_t xhigh = 1024;
   // Int_t nbins = 1024;
   Int_t nbins = 256;
   TH1F* h = new TH1F("h", Form("h;channels;events/%0.0f channels",(xhigh-xlow)/nbins), nbins,xlow,xhigh);

   h->SetMarkerStyle(7);
   h->SetMarkerColor(8);
   h->SetLineColor(8);
   h->SetFillColor(8);

   // background
   TF1* fbkg = new TF1("fbkg", "gaus(0)", 0, 1024);   // gaus(0) is a substitute for [0]*exp(-0.5*((x-[1])/[2])**2)
   fbkg->SetParameter(0, 1);  // this parameter would be ignored, but it should not be 0?
   fbkg->SetParameter(1, 200);
   fbkg->SetParameter(2, 30);

   h->FillRandom("fbkg", 100000);

   new TCanvas;
   h->Draw();

   h->Fit("gaus");

   cout<< "h->GetNbinsX() = " << h->GetNbinsX() <<endl;
   cout<< "h->Integral() = " << h->Integral() <<endl;
   cout<< "h->Integral(\"width\") = " << h->Integral("width") <<endl;

   Double_t gampl = h->GetFunction("gaus")->GetParameter(0);
   Double_t gmean = h->GetFunction("gaus")->GetParameter(1);
   Double_t gsigma = h->GetFunction("gaus")->GetParameter(2);
   Double_t nevt = TMath::Sqrt(2.*TMath::Pi()) * gampl*gsigma;
   cout<< "TMath::Sqrt(2.*TMath::Pi()) * gampl*gsigma = " << nevt <<endl;
   cout<< "h->GetBinWidth(1) = " << h->GetBinWidth(1) <<endl; 
   Double_t nevt_w = TMath::Sqrt(2.*TMath::Pi()) * gampl*gsigma / h->GetBinWidth(1);
   cout<< "TMath::Sqrt(2.*TMath::Pi()) * gampl*gsigma / h->GetBinWidth(1) = " << nevt_w <<endl;

   // gaussian built out from fit parameters
   TF1* fun = new TF1("fun", "gaus(0)", 0, 1024);
   fun->SetLineColor(4);
   fun->SetLineStyle(2);
   fun->SetParameters(gampl, gmean, gsigma);
   fun->Draw("same");

   cout<< "fun->Integral(0,1024) = " << fun->Integral(0,1024) <<endl;

   // just to be sure: use our standalone function fgauss
   TF1* fungauss = new TF1("fungauss", fgauss, 0, 1024, 3);    // the number of pars = 3
   fungauss->SetLineStyle(10);
   fungauss->SetLineColor(2);
   fungauss->SetParameters(gampl, gmean, gsigma);
   fungauss->Draw("same");

   // integrate fungauss using GNU Scientific Library (GSL) support in MathMore
   ROOT::Math::GSLIntegrator gslIntegrator; 
   ROOT::Math::WrappedTF1 wrappedTF1(*fungauss); 
   gslIntegrator.SetFunction(wrappedTF1);
   double val2 = gslIntegrator.Integral(0,1024);
   cout<< "MathMore integral = " << val2 <<endl;

   // integrate "native" GSL function double fgaussGSL(double, void*)
   Double_t parGSL[3];
   parGSL[0] = gampl;
   parGSL[1] = gmean;
   parGSL[2] = gsigma;
   ROOT::Math::GSLIntegrator gslIntegratorGSL;
   gslIntegratorGSL.SetFunction(fgaussGSL, parGSL);
   double val3 = gslIntegratorGSL.Integral(0,1024);
   cout<< "MathMore direct integral = " << val3 <<endl;

   //-- TH1F* h = new TH1F("h", "h", 1024, 0, 1024);
   // root$ gcc `$ROOTSYS/bin/root-config --cflags --glibs` -lMathMore -o gauss gauss.C
   // root$ ./gauss
   // *-- Default rootlogon
   // *-- Local rootlogon
   //  FCN=206.133 FROM MIGRAD    STATUS=CONVERGED      65 CALLS          66 TOTAL
   //                      EDM=1.67583e-10    STRATEGY= 1      ERROR MATRIX ACCURATE 
   //   EXT PARAMETER                                   STEP         FIRST   
   //   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   //    1  Constant     1.33257e+03   5.16753e+00   2.96437e-02   3.95470e-06
   //    2  Mean         1.99969e+02   9.46714e-02   6.65234e-04  -5.31738e-05
   //    3  Sigma        2.98793e+01   6.70342e-02   4.28108e-06   2.34280e-02
   //-- h->GetNbinsX() = 1024
   // h->Integral() = 100000
   // h->Integral("width") = 100000
   //-- TMath::Sqrt(2.*TMath::Pi()) * gampl*gsigma = 99804.7
   // h->GetBinWidth(1) = 1
   //-- TMath::Sqrt(2.*TMath::Pi()) * gampl*gsigma / h->GetBinWidth(1) = 99804.7
   // fun->Integral(0,1024) = 99804.7
   //-- MathMore integral = 99804.7
   // MathMore direct integral = 99804.7

   //-- TH1F* h = new TH1F("h", "h", 256, 0, 1024);
   // root$ gcc `$ROOTSYS/bin/root-config --cflags --glibs` -lMathMore -o gauss gauss.C
   // root$ ./gauss
   // *-- Default rootlogon
   // *-- Local rootlogon
   // FCN=70.628 FROM MIGRAD    STATUS=CONVERGED      68 CALLS          69 TOTAL
   // EDM=2.0748e-15    STRATEGY= 1      ERROR MATRIX ACCURATE 
   // EXT PARAMETER                                   STEP         FIRST   
   // NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   // 1  Constant     5.32341e+03   2.06317e+01   6.95910e-02  -4.13516e-11
   // 2  Mean         1.99983e+02   9.48129e-02   3.91722e-04  -7.02884e-08
   // 3  Sigma        2.99563e+01   6.71117e-02   2.51639e-06  -8.63335e-05
   //-- h->GetNbinsX() = 256
   // h->Integral() = 100000
   // h->Integral("width") = 400000
   //-- TMath::Sqrt(2.*TMath::Pi()) * gampl*gsigma = 399731
   // fun->Integral(0,1024) = 399731
   // h->GetBinWidth(1) = 4
   //-- TMath::Sqrt(2.*TMath::Pi()) * gampl*gsigma / h->GetBinWidth(1) = 99932.7
   //-- MathMore integral = 399731
   // MathMore direct integral = 399731
}

#include <TArrow.h>
#include <TStyle.h>

void triagArea()
{
   TH1F* h3 = new TH1F("h3","N = 600*(40/10) = 2400. Exact = 2500", 12, 0, 120);
   h3->SetTitle("N_{tot} = 600*(40/10) = 2400 counts. Exact = 2500;E, arb. units, bin width #DeltaE = 10;N, counts");
   h3->SetFillColor(8);
   h3->SetBinContent(3, 100);
   h3->SetBinContent(4, 200);
   h3->SetBinContent(5, 400);
   h3->SetBinContent(6, 500);
   h3->SetBinContent(7, 600);
   h3->SetBinContent(8, 400);
   h3->SetBinContent(9, 200);
   h3->SetBinContent(10, 100);

   gStyle->SetOptStat(1001111);

   new TCanvas;
   h3->Draw();

   TArrow* arrow = new TArrow(40,300, 80,300, 0.05, "<>");
   arrow->SetFillStyle(3001);
   arrow->SetLineColor(2);
   arrow->SetLineWidth(3);
   arrow->SetAngle(45);
   arrow->Draw();
}

#if !defined(__CINT__)
/*
   --> to compile with MathMore
g++ `$ROOTSYS/bin/root-config --cflags --glibs` -lMathMore -o gauss gauss.C
*/
int main(int argc, char *argv[])
{
  // cout<< "argc = " << argc <<endl;  for (int i=0; i<argc; ++i) cout<< i <<" "<< argv[i] <<endl;
  // if (argc == 1) {
  //   cout<< "Usage: " << argv[0] << " par1 par2 par3" <<endl;
  //   return 0;
  // }

   TRint* theApp = new TRint("Rint", 0, 0, 0, 0, 1);  // do not show splash screen

   gauss();

   triagArea();

   theApp->Run();
   delete theApp;
   return 0;
}
#endif
