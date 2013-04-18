#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TLinearFitter.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TMarker.h>
#include <TText.h>

#include <iostream>
#include <cstdio>

using std::cout;	using std::endl;

// utils.C stuff
void nbi(Int_t nx=0, Int_t ny=0);
void fitgr(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
Double_t gaus_ampl();
Double_t gaus_mean();
Double_t gaus_sigma();

void peakfit_26_8V_noPZ()
{
   const char* ifname = "ketek50-100Hz-nolf-26.8V-1.dat.root-peaks.root";
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }
   TTree* peaks = (TTree*) gDirectory->Get("peaks");

   nbi(400);   // Set 400 bins by default

   new TCanvas;
   gPad->SetLogy();
   Int_t Nentries = peaks->Draw("pmax", "");

   Double_t w = peaks->GetHistogram()->GetXaxis()->GetXmax() / peaks->GetHistogram()->GetNbinsX();

   Double32_t gampl[100];
   Double32_t gmean[100];
   Double32_t gsigma[100];
   Int_t np = 0;

   fitgr(0.004, 0.006, "L");   
   gampl[np] = gaus_ampl();
   gmean[np] = gaus_mean();
   gsigma[np] = gaus_sigma();
   np++;

   fitgr(0.0087, 0.011, "L+");   
   // fitgr(0.0087, 0.011, "L");   
   gampl[np] = gaus_ampl();
   gmean[np] = gaus_mean();
   gsigma[np] = gaus_sigma();
   np++;

   fitgr(0.0138, 0.016, "L+", "sames");   
   // fitgr(0.0138, 0.016, "L", "sames");   
   gampl[np] = gaus_ampl();
   gmean[np] = gaus_mean();
   gsigma[np] = gaus_sigma();
   np++;

   fitgr(0.018, 0.0218, "L+", "sames");   
   // fitgr(0.018, 0.0218, "L", "sames");   
   gampl[np] = gaus_ampl();
   gmean[np] = gaus_mean();
   gsigma[np] = gaus_sigma();
   np++;

   fitgr(0.023, 0.027, "L+", "sames");   
   // fitgr(0.023, 0.027, "L", "sames");   
   gampl[np] = gaus_ampl();
   gmean[np] = gaus_mean();
   gsigma[np] = gaus_sigma();
   np++;

   fitgr(0.028, 0.032, "L+", "sames");   
   // fitgr(0.028, 0.032, "L", "sames");   
   gampl[np] = gaus_ampl();
   gmean[np] = gaus_mean();
   gsigma[np] = gaus_sigma();
   np++;

   // for nbi(400)
   TText* text = new TText(0.006, 500, "");
   text->DrawText(0.006, 500, "1.6 MHz");
   text->DrawText(0.011, 150, "400 kHz");
   text->DrawText(0.016, 40, "140 kHz");
   text->DrawText(0.022, 15, "60 kHz");
   text->DrawText(0.026, 6, "20 kHz");
   text->DrawText(0.032, 2.5, "5 kHz");

   TGraph* gr_dpmax = new TGraph(np-1);
   gr_dpmax->SetNameTitle("gr_dpmax", "KETEK, pixel 50 #mum, V_{bias} = 26.8 V;#Deltan;#Deltapmax, mV");
   gr_dpmax->SetMarkerStyle(20);
   gr_dpmax->SetMarkerColor(2);

   for (int i=1; i<np; ++i) {
      Double_t dpmax = gmean[i] - gmean[i-1];
      gr_dpmax->SetPoint(i-1, i, 1000.*dpmax);
      cout<< i << "\t dpmax = " << dpmax <<endl;
   }
   gr_dpmax->SetMinimum(4);
   
   new TCanvas;
   gr_dpmax->Draw("ap");

   for (int i=0; i<np; ++i) {
      Double32_t Npeaks = 2.5 * gampl[i] * gsigma[i] / w;
      Double32_t rate = Npeaks / (Nentries * 500e-3);
      cout<< i << "\t rate = " << rate << " MHz" <<endl;
   }

   Double32_t rate_tot = 0;
   for (int i=np-1; i>=0; --i) {
      Double32_t Npeaks = 2.5 * gampl[i] * gsigma[i] / w;
      Double32_t rate = Npeaks / (Nentries * 500e-3);
      rate_tot += rate;
      cout<< i << "\t rate_tot = " << rate_tot << " MHz" <<endl;
   }

   /*
      for nbi(400)
text = new TText(0.006, 500, "")
text->DrawText(0.006, 500, "1.6 MHz")
text->DrawText(0.011, 150, "400 kHz")
text->DrawText(0.016, 40, "140 kHz")
text->DrawText(0.022, 15, "60 kHz")
text->DrawText(0.026, 6, "20 kHz")
text->DrawText(0.032, 2.5, "5 kHz")

      for nbi(100)
text = new TText(0.006, 1000, "")
text->DrawText(0.006, 1000, "1.6 MHz")
text->DrawText(0.011, 300, "400 kHz")
text->DrawText(0.016, 80, "140 kHz")
text->DrawText(0.022, 30, "60 kHz")
text->DrawText(0.026, 12, "20 kHz")
text->DrawText(0.032, 5, "5 kHz")
      */
}

//void peakfit_26_8V_PZ()
void peakfit()
{
   const char* ifname = "ketek50-100Hz-lf_100pF_330ohm-26.8V-1.dat.root-peaks.root";
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }
   TTree* peaks = (TTree*) gDirectory->Get("peaks");

   nbi(400);   // Set 400 bins by default

   new TCanvas;
   gPad->SetLogy();
   Int_t Nentries = peaks->Draw("pmax", "");

   Double_t w = peaks->GetHistogram()->GetXaxis()->GetXmax() / peaks->GetHistogram()->GetNbinsX();

   Double32_t gampl[100];
   Double32_t gmean[100];
   Double32_t gsigma[100];
   Int_t np = 0;

   fitgr(0.0034, 0.0055, "L");   
   // fitgr(0.0034, 0.0055, "L");   
   gampl[np] = gaus_ampl();
   gmean[np] = gaus_mean();
   gsigma[np] = gaus_sigma();
   cout<< "np = " << np << "   gampl = " << gampl[np] << "   gmean = " << gmean[np] << "   gsigma = " << gsigma[np] <<endl;
   np++;

   fitgr(0.0072, 0.0098, "L+");   
   // fitgr(0.0072, 0.0098, "L");   
   gampl[np] = gaus_ampl();
   gmean[np] = gaus_mean();
   gsigma[np] = gaus_sigma();
   cout<< "np = " << np << "   gampl = " << gampl[np] << "   gmean = " << gmean[np] << "   gsigma = " << gsigma[np] <<endl;
   np++;

   fitgr(0.0115, 0.014, "L+", "sames");   
   // fitgr(0.0115, 0.014, "L", "sames");   
   gampl[np] = gaus_ampl();
   gmean[np] = gaus_mean();
   gsigma[np] = gaus_sigma();
   cout<< "np = " << np << "   gampl = " << gampl[np] << "   gmean = " << gmean[np] << "   gsigma = " << gsigma[np] <<endl;
   np++;

   fitgr(0.0156, 0.0184, "L+", "sames");   
   // fitgr(0.0156, 0.0184, "L", "sames");   
   gampl[np] = gaus_ampl();
   gmean[np] = gaus_mean();
   gsigma[np] = gaus_sigma();
   cout<< "np = " << np << "   gampl = " << gampl[np] << "   gmean = " << gmean[np] << "   gsigma = " << gsigma[np] <<endl;
   np++;

   fitgr(0.019, 0.022, "L+", "sames");   
   // fitgr(0.019, 0.022, "L", "sames");   
   gampl[np] = gaus_ampl();
   gmean[np] = gaus_mean();
   gsigma[np] = gaus_sigma();
   cout<< "np = " << np << "   gampl = " << gampl[np] << "   gmean = " << gmean[np] << "   gsigma = " << gsigma[np] <<endl;
   np++;

   fitgr(0.024, 0.027, "L+", "sames");   
   // fitgr(0.024, 0.027, "L", "sames");   
   gampl[np] = gaus_ampl();
   gmean[np] = gaus_mean();
   gsigma[np] = gaus_sigma();
   cout<< "np = " << np << "   gampl = " << gampl[np] << "   gmean = " << gmean[np] << "   gsigma = " << gsigma[np] <<endl;
   np++;

   for (int i=0; i<np; ++i) {
      cout<< i << "   gampl = " << gampl[i] << "   gmean = " << gmean[i] << "   gsigma = " << gsigma[i] <<endl;
   }

   // for nbi(400)
   TText* text = new TText(0.005, 500, "");
   text->DrawText(0.005, 800, "1.8 MHz");
   text->DrawText(0.009, 200, "440 kHz");
   text->DrawText(0.013, 70, "145 kHz");
   text->DrawText(0.018, 20, "53 kHz");
   text->DrawText(0.021, 8, "18 kHz");
   text->DrawText(0.026, 3.5, "7 kHz");

   TGraph* gr_dpmax = new TGraph(np-1);
   gr_dpmax->SetNameTitle("gr_dpmax", "KETEK, pixel 50 #mum, V_{bias} = 26.8 V;#Deltan;#Deltapmax, mV");
   gr_dpmax->SetMarkerStyle(20);
   gr_dpmax->SetMarkerColor(2);

   for (int i=1; i<np; ++i) {
      Double_t dpmax = gmean[i] - gmean[i-1];
      gr_dpmax->SetPoint(i-1, i, 1000.*dpmax);
      cout<< i << "\t dpmax = " << dpmax <<endl;
   }
   gr_dpmax->SetMinimum(4);
   
   new TCanvas;
   gr_dpmax->Draw("ap");

   for (int i=0; i<np; ++i) {
      Double32_t Npeaks = 2.5 * gampl[i] * gsigma[i] / w;
      Double32_t rate = Npeaks / (Nentries * 500e-3);
      cout<< i << "\t rate = " << rate << " MHz" <<endl;
   }

   Double32_t rate_tot = 0;
   for (int i=np-1; i>=0; --i) {
      Double32_t Npeaks = 2.5 * gampl[i] * gsigma[i] / w;
      Double32_t rate = Npeaks / (Nentries * 500e-3);
      rate_tot += rate;
      cout<< i << "\t rate_tot = " << rate_tot << " MHz" <<endl;
   }
}
