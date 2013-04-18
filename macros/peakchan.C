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
#include <TCut.h>
#include <TEventList.h>
#include <TObjArray.h>
#include <TBranch.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>

using std::cout;	using std::endl;

class FunRC {
private:
   static Double_t fRC0(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t sigma = par[npar++];

      Double_t x = *xx - x0;

      Double_t fRC0 = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t arg1 = -x/tau;
         Double_t exp1 = TMath::Abs(arg1) < 50? TMath::Exp(arg1): 0;

         Double_t arg2 = sigma*sigma/(2.*tau*tau);
         Double_t exp2 = TMath::Abs(arg2) < 50? TMath::Exp(arg2): 0;

         Double_t arg3 = sigma/tau/TMath::Sqrt(2.);
         if (TMath::Abs(sigma) > eps) arg3 -= x/sigma/TMath::Sqrt(2.);
         Double_t erfc = TMath::Erfc(arg3);
         fRC0 = A/(2.*tau) * exp1 * exp2 * erfc;
      }

      return fRC0;
   }
   static Double_t fRC0b(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t sigma = par[npar++];
      Double_t b = par[npar++];

      Double_t x = *xx - x0;

      Double_t fRC0b = b;
      if (TMath::Abs(tau) > eps) {
         Double_t arg1 = -x/tau;
         Double_t exp1 = TMath::Abs(arg1) < 50? TMath::Exp(arg1): 0;

         Double_t arg2 = sigma*sigma/(2.*tau*tau);
         Double_t exp2 = TMath::Abs(arg2) < 50? TMath::Exp(arg2): 0;

         Double_t arg3 = sigma/tau/TMath::Sqrt(2.);
         if (TMath::Abs(sigma) > eps) arg3 -= x/sigma/TMath::Sqrt(2.);
         Double_t erfc = TMath::Erfc(arg3);
         fRC0b = b + A/(2.*tau) * exp1 * exp2 * erfc;
      }

      return fRC0b;
   }
   static Double_t fRC1(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];

      Double_t x = *xx - x0;

      if (x <= 0) return 0;

      Double_t fRC1 = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t arg = -x/tau;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         Double_t integral = tau*tau;     // integral_0_inf x*exp(-x/tau)
         fRC1 = A * (x/integral) * exp;
      }

      return fRC1;
   }
   static Double_t fRC1b(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t b = par[npar++];

      Double_t x = *xx - x0;

      if (x <= 0) return b;

      Double_t fRC1b = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = tau*tau;     // integral_0_inf x*exp(-x/tau)
         fRC1b = b + (A/tau) * frac * exp;
      }

      return fRC1b;
   }
   static Double_t fRC2(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];

      Double_t x = *xx - x0;

      if (x <= 0) return 0;

      Double_t fRC2 = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = 2*tau*tau*tau;     // integral_0_inf x*x*exp(-x/tau)
         fRC2 = (A/tau) * frac*frac/2 * exp;
      }

      return fRC2;
   }
   static Double_t fRC2b(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t b = par[npar++];

      Double_t x = *xx - x0;

      if (x <= 0) return b;

      Double_t fRC2b = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = 2*tau*tau*tau;     // integral_0_inf x*x*exp(-x/tau)
         fRC2b = b + (A/tau) * frac*frac/2 * exp;
      }

      return fRC2b;
   }
   static Double_t fRC3(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];

      Double_t x = *xx - x0;

      if (x <= 0) return 0;

      Double_t fRC3 = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = (2*3)*tau*tau*tau*tau;   // integral_0_inf x*x*x*exp(-x/tau)
         fRC3 = (A/tau) * frac*frac*frac/(2*3) * exp;
      }

      return fRC3;
   }
   static Double_t fRC3b(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t b = par[npar++];

      Double_t x = *xx - x0;

      if (x <= 0) return b;

      Double_t fRC3b = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = (2*3)*tau*tau*tau*tau;   // integral_0_inf x*x*x*exp(-x/tau)
         fRC3b = b + (A/tau) * frac*frac*frac/(2*3) * exp;
      }

      return fRC3b;
   }
   static Double_t fRC4(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];

      Double_t x = *xx - x0;

      if (x <= 0) return 0;

      Double_t fRC4 = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = (2*3*4)*tau*tau*tau*tau*tau;   // integral_0_inf x*x*x*x*exp(-x/tau)
         fRC4 = (A/frac) * frac*frac*frac*frac/(2*3*4) * exp;
      }

      return fRC4;
   }
   static Double_t fRC4b(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t b = par[npar++];

      Double_t x = *xx - x0;

      if (x <= 0) return b;

      Double_t fRC4b = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = (2*3*4)*tau*tau*tau*tau*tau;   // integral_0_inf x*x*x*x*exp(-x/tau)
         fRC4b = b + (A/frac) * frac*frac*frac*frac/(2*3*4) * exp;
      }

      return fRC4b;
   }
public:
   static TF1* RC0(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau, Double_t sigma
         , const char* name="fRC0"
         )
   {
      TF1* fps0 = new TF1(name, fRC0, xmin, xmax, 4);
      fps0->SetNpx(1024);
      Int_t npar = 0;
      fps0->SetParName(npar++, "A");
      fps0->SetParName(npar++, "x0");
      fps0->SetParName(npar++, "tau");
      fps0->SetParName(npar++, "sigma");
      npar = 0;
      fps0->SetParameter(npar++, A);
      fps0->SetParameter(npar++, x0);
      fps0->SetParameter(npar++, tau);
      fps0->SetParameter(npar++, sigma);
      return fps0;
   }
   static TF1* RC0b(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau, Double_t sigma
         , Double_t b=0
         , const char* name="fRC0b"
         )
   {
      TF1* fps0b = new TF1(name, fRC0, xmin, xmax, 4);
      fps0b->SetNpx(1024);
      Int_t npar = 0;
      fps0b->SetParName(npar++, "A");
      fps0b->SetParName(npar++, "x0");
      fps0b->SetParName(npar++, "tau");
      fps0b->SetParName(npar++, "sigma");
      fps0b->SetParName(npar++, "b");
      npar = 0;
      fps0b->SetParameter(npar++, A);
      fps0b->SetParameter(npar++, x0);
      fps0b->SetParameter(npar++, tau);
      fps0b->SetParameter(npar++, sigma);
      fps0b->SetParameter(npar++, b);
      return fps0b;
   }
   static TF1* RC1(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , const char* name="fRC1"
         )
   {
      TF1* fps1 = new TF1(name, fRC1, xmin, xmax, 3);
      fps1->SetNpx(1024);
      Int_t npar = 0;
      fps1->SetParName(npar++, "A");
      fps1->SetParName(npar++, "x0");
      fps1->SetParName(npar++, "tau");
      npar = 0;
      fps1->SetParameter(npar++, A);
      fps1->SetParameter(npar++, x0);
      fps1->SetParameter(npar++, tau);
      return fps1;
   }
   static TF1* RC1b(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , Double_t b=0
         , const char* name="fRC1b"
         )
   {
      TF1* fps1b = new TF1(name, fRC1b, xmin, xmax, 4);
      fps1b->SetNpx(1024);
      Int_t npar = 0;
      fps1b->SetParName(npar++, "A");
      fps1b->SetParName(npar++, "x0");
      fps1b->SetParName(npar++, "tau");
      fps1b->SetParName(npar++, "b");
      npar = 0;
      fps1b->SetParameter(npar++, A);
      fps1b->SetParameter(npar++, x0);
      fps1b->SetParameter(npar++, tau);
      fps1b->SetParameter(npar++, b);
      return fps1b;
   }
   static TF1* RC2(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , const char* name="fRC2"
         )
   {
      TF1* fps2 = new TF1(name, fRC2, xmin, xmax, 3);
      fps2->SetNpx(1024);
      Int_t npar = 0;
      fps2->SetParName(npar++, "A");
      fps2->SetParName(npar++, "x0");
      fps2->SetParName(npar++, "tau");
      npar = 0;
      fps2->SetParameter(npar++, A);
      fps2->SetParameter(npar++, x0);
      fps2->SetParameter(npar++, tau);
      return fps2;
   }
   static TF1* RC2b(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , Double_t b=0
         , const char* name="fRC2b"
         )
   {
      TF1* fps2b = new TF1(name, fRC2b, xmin, xmax, 4);
      fps2b->SetNpx(1024);
      Int_t npar = 0;
      fps2b->SetParName(npar++, "A");
      fps2b->SetParName(npar++, "x0");
      fps2b->SetParName(npar++, "tau");
      fps2b->SetParName(npar++, "b");
      npar = 0;
      fps2b->SetParameter(npar++, A);
      fps2b->SetParameter(npar++, x0);
      fps2b->SetParameter(npar++, tau);
      fps2b->SetParameter(npar++, b);
      return fps2b;
   }
   static TF1* RC3(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , const char* name="fRC3"
         )
   {
      TF1* fps3 = new TF1(name, fRC3, xmin, xmax, 3);
      fps3->SetNpx(1024);
      Int_t npar = 0;
      fps3->SetParName(npar++, "A");
      fps3->SetParName(npar++, "x0");
      fps3->SetParName(npar++, "tau");
      npar = 0;
      fps3->SetParameter(npar++, A);
      fps3->SetParameter(npar++, x0);
      fps3->SetParameter(npar++, tau);
      return fps3;
   }
   static TF1* RC3b(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , Double_t b=0
         , const char* name="fRC3b"
         )
   {
      TF1* fps3b = new TF1(name, fRC3b, xmin, xmax, 4);
      fps3b->SetNpx(1024);
      Int_t npar = 0;
      fps3b->SetParName(npar++, "A");
      fps3b->SetParName(npar++, "x0");
      fps3b->SetParName(npar++, "tau");
      fps3b->SetParName(npar++, "b");
      npar = 0;
      fps3b->SetParameter(npar++, A);
      fps3b->SetParameter(npar++, x0);
      fps3b->SetParameter(npar++, tau);
      fps3b->SetParameter(npar++, b);
      return fps3b;
   }
   static TF1* RC4(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , const char* name="fRC4"
         )
   {
      TF1* fps4 = new TF1(name, fRC4, xmin, xmax, 3);
      fps4->SetNpx(1024);
      Int_t npar = 0;
      fps4->SetParName(npar++, "A");
      fps4->SetParName(npar++, "x0");
      fps4->SetParName(npar++, "tau");
      npar = 0;
      fps4->SetParameter(npar++, A);
      fps4->SetParameter(npar++, x0);
      fps4->SetParameter(npar++, tau);
      return fps4;
   }
   static TF1* RC4b(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , Double_t b=0
         , const char* name="fRC4b"
         )
   {
      TF1* fps4b = new TF1(name, fRC4b, xmin, xmax, 4);
      fps4b->SetNpx(1024);
      Int_t npar = 0;
      fps4b->SetParName(npar++, "A");
      fps4b->SetParName(npar++, "x0");
      fps4b->SetParName(npar++, "tau");
      fps4b->SetParName(npar++, "b");
      npar = 0;
      fps4b->SetParameter(npar++, A);
      fps4b->SetParameter(npar++, x0);
      fps4b->SetParameter(npar++, tau);
      fps4b->SetParameter(npar++, b);
      return fps4b;
   }
   //FunRC() {}
};

/*
   root -l ketek50-100Hz-lf_100pF_330ohm-26.0V-1.dat.root 'cpeaks.C+(2)'
*/

class Peak: public TObject {
public:
   Int_t chan;
   Int_t evt;
   Int_t npeak;
   Double_t x0;         // beginning of the pulse
   Double_t halfmax1;
   Double_t halfmax2;
   Double_t pmaxx;
   Double_t pmax;
   Double_t fwhm;
   Double_t charge;     // integral
   Double_t sbkg;
   Double_t baseline;
   Int_t np;
   Double_t* x;         //[np]
   Double_t* y;         //[np]
public:
   void Clear(Option_t*) {
      // cout<< "Peak::Clear(Option_t*)" <<endl;
      evt = -1;
      npeak = -1;
      x0 = 0;
      halfmax1 = 0;
      halfmax2 = 0;
      pmaxx = 0;
      pmax = 0;
      fwhm = 0;
      charge = 0;
      sbkg = 0;
      baseline = 0;
      np = 0;
      delete[] x;    x = 0;
      delete[] y;    y = 0;
   }
   Peak(): TObject(), chan(0), np(0), x(0), y(0) {}
   Peak(const Peak& peak): TObject(peak) {
      // cout<< "Peak copy ctor" <<endl;
      chan = peak.chan;
      evt = peak.evt;
      npeak = peak.npeak;
      x0 = peak.x0;
      halfmax1 = peak.halfmax1;
      halfmax2 = peak.halfmax2;
      pmaxx = peak.pmaxx;
      pmax = peak.pmax;
      fwhm = peak.fwhm;
      charge = peak.charge;
      sbkg = peak.sbkg;
      baseline = peak.baseline;
      np = peak.np;
      x = new Double_t[np];
      y = new Double_t[np];
      for (int i=0; i<np; ++i) {
         x[i] = peak.x[i];
         y[i] = peak.y[i];
      }
   }
   Peak& operator =(const Peak& peak) {
      // cout<< "Peak operator =" <<endl;
      chan = peak.chan;
      evt = peak.evt;
      npeak = peak.npeak;
      if (&peak == this) return *this;
      x0 = peak.x0;
      halfmax1 = peak.halfmax1;
      halfmax2 = peak.halfmax2;
      pmaxx = peak.pmaxx;
      pmax = peak.pmax;
      fwhm = peak.fwhm;
      charge = peak.charge;
      sbkg = peak.sbkg;
      baseline = peak.baseline;
      if (x) delete[] x;   x = 0;
      if (y) delete[] y;   y = 0;
      np = peak.np;
      x = new Double_t[np];
      y = new Double_t[np];
      for (int i=0; i<np; ++i) {
         x[i] = peak.x[i];
         y[i] = peak.y[i];
      }
      return *this;
   }
   ~Peak() {
      // cout<< "Peak dtor" <<endl;
      delete[] x;    x = 0;
      delete[] y;    y = 0;
   }
   void FillData(Int_t i1, Int_t i2, const Double_t* x_, const Double_t* y_) {
      //cout<< "Peak::FillData x=" << x << " y=" << y <<endl;
      if (x) {delete[] x; x = 0;}
      if (y) {delete[] y; y = 0;}
      np = i2 - i1 + 1;
      if (!x) x = new Double_t[np];
      if (!y) y = new Double_t[np];
      //cout<< "Peak::FillData will fill " << np << " points x_[" << i1 << "] = " << x_[i1] << " x_[" << i2 << "] = " << x_[i2] <<endl;
      for (int i=i1,ilocal=0; i<=i2; ++i,++ilocal) {
         x[ilocal] = x_[i];
         y[ilocal] = y_[i];
      }
      //cout<< "Peak::FillData filled " << np << " points" <<endl;
   }
   //-- TGraphErrors* Plot(Int_t nRC=3, Option_t* opt="ap") const
   TGraphErrors* Plot(Int_t nRC=0, Double_t sigma=-1., Option_t* opt="ap") const
   {
      TGraphErrors* gr = new TGraphErrors(np);
      gr->SetNameTitle(Form("peak_evt_%d_npeak_%d",evt,npeak), Form("evt %d chan %d peak %d I = %0.3f a.u.",evt,chan,npeak,charge));
      gr->SetMarkerStyle(20);
      gr->SetMarkerColor(8);
      gr->SetLineColor(8);

      // adjust baseline
      Double_t baseline_corr = 0;
      Double_t sbkg_corr = 0;
      Double_t sumy = 0;
      Double_t sumy2 = 0;
      Int_t npoints = 0;
      for (int i=0; i<np; ++i) {
         if (x[i] >= x0) break;
         npoints++;
         sumy += y[i];
         sumy2 += y[i]*y[i];
      }
      if (npoints >= 12) {
         baseline_corr = sumy / npoints;
         Double_t variance = npoints*sumy2 - sumy*sumy;
         sbkg_corr = (variance > 0)? TMath::Sqrt(variance/(npoints*(npoints-1))): 0;
      }

      for (int i=0; i<np; ++i) {
         gr->SetPoint(i, x[i], y[i]-baseline_corr);
         gr->SetPointError(i, 0, sbkg_corr);
      }
      gr->Draw(opt);

      Int_t i_x0 = 0;
      Int_t i_halfmax1 = 0;
      Int_t i_halfmax2 = 0;
      Int_t i_pmax = 0;
      Double_t delta = 0.1*(x[1] - x[0]);
      for (int i=1; i<np; ++i) {
         if (TMath::Abs(x[i] - x0) < delta) i_x0 = i;
         if (TMath::Abs(x[i] - halfmax1) < delta) i_halfmax1 = i;
         if (TMath::Abs(x[i] - halfmax2) < delta) i_halfmax2 = i;
         if (TMath::Abs(x[i] - pmaxx) < delta) i_pmax = i;
      }
      TMarker* marker_x0 = new TMarker(x[i_x0], y[i_x0], 20);
      marker_x0->SetMarkerColor(1);
      marker_x0->Draw();
      TMarker* marker_halfmax1 = new TMarker(x[i_halfmax1], y[i_halfmax1], 20);
      marker_halfmax1->SetMarkerColor(4);
      marker_halfmax1->Draw();
      TMarker* marker_halfmax2 = new TMarker(x[i_halfmax2], y[i_halfmax2], 20);
      marker_halfmax2->SetMarkerColor(9);
      marker_halfmax2->Draw();
      TMarker* marker_pmax = new TMarker(x[i_pmax], y[i_pmax], 20);
      marker_pmax->SetMarkerColor(2);
      marker_pmax->Draw();
      TMarker* marker_pmax_smooth = new TMarker(pmaxx, pmax, 24);
      marker_pmax_smooth->SetMarkerColor(2);
      marker_pmax_smooth->Draw();

      TF1* fRC = 0;
      switch (nRC) {
         case 0: fRC = FunRC::RC0b(gr->GetX()[0], gr->GetX()[gr->GetN()-1], charge, x0, fwhm/2.35, fwhm/5.70);
                 if (sigma >= 0) fRC->FixParameter(3, sigma);
                 break;
         case 1: fRC = FunRC::RC1b(gr->GetX()[0], gr->GetX()[gr->GetN()-1], charge, x0, fwhm/2.35);
                 break;
         case 2: fRC = FunRC::RC2b(gr->GetX()[0], gr->GetX()[gr->GetN()-1], charge, x0, fwhm/2.35);
                 break;
         case 3: fRC = FunRC::RC3b(gr->GetX()[0], gr->GetX()[gr->GetN()-1], charge, x0, fwhm/2.35);
                 break;
         case 4: fRC = FunRC::RC4b(gr->GetX()[0], gr->GetX()[gr->GetN()-1], charge, x0, fwhm/2.35);
                 break;
         default: cout<< "nRC is out of range" <<endl;
      }
      gr->Fit(fRC, "R", "", gr->GetX()[0],halfmax2);
      
      Double_t t = fRC->GetX(0.5*pmax, x0, pmaxx);
      printf("t = %0.3f ns at the halfmax %0.4f V\n", t, 0.5*pmax);

      gr->SetTitle(Form("%s t = %0.3f ns", gr->GetTitle(), t));
      gPad->Modified();
      gPad->Update();

      return gr;
   }

   ClassDef(Peak, 7);
};

//
// plot y:x with peaks->Draw("peaks[0].y:peaks[0].x","Entry$==2")
//
class PeakChannel: public TObject {
public:
   Int_t chan;
   Double_t fwhm_min;         // minimum peak width
   Double_t sbkg;          // sigma of background is used to find the peak boundary
   bool debug;             //!
   Int_t evt;
   Int_t np;               //!
   Double_t* x;            //!
   Double_t* y;            //!
   TClonesArray* peaks;    //->
   Double_t par[3];        // baseline fit
public:
   PeakChannel(): TObject()
        , chan(0)
        , fwhm_min(0), sbkg(0)
        , debug(false)
        , evt(-1)
        , np(0), x(0), y(0), peaks(0)
   {
      peaks = new TClonesArray("Peak");
      for (int ipar=0; ipar<3; ++ipar) par[ipar] = 0;
   }
   PeakChannel(Int_t chan_, Double_t fwhm_min_, bool debug_=false)
      :chan(chan_)
       , fwhm_min(fwhm_min_)
       , debug(debug_)
       , evt(-1)
       , np(0), x(0), y(0), peaks(0)
   {
      peaks = new TClonesArray("Peak");
      for (int ipar=0; ipar<3; ++ipar) par[ipar] = 0;
   }
   // PeakChannel(Int_t chan_, Double_t fwhm_min_, Float_t* x_float_, Float_t* y_float_, bool debug_=false)
   //    :chan(chan_)
   //     , fwhm_min(fwhm_min_)
   //     , debug(debug_)
   //     , evt(-1)
   //     , np(0), x(0), y(0)
   //     , x_float(x_float_), y_float(y_float_)
   //     , peaks(0)
   // {
   //    peaks = new TClonesArray("Peak");
   //    for (int ipar=0; ipar<3; ++ipar) par[ipar] = 0;
   // }
   PeakChannel(const PeakChannel& peakChannel): TObject(peakChannel)
   {
      //cout<< "PeakChannel copy ctor" <<endl;

      chan = peakChannel.chan;
      fwhm_min = peakChannel.fwhm_min;
      sbkg = peakChannel.sbkg;
      debug = peakChannel.debug;
      evt = peakChannel.evt;
      np = peakChannel.np;
      x = peakChannel.x;
      y = peakChannel.y;
      peaks = (TClonesArray*) peakChannel.peaks->Clone();  // http://root.cern.ch/root/roottalk/roottalk10/0488.html
      // peaks = new TClonesArray(*peakChannel.peaks);     // should work too
      //
      // peaks = new TClonesArray(peakChannel.GetClass());
      // for (int ipeak=0; ipeak<peakChannel.peaks->GetLast()+1; ++ipeak) {
      //    const Peak* peak = (const Peak*) peakChannel.peaks->At(ipeak); 
      //    new ((*peaks)[peaks->GetLast()+1]) Peak(*peak);
      // }
      for (int i=0; i<3; ++i) par[i] = peakChannel.par[i];
   }
   PeakChannel& operator =(const PeakChannel& peakChannel)
   {
      //cout<< "PeakChannel::operator =" <<endl;

      if (&peakChannel == this) return *this;

      chan = peakChannel.chan;
      fwhm_min = peakChannel.fwhm_min;
      sbkg = peakChannel.sbkg;
      debug = peakChannel.debug;
      evt = peakChannel.evt;
      np = peakChannel.np;
      x = peakChannel.x;
      y = peakChannel.y;
      peaks = (TClonesArray*) peakChannel.peaks->Clone();  // http://root.cern.ch/root/roottalk/roottalk10/0488.html
      for (int i=0; i<3; ++i) par[i] = peakChannel.par[i];
      return *this;
   }
   ~PeakChannel() {
      delete peaks;     peaks = 0;
   }

   void SetEventData(Int_t np_, Double_t* x_, Double_t* y_, Int_t evt_=-1) {
      evt = evt_;
      np = np_;
      x = x_;
      y = y_;
      peaks->Clear("C");
   }

   const Peak* GetPeak(Int_t npeak) const {
      return (const Peak*) peaks->At(npeak);
   }

   void PlotPeaks(Int_t iRC=0) const {
      for (int ipeak=0; ipeak<peaks->GetLast()+1; ++ipeak) {
         new TCanvas;
         const Peak* peak = (const Peak*) peaks->At(ipeak);
         peak->Plot(iRC);
      }
   }

   bool find_peak(Int_t i0, const Double_t* ys, Peak& peak, Int_t& iend)
   {
      // cout<< "PeakChannel::find_peak i0 = " << i0 <<endl;

      // peak.Clear("C");

      Int_t istart = i0;
      iend = np-1;

      // set threshold at 2 sigma
      // 5 consequtive points around the pulse maximum should exceed 3 sigma
      //-- const Double_t thres = 2.*sbkg;
      const Double_t thres = 3.*sbkg;
      //-- const Double_t thres_max = 2.*sbkg;

      Int_t ithres = np-1;
      for (int i=i0; i<np; ++i) {
         if (ys[i] < thres) continue;
         ithres = i;
         break;
      }
      if (ithres >= np - 3) return false;

      Int_t imaximum = ithres;
      Double_t threshold_low = sbkg;
      bool below_threshold_low = false;
      Int_t ibelow_threshold_low = ithres;

      // find the pulse end, when one of conditions below is satisfied:
      // --- ys < 0
      // or
      // --- ys fell below thrheshold (at some point) and start to rise

      for (int i=ithres+1; i<np; ++i)
      {
         if (ys[i] > ys[imaximum]) {
            imaximum = i;
            continue;
         }

         if (ys[i] < 0) {
            iend = i;
            break;
         }

         if (ys[i] < threshold_low) {
            below_threshold_low = true;
            if (ys[i] < ys[ibelow_threshold_low]) ibelow_threshold_low = i;
         }

         if (below_threshold_low) {
            if (ys[i] - ys[ibelow_threshold_low] > 0.5*sbkg) {   // beginning of the next peak
               iend = ibelow_threshold_low;
               break;
            }
         }

         //-- if (below_threshold_low and ys[i] > ys[i-1]) {
         //--    iend = i - 1;
         //--    break;
         //-- }
      }

      // find the pulse start. NB: at this point we know imaximum and threshold

      below_threshold_low = false;
      for (int i=ithres-1; i>=istart; --i) {
         if (ys[i] < threshold_low) below_threshold_low = true;

         if (ys[i] < 0) {
            istart = i;      // start from the negative point
            break;
         }

         if (below_threshold_low and ys[i] > ys[i+1]) {
            // istart = i + 1;
            istart = i;      // start from the high point
            break;
         }
      }

      // FWHM

      Int_t ihalfmax1 = istart;
      for (int i=istart; i<imaximum; ++i) {
         if (ys[i] > 0.5 * ys[imaximum]) {
            ihalfmax1 = i;
            break;
         }
      }

      Int_t ihalfmax2 = iend;
      for (int i=iend; i>imaximum; --i) {
         if (ys[i] > 0.5 * ys[imaximum]) {
            ihalfmax2 = i;
            break;
         }
      }

      // calculate sum between ihalfmax1 and ihalfmax2 (use actual values)
      Double_t halfmax_sum = 0;
      Int_t halfmax_n = 0;
      for (int i=ihalfmax1; i<=ihalfmax2; ++i) {
         halfmax_n++;
         halfmax_sum += y[i];
      }
      Double_t halfmax_sum_ratio = halfmax_sum / (sbkg*TMath::Sqrt(halfmax_n));

      peak.x0 = x[istart];
      peak.halfmax1 = x[ihalfmax1];
      peak.halfmax2 = x[ihalfmax2];
      peak.pmaxx = x[imaximum];
      peak.pmax = ys[imaximum];           // take ymax from ys, not from y
      peak.fwhm = x[ihalfmax2] - x[ihalfmax1];

      Int_t ibkg = istart - 2*(ihalfmax2 - ihalfmax1 + 1);  // add some background
      if (ibkg < i0) ibkg = i0;
      peak.FillData(ibkg, iend, x, y);
      
      // calculate integral using y, not ys
      peak.charge = 0;
      for (int i=istart+1; i<=iend; ++i) peak.charge += y[i]*(x[i] - x[i-1]);

      if (true
            && peak.fwhm > fwhm_min
            && halfmax_sum_ratio > 5
         ) {
         if (debug) {
            cout<< "::find_peak: x[istart] = " << x[istart] << " x[iend] = " << x[iend]
               << " fwhm = " << peak.fwhm << " charge = " << peak.charge
               << " halfmax_sum_ratio = " << halfmax_sum_ratio
            <<endl;
         }
         return true;
      }
      else {
         if (debug) {
            cout<< "::find_peak: --- rejected: x[istart] = "
               << x[istart] << " x[iend] = " << x[iend] << " fwhm = " << peak.fwhm << " charge = " << peak.charge
               << " halfmax_sum_ratio = " << halfmax_sum_ratio
            <<endl;
         }
         return false;
      }
   }

   void smooth(Double_t* ys) {
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

   void ProcessEvent(Int_t np_, Double_t* x_, Double_t* y_, Int_t evt_=-1)
   {
      //cout<< "PeakChannel::ProcessEvent" <<endl;

      SetEventData(np_, x_, y_, evt_);

      peaks->Clear("C");

      Peak peak;

      for (int ipar=0; ipar<3; ++ipar) par[ipar] = 0;

      // TLinearFitter* lfit = new TLinearFitter(1, "pol0");
      TLinearFitter lfit(1, "pol0");
      lfit.AssignData(np, 1, x, y);
      // Double_t x_double[1024], y_double[1024];
      // for (int i=0; i<1024; ++i) {
      //    x_double[i] = x[i];
      //    y_double[i] = y[i];
      // }
      // lfit->AssignData(np, 1, x_double, y_double);

      // lfit->EvalRobust(0.75);
      lfit.EvalRobust(0.80);

      // TVectorD params;
      // //-- TVectorD errors;
      // lfit.GetParameters(params);
      // //-- lfit->GetErrors(errors);

      // // TF1* fpol = new TF1("fpol","[0]+[1]*x+[2]*x**2", x[0], x[np-1]);
      // // fpol->SetLineColor(1);
      // // fpol->SetLineWidth(2);
      // // fpol->SetParameters(params[0], params[1], params[2]);
      // // fpol->Draw("same");
      // TF1* fpol = new TF1("fpol","[0]", x[0], x[np-1]);
      // fpol->SetLineColor(1);
      // fpol->SetLineWidth(2);
      // fpol->SetParameter(0, params[0]);

      // for pol0 only
      // Double_t baseline = params[0];
      // Double_t baseline = lfit.GetParameter(0);

      //-- for (int ipar=0; ipar<npar; ++ipar) par[ipar] = lfit.GetParameter(ipar);
      for (int ipar=0; ipar<lfit.GetNumberTotalParameters(); ++ipar) par[ipar] = lfit.GetParameter(ipar);

      // find the noise RMS

      // Double_t sumy = 0;
      // for (int i=0; i<np; ++i) sumy += y[i] - (par[0] + par[1]*x[i] + par[2]*x[i]*x[i]);
      // if (debug) cout<< "sumy / np = " << sumy / np <<endl;

      Double_t sumy_neg = 0;
      Double_t sumy2_tot = 0;
      Double_t sumy2_neg = 0;
      Double_t sumy2_pos = 0;
      Int_t np_neg = 0;
      Int_t np_pos = 0;
      for (int i=0; i<np; ++i) {
         //-- Double_t dy = y[i] - fpol->Eval(x[i]);
         Double_t dy = y[i] - (par[0] + par[1]*x[i] + par[2]*x[i]*x[i]);
         sumy2_tot += dy*dy;
         if (dy <= 0) {
            sumy_neg += dy;
            sumy2_neg += dy*dy;
            np_neg++;
         }
         else {
            sumy2_pos += dy*dy;
            np_pos++;
         }
      }
      // Double_t sbkg_tot = TMath::Sqrt(sumy2_tot/np);
      if (np_neg < 2) np_neg = 2;
      Double_t variance_neg = np_neg*sumy2_neg - sumy_neg*sumy_neg;
      Double_t sbkg_neg = (variance_neg > 0)? TMath::Sqrt(variance_neg/(np_neg*(np_neg-1))): 0;
      // cout<< "--- variance_neg = " << variance_neg <<endl;
      // Double_t sbkg_pos = TMath::Sqrt(sumy2_pos/np_pos);
      // if (debug) {
      //    cout<< "sbkg_tot = " << sbkg_tot <<endl;
      //    cout<< "sbkg_neg = " << sbkg_neg <<endl;
      //    cout<< "sbkg_pos = " << sbkg_pos <<endl;
      // }

      sbkg = sbkg_neg;

      Double_t ey[10240];     // large buffer

      for (int i=0; i<np; ++i) {
         y[i] = y[i] - (par[0] + par[1]*x[i] + par[2]*x[i]*x[i]);
         ey[i] = sbkg;
      }

      // smoothing
      Double_t ys[10240];     // large buffer
      smooth(ys);

      if (debug) {
         TGraphErrors* gr = new TGraphErrors(np, x, y, 0, ey);
         gr->SetNameTitle(Form("gr_%d",evt), Form("chan %d evt %d",chan,evt));
         gr->SetMarkerStyle(6);
         gr->SetMarkerColor(2);
         gr->SetLineColor(2);

         new TCanvas;
         gr->Draw("apl");

         // Double_t ymin = y[0];
         // Double_t ymax = y[0];
         // for (int i=1; i<np; ++i) {
         //    if (y[i] < ymin) ymin = y[i];
         //    if (y[i] > ymax) ymax = y[i];
         // }
         // ymin *= (ymin > 0)? 0.90: 1.10;
         // ymax *= (ymax > 0)? 1.10: 0.90;
         // TH1F* h = new TH1F("h", Form("Entry$==%d", evt), 100, ymin, ymax);
         // // TH1F* h = new TH1F("h", Form("Entry$==%d", evt), 100, -0.005, 0.012);   // good for event #2
         // for (int i=0; i<np; ++i) h->Fill(y[i]);
         // h->Fit("gaus");

         // new TCanvas;
         // h->Draw();

         TGraphErrors* gr_smooth = new TGraphErrors(np, x, ys, 0, ey);
         gr_smooth->SetNameTitle(Form("gr_smooth_%d",evt), Form("Smoothed %s",gr->GetTitle()));
         gr_smooth->SetMarkerStyle(7);
         gr_smooth->SetMarkerColor(4);
         gr_smooth->SetLineColor(4);

         new TCanvas;
         gr_smooth->Draw("apl");
      }

      //Double_t thr = 3*sbkg;
      //Double_t thr = 2*sbkg;
      if (debug) cout<< "sbkg = " << sbkg << " 2*sbkg = " << 2*sbkg << " 3*sbkg = " << 3*sbkg <<endl;

      for (int i=0; i<np; ++i) {
         Int_t iend;
         bool peak_found = find_peak(i, ys, peak, iend);
         // if (debug) cout<< "  i = " << i << " istart = " << istart << " x[istart] = " << x[istart] << "\t iend = " << iend << "\t x[iend] = " << x[iend] <<endl;
         if (peak_found) {
            Peak* p = new((*peaks)[peaks->GetLast()+1]) Peak(peak);
            // add to the peak the baseline data
            p->baseline = (par[0] + par[1]*x[i] + par[2]*x[i]*x[i]);
            p->sbkg = sbkg;
            p->chan = chan;
            p->evt = evt;
            p->npeak = peaks->GetLast();
         }
         i = iend;
         // add dead time
         Double_t xdead = x[i] + 2.*fwhm_min;
         while (x[i] < xdead) i++;
      }
      if (debug) cout<< "ProcessEvent: peaks->GetLast()+1 = " << peaks->GetLast()+1 <<endl;
   }

   ClassDef(PeakChannel, 6);
};

class PeakEvent: public TObject {
public:
   bool debug;                         //!
   //..PeakChannel* peakChannel[8];
   // TClonesArray* peakChannel;          //->
   TClonesArray* pc;          //->
   Double_t fwhm_min[8];
public:
   PeakEvent(): TObject() {
      pc = new TClonesArray("PeakChannel");
      for (int ichan=0; ichan<8; ++ichan) {
         //..peakChannel[ichan] = 0;
         fwhm_min[ichan] = 0;
      }
   }
   ~PeakEvent() {
      delete pc;
   }
   PeakEvent(const PeakEvent& peakEvent): TObject(peakEvent) {
      cout<< "PeakEvent copy ctor" <<endl;
      *this = peakEvent;
   }
   PeakEvent& operator =(const PeakEvent& peakEvent) {
      cout<< "PeakEvent::operator =" <<endl;
      if (&peakEvent == this) return *this;
      debug = peakEvent.debug;

      //pc->Clear("C");
      delete pc;
      pc = (TClonesArray*) peakEvent.pc->Clone();
      // for (int ichan=0; ichan<8; ++ichan) {
      //    //..peakChannel[ichan] = peakEvent.peakChannel[ichan];
      //    // const PeakChannel* peakChan= (const PeakChannel*) peakEvent.pc->At(ichan); 
      //    // new ((*pc)[pc->GetLast()+1]) PeakChannel(*peakChan);

      //    //peakChannel[ichan] = peakEvent.peakChannel[ichan];

      //    fwhm_min[ichan] = peakEvent.fwhm_min[ichan];
      // }
      return *this;
   }
   PeakEvent(const char* cfg, bool debug_=false): TObject(), debug(debug_)
   {
      if (!cfg || *cfg==0) cfg = "peakchan.cfg";

      //..for (int ichan=0; ichan<8; ++ichan) peakChannel[ichan] = 0;
      //---------- pc = new TClonesArray("Peak");

      pc = new TClonesArray("PeakChannel");

      // read the configuration file
      std::ifstream fcfg(cfg);
      if (!fcfg) {
         cout<< "PeakChannel: Could not open configuration file " << cfg <<endl;
         exit(0);
      }
      int chan;
      double fwhm = 0;
      double fwhm_default = 0;
      std::string fwhm_str;
      int iline = 0;
      std::string line;
      while (std::getline(fcfg, line)) {
         iline++;
         while (std::isspace(line[0])) line.erase(0,1);
         std::stringstream input;
         input << line;

         if (input.str().size() == 0) continue;
         if (input.str()[0] == '#') continue;
         if (input.str().size() > 1 && input.str()[0] == '/' && input.str()[1] == '/') continue;

         fwhm = fwhm_default = 0;

         input >> chan >> fwhm_str >> fwhm;

         if (true
               && fwhm_str == "fwhm"
            )
         {
            if (chan == 0) {
               cout<< "default chan fwhm = " << fwhm <<endl;
               fwhm_default = fwhm;
            }
            else {
               cout<< "chan = " << chan << " fwhm = " << fwhm <<endl;
               fwhm_min[chan-1] = fwhm;
            }
         }
         else {
            cout<<endl<< "PeakChannel: Input error in line " << iline << ": unknown keyword " << fwhm_str <<endl;
            exit(0);
         }
      }

      for (int ichan=0; ichan<8; ++ichan) {
         if (fwhm_min[ichan] == 0) {
            cout<< "Set to default channel " << ichan+1 <<endl;
            fwhm_min[ichan] = fwhm_default;
         }
      }

      for (int ichan=0; ichan<8; ++ichan) {
         cout<< "fwhm_min[" << ichan+1 << "] = " << fwhm_min[ichan] <<endl;
      }
   }

   ClassDef(PeakEvent, 3);
};

#ifdef __MAKECINT__
#pragma link C++ class Peak+;
#pragma link C++ class PeakChannel+;
#pragma link C++ class PeakEvent+;
#endif

ClassImp(Peak);
ClassImp(PeakChannel);
ClassImp(PeakEvent);

/*
peaks->Draw("peakEvent.pc[0].evt")
peaks->Draw("peakEvent.pc[0].@peaks.GetLast()+1")
*/

/*
root -l 'peakchan.C+("tb2012_run_24.root.pulse.root", "peakchan.cfg", 0,3, false)'
*/

void peakchan(const char* ifname="tb2012_run_24.root.pulse.root", const char* cfgname="peakchan.cfg"
      , Int_t ientry1=0
      , Int_t ientry2=-1
      , bool debug=false
      )
{
   PeakEvent* peakEvent = new PeakEvent(cfgname, debug);

   Float_t t1[1024];
   Float_t t2[1024];
   Float_t c[8][1024];

   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      exit(0);
   }

   TTree* itree = (TTree*) ifile->Get("p");
   TBranch* b;
   std::vector<int> slots;      // actual data channels (0..7)

   if ((b = itree->GetBranch("t1"))) {
      itree->SetBranchAddress("t1", &t1);
      for (int ichan=0; ichan<4; ++ichan) {
         if ((b = itree->GetBranch(Form("c%d",ichan+1)))) {
            slots.push_back(ichan);
            itree->SetBranchAddress(Form("c%d",ichan+1), c[ichan]);
            //..peakEvent->peakChannel[ichan] = new PeakChannel(ichan+1, peakEvent->fwhm_min[ichan], debug);
            new ((*peakEvent->pc)[slots[ichan]]) PeakChannel(ichan+1, peakEvent->fwhm_min[ichan], debug);
         }
      }
   }
   if ((b = itree->GetBranch("t2"))) {
      itree->SetBranchAddress("t2", &t2);
      for (int ichan=4; ichan<8; ++ichan) {
         if ((b = itree->GetBranch(Form("c%d",ichan+1)))) {
            slots.push_back(ichan);
            itree->SetBranchAddress(Form("c%d",ichan+1), c[ichan]);
            //..peakEvent->peakChannel[ichan] = new PeakChannel(ichan+1, peakEvent->fwhm_min[ichan], debug);
            new ((*peakEvent->pc)[slots[ichan]]) PeakChannel(ichan+1, peakEvent->fwhm_min[ichan], debug);
         }
      }
   }

   // for (unsigned islot=0; islot<slots.size(); ++islot) {
   //    new (*peakEvent->pc[peakEvent->pc->GetLast()+1]) PeakChannel;
   // }

   //..for (int i=0; i<8; ++i) peakEvent->peakChannel[i] = 0;
   // peakEvent->pc->Clear("C");

   TFile* ofile = TFile::Open(Form("%s.peaks.root",ifname), "recreate");
   TTree* otree = new TTree("peaks", "peaks tree");
   otree->Branch("peakEvent", "PeakEvent", &peakEvent);

   if (ientry2 < ientry1) ientry2 = itree->GetEntries()-1;
   Int_t nprocess = ientry2 - ientry1 + 1;

   cout<< nprocess << " entries to process" <<endl;

   for (int ientry=ientry1; ientry<=ientry2; ++ientry)
   {
      if (itree->LoadTree(ientry) < 0) break;
      itree->GetEntry(ientry);

      if (nprocess < 100) cout<< "--> processing entry " << ientry << "\n";
      else if (ientry % 100 == 0) cout<< "--> processing entry " << ientry << "\n";

      for (unsigned islot=0; islot<slots.size(); ++islot)
      {
         int ichan = slots[islot];
         //..if (peakEvent->peakChannel[ichan] != 0)      // NB: does not work: !peakEvent->peakChannel[ichan]
         if (debug) cout<< "--- channel " << ichan+1 << "   evt " << ientry <<endl;

         Float_t* x_float = (ichan < 4)? t1: t2;
         Float_t* y_float = c[ichan];
         Double_t x[1024];
         Double_t y[1024];

         for (int i=0; i<1024; ++i) {
            x[i] = x_float[i];
            y[i] = -1.*y_float[i];
         }
         //..peakEvent->peakChannel[ichan]->ProcessEvent(1024, x, y, ientry);
         PeakChannel* peakChannel = (PeakChannel*) peakEvent->pc->At(slots[ichan]);
         peakChannel->ProcessEvent(1024, x, y, ientry);
      }

      //cout<< "peakchan: otree->Fill()" <<endl;
      otree->Fill();
      //cout<< "peakchan: otree was filled successfully" <<endl;
   }

   cout<< "Writing " << otree->GetEntries() << " entries into file " << otree->GetCurrentFile()->GetName() <<endl;
   ofile->Write();

   // cout<< "last entry" <<endl;
   // for (int ichan=0; ichan<8; ++ichan) {
   //    TClonesArray* peaks = peakEvent->peakChannel[ichan]->peaks;
   //    for (int ipeak=0; ipeak<peaks->GetLast()+1; ++ipeak) {
   //       Peak* peak = (Peak*) peaks->At(ipeak);
   //       cout<< "ichan = " << ichan << " ipeak = " << ipeak << " peak->pmax = " << peak->pmax <<endl;
   //    }
   // }

   itree->ResetBranchAddresses();
}

void read_peakEvent(const char* ifname="tb2012_run_24.root.pulse.root.peaks.root", Int_t nRC=3, Int_t chan=-1, Int_t evt=-1)
{
   PeakEvent* peakEvent = 0;

   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }
   TTree* tree = (TTree*) ifile->Get("peaks");
   if (!tree) {
      cout<< "Could not found tree 'peaks'" <<endl;
      return;
   }

   tree->SetBranchAddress("peakEvent", &peakEvent);
   cout<< "tree->GetEntries() = " << tree->GetEntries() <<endl;

   // std::vector<Float_t> buffer;

   for (int ientry=0; ientry<tree->GetEntries(); ++ientry) {
      if (tree->LoadTree(ientry) < 0) break;
      tree->GetEntry(ientry);

      if (ientry < 5) {
         cout<< "plot peaks from the event " << ientry <<endl;
         for (int ichan=0; ichan<8; ++ichan) {
            //..TClonesArray* peaks = peakEvent->peakChannel[ichan]->peaks;
            PeakChannel* peakChannel = (PeakChannel*) peakEvent->pc->At(ichan);
            TClonesArray* peaks = peakChannel->peaks;
            for (int ipeak=0; ipeak<peaks->GetLast()+1; ++ipeak) {
               Peak* peak = (Peak*) peaks->At(ipeak);
               if (chan >= 0 && chan != peak->chan) continue;
               if (evt >= 0 && evt != peak->evt) continue;
               cout<< "ichan = " << ichan << " ipeak = " << ipeak << " peak->pmax = " << peak->pmax <<endl;
               new TCanvas;
               peak->Plot(nRC);

               // if (true
               //       && peak->pmax < 20
               //    )
               // {
               //    buffer.push_back(peak->pmax);
               // }
            }
         }
      }
   }

   // Double_t buffer_min = buffer[0];
   // Double_t buffer_max = buffer[0];
   // for (unsigned i=1; i<buffer.size(); ++i) {
   //    if (buffer[i] < buffer_min) buffer_min = buffer[i];
   //    if (buffer[i] > buffer_max) buffer_max = buffer[i];
   // }

   // Int_t nbins = 100;
   // if (buffer_min > 0) {
   //    if (0.01*(buffer_max - buffer_min) > buffer_min) buffer_min = 0;
   // }
   // else buffer_min *= 1.10;
   // buffer_max *= 1.10;           // assume positive

   // TH1F* h_pmax = new TH1F("h_pmax", "h_pmax", nbins, buffer_min, buffer_max);
   // for (unsigned i=0; i<buffer.size(); ++i) h_pmax->Fill(buffer[i]);

   // new TCanvas;
   // h_pmax->Draw();
}

TTree* peakfinder(TTree* tree, Int_t evt1=0, Int_t evt2=-1, Double_t fwhm=3., bool debug=false);
TTree* peakfinder(TTree* tree, Int_t evt1, Int_t evt2, Double_t fwhm, bool debug)
{
   Double_t x[1024];
   Double_t y[1024];

   if (evt2 < evt1) evt2 = tree->GetEntries()-1;
   Int_t nprocess = evt2 - evt1 + 1;
   cout<< nprocess << " entries to process" <<endl;

   if (debug && nprocess > 25) {
      cout<< "\n***WARNING peakfinder: Turn off debug mode\n" <<endl;
      debug = false;
   }

   PeakChannel* peakChannel = new PeakChannel(fwhm, debug);

   TTree* otree = new TTree("peaks", Form("peak tree for %s", tree->GetCurrentFile()->GetName()));
   otree->SetMarkerStyle(6);
   otree->SetMarkerColor(46);
   otree->SetLineColor(46);
   otree->Branch("peakChannel", "PeakChannel", &peakChannel);

   for (int ientry=evt1; ientry<=evt2; ++ientry) {
      if (tree->LoadTree(ientry) >= 0) tree->GetEntry(ientry);

      if (nprocess < 100) cout<< "--> processing entry " << ientry << "\n";
      else if (ientry % 100 == 0) cout<< "--> processing entry " << ientry << "\n";

      Float_t* x_float = (Float_t*) tree->GetLeaf("t1")->GetValuePointer();
      Float_t* y_float = (Float_t*) tree->GetLeaf("c1")->GetValuePointer();
      for (int i=0; i<1024; ++i) {
         x[i] = x_float[i];
         y[i] = -1.*y_float[i];
      }
      peakChannel->ProcessEvent(1024, x, y, ientry);
      otree->Fill();
   }

   if (debug) cout<< "otree->GetEntries() = " << otree->GetEntries() <<endl;

   return otree;
}

TTree* peakfinder(const char* ifname, Int_t evt1=0, Int_t evt2=-1, Double_t fwhm=3., bool debug=false);
TTree* peakfinder(const char* ifname, Int_t evt1, Int_t evt2, Double_t fwhm, bool debug)
{
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return 0;
   }

   TTree* tree = (TTree*) gDirectory->Get("p");
   if (!tree) {
      cout<< "Could not find tree t" <<endl;
      return 0;
   }

   TFile* ofile = TFile::Open(Form("%s-peaks.root",ifname), "recreate");

   TTree* otree = peakfinder(tree, evt1,evt2, fwhm, debug);

   cout<< "Write " << otree->GetEntries() << " entries into output file " << ofile->GetName() <<endl;
   ofile->Write();

   return otree;
}

void readpeak(Int_t evt, Int_t npeak, TTree* tree=0);
void readpeak(Int_t evt, Int_t npeak, TTree* tree)
{
   if (!tree) tree = (TTree*) gDirectory->Get("peaks");
   if (!tree) {
      cout<< "Could not find tree \"peaks\"" <<endl;
      return;
   }

   PeakChannel* peakChannel = 0;
   tree->SetBranchAddress("peakChannel", &peakChannel);

   if (tree->LoadTree(evt) < 0) {
      cout<< "Could not find entry " << evt << " in the tree \"peaks\"" <<endl;
      return;
   }
   tree->GetEntry(evt);

   Peak* peak = (Peak*) peakChannel->peaks->At(npeak);
   if (!peak) {
      cout<< "npeak " << npeak << " is out of range: peaks->GetLast()+1 = " << peakChannel->peaks->GetLast()+1 <<endl;
      return;
   }

   TGraphErrors* gr = new TGraphErrors(peak->np);
   gr->SetNameTitle(Form("gr_%d_%d", evt,npeak), Form("Entry %d peak %d", evt,npeak));
   gr->SetMarkerStyle(7);
   gr->SetMarkerColor(46);
   gr->SetLineColor(46);
   for (int i=0; i<peak->np; ++i) {
      gr->SetPoint(i, peak->x[i], peak->y[i]);
      gr->SetPointError(i, 0, peak->sbkg);
   }

   new TCanvas;
   gr->Draw("ap");
}

/*
speaks->Draw("y:x","Entry$==0","lp")
*/
void select_peaks(const char* ifname, Int_t chan, Double_t pmax1, Double_t pmax2, Int_t entry1=0, Int_t entry2=-1);
void select_peaks(const char* ifname, Int_t chan, Double_t pmax1, Double_t pmax2, Int_t entry1, Int_t entry2)
{
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }

   TTree* itree = (TTree*) gDirectory->Get("peaks");
   if (!itree) {
      cout<< "Could not find input tree \"peaks\"" <<endl;
      return;
   }

   PeakEvent* peakEvent = 0;
   itree->SetBranchAddress("peakEvent", &peakEvent);

   TFile* ofile = TFile::Open(Form("%s-sel_%0.4f-%0.4f.root", ifname,pmax1,pmax2), "recreate");

   Peak* opeak = new Peak;
   TTree* otree = new TTree("speaks", Form("selected Peak samples"));
   otree->Branch("speak", "Peak", &opeak);
   otree->SetMarkerStyle(7);
   otree->SetMarkerColor(8);
   otree->SetLineColor(8);

   if (entry2 < entry1) entry2 = itree->GetEntries();
   Int_t nentries = entry2 - entry1 + 1;
   cout<< nentries << " entries to process" <<endl;

   for (int ientry=entry1; ientry<=entry2; ++ientry)
   {
      if (itree->LoadTree(ientry) < 0) break;
      itree->GetEntry(ientry);

      if (ientry % 100 == 0) cout<< "processing entry " << ientry <<endl;

      PeakChannel* peakChannel = (PeakChannel*) peakEvent->pc->At(chan-1);
      if (!peakChannel) {
         cout<< "Could not find channel " << chan << " in the data file " << itree->GetCurrentFile()->GetName() <<endl;
         return;
      }

      for (int ipeak=0; ipeak<peakChannel->peaks->GetLast()+1; ++ipeak) {
         Peak* peak = (Peak*) peakChannel->peaks->At(ipeak);
         if (true
               and peak->pmax > pmax1
               and peak->pmax < pmax2
            )
         {
            *opeak = *peak;
            otree->Fill();
         }
      }
   }

   new TCanvas;
   otree->Draw("pmax", "");

   cout<< "\nWriting " << otree->GetEntries() << " events into file " << ofile->GetName() <<endl;
   ofile->Write();
}

void select_peaks_gate(const char* ifname, Double_t trig, Double_t gate, Int_t entry1=0, Int_t entry2=-1);
void select_peaks_gate(const char* ifname, Double_t trig, Double_t gate, Int_t entry1, Int_t entry2)
{
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }

   TTree* itree = (TTree*) gDirectory->Get("peaks");
   if (!itree) {
      cout<< "Could not find input tree \"peaks\"" <<endl;
      return;
   }

   PeakChannel* peakChannel = 0;
   itree->SetBranchAddress("peakChannel", &peakChannel);

   TFile* ofile = TFile::Open(Form("%s-sel_%0.0fns.root", ifname,trig), "recreate");

   Peak* opeak = new Peak;
   TTree* otree = new TTree("speaks", Form("selected Peak samples"));
   otree->Branch("speak", "Peak", &opeak);
   otree->SetMarkerStyle(7);
   otree->SetMarkerColor(8);
   otree->SetLineColor(8);

   if (entry2 < entry1) entry2 = itree->GetEntries();
   Int_t nentries = entry2 - entry1 + 1;
   cout<< nentries << " entries to process" <<endl;

   for (int ientry=entry1; ientry<=entry2; ++ientry)
   {
      if (itree->LoadTree(ientry) < 0) break;
      itree->GetEntry(ientry);

      if (ientry % 100 == 0) cout<< "processing entry " << ientry <<endl;

      for (int ipeak=0; ipeak<peakChannel->peaks->GetLast()+1; ++ipeak) {
         Peak* peak = (Peak*) peakChannel->peaks->At(ipeak);
         if (true
               and peak->x0 > trig
               and peak->x0 < trig+gate
            )
         {
            *opeak = *peak;
            otree->Fill();
         }
      }
   }

   new TCanvas;
   otree->Draw("pmax", "");

   cout<< "\nWriting " << otree->GetEntries() << " events into file " << ofile->GetName() <<endl;
   ofile->Write();
}

void select_peaks_custom(const char* ifname, const char* suffix="", Int_t entry1=0, Int_t entry2=-1);
void select_peaks_custom(const char* ifname, const char* suffix, Int_t entry1, Int_t entry2)
{
   Double_t trig = 50;     // ns
   Double_t gate = 10;     // ns
   Double_t pmax_min = 2;  // mV

   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }

   TTree* itree = (TTree*) gDirectory->Get("peaks");
   if (!itree) {
      cout<< "Could not find input tree \"peaks\"" <<endl;
      return;
   }

   PeakChannel* peakChannel = 0;
   itree->SetBranchAddress("peakChannel", &peakChannel);

   TFile* ofile = TFile::Open(Form("%s-sel_%s.root", ifname,suffix), "recreate");

   Peak* opeak = new Peak;
   TTree* otree = new TTree("speaks", Form("selected Peak samples"));
   otree->Branch("speak", "Peak", &opeak);
   otree->SetMarkerStyle(7);
   otree->SetMarkerColor(8);
   otree->SetLineColor(8);

   if (entry2 < entry1) entry2 = itree->GetEntries();
   Int_t nentries = entry2 - entry1 + 1;
   cout<< nentries << " entries to process" <<endl;

   for (int ientry=entry1; ientry<=entry2; ++ientry)
   {
      if (itree->LoadTree(ientry) < 0) break;
      itree->GetEntry(ientry);

      if (ientry % 100 == 0) cout<< "processing entry " << ientry <<endl;

      for (int ipeak=0; ipeak<peakChannel->peaks->GetLast()+1; ++ipeak) {
         Peak* peak = (Peak*) peakChannel->peaks->At(ipeak);
         if (true
               and peak->x0 > trig
               and peak->x0 < trig+gate
               and peak->pmax > pmax_min
            )
         {
            *opeak = *peak;
            otree->Fill();
         }
      }
   }

   new TCanvas;
   otree->Draw("pmax", "");

   cout<< "\nWriting " << otree->GetEntries() << " events into file " << ofile->GetName() <<endl;
   ofile->Write();
}

PeakChannel* getPeakChannel(Int_t chan, Int_t evt, TTree* peaks=0);
PeakChannel* getPeakChannel(Int_t chan, Int_t evt, TTree* peaks)
{
   if (!peaks) peaks = (TTree*) gDirectory->Get("peaks");
   if (!peaks) {
      cout<< "Could not find tree 'peaks'" <<endl;
      return 0;
   }

   PeakEvent* peakEvent = 0;
   peaks->SetBranchAddress("peakEvent", &peakEvent);

   if (peaks->LoadTree(evt) >= 0) {
      peaks->GetEntry(evt);

      for (int ichan=0; ichan<peakEvent->pc->GetLast()+1; ++ichan) {
         PeakChannel* peakChannel = (PeakChannel*) peakEvent->pc->At(ichan);
         //cout<< "peakChannel->chan = " << peakChannel->chan <<endl;
         if (peakChannel->chan == chan) {
            cout<< "peakChannel->peaks->GetLast()+1 = " << peakChannel->peaks->GetLast()+1 <<endl;
            // peakChannel->PlotPeaks(iRC);
            return peakChannel;
         }
      }
   }
   
   return 0;
}

PeakEvent* getPeakEvent(Int_t evt, TTree* peaks=0);
PeakEvent* getPeakEvent(Int_t evt, TTree* peaks)
{
   if (!peaks) peaks = (TTree*) gDirectory->Get("peaks");
   if (!peaks) {
      cout<< "Could not find tree 'peaks'" <<endl;
      return 0;
   }

   PeakEvent* peakEvent = 0;
   peaks->SetBranchAddress("peakEvent", &peakEvent);

   if (peaks->LoadTree(evt) >= 0) {
      peaks->GetEntry(evt);
      return peakEvent;
   }
   
   return 0;
}

void plotpeaks(Int_t chan, Int_t evt, Int_t iRC=0, TTree* peaks=0);
void plotpeaks(Int_t chan, Int_t evt, Int_t iRC, TTree* peaks)
{
   if (!peaks) peaks = (TTree*) gDirectory->Get("peaks");
   if (!peaks) {
      cout<< "Could not find tree 'peaks'" <<endl;
      return;
   }

   PeakEvent* peakEvent = 0;
   peaks->SetBranchAddress("peakEvent", &peakEvent);

   if (peaks->LoadTree(evt) >= 0) {
      peaks->GetEntry(evt);

      for (int ichan=0; ichan<peakEvent->pc->GetLast()+1; ++ichan) {
         PeakChannel* peakChannel = (PeakChannel*) peakEvent->pc->At(ichan);
         // cout<< "peakChannel->chan = " << peakChannel->chan <<endl;
         if (peakChannel->chan == chan) {
            cout<< "peakChannel->peaks->GetLast()+1 = " << peakChannel->peaks->GetLast()+1 <<endl;
            peakChannel->PlotPeaks(iRC);
         }
      }
   }
}

Peak* getpeak(Int_t evt, TTree* speaks=0);
Peak* getpeak(Int_t evt, TTree* speaks)
{
   if (!speaks) speaks = (TTree*) gDirectory->Get("speaks");
   if (!speaks) {
      cout<< "Could not find tree 'speaks'" <<endl;
      return 0;
   }

   Peak* peak = 0;
   speaks->SetBranchAddress("speak", &peak);

   if (speaks->LoadTree(evt) >= 0) speaks->GetEntry(evt);

   //new TCanvas;
   //peak->plot();
   return new Peak(*peak);
}

TGraphErrors* peakplot(Int_t evt, Int_t nRC=3, TTree* speaks=0);
TGraphErrors* peakplot(Int_t evt, Int_t nRC, TTree* speaks)
{
   if (!speaks) speaks = (TTree*) gDirectory->Get("speaks");
   if (!speaks) {
      cout<< "Could not find tree 'speaks'" <<endl;
      return 0;
   }

   Peak* peak = 0;
   speaks->SetBranchAddress("speak", &peak);

   if (speaks->LoadTree(evt) >= 0) speaks->GetEntry(evt);

   new TCanvas;
   return peak->Plot(nRC);
}
