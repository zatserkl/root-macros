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

#include <iostream>
#include <cstdio>

using std::cout;	using std::endl;

class FunRC {
private:
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
   Peak(): TObject(), np(0), x(0), y(0) {}
   Peak(const Peak& peak): TObject(peak) {
      // cout<< "Peak copy ctor" <<endl;
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
   TGraphErrors* Plot(Int_t nRC=3, Option_t* opt="ap") const {
      TGraphErrors* gr = new TGraphErrors(np);
      gr->SetNameTitle(Form("peak_evt_%d_npeak_%d",evt,npeak), Form("evt %d peak %d I = %0.3f a.u.",evt,npeak,charge));
      gr->SetMarkerStyle(20);
      gr->SetMarkerColor(8);
      gr->SetLineColor(8);
      for (int i=0; i<np; ++i) {
         gr->SetPoint(i, x[i], y[i]);
         gr->SetPointError(i, 0, sbkg);
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
class PeakEvent: public TObject {
public:
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
   PeakEvent(): TObject()
        , fwhm_min(0), sbkg(0)
        , debug(false)
        , evt(-1)
        , np(0), x(0), y(0), peaks(0)
   {
      peaks = new TClonesArray("Peak");
      for (int ipar=0; ipar<3; ++ipar) par[ipar] = 0;
   }
   PeakEvent(Double_t fwhm_min_, bool debug_=false)
      : fwhm_min(fwhm_min_)
        , debug(debug_)
        , evt(-1)
        , np(0), x(0), y(0), peaks(0)
   {
      peaks = new TClonesArray("Peak");
      for (int ipar=0; ipar<3; ++ipar) par[ipar] = 0;
   }
   PeakEvent(const PeakEvent& peakEvent): TObject(peakEvent) {
      fwhm_min = peakEvent.fwhm_min;
      sbkg = peakEvent.sbkg;
      debug = peakEvent.debug;
      evt = peakEvent.evt;
      np = peakEvent.np;
      x = peakEvent.x;
      y = peakEvent.y;
      peaks = (TClonesArray*) peakEvent.peaks->Clone();  // http://root.cern.ch/root/roottalk/roottalk10/0488.html
      // peaks = new TClonesArray(*peakEvent.peaks);     // should work too
      //
      // peaks = new TClonesArray(peakEvent.GetClass());
      // for (int ipeak=0; ipeak<peakEvent.peaks->GetLast()+1; ++ipeak) {
      //    const Peak* peak = (const Peak*) peakEvent.peaks->At(ipeak); 
      //    new ((*peaks)[peaks->GetLast()+1]) Peak(*peak);
      // }
      for (int i=0; i<3; ++i) par[i] = peakEvent.par[i];
   }
   ~PeakEvent() {
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

   void PlotPeaks() const {
      for (int ipeak=0; ipeak<peaks->GetLast()+1; ++ipeak) {
         new TCanvas;
         const Peak* peak = (const Peak*) peaks->At(ipeak);
         peak->Plot();
      }
   }

   bool find_peak(Int_t i0, const Double_t* ys, Peak& peak, Int_t& iend)
   {
      // cout<< "PeakEvent::find_peak i0 = " << i0 <<endl;

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
            sumy2_neg += dy*dy;
            np_neg++;
         }
         else {
            sumy2_pos += dy*dy;
            np_pos++;
         }
      }
      // Double_t sbkg_tot = TMath::Sqrt(sumy2_tot/np);
      Double_t sbkg_neg = TMath::Sqrt(sumy2_neg/np_neg);
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
         gr->SetNameTitle(Form("gr_%d",evt), Form("Entry$==%d",evt));
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
         gr_smooth->SetNameTitle(Form("gr_smooth_%d",evt), Form("Smoothed Entry$==%d",evt));
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

   ClassDef(PeakEvent, 6);
};

#ifdef __MAKECINT__
#pragma link C++ class Peak+;
#pragma link C++ class PeakEvent+;
#endif

/*
root -l ketek50-100Hz-lf_100pF_330ohm-26.0V-1.dat.root
peakfinder(p, 2,2, true)

plot y:x with peaks->Draw("peaks[0].y:peaks[0].x","Entry$==2")

--> works
root -l ketek50-100Hz-lf_100pF_330ohm-26.0V-1.dat.root 'peakfinder.C+(p, 20,29, true)'
--> but not: root -l ketek50-100Hz-lf_100pF_330ohm-26.0V-1.dat.root peakfinder.C+\(p, 20,29, true\)
*/

ClassImp(Peak);
ClassImp(PeakEvent);

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

   PeakEvent* peakEvent = new PeakEvent(fwhm, debug);

   TTree* otree = new TTree("peaks", Form("peak tree for %s", tree->GetCurrentFile()->GetName()));
   otree->SetMarkerStyle(6);
   otree->SetMarkerColor(46);
   otree->SetLineColor(46);
   otree->Branch("peakEvent", "PeakEvent", &peakEvent);

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
      peakEvent->ProcessEvent(1024, x, y, ientry);
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

   PeakEvent* peakEvent = 0;
   tree->SetBranchAddress("peakEvent", &peakEvent);

   if (tree->LoadTree(evt) < 0) {
      cout<< "Could not find entry " << evt << " in the tree \"peaks\"" <<endl;
      return;
   }
   tree->GetEntry(evt);

   Peak* peak = (Peak*) peakEvent->peaks->At(npeak);
   if (!peak) {
      cout<< "npeak " << npeak << " is out of range: peaks->GetLast()+1 = " << peakEvent->peaks->GetLast()+1 <<endl;
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
void select_peaks(const char* ifname, Double_t pmax1, Double_t pmax2, Int_t entry1=0, Int_t entry2=-1);
void select_peaks(const char* ifname, Double_t pmax1, Double_t pmax2, Int_t entry1, Int_t entry2)
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

      for (int ipeak=0; ipeak<peakEvent->peaks->GetLast()+1; ++ipeak) {
         Peak* peak = (Peak*) peakEvent->peaks->At(ipeak);
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

   PeakEvent* peakEvent = 0;
   itree->SetBranchAddress("peakEvent", &peakEvent);

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

      for (int ipeak=0; ipeak<peakEvent->peaks->GetLast()+1; ++ipeak) {
         Peak* peak = (Peak*) peakEvent->peaks->At(ipeak);
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

   PeakEvent* peakEvent = 0;
   itree->SetBranchAddress("peakEvent", &peakEvent);

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

      for (int ipeak=0; ipeak<peakEvent->peaks->GetLast()+1; ++ipeak) {
         Peak* peak = (Peak*) peakEvent->peaks->At(ipeak);
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
      return new PeakEvent(*peakEvent);
   }
   
   return 0;
}

void plotpeaks(Int_t evt, TTree* peaks=0);
void plotpeaks(Int_t evt, TTree* peaks)
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
      peakEvent->PlotPeaks();
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
