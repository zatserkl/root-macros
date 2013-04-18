#include "drs.C"

#include <TROOT.h>
#include <TEnv.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TCanvas.h>
//-- #include <Math/GSLIntegrator.h>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2.h>

#include <iostream>
#include <sstream>
#include <cmath>

using std::cout;     using std::endl;

// utils.C stuff
void pic(TString pathname="", TString suffix="");
TH1* htemp(const char* name=0, const char* title=0);
TF1* ftemp(const char* name=0, const char* title=0);
TGraph* gtemp(const char* name=0, const char* title=0);
void fitgaus(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void fitpol(Double_t xmin=0, Double_t xmax=0, Int_t power=1, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void rightgaus();
void leftgaus();
void left();
void right();
Int_t nbi(Int_t n=0);
const char* addtit(const char* title, TCanvas* can=0);

bool debug = false;
bool gdebug = false;

class OscFit: public TObject
{
public:
   Int_t evt;
   Float_t adc[8];
   Float_t adcf[8];        // ADC from fit
   Float_t t[8];
   Float_t d[8];
   Float_t tau1[8];
   Float_t tau2[8];
   Float_t T[8];
   Float_t sigma[8];
   Float_t bkg[8];         // flat background before the signal
   Float_t sbkg[8];        // sigma of the flat background
   Float_t v[8];
   Float_t chi2[8];
public:
   void clear() {
      evt = 0;
      for (int i=0; i<8; ++i) {
         adc[i] = 0;
         adcf[i] = 0;
         t[i] = 0;
         d[i] = 0;
         tau1[i] = 0;
         tau2[i] = 0;
         T[i] = 0;
         sigma[i] = 0;
         bkg[i] = 0;
         sbkg[i] = 0;
         v[i] = 0;
         chi2[i] = 0;
      }
   }
   OscFit(): TObject() {clear();}
   ClassDef(OscFit,8);
};

#ifdef __MAKECINT__
#pragma link C++ class OscFit;
#endif

class MathLim {
public:
   static Double_t pow_(Double_t x, Int_t n) {
      Double_t pow_ = std::pow(x,n);               // from <cmath>
      return pow_;
   }
   static Double_t sqrt_(Double_t x) {
      Double_t sqrt_ = 0;
      if (x > 0) sqrt_ = TMath::Sqrt(x);
      return sqrt_;
   }
   static Double_t exp_(Double_t x) {
      Double_t exp_ = 0;
      if (x < 500.) exp_ = TMath::Exp(x);
      return exp_;
   }
   static Double_t erfc_(Double_t x) {
      Double_t erfc_ = 0;
      if (TMath::Abs(x) < 100.) erfc_ = TMath::Erfc(x);
      else erfc_ = x < 0? 2.: 0;
      return erfc_;
   }
};

//bool debug_ = 0;

class PulseFunction: public MathLim {
public:
   enum {ipar_A=0, ipar_x0, ipar_tau1, ipar_tau2, ipar_T, ipar_sigma};
public:
   Double_t par_[10];
   Double_t A_, x0_, tau1_, tau2_, T_, sigma_;
   Double_t xmin_, xmax_;
   Int_t color_;
   Int_t npar_;
   TF1* fpulse_;
public:
   void Set_A(Double_t A) {par_[ipar_A] = A_ = A;}
   void Set_x0(Double_t x0) {par_[ipar_x0] = x0_ = x0;}
   void Set_tau1(Double_t tau1) {par_[ipar_tau1] = tau1_ = tau1;}
   void Set_tau2(Double_t tau2) {par_[ipar_tau2] = tau2_ = tau2;}
   void Set_T(Double_t T) {par_[ipar_T] = T_ = T;}
   void Set_sigma(Double_t sigma) {par_[ipar_sigma] = sigma_ = sigma;}
public:
   PulseFunction(Double_t A, Double_t x0, Double_t tau1, Double_t tau2, Double_t T, Double_t sigma,
         Double_t xmin, Double_t xmax, Int_t color=2)
      : A_(A), x0_(x0), tau1_(tau1), tau2_(tau2), T_(T), sigma_(sigma)
      , xmin_(xmin), xmax_(xmax)
      , color_(color)
      , npar_(6)
      , fpulse_(0)
   {
      createFunction();
   }
   // PulseFunction(const PulseFunction& pf) {
   //    A_ = pf.A_;
   //    x0_ = pf.x0_;
   //    tau1_ = pf.tau1_;
   //    tau2_ = pf.tau2_;
   //    T_ = pf.T_;
   //    sigma_ = pf.sigma_;
   //    xmin_ = pf.xmin_;
   //    xmax = pf.xmax;
   //    npar_ = pf.npar_;
   //    for (int ipar=0; i<npar_; ++ipar_) par_[ipar] = pf.par_[i];
   // }
  ~PulseFunction() {
     delete fpulse_;    fpulse_ = 0;
   }
   TF1* fun() {return fpulse_;}
   void createFunction() {
      if (fpulse_) delete fpulse_;
      npar_ = 6;
      fpulse_ = new TF1("fpulse", fPsigma, xmin_, xmax_, npar_);
      fpulse_->SetLineColor(color_);
      fpulse_->SetNpx(1024);
      par_[ipar_A] = A_;
      par_[ipar_x0] = x0_;
      par_[ipar_tau1] = tau1_;
      par_[ipar_tau2] = tau2_;
      par_[ipar_T] = T_;
      par_[ipar_sigma] = sigma_;
      fpulse_->SetParameters(par_);

      std::stringstream sstit;
      sstit << fpulse_->GetName() << " #tau_{1}=" << tau1_ << " #tau_{2}=" << tau2_ << " T=" << T_ << " #sigma=" << sigma_ <<endl;
      fpulse_->SetTitle(sstit.str().c_str());

      fpulse_->SetParName(ipar_A, "A");
      fpulse_->SetParName(ipar_x0, "x0");
      fpulse_->SetParName(ipar_tau1, "#tau_{1}");
      fpulse_->SetParName(ipar_tau2, "#tau_{2}");
      fpulse_->SetParName(ipar_T, "T");
      fpulse_->SetParName(ipar_sigma, "#sigma");
   }
   static Double_t ITtausigma(Double_t x, Double_t tau, Double_t T, Double_t sigma)
   {
      using TMath::Pi;
      using TMath::Exp;
      using TMath::Sqrt;
      using TMath::Abs;
      using TMath::Erfc;

      const Double_t eps = 1e-12;

      Double_t ITtausigma = 0;

      Double_t arg;

      if (T < 0 || tau < 0) return 0;     // sanity check

      if (tau < eps) return 0;            // tau is a factor in front of ITtausigma

      if (Abs(sigma) < eps && x < eps) return 0;   // range of P(x)

      // tau --> T
      if (Abs(T-tau) < eps) {
         Double_t tt = T > tau? T: tau;
         if (tt < eps) return 0;          // both tau and T are zero

         if (Abs(sigma) < eps)
         {
            // no smearing: return ITtau(x,tau,T);
            /// ITtausigma = (x/tt)*Exp(-x/tt);
            ITtausigma = (x/tt)*exp_(-x/tt);
            return ITtausigma;
         }
         else {
            /// // the first term
            /// Double_t exp_x = Exp(-x/tt);
            /// Double_t exp_sigma = Exp(sigma*sigma/2/tt/tt);
            /// arg = (x - sigma*sigma/tau);
            /// Double_t exp_x2 = Exp(-arg*arg);
            /// Double_t term1 = 1./(2.*Pi()) * exp_x * exp_sigma * (sigma/tt) * exp_x2;
            /// // the second term
            /// Double_t term2 = (x - sigma*sigma/tau)/tau * exp_x;
            /// // result for tau --> T
            /// ITtausigma = term1 + term2;
            /// return ITtausigma;
            // the first term
            Double_t exp_x = exp_(-x/tt);
            Double_t exp_sigma = exp_(sigma*sigma/2/tt/tt);
            arg = (x - sigma*sigma/tau);
            Double_t exp_x2 = exp_(-arg*arg);
            Double_t term1 = 1./(2.*Pi()) * exp_x * exp_sigma * (sigma/tt) * exp_x2;
            // the second term
            Double_t term2 = (x - sigma*sigma/tau)/tau * exp_x;
            // result for tau --> T
            ITtausigma = term1 + term2;
            if (ITtausigma < eps) ITtausigma = 0;
            return ITtausigma;
         }
      }

      if (T < eps)
      {
         if (sigma < eps)
         {
            // no smearing
            ITtausigma = exp_(-x/tau);
         }
         else {
            /// ITtausigma = 0.5*Exp(-x/tau)*Exp(sigma*sigma/2./tau/tau)*Erfc(-(x-sigma*sigma/tau)/sigma/Sqrt(2.));
            ITtausigma = 0.5*exp_(-x/tau)*exp_(sigma*sigma/2./tau/tau)*erfc_(-(x-sigma*sigma/tau)/sigma/sqrt_(2.));
         }
         return ITtausigma;
      }

      // regular case

      if (sigma < eps)
      {
         // no smearing
         ITtausigma = tau * (exp_(-x/T) - exp_(-x/tau)) / (T - tau);
      }
      else {
         /// Double_t exp_x_T = 0;
         /// arg = -x/T;
         /// if (arg < 500.) exp_x_T = Exp(arg);
         /// Double_t exp_sigma_T = 0;
         /// arg = sigma*sigma/2./T/T;
         /// if (arg < 500.) exp_sigma_T = Exp(arg);
         /// Double_t erfc_T = 0;
         /// arg = -(x-sigma*sigma/T)/sigma/Sqrt(2.);
         /// if (Abs(arg) < 10.) erfc_T = Erfc(arg);
         /// else erfc_T = (arg < 0)? 2.: 0;
         /// if (debug_) cout<< "exp_x_T " << exp_x_T << " exp_sigma_T " << exp_sigma_T << " erfc_T " << erfc_T <<endl;
         /// // Double_t ITsigma = 0.5 * Exp(-x/T) * Exp(sigma*sigma/2./T/T) * Erfc(-(x-sigma*sigma/T)/sigma/Sqrt(2.));
         /// Double_t ITsigma = 0.5 * exp_x_T * exp_sigma_T * erfc_T;
         Double_t ITsigma = 0.5 * exp_(-x/T) * exp_(sigma*sigma/2./T/T) * erfc_(-(x-sigma*sigma/T)/sigma/sqrt_(2.));

         /// Double_t exp_x_tau = 0;
         /// arg = -x/tau;
         /// if (arg < 500.) exp_x_tau = Exp(arg);
         /// Double_t exp_sigma_tau = 0;
         /// arg = sigma*sigma/2./tau/tau;
         /// if (arg < 500.) exp_sigma_tau = Exp(arg);
         /// Double_t erfc_tau = 0;
         /// arg = -(x-sigma*sigma/tau)/sigma/Sqrt(2.);
         /// if (Abs(arg) < 10.) erfc_tau = Erfc(arg);
         /// else erfc_tau = (arg < 0)? 2.: 0;
         /// if (debug_) cout<< "exp_x_tau " << exp_x_tau << " exp_sigma_tau " << exp_sigma_tau << " erfc_tau " << erfc_tau <<endl;
         /// // Double_t Itausigma = 0.5 * Exp(-x/tau) * Exp(sigma*sigma/2./tau/tau) * Erfc(-(x-sigma*sigma/tau)/sigma/Sqrt(2.));
         /// Double_t Itausigma = 0.5 * exp_x_tau * exp_sigma_tau * erfc_tau;
         Double_t Itausigma = 0.5 * exp_(-x/tau) * exp_(sigma*sigma/2./tau/tau) * erfc_(-(x-sigma*sigma/tau)/sigma/sqrt_(2.));

         ITtausigma = tau/(T-tau) * (ITsigma - Itausigma);
         //if (debug_) cout<< "ITsigma = " << ITsigma << " Itausigma = " << Itausigma << " ITtausigma = " << ITtausigma <<endl;
         if (ITtausigma < 0) ITtausigma = 0;
      }

      return ITtausigma;
   }

   static Double_t Psigma(Double_t x, Double_t tau1, Double_t tau2, Double_t T, Double_t sigma)
   {
      const Double_t eps = 1.e-12;

      if (tau1<0 || tau2<0 || T<0) return 0;       // sanity check

      if (tau2 < eps) return 0;                    // by definition tau2 should be finite

      if (TMath::Abs(sigma) < eps && x < eps) return 0;  // range of P(x)

      Double_t tau12 = tau1*tau2/(tau1+tau2);
      Double_t ITtausigma2 = ITtausigma(x,tau2,T,sigma);
      Double_t ITtausigma12 = ITtausigma(x,tau12,T,sigma);
      //-- Double_t ITtausigma12 = 0;
      Double_t norm = (tau1+tau2)/tau2/tau2;
      Double_t Psigma = norm * (ITtausigma2 - ITtausigma12);
      if (Psigma < eps) Psigma = 0;
      return Psigma;
   }

   static Double_t fPsigma(Double_t* xx, Double_t* par)
   {
      const Double_t eps = 1.e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau1 = par[npar++];
      Double_t tau2 = par[npar++];
      Double_t T = par[npar++];
      Double_t sigma = par[npar++];

      Double_t x = *xx - x0;

      if (tau1 < 0 || tau2 < 0 || T < 0) return 0;    // sanity check
      if (tau2 < eps) return 0;                       // discharge time physical limit

      if (TMath::Abs(sigma) < eps && x < eps) return 0;  // range of P(x)

      return A*Psigma(x,tau1,tau2,T,sigma);
   }
}; //------------------------- end of class PulseFunction -----------------------------

void pulseclass(
        Double_t A = 100.
      , Double_t x0 = 80.
      , Double_t tau1 = 0.100
      , Double_t tau2 = 100.
      , Double_t T = 40.
      , Double_t sigma = 0.100
      , Double_t xmin=0.
      , Double_t xmax=512.
      )
{
   const char ifname[] = "Co60_STM_LSO2x2_NOCAP_NOAMP_split.xml.root";
   // const char ifname[] = "Co60_STM_LSO2x2_nocap_noamp_50ns.xml.root";
   Long64_t entry = 1;
   Int_t ich = 1;             // counting from 0
   Double_t bkgmin = 0;
   Double_t bkgmax = 80;
   Double_t sigmin = 80;
   // Double_t sigmax = 512;
   Double_t sigmax = 360;

   // PulseFunction* pfexample = new PulseFunction(A,x0,tau1,tau2,T,sigma, xmin,xmax);
   // pfexample->fpulse_->Draw();
   PulseFunction pfexample(A,x0,tau1,tau2,T,sigma, xmin,xmax, 2);
   pfexample.fpulse_->DrawCopy();                               // NB: DrawCopy instead of just Draw

   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return;
   }
   TTree* tree = (TTree*) ifile->Get("t");
   if (!tree) {
      cout<< "No tree \"t\" has been found in " << ifname <<endl;
      return;
   }
   cout<< "tree->GetEntries() = " << tree->GetEntries() <<endl;

   OscEvent* oscEvent = 0;
   tree->SetBranchAddress("oscEvent",&oscEvent);

   // histogram which will be reused
   TH1F* hbkg = new TH1F("hbkg","hbkg", 800,-400,400);

   // get entry and obtain pointer to channel ich
   tree->GetEntry(entry);
   const OscChannel* oscChannel = (OscChannel*) oscEvent->oscChannels->At(ich);

   // find flat background

   Double_t x[1024], y[1024], ex[1024], ey[1024];
   Int_t np = 0;

   np = 0;
   for (int i=0; i<1024; ++i) {
      if (oscChannel->x[i] < bkgmin) continue;
      if (oscChannel->x[i] > bkgmax) break;
      x[np] = oscChannel->x[i];
      y[np] = oscChannel->y[i];
      ++np;
   }

   // sort bkg array to eliminate possible USB spikes: three highest points
   Int_t index[1024];
   TMath::Sort(np, y, index, kFALSE);

   hbkg->Reset();
   hbkg->SetNameTitle(Form("hbkg_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch),Form("hbkg_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch));
   for (int i=0; i<np-3; ++i) hbkg->Fill(y[index[i]]);
   new TCanvas;
   // hbkg->Fit("gaus", "L0Q", "goff");                            // NB: Loglikelihood option
   hbkg->Fit("gaus", "L", "");                            // NB: Loglikelihood option

   Double_t bkg_mean = hbkg->GetFunction("gaus")->GetParameter("Mean");
   Double_t bkg_sigma = hbkg->GetFunction("gaus")->GetParameter("Sigma");

   np = 0;
   Double_t ymax = 0;
   for (int i=0; i<1024; ++i) {
      if (oscChannel->x[i] < sigmin) continue;
      if (oscChannel->x[i] > sigmax) break;
      x[np] = oscChannel->x[i];
      y[np] = -1.*(oscChannel->y[i] - bkg_mean);
      if (y[np] > ymax) ymax = y[np];
      ex[np] = 0;
      ey[np] = bkg_sigma;
      ++np;
   }

   Double_t thres = 20;
   if (ymax > thres) cout<< "ymax = " << ymax << " exceeds threshold" <<endl;

   TGraphErrors* gr = new TGraphErrors(np,x,y,ex,ey);
   gr->SetNameTitle(Form("gr_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch), Form("gr_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch));
   gr->SetMarkerStyle(24);
   gr->SetMarkerColor(8);
   gr->SetLineColor(8);

   new TCanvas;
   gr->Draw("ap");

   Double_t integral = 0;
   for (int i=0; i<np-1; ++i) integral += 0.5*(gr->GetY()[i] + gr->GetY()[i+1]) * (gr->GetX()[i+1]-gr->GetX()[i]);

   // normalize to pC by division to R = 50 Hom
   integral /= 50.;
   // cout<< ".. integral = " << integral << " pC" <<endl;

   // -------------------------------------------------------------
   //    Fit
   // -------------------------------------------------------------

   PulseFunction pf(A,x0,tau1,tau2,T,sigma, sigmin,sigmax, 2);

   //pf.fun()->FixParameter(pf.fun()->GetParNumber("tau1"), pf.fun->GetParameter("tau1"));
   //pf.fun()->FixParameter(pf.fun()->GetParNumber("sigma"), pf.fun->GetParameter("sigma"));
   //pf.fun()->FixParameter(pf.fun()->GetParNumber("T"), pf.fun->GetParameter("T"));
   pf.fun()->FixParameter(PulseFunction::ipar_tau1, 0.050);
   pf.fun()->FixParameter(PulseFunction::ipar_sigma, 0.050);
   pf.fun()->FixParameter(PulseFunction::ipar_T, 40.);
   gr->Fit(pf.fpulse_, "R", "", sigmin,sigmax);
   cout<< ".. integral = " << integral << " pC" << " pf.fpulse_->Integral(0,1000)/50. = " << pf.fpulse_->Integral(0,1000)/50. <<endl;
}

// class DRS4Fitter {
// public:
//    TGraphErrors* gr_;
//    PulseFunction* pf_fit_;
//    PulseFunction* pf_refit_;
//    Double_t sigmin_, sigmax_;
//    DRS4Fitter(TGraphErrors* gr, Double_t sigmin, Double_t sigmax)
//       : gr_(gr)
//       , sigmin_(sigmin)
//       , sigmax_(sigmax)
//    {
//    }
//    void SetParameters(Double_t A, Double_t x0, Double_t tau1, Double_t tau2, Double_t T, Double_t sigma) {
//    }
//    void Fit() {
//    }
// };

Double_t ITtausigma(Double_t x, Double_t tau, Double_t T, Double_t sigma)
{
   using TMath::Pi;
   using TMath::Exp;
   using TMath::Sqrt;
   using TMath::Abs;
   using TMath::Erfc;

   const Double_t eps = 1e-12;
   Double_t ITtausigma = 0;

   if (T < 0 || tau < 0) return 0;     // sanity check

   if (tau < eps) return 0;            // tau is a factor in front of ITtausigma

   if (Abs(sigma) < eps && x < eps) return 0;   // range of P(x)

   // tau --> T
   if (Abs(T-tau) < eps) {
      Double_t tt = T > tau? T: tau;
      if (tt < eps) return 0;          // both tau and T are zero

      if (Abs(sigma) < eps)
      {
         // no smearing: return ITtau(x,tau,T);
         ITtausigma = (x/tt)*Exp(-x/tt);
         return ITtausigma;
      }
      else {
         // the first term
         Double_t exp_x = Exp(-x/tt);
         Double_t exp_sigma = Exp(sigma*sigma/2/tt/tt);
         Double_t arg = (x - sigma*sigma/tau);
         Double_t exp_x2 = Exp(-arg*arg);
         Double_t term1 = 1./(2.*Pi()) * exp_x * exp_sigma * (sigma/tt) * exp_x2;
         // the second term
         Double_t term2 = (x - sigma*sigma/tau)/tau * exp_x;
         // result for tau --> T
         ITtausigma = term1 + term2;
         return ITtausigma;
      }
   }

   if (T < eps)
   {
      if (sigma < eps)
      {
         // no smearing
         ITtausigma = Exp(-x/tau);
      }
      else {
         ITtausigma = 0.5*Exp(-x/tau)*Exp(sigma*sigma/2./tau/tau)*Erfc(-(x-sigma*sigma/tau)/sigma/Sqrt(2.));
      }
      return ITtausigma;
   }

   // regular case

   if (sigma < eps)
   {
      // no smearing
      ITtausigma = tau * (Exp(-x/T) - Exp(-x/tau)) / (T - tau);
   }
   else {
      Double_t ITsigma = 0.5 * Exp(-x/T) * Exp(sigma*sigma/2./T/T) * Erfc(-(x-sigma*sigma/T)/sigma/Sqrt(2.));
      Double_t Itausigma = 0.5 * Exp(-x/tau) * Exp(sigma*sigma/2./tau/tau) * Erfc(-(x-sigma*sigma/tau)/sigma/Sqrt(2.));
      ITtausigma = tau/(T-tau) * (ITsigma - Itausigma);
   }

   return ITtausigma;
}

Double_t Psigma(Double_t x, Double_t tau1, Double_t tau2, Double_t T, Double_t sigma)
{
   const Double_t eps = 1.e-12;

   if (tau1<0 || tau2<0 || T<0) return 0;       // sanity check

   if (tau2 < eps) return 0;                    // by definition tau2 should be finite

   if (TMath::Abs(sigma) < eps && x < eps) return 0;  // range of P(x)

   Double_t tau12 = tau1*tau2/(tau1+tau2);
   Double_t ITtausigma2 = ITtausigma(x,tau2,T,sigma);
   Double_t ITtausigma12 = ITtausigma(x,tau12,T,sigma);
   Double_t norm = (tau1+tau2)/tau2/tau2;
   Double_t Psigma = norm * (ITtausigma2 - ITtausigma12);
   return Psigma;
}

Double_t fPsigma(Double_t* xx, Double_t* par)
{
   const Double_t eps = 1.e-12;

   Int_t npar = 0;
   Double_t A = par[npar++];
   Double_t x0 = par[npar++];
   Double_t tau1 = par[npar++];
   Double_t tau2 = par[npar++];
   Double_t T = par[npar++];
   Double_t sigma = par[npar++];

   Double_t x = *xx - x0;

   if (tau1 < 0 || tau2 < 0 || T < 0) return 0;    // sanity check
   if (tau2 < eps) return 0;                       // discharge time physical limit

   if (TMath::Abs(sigma) < eps && x < eps) return 0;  // range of P(x)

   return A*Psigma(x,tau1,tau2,T,sigma);
}

TTree* pulse(TTree *tree
      , Double_t bkgmin=0, Double_t bkgmax=40, Double_t sigmin=40, Double_t sigmax=120
      , Int_t entry_first=0, Int_t entry_last=-1
      , bool setdebug=false, bool setgdebug=false
      , Int_t channel=-1
      )
{
   if (setdebug) debug = true;
   if (setgdebug) gdebug = true;

   //-- Double_t thres = 20.;
   Double_t thres = 10.;
   Int_t nthres_min = 5;                  // data are required to contain at least nthres_min points over the thres value
   Double_t ysaturation = 499;

   tree->SetMarkerStyle(7);
   tree->SetMarkerColor(2);

   OscEvent* oscEvent = 0;
   tree->SetBranchAddress("oscEvent",&oscEvent);

   cout<< "tree->GetEntries() = " << tree->GetEntries() <<endl;

   // output (fit results) tree
   TTree* otree = new TTree("ft", "Fit result tree");
   OscFit* oscFit = new OscFit;
   otree->Branch("oscFit", "OscFit", &oscFit);

   // the number of channels in the data
   tree->GetEntry(entry_first);
   Int_t nchannels = oscEvent->oscChannels->GetEntries();
   
   // histogram for background fit
   TH1F *hbkg = new TH1F("hbkg", "hbkg", 400, -100, 100);

   Double_t A = 1000;
   Double_t x0 = 50;
   Double_t tau1 = 1;
   Double_t tau2 = 5;
   Double_t T = 10;
   Double_t sigma = 1;

   Double_t par[10];
   Int_t npar = 0;

   enum {ipar_A=0, ipar_x0, ipar_tau1, ipar_tau2, ipar_T, ipar_sigma};

   par[npar++] = A;
   par[npar++] = x0;
   par[npar++] = tau1;
   par[npar++] = tau2;
   par[npar++] = T;
   par[npar++] = sigma;
   cout<< "\t npar = " << npar <<endl;

   // function for fit
   TF1* fun_Psigma = new TF1("fun_Psigma", fPsigma, sigmin, sigmax, npar);
   fun_Psigma->SetNpx(1024);
   fun_Psigma->SetName("fun_Psigma");
   std::stringstream sstit;
   sstit << fun_Psigma->GetName() << " #tau_{1}=" << tau1 << " #tau_{2}=" << tau2 << " T=" << T << " #sigma=" << sigma;
   fun_Psigma->SetTitle(sstit.str().c_str());
   fun_Psigma->SetLineColor(1);

   fun_Psigma->SetParameters(par);
   fun_Psigma->SetParName(ipar_A, "A");
   fun_Psigma->SetParName(ipar_x0, "x0");
   fun_Psigma->SetParName(ipar_tau1, "tau1");
   fun_Psigma->SetParName(ipar_tau2, "tau2");
   fun_Psigma->SetParName(ipar_T, "T");
   fun_Psigma->SetParName(ipar_sigma, "sigma");

   fun_Psigma->FixParameter(ipar_tau1,1);
   fun_Psigma->FixParameter(ipar_sigma,1);

   // fun_Psigma->SetParLimits(ipar_tau2, 1., 20.);
   // fun_Psigma->SetParLimits(ipar_T, 1., 60.);

   // function for refit
   TF1* fun_Psigma_refit = new TF1("fun_Psigma_refit", fPsigma, sigmin, sigmax, npar);
   fun_Psigma_refit->SetNpx(1024);
   fun_Psigma_refit->SetName("fun_Psigma_refit");
   std::stringstream sstit_refit;
   sstit_refit << fun_Psigma_refit->GetName() << " #tau_{1}=" << tau1 << " #tau_{2}=" << tau2 << " T=" << T << " #sigma=" << sigma;
   fun_Psigma_refit->SetTitle(sstit_refit.str().c_str());
   fun_Psigma_refit->SetLineColor(2);

   fun_Psigma_refit->SetParameters(par);
   fun_Psigma_refit->SetParName(ipar_A, "A");
   fun_Psigma_refit->SetParName(ipar_x0, "x0");
   fun_Psigma_refit->SetParName(ipar_tau1, "tau1");
   fun_Psigma_refit->SetParName(ipar_tau2, "tau2");
   fun_Psigma_refit->SetParName(ipar_T, "T");
   fun_Psigma_refit->SetParName(ipar_sigma, "sigma");

   if (entry_last < 0) entry_last = tree->GetEntries() - 1;

   for (int jentry=entry_first; jentry<=entry_last; ++jentry)
   {
      tree->GetEntry(jentry);
      cout<< "\n---> jentry = " << jentry << " evt = " << oscEvent->evt <<endl;

      oscFit->clear();
      oscFit->evt = oscEvent->evt;

      cout<< ".. loop over osc channels. oscEvent->oscChannels = " << oscEvent->oscChannels->GetEntries() <<endl;
      for (int ich=0; ich<nchannels; ++ich)
      {
         const OscChannel* oscChannel = (OscChannel*) oscEvent->oscChannels->At(ich);

         if (oscChannel->ch < 3 || oscChannel->ch >4) continue;

         if (channel > -1 && oscChannel->ch != channel) continue;
         //-- cout<< "oscChannel->ch = " << oscChannel->ch <<endl;
         //-- if (oscChannel->ch == 1) {
         //--    cout<< "------------------> skip oscilloscope channel " << oscChannel->ch <<endl;
         //--    continue;
         //-- }

         //
         // find the level of flat background
         //

         Double_t x[1024], y[1024], ex[1024], ey[1024];
         Int_t np = 0;

         Double_t ybkg[1024];
         Int_t np_bkg = 0;
         for (int i=0; i<1024; ++i) {
            if (oscChannel->x[i] < bkgmin) continue;
            if (oscChannel->x[i] > bkgmax) break;
            ybkg[np_bkg] = -1.*oscChannel->y[i];
            ++np_bkg;
         }
         // sort bkg array to eliminate possible USB spikes: three highest points
         Int_t index[1024];
         TMath::Sort(np_bkg, ybkg, index, kFALSE);

         hbkg->Reset();
         for (int i=0; i<np_bkg-3; ++i) hbkg->Fill(ybkg[index[i]]);
         hbkg->Fit("gaus", "L0Q", "goff");                            // NB: Loglikelihood option

         Double_t bkg_mean = hbkg->GetFunction("gaus")->GetParameter("Mean");
         Double_t bkg_sigma = hbkg->GetFunction("gaus")->GetParameter("Sigma");
         // if (gdebug) {
         //    hbkg->SetNameTitle(Form("hbkg_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch),Form("hbkg_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch));
         //    new TCanvas;
         //    hbkg->DrawCopy();
         // }

         Int_t nthres = 0;                // the number of channels which exceed threshold
         Float_t ampl_max = 0;
         np = 0;
         for (int i=0; i<1024; ++i)
         {
            if (oscChannel->x[i] < sigmin) continue;
            if (oscChannel->x[i] > sigmax) break;

            // maximum in signal window NB: (-1)*y
            if (oscChannel->x[i] > sigmin && oscChannel->x[i] < sigmax) if (-1.*oscChannel->y[i] > thres) ++nthres;
            //. if (oscChannel->x[i] >= bkgmin && oscChannel->x[i] <= bkgmax) h->Fill(-oscChannel->y[ic]);
            if (oscChannel->x[i] > sigmin && oscChannel->x[i] < sigmax) if (-1.*oscChannel->y[i] > ampl_max) ampl_max = -1.*oscChannel->y[i];

            // fill the graph
            x[np] = oscChannel->x[i];
            y[np] = -1.*oscChannel->y[i] - bkg_mean;
            ex[np] = 0;
            ey[np] = bkg_sigma;
            ++np;
         }

         // pointers to graphs
         TGraphErrors *grsig = 0;
         TGraphErrors *grsig_refit = 0;

         grsig = new TGraphErrors(np, x, y, ex, ey);
         grsig->SetNameTitle(Form("grsig_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch), Form("grsig_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch));
         grsig->SetMarkerStyle(24);
         grsig->SetMarkerColor(8);
         grsig->SetLineColor(8);
         if (gdebug) {
            new TCanvas;
            grsig->Draw("ap");
         }

         Double_t integral = 0;        // integral from direct sum
         Double_t integralf = 0;       // integral from fit (will be computed later)

         for (int i=0; i<np-1; ++i) integral += 0.5*(grsig->GetY()[i] + grsig->GetY()[i+1]) * (grsig->GetX()[i+1]-grsig->GetX()[i]);

         // normalize to pC by division to R = 50 Hom
         integral /= 50.;
         cout<< ".. sum integral = " << integral << " pC" <<endl;
         oscFit->adc[ich] = integral;

         // apply const threshold
         if (nthres < nthres_min) {
            cout<< "skip event evt " << oscEvent->evt << " channel " << oscChannel->ch << " the number of points which exceed threshold " << thres << " is too small: " << nthres <<endl;
            continue;
         }

         // skip events with saturation
         if (ampl_max > ysaturation) {
            cout<< "skip saturated event evt " << oscEvent->evt << " channel " << oscChannel->ch <<endl;
            continue;
         }

         //
         // fit pulse
         //

         Double_t pulse_max = grsig->GetY()[0];
         // find value and position of the maximum
         for (int i=0; i<np; ++i) {
            if (y[i] > pulse_max) {
               pulse_max = y[i];
            }
         }

         Double_t pulse_halfmax =   0.5*pulse_max;
         Double_t pulse_xhalfmax = 0;
         // Double_t pulse_y1 =        0.25*pulse_max;
         //-- Double_t pulse_y1 =        0.20*pulse_max;
         // Double_t pulse_y1 =        0.10*pulse_max;
         Double_t pulse_y1 =        0.05*pulse_max;
         //////////////////////////////////////////////// Double_t pulse_y2 = 0.75*pulse_max;
         Double_t pulse_y2 = 0.85*pulse_max;
         Double_t pulse_x1 = 0;
         Double_t pulse_x2 = 0;
         Int_t ipulse_halfmax = 0;
         Int_t ipulse_x1 = 0;
         Int_t ipulse_x2 = 0;
         Int_t npoint_fit = 0;

         // find the middle of the leading edge
         for (int i=0; i<np; ++i) {
            if (y[i] > pulse_halfmax) {
               ipulse_halfmax = i;
               pulse_xhalfmax = x[i];
               break;
            }
         }
         // find high limit
         for (int i=ipulse_halfmax; i<np; ++i) {
            if (y[i] > pulse_y2) {
               ipulse_x2 = i;
               pulse_x2 = x[i];
               break;
            }
         }
         // find lower limit
         for (int i=ipulse_halfmax; i>=0; --i) {
            if (y[i] < pulse_y1) {
               ipulse_x1 = i;
               pulse_x1 = x[i];
               break;
            }
         }

         //------------------------ test ------------------------

         pulse_x1 = sigmin;   // approx. see loop below
         ipulse_x1 = -1;
         for (int i=0; i<1024; ++i) {
            if (x[i] < sigmin) continue;
            if (x[i] >= sigmin) {
               ipulse_x1 = i;
               pulse_x1 = x[i];
               break;
            }
         }

         // optimal setup for Photek-240
         A = 1000.;
         //-- x0 = 200.;
         x0 = 60.;
         tau1 = 0.050;
         tau2 = 2.0;
         T = 20;
         sigma = 0.200;

         // Int_t ipar = 0;
         // par[ipar++] = A;
         // par[ipar++] = x0;
         // par[ipar++] = tau1;
         // par[ipar++] = tau2;
         // par[ipar++] = T;
         // par[ipar++] = sigma;

         par[ipar_A] = A;
         par[ipar_x0] = x0;
         par[ipar_tau1] = tau1;
         par[ipar_tau2] = tau2;
         par[ipar_T] = T;
         par[ipar_sigma] = sigma;

         fun_Psigma->SetParameters(par);

         for (int i=0; i<fun_Psigma->GetNpar(); ++i) fun_Psigma->ReleaseParameter(i);

         // fun_Psigma->FixParameter(ipar_tau1, 1.);
         fun_Psigma->FixParameter(ipar_tau1, fun_Psigma->GetParameter(ipar_tau1));
         // fun_Psigma->FixParameter(ipar_sigma, fun_Psigma->GetParameter(ipar_sigma));
         fun_Psigma->FixParameter(ipar_T, fun_Psigma->GetParameter(ipar_T));

         TF1* fitted_function = 0;
         TF1* fitted_function_refit = 0;
         npoint_fit = ipulse_x2 - ipulse_x1 + 1;
         // cout<< ".. fit pulse: pulse_max " << pulse_max << " pulse_xhalmax " << pulse_xhalfmax << " pulse_x1 = " << pulse_x1 << " pulse_x2 = " << pulse_x2 << " npoint_fit = " << npoint_fit <<endl;
         Double_t y0cross = 0;
         Double_t y05cross = 0;
         if (npoint_fit > 2)
         {
            if (gdebug) grsig->Fit(fun_Psigma, "R", "", sigmin,sigmax);
            else grsig->Fit(fun_Psigma,"R0", "goff", sigmin,sigmax);

            fitted_function = grsig->GetFunction(fun_Psigma->GetName());

            // fitted integral
            integralf = fitted_function->Integral(sigmin,1024.) / 50.;
            cout<< ".. function integral = " << integralf << " pC" <<endl;

            // refit
            cout<< "-- refit" <<endl;

            //Double_t sigmax_refit = grsig->GetFunction(fun_Psigma->GetName())->GetMaximumX();

            //-- grsig_refit = (TGraphErrors*) grsig->Clone("grsig_refit");
            //-- grsig_refit->SetTitle(Form("%s refit",grsig_refit->GetTitle()));
            //-- if (gdebug) {
            //--    new TCanvas;
            //--    grsig_refit->Draw("ap");
            //-- }
            grsig_refit = grsig;
            fun_Psigma_refit->SetParameters(grsig->GetFunction(fun_Psigma->GetName())->GetParameters());
            // -- for (int i=0; i<fun_Psigma->GetNpar(); ++i) cout<< i << "\t " << fun_Psigma->GetParName(i) << "\t " << fun_Psigma->GetParameter(i) <<endl;
            
            // fix T for refit
            fun_Psigma_refit->FixParameter(ipar_tau1, grsig->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_tau1));
            fun_Psigma_refit->FixParameter(ipar_tau2, grsig->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_tau2));
            fun_Psigma_refit->FixParameter(ipar_T, grsig->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_T));
            fun_Psigma_refit->FixParameter(ipar_sigma, grsig->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_sigma));

            // fun_Psigma->FixParameter(ipar_tau2, grsig->GetFunction(fun_Psigma->GetName())->GetParameter(3));
            // //fun_Psigma->SetParameter(ipar_tau1, 1.);
            // fun_Psigma->ReleaseParameter(ipar_tau1);
            // fun_Psigma->SetParLimits(ipar_tau1, 0., 10.);

            if (gdebug) grsig_refit->Fit(fun_Psigma_refit, "+R", "", pulse_x1,pulse_x2);
            else grsig_refit->Fit(fun_Psigma_refit, "+R0", "goff", pulse_x1,pulse_x2);

            // swap functions such that the stats box shows the last fitted function
            TF1* frefit = (TF1*) grsig->GetListOfFunctions()->Last();
            grsig->GetListOfFunctions()->Remove(frefit);
            grsig->GetListOfFunctions()->AddFirst(frefit);

            fitted_function_refit = grsig_refit->GetFunction(fun_Psigma_refit->GetName());
            for (int i=0; i<fitted_function_refit->GetNpar(); ++i) cout<< i << "\t " << fitted_function_refit->GetParName(i) << "\t " << fitted_function_refit->GetParameter(i) <<endl;

            y0cross = fitted_function_refit->GetParameter("x0");
            y05cross = 0;
         }
         else {
            cout<< "--- no fit for this channel: noint_fit = " << npoint_fit <<endl;
         }

         // assign fit result tree
         Int_t ichannel = oscChannel->ch - 1;                     // to count array element from 0
         oscFit->adcf[ichannel] = integralf;
         // oscFit->t[ichannel] = y0cross;
         // oscFit->T[ichannel] = y05cross;

         if (fitted_function)
         {
            //. take T and tau2 from the first fit
            //oscFit->T[ichannel] = fitted_function->GetParameter("T");
            //oscFit->tau2[ichannel] = fitted_function->GetParameter("tau2");
            //. take tau1, sigma from the refit
            //oscFit->tau1[ichannel] = fitted_function_refit->GetParameter("tau1");
            //oscFit->sigma[ichannel] = fitted_function_refit->GetParameter("sigma");
            oscFit->v[ichannel] = pulse_max;
            oscFit->bkg[ichannel] = bkg_mean;
            oscFit->sbkg[ichannel] = bkg_sigma;
            if (fitted_function_refit)
            {
               // calculate intersection of tangent to the middle of the leading edge with axis x
               Double_t fval = fitted_function_refit->Eval(pulse_xhalfmax);
               Double_t dval = fitted_function_refit->Derivative(pulse_xhalfmax);
               Double_t dintersect = 0;
               const Double_t eps = 1e-12;
               if (TMath::Abs(dval) > eps) dintersect = pulse_xhalfmax - fval/dval;
               oscFit->d[ichannel] = dintersect;
               cout<< "--> dintersect = " << dintersect << " pulse_xhalfmax = " << pulse_xhalfmax << " fval = " << fval << " dval = " << dval <<endl;

               oscFit->t[ichannel] = fitted_function_refit->GetParameter("x0");
               oscFit->chi2[ichannel] = (fitted_function_refit->GetNDF() > 0)? fitted_function_refit->GetChisquare() / fitted_function_refit->GetNDF(): 0;
               oscFit->tau1[ichannel] = fitted_function_refit->GetParameter("tau1");
               oscFit->tau2[ichannel] = fitted_function->GetParameter("tau2");
               oscFit->T[ichannel] = fitted_function_refit->GetParameter("T");
               oscFit->sigma[ichannel] = fitted_function_refit->GetParameter("sigma");
            }
         }

         // cleanup
         if (!gdebug) {
            if (grsig) delete grsig;
            //-- if (grsig_refit) delete grsig_refit;
         }
      }
      otree->Fill();
   }

   for (int ich=0; ich<4; ++ich) {
      if (ich < 2) continue;
      new TCanvas;
      otree->Draw(Form("adc[%d]",ich), Form("adc[%d]>0&&adc[%d]<1000",ich,ich));
      new TCanvas;
      otree->Draw(Form("chi2[%d]",ich), Form("chi2[%d]>0&&chi2[%d]<10",ich,ich));
   }

   // turn off debug
   debug = false;
   gdebug = false;

   return otree;
}

TTree* pulse(const char *ifname="Co60_STM_LSO2x2_NOCAP_NOAMP_split.xml.root"
      , Double_t bkgmin=0, Double_t bkgmax=80, Double_t sigmin=80, Double_t sigmax=512
      , Int_t entry_first=0, Int_t entry_last=-1
      , bool setdebug=false, bool setgdebug=false
      )
{
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return 0;
   }
   cout<< "Processing file " << ifname <<endl;

   TTree* tree = (TTree*) ifile->Get("t");
   if (!tree) {
      cout<< "tree \"t\" was not found in file " << ifname <<endl;
      return 0;
   }

   // output file with tree "ft" ("fit tree")
   TFile* ofile = TFile::Open(Form("%s-ft.root",ifname),"recreate");

   TTree* otree = pulse(tree, bkgmin,bkgmax, sigmin,sigmax, entry_first,entry_last, setdebug,setgdebug);
   
   ofile->Write();
   return otree;
}
