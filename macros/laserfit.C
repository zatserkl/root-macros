#include "drs.C"

#include <TROOT.h>
#include <TEnv.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TMarker.h>
#include <TCanvas.h>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2.h>

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <map>

using std::cout;     using std::endl;

// utils.C stuff
std::ostream& exitl(std::ostream& os);
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

class PulseParameters: public TNamed {
public:
   Double_t A;
   Double_t x0;
   Double_t tau1;
   Double_t tau2;
   Double_t T;
   Double_t sigma;
public:
   void clear() {A=x0=tau1=tau2=T=sigma=0;}
   PulseParameters(): TNamed() {
      SetNameTitle("pulseParameters", "initParameters");
      clear();
   }
   void SetParameters(Double_t* par) {
      Int_t npar = 0;
      A =      par[npar++];
      x0 =     par[npar++];
      tau1 =   par[npar++];
      tau2 =   par[npar++];
      T =      par[npar++];
      sigma =  par[npar++];
   }
   void GetParameters(Double_t* par) const {
      Int_t npar = 0;
      par[npar++] = A;
      par[npar++] = x0;
      par[npar++] = tau1;
      par[npar++] = tau2;
      par[npar++] = T;
      par[npar++] = sigma;
   }
   Double_t Get_A() const {return A;}
   Double_t Get_T() const {return T;}
   Int_t GetNpar() const {
      Int_t npar = 6;
      return npar;
   }
   ClassDef(PulseParameters, 2);
};

class OscFit: public TNamed
{
public:
   Int_t evt;
   Float_t adc[8];
   Float_t adcf[8];        // ADC from fit
   Float_t t[8];
   Float_t d[8];
   Float_t dydx[8];        // derivative
   Float_t tau1[8];
   Float_t tau2[8];
   Float_t T[8];
   Float_t sigma[8];
   Float_t bkg[8];         // flat background before the signal
   Float_t sbkg[8];        // sigma of the flat background
   Float_t v[8];
   Float_t xmaxi[8];       // beginning of the local maximum
   Float_t ymaxi[8];       // beginning of the local maximum
   Float_t xt[8];          // x-coordinate of point to draw tangent
   Float_t yt[8];          // y-coordinate of point to draw tangent
   Float_t chi2[8];
   Float_t np[8];
   // line fit
   Float_t p0[8];          // slope in y = a*x + b
   Float_t p1[8];          // intercept
   Float_t dline[8];       // intersection of the fitted line with y-axis
   Float_t chi2line[8];    // chi2 of line fit
   Float_t npline[8];
public:
   void clear() {
      evt = 0;
      for (int i=0; i<8; ++i) {
         adc[i] = 0;
         adcf[i] = 0;
         t[i] = 0;
         d[i] = 0;
         dydx[i] = 0;
         tau1[i] = 0;
         tau2[i] = 0;
         T[i] = 0;
         sigma[i] = 0;
         bkg[i] = 0;
         sbkg[i] = 0;
         v[i] = 0;
         xmaxi[i] = 0;
         ymaxi[i] = 0;
         xt[i] = 0;
         yt[i] = 0;
         chi2[i] = 10;
         np[i] = 0;
         p0[i] = 0;
         p1[i] = 0;
         dline[i] = 0;
         chi2line[i] = 10;
         npline[i] = 0;
      }
   }
   OscFit(): TNamed() {
      SetNameTitle("oscFit", "oscFit");
      clear();
   }
   ClassDef(OscFit,14);
};

#ifdef __MAKECINT__
#pragma link C++ class OscFit;
#pragma link C++ class PulseParameters;
#endif

ClassImp(OscFit);
ClassImp(PulseParameters);

void PrintPulseParameters(TTree* ft)
{
   OscFit* oscFit = 0;
   PulseParameters* initParameters = 0;

   ft->SetBranchAddress("oscFit", oscFit);
   cout<< "ft->GetEntries() = " << ft->GetEntries() <<endl;

   //ft->GetEntry(0);
   cout<< "GetUserInfo" <<endl;

   cout<< "ft->GetUserInfo()->GetEntries() = " << ft->GetUserInfo()->GetEntries() <<endl;

   if (ft->GetUserInfo()->GetEntries() > 0)
   {
      initParameters = (PulseParameters*) ft->GetUserInfo()->At(0);
      cout<< "A = " << initParameters->A <<endl;
      cout<< "x0 = " << initParameters->x0 <<endl;
      cout<< "tau1 = " << initParameters->tau1 <<endl;
      cout<< "tau2 = " << initParameters->tau2 <<endl;
      cout<< "T = " << initParameters->T <<endl;
      cout<< "sigma = " << initParameters->sigma <<endl;
   }
}

//--------------------------------------------------------------------------------

Double_t fline(Double_t xx[], Double_t par[])
{
   Double_t& x = *xx;
   Double_t& intercept = par[0];
   Double_t& slope = par[1];
   return intercept + slope*x;
}

/////////////////// pulse function begin ///////////////////////////

Double_t ITtau(Double_t x, Double_t tau, Double_t T)
{
   const Double_t eps = 1e-12;
   const Double_t explim = 100.;
   Double_t arg;

   if (x < 0) return 0;                // function range
   if (tau < 0 || T < 0) return 0;     // sanity check

   Double_t ITtau = 0;

   // case of tau --> 0
   if (tau < eps) return 0;
   
   // case of T --> 0
   if (T < eps) {
      arg = -x/tau;
      if (TMath::Abs(arg) < explim) ITtau = TMath::Exp(arg);
      return ITtau;
   }

   // case of T --> tau
   if (TMath::Abs(tau - T) < eps) {
      // NB: tau > eps here
      arg = -x/tau;
      if (TMath::Abs(arg) < explim) ITtau = (x/tau)*TMath::Exp(arg);
      return ITtau;
   }

   // general case
   Double_t exp_T = 0;
   arg = -x/T;
   if (TMath::Abs(arg) < explim) exp_T = TMath::Exp(arg);
   Double_t exp_tau = 0;
   arg = -x/tau;
   if (TMath::Abs(arg) < explim) exp_tau = TMath::Exp(arg);
   // NB: T-tau is finite here
   ITtau = (tau/(T-tau)) * (exp_T - exp_tau);

   return ITtau;
}

Double_t ITtausigma(Double_t x, Double_t tau, Double_t T, Double_t sigma)
{
   const Double_t eps = 1e-12;
   const Double_t explim = 100.;
   const Double_t erfclim = 100.;
   Double_t arg;

   if (tau < 0 || T < 0) return 0;                 // sanity check
   if (sigma < 0) return 0;                        // is that possible? -- Yes!

   if (sigma == 0) return ITtau(x,tau,T);    // exect solution for sigma = 0

   // case tau --> 0
   if (tau < eps) return 0;

   // case T --> 0
   if (T < eps)
   {
      Double_t exp_x = 0;
      arg = -x/tau;
      if (TMath::Abs(arg) < explim) exp_x = TMath::Exp(arg);
      Double_t exp_sigma2 = 0;
      arg = sigma*sigma/2/tau/tau;
      if (TMath::Abs(arg) < explim) exp_sigma2 = TMath::Exp(arg);

      Double_t erfc_xsigmatau = 0;
      if (sigma*erfclim > TMath::Abs((x-sigma*sigma/tau)/sqrt(2.)))
      {
         arg = -(x-sigma*sigma/tau)/(sigma*sqrt(2.));
         erfc_xsigmatau = TMath::Erfc(arg);
      }
      else if (x-sigma*sigma/tau > 0) erfc_xsigmatau = 2;   // negative argument
      else erfc_xsigmatau = 0;                              // positive argument

      Double_t ITtausigma = 0.5*exp_x*exp_sigma2*erfc_xsigmatau;
      return ITtausigma;
   }

   // case tau --> T
   if (TMath::Abs(T-tau) < eps)
   {
      Double_t tt = tau > T? tau: T;

      Double_t erfc_xsigmatt = 0;
      if (sigma*erfclim > TMath::Abs((x-sigma*sigma/tt)/sqrt(2.)))
      {
         arg = -(x-sigma*sigma/tt)/(sigma*sqrt(2.));
         erfc_xsigmatt = TMath::Erfc(arg);
      }
      else if (x-sigma*sigma/tt > 0) erfc_xsigmatt = 2;  // negative argument
      else erfc_xsigmatt = 0;                            // positive argument

      Double_t exp_xsigmatt2 = 0;
      if (sigma*TMath::Sqrt(erfclim) > (x-sigma*sigma/tt)/sqrt(2.))
      {
         arg = (x-sigma*sigma/tt)/(sigma*sqrt(2.));
         exp_xsigmatt2 = TMath::Exp(-arg*arg);
      }

      Double_t term1 = sigma*exp_xsigmatt2;
      Double_t sqrt_pi = TMath::Sqrt(TMath::Pi());
      Double_t sqrt_2 = TMath::Sqrt(2.);
      Double_t term2 = sqrt_pi*(x-sigma*sigma/tt)/sqrt_2*erfc_xsigmatt;
      Double_t exp_x = 0;
      arg = -x/tt;
      if (TMath::Abs(arg) < explim) exp_x = TMath::Exp(arg);
      Double_t exp_sigma2 = 0;
      arg = sigma*sigma/2/tt/tt;
      if (TMath::Abs(arg) < explim) exp_sigma2 = TMath::Exp(arg);
      Double_t inv_sqrt_2pi = (1./TMath::Sqrt(2.*TMath::Pi()));
      Double_t ITtausigma = inv_sqrt_2pi*exp_x*exp_sigma2*(1./tt)*(term1+term2);
      return ITtausigma;
   }

   // regular case

   // T term
   Double_t exp_xT = 0;
   arg = -x/T;
   if (TMath::Abs(arg) < explim) exp_xT = TMath::Exp(arg);
   //
   Double_t exp_sigma2T2 = 0;
   arg = sigma*sigma / (2*T*T);
   if (TMath::Abs(arg) < explim) exp_sigma2T2 = TMath::Exp(arg);
   //
   Double_t erfc_xsigmaT = 0;
   if (sigma*erfclim > TMath::Abs((x-sigma*sigma/T)/sqrt(2.)))
   {
      arg = -(x-sigma*sigma/T)/(sigma*sqrt(2.));
      erfc_xsigmaT = TMath::Erfc(arg);
   }
   else if (x-sigma*sigma/T > 0) erfc_xsigmaT = 2;       // negative argument
   else erfc_xsigmaT = 0;                                // positive argument
   //
   Double_t term_T = 0.5*exp_xT*exp_sigma2T2*erfc_xsigmaT;

   // tau term
   Double_t exp_xtau = 0;
   arg = -x/tau;
   if (TMath::Abs(arg) < explim) exp_xtau = TMath::Exp(arg);
   //
   Double_t exp_sigma2tau2 = 0;
   arg = sigma*sigma / (2*tau*tau);
   if (TMath::Abs(arg) < explim) exp_sigma2tau2 = TMath::Exp(arg);
   //
   Double_t erfc_xsigmatau = 0;
   if (sigma*erfclim > TMath::Abs((x-sigma*sigma/tau)/sqrt(2.)))
   {
      arg = -(x-sigma*sigma/tau)/(sigma*sqrt(2.));
      erfc_xsigmatau = TMath::Erfc(arg);
   }
   else if (x-sigma*sigma/tau > 0) erfc_xsigmatau = 2;   // negative argument
   else erfc_xsigmatau = 0;                              // positive argument
   //
   Double_t term_tau = 0.5*exp_xtau*exp_sigma2tau2*erfc_xsigmatau;

   Double_t ITtausigma = (tau/(T-tau))*(term_T - term_tau);
   return ITtausigma;
}

Double_t fPsigma(Double_t *xx, Double_t *par)
{
   const Double_t eps = 1e-12;

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
   if (sigma < 0) return 0;                        // is that possible? -- Yes!

   if (x <= -10.*TMath::Abs(sigma)) return 0;
   
   Double_t norm = A*(tau1+tau2)/tau2/tau2;
   Double_t tau12 = tau1*tau2/(tau1+tau2);

   Double_t fPsigma = norm * (ITtausigma(x,tau2,T,sigma) - ITtausigma(x,tau12,T,sigma));
   return fPsigma;
}

/////////////////// pulse function end ///////////////////////////

/////////////////// FitFun begin ///////////////////////

// class FitFun {
// public:
//    TF1* ffit_;
//    TF1* frefit_;
//    Double_t par_[10];
//    Int_t npar_;
//    std::vector<int> fix_fit;
//    std::vector<int> fix_refit;
//    void clear() {
//       for (int i=0; i<10; ++i) par[i] = 0;
//    }
//    FitFun(): ffit_(0), frefit_(0)
//              , npar_(0)
//    {
//       clear();
//    }
//    FitFun(TF1* ffit_, TF1* frefit): ffit_(f), frefit_(0)
//                    , npar_(0)
//    {
//       clear();
//    }
//    void FixParFit(Int_t ipar) {
//       fix_fit.push_back(ipar);
//    }
//    void FixParFit(const char* parname) {
//       cout<< "FixFit " << parname <<endl;
//    }
//    void FixParRefit(Int_t ipar) {
//       fix_refit.push_back(ipar);
//    }
//    void FixParRefit(const char* parname) {
//       cout<< "FixRefit " << parname <<endl;
//    }
//    void PrepareFit() {
//       for (unsigned ipar=0; ipar<npar_; ++ipar) ffit_->ReleaseParameter(ipar);
//       for (unsigned i=0; i<fix_fit.size(); ++i) {
//          ffit_->FixParameter(fix_fit[i], ffit_->GetParameter(fix_fit[i]));
//       }
//    }
// };

class Fit {
private:
   std::map<std::string, int> index_;
   TF1* ffit_;
   Double_t par_[10];
   Int_t npar_;
   std::vector<std::string> fix_;
public:
   //
   // redefine input stream operator in my non-standard way to use like:
   // fit[0] << fun_Psigma << "A" << "x0" << "tau1" << "tau2" << "T" << "sigma";
   //
   Fit& operator <<(TF1* f) {
      if (!ffit_) ffit_ = f;
      else cout<< "***Warning Fit: fit function (ffit_) is alreaky filled" <<endl;
      return *this;
   }
   Fit& operator <<(const std::string& parname) {
      index_[parname] = index_.size() - 1;      // NB: size() - 1
      return *this;
   }
   Int_t Index(const std::string& parname) const {
      std::map<std::string, int>::const_iterator it = index_.find(parname);
      return it->second;
   }
   void Set(const std::string& parname, Double_t parvalue, Bool_t fix=0) {
      if (fix) ffit_->FixParameter(Index(parname), parvalue);
      else ffit_->SetParameter(Index(parname),parvalue);
   }
   void Fix(const std::string& parname) {
      //ffit_->FixParameter(Index(parname))
      fix_.push_back(parname);
   }
   void clear() {
      npar_ = 0;
      for (int i=0; i<10; ++i) par_[i] = 0;
   }
   Fit(): ffit_(0)
          , npar_(0)
   {
      clear();
   }
   void PrepareFit(Fit* fitter=0)
   {
      TF1* fun = 0;
      if (fitter) {
         fun = fitter->ffit_;
         for (int i=0; i<fun->GetNpar(); ++i) par_[i] = fun->GetParameter(i);
         ffit_->SetParameters(par_);
      }
      for (int ipar=0; ipar<npar_; ++ipar) ffit_->ReleaseParameter(ipar);
      // for (int i=0; i<fix_.size(); ++i) {
      //    ffit_->FixParameter(index_[i], ffit_->GetParameter(fix_[i]));
      // }
      for (std::map<std::string, int>::const_iterator it=index_.begin(); it!=index_.end(); ++it) {
         Int_t ipar = it->second;
         ffit_->FixParameter(ipar, ffit_->GetParameter(ipar));
      }
   }
};

// class FitBoard0 {
// private:
//    static const Int_t nchan_max_ = 8;
//    unsigned nchan_;
//    Fit fit_[nchan_max_+1];                // max number of channels + default channel
//    Fit fit_out_of_range_;
// public:
//    FitBoard0(Int_t nchan): nchan_(nchan) {}
//    Fit& operator [](unsigned chan) {
//       if (chan < nchan_) return fit_[chan];
//       else {
//          cout<< "***ERROR: channel " << chan << " is out of range. The maximum channel is " << nchan_ <<endl;
//          return fit_out_of_range_;
//       }
//    }
//    const Fit& operator [](unsigned chan) const {
//       if (chan < nchan_) return fit_[chan];
//       else {
//          cout<< "***ERROR: channel " << chan << " is out of range. The maximum channel is " << nchan_ <<endl;
//          return fit_out_of_range_;
//       }
//    }
//    void PrepareFit() {
//       cout<< "FitBoard0::Preapare" <<endl;
//    }
//    void Set(Int_t chan, const std::string& parname, Double_t parvalue, Bool_t fix=0) {
//       fit_[chan].Set(parname, parvalue, fix);
//    }
// };

class Fitter {
protected:
   TF1* tf1_;                             //!
   std::map<std::string, int> index_;
   std::vector<std::string> fix_;
   Double_t par_[10];
   Double_t epar_[10];
   Int_t npar_;
   Double_t xmin_;
   Double_t xmax_;
   // fit results
   Double_t chi2_;
public:
   //
   // Finally I will have class FitBoard with array of pointer to Fitter.
   // The first element will be default fitter: child class of the Fitter.
   // I'd like to clone default class to the rest of the channels.
   // I cannot do that using a[i] = new FitX((*a)[0]) because
   //    -- I do not know type FitX: it can be FitPulse of FitLine, etc.
   //    -- I cannot use Fitter((*a)[0]) because it will create object just of type Fitter
   //       (because virtual constructors do not exist)
   // Therfore I need a virtual function e.g. Clone() which will return pointer on cloned child object
   // To prevent creation of the object of type Fitter method Fitter::Clone can be made pure virtual.
   // Thus I need:
   // 2) copy constructor
   // 3) virtual method Clone which creates new Fitter object by copy constructor
   //
   void clear() {
      for (int i=0; i<10; ++i) {
         par_[i] = 0;
         epar_[i] = 0;
      }
   }
   Fitter(): tf1_(0)
             , npar_(0)
             , xmin_(0)
             , xmax_(0)
   {
      clear();
   }
   Fitter(Double_t xmin, Double_t xmax): tf1_(0)
                                         , npar_(0)
                                         , xmin_(xmin)
                                         , xmax_(xmax)
   {
      clear();
   }
   Fitter(const Fitter& fitter): tf1_(fitter.tf1_)
                                 , npar_(fitter.npar_)
                                 , xmin_(fitter.xmin_)
                                 , xmax_(fitter.xmax_)
   {
      for (int i=0; i<fitter.npar_; ++i) {
         par_[i] = fitter.par_[i];
         epar_[i] = fitter.epar_[i];
      }
   }
   virtual Fitter* Clone() const = 0;
   Fitter& operator <<(const std::string& parname) {
      index_[parname] = index_.size() - 1;      // NB: size() - 1
      return *this;
   }
   // virtual Fitter& operator =(const Fitter& fitter) {
   //    tf1_ = fitter.tf1_;
   //    return *this;
   // }
   Int_t Index(const std::string& parname) const {
      std::map<std::string, int>::const_iterator it = index_.find(parname);
      return it->second;
   }
   Bool_t Exist(const std::string parname) const {
      // look for key
      return index_.find(parname) != index_.end();
   }
   Bool_t Exist(Int_t i) const {
      // look for value
      for (std::map<std::string, int>::const_iterator it=index_.begin(); it!=index_.end(); ++it) {
         if (it->second == i) return kTRUE;
      }
      return kFALSE;
   }
   TF1* FitFun()  // use this for fit
   {
      // ensure that the fit function is ready to use: parameter have been set, etc.
      if (true
            && index_.size() != 0
         )
         return tf1_;
      else {
         cout<< "***Error Fitter::FitFun: the function is not initialized" <<endl;
         return 0;
      }
   }
   void Fit(TH1* h, Option_t* option, Option_t* goption, Double_t xmin, Double_t xmax) {
      h->Fit(tf1_, option, goption, xmin, xmax);
      for (int i=0; i<npar_; ++i) {
         par_[i] = tf1_->GetParameter(i);
         epar_[i] = tf1_->GetParError(i);
      }
   }
   void Fit(TGraph* g, Option_t* option, Option_t* goption, Double_t xmin, Double_t xmax) {
      g->Fit(tf1_, option, goption, xmin, xmax);
      for (int i=0; i<npar_; ++i) {
         par_[i] = tf1_->GetParameter(i);
         epar_[i] = tf1_->GetParError(i);
      }
   }
   void Set(const std::string& parname, Double_t parvalue, Bool_t fix=0) {
      if (!Exist(parname)) {
         cout<< "***Error Fitter::Fix: unknown parameter " << parname <<endl;
         exit(0);
      }
      par_[Index(parname)] = parvalue;
      if (fix) fix_.push_back(parname);
      else Release(parname);
   }
   void Fix(const std::string& parname) {
      if (!Exist(parname)) {
         cout<< "***Error Fitter::Fix: unknown parameter " << parname <<endl;
         exit(0);
      }
      fix_.push_back(parname);
   }
   void FixAll() {
      fix_.clear();
      for (std::map<std::string, int>::const_iterator it=index_.begin(); it!=index_.end(); ++it) {
         fix_.push_back(it->first);
      }
   }
   void Release(const std::string& parname) {
      cout<< "Fitter::Release" <<endl;
      if (!Exist(parname)) {
         cout<< "***Error Fitter::Fix: unknown parameter " << parname <<endl;
         exit(0);
      }
      for (std::vector<std::string>::iterator it=fix_.begin(); it!=fix_.end(); ++it) {
         if (*it == parname) {
            it = fix_.erase(it);
            if (it == fix_.end()) break;     // prevents crash at erasing of the last element
         }
      }
   }
   void ReleaseAll()
   {
      fix_.clear();
   }
   virtual void PrepareFit(Fitter* fitter=0)
   {
      TF1* fun = 0;
      if (fitter) {
         fun = fitter->tf1_;
         if (fun) {
            for (int i=0; i<fun->GetNpar(); ++i) par_[i] = fun->GetParameter(i);
            tf1_->SetParameters(par_);
         }
         else {
            cout<< "***Error Fitter::PrepareFit: tf1_ == 0" <<endl;
            exit(0);
         }
      }
      for (std::map<std::string, int>::const_iterator it=index_.begin(); it!=index_.end(); ++it) {
         tf1_->SetParName(it->second, it->first.c_str());
      }
      tf1_->SetParameters(par_);
      for (int ipar=0; ipar<npar_; ++ipar) tf1_->ReleaseParameter(ipar);
      for (unsigned i=0; i<fix_.size(); ++i) {
         tf1_->FixParameter(Index(fix_[i]), tf1_->GetParameter(Index(fix_[i])));
      }
      Complete();
   }
   virtual Bool_t Complete() const        // NB: virtual
   {
      //cout<< "Fitter::Complete" <<endl;
      if (!tf1_) cout<< "***Error Fitter::Complete: tf1_ == 0" <<endl<<exitl;
      if (xmin_ >= xmax_) cout<< "***Error Fitter::Complete: xmin_ >= xmax_" <<endl<<exitl;
      if (index_.size()==0) cout<< "*Error Fitter::Complete: index_.size()==0" <<endl<<exitl;
      return kTRUE;
   }
};

class FitPulse: public Fitter {
public:
   FitPulse(Double_t xmin, Double_t xmax): Fitter(xmin,xmax)
   {
      npar_ = 6;
      *this << "A" << "x0" << "tau1" << "tau2" << "T" << "sigma";
      cout<< "FitPulse::FitPulse: index_.size() = " << index_.size() <<endl;
      tf1_ = new TF1("fun_Psigma", fPsigma, xmin_, xmax_, npar_);
      tf1_->SetNpx(1000);
   }
   FitPulse(const FitPulse& fitter): Fitter(fitter)
   {
      npar_ = fitter.npar_;
      *this << "A" << "x0" << "tau1" << "tau2" << "T" << "sigma";
      cout<< "FitPulse::FitPulse: index_.size() = " << index_.size() <<endl;
      tf1_ = fitter.tf1_;
      tf1_->SetNpx(1000);

      index_ = fitter.index_;
      fix_ = fitter.fix_;
      for (int i=0; i<npar_; ++i) {
         par_[i] = fitter.par_[i];
         epar_[i] = fitter.epar_[i];
      }
   }
   FitPulse* Clone() const {
      FitPulse* fitPulse = new FitPulse(*this);
      return fitPulse;
   }
   Bool_t Complete() const
   {
      Fitter::Complete();
      //cout<< "FitPulse::Complete" <<endl;
      if (!Exist("A")) cout<< "***Error FitPulse::Complete: missing parameter A" <<endl<<exitl;
      if (!Exist("x0")) cout<< "***Error FitPulse::Complete: missing parameter x0" <<endl<<exitl;
      if (!Exist("tau1")) cout<< "***Error FitPulse::Complete: missing parameter tau1" <<endl<<exitl;
      if (!Exist("tau2")) cout<< "***Error FitPulse::Complete: missing parameter tau2" <<endl<<exitl;
      if (!Exist("T")) cout<< "***Error FitPulse::Complete: missing parameter T" <<endl<<exitl;
      if (!Exist("sigma")) cout<< "***Error FitPulse::Complete: missing parameter sigma" <<endl<<exitl;
      return kTRUE;
   }
   // FitPulse& operator =(const FitPulse& fitter) {
   //    Fitter::operator=(fitter);    // that is: call operator = from the parent class Fitter
   //    return *this;
   // }
};

class FitBoard {
private:
   static const Int_t nchan_max_ = 8;
   unsigned nchan_;
   //Fitter fit_out_of_range_;
   // control flags
   Bool_t defaultComplete_;
   Bool_t defaultClone_;          // required to be set before setting up of some particular channel
   Bool_t defaultLock_;
public:
   Fitter* fit_[nchan_max_+1];      // max number of channels + default channel
public:
   FitBoard(Int_t nchan): nchan_(nchan)
                          , defaultComplete_(kFALSE)
   {
      for (int i=0; i<nchan_max_+1; ++i) fit_[i] = 0;
   }
   Bool_t CompleteDefault() const {
      return fit_[0]->Complete();
   }
   void PrepareFit() {
      cout<< "FitBoard::Prepare" <<endl;
      if (!defaultClone_) cout<< "FitBoard::PrepareFit: default did not applied yet" <<endl<<exitl;
   }
   Fitter* operator [](Int_t chan) {
      return fit_[chan];
   }
   // void Set(Int_t chan, Fitter* fitter) {
   //    fit_[chan] = new Fitter(*fitter);
   // }
   void Set(Int_t chan, const std::string& parname, Double_t parvalue, Bool_t fix=0) {
      assert(fit_[chan] != 0);
      // if (!fit_[chan]) cout<< "pointer to channel " << chan << " == 0" <<endl<<exitl;
      if (!CompleteDefault() && chan != 0) cout<< "***Warning: you must complete default fitter first" <<endl<<exitl;
      fit_[chan]->Set(parname, parvalue, fix);
   }
   void DefaultClone() {
      cout<< "DefaultClone" <<endl;
      fit_[0]->PrepareFit();
      for (unsigned i=1; i<nchan_; ++i) fit_[i] = fit_[0]->Clone();
      defaultClone_ = kTRUE;
      defaultLock_ = kTRUE;
   }
};

void laserfit_alone() {
   cout<< "test of laserfit" <<endl;
   const char ifname[] = "stm_noC_ch123_trigext_5GS.root.osc.root";
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return;
   }
   TTree* t = (TTree*) ifile->Get("t");
   if (!t) {
      cout<< "No tree t was found in file " << ifname <<endl;
      return;
   }

   Int_t evt = 0;
   Int_t channel = 1;

   OscEvent* oscEvent = 0;
   t->SetBranchAddress("oscEvent",&oscEvent);

   t->GetEntry(evt);

   const OscChannel* oscChannel = (OscChannel*) oscEvent->oscChannels->At(channel-1);
   cout<< "oscChannel->ch = " << oscChannel->ch <<endl;

   Double_t sigmin = 45;
   Double_t sigmax = 200;
   // Double_t sigmax = 60;
   // Double_t sigmax = 50;
   
   Double_t bkg_mean_arr[] = {   1.60,    1.80,    1.96};
   Double_t bkg_sigma_arr[] = {  0.58,    0.54,    0.54};
   Double_t bkg_mean = bkg_mean_arr[channel-1];
   Double_t bkg_sigma = bkg_sigma_arr[channel-1];

   Double_t x[1024], y[1024], ey[1024];
   Int_t np = 0;

   for (int i=0; i<1024; ++i) {
      if (oscChannel->x[i] < sigmin) continue;
      if (oscChannel->x[i] > sigmax) break;
      x[np] = oscChannel->x[i];
      y[np] = (-1.)*oscChannel->y[i] - bkg_mean;
      ey[np] = bkg_sigma;
      ++np;
   }
   TGraphErrors* gr = new TGraphErrors(np, x, y, 0, ey);
   gr->SetNameTitle("gr", Form("evt %d, ch %d", evt, channel));
   gr->SetMarkerStyle(24);
   gr->SetMarkerColor(8);
   gr->SetLineColor(8);
   new TCanvas;
   gr->Draw("ap");

   Bool_t fix = kTRUE;
   Bool_t nofix = kFALSE;
   nofix = !fix;           // to avoid warning: unused variable ‘nofix’

   FitPulse* fitPulse = new FitPulse(sigmin,sigmax);
   // fitPulse->Set("A",       1.37309e+04,   fix);
   // fitPulse->Set("x0",      4.86460e+01,   fix);
   // fitPulse->Set("tau1",    3.43145e-01,   fix);
   // fitPulse->Set("tau2",    7.46879e+01,   fix);
   // fitPulse->Set("T",       0.00000e+00,   fix);
   // fitPulse->Set("sigma",   2.00000e-01,   fix);
   //
   fitPulse->Set("A",       10000.);
   fitPulse->Set("x0",      50.);
   //fitPulse->Set("tau1",    0.050,   fix);
   // fitPulse->Set("tau1",    0.050);
   fitPulse->Set("tau1",    0.100,   fix);
   fitPulse->Set("tau2",    75.);
   fitPulse->Set("T",       0.017,   fix);
   fitPulse->Set("sigma",   0.200,   fix);
   //fitPulse->Set("sigma",   0.000,   fix);
   //fitPulse->Set("sigma",   2e-3,    nofix);
   //fitPulse->Fix("tau2");

   fitPulse->PrepareFit();
   fitPulse->Fit(gr, "", "", sigmin,sigmax);

   // cout<<endl<< "-- FitBoard test" <<endl;

   // Int_t nchan = 4;
   // FitBoard fit(nchan);

   // channel = 0;
   // fit.fit_[channel] = new FitPulse(sigmin,sigmax);
   // cout<< "call fit.Set for A" <<endl;
   // fit.Set(channel, "A",       10000.);
   // cout<< "call fit.Set for x0" <<endl;
   // fit.Set(channel, "x0",      50.);
   // //fit.Set(channel, "tau1",    0.050,   fix);
   // // fit.Set(channel, "tau1",    0.050);
   // fit.Set(channel, "tau1",    0.100,   fix);
   // fit.Set(channel, "tau2",    75.);
   // fit.Set(channel, "T",       0.017,   fix);
   // fit.Set(channel, "sigma",   0.200,   fix);
   // //fit.Set(channel, "sigma",   0.000,   fix);
   // //fit.Set(channel, "sigma",   2e-3,    nofix);
   // //fitPulse.Fix("tau2");
   // fit.DefaultClone();

   // cout<< "call fit.PrepareFit()" <<endl;
   // fit.PrepareFit();

   // // fit[0]->PrepareFit();
   // // fit[0]->Fit(gr, "", "", sigmin, sigmax);
   // fit[1]->PrepareFit();
   // fit[1]->Fit(gr, "", "", sigmin, sigmax);
}

void laserfit() {
   cout<< "test of laserfit" <<endl;
   const char ifname[] = "stm_noC_ch123_trigext_5GS.root.osc.root";
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return;
   }
   TTree* t = (TTree*) ifile->Get("t");
   if (!t) {
      cout<< "No tree t was found in file " << ifname <<endl;
      return;
   }

   Int_t evt = 0;
   Int_t channel = 1;

   OscEvent* oscEvent = 0;
   t->SetBranchAddress("oscEvent",&oscEvent);

   t->GetEntry(evt);

   const OscChannel* oscChannel = (OscChannel*) oscEvent->oscChannels->At(channel-1);
   cout<< "oscChannel->ch = " << oscChannel->ch <<endl;

   Double_t sigmin = 45;
   Double_t sigmax = 200;
   // Double_t sigmax = 60;
   // Double_t sigmax = 50;
   
   Double_t bkg_mean_arr[] = {   1.60,    1.80,    1.96};
   Double_t bkg_sigma_arr[] = {  0.58,    0.54,    0.54};
   Double_t bkg_mean = bkg_mean_arr[channel-1];
   Double_t bkg_sigma = bkg_sigma_arr[channel-1];

   Double_t x[1024], y[1024], ey[1024];
   Int_t np = 0;

   for (int i=0; i<1024; ++i) {
      if (oscChannel->x[i] < sigmin) continue;
      if (oscChannel->x[i] > sigmax) break;
      x[np] = oscChannel->x[i];
      y[np] = (-1.)*oscChannel->y[i] - bkg_mean;
      ey[np] = bkg_sigma;
      ++np;
   }
   TGraphErrors* gr = new TGraphErrors(np, x, y, 0, ey);
   gr->SetNameTitle("gr", Form("evt %d, ch %d", evt, channel));
   gr->SetMarkerStyle(24);
   gr->SetMarkerColor(8);
   gr->SetLineColor(8);
   new TCanvas;
   gr->Draw("ap");

   Bool_t fix = kTRUE;
   Bool_t nofix = kFALSE;
   nofix = !fix;           // to avoid warning: unused variable ‘nofix’

   // FitPulse* fitPulse = new FitPulse(sigmin,sigmax);
   // // fitPulse->Set("A",       1.37309e+04,   fix);
   // // fitPulse->Set("x0",      4.86460e+01,   fix);
   // // fitPulse->Set("tau1",    3.43145e-01,   fix);
   // // fitPulse->Set("tau2",    7.46879e+01,   fix);
   // // fitPulse->Set("T",       0.00000e+00,   fix);
   // // fitPulse->Set("sigma",   2.00000e-01,   fix);
   // //
   // fitPulse->Set("A",       10000.);
   // fitPulse->Set("x0",      50.);
   // //fitPulse->Set("tau1",    0.050,   fix);
   // // fitPulse->Set("tau1",    0.050);
   // fitPulse->Set("tau1",    0.100,   fix);
   // fitPulse->Set("tau2",    75.);
   // fitPulse->Set("T",       0.017,   fix);
   // fitPulse->Set("sigma",   0.200,   fix);
   // //fitPulse->Set("sigma",   0.000,   fix);
   // //fitPulse->Set("sigma",   2e-3,    nofix);
   // //fitPulse->Fix("tau2");

   // fitPulse->PrepareFit();
   // fitPulse->Fit(gr, "", "", sigmin,sigmax);

   cout<<endl<< "-- FitBoard test" <<endl;

   Int_t nchan = 4;
   FitBoard fit(nchan);

   channel = 0;
   fit.fit_[channel] = new FitPulse(sigmin,sigmax);
   cout<< "call fit.Set for A" <<endl;
   fit.Set(channel, "A",       10000.);
   cout<< "call fit.Set for x0" <<endl;
   fit.Set(channel, "x0",      50.);
   //fit.Set(channel, "tau1",    0.050,   fix);
   // fit.Set(channel, "tau1",    0.050);
   fit.Set(channel, "tau1",    0.100,   fix);
   fit.Set(channel, "tau2",    75.);
   fit.Set(channel, "T",       0.014,   fix);
   // fit.Set(channel, "T",       40.,   fix);
   fit.Set(channel, "sigma",   0.200,   fix);
   //fit.Set(channel, "sigma",   0.000,   fix);
   //fit.Set(channel, "sigma",   2e-3,    nofix);
   //fitPulse.Fix("tau2");
   fit.DefaultClone();

   // some specific requirements for the channel 2
   channel = 1;
   fit.Set(channel, "sigma", 0.300, nofix);
   fit[channel]->Release("tau1");
   //fit[channel]->Release("T");

   cout<< "call fit.PrepareFit()" <<endl;
   fit.PrepareFit();

   // fit[0]->PrepareFit();
   // fit[0]->Fit(gr, "", "", sigmin, sigmax);
   fit[1]->PrepareFit();
   fit[1]->Fit(gr, "", "", sigmin, sigmax);
}

void laserfit_works() {
   cout<< "test of laserfit" <<endl;
   const char ifname[] = "stm_noC_ch123_trigext_5GS.root.osc.root";
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return;
   }
   TTree* t = (TTree*) ifile->Get("t");
   if (!t) {
      cout<< "No tree t was found in file " << ifname <<endl;
      return;
   }

   Int_t evt = 0;
   Int_t channel = 1;

   OscEvent* oscEvent = 0;
   t->SetBranchAddress("oscEvent",&oscEvent);

   t->GetEntry(evt);

   const OscChannel* oscChannel = (OscChannel*) oscEvent->oscChannels->At(channel-1);
   cout<< "oscChannel->ch = " << oscChannel->ch <<endl;

   Double_t sigmin = 45;
   Double_t sigmax = 200;
   // Double_t sigmax = 60;
   // Double_t sigmax = 50;
   
   Double_t bkg_mean_arr[] = {   1.60,    1.80,    1.96};
   Double_t bkg_sigma_arr[] = {  0.58,    0.54,    0.54};
   Double_t bkg_mean = bkg_mean_arr[channel-1];
   Double_t bkg_sigma = bkg_sigma_arr[channel-1];

   Double_t x[1024], y[1024], ey[1024];
   Int_t np = 0;

   for (int i=0; i<1024; ++i) {
      if (oscChannel->x[i] < sigmin) continue;
      if (oscChannel->x[i] > sigmax) break;
      x[np] = oscChannel->x[i];
      y[np] = (-1.)*oscChannel->y[i] - bkg_mean;
      ey[np] = bkg_sigma;
      ++np;
   }
   TGraphErrors* gr = new TGraphErrors(np, x, y, 0, ey);
   gr->SetNameTitle("gr", Form("evt %d, ch %d", evt, channel));
   gr->SetMarkerStyle(24);
   gr->SetMarkerColor(8);
   gr->SetLineColor(8);
   new TCanvas;
   gr->Draw("ap");

   Bool_t fix = kTRUE;
   Bool_t nofix = kFALSE;
   fix = !nofix;           // to avoid warning: unused variable ‘nofix’

   FitPulse fitPulse(sigmin,sigmax);
   // fitPulse.Set("A",       1.37309e+04,   fix);
   // fitPulse.Set("x0",      4.86460e+01,   fix);
   // fitPulse.Set("tau1",    3.43145e-01,   fix);
   // fitPulse.Set("tau2",    7.46879e+01,   fix);
   // fitPulse.Set("T",       0.00000e+00,   fix);
   // fitPulse.Set("sigma",   2.00000e-01,   fix);
   //
   fitPulse.Set("A",       10000.);
   fitPulse.Set("x0",      50.);
   //fitPulse.Set("tau1",    0.050,   fix);
   // fitPulse.Set("tau1",    0.050);
   fitPulse.Set("tau1",    0.100,   fix);
   fitPulse.Set("tau2",    75.);
   fitPulse.Set("T",       0.017,   fix);
   fitPulse.Set("sigma",   0.200,   fix);
   //fitPulse.Set("sigma",   0.000,   fix);
   //fitPulse.Set("sigma",   2e-3,    nofix);
   //fitPulse.Fix("tau2");

   fitPulse.PrepareFit();
   fitPulse.Fit(gr, "", "", sigmin,sigmax);
}

// void fitfun()
// {
//    // example of using of FitFun class in main routine
// 
//    Double_t sigmin = 40;
//    Double_t sigmax = 200;
//    Int_t npar = 6;
// 
//    TF1* fun_Psigma = new TF1("fun_Psigma", fPsigma, sigmin, sigmax, npar);
//    TF1* fun_Psigma_refit = new TF1("fun_Psigma_refit", fPsigma, sigmin, sigmax, npar);
// 
//    const Bool_t fix = kTRUE;
// 
//    //const int NCH = 9;      // chan 0: default, chan N: for DRS4 channel N (N=4 for one board, N=8 for two boards)
// 
//    //FitBoard0 fit(4);                // uses operator [] to work with particular channel
//    //FitBoard0 refit(4);                // uses operator [] to work with particular channel
// 
//    FitBoard fit(4);
//    FitBoard refit(4);
// 
//    // prepare fit
// 
//    // default: need to be filled anyway
//    fit[0] << fun_Psigma << "A" << "x0" << "tau1" << "tau2" << "T" << "sigma";
//    fit[0].Set("A",      1000.);
//    fit[0].Set("x0",     50.);
//    fit[0].Set("tau1",   0.050,   fix);
//    fit[0].Set("tau2",   80.);
//    fit[0].Set("T",      40.,     fix);
//    fit[0].Set("sigma",  0.200,   fix);
//    refit[0] << fun_Psigma_refit << "A" << "x0" << "tau1" << "tau2" << "T" << "sigma";
//    refit[0].Fix("tau2");
//    refit[0].Fix("T");
//    refit[0].Fix("sigma");
// 
//    fit[1].Set("sigma",  0.300,   fix);
//    refit[1].Fix("sigma");
// 
//    // specific for channel 3
//    fit[3] << fun_Psigma << "A" << "x0" << "tau1" << "tau2" << "T" << "sigma";
//    fit[3].Set("A",      1000.);
//    fit[3].Set("x0",     50.);
//    fit[3].Set("tau1",   0.050,   fix);
//    fit[3].Set("tau2",   80.);
//    fit[3].Set("T",      40.,     fix);
//    fit[3].Set("sigma",  0.200,   fix);
//    refit[3] << fun_Psigma_refit << "A" << "x0" << "tau1" << "tau2" << "T" << "sigma";
//    refit[3].Fix("tau2");
//    refit[3].Fix("T");
//    refit[3].Fix("sigma");
// 
//    fit.PrepareFit();
//    refit.PrepareFit();
// 
//    // // loop over channels
// 
//    //    fit[ich].PrepareFit();
//    //    // doing fit ...
//    //    refit[ich].PrepareFit(fit[ich]);  // <--------- how to link with fit?
//    //    // doing refit ...
// }

/////////////////// FitFun end ///////////////////////

bool debug = false;
bool gdebug = false;

TTree* pulse(TTree *tree
      , Double_t bkgmin=0, Double_t bkgmax=40, Double_t sigmin=40, Double_t sigmax=120
      , Int_t entry_first=0, Int_t entry_last=-1
      , bool setdebug=false, bool setgdebug=false
      , Int_t ignorech1=-1
      , Int_t ignorech2=-1
      , Int_t ignorech3=-1
      )
{
   if (setdebug) debug = true;
   if (setgdebug) gdebug = true;

   //-- Double_t thres = 20.;
   //-- Double_t thres = 10.;
   Double_t thres = 1.;
   Int_t nthres_min = 5;                  // data are required to contain at least nthres_min points over the thres value
   //--------------Double_t ysaturation = 499;

   tree->SetMarkerStyle(7);
   tree->SetMarkerColor(2);

   OscEvent* oscEvent = 0;
   tree->SetBranchAddress("oscEvent",&oscEvent);

   cout<< "tree->GetEntries() = " << tree->GetEntries() <<endl;

   // output (fit results) tree
   TTree* otree = new TTree("ft", "Fit result tree");
   OscFit* oscFit = new OscFit;
   otree->Branch("oscFit", "OscFit", &oscFit);
   otree->SetMarkerStyle(6);

   // UserInfo: will be filled befor the first fit
   PulseParameters* initParameters = 0;

   // the number of channels in the data
   tree->GetEntry(entry_first);
   Int_t nchannels = oscEvent->oscChannels->GetEntries();
   
   // histogram for background fit
   TH1F *hbkg = new TH1F("hbkg", "hbkg", 400, -100, 100);

   enum {ipar_A=0, ipar_x0, ipar_tau1, ipar_tau2, ipar_T, ipar_sigma};

   Double_t par[10];
   Int_t npar = 0;
   Double_t& A = par[npar++];
   Double_t& x0 = par[npar++];
   Double_t& tau1 = par[npar++];
   Double_t& tau2 = par[npar++];
   Double_t& T = par[npar++];
   Double_t& sigma = par[npar++];

   // function for fit
   TF1* fun_Psigma = new TF1("fun_Psigma", fPsigma, sigmin, sigmax, npar);
   //-- TF1* fun_Psigma = new TF1("fun_Psigma", fPsigma_std, sigmin, sigmax, npar);
   fun_Psigma->SetNpx(1024);
   fun_Psigma->SetName("fun_Psigma");
   fun_Psigma->SetLineColor(1);
   fun_Psigma->SetLineWidth(1);

   fun_Psigma->SetParameters(par);
   fun_Psigma->SetParName(ipar_A, "A");
   fun_Psigma->SetParName(ipar_x0, "x0");
   fun_Psigma->SetParName(ipar_tau1, "tau1");
   fun_Psigma->SetParName(ipar_tau2, "tau2");
   fun_Psigma->SetParName(ipar_T, "T");
   fun_Psigma->SetParName(ipar_sigma, "sigma");

   // function for refit
   TF1* fun_Psigma_refit = new TF1("fun_Psigma_refit", fPsigma, sigmin, sigmax, npar);
   //-- TF1* fun_Psigma_refit = new TF1("fun_Psigma_refit", fPsigma_std, sigmin, sigmax, npar);
   fun_Psigma_refit->SetNpx(1024);
   fun_Psigma_refit->SetName("fun_Psigma_refit");
   fun_Psigma_refit->SetLineColor(2);

   fun_Psigma_refit->SetParameters(par);
   fun_Psigma_refit->SetParName(ipar_A, "A");
   fun_Psigma_refit->SetParName(ipar_x0, "x0");
   fun_Psigma_refit->SetParName(ipar_tau1, "tau1");
   fun_Psigma_refit->SetParName(ipar_tau2, "tau2");
   fun_Psigma_refit->SetParName(ipar_T, "T");
   fun_Psigma_refit->SetParName(ipar_sigma, "sigma");

   //------------ beginning of line fit function definition -------------

   Double_t parline[10];
   Int_t nparline = 0;

   enum {iparline_a0=0, iparline_a1};
   const char* parlinename[] = {"a0", "a1"};

   Double_t a0 = -1000;
   Double_t a1 = 10;
   parline[iparline_a0] = a0;
   parline[iparline_a1] = a1;
   nparline = 2;

   TF1* fun_line = new TF1("fun_line", fline, 0, 1000, nparline);
   fun_line->SetNpx(1024);
   fun_line->SetName("fun_line");
   std::stringstream sstitline;
   sstitline << fun_line->GetName();
   //sstit << fun_line->GetName() << " intercept=" << intercept << " slope=" << slope;
   fun_line->SetTitle(sstitline.str().c_str());
   fun_line->SetLineColor(4);

   fun_line->SetParameters(parline);
   fun_line->SetParName(iparline_a0, parlinename[iparline_a0]);
   fun_line->SetParName(iparline_a1, parlinename[iparline_a1]);

   //------------ end of line fit function definition ------------

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
         cout<< "oscChannel->ch = " << oscChannel->ch <<endl;

         if (oscChannel->ch == ignorech1) continue;
         if (oscChannel->ch == ignorech2) continue;
         if (oscChannel->ch == ignorech3) continue;

         // prepare the data: invert and squeeze

         Double_t x[1024], y[1024], ex[1024], ey[1024];

         for (int i=0; i<1024; ++i) {
            x[i] = oscChannel->x[i];
            y[i] = (-1.)*oscChannel->y[i];
         }

         Double_t ybkg[1024];
         Int_t np_bkg = 0;
         for (int i=0; i<1024; ++i) {
            if (x[i] < bkgmin) continue;
            if (x[i] > bkgmax) break;
            ybkg[np_bkg] = y[i];
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
         if (debug) cout<< "bkg_mean = " << bkg_mean << " bkg_sigma = " << bkg_sigma <<endl;

         Int_t isigmin = 0;
         Int_t isigmax = 0;

         for (int i=0; i<1024; ++i)
         {
            if (x[i] > sigmin && isigmin == 0) isigmin = i;
            if (x[i] < sigmax) isigmax = i;

            y[i] -= bkg_mean;
            ex[i] = 0;
            ey[i] = bkg_sigma;
         }
         if (debug) cout<< "isigmin = " << isigmin << " isigmax = " << isigmax <<endl;

         Int_t nthres = 0;                         // the number of channels which exceed threshold
         Double_t ymax = 0;
         Double_t integral = 0;                    // integral from direct sum
         Double_t dx = x[isigmin+1] - x[isigmin];
         for (int i=isigmin; i<=isigmax; ++i)
         {
            if (y[i] > ymax) ymax = y[i];
            if (y[i] > thres) ++nthres;
            integral += 0.5*(y[i-1] + y[i]) * dx;
         }
         // if (nthres < nthres_min) {
         //    cout<< "--> oscChannel->ch " << oscChannel->ch << ": the number of channel above threshold " << thres << " is less than " << nthres_min <<endl;
         //    continue;
         // }

         // normalize integral to pC by division to R = 50 Ohm
         integral /= 50.;
         cout<< ".. sum integral = " << integral << " pC" <<endl;
         oscFit->adc[ich] = integral;

         if (nthres < nthres_min) {
            cout<< "--> oscChannel->ch " << oscChannel->ch << ": the number of channel above threshold " << thres << " is less than " << nthres_min <<endl;
            continue;
         }

         TGraphErrors* grsig = new TGraphErrors(isigmax-isigmin+1, &x[isigmin], &y[isigmin], &ex[isigmin], &ey[isigmin]);
         grsig->SetNameTitle(Form("grsig_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch), Form("grsig_evt_%d_ch_%d I=%0.1f",oscEvent->evt,oscChannel->ch,integral));
         grsig->SetMarkerStyle(24);
         grsig->SetMarkerColor(8);
         grsig->SetLineColor(8);
         if (gdebug) {
            new TCanvas;
            grsig->Draw("ap");
         }

         //-- find beginning of the edge and beginning of the local maximum

         TMarker* marker_edge_halfmax = 0;
         TMarker* marker_edge_max = 0;
         TMarker* marker_edge_min = 0;

         // find start of the nearest maximum

         //-- find a channel of the half-max
         // Double_t edge_halfmax = ymax/2;           // no clipping C
         Double_t edge_halfmax = (2./3.)*ymax;           // no clipping C
         //-- Double_t edge_halfmax = 0.3*ymax;    // clipping C
         Int_t iedge_halfmax = 0;
         //------------------------------------for (int i=0; i<np; ++i) {
         for (int i=isigmin; i<=isigmax; ++i) {
            if (y[i] >= edge_halfmax) {
               iedge_halfmax = i;
               break;
            }
         }
         cout<< "ymax = " << ymax << " iedge_halfmax = " << iedge_halfmax << " x[iedge_halfmax] = " << x[iedge_halfmax] << " y[iedge_halfmax] = " << y[iedge_halfmax] <<endl;
         if (gdebug) {
            marker_edge_halfmax = new TMarker(x[iedge_halfmax],y[iedge_halfmax],20);
            marker_edge_halfmax->SetMarkerColor(1);
            marker_edge_halfmax->Draw("same");
         }

         // find beginning of the maximum
         Int_t edge_max_check = 5;
         Int_t iedge_max = 0;
         for (int i=iedge_halfmax-1; i<1024; ++i) {
            // Double_t ymax_curr = y[i] + 1.*ey[i];
            // Double_t ymax_curr = y[i] + 2.*ey[i];
            Double_t ymax_curr = y[i] + 3.*ey[i];
            // Double_t ymax_curr = y[i] + 5.*ey[i];

            // count the number of channels with value less than ymax_curr
            Int_t nless = 0;
            for (int j=i+1; j<1024; ++j) {
               if (y[j] < ymax_curr) ++nless;
               else break;
               if (nless > edge_max_check) break;              // do not need to test more
            }
            if (nless > edge_max_check) {
               iedge_max = i;
               break;
            }
         }
         cout<< "iedge_max = " << iedge_max << " x[iedge_max] = " << x[iedge_max] << " y[iedge_max] = " << y[iedge_max] <<endl;
         if (gdebug) {
            marker_edge_max = new TMarker(x[iedge_max],y[iedge_max],20);
            marker_edge_max->SetMarkerColor(2);
            marker_edge_max->Draw("same");
         }

         // find lower limit
         Int_t iedge_min = 0;
         Double_t low_limit = 3.*bkg_sigma;
         for (int i=iedge_halfmax-1; i>=0; --i) {
            if (y[i] < low_limit) {
               iedge_min = i;
               break;
            }
         }
         cout<< "iedge_min = " << iedge_min << " x[iedge_min] = " << x[iedge_min] << " y[iedge_min] = " << y[iedge_min] <<endl;
         if (iedge_max - iedge_min < 4) {
            iedge_max = iedge_min + 4;
            cout<< "--> forced iedge_max = iedge_min + 4 = " << iedge_max <<endl;
         }
         if (gdebug) {
            marker_edge_min = new TMarker(x[iedge_min],y[iedge_min],20);
            marker_edge_min->SetMarkerColor(4);
            marker_edge_min->Draw("same");
            // // redraw edge_halfmax marker
            // marker_edge_halfmax->Draw("same");
         }

         // put into tree the beginning of the local maximum
         oscFit->ymaxi[ich] = y[iedge_max];
         oscFit->xmaxi[ich] = x[iedge_max];

         if (x[iedge_min] < sigmin || x[iedge_max] > sigmax)   // FIXME
         {
            // the lower channel or upper channel is out of range (sigmin,sigmax): fill the tree and goto next DRS4 channel
            tree->Fill();
            continue;
         }

         // fill some data into the tree
         oscFit->v[ich] = ymax;
         oscFit->bkg[ich] = bkg_mean;
         oscFit->sbkg[ich] = bkg_sigma;

         //---------------------- line fit begin ----------------------------

         fun_line->SetParameters(parline);

         Double_t epsline = 1e-7;
         Int_t ifit_first_line = iedge_min+1;                        // do not include min
         Int_t ifit_last_line = iedge_max-1;                         // do not include max
         Int_t np_edge_line = ifit_last_line - ifit_first_line + 1;  // the number of points on the edge
         oscFit->npline[ich] = np_edge_line;                         // store in the tree
         Double_t fit_edge_xmin = x[ifit_first_line] - epsline;
         Double_t fit_edge_xmax = x[ifit_last_line] + epsline;
         Int_t npline = ifit_last_line - ifit_first_line + 1;
         oscFit->npline[ich] = npline;
         cout<< "line fit from fit_edge_xmin " << fit_edge_xmin << " to fit_edge_xmax " << fit_edge_xmax <<endl;
         if (ifit_last_line - ifit_first_line + 1 >= nparline)
         {
            grsig->Fit(fun_line, "", "", fit_edge_xmin, fit_edge_xmax);

            // results
            Double_t res_intercept = grsig->GetFunction(fun_line->GetName())->GetParameter(0);
            Double_t res_slope = grsig->GetFunction(fun_line->GetName())->GetParameter(1);
            Double_t intersection = 0;
            Double_t chi2line = grsig->GetFunction(fun_line->GetName())->GetChisquare();
            Double_t ndfline = grsig->GetFunction(fun_line->GetName())->GetNDF();
            Double_t chi2lineNDF = ndfline > 0? chi2line/ndfline: 0;
            if (TMath::Abs(res_slope) > epsline) intersection = (-1)*res_intercept/res_slope;
            cout<< "intersection = " << intersection << " chi2line = " << chi2line << " ndfline = " << ndfline <<endl;

            // fill the tree variables
            oscFit->p0[ich] = res_intercept;
            oscFit->p1[ich] = res_slope;
            oscFit->chi2line[ich] = chi2lineNDF;
            oscFit->dline[ich] = intersection;
         }

         //---------------------- line fit end --------------------------

         //---------------------- pulse fit begin -------------------------

         Int_t ipulse_x1 = isigmin;
         // Int_t ipulse_x2 = iedge_max+1;       // include point after the start of local max: important for DRS4 v.2
         Int_t ipulse_x2 = iedge_max;
         Double_t pulse_x1 = x[ipulse_x1] - epsline;      // include this point
         Double_t pulse_x2 = x[ipulse_x2] + epsline;      // include this point

         TF1* fitted_function = 0;
         TF1* fitted_function_refit = 0;

         A = 1000;
         x0 = 50;
         // tau1 = 0.050;
         tau1 = 0.750;
         tau2 = 1;      // clipping C
         // tau2 = 80;         // noC
         T = 0;
         sigma = 0.200;
         fun_Psigma->SetParameters(par);

         // fill UserInfo
         if (initParameters == 0) {
            initParameters = new PulseParameters;
            initParameters->SetParameters(par);
         }

         // release all parameters
         for (int i=0; i<fun_Psigma->GetNpar(); ++i) fun_Psigma->ReleaseParameter(i);

         // fix some of parameters
         fun_Psigma->FixParameter(ipar_tau1, fun_Psigma->GetParameter(ipar_tau1));
         // fun_Psigma->FixParameter(ipar_tau2, fun_Psigma->GetParameter(ipar_tau2));
         // fun_Psigma->FixParameter(ipar_sigma, fun_Psigma->GetParameter(ipar_sigma));
         fun_Psigma->FixParameter(ipar_T, fun_Psigma->GetParameter(ipar_T));
         fun_Psigma->FixParameter(ipar_sigma,fun_Psigma->GetParameter(ipar_sigma));

         // fun_Psigma->SetParLimits(ipar_tau2, 1., 20.);
         // fun_Psigma->SetParLimits(ipar_T, 1., 60.);

         Int_t np = isigmax - isigmin + 1;
         oscFit->np[ich] = np;

         if (np >= npar)
         {
            if (gdebug) grsig->Fit(fun_Psigma, "+R", "", sigmin,sigmax);
            else grsig->Fit(fun_Psigma,"+R0", "goff", sigmin,sigmax);

            fitted_function = grsig->GetFunction(fun_Psigma->GetName()); 

            //------------- refit -------------

            cout<< "-- refit" <<endl;
            
            fun_Psigma_refit->SetParameters(grsig->GetFunction(fun_Psigma->GetName())->GetParameters());

            // release all parameters
            for (int i=0; i<fun_Psigma->GetNpar(); ++i) fun_Psigma->ReleaseParameter(i);

            // fix T for refit
            // fun_Psigma_refit->FixParameter(ipar_tau1, grsig->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_tau1));
            fun_Psigma_refit->FixParameter(ipar_tau2, grsig->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_tau2));
            fun_Psigma_refit->FixParameter(ipar_T, grsig->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_T));
            fun_Psigma_refit->FixParameter(ipar_sigma, grsig->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_sigma));

            // fun_Psigma->FixParameter(ipar_tau2, grsig->GetFunction(fun_Psigma->GetName())->GetParameter(3));
            // //fun_Psigma->SetParameter(ipar_tau1, 1.);
            // fun_Psigma->ReleaseParameter(ipar_tau1);
            // fun_Psigma->SetParLimits(ipar_tau1, 0., 10.);

            if (gdebug) grsig->Fit(fun_Psigma_refit, "+R", "", pulse_x1,pulse_x2);
            else grsig->Fit(fun_Psigma_refit, "+R0", "goff", pulse_x1,pulse_x2);
            //cout<< "--> refit is done" <<endl;

            fitted_function_refit = grsig->GetFunction(fun_Psigma_refit->GetName());
            for (int i=0; i<fitted_function_refit->GetNpar(); ++i) cout<< i << "\t " << fitted_function_refit->GetParName(i) << "\t " << fitted_function_refit->GetParameter(i) <<endl;

            // swap functions such that the stats box shows the last fitted function
            TF1* frefit = (TF1*) grsig->GetListOfFunctions()->Last();
            grsig->GetListOfFunctions()->Remove(frefit);
            grsig->GetListOfFunctions()->AddFirst(frefit);

            //-- put line fit parameters first
            // TF1* flinefit = (TF1*) grsig->GetListOfFunctions()->FindObject(fun_line->GetName());
            // if (flinefit) {
            //    grsig->GetListOfFunctions()->Remove(flinefit);
            //    grsig->GetListOfFunctions()->AddFirst(flinefit);
            // }

            // results
            if (fitted_function_refit)
            {
               Double_t accuracy_x = 1e-3;   // optional parameter
               Double_t middle_point_y = 0.5*(y[iedge_min]+y[iedge_max]);
               Double_t middle_point_x = fitted_function_refit->GetX(middle_point_y, x[iedge_min],x[iedge_max],accuracy_x);
               cout<< "--> middle_point_y = " << middle_point_y << " middle_point_x = " << middle_point_x <<endl;
               //cout<< "--> before Eval" <<endl;
               Double_t fval = fitted_function_refit->Eval(middle_point_x);
               //cout<< "--> after Eval" <<endl;
               oscFit->xt[ich] = middle_point_x;
               oscFit->yt[ich] = fval;
               if (gdebug) {
                  TMarker* marker_tang = new TMarker(middle_point_x,fval,3);
                  marker_tang->SetMarkerColor(7);
                  marker_tang->SetMarkerSize(2);
                  marker_tang->Draw("same");
               }
               Double_t dval = fitted_function_refit->Derivative(middle_point_x);
               oscFit->dydx[ich] = dval;
               Double_t dintersect = 0;
               const Double_t eps = 1e-12;
               if (TMath::Abs(dval) > eps) dintersect = middle_point_x - fval/dval;
               oscFit->d[ich] = dintersect;
               cout<< "--> dintersect = " << dintersect << " middle_point_x = " << middle_point_x << " fval = " << fval << " dval = " << dval <<endl;

               oscFit->t[ich] = fitted_function_refit->GetParameter("x0");
               Double_t chi2 = fitted_function_refit->GetChisquare();
               Double_t ndf = fitted_function_refit->GetNDF();
               Double_t chi2ndf = ndf > 0? chi2 / ndf: 0;
               oscFit->chi2[ich] = chi2ndf;
               oscFit->tau1[ich] = fitted_function_refit->GetParameter("tau1");
               oscFit->tau2[ich] = fitted_function_refit->GetParameter("tau2");
               oscFit->T[ich] = fitted_function_refit->GetParameter("T");
               oscFit->sigma[ich] = fitted_function_refit->GetParameter("sigma");

               grsig->SetTitle(Form("%s d=%0.1f #chi^{2}=%0.1f",grsig->GetTitle(), dintersect, chi2ndf));
            }
         }
         else {
            cout<< "--- no fit for this channel: isigmax - isigmin + 1 = " << isigmax - isigmin + 1 <<endl;
         }

         //---------------------- pulse fit end -------------------------

         // cleanup
         if (!gdebug) {
            if (grsig) delete grsig;
            //-- if (grsig_refit) delete grsig_refit;
         }
      }  // loop over ich
      otree->Fill();
   }

   otree->GetUserInfo()->Add(initParameters);
   cout<< "otree->GetUserInfo()->GetEntries() = " << otree->GetUserInfo()->GetEntries() <<endl;

   // for (int ich=0; ich<4; ++ich) {
   //    //if (ich < 2) continue;
   //    new TCanvas;
   //    otree->Draw(Form("adc[%d]",ich), Form("adc[%d]>0&&adc[%d]<1000",ich,ich));
   //    new TCanvas;
   //    otree->Draw(Form("chi2[%d]",ich), Form("chi2[%d]>0&&chi2[%d]<10",ich,ich));
   // }

   // turn off debug
   debug = false;
   gdebug = false;

   return otree;
}

TTree* pulse(const char *ifname="Co60_STM_LSO2x2_NOCAP_NOAMP_split.xml.root"
      , Double_t bkgmin=0, Double_t bkgmax=80, Double_t sigmin=80, Double_t sigmax=512
      , Int_t entry_first=0, Int_t entry_last=-1
      , bool setdebug=false, bool setgdebug=false
      , Int_t ignorech1=-1
      , Int_t ignorech2=-1
      , Int_t ignorech3=-1
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

   TTree* otree = pulse(tree, bkgmin,bkgmax, sigmin,sigmax, entry_first,entry_last, setdebug,setgdebug, ignorech1,ignorech2,ignorech3);
   
   ofile->Write();
   return otree;
}

