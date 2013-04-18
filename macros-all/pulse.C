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
TF1* intfun(TF1* fun);

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
   // flag
   Bool_t ok[8];
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
         //
         ok[i] = true;
      }
   }
   OscFit(): TNamed() {
      SetNameTitle("oscFit", "oscFit");
      clear();
   }
   ClassDef(OscFit,15);
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

class Fitter {
protected:
   TF1* tf1_;                             //!
   std::map<std::string, int> index_;
   std::vector<std::string> fix_;
   Double_t par_[10];
   Double_t par_init_[10];
   Double_t epar_[10];
   Int_t npar_;
   Double_t xmin_;
   Double_t xmax_;
   // fit results
   Double_t chi2_;
   Int_t color_;
   Int_t lwidth_;
   Int_t Npx_;
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
             , color_(1)
             , lwidth_(1)
             // , Npx_(1024)
             , Npx_(10000)
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
   //Fitter(const Fitter& fitter): tf1_(fitter.tf1_)
   Fitter(const Fitter& fitter): tf1_(new TF1(*fitter.tf1_))
                                 , npar_(fitter.npar_)
                                 , xmin_(fitter.xmin_)
                                 , xmax_(fitter.xmax_)
                                 , color_(fitter.color_)
                                 , lwidth_(fitter.lwidth_)
                                 , Npx_(fitter.Npx_)
   {
      cout<< "Fitter::Fitter(const Fitter&)" <<endl;
      tf1_ = new TF1(*fitter.tf1_);
      for (int i=0; i<fitter.npar_; ++i) {
         par_[i] = fitter.par_[i];
         epar_[i] = fitter.epar_[i];
         par_init_[i] = fitter.par_init_[i];
         //par_[i] = fitter.tf1_->GetParameter(i);
         //epar_[i] = fitter.tf1_->GetParError(i);
      }
      tf1_->SetParameters(par_);
      tf1_->SetNpx(Npx_);
      tf1_->SetLineColor(color_);
      tf1_->SetLineWidth(lwidth_);
   }
   virtual ~Fitter() {
      if (tf1_) delete tf1_; tf1_=0;
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
   void SetLineColor(Int_t color) {
      color_ = color;
      if (tf1_) tf1_->SetLineColor(color_);
   }
   void SetLineWidth(Int_t lwidth) {
      lwidth_ = lwidth;
      if (tf1_) tf1_->SetLineWidth(lwidth_);
   }
   void SetNpx(Int_t Npx) {
      Npx_ = Npx;
      if (tf1_) tf1_->SetNpx(Npx_);
   }
   Int_t GetLineColor() const {
      cout<< "GetLineColor: color_ = " << color_ << " tf1_->GetLineColor() = " << tf1_->GetLineColor() <<endl;
      return color_;
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
      cout<< "-- Fitter::Fit: "; GetLineColor();
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
   void SetInit(const std::string& parname, Double_t parvalue, Bool_t fix=0) {
      if (!Exist(parname)) {
         cout<< "***Error Fitter::Fix: unknown parameter " << parname <<endl;
         exit(0);
      }
      par_init_[Index(parname)] = parvalue;
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
      //cout<< "Fitter::Release" <<endl;
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
   void PrintPars() const {for (int ipar=0; ipar<npar_; ++ipar) cout<< par_[ipar] << " ";}
   void PrintParsInit() const {for (int ipar=0; ipar<npar_; ++ipar) cout<< par_init_[ipar] << " ";}
   void SaveInit() {
      for (int ipar=0; ipar<npar_; ++ipar) par_init_[ipar] = par_[ipar];
   }
   virtual void PrepareFit(Fitter* fitter=0)
   {
      //cout<< "Fitter::PrepareFit: color_ = " << color_ <<endl;
      if (fitter)
      {
         // use parameters from the fitter
         TF1* fun = 0;     // just pointer
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
      else {
         // set initilal parameters
         for (int ipar=0; ipar<npar_; ++ipar) par_[ipar] = par_init_[ipar];
      }

      for (std::map<std::string, int>::const_iterator it=index_.begin(); it!=index_.end(); ++it) {
         tf1_->SetParName(it->second, it->first.c_str());
      }
      tf1_->SetParameters(par_);
      tf1_->SetNpx(Npx_);
      tf1_->SetLineColor(color_);
      tf1_->SetLineWidth(lwidth_);
      for (int ipar=0; ipar<npar_; ++ipar) tf1_->ReleaseParameter(ipar);
      for (unsigned i=0; i<fix_.size(); ++i) {
         tf1_->FixParameter(Index(fix_[i]), tf1_->GetParameter(Index(fix_[i])));
      }
      Complete();
      //cout<< "at the end of Fitter::PrepareFit: color_ = " << color_ << " tf1_->GetLineColor() = " << tf1_->GetLineColor() <<endl;
      //cout<< "end of Fitter::PrepareFit" <<endl;
   }
   virtual Bool_t Complete() const        // NB: virtual
   {
      //cout<< "Fitter::Complete" <<endl;
      if (!tf1_) cout<< "***Error Fitter::Complete: tf1_ == 0" <<endl<<exitl;
      if (xmin_ >= xmax_) cout<< "***Error Fitter::Complete: xmin_ >= xmax_" <<endl<<exitl;
      if (index_.size()==0) cout<< "*Error Fitter::Complete: index_.size()==0" <<endl<<exitl;
      return kTRUE;
   }
   Double_t GetX(Double_t y, Double_t xmin, Double_t xmax, Double_t accuracy) const {
      return tf1_->GetX(y, xmin, xmax, accuracy);
   }
   Double_t GetMaximum(Double_t xmin, Double_t xmax, Double_t epsilon=1.E-10, Int_t maxiter=100) const {
      // return tf1_->GetMaximum(xmin, xmax, epsilon, maxiter, logx);
      return tf1_->GetMaximum(xmin, xmax, epsilon, maxiter);
   }
   Double_t GetMaximumX(Double_t xmin=0, Double_t xmax=0, Double_t epsilon=1.E-10, Int_t maxiter=100) const {
      // return tf1_->GetMaximumX(xmin, xmax, epsilon, maxiter, logx);
      return tf1_->GetMaximumX(xmin, xmax, epsilon, maxiter);
   }
   Double_t Eval(Double_t x) const {
      return tf1_->Eval(x);
   }
   Double_t Derivative(Double_t x, Double_t* params=0, Double_t epsilon=0.001) const {
      return tf1_->Derivative(x, params, epsilon);
   }
   Double_t GetParameter(const char* name) const {
      return tf1_->GetParameter(name);
   }
   Double_t GetChisquare() const {
      return tf1_->GetChisquare();
   }
   Double_t GetNDF() const {
      return tf1_->GetNDF();
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
      tf1_->SetNpx(1024);
   }
   FitPulse(const FitPulse& fitter): Fitter(fitter)
   {
      cout<< "FitPulse(const FitPulse&)" <<endl;
      npar_ = fitter.npar_;
      *this << "A" << "x0" << "tau1" << "tau2" << "T" << "sigma";
      cout<< "FitPulse::FitPulse: index_.size() = " << index_.size() <<endl;
      //tf1_ = new TF1(*fitter.tf1_);
      ////tf1_->SetNpx(10240);

      index_ = fitter.index_;
      fix_ = fitter.fix_;
      for (int i=0; i<npar_; ++i) {
         par_[i] = fitter.par_[i];
         epar_[i] = fitter.epar_[i];
         par_init_[i] = fitter.par_init_[i];
      }
      tf1_->SetParameters(par_);
      tf1_->SetNpx(Npx_);
      tf1_->SetLineColor(color_);
      tf1_->SetLineWidth(lwidth_);
   }
   FitPulse* Clone() const {
      cout<< "FitPulse::Clone" <<endl;
      FitPulse* fitPulse = new FitPulse(*this);

      cout<< "FitPulse::Clone: my PrintPars --> ";
      fitPulse->PrintPars();
      cout<<endl;
      cout<< "FitPulse::Clone: cloned PrintPars: ";
      fitPulse->PrintPars();
      cout<<endl;
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
   FitBoard(const FitBoard& fit): nchan_(fit.nchan_)
                                  , defaultComplete_(fit.defaultComplete_)
                                  , defaultClone_(fit.defaultClone_)
                                  , defaultLock_(fit.defaultLock_)
   {
      for (unsigned i=0; i<nchan_+1; ++i) {
         fit_[i] = fit.fit_[i]->Clone();
      }
      //cout<< "FitBoard copy constructor: after Clone" <<endl;
   }
   Bool_t CompleteDefault() const {
      return fit_[0]->Complete();
   }
   void PrepareFit() {
      cout<< "FitBoard::PrepareFit()" <<endl;
      if (!defaultClone_) cout<< "FitBoard::PrepareFit: default did not applied yet" <<endl<<exitl;
   }
   void PrepareFit(int chan, Fitter* fitter=0) {
      cout<< ".. begin of FitBoard::PrepareFit: chan = " << chan <<endl;
      if (!defaultClone_) cout<< "FitBoard::PrepareFit: default did not applied yet" <<endl<<exitl;
      fit_[chan]->PrepareFit(fitter);
      //cout<< ".. FitBoard::PrepareFit: fit_[chan]->GetLineColor: "; fit_[chan]->GetLineColor();
      //cout<< ".. end of FitBoard::PrepareFit" <<endl;
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
   void SetInit(Int_t chan, const std::string& parname, Double_t parvalue, Bool_t fix=0) {
      assert(fit_[chan] != 0);
      // if (!fit_[chan]) cout<< "pointer to channel " << chan << " == 0" <<endl<<exitl;
      if (!CompleteDefault() && chan != 0) cout<< "***Warning: you must complete default fitter first" <<endl<<exitl;
      fit_[chan]->SetInit(parname, parvalue, fix);
   }
   void DefaultClone() {
      cout<< "DefaultClone" <<endl;
      //////////////////////////fit_[0]->PrepareFit();
      for (unsigned i=1; i<nchan_+1; ++i) fit_[i] = fit_[0]->Clone();
      defaultClone_ = kTRUE;
      defaultLock_ = kTRUE;
   }
   void PrintPars() const {
      for (unsigned ich=0; ich<nchan_+1; ++ich) {
         cout<< "channel " << ich << " par_: ";
         fit_[ich]->PrintPars();
         cout<<endl;
      }
   }
   void PrintParsInit() const {
      for (unsigned ich=0; ich<nchan_+1; ++ich) {
         cout<< "channel " << ich << " par_init_: ";
         fit_[ich]->PrintParsInit();
         cout<<endl;
      }
   }
   void SaveInit() {
      for (unsigned i=1; i<nchan_+1; ++i) fit_[i]->SaveInit();
   }
   void SetLineColor(Int_t chan, Int_t color) {
      assert(fit_[chan] != 0);
      if (!CompleteDefault() && chan != 0) cout<< "***Warning: you must complete default fitter first" <<endl<<exitl;
      fit_[chan]->SetLineColor(color);
   }
   void SetLineWidth(Int_t chan, Int_t lwidth) {
      assert(fit_[chan] != 0);
      if (!CompleteDefault() && chan != 0) cout<< "***Warning: you must complete default fitter first" <<endl<<exitl;
      fit_[chan]->SetLineWidth(lwidth);
   }
   void SetNpx(Int_t chan, Int_t Npx) {
      assert(fit_[chan] != 0);
      if (!CompleteDefault() && chan != 0) cout<< "***Warning: you must complete default fitter first" <<endl<<exitl;
      fit_[chan]->SetNpx(Npx);
   }
   void Fix(const std::string& parname) {
      bool found = false;
      for (unsigned i=0; i<nchan_+1; ++i) {
         if (fit_[i]->Exist(parname)) {
            found = true;
            fit_[i]->Fix(parname);
         }
      }
      if (!found) {
         cout<< "***Error FixBoard::Fix: no such parname: " << parname <<endl<<exitl;
      }
   }
   void Release(const std::string& parname) {
      bool found = false;
      for (unsigned i=0; i<nchan_+1; ++i) {
         if (fit_[i]->Exist(parname)) {
            found = true;
            fit_[i]->Release(parname);
         }
      }
      if (!found) {
         cout<< "***Error FixBoard::Release: no such parname: " << parname <<endl<<exitl;
      }
   }
};

void beamfit()
{
   const char ifname[] = "mppc16_18.root.osc.root";
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) cout<< "File not found: " << ifname <<endl<<exitl;
   TTree* t = (TTree*) ifile->Get("t");
   if (!t) cout<< "No tree t was found in file " << ifname <<endl<<exitl;

   Int_t evt = 0;
   Int_t channel = 1;

   OscEvent* oscEvent = 0;
   t->SetBranchAddress("oscEvent",&oscEvent);

   t->GetEntry(evt);

   const OscChannel* oscChannel = (OscChannel*) oscEvent->oscChannels->At(channel-1);
   cout<< "oscChannel->ch = " << oscChannel->ch <<endl;

   const Double_t ythres = 10;
   Double_t xthres = 0;

   // find approx value of the bkg
   Int_t nbkg = 20;
   Double_t ybkg_approx = 0;
   for (int i=0; i<nbkg; ++i) ybkg_approx += -1.*oscChannel->y[i];
   ybkg_approx /= nbkg;
   cout<< "ybkg_approx = " << ybkg_approx <<endl;

   Int_t ithres = 0;
   for (int i=0; i<1024; ++i) {
      Double_t y = -1.*oscChannel->y[i] - ybkg_approx;
      if (y > ythres) {
         xthres = oscChannel->x[i];
         ithres = i;
         break;
      }
   }
   if (ithres == 0) cout<< "no pulse found" <<endl<<exitl;
   else cout<< "ithres = " << ithres << " xthres = " << xthres <<endl;

   // define as a start of the signal noffset data points
   const Int_t ioffset = 20;
   Int_t isigmin = ithres - ioffset;

   // find bkg level

   // histogram for background fit
   TH1F *hbkg = new TH1F("hbkg", "hbkg", 400, -100, 100);

   Double_t ybkg[1024];
   Int_t np_bkg = 0;
   for (int i=0; i<isigmin; ++i) {
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
   cout<< "bkg_mean = " << bkg_mean << " bkg_sigma = " << bkg_sigma <<endl;

   // find the maximum

   Double_t sigmin = oscChannel->x[isigmin];
   cout<< "isigmin = " << isigmin << " sigmin = " << sigmin <<endl;

   Double_t sigmax = 0;
   Int_t isigmax = 0;

   Double_t ymaximum = 0;
   Double_t xmaximum = 0;
   Int_t imaximum = 0;
   Double_t yhalfmaximum = 0;
   Double_t xhalfmaximum = 0;
   Int_t ihalfmaximum = 0;
   for (int i=ithres; i<1024; ++i) {
      Double_t y = -1.*oscChannel->y[i] - bkg_mean;
      //cout<< i << " y = " << y <<endl;
      if (y > ymaximum) {
         ymaximum = y;
         imaximum = i;
         xmaximum = oscChannel->x[i];
         //cout<< "i = " << i << " y = " << y << " ymaximum/8. = " << ymaximum/8. <<endl;
      }
      if (ihalfmaximum == 0 && y < ymaximum / 2.) {
         ihalfmaximum = i;
         yhalfmaximum = y;
         xhalfmaximum = oscChannel->x[i];
      }
      if (y < ymaximum / 8.) {
         //cout<< "find isigmax: i = " << i << " y = " << y << " ymaximum/8. = " << ymaximum/8. <<endl;
         isigmax = i;
         sigmax = oscChannel->x[i];
         break;
      }
   }
   cout<< "imaximum = " << imaximum << " ymaximum = " << ymaximum <<endl;
   cout<< "ihalfmaximum = " << ihalfmaximum << " xhalfmaximum = " << xhalfmaximum << " yhalfmaximum = " << yhalfmaximum <<endl;

   cout<< "isigmin " << isigmin << " isigmax " << isigmax << " sigmin " << sigmin << " sigmax " << sigmax <<endl;
   
   Double_t x[1024], y[1024], ey[1024];
   Int_t np = 0;

   for (int i=isigmin; i<=isigmax; ++i) {
      x[np] = oscChannel->x[i];
      y[np] = -1.*oscChannel->y[i] - bkg_mean;
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

   cout<<endl<< "-- FitBoard test" <<endl;

   Int_t nchan = 4;
   //
   // define fit
   //
   FitBoard fit(nchan);
   cout<< "before SetLineWidth" <<endl;

   cout<< "define default" <<endl;
   channel = 0;
   fit.fit_[channel] = new FitPulse(sigmin,sigmax);
   fit.SetLineColor(channel, 1);
   fit.SetLineWidth(channel, 1);
   fit.SetNpx(channel, 1024);
   fit.Set(channel, "A",       100.);
   fit.Set(channel, "x0",      xthres);
   fit.Set(channel, "tau1",    0.,   fix);
   fit.Set(channel, "tau2",    2.);
   fit.Set(channel, "T",       0.014,   fix);
   // fit.Set(channel, "T",       40.,   fix);
   fit.Set(channel, "sigma",   0.300);
   fit.DefaultClone();

   // some specific requirements for the channel 2
   channel = 2;
   fit.Set(channel, "sigma", 0.300, nofix);
   // fit.Set(channel, "sigma", 0., fix);
   // fit[channel]->Release("tau1");
   // fit[channel]->Release("T");

   cout<< "refit" <<endl;
   //
   // define refit
   //
   FitBoard refit(fit);
   for (int i=0; i<nchan+1; ++i) {
      refit.SetLineWidth(i,2);
      refit.SetLineColor(i,2);
   }
   refit.Release("sigma");
   // refit.Release("tau1");
   refit.Fix("tau1");
   refit.Fix("tau2");

   // This is suppose to be a loop over channels

   int ich = 1;

   fit.PrepareFit(ich);
   fit[ich]->Fit(gr, "", "", sigmin, xhalfmaximum);

   cout<< "\nrefit\n" <<endl;

   Double_t xmin_refit = sigmin;
   Double_t xmax_refit = xmaximum;

   refit.PrepareFit(ich,fit[ich]);
   refit[ich]->Fit(gr, "+", "", xmin_refit, xmax_refit);
}

/////////////////////////////////////////////////////////////////////////////////

bool debug = false;
bool gdebug = false;

TTree* pulse(TTree *tree
      //, Double_t bkgmin=0, Double_t bkgmax=40, Double_t sigmin=40, Double_t sigmax=120
      , Int_t entry_first=0, Int_t entry_last=-1
      , bool setdebug=false, bool setgdebug=false
      , Int_t ignorech1=-1
      , Int_t ignorech2=-1
      , Int_t ignorech3=-1
      , Int_t ignorech4=-1
      , Int_t ignorech5=-1
      , Int_t ignorech6=-1
      , Int_t ignorech7=-1
      , Int_t ignorech8=-1
      )
{
   if (setdebug) debug = true;
   if (setgdebug) gdebug = true;

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

   // UserInfo: will be filled before the first fit
   //-- PulseParameters* initParameters = 0;
   PulseParameters* initParameters = new PulseParameters;

   // the number of channels in the data
   tree->GetEntry(entry_first);
   Int_t nchannels = oscEvent->oscChannels->GetEntries();

   // hisstogram for kkg
   TH1F *hbkg = new TH1F("hbkg", "hbkg", 400, -100, 100);

   /////////////////////////////////////////////////////////////////////////

   Bool_t fix = kTRUE;
   Bool_t nofix = kFALSE;
   nofix = !fix;           // to avoid warning: unused variable ‘nofix’

   Int_t nchan = 8;
   const Double_t fun_xmin = 0;
   const Double_t fun_xmax = 200.;
   //
   // define fit
   //
   Int_t channel = 0;

   FitBoard fit(nchan);

   channel = 0;
   fit.fit_[channel] = new FitPulse(fun_xmin,fun_xmax);
   fit.SetLineColor(channel, 1);
   fit.SetLineWidth(channel, 1);
   fit.SetNpx(channel, 1024);
   fit.Set(channel, "A",       100.);
   fit.Set(channel, "x0",      0.5*(fun_xmin+fun_xmax));
   fit.Set(channel, "tau1",    0.,   fix);
   fit.Set(channel, "tau2",    2.);
   fit.Set(channel, "T",       0.014,   fix);
   // fit.Set(channel, "T",       40.,   fix);
   fit.Set(channel, "sigma",   0.300);
   cout<< "default PrintPars: ";
   fit.fit_[channel]->PrintPars();
   cout<<endl;
   cout<< "\ncall fit.DefaultClone" <<endl;
   fit.DefaultClone();
   cout<< "\ncall fit.PrintPars" <<endl;
   fit.PrintPars();

   //return 0;

   // some specific requirements for the channel 2
   channel = 2;
   fit.Set(channel, "sigma", 0.300, nofix);
   // fit.Set(channel, "sigma", 0., fix);
   // fit[channel]->Release("tau1");
   // fit[channel]->Release("T");

   cout<< "call fit.SaveInit" <<endl;
   fit.SaveInit();
   cout<< "call fit.PrintPars" <<endl;
   fit.PrintPars();
   cout<< "call fit.PrintParsInit" <<endl;
   fit.PrintParsInit();

   // //
   // // define refit
   // //
   // cout<< "--- define refit" <<endl;
   // FitBoard refit(fit);
   // for (int i=0; i<nchan+1; ++i) {
   //    refit.SetLineWidth(i,2);
   //    refit.SetLineColor(i,2);
   // }
   // refit.Release("sigma");
   // // refit.Release("tau1");
   // refit.Fix("tau1");
   // refit.Fix("tau2");

   // refit.SaveInit();

   /////////////////////////////////////////////////////////////////////////

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

      // this arrays will be used to set event OK flag before the otree->Fill()
      Double_t fit_xmin[8];
      Double_t fit_xmax[8];

      cout<< ".. loop over osc channels. oscEvent->oscChannels = " << oscEvent->oscChannels->GetEntries() <<endl;
      for (int ich=0; ich<nchannels; ++ich)
      {
         channel = ich + 1;
         const OscChannel* oscChannel = (OscChannel*) oscEvent->oscChannels->At(ich);
         cout<< "oscChannel->ch = " << oscChannel->ch <<endl;

         if (oscChannel->ch == ignorech1) continue;
         if (oscChannel->ch == ignorech2) continue;
         if (oscChannel->ch == ignorech3) continue;
         if (oscChannel->ch == ignorech4) continue;
         if (oscChannel->ch == ignorech5) continue;
         if (oscChannel->ch == ignorech6) continue;
         if (oscChannel->ch == ignorech7) continue;
         if (oscChannel->ch == ignorech8) continue;

         const Double_t ythres = 10;
         Double_t xthres = 0;

         // find approx value of the bkg
         Int_t nbkg = 20;
         Double_t ybkg_approx = 0;
         for (int i=0; i<nbkg; ++i) ybkg_approx += -1.*oscChannel->y[i];
         ybkg_approx /= nbkg;
         cout<< "ybkg_approx = " << ybkg_approx <<endl;

         Int_t ithres = 0;
         for (int i=0; i<1024; ++i) {
            Double_t y = -1.*oscChannel->y[i] - ybkg_approx;
            if (y > ythres) {
               xthres = oscChannel->x[i];
               ithres = i;
               break;
            }
         }
         if (ithres == 0) {
            cout<< "no pulse found" <<endl;
            continue;
         }
         else cout<< "ithres = " << ithres << " xthres = " << xthres <<endl;

         // define as a start of the signal noffset data points
         const Int_t ioffset = 20;
         Int_t isigmin = ithres - ioffset;
         if (isigmin < 2) isigmin = 2;

         // find bkg level

         Double_t ybkg[1024];
         Int_t np_bkg = 0;
         for (int i=0; i<isigmin; ++i) {
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
         cout<< "bkg_mean = " << bkg_mean << " bkg_sigma = " << bkg_sigma <<endl;

         // find the maximum

         Double_t sigmin = oscChannel->x[isigmin];
         cout<< "isigmin = " << isigmin << " sigmin = " << sigmin <<endl;

         Double_t sigmax = 0;
         Int_t isigmax = 0;

         Double_t ymaximum = 0;
         Double_t xmaximum = 0;
         Int_t    imaximum = 0;
         Double_t ytail = 0;
         Double_t xtail = 0;
         Int_t    itail = 0;
         for (int i=ithres; i<1024; ++i) {
            Double_t y = -1.*oscChannel->y[i] - bkg_mean;
            //cout<< i << " y = " << y <<endl;
            if (y > ymaximum) {
               ymaximum = y;
               imaximum = i;
               xmaximum = oscChannel->x[i];
               //cout<< "i = " << i << " y = " << y << " ymaximum/8. = " << ymaximum/8. <<endl;
            }
            if (itail == 0 && y < ymaximum / 2.) {
               itail = i;
               ytail = y;
               xtail = oscChannel->x[i];
            }
            if (y < ymaximum / 8.) {
               //cout<< "find isigmax: i = " << i << " y = " << y << " ymaximum/8. = " << ymaximum/8. <<endl;
               isigmax = i;
               sigmax = oscChannel->x[i];
               break;
            }
         }
         cout<< "imaximum = " << imaximum << " ymaximum = " << ymaximum << " xmaximum " << xmaximum <<endl;
         cout<< "itail = " << itail << " xtail = " << xtail << " ytail = " << ytail <<endl;

         cout<< "isigmin " << isigmin << " isigmax " << isigmax << " sigmin " << sigmin << " sigmax " << sigmax <<endl;

         // Double_t ymaximum_test = 0;
         // Double_t xmaximum_test = 0;
         // Int_t    imaximum_test = 0;
         // for (int i=ithres; i<1024; ++i) {
         //    Double_t y = -1.*oscChannel->y[i] - bkg_mean;
         //    //cout<< i << " y = " << y <<endl;
         //    if (y > ymaximum_test) {
         //       ymaximum_test = y;
         //       xmaximum_test = oscChannel->x[i];
         //       imaximum_test = i;
         //    }
         // }
         // cout<< "--- ymaximum_test " << ymaximum_test << " xmaximum_test " << xmaximum_test <<endl;

         Double_t x[1024], y[1024], ey[1024];
         Int_t np = 0;

         Double_t integral = 0;                    // integral from direct sum
         Double_t dx = oscChannel->x[isigmin+1] - oscChannel->x[isigmin];

         Double_t yprev = -1*oscChannel->y[isigmin-1] - bkg_mean;
         for (int i=isigmin; i<=isigmax; ++i)
         {
            x[np] = oscChannel->x[i];
            y[np] = -1.*oscChannel->y[i] - bkg_mean;
            //integral += 0.5*(y[i-1] + y[i]) * dx;
            integral += 0.5*(yprev + y[np]) * dx;
            yprev = y[np];
            ey[np] = bkg_sigma;
            ++np;
         }
         // normalize integral to pC by division to R = 50 Ohm
         integral /= 50.;
         cout<< ".. sum integral = " << integral << " pC" <<endl;
         oscFit->adc[ich] = integral;

         TGraphErrors* grsig = new TGraphErrors(np, x, y, 0, ey);
         grsig->SetNameTitle(Form("grsig_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch), Form("grsig_evt_%d_ch_%d I=%0.1f",oscEvent->evt,oscChannel->ch,integral));
         grsig->SetMarkerStyle(24);
         grsig->SetMarkerColor(8);
         grsig->SetLineColor(8);

         if (gdebug) {
            new TCanvas;
            grsig->Draw("ap");
         }

         // TMarker* marker_edge_halfmax = 0;
         // TMarker* marker_edge_max = 0;
         // TMarker* marker_edge_min = 0;

         if (np >= 6)
         {
            // set initial value for the start time
            fit.Set(channel, "x0", xthres);

            //cout<< "  before fitting: call fit.PrepareFit(channel)" <<endl;
            //cout<< "  par_: "; fit[channel]->PrintPars(); cout<<endl;
            //cout<< "  par_init_: "; fit[channel]->PrintParsInit(); cout<<endl;
            //cout<< "  fitting: call fit.PrepareFit(channel)" <<endl;
            //-- fit.Set(channel, "x0", xthres);
            fit.SetInit(channel, "x0", xthres);    //--FIXIT
            fit.PrepareFit(channel);
            cout<< "  after fit.PrepareFit(channel) par_: "; fit[channel]->PrintPars(); cout<<endl;
            //Double_t fit_xmin = sigmin;
            //Double_t fit_xmax = xtail;
            fit_xmin[ich] = sigmin;
            fit_xmax[ich] = xtail;
            if (gdebug) fit[channel]->Fit(grsig, "+R", "", fit_xmin[ich], fit_xmax[ich]);
            else fit[channel]->Fit(grsig, "+R0", "goff", fit_xmin[ich], fit_xmax[ich]);

            //------------- refit -------------

            /// cout<< "\nrefit\n" <<endl;

            /// Double_t xmin_refit = sigmin;
            /// Double_t xmax_refit = xmaximum;

            /// refit.PrepareFit(channel,fit[channel]);
            /// refit[channel]->Fit(gr, "+", "", xmin_refit, xmax_refit);

            // results
            if (true)
            {
               // Double_t accuracy_x = 1e-3;   // optional parameter
               // Double_t middle_point_y = 0.5*ymaximum;
               // Double_t middle_point_x = fit[channel]->GetX(middle_point_y, xthres,xmaximum,accuracy_x);
               Double_t accuracy_x = 1e-3;   // optional parameter
               Double_t maximum_x = fit[channel]->GetMaximumX(sigmin,xtail);
               Double_t maximum_y = fit[channel]->Eval(maximum_x);
               Double_t middle_point_y = 0.5*maximum_y;

               Double_t middle_point_x = fit[channel]->GetX(middle_point_y, sigmin,maximum_x,accuracy_x);
               cout<< "--> middle_point_y = " << middle_point_y << " middle_point_x = " << middle_point_x <<endl;
               Double_t fval = fit[channel]->Eval(middle_point_x);
               oscFit->xt[ich] = middle_point_x;
               oscFit->yt[ich] = fval;
               if (gdebug) {
                  TMarker* marker_tang = new TMarker(middle_point_x,fval,3);
                  marker_tang->SetMarkerColor(7);
                  marker_tang->SetMarkerSize(2);
                  marker_tang->Draw("same");
               }
               Double_t dval = fit[channel]->Derivative(middle_point_x);
               oscFit->dydx[ich] = dval;
               Double_t dintersect = 0;
               const Double_t eps = 1e-12;
               if (TMath::Abs(dval) > eps) dintersect = middle_point_x - fval/dval;
               oscFit->d[ich] = dintersect;
               cout<< "--> dintersect = " << dintersect << " middle_point_x = " << middle_point_x << " fval = " << fval << " dval = " << dval <<endl;

               oscFit->t[ich] = fit[channel]->GetParameter("x0");
               Double_t chi2 = fit[channel]->GetChisquare();
               Double_t ndf = fit[channel]->GetNDF();
               Double_t chi2ndf = ndf > 0? chi2 / ndf: 0;
               oscFit->chi2[ich] = chi2ndf;
               oscFit->tau1[ich] = fit[channel]->GetParameter("tau1");
               oscFit->tau2[ich] = fit[channel]->GetParameter("tau2");
               oscFit->T[ich] = fit[channel]->GetParameter("T");
               oscFit->sigma[ich] = fit[channel]->GetParameter("sigma");

               if (gdebug) grsig->SetTitle(Form("%s d=%0.1f #chi^{2}=%0.1f",grsig->GetTitle(), dintersect, chi2ndf));
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
      // set channels ok flag
      for (int i=0; i<nchan; ++i)
      {
         if (oscFit->t[i] < fit_xmin[i] || oscFit->t[i] > fit_xmax[i]) oscFit->ok[i] = kFALSE;
         if (oscFit->d[i] < fit_xmin[i] || oscFit->d[i] > fit_xmax[i]) oscFit->ok[i] = kFALSE;
         if (oscFit->sigma[i] < 0 || oscFit->sigma[i] > 10) oscFit->ok[i] = kFALSE;
         if (oscFit->tau2[i] > 1000) oscFit->ok[i] = kFALSE;
         if (oscFit->adc[i] <= 0) oscFit->ok[i] = kFALSE;
      }
      otree->Fill();
   }

   cout<< "add UserInfo" <<endl;
   otree->GetUserInfo()->Add(initParameters);
   cout<< "otree->GetUserInfo()->GetEntries() = " << otree->GetUserInfo()->GetEntries() <<endl;

   cout<< "gDirectory->GetName() = " << gDirectory->GetName() <<endl;

   // turn off debug
   debug = false;
   gdebug = false;

   return otree;
}

TTree* pulse(const char *ifname
      //, Double_t bkgmin=0, Double_t bkgmax=80, Double_t sigmin=80, Double_t sigmax=512
      , Int_t entry_first=0, Int_t entry_last=-1
      , bool setdebug=false, bool setgdebug=false
      , Int_t ignorech1=-1
      , Int_t ignorech2=-1
      , Int_t ignorech3=-1
      , Int_t ignorech4=-1
      , Int_t ignorech5=-1
      , Int_t ignorech6=-1
      , Int_t ignorech7=-1
      , Int_t ignorech8=-1
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

   TTree* otree = pulse(tree, entry_first,entry_last, setdebug,setgdebug, ignorech1,ignorech2,ignorech3,ignorech4,ignorech5,ignorech6,ignorech7,ignorech8);
   
   ofile->Write();
   return otree;
}
