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
   // line fit
   Float_t p0[8];          // slope in y = a*x + b
   Float_t p1[8];          // intercept
   Float_t dline[8];       // intersection of the fitted line with y-axis
   Float_t chi2line[8];    // chi2 of line fit
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
         p0[i] = 0;
         p1[i] = 0;
         dline[i] = 0;
         chi2line[i] = 0;
      }
   }
   OscFit(): TObject() {clear();}
   ClassDef(OscFit,9);
};

#ifdef __MAKECINT__
#pragma link C++ class OscFit;
#endif

ClassImp(OscFit);

//--------------------------------------------------------------------------------

Double_t fline(Double_t xx[], Double_t par[])
{
   Double_t& x = *xx;
   Double_t& intercept = par[0];
   Double_t& slope = par[1];
   return intercept + slope*x;
}

void pulsemax()
{
   const char ifname[] = "UCNa22stm2_70mm_0.root.osc.root";
   // Int_t entry = 9; Int_t channel = 3;
   // Int_t entry = 8; Int_t channel = 3;
   // Int_t entry = 7; Int_t channel = 3;
   // Int_t entry = 6; Int_t channel = 3;       // bad data
   // Int_t entry = 5; Int_t channel = 3;
   // Int_t entry = 3; Int_t channel = 3;
   // Int_t entry = 1; Int_t channel = 3;
   Int_t entry = 0; Int_t channel = 3;
   Double_t bkgmin = 0.;
   Double_t bkgmax = 50.;
   Double_t sigmin = 50.;
   Double_t sigmax = 120.;
   Double_t thres = 10;

   bool gdebug = true;

   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "***ERROR: File not found: " << ifname <<endl;
      exit(1);
   }

   const char tname[] = "t";
   TTree* tree = (TTree*) ifile->Get(tname);
   if (!tree) {
      cout<< "There is no tree \"" << tname << "\" in the file " << ifname <<endl;
      exit(1);
   }

   tree->Draw("y:x", Form("ch==%d&&Entry$==%d",channel,entry), "goff");
   // new TCanvas;
   // tree->Draw("y:x", Form("ch==%d&&Entry$==%d",channel,entry));

   //-------------- fit function definition start -----------------

   Int_t npar = 2;
   TF1* fun_line = new TF1("fun_line", fline, 0, 1000, npar);
   fun_line->SetNpx(1024);
   fun_line->SetName("fun_line");
   std::stringstream sstit;
   sstit << fun_line->GetName();
   //sstit << fun_line->GetName() << " intercept=" << intercept << " slope=" << slope;
   fun_line->SetTitle(sstit.str().c_str());
   fun_line->SetLineColor(1);

   Double_t intercept = 0;
   Double_t slope = 1;

   Double_t par[10];
   fun_line->SetParameters(par);

   enum {ipar_intercept=0, ipar_slope};
   fun_line->SetParName(ipar_intercept, "intercept");
   fun_line->SetParName(ipar_slope, "slope");

   par[ipar_intercept] = intercept;
   par[ipar_slope] = slope;

   //-------------- fit function definition end -----------------

   Double_t x[1024], y[1024], ex[1024], ey[1024];

   // usual names
   Int_t evt = entry;
   Int_t ch = channel;

   for (int i=0; i<tree->GetSelectedRows(); ++i) {
      y[i] = (-1.)*tree->GetV1()[i];
      ey[i] = 0;
      x[i] = tree->GetV2()[i];
      ex[i] = 0;
   }

   // scale factor
   Double_t y1y2 = 0;
   Double_t y2y2 = 0;
   for (int i=0; i<1024; i+=2) {
      Double_t y1 = y[i];
      Double_t y2 = y[i+1];
      y1y2 += y1*y2;
      y2y2 += y2*y2;
   }
   Double_t alpha = y2y2/y1y2;
   cout<< "alpha = " << alpha <<endl;
   for (int i=0; i<1024; i+=2) y[i] *= alpha;

   // find the background
   TH1F* hbkg = new TH1F("hbkg", "hbkg", 800,-400,400);

   Double_t ybkg[1024];
   Int_t np_bkg = 0;
   for (int i=0; i<1024; ++i)
   {
      if (x[i] < bkgmin) continue;
      if (x[i] > bkgmax) break;
      ybkg[np_bkg++] = y[i];
   }

   // sort bkg array to eliminate possible USB spikes: three highest points
   Int_t index[1024];
   TMath::Sort(np_bkg, ybkg, index, kFALSE);

   hbkg->Reset();
   for (int i=0; i<np_bkg-3; ++i) hbkg->Fill(ybkg[index[i]]);
   hbkg->Fit("gaus", "L0Q", "goff");                            // NB: Loglikelihood option
   // new TCanvas;
   // hbkg->Draw();
   // hbkg->Fit("gaus", "LQ", "");                            // NB: Loglikelihood option

   Double_t bkg_mean = hbkg->GetFunction("gaus")->GetParameter("Mean");
   Double_t bkg_sigma = hbkg->GetFunction("gaus")->GetParameter("Sigma");
   cout<< "bkg_mean = " << bkg_mean << " bkg_sigma = " << bkg_sigma <<endl;

   // remove bkg and fill errors
   for (int i=0; i<1024; ++i)
   {
      y[i] -= bkg_mean;
      ey[i] = bkg_sigma;
      ex[i] = 0;
   }

   Int_t nthres = 0;                // the number of channels which exceed threshold
   Float_t ampl_max = 0;

   Double_t x_gr[1024], y_gr[1024], ex_gr[1024], ey_gr[1024];
   Int_t np_gr = 0;
   for (int i=0; i<1024; ++i)
   {
      if (x[i] < sigmin) continue;
      if (x[i] > sigmax) break;

      y_gr[np_gr] = y[i];
      ey_gr[np_gr] = ey[i];
      x_gr[np_gr] = x[i];
      ex_gr[np_gr] = ex[i];

      // maximum in signal window NB: (-1)*y
      if (y_gr[np_gr] > thres) ++nthres;
      if (y_gr[np_gr] > ampl_max) ampl_max = y_gr[np_gr];

      ++np_gr;
   }

   // pointers to graphs
   TGraphErrors *grsig = 0;
   //-- TGraphErrors *grsig_refit = 0;

   grsig = new TGraphErrors(np_gr, x_gr, y_gr, ex_gr, ey_gr);
   grsig->SetNameTitle(Form("grsig_evt_%d_ch_%d",evt,ch), Form("grsig_evt_%d_ch_%d",evt,ch));
   grsig->SetMarkerStyle(24);
   grsig->SetMarkerColor(8);
   grsig->SetLineColor(8);
   if (gdebug) {
      new TCanvas;
      grsig->Draw("ap");
   }

   Double_t integral = 0;        // integral from direct sum
   //-- Double_t integralf = 0;       // integral from fit (will be computed later)

   for (int i=0; i<np_gr-1; ++i) integral += 0.5*(grsig->GetY()[i] + grsig->GetY()[i+1]) * (grsig->GetX()[i+1]-grsig->GetX()[i]);

   // normalize to pC by division to R = 50 Hom
   integral /= 50.;
   cout<< ".. sum integral = " << integral << " pC" <<endl;

   // find start of the nearest maximum

   // use complete arrays

   Double_t ymax = 0;      // initialyze by 0
   for (int i=0; i<1024; ++i) if (y[i] > ymax) ymax = y[i];

   // find a channel of the half-max
   Double_t halfmax = ymax/2;
   Int_t ihalfmax = 0;
   for (int i=0; i<1024; ++i) {
      if (y[i] >= halfmax) {
         ihalfmax = i;
         break;
      }
   }
   cout<< "ihalfmax = " << ihalfmax << " x[ihalfmax] = " << x[ihalfmax] << " y[ihalfmax] = " << y[ihalfmax] <<endl;
   TMarker* marker_halfmax = new TMarker(x[ihalfmax],y[ihalfmax],20);
   marker_halfmax->SetMarkerColor(1);
   marker_halfmax->Draw("same");

   // find beginning of the maximum
   Int_t max_check = 5;
   Int_t imax = 0;
   for (int i=ihalfmax-1; i<1024; ++i) {
      Int_t nless = 0;
      // Double_t ymax_curr = y[i] + 1.*ey[i];
      Double_t ymax_curr = y[i] + 2.*ey[i];
      // Double_t ymax_curr = y[i] + 3.*ey[i];
      //if (x[i] > 65) cout<< "x[i] = " << x[i] << " y[i] = " << y[i] << " ymax_curr = " << ymax_curr <<endl;
      for (int j=i+1; j<1024; ++j) {
         //if (x[i] > 65) cout<< "x[j] = " << x[j] << " y[i] = " << y[i] <<endl;
         if (y[j] < ymax_curr) ++nless;
         else break;
         if (nless > max_check) break;              // do not need to test more
      }
      if (nless > max_check) {
         imax = i;
         break;
      }
   }
   cout<< "imax = " << imax << " x[imax] = " << x[imax] << " y[imax] = " << y[imax] <<endl;
   TMarker* marker_max = new TMarker(x[imax],y[imax],20);
   marker_max->SetMarkerColor(2);
   marker_max->Draw("same");

   // find lower limit
   Int_t imin = 0;
   Double_t low_limit = 3.*bkg_sigma;
   for (int i=ihalfmax-1; i>=0; --i) {
      if (y[i] < low_limit) {
         imin = i;
         break;
      }
   }
   cout<< "imin = " << imin << " x[imin] = " << x[imin] << " y[imin] = " << y[imin] <<endl;
   TMarker* marker_min = new TMarker(x[imin],y[imin],20);
   marker_min->SetMarkerColor(4);
   marker_min->Draw("same");
   // redraw halfmax marker
   marker_halfmax->Draw("same");

   // fit

   par[ipar_intercept] = -1000.;
   par[ipar_slope] = 10.;
   fun_line->SetParameters(par);

   Double_t eps = 1e-7;
   Double_t fit_xmin = x[imin] - eps;                 // include min
   //Double_t fit_xmin = x[imin+1] - eps;                 // do not include min
   Double_t fit_xmax = x[imax-1] + eps;                 // do not include max
   grsig->Fit(fun_line, "", "", fit_xmin, fit_xmax);

   // results
   Double_t res_intercept = grsig->GetFunction(fun_line->GetName())->GetParameter(0);
   Double_t res_slope = grsig->GetFunction(fun_line->GetName())->GetParameter(1);
   Double_t intersection = 0;
   if (TMath::Abs(res_slope) > eps) intersection = (-1)*res_intercept/res_slope;
   cout<< "intersection = " << intersection <<endl;
}

//--
//---------------------------------------------------------------------------
//--

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

/////////////////////////////////////////////////////////////////////////////////

bool debug = false;
bool gdebug = false;

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

         //-- integral over the grsig

         Double_t integral = 0;        // integral from direct sum
         Double_t integralf = 0;       // integral from fit (will be computed later)

         for (int i=0; i<np-1; ++i) integral += 0.5*(grsig->GetY()[i] + grsig->GetY()[i+1]) * (grsig->GetX()[i+1]-grsig->GetX()[i]);

         // normalize to pC by division to R = 50 Hom
         integral /= 50.;
         cout<< ".. sum integral = " << integral << " pC" <<endl;
         oscFit->adc[ich] = integral;

         grsig->SetTitle(Form("%s I=%0.1f",grsig->GetTitle(), integral));

         //-- debug mode selection
         if (ich == 2 && integral < 25) continue;
         if (ich == 3 && integral < 45) continue;

         //---------------------- line fit begin ----------------------------

         // find start of the nearest maximum

         Double_t ymax = 0;      // initialyze by 0
         for (int i=0; i<np; ++i) if (y[i] > ymax) ymax = y[i];

         TMarker* marker_edge_halfmax = 0;
         TMarker* marker_edge_max = 0;
         TMarker* marker_edge_min = 0;

         // find a channel of the half-max
         //-- Double_t edge_halfmax = ymax/2;
         Double_t edge_halfmax = 0.3*ymax;
         Int_t iedge_halfmax = 0;
         for (int i=0; i<np; ++i) {
            if (y[i] >= edge_halfmax) {
               iedge_halfmax = i;
               break;
            }
         }
         cout<< "iedge_halfmax = " << iedge_halfmax << " x[iedge_halfmax] = " << x[iedge_halfmax] << " y[iedge_halfmax] = " << y[iedge_halfmax] <<endl;
         if (gdebug) {
            marker_edge_halfmax = new TMarker(x[iedge_halfmax],y[iedge_halfmax],20);
            marker_edge_halfmax->SetMarkerColor(1);
            marker_edge_halfmax->Draw("same");
         }

         // find beginning of the maximum
         Int_t edge_max_check = 5;
         Int_t iedge_max = 0;
         for (int i=iedge_halfmax-1; i<1024; ++i) {
            Int_t nless = 0;
            // Double_t ymax_curr = y[i] + 1.*ey[i];
            // Double_t ymax_curr = y[i] + 2.*ey[i];
            Double_t ymax_curr = y[i] + 3.*ey[i];
            // Double_t ymax_curr = y[i] + 5.*ey[i];
            //if (x[i] > 65) cout<< "x[i] = " << x[i] << " y[i] = " << y[i] << " ymax_curr = " << ymax_curr <<endl;
            for (int j=i+1; j<1024; ++j) {
               //if (x[i] > 65) cout<< "x[j] = " << x[j] << " y[i] = " << y[i] <<endl;
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
         if (gdebug) {
            marker_edge_min = new TMarker(x[iedge_min],y[iedge_min],20);
            marker_edge_min->SetMarkerColor(4);
            marker_edge_min->Draw("same");
            // // redraw edge_halfmax marker
            // marker_edge_halfmax->Draw("same");
         }

         // line fit

         parline[iparline_a0] = -1000.;
         parline[iparline_a1] = 10.;
         fun_line->SetParameters(parline);

         Double_t epsline = 1e-7;
         // Double_t fit_edge_xmin = x[imin] - epsline;                 // include min
         Double_t fit_edge_xmin = x[iedge_min+1] - epsline;                 // do not include min
         Double_t fit_edge_xmax = x[iedge_max-1] + epsline;                 // do not include max
         grsig->Fit(fun_line, "", "", fit_edge_xmin, fit_edge_xmax);

         // results
         Double_t res_intercept = grsig->GetFunction(fun_line->GetName())->GetParameter(0);
         Double_t res_slope = grsig->GetFunction(fun_line->GetName())->GetParameter(1);
         Double_t intersection = 0;
         Double_t chi2line = grsig->GetFunction(fun_line->GetName())->GetChisquare();
         Double_t ndfline = grsig->GetFunction(fun_line->GetName())->GetNDF();
         Double_t chi2lineNDF = (TMath::Abs(ndfline) != 0)? chi2line/ndfline: 0;
         if (TMath::Abs(res_slope) > epsline) intersection = (-1)*res_intercept/res_slope;
         cout<< "intersection = " << intersection << " chi2line = " << chi2line << " ndfline = " << ndfline <<endl;

         // fill the tree variables
         oscFit->p0[ich] = res_intercept;
         oscFit->p1[ich] = res_slope;
         oscFit->chi2line[ich] = chi2lineNDF;
         oscFit->dline[ich] = intersection;

         //---------------------- line fit end --------------------------

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
         // Double_t pulse_y2 = 0.75*pulse_max;
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

         //-- correct high limit using result from the linear fit section
         // pulse_x2 = fit_edge_xmax;
         pulse_x2 = x[iedge_max] + epsline;

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
         // A = 3000.;
         //-- x0 = 200.;
         x0 = 60.;
         tau1 = 0.050;
         tau2 = 2.0;
         // T = 20;
         T = 40;
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
         // fun_Psigma->FixParameter(ipar_tau1, fun_Psigma->GetParameter(ipar_tau1));
         // fun_Psigma->FixParameter(ipar_sigma, fun_Psigma->GetParameter(ipar_sigma));
         // fun_Psigma->FixParameter(ipar_T, fun_Psigma->GetParameter(ipar_T));

         TF1* fitted_function = 0;
         TF1* fitted_function_refit = 0;
         npoint_fit = ipulse_x2 - ipulse_x1 + 1;
         // cout<< ".. fit pulse: pulse_max " << pulse_max << " pulse_xhalmax " << pulse_xhalfmax << " pulse_x1 = " << pulse_x1 << " pulse_x2 = " << pulse_x2 << " npoint_fit = " << npoint_fit <<endl;
         Double_t y0cross = 0;
         Double_t y05cross = 0;
         if (npoint_fit > 2)
         {
            if (gdebug) grsig->Fit(fun_Psigma, "+R", "", sigmin,sigmax);
            else grsig->Fit(fun_Psigma,"+R0", "goff", sigmin,sigmax);

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
            // fun_Psigma_refit->FixParameter(ipar_tau2, grsig->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_tau2));
            fun_Psigma_refit->FixParameter(ipar_T, grsig->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_T));
            // fun_Psigma_refit->FixParameter(ipar_sigma, grsig->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_sigma));

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

            //-- put line fit parameters first
            TF1* flinefit = (TF1*) grsig->GetListOfFunctions()->FindObject(fun_line->GetName());
            if (flinefit) {
               grsig->GetListOfFunctions()->Remove(flinefit);
               grsig->GetListOfFunctions()->AddFirst(flinefit);
            }

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
               // Double_t middle_point = pulse_xhalfmax;
               Double_t middle_point = (fit_edge_xmin + fit_edge_xmax) / 2.;
               Double_t fval = fitted_function_refit->Eval(middle_point);
               Double_t dval = fitted_function_refit->Derivative(middle_point);
               Double_t dintersect = 0;
               const Double_t eps = 1e-12;
               if (TMath::Abs(dval) > eps) dintersect = middle_point - fval/dval;
               oscFit->d[ichannel] = dintersect;
               cout<< "--> dintersect = " << dintersect << " middle_point = " << middle_point << " fval = " << fval << " dval = " << dval <<endl;

               oscFit->t[ichannel] = fitted_function_refit->GetParameter("x0");
               Double_t chi2ndf = (fitted_function_refit->GetNDF() > 0)? fitted_function_refit->GetChisquare() / fitted_function_refit->GetNDF(): 0;
               oscFit->chi2[ichannel] = chi2ndf;
               oscFit->tau1[ichannel] = fitted_function_refit->GetParameter("tau1");
               oscFit->tau2[ichannel] = fitted_function->GetParameter("tau2");
               oscFit->T[ichannel] = fitted_function_refit->GetParameter("T");
               oscFit->sigma[ichannel] = fitted_function_refit->GetParameter("sigma");

               grsig_refit->SetTitle(Form("%s #chi^{2}=%0.1f",grsig_refit->GetTitle(), chi2ndf));
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
