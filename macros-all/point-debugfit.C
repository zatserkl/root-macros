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
   // for point algorithm
   Float_t yscal[8];       // average of three points around the maximum
   Float_t pti[8];         // integral which use pt algorithm
   Float_t pty[8];
   Float_t ptx[8];
   Float_t pt3x[8];
   Float_t pt3c[8];        // chi2/NDF for 3-point fit
   Float_t pt4x[8];
   Float_t pt4c[8];        // chi2/NDF for 4-point fit
   Float_t t3x[8];
   Float_t t4x[8];
   Bool_t ptok[8];
   // flag
   Bool_t sat[8];
   Bool_t ok[8];
   Bool_t dok[8];
   Bool_t tok[8];
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
         yscal[i] = 0;
         pti[i] = 0;
         pty[i] = 0;
         ptx[i] = 0;
         pt3x[i] = 0;
         pt3c[i] = 0;
         pt4x[i] = 0;
         pt4c[i] = 0;
         t3x[i] = 0;
         t4x[i] = 0;
         ptok[i] = true;
         //
         sat[i] = false;
         ok[i] = true;
         dok[i] = true;
         tok[i] = true;
      }
   }
   OscFit(): TNamed() {
      SetNameTitle("oscFit", "oscFit");
      clear();
   }
   ClassDef(OscFit,20);
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

TF1* fpulse(Double_t xmin, Double_t xmax
      , Double_t A, Double_t x0, Double_t tau1, Double_t tau2, Double_t T, Double_t sigma
      , const char* name="fpulse"
      )
{
   TF1* fpulse = new TF1(name, fPsigma, xmin, xmax, 6);
   fpulse->SetNpx(1024);
   fpulse->SetParName(0, "A");
   fpulse->SetParName(1, "x0");
   fpulse->SetParName(2, "tau1");
   fpulse->SetParName(3, "tau2");
   fpulse->SetParName(4, "T");
   fpulse->SetParName(5, "sigma");
   fpulse->SetParameter(0, A);
   fpulse->SetParameter(1, x0);
   fpulse->SetParameter(2, tau1);
   fpulse->SetParameter(3, tau2);
   fpulse->SetParameter(4, T);
   fpulse->SetParameter(5, sigma);
   return fpulse;
}

/////////////////// pulse function end ///////////////////////////

Double_t lfit(const Float_t x[], const Float_t y[]
      , Int_t x1_i, Int_t npoints
      , Double_t& p0, Double_t& p1
      )
{
   Double_t chi2 = 0;
   p0 = p1 = 0;
   const Double_t eps = 1e-7;
   TGraph* gr_lfit = new TGraph(npoints, &x[x1_i], &y[x1_i]);
   // TGraph* gr_lfit = new TGraph(npoints+10, &x[x1_i-5], &y[x1_i-5]);
   gr_lfit->SetNameTitle("gr_lfit", "gr_lfit");
   gr_lfit->SetMarkerStyle(20);
   // new TCanvas;
   // gr_lfit->Draw("ap");
   Double_t xmin = x[x1_i] - eps;
   Double_t xmax = x[x1_i+npoints-1] + eps;
   gr_lfit->Fit("pol1", "Q0", "goff", xmin, xmax);
   // gr_lfit->Fit("pol1", "", "", xmin, xmax);
   TF1* fpol1 = gr_lfit->GetFunction("pol1");
   if (fpol1) {
      p0 = fpol1->GetParameter("p0");
      p1 = fpol1->GetParameter("p1");
      if (fpol1->GetNDF() > 0) chi2 = fpol1->GetChisquare() / fpol1->GetNDF();
   }
   delete gr_lfit;
   return chi2;
}

///////////////////////////////////////////////////////////////////
//
// main, the most general version
//
///////////////////////////////////////////////////////////////////

Double_t pol1fast(const Float_t x[], const Float_t y[], const Float_t ey[], Int_t ifirst, Int_t np, Double_t& am, Double_t& bm) 
{
   // if ey == 0 the weights are considered equal

   Double_t eps = 1e-7;
   Double_t huge = 1e10;

   Double_t chi2 = huge;
   am = 0;
   bm = 0;

   // fit of the straight line with equal weights
   Double_t wx = 0;
   Double_t wy = 0;
   Double_t wxy = 0;
   Double_t wx2 = 0;
   Double_t W = 0;

   for (int i=0; i<np; ++i)
   {
      Int_t icurr = ifirst+i;
      Float_t xi = x[icurr];
      Float_t yi = y[icurr];
      Float_t eyi = (ey)? ey[icurr]: 0;
      Float_t wi = 1.;
      if (eyi) wi = TMath::Abs(eyi) > eps? 1./(eyi*eyi): 1.;
      W += wi;
      wx += wi*xi;
      wx2 += wi*xi*xi;
      wy += wi*yi;
      wxy += wi*xi*yi;
   }

   Double_t Discriminant = W*wx2 - wx*wx;
   if (TMath::Abs(Discriminant) < eps) return huge;   // should not be happen, actually

   am = (wy*wx2 - wxy*wx) / Discriminant;    // y = am + bm*x
   bm = (W*wxy - wx*wy) / Discriminant;

   chi2 = 0;
   for (int i=0; i<np; ++i)
   {
      Int_t icurr = ifirst+i;
      Double_t yfit = am + bm*x[icurr];
      Double_t r = y[icurr] - yfit;
      if (ey) r = TMath::Abs(ey[icurr]) > eps? r / ey[icurr]: 1.;
      chi2 += r*r;
   }
   //cout<< ".. pol1fast (complete version): chi2 = " << chi2 << " am = " << am << " bm = " << bm <<endl;

   Double_t chi2ndf = np > 2? chi2/(np - 2): 0;
   return chi2ndf;
}

/////////////////////////////////////////////////////////////////////////////////
// inline ysignal(Double_t yraw, Double_t bkg_mean) {return yraw - bkg_mean);}

///////////////////////////////////////////////////////////////////////////////

void point(const Float_t xraw[], const Float_t yraw[]
      , OscFit* oscFit
      , Double_t thres = 10.
      , Bool_t debug = kFALSE
      , Int_t evt = 0
      , Int_t channel = 0        // actual channel number starts from 1
      )
{
   /*
   Good events to test
   root -l run_1_008.root.osc.root
   .L par.C+
   pulse(t, 0,0, 1,1, 1,2,3,4,5)

   some events for debug
   local_maximum:
   root -l run_1_009.root.osc.root     pulse(t, 127,127, 1,1, 1,2,3,5,6,7,8)
   problem with high bkg after subtracting (two close particles):
   root -l run_1_009.root.osc.root     pulse(t, 1817,1817, 1,1, 1,2,3,4,5,6,8)
   spike
   root -l run_1_060.root.osc.root     pulse(t, 4,4, 1,1, 1,2,3,4,5,6,8)
   example of the single spike in the data point i=1
   root -l run_1_060.root.osc.root     pulse(t, 2,2, 1,1, 1,2,3,4,5,6,7)
   */

   cout<< "------> evt = " << evt << " channel = " << channel <<endl;

   oscFit->pti[channel-1] = 0;
   oscFit->pty[channel-1] = 0;
   oscFit->ptx[channel-1] = 0;
   oscFit->pt3x[channel-1] = 0;
   oscFit->pt4x[channel-1] = 0;
   oscFit->t3x[channel-1] = 0;
   oscFit->t4x[channel-1] = 0;

   const Float_t* x = xraw;
   Float_t y[1024], ey[1024];

   // histogram for bkg
   static TH1F *hbkg;
   if (!hbkg) {
      hbkg = new TH1F("hbkg_partime", "hbkg_partime", 400, -100, 100);
   }
   hbkg->SetDirectory(0);

   // find approx value of the bkg

   // DRS4 version 3 sometime has
   // (formerly known as USB) spikes: 2 spikes of two neiboring channels, low in yraw
   // sometime I observed single-channel spike (channel 2): high in yraw

   const Int_t nbkg_approx = 20;

   // sort bkg array to eliminate possible spikes: three highest points
   Int_t index[1024];
   if (debug) cout<< "start sorting" <<endl;
   TMath::Sort(nbkg_approx, yraw, index, kFALSE);   // sort in the increasing order: spikes are high in yraw
   //for (int i=0; i<nbkg_approx; ++i) cout<< i << "\t " << yraw[index[i]] <<endl;
   if (debug) cout<< "finished sorting" <<endl;

   // look at the lowest of the data point: this can be single spike (low in yraw)
   const Double_t thres_single_spike = -100.;                              // NB high threshold value!
   Int_t single_spike_i = yraw[index[0]] < thres_single_spike? index[0]: -1;

   Float_t bkg_mean_approx = 0;
   //Float_t bkg_sigma_approx = 0;
   for (int i=1; i<nbkg_approx-4; ++i) bkg_mean_approx += yraw[index[i]];  // omit first 1 and last 4
   bkg_mean_approx /= (nbkg_approx-1-4);
   bkg_mean_approx /= -1.;                          // invert bkg_mean built out of yraw
   cout<< "approx value of the bkg is " << bkg_mean_approx <<endl;

   Int_t maximum_i = nbkg_approx;
   Float_t maximum_x = x[maximum_i];
   Float_t maximum_y = yraw[maximum_i];

   bool passed_thres_rise = false;     // exceeded thres at the rising (left) slope
   bool passed_peak = false;           // passed the peak

   std::vector<TMarker*> marker_spikes;

   // in the following loop over all the points:
   //
   // 1) convert yraw to y using approx value of the bkg
   // 2) correct spikes (if any)
   // 3) find the maximum (consider data points over the thres only)
   // 4) find the right peak boundaty
   //    -- 1/8 of the maximum
   //    -- process the first peak only (in case of event with multiple peaks)

   //cout<< "apply approx threshold: start" <<endl;
   Bool_t saturation = kFALSE;
   Double_t eof_peak = yraw[maximum_i]/8.;
   for (int i=0; i<1024; ++i)
   {
      if (y[i] < -499.) saturation = kTRUE;

      y[i] = -1.*yraw[i] - bkg_mean_approx;

      // if (i == single_spike_i) {
      //    y[i] = 0;
      //    continue;
      // }

      if (i >= 3)
      {
         // identify spikes. NB: spikes are down ==> do not affect the maximum
         const Double_t delta_x = x[1] - x[0];
         Double_t tangent = (y[i] - y[i-3]) / (3.*delta_x);
         Double_t delta_y = tangent*delta_x;
         Double_t interpol_y1 = y[i-3] + delta_y;
         Double_t dy1 = interpol_y1 - y[i-2];
         if (dy1 > 10. && dy1 < 15.) {
            Double_t interpol_y2 = interpol_y1 + delta_y;
            Double_t dy2 = interpol_y2 - y[i-1];
            if (dy2 > 10. && dy2 < 15.) {
               if (TMath::Abs(dy1 - dy2) < 5.) {
                  TMarker* marker_spike1 = new TMarker(x[i-2], y[i-2], 20);
                  marker_spike1->SetMarkerColor(8);
                  marker_spikes.push_back(marker_spike1);
                  TMarker* marker_spike2 = new TMarker(x[i-1], y[i-1], 20);
                  marker_spike2->SetMarkerColor(8);
                  marker_spikes.push_back(marker_spike2);
                  y[i-2] = interpol_y1;
                  y[i-1] = interpol_y2;
                  if (debug) cout<< "... corrected spikes at x[" << i-2 << "] = " << x[i-2] <<endl;
               }
            }
         }
      }
      if (!passed_peak
            && i > nbkg_approx
            && y[i] > thres
            && y[i] > maximum_y
            )
      {
         passed_thres_rise = true;
         maximum_i = i;
         maximum_x = x[i];
         maximum_y = y[i];
         eof_peak = maximum_y / 8.;
      }
      if (!passed_peak
            && passed_thres_rise
            && y[i] < eof_peak
            )
      {
         if (debug) cout<< "... passed the peak: x[" << i << "] " << x[i] << " y " << y[i] << " < eof_peak " << eof_peak <<endl;
         passed_peak = true;
      }
   }
   if (!passed_thres_rise) {
      cout<< "point: no pulse found" <<endl;
      return;
   }
   else cout<< "maximum_x " << maximum_x << " maximum_y " << maximum_y << " maximum_i " << maximum_i <<endl;

   // find width at the half maximum

   Double_t halfmax = maximum_y / 2.;
   Double_t quotermax = maximum_y / 4.;
   Int_t halfmax_x1_i = 0;
   Int_t halfmax_x2_i = 0;

   for (int i=0; i<1024; ++i)
   {
      if (i == single_spike_i) continue;

      if (true                // find the first point which exceeds half-max at the leading edge
            && halfmax_x1_i == 0
            && y[i] > halfmax
            )
      {
         halfmax_x1_i = i;
      }
      if (true                // find the last point which exceeds half-max at the tail
            && i > maximum_i
            && y[i] > halfmax
            )
      {
         halfmax_x2_i = i;
      }
      if (true                // break when down the quoter-max to avoid the second peak (if any)
            && i > maximum_i
            && y[i] < quotermax
         )
      {
         break;
      }
   }
   // adjust indices
   //halfmax_x1_i--;
   halfmax_x2_i++;
   if (halfmax_x2_i - halfmax_x1_i < 3) {
      ++halfmax_x2_i;
      halfmax_x1_i = halfmax_x2_i - 3;
   }
   if (debug) cout<< "left peak slope: " << x[halfmax_x1_i] << " right peak slope: " << x[halfmax_x2_i] <<endl;;

   // three sigma to the left
   // Exactly: sigma = FWHM/2.35 = delta_x * 2*(maximum_i - halfmax_x1_i + 1) / 2.35
   // Approx: sigma = delta_x * (maximum_i - halfmax_x1_i + 1)

   Int_t i3sigma1 = 3 * (maximum_i - halfmax_x1_i);   // 3 sigma1 to the left
   Int_t i3sigma2 = 3 * (halfmax_x2_i - maximum_i);   // 3 sigma2 to the right
   Int_t isigmin = maximum_i - i3sigma1;
   if (isigmin < 10) isigmin = 10;                    // isigmin >= 10
   Int_t isigmax = maximum_i + i3sigma2;
   if (isigmax > 1023) isigmax = 1023;
   Int_t ibkgmax = isigmin - 1;
   Int_t ibkgmin = isigmin < 40? 0: isigmin - 40;        // will be used for the plot
   if (debug) cout<< "i3sigma1 " << i3sigma1 << " i3sigma2 " << i3sigma2 << " x[isigmin] " << x[isigmin] << " x[isigmax] " << x[isigmax] <<endl;
   if (debug) cout<< "left edge: " << x[isigmin] << " right edge: " << x[isigmax] <<endl;

   // find the exact bkg
   Int_t nbkg = ibkgmax;
   // /* NB: we already corrected double spikes (formerly known as USB) */ TMath::Sort(nbkg, y, index, kTRUE); // sort in the decreasing order: spikes are high in yraw
   TMath::Sort(nbkg, y, index, kFALSE); // sort in the increasing order: the single spike is low in yraw
   hbkg->Reset();
   //for (int i=1; i<nbkg-4; ++i) hbkg->Fill(y[index[i]]);    // omit first 1 and last 4
   for (int i=1; i<nbkg-1; ++i) hbkg->Fill(y[index[i]]);    // omit last 1 -- protection against the single spike
   hbkg->Fit("gaus", "L0Q", "goff");                        // NB: Loglikelihood option

   Double_t bkg_mean = hbkg->GetFunction("gaus")->GetParameter("Mean");
   Double_t bkg_sigma = hbkg->GetFunction("gaus")->GetParameter("Sigma");
   cout<< "bkg_mean = " << bkg_mean << " bkg_sigma = " << bkg_sigma <<endl;   

   // correct the data with the updated bkg level
   Float_t integral = 0;
   Float_t yprev = y[isigmin] - bkg_mean;
   for (int i=0; i<1024; ++i) {
      y[i] -= bkg_mean;
      ey[i] = bkg_sigma;
      if (i > isigmin && i <= isigmax) {
         integral += 0.5*(yprev + y[i]);
         yprev = y[i];
      }
   }
   maximum_y -= bkg_mean;
   // complete integral
   integral *= x[1] - x[0];
   // normalize integral to pC by division to R = 50 Ohm
   integral /= 50.;
   cout<< "pulse integral = " << integral << " pC" <<endl;

   Int_t np = isigmax-ibkgmin+1;
   TGraphErrors* gr = new TGraphErrors(np, &x[ibkgmin], &y[ibkgmin], 0, &ey[ibkgmin]);
   gr->SetNameTitle("gr", Form("evt %d, ch %d", evt, channel));
   gr->SetMarkerStyle(24);
   gr->SetMarkerColor(8);
   gr->SetLineColor(8);
   if (debug) {
      new TCanvas;
      gr->Draw("ap");
      // spikes
      for (unsigned imarker=0; imarker<marker_spikes.size(); ++imarker) marker_spikes[imarker]->Draw("same");
   }

   // time stamp

   // find the local maximum

   // Double_t point_first_y = maximum_y/4.;
   Double_t point_first_y = maximum_y/8.;
   Int_t point_first_i = maximum_i;
   do --point_first_i; while (point_first_i > isigmin && y[point_first_i] > point_first_y);

   if (debug) cout<< "point_first_i " << point_first_i << " x[point_first_i] " << x[point_first_i] <<endl; 
   Int_t local_maximum_i = maximum_i;

   const Float_t eps = 1e-7;  // used for fit limits and against div by 0 at ptx

   // gaussian fit -- good for test, but slow
   //
   // for (int i=point_first_i; i<isigmax-6; ++i)
   // {
   //    Int_t point_last_i = i + 5;
   //    if (point_last_i == maximum_i) {
   //       cout<< "... local_maximum_i == maximum_i" <<endl;
   //       local_maximum_i = maximum_i;
   //       break;
   //    }
   //    Double_t fit_xmin = x[i] - eps;
   //    Double_t fit_xmax = x[point_last_i] + eps;
   //    cout<< "... fit from xmin " << fit_xmin << " fit_xmax " << fit_xmax <<endl;
   //    //-- gr->Fit("gaus", "0Q", "goff", fit_xmin, fit_xmax);
   //    // gr->Fit("gaus", "", "", fit_xmin, fit_xmax);
   //    gr->Fit("gaus", "Q", "", fit_xmin, fit_xmax);
   //    Double_t gmean = gr->GetFunction("gaus")->GetParameter("Mean");

   //    Int_t point_third_i = i + 2;
   //    if (point_third_i >= maximum_i) {
   //       cout<< "... point_third_i " << point_third_i << " <= maximum_i " << maximum_i <<endl;
   //       local_maximum_i = maximum_i;
   //       break;
   //    }
   //    else if (gmean < x[point_third_i]) {
   //       cout<< "... x[point_third_i] " << x[point_third_i] << " > gmean " << gmean <<endl;
   //       local_maximum_i = point_third_i;
   //       break;
   //    }
   // }

   // set local_maximum as the maximum of current 5 points -- I think, this is extra
   //
   // for (int i=point_first_i; i<isigmax-6; ++i)
   // {
   //    Int_t point_last_i = i + 5;
   //    if (point_last_i == maximum_i) {
   //       if (debug) cout<< "... local_maximum_i == maximum_i" <<endl;
   //       local_maximum_i = maximum_i;
   //       break;
   //    }

   //    Double_t ycurr = y[i] + ey[i];
   //    Int_t nlow = 0;
   //    Int_t maximum_current_i = i;
   //    for (int ii=i+1; ii<point_last_i; ++ii)
   //    {
   //       if (y[ii] > y[maximum_current_i]) maximum_current_i = ii;
   //       if (y[ii] < ycurr) ++nlow;
   //    }
   //    if (nlow > 3) {
   //       // this is a local maximum
   //       local_maximum_i = maximum_current_i;
   //       break;
   //    }
   // }

   // for (int i=point_first_i; i<isigmax-6; ++i)
   // {
   //    Int_t point_last_i = i + 5;
   //    if (point_last_i == maximum_i) {
   //       if (debug) cout<< "... local_maximum_i == maximum_i" <<endl;
   //       local_maximum_i = maximum_i;
   //       break;
   //    }

   //    Double_t ycurr = y[i] + ey[i];
   //    Int_t nlow = 0;
   //    for (int ii=i+1; ii<point_last_i; ++ii)
   //    {
   //       if (y[ii] < ycurr) ++nlow;
   //    }
   //    if (nlow > 3) {
   //       // this is a local maximum
   //       local_maximum_i = i;
   //       break;
   //    }
   // }
   for (int i=point_first_i; i<isigmax-6; ++i)
   {
      Int_t ncheck = 5;
      if (i+ncheck == maximum_i) {
         if (debug) cout<< "... local_maximum_i == maximum_i" <<endl;
         local_maximum_i = maximum_i;
         break;
      }

      Double_t ycurr = y[i] + ey[i];         // current level to compare with
      Int_t nlow = 0;
      for (int ii=i+1; ii<i+ncheck; ++ii)
      {
         if (y[ii] < ycurr) ++nlow;
      }
      if (nlow > 3) {
         // 3 points out of ncheck are lower: this is a local maximum
         local_maximum_i = i;
         break;
      }
   }
   Double_t local_maximum_x = x[local_maximum_i];
   Double_t local_maximum_y = y[local_maximum_i];

   if (local_maximum_i > maximum_i) {
      if (debug) cout<< "... ---> local_maximum_i > maximum_i. Set to maximum" <<endl;
      local_maximum_i = maximum_i;
      local_maximum_x = maximum_x;
      local_maximum_y = maximum_y;
   }

   if (debug) {
      TMarker* marker_local_maximum = new TMarker(local_maximum_x, local_maximum_y, 20);
      marker_local_maximum->SetMarkerColor(4);
      marker_local_maximum->SetMarkerSize(1);
      marker_local_maximum->Draw("same");
   }

   cout<< "local_maximum_x " << local_maximum_x << " local_maximum_y " << local_maximum_y <<endl;

   // find half-max
   Double_t pty = local_maximum_y / 2.;
   Int_t x1_i = local_maximum_i;
   do x1_i--; while (x1_i > isigmin && y[x1_i] > pty);
   Int_t x2_i = x1_i + 1;
   if (debug) cout<< "x[x1_i] " << x[x1_i] << " x[x2_i] " << x[x2_i] <<endl;
   Double_t delta_x = x[x2_i] - x[x1_i];
   Double_t delta_y = y[x2_i] - y[x1_i];           // NB: delta_y > 0 in regular case
   Double_t ptx = 0;
   if (delta_y > eps) ptx = x[x2_i] - delta_x * (y[x2_i] - pty) / delta_y;
   if (ptx < x[x1_i] || ptx > x[x2_i]) {
      cout<< "--> warning: evt " << evt << " channel " << channel << ": ptx is out of limits" <<endl;
      ptx = 0;    // something wrong
   }
   cout<< "ptx " << ptx << " pty " << pty <<endl;

   // fit straight line

   Double_t chi2_1, chi2_2;
   Double_t chi2_min3 = 0;
   Double_t p0_1, p0_2, p1_1, p1_2;
   Double_t p0_best = 0, p1_best = 0;
   Int_t npfit = 0;
   // fit straight line through 3 points
   npfit = 3;
   Int_t i1 = x1_i-1;
   Int_t i2 = x1_i;
   // cout<< "fit1: x[i1] " << x[i1] << " to x[i1+npfit-1] " << x[i1+npfit-1] <<endl;
   //chi2_1 = lfit(x,y, i1, npfit, p0_1, p1_1);
   chi2_1 = pol1fast(x,y,ey, i1, npfit, p0_1, p1_1);
   // cout<< "fit2: x[i2] " << x[i2] << " to x[i2+npfit-1] " << x[i2+npfit-1] <<endl;
   //chi2_2 = lfit(x,y, i2, npfit, p0_2, p1_2);
   chi2_2 = pol1fast(x,y,ey, i2, npfit, p0_2, p1_2);
   // cout<< "chi2_1 " << chi2_1 << " chi2_2 " << chi2_2 <<endl;
   // Double_t pt3x = chi2_1 < chi2_2? (pty - p0_1)/p1_1: (pty - p0_2)/p1_2;

   Double_t pt3x = 0;            // x-coordinate of the middle slope point
   Double_t t3x = 0;             // x of intersection of the straight line with x-axis
   if (chi2_1 < chi2_2) {
      chi2_min3 = chi2_1;
      pt3x = (pty - p0_1)/p1_1;
      t3x = -p0_1/p1_1;
      i1 = i1;
      i2 = i1 - 1;
   }
   else {
      chi2_min3 = chi2_2;
      pt3x = (pty - p0_2)/p1_2;
      t3x = -p0_2/p1_2;
      i1 = i2;
      i2 = i1 - 1;
   }
   cout<< "pt3x = " << pt3x << " pty = " << pty << " x[i1] = " << x[i1] <<endl;

   // fit two points above and two points below
   npfit = 4;
   chi2_min3 = pol1fast(x,y,ey, x1_i-1, npfit, p0_best, p1_best);
   pt3x = (pty - p0_2)/p1_2;
   t3x = -p0_2/p1_2;
   //
   p0_best = p0_best;
   p1_best = p1_best;

   // fit straight line through 4 points
   npfit = 4;
   // cout<< "fit1: x[i1] " << x[i1] << " to x[i1+npfit-1] " << x[i1+npfit-1] <<endl;
   //chi2_1 = lfit(x,y, i1, npfit, p0_1, p1_1);
   chi2_1 = pol1fast(x,y,ey, i1, npfit, p0_1, p1_1);
   // cout<< "fit2: x[i2] " << x[i2] << " to x[i2+npfit-1] " << x[i2+npfit-1] <<endl;
   //chi2_2 = lfit(x,y, i2, npfit, p0_2, p1_2);
   chi2_2 = pol1fast(x,y,ey, i2, npfit, p0_2, p1_2);
   // cout<< "chi2_1 " << chi2_1 << " chi2_2 " << chi2_2 <<endl;
   // Double_t pt3x = chi2_1 < chi2_2? (pty - p0_1)/p1_1: (pty - p0_2)/p1_2;

   TF1* pol1_4 = 0;

   Double_t pt4x = 0;            // x-coordinate of the middle slope point
   Double_t t4x = 0;             // x of intersection of the straight line with x-axis
   Double_t chi2_min4 = 0;
   if (chi2_1 < chi2_2) {
      chi2_min4 = chi2_1;
      pt4x = (pty - p0_1)/p1_1;
      t4x = -p0_1/p1_1;
      i1 = i1;
      //i2 = i1 - 1;
      i2 = i1+npfit - 1;
      // //
      // p0_best = p0_1;
      // p1_best = p1_1;
   }
   else {
      chi2_min4 = chi2_2;
      pt4x = (pty - p0_2)/p1_2;
      t4x = -p0_2/p1_2;
      i1 = i2;
      //i2 = i1 - 1;
      i2 = i1+npfit - 1;
      // //
      // p0_best = p0_2;
      // p1_best = p1_2;
   }
   cout<< "pt4x = " << pt4x << " pty = " << pty << " x[i1] = " << x[i1] <<endl;

   if (debug) {
      // pol1_4 = new TF1("pol1_4","pol1", x[i1],x[i2]);
      pol1_4 = new TF1("pol1_4","pol1", x[x1_i-1],x[x2_i+1]);
      pol1_4->SetNpx(1024);
      pol1_4->SetLineWidth(1);
      pol1_4->SetLineColor(2);
      pol1_4->SetParameter(0, p0_best);
      pol1_4->SetParameter(1, p1_best);
      pol1_4->Draw("same");
   }

   if (debug) {
      TMarker* marker_pt = new TMarker(ptx, pty, 5);
      marker_pt->SetMarkerColor(2);
      marker_pt->SetMarkerSize(2);
      marker_pt->Draw("same");
   }

   // assign output values
   oscFit->yscal[channel-1] = (y[maximum_i-1] + y[maximum_i] + y[maximum_i+1]) / 3.;
   oscFit->pti[channel-1] = integral;
   oscFit->sat[channel-1] = saturation;
   oscFit->pty[channel-1] = pty;
   oscFit->ptx[channel-1] = ptx;
   oscFit->t3x[channel-1] = t3x;
   oscFit->pt3c[channel-1] = chi2_min3;
   oscFit->pt3x[channel-1] = pt3x;
   oscFit->t4x[channel-1] = t4x;
   oscFit->pt4x[channel-1] = pt4x;
   oscFit->pt4c[channel-1] = chi2_min4;
}
/////////////////////////////////////////////////////////////////////////////////

bool debug = false;
bool gdebug = false;

TTree* pulseosc(TTree *tree
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

   // the number of channels in the data
   tree->GetEntry(entry_first);
   Int_t nchannels = oscEvent->oscChannels->GetEntries();

   if (entry_last < 0) entry_last = tree->GetEntries() - 1;

   for (int jentry=entry_first; jentry<=entry_last; ++jentry)
   {
      tree->LoadTree(jentry);
      tree->GetEntry(jentry);
      cout<< "\n---------------------------------------------------> jentry = " << jentry << " evt = " << oscEvent->evt <<endl;

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
         if (oscChannel->ch == ignorech4) continue;
         if (oscChannel->ch == ignorech5) continue;
         if (oscChannel->ch == ignorech6) continue;
         if (oscChannel->ch == ignorech7) continue;
         if (oscChannel->ch == ignorech8) continue;

         const Double_t thres = 10;
         //-- const Double_t thres = 5;
         
         //
         // pt algorithm
         //
         point(oscChannel->x, oscChannel->y
               , oscFit
               , thres, debug
               , jentry, oscChannel->ch
               );
         oscFit->ptok[ich] = true
            && oscFit->pti[ich] > 0
            && oscFit->ptx[ich] > 0
            && oscFit->sat[ich] == kFALSE
         ;
      }
      otree->Fill();
   }

   cout<< "\ngDirectory->GetName() = " << gDirectory->GetName() <<endl;

   // turn off debug
   debug = false;
   gdebug = false;

   return otree;
}

TTree* pulseosc(const char *ifname
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

   TTree* otree = pulseosc(tree, entry_first,entry_last, setdebug,setgdebug, ignorech1,ignorech2,ignorech3,ignorech4,ignorech5,ignorech6,ignorech7,ignorech8);
   
   ofile->Write();
   return otree;
}

TTree* pulse3(TTree *tree
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

   Float_t b1_t[1024], b1_c1[1024], b1_c2[1024], b1_c3[1024], b1_c4[1024];
   Float_t b2_t[1024], b2_c1[1024], b2_c2[1024], b2_c3[1024], b2_c4[1024];
   Int_t event, tc1, tc2;

   tree->SetMarkerStyle(7);
   tree->SetMarkerColor(2);

   //OscEvent* oscEvent = 0;
   //tree->SetBranchAddress("oscEvent",&oscEvent);

   TBranch* t_event = tree->GetBranch("event");
   TBranch* t_tc1   = tree->GetBranch("tc1");
   TBranch* t_b1_t  = tree->GetBranch("b1_t");
   TBranch* t_b1_c1 = tree->GetBranch("b1_c1");
   TBranch* t_b1_c2 = tree->GetBranch("b1_c2");
   TBranch* t_b1_c3 = tree->GetBranch("b1_c3");
   TBranch* t_b1_c4 = tree->GetBranch("b1_c4");

   TBranch* t_tc2   = tree->GetBranch("tc2");
   TBranch* t_b2_t  = tree->GetBranch("b2_t");
   TBranch* t_b2_c1 = tree->GetBranch("b2_c1");
   TBranch* t_b2_c2 = tree->GetBranch("b2_c2");
   TBranch* t_b2_c3 = tree->GetBranch("b2_c3");
   TBranch* t_b2_c4 = tree->GetBranch("b2_c4");

   t_event->SetAddress( &event);
   t_tc1->SetAddress( &tc1);
   t_b1_t->SetAddress( b1_t);
   t_b1_c1->SetAddress( b1_c1);
   t_b1_c2->SetAddress( b1_c2);
   t_b1_c3->SetAddress( b1_c3);
   t_b1_c4->SetAddress( b1_c4);

   if (t_b2_t) {
      cout<< "--> Found board #2" <<endl;
      t_tc2->SetAddress( &tc2);
      t_b2_t->SetAddress( b2_t);
      t_b2_c1->SetAddress( b2_c1);
      t_b2_c2->SetAddress( b2_c2);
      t_b2_c3->SetAddress( b2_c3);
      t_b2_c4->SetAddress( b2_c4);
   }

   cout<< "tree->GetEntries() = " << tree->GetEntries() <<endl;

   // output (fit results) tree
   TTree* otree = new TTree("ft", "Fit result tree");
   OscFit* oscFit = new OscFit;
   otree->Branch("oscFit", "OscFit", &oscFit);
   otree->SetMarkerStyle(6);

   // the number of channels in the data
   tree->GetEntry(entry_first);
   //Int_t nchannels = oscEvent->oscChannels->GetEntries();
   Int_t nchannels = 4;
   if(t_b2_t) nchannels = 8;

   if (entry_last < 0) entry_last = tree->GetEntries() - 1;

   for (int jentry=entry_first; jentry<=entry_last; ++jentry)
   {
      tree->LoadTree(jentry);
      tree->GetEntry(jentry);
      cout<< "\n---------------------------------------------------> jentry = " << jentry <<endl;

      oscFit->clear();
      oscFit->evt = jentry;

      //cout<< ".. loop over osc channels. oscEvent->oscChannels = " << oscEvent->oscChannels->GetEntries() <<endl;
      for (int ich=0; ich<nchannels; ++ich)
      {
         //const OscChannel* oscChannel = (OscChannel*) oscEvent->oscChannels->At(ich);
         //cout<< "oscChannel->ch = " << oscChannel->ch <<endl;

         if (ich+1 == ignorech1) continue;
         if (ich+1 == ignorech2) continue;
         if (ich+1 == ignorech3) continue;
         if (ich+1 == ignorech4) continue;
         if (ich+1 == ignorech5) continue;
         if (ich+1 == ignorech6) continue;
         if (ich+1 == ignorech7) continue;
         if (ich+1 == ignorech8) continue;

         const Double_t thres = 10;
         //-- const Double_t thres = 5;
         
         //
         // pt algorithm
         //
         if (ich == 0) point(b1_t, b1_c1, oscFit, thres, debug, jentry, ich+1);
         if (ich == 1) point(b1_t, b1_c2, oscFit, thres, debug, jentry, ich+1);
         if (ich == 2) point(b1_t, b1_c3, oscFit, thres, debug, jentry, ich+1);
         if (ich == 3) point(b1_t, b1_c4, oscFit, thres, debug, jentry, ich+1);

         if (ich > 3 && t_b2_t == 0) continue;
         if (ich == 4) point(b2_t, b2_c1, oscFit, thres, debug, jentry, ich+1);
         if (ich == 5) point(b2_t, b2_c2, oscFit, thres, debug, jentry, ich+1);
         if (ich == 6) point(b2_t, b2_c3, oscFit, thres, debug, jentry, ich+1);
         if (ich == 7) point(b2_t, b2_c4, oscFit, thres, debug, jentry, ich+1);

         oscFit->ptok[ich] = true
            && oscFit->pti[ich] > 0
            && oscFit->ptx[ich] > 0
            && oscFit->sat[ich] == kFALSE
         ;
      }
      otree->Fill();
   }

   cout<< "\ngDirectory->GetName() = " << gDirectory->GetName() <<endl;

   // turn off debug
   debug = false;
   gdebug = false;

   return otree;
}

TTree* pulse3(const char *ifname
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

   TTree* tree = (TTree*) ifile->Get("pulse");
   if (!tree) {
      cout<< "tree \"pulse\" was not found in file " << ifname <<endl;
      return 0;
   }

   // output file with tree "ft" ("fit tree")
   TFile* ofile = TFile::Open(Form("%s-ft.root",ifname),"recreate");

   TTree* otree = pulse3(tree, entry_first,entry_last, setdebug,setgdebug, ignorech1,ignorech2,ignorech3,ignorech4,ignorech5,ignorech6,ignorech7,ignorech8);
   
   ofile->Write();
   return otree;
}
