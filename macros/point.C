#include <TROOT.h>
#include <TStyle.h>
#include <TEnv.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TDirectory.h>
#include <TMarker.h>
#include <TClonesArray.h>
#include <TCanvas.h>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TEventList.h>
#include <TTimer.h>

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <map>

using std::cout;     using std::endl;

const char* nextname(const char* base);
const char* nextcan(const char* base);

class OscFit: public TNamed
{
public:
   Int_t evt;
   Float_t thres[8];
   Float_t trig[8];
   Float_t gate[8];
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
   Float_t tq[8][30];      // time when the area under the fpulse reached 1, 2, 3, ...
   // line fit
   Float_t p0[8];          // slope in y = a*x + b
   Float_t p1[8];          // intercept
   Float_t dline[8];       // intersection of the fitted line with y-axis
   Float_t chi2line[8];    // chi2 of line fit
   Float_t npline[8];
   // for point algorithm
   Float_t sigmin[8];
   Float_t sigmax[8];
   Float_t yscal[8];       // average of three points around the maximum
   Float_t lmaxx[8];       // x position of the local maximum
   Float_t lmaxy[8];       // y position of the local maximum
   Float_t pmax[8];        // pulse maximum
   Float_t pti[8];         // integral which use pt algorithm
   Float_t pty[8];
   Float_t ptx[8];
   Float_t plx[8];
   Float_t plx0[8];        // intersection with the x-axis
   Float_t plc[8];         // chi2/NDF for 4-point fit
   Bool_t ptok[8];
   // fit
   Double_t plxmin[8];     // fit limits
   Double_t plxmax[8];     // fit limits
   Double_t plp0[8];       // y = p0 + p1*x
   Double_t plp1[8];       // y = plp0[i] + plp1[i]*x
   // flag
   Bool_t sat[8];
   Bool_t ok[8];
   Bool_t dok[8];
   Bool_t tok[8];
   //
   Float_t itot[8];
   Float_t ibkg[8];
   Float_t itail[8];
   // rise time
   Float_t trisef[8];      // pulse rise time from 10% to 90% of pmax from pulse function
public:
   void clear() {
      evt = 0;
      for (int i=0; i<8; ++i) {
         thres[i] = 0;
         trig[i] = 0;
         gate[i] = 0;
         trig[i] = 0;
         gate[i] = 0;
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
         for (int j=0; j<30; ++j) tq[i][j] = 0;

         p0[i] = 0;
         p1[i] = 0;
         dline[i] = 0;
         chi2line[i] = 10;
         npline[i] = 0;
         //
         sigmin[i] = 0;
         sigmax[i] = 0;
         yscal[i] = 0;
         lmaxx[i] = 0;
         lmaxy[i] = 0;
         pmax[i] = 0;
         pti[i] = 0;
         pty[i] = 0;
         ptx[i] = 0;
         plx[i] = 0;
         plx0[i] = 0;
         plc[i] = 0;
         ptok[i] = kTRUE;
         //
         plxmin[i] = 0;
         plxmax[i] = 0;
         plp0[i] = 0;
         plp1[i] = 0;
         //
         sat[i] = kFALSE;
         ok[i] = kTRUE;
         dok[i] = kTRUE;
         tok[i] = kTRUE;
         //
         itot[i] = 0;
         ibkg[i] = 0;
         itail[i] = 0;
         // rise time
         trisef[i] = 0;
      }
   }
   OscFit(): TNamed() {
      SetNameTitle("oscFit", "oscFit");
      clear();
   }
   ClassDef(OscFit,34);
};

class OscChannel: public TObject
{
   public:
      Int_t ch;
      Float_t x[1024];
      Float_t y[1024];

   public:
      void clear() {
         ch = 0;
         for (int i=0; i<1024; ++i) {
            x[i] = 0;
            y[i] = 0;
         }
      }
      OscChannel(): TObject() {
         clear();
      }
      ~OscChannel() {}
      ClassDef(OscChannel, 1)
};

class OscEvent: public TObject
{
   public:
      Int_t evt;
      TClonesArray* oscChannels;          //-> array of OscChannel

   public:
      void clear() {
         evt = 0;
         oscChannels->Clear();
      }
      OscEvent(): TObject() {
         oscChannels = new TClonesArray("OscChannel");
         clear();
      }
      ~OscEvent() {
         delete oscChannels;   oscChannels = 0;
      }
      const OscChannel* oscChannel(int channelNo) const {
         for (int ich=0; ich<oscChannels->GetEntries(); ++ich) {
            const OscChannel* channel = (const OscChannel*) oscChannels->At(ich);
            if (channel->ch == channelNo) return channel;
         }
         return 0;
      }
      ClassDef(OscEvent, 2)
};

#ifdef __MAKECINT__
#pragma link C++ class OscChannel;
#pragma link C++ class OscEvent;
#pragma link C++ class OscFit;
#endif

ClassImp(OscChannel);
ClassImp(OscEvent);
ClassImp(OscFit);

struct PulseBuffer
{
   Int_t event;
   Int_t year;
   Int_t month;
   Int_t day;
   Int_t hour;
   Int_t minute;
   Int_t second;
   Int_t millisecond;
   Int_t tc1;              // not in use, for compatibility with pulse tree only
   Float_t b1_t[1024];
   Float_t b1_c1[1024];
   Float_t b1_c2[1024];
   Float_t b1_c3[1024];
   Float_t b1_c4[1024];
   Int_t tc2;              // not in use, for compatibility with pulse tree only
   Float_t b2_t[1024];
   Float_t b2_c1[1024];
   Float_t b2_c2[1024];
   Float_t b2_c3[1024];
   Float_t b2_c4[1024];
   Int_t usedchan[4];
   TBranch* b_event;
   TBranch* b_year;
   TBranch* b_month;
   TBranch* b_day;
   TBranch* b_hour;
   TBranch* b_minute;
   TBranch* b_second;
   TBranch* b_millisecond;
   TBranch* b_tc1;              // not in use, for compatibility with pulse tree only
   TBranch* b_b1_t;
   TBranch* b_b1_c1;
   TBranch* b_b1_c2;
   TBranch* b_b1_c3;
   TBranch* b_b1_c4;
   TBranch* b_tc2;              // not in use, for compatibility with pulse tree only
   TBranch* b_b2_t;
   TBranch* b_b2_c1;
   TBranch* b_b2_c2;
   TBranch* b_b2_c3;
   TBranch* b_b2_c4;
   TBranch* b_usedchan;
   std::map<TBranch*, void*> baddr;
   void clear() {
      event = 0;
      year = 0;
      month = 0;
      day = 0;
      hour = 0;
      minute = 0;
      second = 0;
      millisecond = 0;
      tc1 = 0;
      for (int ich=0; ich<4; ++ich) usedchan[ich] = 0;
      for (int ipoint=0; ipoint<1024; ++ipoint) {
         b1_t[ipoint] = 0;
         b1_c1[ipoint] = 0;
         b1_c2[ipoint] = 0;
         b1_c3[ipoint] = 0;
         b1_c4[ipoint] = 0;
         b2_t[ipoint] = 0;
         b2_c1[ipoint] = 0;
         b2_c2[ipoint] = 0;
         b2_c3[ipoint] = 0;
         b2_c4[ipoint] = 0;
      }
   }
   ~PulseBuffer() {
   }
   PulseBuffer():
      b_event(0)
      , b_year(0)
      , b_month(0)
      , b_day(0)
      , b_hour(0)
      , b_minute(0)
      , b_second(0)
      , b_millisecond(0)
      , b_tc1(0)              // not in use, for compatibility with pulse tree only
      , b_b1_t(0)
      , b_b1_c1(0)
      , b_b1_c2(0)
      , b_b1_c3(0)
      , b_b1_c4(0)
      , b_tc2(0)              // not in use, for compatibility with pulse tree only
      , b_b2_t(0)
      , b_b2_c1(0)
      , b_b2_c2(0)
      , b_b2_c3(0)
      , b_b2_c4(0)
      , b_usedchan(0)
   {
      clear();
   }
   PulseBuffer(TTree* tree) {
      clear();
      Connect(tree);
   }
   PulseBuffer(const PulseBuffer& buf) {
      operator =(buf);
   }
   PulseBuffer& operator =(const PulseBuffer& buf) {
      if (&buf == this) return *this;
      event = buf.event;
      year = buf.year;
      month = buf.month;
      day = buf.day;
      hour = buf.hour;
      minute = buf.minute;
      second = buf.second;
      millisecond = buf.millisecond;
      for (int ich=0; ich<4; ++ich) usedchan[ich] = buf.usedchan[ich];
      for (int ipoint=0; ipoint<1024; ++ipoint) {
         b1_t[ipoint] = buf.b1_t[ipoint];
         b1_c1[ipoint] = buf.b1_c1[ipoint];
         b1_c2[ipoint] = buf.b1_c2[ipoint];
         b1_c3[ipoint] = buf.b1_c3[ipoint];
         b1_c4[ipoint] = buf.b1_c4[ipoint];
         b2_t[ipoint] = buf.b2_t[ipoint];
         b2_c1[ipoint] = buf.b2_c1[ipoint];
         b2_c2[ipoint] = buf.b2_c2[ipoint];
         b2_c3[ipoint] = buf.b2_c3[ipoint];
         b2_c4[ipoint] = buf.b2_c4[ipoint];
      }
      b_event = buf.b_event;
      b_year = buf.b_year;
      b_month = buf.b_month;
      b_day = buf.b_day;
      b_hour = buf.b_hour;
      b_minute = buf.b_minute;
      b_second = buf.b_second;
      b_millisecond = buf.b_millisecond;
      b_tc1 = buf.b_tc1;
      b_b1_t = buf.b_b1_t;
      b_b1_c1 = buf.b_b1_c1;
      b_b1_c2 = buf.b_b1_c2;
      b_b1_c3 = buf.b_b1_c3;
      b_b1_c4 = buf.b_b1_c4;
      b_tc2 = buf.b_tc2;
      b_b2_t = buf.b_b2_t;
      b_b2_c1 = buf.b_b2_c1;
      b_b2_c2 = buf.b_b2_c2;
      b_b2_c3 = buf.b_b2_c3;
      b_b2_c4 = buf.b_b2_c4;
      b_usedchan = buf.b_usedchan;
      return *this;
   }
   void UsedChan() const {
      cout<< "DRS4 channels: ";
      for (int ich=0; ich<4; ++ich) {
         if (usedchan[ich]) cout<< usedchan[ich] << "  ";
      }
      cout<<endl;
   }
   void Book(TTree* tree)
   {
      // books branches of the new tree
      b_event = tree->Branch("event", &event, "event/I");
      b_year = tree->Branch("year", &year, "year/I");
      b_month = tree->Branch("month", &month, "month/I");
      b_day = tree->Branch("day", &day, "day/I");
      b_hour = tree->Branch("hour", &hour, "hour/I");
      b_minute = tree->Branch("minute", &minute, "minute/I");
      b_second = tree->Branch("second", &second, "second/I");
      b_millisecond = tree->Branch("millisecond", &millisecond, "millisecond/I");
      b_tc1 = tree->Branch("tc1", &tc1, "tc1/I");
      b_b1_t = tree->Branch("b1_t", &b1_t, "b1_t[1024]/F");
      // channel in use: from 1 to 4 (its number), not in use: 0
      b_usedchan = tree->Branch("usedchan", &usedchan, "usedchan[4]/I");
      // book used channels only
      if (usedchan[0] > 0) b_b1_c1 = tree->Branch("b1_c1", &b1_c1, "b1_c1[1024]/F");
      if (usedchan[1] > 0) b_b1_c2 = tree->Branch("b1_c2", &b1_c2, "b1_c2[1024]/F");
      if (usedchan[2] > 0) b_b1_c3 = tree->Branch("b1_c3", &b1_c3, "b1_c3[1024]/F");
      if (usedchan[3] > 0) b_b1_c4 = tree->Branch("b1_c4", &b1_c4, "b1_c4[1024]/F");
   }
   void Connect(TTree* tree)
   {
      clear();
      //SaveBranchAddresses(tree);
      // connects tree buffers with variables to use for event-by-event analysis
      // rest of the event
      cout<< "&event = " << &event <<endl;
      if ((b_event = tree->GetBranch("event"))) b_event->SetAddress(&event);
      if ((b_year = tree->GetBranch("year"))) b_year->SetAddress(&year);
      if ((b_month = tree->GetBranch("month"))) b_month->SetAddress(&month);
      if ((b_day = tree->GetBranch("day"))) b_day->SetAddress(&day);
      if ((b_hour = tree->GetBranch("hour"))) b_hour->SetAddress(&hour);
      if ((b_minute = tree->GetBranch("minute"))) b_minute->SetAddress(&minute);
      if ((b_second = tree->GetBranch("second"))) b_second->SetAddress(&second);
      if ((b_millisecond = tree->GetBranch("millisecond"))) b_millisecond->SetAddress(&millisecond);
      if ((b_tc1 = tree->GetBranch("tc1"))) b_tc1->SetAddress(&tc1);
      if ((b_b1_t = tree->GetBranch("b1_t"))) b_b1_t->SetAddress(&b1_t);
      // channel in use: from 1 to 4 (its number), not in use: 0
      if ((b_usedchan = tree->GetBranch("usedchan"))) b_usedchan->SetAddress(&usedchan);
      // connect used channels only
      if ((b_b1_c1 = tree->GetBranch("b1_c1"))) b_b1_c1->SetAddress(&b1_c1);
      if ((b_b1_c2 = tree->GetBranch("b1_c2"))) b_b1_c2->SetAddress(&b1_c2);
      if ((b_b1_c3 = tree->GetBranch("b1_c3"))) b_b1_c3->SetAddress(&b1_c3);
      if ((b_b1_c4 = tree->GetBranch("b1_c4"))) b_b1_c4->SetAddress(&b1_c4);
      if ((b_tc1 = tree->GetBranch("tc1"))) b_tc1->SetAddress(&tc1);
      if ((b_b2_t = tree->GetBranch("b2_t"))) b_b2_t->SetAddress(&b2_t);
      // channel in use: from 1 to 4 (its number), not in use: 0
      if ((b_usedchan = tree->GetBranch("usedchan"))) b_usedchan->SetAddress(&usedchan);
      // connect used channels only
      if ((b_b2_c1 = tree->GetBranch("b2_c1"))) b_b2_c1->SetAddress(&b2_c1);
      if ((b_b2_c2 = tree->GetBranch("b2_c2"))) b_b2_c2->SetAddress(&b2_c2);
      if ((b_b2_c3 = tree->GetBranch("b2_c3"))) b_b2_c3->SetAddress(&b2_c3);
      if ((b_b2_c4 = tree->GetBranch("b2_c4"))) b_b2_c4->SetAddress(&b2_c4);
   }
   void Disconnect(TTree* tree) {
      tree->ResetBranchAddresses();
   }
   TGraph* plot(Int_t ch=0) const      // NB: init to non-existant channel
   {
      if (ch == 0) {
         cout<< "Usage: plot(channel)" <<endl;
         UsedChan();
         return 0;
      }

      const Int_t color[] = {2, 4, 6, 8};

      const Float_t* v = 0;
      if (ch == 1 and usedchan[0]) v = b1_c1;
      if (ch == 2 and usedchan[1]) v = b1_c2;
      if (ch == 3 and usedchan[2]) v = b1_c3;
      if (ch == 4 and usedchan[3]) v = b1_c4;
      if (!v) {
         cout<< "Channel not found: " << ch <<endl;
         UsedChan();
         return 0;
      }
      TGraph* gr = new TGraph(1024,b1_t,v);
      gr->SetNameTitle(Form("gr_evt_%d_chan_%d",event,ch), Form("gr_evt_%d_chan_%d;time, ns;V",event,ch));
      gr->SetMarkerStyle(6);
      gr->SetMarkerColor(color[ch]);
      gr->Draw("ap");
      return gr;
   }
};
// //--------------------- test function ------------------------
// void pulsebuffer(TTree* tree) {
//    PulseBuffer buf(tree);
// }

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

//--------------------------------------------------------------------------------

TGraphErrors* ftplot(Int_t evt, Int_t chan, TTree* tree_ft=0, TTree* tree_pulse=0, bool allchan=false, bool samecan=false)
{
   // bool implicit_ft = tree_ft == 0;
   // bool implicit_pulse = tree_pulse == 0;
   // if (implicit_ft or implicit_pulse) cout<< "--> ftplot: tree_ft = " << tree_ft << " tree_pulse = " << tree_pulse <<endl;

   if (!tree_pulse) tree_pulse = (TTree*) gDirectory->Get("pulse");
   if (!tree_pulse)
   {
      // tree tree_pulse may be friend of the tree ft
      // cout<< "find tree ft in search for tree_pulse" <<endl;
      if (!tree_ft) {
         tree_ft = (TTree*) gDirectory->Get("ft");   // look for the tree ft
         if (!tree_ft) {
            cout<< "***Error ftplot: no trees" <<endl;
            return 0;
         }
      }
      // cout<< "trying to find tree_pulse as a friend of the tree_ft" <<endl;
      tree_pulse = tree_ft->GetFriend("pulse");
      if (!tree_pulse) {
         tree_pulse = tree_ft->GetFriend("p");
         if (!tree_pulse) {
            cout<< "***Warning ftplot: tree \"pulse\" was found neither in the current dir not as a friend" <<endl;
            return 0;
         }
      }
   }
   // we have tree tree_pulse at this point (both trees if we got lucky)

   // if (implicit_ft or implicit_pulse) cout<< "at this point: ftplot: tree_ft = " << tree_ft << " tree_pulse = " << tree_pulse <<endl;

   if (!tree_ft) tree_ft = (TTree*) gDirectory->Get("ft");
   if (!tree_ft) {
      // tree ft may be friend of the tree pulse
      tree_ft = tree_pulse->GetFriend("ft");
      // if (!tree_ft) {
      //    cout<< "***Warning ftplot: no friend tree \"ft\" has been found" <<endl;
      // }
   }

   // cout<< "before the start: ftplot: tree_ft = " << tree_ft << " tree_pulse = " << tree_pulse <<endl;
   // cout<< "Processing tree_ft->GetName() = " << tree_ft->GetName() << " and tree_pulse->GetName() = " << tree_pulse->GetName() <<endl;

   //-- const Int_t nsample = 1024;

   Int_t board = (chan-1)/4 + 1;
   Int_t channel = (chan-1) % 4 + 1;
   // cout<< "--> board = " << board << " channel = " << channel <<endl;

   std::string tleaf = Form("t%d",board);
   std::string vleaf = Form("c%d",chan);
   if (tree_pulse->GetLeaf(tleaf.c_str()))
   {
      // found time brach (new name), look for the voltage branch
      if (!tree_pulse->GetLeaf(vleaf.c_str())) {
         cout<< "Could not find voltage branch" <<endl;
         return 0;
      }
   }
   else {
      // look for the old branch names
      Int_t chan_1_4 = (chan-1) % 4 + 1;
      tleaf = Form("b%d_t", board);
      vleaf = Form("b%d_c%d",board,chan_1_4);
      if (!tree_pulse->GetLeaf(tleaf.c_str())) {
         cout<< "Could not find time branch" <<endl;
         return 0;
      }
      if (!tree_pulse->GetLeaf(tleaf.c_str())) {
         cout<< "Could not find voltage branch" <<endl;
         return 0;
      }
   }

   // std::string lname = Form("b%d_c%d",board,channel);
   // if (!tree_pulse->GetLeaf(lname.c_str())) {
   //    cout<< "No such leaf: " << lname.c_str() <<endl;
   //    cout<< "TTree " << tree_pulse->GetName() << " has leaves: " <<endl;
   //    for (int i=0; i<tree_pulse->GetListOfLeaves()->GetEntries(); ++i) {
   //       TLeaf* leaf = (TLeaf*) tree_pulse->GetListOfLeaves()->At(i);
   //       cout<< i << "\t " << leaf->GetName() <<endl;
   //    }
   //    return 0;
   // }
   Float_t* x = (Float_t*) tree_pulse->GetLeaf(tleaf.c_str())->GetValuePointer();
   Float_t* y = (Float_t*) tree_pulse->GetLeaf(vleaf.c_str())->GetValuePointer();
   //-- Int_t* evt_pulse_ptr = (Int_t*) tree_pulse->GetLeaf("event")->GetValuePointer();

   // fit tree

   OscFit* oscFit = 0;
   if (tree_ft) tree_ft->SetBranchAddress("oscFit", &oscFit);

   // plot

   tree_pulse->LoadTree(evt);
   tree_pulse->GetEntry(evt);
   if (tree_ft) {
      tree_ft->LoadTree(evt);
      tree_ft->GetEntry(evt);
      //cout<< "oscFit->evt = " << oscFit->evt << " oscFit->xt[0] = " << oscFit->xt[0] <<endl;
      //-- if (*evt_pulse_ptr != oscFit->evt) {
      //--       cout<< "Event mismatch:\nevt_pulse = " << *evt_pulse_ptr << "\noscFit->evt = " << oscFit->evt <<endl;
      //--       return 0;
      //-- }
   }

   // cout<< "--> evt = " << evt <<endl;
   // cout<< "sigmin = " << oscFit->sigmin[channel-1] <<endl;
   // cout<< "sigmax = " << oscFit->sigmax[channel-1] <<endl;
   // cout<< "ptx = " << oscFit->ptx[channel-1] <<endl;
   // cout<< "pty = " << oscFit->pty[channel-1] <<endl;
   // cout<< "lmaxx = " << oscFit->lmaxx[channel-1] <<endl;
   // cout<< "lmaxy = " << oscFit->lmaxy[channel-1] <<endl;

   Float_t gr_y[1024], gr_ey[1024];
   for (int i=0; i<1024; ++i) {
      Float_t bkg = tree_ft? oscFit->bkg[channel-1]: 0;
      gr_y[i] = -y[i] - bkg;
      //-------- gr_y[i] = y[i] - bkg;
      Float_t sbkg = tree_ft? oscFit->sbkg[channel-1]: 0;
      gr_ey[i] = sbkg;
   }

   // find limits for the plot
   Float_t sigmin = tree_ft? oscFit->sigmin[channel-1]: x[0];
   Float_t sigmax = tree_ft? oscFit->sigmax[channel-1]: x[1023];
   Int_t isigmin = 0;
   Int_t isigmax = 0;
   const Float_t eps = 1e-7;
   for (int i=0; i<1024; ++i)
   {
      if (TMath::Abs(x[i] - sigmin) < eps) isigmin = i;
      if (TMath::Abs(x[i] - sigmax) < eps) {
         isigmax = i;
         break;
      }
   }

   Int_t ibkgmin = isigmin > 40? isigmin - 40: 0;     // was used for debug plot
   Int_t np = isigmax-ibkgmin+1;
   if (allchan)
   {
      ibkgmin = 0;
      np = 1024;
   }

   Int_t markerStyle = tree_ft? 24: 7;   // set different marker style w/ and w/o tree ft
   Int_t markerColor = tree_ft? 8: 2;    // set different marker color w/ and w/o tree ft

   TGraphErrors* gr = new TGraphErrors(np, &x[ibkgmin], &gr_y[ibkgmin], 0, &gr_ey[ibkgmin]);
   // gr->SetNameTitle(Form("gr_evt%d_ch%d",evt,channel), Form("gr_evt%d_ch%d I=%f",evt,channel,oscFit->pti[channel-1]));
   gr->SetNameTitle(Form("gr_evt%d_ch%d",evt,channel), Form("gr_evt%d_ch%d I=%f Itot=%f",evt,channel,oscFit->pti[channel-1],oscFit->itot[channel-1]));
   gr->SetMarkerStyle(markerStyle);
   gr->SetMarkerColor(markerColor);
   gr->SetLineColor(markerColor);
   // fit function
   if (tree_ft) {
      TF1* fpl = new TF1("fpl", "pol1", oscFit->plxmin[channel-1],oscFit->plxmax[channel-1]);
      // fpl->SetNpx(1024);
      fpl->SetLineWidth(1);
      fpl->SetLineColor(2);
      fpl->SetParameter(0, oscFit->plp0[channel-1]);
      fpl->SetParameter(1, oscFit->plp1[channel-1]);
      Int_t ndf = 4 - 2;
      fpl->SetNDF(ndf);
      fpl->SetChisquare(oscFit->plc[channel-1] * ndf);
      gr->GetListOfFunctions()->Add(fpl);
   }

   if (!samecan) {
      if (gROOT->GetListOfCanvases()->GetEntries() < 2) new TCanvas;
      else {
         std::string cname = nextcan(Form("chan_%d_event_%d", channel,evt));
         new TCanvas(cname.c_str(),cname.c_str(), gStyle->GetCanvasDefX(), gStyle->GetCanvasDefY(), gStyle->GetCanvasDefW(), gStyle->GetCanvasDefH());
         //cout<< "created new TCanvas " << gPad->GetName() <<endl;
      }
   }

   if (gr->GetN() > 0) gr->Draw("ap");
   else {
      cout<< gr->GetTitle() << ": n points == 0" <<endl;
      gPad->DrawFrame(0,0,1,1, gr->GetTitle());
   }
   if (tree_ft) {
      TMarker* marker_pt = new TMarker(oscFit->ptx[channel-1], oscFit->pty[channel-1], 5);
      marker_pt->SetMarkerColor(2);
      marker_pt->SetMarkerSize(2);
      marker_pt->Draw("same");

      TMarker* marker_local_maximum = new TMarker(oscFit->lmaxx[channel-1], oscFit->lmaxy[channel-1], 20);
      marker_local_maximum->SetMarkerColor(4);
      marker_local_maximum->SetMarkerSize(1);
      marker_local_maximum->Draw("same");
   }

   tree_pulse->ResetBranchAddresses();
   tree_ft->ResetBranchAddresses();
   return gr;
}

/*
root -l run_1_009.root run_1_009.root-ft.root 
.L point.C+
ft->Draw(">>elist","plx[7]>50&&plx[7]<60 &&sat[7] &&plc[7]<200")
elist=elist
ft=ft
cd("run_1_009.root")
pulse=pulse
ftplot(elist,7,pulse,ft)
// NB event 1014 with 3 protons
// NB2: fit for the ch8 (evt 1014) seems weird on ftplot because the first two points are spikes
ftplot(elist,8,pulse,ft,true)
*/

/*
--> the cuts.C should contain at least this definition:
TCut run35[8] = {
   "ptok[0]&&ptx[0]>50&&ptx[0]<60",
   "ptok[1]&&ptx[1]>50&&ptx[1]<60",
   "ptok[2]&&ptx[2]>50&&ptx[2]<60",
   "ptok[3]&&ptx[3]>40&&ptx[3]<60 &&pti[3]>1&&pti[3]<2",
   "ptok[4]&&ptx[4]>40&&ptx[4]<60",
   "ptok[5]&&ptx[5]>45&&ptx[5]<55 &&pti[5]>7&&pti[5]<14",
   "ptok[6]&&ptx[6]>45&&ptx[6]<55 &&pti[6]>13&&pti[6]<20",
   "ptok[7]&&ptx[7]>40&&ptx[7]<60 &&pti[7]>1&&pti[7]<2"
};

root -l tb2012_run_34.root.pulse.root-ft_10.000.root
.L point.C+
.x cuts.C
ft->AddFriend("pulse","tb2012_run_34.root.pulse.root")
ft->Draw(">>elist",run35[0]+run35[3]+"pti[0]>2.5 &&plx[0]-plx[3]>2")
elist=elist
ftplot(elist,4)
// ftplot(elist,1)
*/
//--old sequence of ft and pulse parameters: void ftplot(TEventList* eventlist, Int_t channel, TTree* tree_pulse=0, TTree* tree_ft=0, bool allchan=false)
void ftplot(TEventList* eventlist, Int_t channel, TTree* tree_ft=0, TTree* tree_pulse=0, bool allchan=false)
{
   TTimer timer("gSystem->ProcessEvents();",50,kFALSE);  // process mouse events every 50 ms

   // read event number
   Int_t nentry_tlist = eventlist->GetN();
   for (Int_t i=0; i<nentry_tlist; ++i)
   {
      if (i>0 && i % 10 == 0)
      {
         cout<< "<CR> = next 10 plots. Q = Quit: ";
         timer.Start();    // start processing of mouse events on TCanvas

         std::string line;
         std::getline(std::cin, line);

         timer.Stop();     // disable timer after <CR>

         if (line.size() == 0) continue;
         if (line.find("Q") != std::string::npos || line.find("q") != std::string::npos) break;

         timer.Stop();     // disable timer after <CR>
      }

      Int_t jentry = eventlist->GetEntry(i);
      // cout<< "--> jentry = " << jentry <<endl;
      ftplot(jentry, channel, tree_pulse, tree_ft, allchan, false);
   }
}

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

/// /*
/// good examples:
/// run_1_009.root, entry 1817, two particles, ~17ns: b1_c4:b1_t clean, b2_c2_b2_t noise, b2_c3:b2_t noC
/// */
/// TH1F* hevt(TTree* tree, Int_t chan, Int_t entry, Int_t nbkg_approx=20)
/// {
///    PulseBuffer buf(tree);
/// 
///    tree->LoadTree(entry);
///    tree->GetEntry(entry);
/// 
///    Float_t* t = 0;
///    Float_t* v = 0;
///    if (chan <= 4) {
///       t = buf.b1_t;
///       switch (chan) {
///          case 1: v = buf.b1_c1; break;
///          case 2: v = buf.b1_c2; break;
///          case 3: v = buf.b1_c3; break;
///          case 4: v = buf.b1_c4; break;
///       }
///    }
///    else {
///       t = buf.b2_t;
///       switch (chan % 4) {
///          case 1: v = buf.b2_c1; break;
///          case 2: v = buf.b2_c2; break;
///          case 3: v = buf.b2_c3; break;
///          case 4: v = buf.b2_c4; break;
///       }
///    }
///    if (v == 0) {
///       cout<< "Data was not found" <<endl;
///       buf.Disconnect(tree);
///       return 0;
///    }
/// 
///    // sort bkg array to eliminate possible spikes: three highest points
///    Int_t index[1024];
///    TMath::Sort(nbkg_approx, v, index, kFALSE);   // sort in the increasing order
/// 
///    Float_t bkg_mean_approx = 0;
///    // omit the first and the last two values
///    for (int i=2; i<nbkg_approx-2; ++i) bkg_mean_approx += v[index[i]];  // omit first 2 and last 2
///    bkg_mean_approx /= (nbkg_approx-2-2);
/// 
///    // noise RMS
///    Float_t sumx2 = 0;
///    for (int i=2; i<nbkg_approx-2; ++i) {  // omit first 2 and last 2
///       Float_t dy = v[index[i]] - bkg_mean_approx;
///       sumx2 += dy*dy;
///    }
///    Float_t bkg_sigma_approx = TMath::Sqrt(sumx2 / (nbkg_approx-2-2));
///    cout<< "approx value of the bkg is " << bkg_mean_approx << " bkg_sigma_approx = " << bkg_sigma_approx <<endl;
/// 
///    Int_t nsample = 1024;
/// 
///    // find the pulse sign: the negative pulses will be inverted
///    Double_t y_average = 0;
///    for (int i=0; i<nsample; ++i) {
///       y_average += v[i];
///    }
///    y_average = y_average/nsample - bkg_mean_approx;
///    Double_t ysign = y_average < 0? -1.: 1.;
///    cout<< "y_average = " << y_average << " ysign = " << ysign <<endl;
/// 
///    Float_t xlow = t[0];
///    Float_t dt = t[1023] - t[1023-1];
///    Float_t xup = t[1023] + dt;
///    TH1F* hchan = new TH1F(nextname(Form("ch_%d_evt_%d",chan,entry)),Form("ch %d evt %d",chan,entry), 1024,xlow,xup);
///    for (int i=0; i<1024; ++i) {
///       hchan->SetBinContent(i+1, ysign*(v[i]-bkg_mean_approx));
///       hchan->SetBinError(i+1, bkg_sigma_approx);
///    }
/// 
///    buf.Disconnect(tree);
///    return hchan;
/// }

std::ostream& endn(std::ostream& os)
{
   // example of usage:
   //
   // bool debug = true;
   // std::ostream& (*end)(std::ostream&);   // end is a pointer to a function
   // if (debug) end = std::endl;
   // else end = endn;

   os << "\n";
   return os;
}

void point(const Float_t xraw[], const Float_t yraw[]
      , OscFit* oscFit
      , Double_t thres
      , Bool_t debug = kFALSE
      , Int_t evt = 0
      , Int_t channel = 0        // actual channel number starts from 1
      , Float_t trig = 0
      , Float_t gate = 0
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

   std::ostream& (*end)(std::ostream&);   // end is a pointer to the function
   if (debug) end = std::endl;
   else end = endn;

   cout<< "------> evt = " << evt << " channel = " << channel <<end;

   Int_t nsample = 1024;

   const Int_t nbkg_approx = 20;

   Int_t istart = nbkg_approx - 1;                 // use the first nbkg_approx points for bkg only
   if (trig < xraw[istart]) {
      trig = xraw[istart];                         // adjust the trig if needed
      if (debug) cout<< "***Warning: point: trig was set to " << trig <<endl;
   }
   for (int i=istart; i<nsample; ++i) {
      if (xraw[i] < trig) istart = i;
      else break;
   }
   if (gate == 0) gate = xraw[nsample-1] - trig;
   Int_t iend = 0;
   for (int i=istart+1; i<nsample; ++i) {
      if (xraw[i] > trig + gate) break;
      iend = i;
   }
   if (debug) cout<< "istart = " << istart << " xraw[istart] = " << xraw[istart] << " iend = " << iend << " xraw[iend] = " << xraw[iend] <<endl;

   /// // find the pulse sign: the negative pulses will be inverted
   /// Double_t y_average = 0;
   /// Double_t y_background = 0;
   /// for (int i=0; i<nsample; ++i) {
   ///    y_average += yraw[i];
   ///    if (i < 20) y_background = y_average;
   /// }
   /// y_average /= nsample;
   /// y_background /= 20;
   /// Double_t ysign = y_average - y_background < 0? -1.: 1.;
   /// if (debug) cout<< "y_average - y_background = " << y_average - y_background <<end;

   oscFit->thres[channel-1] = thres;
   oscFit->trig[channel-1] = trig;
   oscFit->gate[channel-1] = gate;
   oscFit->trig[channel-1] = 0;
   oscFit->gate[channel-1] = 0;
   for (int j=0; j<30; ++j) oscFit->tq[channel-1][j] = 0;
   // assign output values
   oscFit->bkg[channel-1] = 0;
   oscFit->sbkg[channel-1] = 0;
   // algorithms
   oscFit->yscal[channel-1] = 0;
   oscFit->sigmin[channel-1] = 0;
   oscFit->sigmax[channel-1] = 0;
   oscFit->lmaxx[channel-1] = 0;
   oscFit->lmaxy[channel-1] = 0;
   oscFit->pti[channel-1] = 0;
   oscFit->sat[channel-1] = 0;
   oscFit->pty[channel-1] = 0;
   oscFit->ptx[channel-1] = 0;
   oscFit->plx[channel-1] = 0;
   oscFit->plx0[channel-1] = 0;
   oscFit->plc[channel-1] = 0;
   // save fit in the tree
   oscFit->plxmin[channel-1] = 0;
   oscFit->plxmax[channel-1] = 0;
   oscFit->plp0[channel-1] = 0;
   oscFit->plp1[channel-1] = 0;

   oscFit->ptok[channel-1] = kFALSE;
   oscFit->itot[channel-1] = 0;
   oscFit->ibkg[channel-1] = 0;
   oscFit->itail[channel-1] = 0;

   const Float_t* x = xraw;
   Float_t y[1024], ey[1024];       // NB: y gives positive signal: negative signal will be inverted

   // histogram for bkg
   static TH1F *hbkg;
   if (!hbkg) {
      //-- hbkg = new TH1F("hbkg_partime", "hbkg_partime", 400, -100, 100);
      // hbkg = new TH1F("hbkg_partime", "hbkg_partime", 400, -10, 10);
      hbkg = new TH1F("hbkg_partime", "hbkg_partime", 400, -0.01, 0.01);
   }
   hbkg->SetDirectory(0);

   ///////////////////////////////////////// pulse function ///////////////////////////////////////////////

   Double_t A = 1000;
   Double_t x0 = 50;
   Double_t tau1 = 0.1;
   Double_t tau2 = 1;
   Double_t T = 0;
   Double_t sigma = 0.1;

   Double_t par[10];
   Int_t npar = 0;

   enum {ipar_A=0, ipar_x0, ipar_tau1, ipar_tau2, ipar_T, ipar_sigma};

   par[npar++] = A;
   par[npar++] = x0;
   par[npar++] = tau1;
   par[npar++] = tau2;
   par[npar++] = T;
   par[npar++] = sigma;
   //if (debug) cout<< "\t npar = " << npar <<endl;

   // function for fit
   TF1* fun_Psigma = new TF1("fun_Psigma", fPsigma, 0, 500, npar);
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

   // function for refit
   TF1* fun_Psigma_refit = new TF1("fun_Psigma_refit", fPsigma, 0, 500, npar);
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

   ////////////////////////////////////////////////////////////////////////////////////////////////////////

   // find approx value of the bkg

   // DRS4 version 3 sometime has
   // (formerly known as USB) spikes: 2 spikes of two neiboring channels, low in yraw
   // sometime I observed single-channel spike (channel 2): high in yraw

   //-- const Int_t nbkg_approx = 20;

   // sort bkg array to eliminate possible spikes: three highest points
   Int_t index[1024];
   TMath::Sort(nbkg_approx, yraw, index, kFALSE);   // sort in the increasing order: spikes are high in yraw

   // look at the lowest of the data point: this can be single spike (low in yraw)
   const Double_t thres_single_spike = -100.;                              // NB high threshold value!
   Int_t single_spike_i = yraw[index[0]] < thres_single_spike? index[0]: -1;

   Float_t bkg_mean_approx = 0;
   //Float_t bkg_sigma_approx = 0;
   for (int i=1; i<nbkg_approx-4; ++i) bkg_mean_approx += yraw[index[i]];  // omit first 1 and last 4
   bkg_mean_approx /= (nbkg_approx-1-4);
   bkg_mean_approx /= -1.;                          // invert bkg_mean built out of yraw
   cout<< "approx value of the bkg is " << bkg_mean_approx <<end;

   // integral over the total pulse with approximate background
   Float_t itot = 0;
   Int_t ntot = 0;
   for (int is=20; is<1020; ++is) {
      itot += (-1.*0.5*(yraw[is] + yraw[is-1]) - bkg_mean_approx) * (x[is] - x[is-1]);
      ++ntot;
   }
   // itot *= x[21] - x[20];
   // normalize integral to pC by division to R = 50 Ohm
   // //--50 Ohm: itot /= 50.;
   // oscFit->itot[channel-1] = itot;

   std::vector<TMarker*> marker_spikes;

   // in the following loop over all the points:
   //
   // 1) convert yraw to y using approx value of the bkg
   // 2) correct spikes (if any)
   // 3) find the maximum (consider data points over the thres only)
   // 4) find the right peak boundaty
   //    -- 1/8 of the maximum
   //    -- process the first peak only (in case of event with multiple peaks)

   Bool_t saturation = kFALSE;
   for (int i=0; i<nsample; ++i)
   {
      if (yraw[i] < -499.) saturation = kTRUE;

      y[i] = -1.*yraw[i] - bkg_mean_approx;           // will use y instead of yraw from this point

      bool remove_spikes = false;                     //-- do not bother with spikes
      if (remove_spikes and i >= 3)
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
                  if (debug) cout<< "... corrected spikes at x[" << i-2 << "] = " << x[i-2] <<end;
               }
            }
         }
      }
   }

   Int_t maximum_i = nbkg_approx;
   Float_t maximum_x = x[maximum_i];
   Float_t maximum_y = yraw[maximum_i];

   bool passed_thres_rise = false;     // exceeded thres at the rising (left) slope
   bool passed_peak = false;           // passed the peak

   Double_t eof_peak = y[maximum_i]/8.;
   for (int i=istart; i<=iend; ++i)
   {
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
         if (debug) cout<< "... passed the peak: x[" << i << "] " << x[i] << " y " << y[i] << " < eof_peak " << eof_peak <<end;
         passed_peak = true;
      }
   }
   if (!passed_thres_rise) {
      // store in the tree itot with approximate background
      cout<< "point: no pulse found" <<end;
      //--> already done: itot *= x[21] - x[20];
      // normalize integral to pC by division to R = 50 Ohm
      //--50 Ohm: itot /= 50.;
      cout<< "point: store itot (with approximate background) " << itot <<end;
      oscFit->itot[channel-1] = itot;
      return;
   }
   else cout<< "maximum_x " << maximum_x << " maximum_y " << maximum_y << " maximum_i " << maximum_i <<end;

   // find width at the half maximum

   Double_t halfmax = maximum_y / 2.;
   Double_t quotermax = maximum_y / 4.;
   Int_t halfmax_x1_i = 0;
   Int_t halfmax_x2_i = 0;

   // if (istart < nbkg_approx) istart = nbkg_approx;    //----- TODO

   for (int i=istart; i<=iend; ++i)
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
   if (debug) cout<< "left peak slope: " << x[halfmax_x1_i] << " right peak slope: " << x[halfmax_x2_i] <<end;;

   // three sigma to the left
   // Exactly: sigma = FWHM/2.35 = delta_x * 2*(maximum_i - halfmax_x1_i + 1) / 2.35
   // Approx: sigma = delta_x * (maximum_i - halfmax_x1_i + 1)

   Int_t i3sigma1 = 3 * (maximum_i - halfmax_x1_i);   // 3 sigma1 to the left
   //---- Int_t i3sigma1 = 5 * (maximum_i - halfmax_x1_i);   // 3 sigma1 to the left
   Int_t i3sigma2 = 3 * (halfmax_x2_i - maximum_i);   // 3 sigma2 to the right
   //---- Int_t i3sigma2 = 5 * (halfmax_x2_i - maximum_i);   // 3 sigma2 to the right

   if (i3sigma1 < 12) i3sigma1 = 12;
   if (i3sigma2 < 12) i3sigma2 = 12;

   Int_t isigmin = maximum_i - i3sigma1;
   //-- if (isigmin < 10) isigmin = 10;               // isigmin >= 10
   if (isigmin < istart) isigmin = istart;         // isigmin >= istart
   Int_t isigmax = maximum_i + i3sigma2;
   //-- if (isigmax > 1023) isigmax = 1023;
   if (isigmax > iend) isigmax = iend;
   Int_t ibkgmax = isigmin - 1;
   if (ibkgmax < nbkg_approx) ibkgmax = nbkg_approx;
   Int_t ibkgmin = isigmin > 40? isigmin - 40: 0;     // will be used for the plot
   if (debug) cout<< "isigmin = " << isigmin << " isigmax = " << isigmax << " ibkgmin = " << ibkgmin << " ibkgmax = " << ibkgmax <<endl;
   if (debug) cout<< "i3sigma1 " << i3sigma1 << " i3sigma2 " << i3sigma2 << " x[isigmin] " << x[isigmin] << " x[isigmax] " << x[isigmax] <<end;
   if (debug) cout<< "left edge: " << x[isigmin] << " right edge: " << x[isigmax] <<end;

   // find the exact bkg
   //-- standard background --// Int_t nbkg = ibkgmax;
   Int_t nbkg = ibkgmax-ibkgmin+1;
   // /* NB: we already corrected double spikes (formerly known as USB) */ TMath::Sort(nbkg, y, index, kTRUE); // sort in the decreasing order: spikes are high in yraw
   //-- standard background range --// TMath::Sort(nbkg, y, index, kFALSE); // sort in the increasing order: the single spike is low in yraw
   TMath::Sort(nbkg, &y[ibkgmin], index, kFALSE); // sort in the increasing order: the single spike is low in yraw
   //--------- if (debug) for (int i=0; i<nbkg; ++i) {cout<< i << "\t x[ibkgmin+index[i]] = " << x[ibkgmin+index[i]] << " y[ibkgmin+index[i]] = " << y[ibkgmin+index[i]] <<endl;}
   /// hbkg->GetXaxis()->SetLimits(y[ibkgmin+index[1]]-0.001, y[ibkgmin+index[nbkg-1-1]]+0.001);
   /// hbkg->Reset();
   /// if (debug) cout<< "y[ibkgmin+index[1]]-0.001 = " << y[ibkgmin+index[1]]-0.001 << " y[ibkgmin+index[nbkg-1-1]]+0.001 = " << y[ibkgmin+index[nbkg-1-1]]+0.001 << " hbkg->GetXaxis()->GetXmin() = " << hbkg->GetXaxis()->GetXmin() << " hbkg->GetXaxis()->GetXmax() = " << hbkg->GetXaxis()->GetXmax() <<endl;
   /// for (int i=1; i<nbkg-1; ++i) hbkg->Fill(y[ibkgmin+index[i]]);    // omit last 1 -- protection against the single spike
   /// //-------- new TCanvas;        hbkg->DrawCopy();
   /// hbkg->Fit("gaus", "L0Q", "goff");                        // NB: Loglikelihood option
   /// if (debug) {
   ///    new TCanvas;
   ///    hbkg->DrawCopy();
   /// }

   /// Double_t bkg_mean_corr = hbkg->GetFunction("gaus")->GetParameter("Mean");
   /// Double_t bkg_sigma = hbkg->GetFunction("gaus")->GetParameter("Sigma");
   /// Double_t bkg_mean = bkg_mean_approx + bkg_mean_corr;        // to be stored in the tree

   Double_t bkg_sum = 0;
   Int_t bkg_n = 0;
   for (int i=ibkgmin; i<=ibkgmax; ++i) {
      bkg_sum += y[i];
      ++bkg_n;
   }
   Double_t bkg_mean_corr = bkg_sum / bkg_n;
   Double_t bkg_sigma = 0;
   Double_t bkg_mean = bkg_mean_approx + bkg_mean_corr;        // to be stored in the tree
   if (debug) cout<< "x[ibkgmin] = " << x[ibkgmin] << " x[ibkgmax] = " << x[ibkgmax] << " nbkg = " << nbkg <<endl;
   cout<< "bkg_mean_corr = " << bkg_mean_corr << " bkg_mean = " << bkg_mean << " bkg_sigma = " << bkg_sigma <<end;

   //if (debug) {new TCanvas; hbkg->DrawCopy();}

   // // correct the data with the updated bkg level
   // Float_t integral = 0;
   // Float_t yprev = y[isigmin] - bkg_mean_corr;
   // Float_t integral_bkg = 0;
   // Float_t integral_tail = 0;
   // for (int i=0; i<nsample; ++i) {
   //    y[i] -= bkg_mean_corr;
   //    ey[i] = bkg_sigma;
   //    if (i > isigmin && i <= isigmax) {
   //       integral += 0.5*(yprev + y[i]);
   //    }
   //    if (i > 0 && i < isigmin) {
   //       integral_bkg += 0.5*(yprev + y[i]);
   //    }
   //    if (i > isigmin) {
   //       integral_tail += 0.5*(yprev + y[i]);
   //    }
   //    yprev = y[i];              // update yprev
   // }
   // maximum_y -= bkg_mean_corr;
   // // complete integral
   // integral *= x[1] - x[0];
   // integral_bkg *= x[1] - x[0];
   // integral_tail *= x[1] - x[0];
   // // normalize integral to pC by division to R = 50 Ohm
   // integral /= 50.;
   // integral_bkg /= 50.;
   // integral_tail /= 50.;
   // // integral over the all samples
   // itot -= ntot*bkg_mean_corr;
   // //--> already done: itot *= x[21] - x[20];
   // // normalize integral to pC by division to R = 50 Ohm
   // itot /= 50.;
   // cout<< "pulse integral = " << integral << " pC   itot = " << itot  << " pC" <<end;

   // correct the data with the updated bkg level
   Float_t integral = 0;
   Float_t integral_bkg = 0;
   Float_t integral_tail = 0;
   Float_t pulse_maximum_y = 0;
   //for (int i=istart; i<=iend; ++i) {
   for (int i=0; i<nsample; ++i) {
      //-- what the heck? --// y[i] -= bkg_mean_corr;
      y[i] -= bkg_mean_corr;
      ey[i] = bkg_sigma;
      Float_t ai = i > 0? 0.5*(y[i] + y[i-1])*(x[i] - x[i-1]): 0;
      if (i > isigmin && i <= isigmax) {
         integral += ai;
         if (y[i] > pulse_maximum_y) pulse_maximum_y = y[i];
      }
      if (i < isigmin) {
         integral_bkg += ai;
      }
      if (i > isigmin) {
         integral_tail += ai;
      }
   }
   maximum_y -= bkg_mean_corr;
   // normalize integral to pC by division to R = 50 Ohm
   //--50 Ohm: integral /= 50.;
   //--50 Ohm: integral_bkg /= 50.;
   //--50 Ohm: integral_tail /= 50.;
   // integral over the all samples
   itot -= ntot*bkg_mean_corr;
   //--> already done: itot *= x[21] - x[20];
   // normalize integral to pC by division to R = 50 Ohm
   //--50 Ohm: itot /= 50.;
   cout<< "pulse integral = " << integral << " pC   itot = " << itot  << " pC" <<end;

   TGraphErrors* gr = 0;
   if (debug) {
      Int_t np = isigmax-ibkgmin+1;
      gr = new TGraphErrors(np, &x[ibkgmin], &y[ibkgmin], 0, &ey[ibkgmin]);
      // gr->SetNameTitle("gr", Form("evt %d, ch %d", evt, channel));
      gr->SetNameTitle("gr", Form("evt %d, ch %d I %0.3g I_{TOT} %0.3g", evt, channel, integral, itot));
      gr->SetMarkerStyle(24);
      gr->SetMarkerColor(8);
      gr->SetLineColor(8);
      new TCanvas;
      if (gr->GetN() > 0) gr->Draw("ap");
      else {
         cout<< gr->GetTitle() << ": n points == 0" <<endl;
         gPad->DrawFrame(0,0,1,1, gr->GetTitle());
      }
      // spikes
      for (unsigned imarker=0; imarker<marker_spikes.size(); ++imarker) marker_spikes[imarker]->Draw("same");
   }

   // time stamp

   // find the local maximum

   // Double_t point_first_y = maximum_y/4.;
   Double_t point_first_y = maximum_y/8.;

   Int_t point_first_i = maximum_i;
   do --point_first_i; while (point_first_i > isigmin && y[point_first_i] > point_first_y);

   if (debug) cout<< "point_first_i " << point_first_i << " x[point_first_i] " << x[point_first_i] <<end; 
   Int_t local_maximum_i = maximum_i;

   const Float_t eps = 1e-7;  // used for fit limits and against dividing by 0 at ptx

   for (int i=point_first_i; i<isigmax-6; ++i)
   {
      //-- Int_t ncheck = 5;
      Int_t ncheck = 10;
      if (i+ncheck == maximum_i) {
         if (debug) cout<< "... local_maximum_i == maximum_i" <<end;
         local_maximum_i = maximum_i;
         break;
      }

      Int_t nlow_max = 3;                    // the number of the next points lower than maximum candidate
      Double_t ycurr = y[i] + ey[i];         // current level to compare with
      Int_t nlow = 0;
      for (int ii=i+1; ii<i+ncheck; ++ii)
      {
         if (y[ii] < ycurr) ++nlow;
      }
      if (nlow > nlow_max) {
         // nlow_max points out of ncheck are lower: this is a local maximum
         local_maximum_i = i+1;
         // find a maximum in this range
         for (int ii=i+1+1; ii<i+ncheck; ++ii) {
            if (y[ii] > y[local_maximum_i]) {
               local_maximum_i = ii;
            }
         }
         break;
      }
   }
   Double_t local_maximum_x = x[local_maximum_i];
   Double_t local_maximum_y = y[local_maximum_i];

   if (local_maximum_i > maximum_i) {
      if (debug) cout<< "... ---> local_maximum_i > maximum_i. Set to maximum" <<end;
      local_maximum_i = maximum_i;
      local_maximum_x = maximum_x;
      local_maximum_y = maximum_y;
   }

   cout<< "local_maximum_x " << local_maximum_x << " local_maximum_y " << local_maximum_y <<end;

   if (debug) {
      TMarker* marker_local_maximum = new TMarker(local_maximum_x, local_maximum_y, 20);
      marker_local_maximum->SetMarkerColor(4);
      marker_local_maximum->SetMarkerSize(1);
      marker_local_maximum->Draw("same");
   }

   // find half-max
   Double_t pty = local_maximum_y / 2.;
   Int_t x1_i = local_maximum_i;
   do x1_i--; while (x1_i > isigmin && y[x1_i] > pty);
   if (x1_i > 1023 - 2) x1_i = 1023 - 2;                 // to avoid out of range in 4-point fit
   Int_t x2_i = x1_i + 1;
   if (debug) cout<< "x[x1_i] " << x[x1_i] << " x[x2_i] " << x[x2_i] <<end;
   Double_t delta_x = x[x2_i] - x[x1_i];
   Double_t delta_y = y[x2_i] - y[x1_i];           // NB: delta_y > 0 in regular case
   Double_t ptx = 0;
   if (delta_y > eps) ptx = x[x2_i] - delta_x * (y[x2_i] - pty) / delta_y;
   if (ptx < x[x1_i] || ptx > x[x2_i]) {
      cout<< "--> warning: evt " << evt << " channel " << channel << ": ptx is out of limits" <<end;
      ptx = 0;    // something wrong
   }

   // fit straight line to two points below and two points above

   Int_t npfit = 4;
   Double_t p0 = 0, p1 = 0;
   Double_t chi2ndf = pol1fast(x,y,ey, x1_i-1, npfit, p0, p1);
   Double_t plx = TMath::Abs(p1) > eps? (pty - p0)/p1: 0;
   Double_t plx0 = TMath::Abs(p1) > eps? -p0/p1: 0;
   // to store in the tree
   Double_t plxmin = x[x1_i-1];
   Double_t plxmax = x[x1_i-1 + npfit - 1];

   cout<< "ptx " << ptx << " plx = " << plx << " pty = " << pty <<end;

   // fit exp to the pulse back slope to find exect pulse width at the half of maximum

   Int_t back_half_above_i = local_maximum_i;
   for (int i=local_maximum_i+1; i<=isigmax; ++i) {
      if (y[i] < local_maximum_y / 2.) break;
      back_half_above_i = i;
   }
   Int_t back_half_below_i = isigmax;
   for (int i=isigmax; i>back_half_above_i; --i) {
      if (y[i] > local_maximum_y / 2.) break;
      back_half_below_i = i;
   }
   // Int_t back_half_below_i = back_half_above_i;
   // for (int i=back_half_above_i; i<=isigmax; ++i) {
   //    if (y[i] > local_maximum_y / 2.) break;
   //    back_half_below_i = i;
   // }
   if (debug) cout<< "x[back_half_above_i] = " << x[back_half_above_i] << " x[back_half_below_i] = " << x[back_half_below_i] <<endl;

   // fit pulse function

   Double_t tau1_est = p1 > 0? 0.5*(local_maximum_y / p1): 0;
   Double_t x0_est = ptx;

   for (int i=0; i<fun_Psigma->GetNpar(); ++i) fun_Psigma->ReleaseParameter(i);
   fun_Psigma->SetParameter(ipar_T, 0);
   fun_Psigma->SetParameter(ipar_A, local_maximum_y/tau1_est);
   fun_Psigma->SetParameter(ipar_tau1, tau1_est);
   fun_Psigma->SetParameter(ipar_tau2, tau1_est);
   fun_Psigma->SetParameter(ipar_sigma, tau1_est);
   fun_Psigma->SetParameter(ipar_x0, x0_est);

   fun_Psigma->FixParameter(ipar_T, fun_Psigma->GetParameter(ipar_T));
   fun_Psigma->FixParameter(ipar_x0, fun_Psigma->GetParameter(ipar_x0));
   fun_Psigma->FixParameter(ipar_tau1, fun_Psigma->GetParameter(ipar_tau1));
   fun_Psigma->FixParameter(ipar_sigma, fun_Psigma->GetParameter(ipar_sigma));

   if (!gr) {
      Int_t np = isigmax-ibkgmin+1;
      gr = new TGraphErrors(np, &x[ibkgmin], &y[ibkgmin], 0, &ey[ibkgmin]);
      // gr->SetNameTitle("gr", Form("evt %d, ch %d", evt, channel));
      gr->SetNameTitle("gr", Form("evt %d, ch %d I %0.3g I_{TOT} %0.3g", evt, channel, integral, itot));
      gr->SetMarkerStyle(24);
      gr->SetMarkerColor(8);
      gr->SetLineColor(8);
      if (gr->GetN() == 0) delete gr;
   }

   TF1* fitted_function_refit = 0;
   Double_t y0cross = 0;
   Double_t y05cross = 0;

   Double_t integral_pulse_fitted = 0;

   if (gr) {
      Double_t xmin_fit, xmax_fit;
      xmin_fit = x[isigmin] - eps;
      xmax_fit = x[local_maximum_i+1] + eps;

      // if (debug) gr->Fit(fun_Psigma, "+R", "", x[isigmin]-eps, local_maximum_x+eps);
      // else gr->Fit(fun_Psigma, "+R0", "goff", x[isigmin]-eps, local_maximum_x+eps);

      // if (debug) gr->Fit(fun_Psigma, "+R", "", x[isigmin]-eps, x[isigmax]+eps);
      // else gr->Fit(fun_Psigma, "+R0", "goff", x[isigmin]-eps, x[isigmax]+eps);

      if (debug) gr->Fit(fun_Psigma, "+R", "", xmin_fit, xmin_fit);
      else gr->Fit(fun_Psigma, "+R0", "goff", xmin_fit, xmin_fit);

      // fitted integral
      integral_pulse_fitted = gr->GetFunction(fun_Psigma->GetName())->Integral(x[isigmin],1024.);

      // refit
      fitted_function_refit = gr->GetFunction(fun_Psigma->GetName());
      fun_Psigma_refit->SetParameters(gr->GetFunction(fun_Psigma->GetName())->GetParameters());
      for (int i=0; i<fun_Psigma_refit->GetNpar(); ++i) fun_Psigma_refit->ReleaseParameter(i);

      fun_Psigma_refit->FixParameter(ipar_tau2, gr->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_tau2));
      fun_Psigma_refit->FixParameter(ipar_T, gr->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_T));
      //
      // fun_Psigma_refit->FixParameter(ipar_A, gr->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_A));
      //fun_Psigma_refit->FixParameter(ipar_sigma, gr->GetFunction(fun_Psigma->GetName())->GetParameter(ipar_sigma));

      gr->GetListOfFunctions()->Remove(gr->GetListOfFunctions()->FindObject(fun_Psigma->GetName()));

      Double_t xmin_refit, xmax_refit;
      xmin_refit = x[isigmin] - eps;
      // xmax_refit = x[local_maximum_i - 1] + eps;
      xmax_refit = x[local_maximum_i] + eps;

      // if (debug) gr->Fit(fun_Psigma_refit, "+R", "", x[isigmin]-eps, local_maximum_x+eps);
      // else gr->Fit(fun_Psigma_refit, "+R0", "goff", x[isigmin]-eps, local_maximum_x+eps);

      if (debug) gr->Fit(fun_Psigma_refit, "+R", "", xmin_refit, xmax_refit);
      else gr->Fit(fun_Psigma_refit, "+R0", "goff", xmin_refit, xmax_refit);

      fitted_function_refit = gr->GetFunction(fun_Psigma_refit->GetName());
      for (int i=0; i<fitted_function_refit->GetNpar(); ++i) cout<< i << "\t " << fitted_function_refit->GetParName(i) << "\t " << fitted_function_refit->GetParameter(i) <<endl;

      y0cross = fitted_function_refit->GetParameter("x0");
      y05cross = 0;
      if (fitted_function_refit)
      {
         Int_t ichannel = channel - 1;

         // calculate intersection of tangent to the middle of the leading edge with axis x
         //Double_t fval = fitted_function_refit->Eval(pulse_xhalfmax);

         Double_t xt = fitted_function_refit->GetX(local_maximum_y/2., x[isigmin], local_maximum_x, eps);
         Double_t fval = fitted_function_refit->Eval(xt);
         Double_t dval = fitted_function_refit->Derivative(xt);
         Double_t dintersect = 0;
         //const Double_t eps = 1e-12;
         if (TMath::Abs(dval) > eps) dintersect = xt - fval/dval;
         oscFit->xt[ichannel] = xt;
         oscFit->yt[ichannel] = fval;
         oscFit->d[ichannel] = dintersect;
         oscFit->dydx[ichannel] = dval;
         cout<< "--> dintersect = " << dintersect << " xt = " << xt << " fval = " << fval << " dval = " << dval <<endl;

         if (gr) gr->SetTitle(Form("%s xt=%0.3f", gr->GetTitle(), xt));

         oscFit->t[ichannel] = fitted_function_refit->GetParameter("x0");
         oscFit->chi2[ichannel] = (fitted_function_refit->GetNDF() > 0)? fitted_function_refit->GetChisquare() / fitted_function_refit->GetNDF(): 0;
         oscFit->tau1[ichannel] = fitted_function_refit->GetParameter("tau1");
         oscFit->tau2[ichannel] = fitted_function_refit->GetParameter("tau2");
         oscFit->T[ichannel] = fitted_function_refit->GetParameter("T");
         oscFit->sigma[ichannel] = fitted_function_refit->GetParameter("sigma");

         // find time for charge 0 - 20
         Double_t time0 = x[isigmin];
         Double_t time = time0;
         Double_t dtime = 0.010;           // 10 ps step
         Double_t time_max = time0 + 10;   // 10 ns window
         Int_t iq[30] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};
         Int_t iq_curr = 0;
         while (time < time_max) {
            time += dtime;
            Double_t* param = 0;
            Double_t q = fitted_function_refit->Integral(time0, time, param, 1e-6);
            if (q > iq[iq_curr]) {
               oscFit->tq[ichannel][iq_curr] = time;
               if (debug) cout<< "oscFit->tq[" << ichannel << "][" << iq_curr << "] = " << oscFit->tq[ichannel][iq_curr] <<endl;
               ++iq_curr;
               if (iq_curr == 30) break;
            }
         }

         // rise time
         Double_t tf10 = fitted_function_refit->GetX(0.10*local_maximum_y, x[isigmin], local_maximum_x, eps);
         Double_t tf90 = fitted_function_refit->GetX(0.90*local_maximum_y, x[isigmin], local_maximum_x, eps);
         Double_t trisef = tf90 - tf10;
         if (debug) cout<< "tf10 = " << tf10 << " tf90 = " << tf90 << " ==> rise time (pulse function) = " << trisef <<endl;
         oscFit->trisef[channel-1] = trisef;
      }
   }

   if (debug) {
      TF1* fpl = 0;
      fpl = new TF1("fpl","pol1", x[x1_i-1]-eps,x[x2_i+1]+eps);
      // fpl->SetNpx(nsample);
      fpl->SetLineWidth(1);
      fpl->SetLineColor(2);
      fpl->SetParameter(0, p0);
      fpl->SetParameter(1, p1);
      Int_t ndf = 4 - 2;
      fpl->SetNDF(ndf);
      fpl->SetChisquare(chi2ndf*ndf);
      if (gr) gr->GetListOfFunctions()->Add(fpl);
      // fpl->Draw("same");
   }

   if (debug) {
      TMarker* marker_pt = new TMarker(ptx, pty, 5);
      marker_pt->SetMarkerColor(2);
      marker_pt->SetMarkerSize(2);
      marker_pt->Draw("same");
   }

   if (debug) {
      TMarker* marker_sigmin = new TMarker(x[isigmin], y[isigmin], 20);
      marker_sigmin->SetMarkerColor(1);
      marker_sigmin->Draw("same");
   }

   // assign output values
   oscFit->bkg[channel-1] = bkg_mean;
   oscFit->sbkg[channel-1] = bkg_sigma;
   // algorithms
   oscFit->yscal[channel-1] = (y[maximum_i-1] + y[maximum_i] + y[maximum_i+1]) / 3.;
   oscFit->sigmin[channel-1] = x[isigmin];
   oscFit->sigmax[channel-1] = x[isigmax];
   oscFit->lmaxx[channel-1] = local_maximum_x;
   oscFit->lmaxy[channel-1] = local_maximum_y;
   oscFit->pmax[channel-1] = pulse_maximum_y;
   oscFit->pti[channel-1] = integral;
   oscFit->sat[channel-1] = saturation;
   oscFit->pty[channel-1] = pty;
   oscFit->ptx[channel-1] = ptx;
   oscFit->plx[channel-1] = plx;
   oscFit->plx0[channel-1] = plx0;
   oscFit->plc[channel-1] = chi2ndf;
   oscFit->p0[channel-1] = p0;
   oscFit->p1[channel-1] = p1;
   // save fit in the tree
   oscFit->plxmin[channel-1] = plxmin;
   oscFit->plxmax[channel-1] = plxmax;
   oscFit->plp0[channel-1] = p0;
   oscFit->plp1[channel-1] = p1;

   oscFit->ptok[channel-1] = kTRUE
      && oscFit->pti[channel-1] > 0
      && oscFit->ptx[channel-1] > 0
      && oscFit->sat[channel-1] == kFALSE
   ;
   oscFit->itot[channel-1] = itot;
   oscFit->ibkg[channel-1] = integral_bkg;
   oscFit->itail[channel-1] = integral_tail;
}
/////////////////////////////////////////////////////////////////////////////////

bool debug = false;
bool gdebug = false;

// TTree* pulseosc(TTree *tree
//       , Int_t entry_first=0, Int_t entry_last=-1
//       , bool setdebug=false, bool setgdebug=false
//       , Double_t threshold = 0
//       , Int_t ignorech1=-1
//       , Int_t ignorech2=-1
//       , Int_t ignorech3=-1
//       , Int_t ignorech4=-1
//       , Int_t ignorech5=-1
//       , Int_t ignorech6=-1
//       , Int_t ignorech7=-1
//       , Int_t ignorech8=-1
//       )
// {
//    if (setdebug) debug = true;
//    if (setgdebug) gdebug = true;
// 
//    tree->SetMarkerStyle(7);
//    tree->SetMarkerColor(2);
// 
//    OscEvent* oscEvent = 0;
//    tree->SetBranchAddress("oscEvent",&oscEvent);
// 
//    cout<< "tree->GetEntries() = " << tree->GetEntries() <<endl;
// 
//    // output (fit results) tree
//    TTree* otree = new TTree("ft", "Fit result tree");
//    //-- otree->SetEstimate(1000000000L);
//    OscFit* oscFit = new OscFit;
//    otree->Branch("oscFit", "OscFit", &oscFit);
//    otree->SetFillStyle(3001);
//    otree->SetMarkerStyle(6);
//    otree->SetMarkerColor(2);
// 
//    // the number of channels in the data
//    tree->LoadTree(entry_first);
//    tree->GetEntry(entry_first);
//    Int_t nchannels = oscEvent->oscChannels->GetEntries();
// 
//    if (entry_last < 0) entry_last = tree->GetEntries() - 1;
// 
//    for (int jentry=entry_first; jentry<=entry_last; ++jentry)
//    {
//       tree->LoadTree(jentry);
//       tree->GetEntry(jentry);
//       cout<< "\n---------------------------------------------------> jentry = " << jentry << " evt = " << oscEvent->evt <<endl;
// 
//       oscFit->clear();
//       oscFit->evt = oscEvent->evt;
// 
//       cout<< ".. loop over osc channels. oscEvent->oscChannels = " << oscEvent->oscChannels->GetEntries() <<endl;
//       for (int ich=0; ich<nchannels; ++ich)
//       {
//       const OscChannel* oscChannel = (OscChannel*) oscEvent->oscChannels->At(ich);
//       cout<< "oscChannel->ch = " << oscChannel->ch <<endl;
// 
//       if (oscChannel->ch == ignorech1) continue;
//       if (oscChannel->ch == ignorech2) continue;
//       if (oscChannel->ch == ignorech3) continue;
//       if (oscChannel->ch == ignorech4) continue;
//       if (oscChannel->ch == ignorech5) continue;
//       if (oscChannel->ch == ignorech6) continue;
//       if (oscChannel->ch == ignorech7) continue;
//       if (oscChannel->ch == ignorech8) continue;
// 
//       //-- Double_t thres = 0.005;
//       // Double_t thres = 0.010;
//       //-- Double_t thres = 0.015;
//       Double_t thres = 0.020;
//       // Double_t thres = 10;
//       //-- Double_t thres = 5;
//       if (threshold > 0) thres = threshold;
//       
//       //
//       // pt algorithm
//       //
//       point(oscChannel->x, oscChannel->y
//             , oscFit
//             , thres, debug
//             , jentry, oscChannel->ch
//             );
//       oscFit->ptok[ich] = true
//          && oscFit->pti[ich] > 0
//          && oscFit->ptx[ich] > 0
//          && oscFit->sat[ich] == kFALSE
//       ;
//       }
//       otree->Fill();
//    }
// 
//    cout<< "\ngDirectory->GetName() = " << gDirectory->GetName() <<endl;
// 
//    // turn off debug
//    debug = false;
//    gdebug = false;
// 
//    tree->ResetBranchAddresses();
//    return otree;
// }
// 
// TTree* pulseosc(const char *ifname
//       , Int_t entry_first=0, Int_t entry_last=-1
//       , bool setdebug=false, bool setgdebug=false
//       , Double_t threshold = 0
//       , Int_t ignorech1=-1
//       , Int_t ignorech2=-1
//       , Int_t ignorech3=-1
//       , Int_t ignorech4=-1
//       , Int_t ignorech5=-1
//       , Int_t ignorech6=-1
//       , Int_t ignorech7=-1
//       , Int_t ignorech8=-1
//       )
// {
//    TFile* ifile = TFile::Open(ifname);
//    if (!ifile) {
//       cout<< "File not found: " << ifname <<endl;
//       return 0;
//    }
//    cout<< "Processing file " << ifname <<endl;
// 
//    TTree* tree = (TTree*) ifile->Get("t");
//    if (!tree) {
//       cout<< "tree \"t\" was not found in file " << ifname <<endl;
//       return 0;
//    }
// 
//    //-- Double_t thres = 0.005;
//    // Double_t thres = 0.010;
//    //-- Double_t thres = 0.015;
//    Double_t thres = 0.020;
//    // Double_t thres = 10;
//    //-- Double_t thres = 5;
//    if (threshold > 0) thres = threshold;
// 
//    // output file with tree "ft" ("fit tree")
//    // TFile* ofile = TFile::Open(Form("%s-ft.root",ifname),"recreate");
//    TFile* ofile = TFile::Open(Form("%s-ft_%0.3lf.root",ifname,threshold),"recreate");
// 
//    TTree* otree = pulseosc(tree, entry_first,entry_last, setdebug,setgdebug, threshold, ignorech1,ignorech2,ignorech3,ignorech4,ignorech5,ignorech6,ignorech7,ignorech8);
//    
//    cout<< "Writing " << otree->GetEntries() << " entries into file " << ofile->GetName() <<endl;
//    ofile->Write();
//    return otree;
// }

TTree* pulse3(TTree *tree
      , Int_t entry_first=0, Int_t entry_last=-1
      , bool setdebug=false, bool setgdebug=false
      , Double_t trig=0, Double_t gate=0
      , Double_t threshold1=0.020
      , Double_t threshold2=0
      , Double_t threshold3=0
      , Double_t threshold4=0
      , Double_t threshold5=0
      , Double_t threshold6=0
      , Double_t threshold7=0
      , Double_t threshold8=0
      )
{
   if (setdebug) debug = true;
   if (setgdebug) gdebug = true;

   Double_t threshold[8];
   threshold[0] = threshold1;
   threshold[1] = threshold2;
   threshold[2] = threshold3;
   threshold[3] = threshold4;
   threshold[4] = threshold5;
   threshold[5] = threshold6;
   threshold[6] = threshold7;
   threshold[7] = threshold8;
   Double_t threshold_current = -1.;
   for (int i=0; i<8; ++i) {
      if (threshold[i] > 0) {
         threshold_current = threshold[i];
         break;
      }
   }
   if (threshold_current <= 0) {
      cout<< "At least one threshold should be > 0" <<endl;
      return 0;
   }
   for (int i=0; i<8; ++i) {
      if (threshold[i] > 0) threshold_current = threshold[i];
      if (threshold[i] == 0) threshold[i] = threshold_current;
   }

   Float_t b1_t[1024], b1_c1[1024], b1_c2[1024], b1_c3[1024], b1_c4[1024];
   Float_t b2_t[1024], b2_c1[1024], b2_c2[1024], b2_c3[1024], b2_c4[1024];
   Int_t event, tc1, tc2;

   // clear buffers before the connection
   event = 0;
   tc1 = tc2 = 0;
   for (int i=0; i<1024; ++i) {
      b1_t[i] = 0;
      b1_c1[i] = 0;
      b1_c2[i] = 0;
      b1_c3[i] = 0;
      b1_c4[i] = 0;
      b2_t[i] = 0;
      b2_c1[i] = 0;
      b2_c2[i] = 0;
      b2_c3[i] = 0;
      b2_c4[i] = 0;
   }

   TBranch* t_event = tree->GetBranch("event");
   TBranch* t_tc1   = tree->GetBranch("tc1");
   TBranch* t_b1_t  = tree->GetBranch("b1_t")? tree->GetBranch("b1_t"): tree->GetBranch("t1");
   TBranch* t_b1_c1 = tree->GetBranch("b1_c1")? tree->GetBranch("b1_c1"): tree->GetBranch("c1");
   TBranch* t_b1_c2 = tree->GetBranch("b1_c2")? tree->GetBranch("b1_c2"): tree->GetBranch("c2");
   TBranch* t_b1_c3 = tree->GetBranch("b1_c3")? tree->GetBranch("b1_c3"): tree->GetBranch("c3");
   TBranch* t_b1_c4 = tree->GetBranch("b1_c4")? tree->GetBranch("b1_c4"): tree->GetBranch("c4");

   TBranch* t_tc2   = tree->GetBranch("tc2");
   TBranch* t_b2_t  = tree->GetBranch("b2_t")? tree->GetBranch("b2_t"): tree->GetBranch("t2");
   TBranch* t_b2_c1 = tree->GetBranch("b2_c1")? tree->GetBranch("b2_c1"): tree->GetBranch("c5");
   TBranch* t_b2_c2 = tree->GetBranch("b2_c2")? tree->GetBranch("b2_c2"): tree->GetBranch("c6");
   TBranch* t_b2_c3 = tree->GetBranch("b2_c3")? tree->GetBranch("b2_c3"): tree->GetBranch("c7");
   TBranch* t_b2_c4 = tree->GetBranch("b2_c4")? tree->GetBranch("b2_c4"): tree->GetBranch("c8");

   if (t_event) t_event->SetAddress(&event);
   if (t_tc1) t_tc1->SetAddress(&tc1);
   if (t_b1_t) t_b1_t->SetAddress(b1_t);
   if (t_b1_c1) t_b1_c1->SetAddress(b1_c1);
   if (t_b1_c2) t_b1_c2->SetAddress(b1_c2);
   if (t_b1_c3) t_b1_c3->SetAddress(b1_c3);
   if (t_b1_c4) t_b1_c4->SetAddress(b1_c4);

   if (t_tc2) t_tc2->SetAddress(&tc2);
   if (t_b2_t) t_b2_t->SetAddress(b2_t);
   if (t_b2_c1) t_b2_c1->SetAddress(b2_c1);
   if (t_b2_c2) t_b2_c2->SetAddress(b2_c2);
   if (t_b2_c3) t_b2_c3->SetAddress(b2_c3);
   if (t_b2_c4) t_b2_c4->SetAddress(b2_c4);

   if (t_b2_t) {
      cout<< "--> Found board #2" <<endl;
   }

   tree->SetMarkerStyle(7);
   tree->SetMarkerColor(2);

   cout<< "tree->GetEntries() = " << tree->GetEntries() <<endl;

   // output (fit results) tree
   TTree* otree = new TTree("ft", "Fit result tree");
   //-- otree->SetEstimate(1000000000L);
   OscFit* oscFit = new OscFit;
   otree->Branch("oscFit", "OscFit", &oscFit);
   otree->SetFillStyle(3001);
   otree->SetMarkerStyle(6);
   otree->SetMarkerColor(2);

   // the number of channels in the data
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

      for (int ich=0; ich<nchannels; ++ich)
      {
         //if (ich+1 == ignorech1) continue;
         //if (ich+1 == ignorech2) continue;
         //if (ich+1 == ignorech3) continue;
         //if (ich+1 == ignorech4) continue;
         //if (ich+1 == ignorech5) continue;
         //if (ich+1 == ignorech6) continue;
         //if (ich+1 == ignorech7) continue;
         //if (ich+1 == ignorech8) continue;
         if (threshold[ich] < 0) continue;

         // // Double_t thres = 0.005;
         // // Double_t thres = 0.010;
         // // Double_t thres = 0.015;
         // Double_t thres = 0.020;
         // // Double_t thres = 10.;                    // for tb2012_run_34.root.pulse.root
         // // Double_t thres = 5.;                     // for tb2012_run_34.root.pulse.root
         // if (threshold > 0) thres = threshold;
         Double_t thres = threshold[ich];
         
         //
         // pt algorithm
         //
         if (ich == 0 && t_b1_c1) point(b1_t, b1_c1, oscFit, thres, debug, jentry, ich+1, trig,gate);
         if (ich == 1 && t_b1_c2) point(b1_t, b1_c2, oscFit, thres, debug, jentry, ich+1, trig,gate);
         if (ich == 2 && t_b1_c3) point(b1_t, b1_c3, oscFit, thres, debug, jentry, ich+1, trig,gate);
         if (ich == 3 && t_b1_c4) point(b1_t, b1_c4, oscFit, thres, debug, jentry, ich+1, trig,gate);

         if (ich > 3 && t_b2_t == 0) continue;
         if (ich == 4 && t_b1_c1) point(b2_t, b2_c1, oscFit, thres, debug, jentry, ich+1, trig,gate);
         if (ich == 5 && t_b1_c2) point(b2_t, b2_c2, oscFit, thres, debug, jentry, ich+1, trig,gate);
         if (ich == 6 && t_b1_c3) point(b2_t, b2_c3, oscFit, thres, debug, jentry, ich+1, trig,gate);
         if (ich == 7 && t_b1_c4) point(b2_t, b2_c4, oscFit, thres, debug, jentry, ich+1, trig,gate);
      }
      otree->Fill();
   }

   cout<< "\ngDirectory->GetName() = " << gDirectory->GetName() <<endl;

   // turn off debug
   debug = false;
   gdebug = false;

   tree->ResetBranchAddresses();
   return otree;
}

TTree* pulse3(const char *ifname
      , Int_t entry_first=0, Int_t entry_last=-1
      , bool setdebug=false, bool setgdebug=false
      , Double_t trig=0, Double_t gate=0
      , Double_t threshold1=0.020
      , Double_t threshold2=0
      , Double_t threshold3=0
      , Double_t threshold4=0
      , Double_t threshold5=0
      , Double_t threshold6=0
      , Double_t threshold7=0
      , Double_t threshold8=0
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
      tree = (TTree*) ifile->Get("p");
      if (!tree) {
         cout<< "tree \"pulse\" was not found in file " << ifname <<endl;
         return 0;
      }
   }

   Double_t threshold[8];
   threshold[0] = threshold1;
   threshold[1] = threshold2;
   threshold[2] = threshold3;
   threshold[3] = threshold4;
   threshold[4] = threshold5;
   threshold[5] = threshold6;
   threshold[6] = threshold7;
   threshold[7] = threshold8;
   Double_t threshold_current = -1.;
   for (int i=0; i<8; ++i) {
      if (threshold[i] > 0) {
         threshold_current = threshold[i];
         break;
      }
   }
   if (threshold_current <= 0) {
      cout<< "At least one threshold should be > 0" <<endl;
      return 0;
   }
   for (int i=0; i<8; ++i) {
      if (threshold[i] > 0) threshold_current = threshold[i];
      if (threshold[i] == 0) threshold[i] = threshold_current;
   }

   Int_t maxbranch = 0;
   for (int i=0; i<8; ++i) {
      Int_t board = (i+1) / 4;
      if (tree->GetBranch(Form("c%d",i+1)) or tree->GetBranch(Form("c%d_b%d",i+1,board))) maxbranch = i;
   }

   std::string ofname = Form("%s-ft",ifname);
   for (int i=0; i<=maxbranch; ++i) {
      ofname = Form("%s_%0.3f",ofname.c_str(),threshold[i]);
   }
   ofname = Form("%s.root",ofname.c_str());

   TFile* ofile = TFile::Open(ofname.c_str(),"recreate");

   TTree* otree = pulse3(tree, entry_first,entry_last, setdebug,setgdebug, trig,gate, threshold1,threshold2,threshold3,threshold4,threshold5,threshold6,threshold7,threshold8);
   
   cout<< "Writing " << otree->GetEntries() << " entries into file " << ofile->GetName() <<endl;
   ofile->Write();
   return otree;
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using std::cout;     using std::endl;

class Ray: public TObject {
   Double_t radiant[3];
   Double_t dircos[3];
   Ray(): TObject() {
      for (int i=0; i<3; ++i) {
         radiant[i] = 0;
         dircos[i] = 0;
      }
   }
   Ray(Double_t x, Double_t y, Double_t z, Double_t cx, Double_t cy, Double_t cz): TObject() {
      radiant[0] = x;
      radiant[1] = y;
      radiant[2] = z;
      dircos[0] = cx;
      dircos[1] = cy;
      dircos[2] = cz;
      Double_t length = TMath::Sqrt(cx*cx + cy*cy + cz*cz);
      if (length > 0) for (int i=0; i<3; ++i) dircos[i] /= length;
   }

   ClassDef(Ray, 1);
};

class BulkyRangeDetector: public TObject
{
   static const Int_t ntiles = 64;
   Float_t d;                             // degrader thickness
   Float_t x1,y1, x2,y2, x3,y3, x4,y4;
   Float_t etot;                          // energy deposited in the degrader and bulky
   Float_t edeg;                          // energy deposited in the degrader
   Float_t ebulky;                        // energy deposited in the bulky
   Int_t last;
   Float_t e[ntiles];
   // calculated quantities
   Float_t ri[6], ro[6];         // components 0..2: radiant, 3..5: direction cosines
   Float_t rdegi[6], rdego[6];   // components 0..2: radiant, 3..5: direction cosines
   Float_t dr;                   // shift in r after the degrader
   Float_t dtheta;               // shift in angle after the degrader
public:
   void clear() {
      d = 0;
      x1 = y1 = x2 = y2 = x3 = y3 = x4 = y4 = 0;
      etot = 0;
      last = 0;
      for (int i=0; i<ntiles; ++i) e[i] = 0;
      for (int i=0; i<3; ++i) {
         ri[i] = ri[3+i] = 0;
         ro[i] = ro[3+i] = 0;
         rdegi[i] = rdegi[3+i] = 0;
         rdego[i] = rdego[3+i] = 0;
      }
      dr = 0;
      dtheta = 0;
   }
   BulkyRangeDetector(): TObject() {
      clear();
   }
   Float_t Sum(Int_t n2=0, Int_t n1=0) {
      if (n2 == 0) n2 = ntiles;
      if (n2 > ntiles) n2 = ntiles;
      Float_t sum = 0;
      for (int itile=n1; itile<n2; ++itile) sum += e[itile];
      return sum;
   }
   Int_t Last(Float_t thres=0) {
      Int_t last_tile = ntiles - 1;
      while (last_tile > 0 and e[last_tile] <= thres) --last_tile;
      return last_tile;
   }
   void SetDegrader(Float_t degrader) {d = degrader;}
   void Init(OscFit* oscFit, bool debug1=false)
   {
      // std::stringstream ss(line);
      // Int_t number;                 // number of the sensitive element in the geant4 txt data file: will be ignored
      // while (ss
      //       >> x1
      //       >> y1
      //       >> x2
      //       >> y2
      //       >> x3
      //       >> y3
      //       >> x4
      //       >> y4
      //       >> etot     // NB: degrader + bulky
      //       >> last
      //       >> number
      //       >> edeg
      //       >> number
      //       >> ebulky
      //       )
      // {
      //    // read the rest of the line
      //    for (int itile=0; itile<ntiles; ++itile) {
      //       ss >> number >> e[itile];
      //    }
      // }

      x1 = y1 = x2 = y2 = x3 = y3 = x4 = y4 = 0;
      etot = 0;
      last = 0;
      edeg = 0;
      ebulky = oscFit->pmax[4];
      for (int itile=0; itile<ntiles; ++itile) e[itile] = 0;
      e[0] = oscFit->pmax[0];
      e[1] = oscFit->pmax[1];
      e[2] = oscFit->pmax[2];

      // const Double_t z1 = -(15.225 + 10.3);  // cm
      // const Double_t z2 = -15.225;        // cm
      // const Double_t z3 = 15.225;         // cm
      // const Double_t z4 = 15.225 + 10.3;     // cm
      // Double_t eps = 1e-7;
      // Double_t t;

      // if (debug1) cout<< "x1 " << x1 << " y1 " << y1 << " z1 " << z1 << " x2 " << x2 << " y2 " << y2 << " z2 " << z2 <<endl;
      // if (debug1) cout<< "x3 " << x3 << " y3 " << y3 << " z3 " << z3 << " x4 " << x4 << " y4 " << y4 << " z4 " << z4 <<endl;
      // if (debug1) cout<< "etot " << etot << " last " << last << " ebulky " << ebulky <<endl;

      // // degrader input
      // ri[0] = x1;
      // ri[1] = y1;
      // ri[2] = z1;
      // // compute the direction cosines
      // ri[3+0] = x2 - x1;
      // ri[3+1] = y2 - y1;
      // ri[3+2] = z2 - z1;
      // Double_t ri_length = TMath::Sqrt(ri[3+0]*ri[3+0] + ri[3+1]*ri[3+1] + ri[3+2]*ri[3+2]);
      // for (int i=0; i<3; ++i) ri[3+i] /= ri_length;
      // if (debug1) cout<< "ri: " << ri[0] <<" "<< ri[1] <<" "<< ri[2] <<" cos: "<< ri[3] <<" "<< ri[4] <<" "<< ri[5] <<endl;

      // // coordinates of the degrader entrance
      // t = TMath::Abs(ri[3+2]) > eps? (-0.5*d - z1) / ri[3+2]: 0;   // parameter
      // for (int i=0; i<3; ++i) {
      //         rdegi[i] = ri[i] + ri[3+i]*t;
      //         rdegi[3+i] = ri[3+i];
      // }
      // if (debug1) cout<< "rdegi: t " << t <<" "<< rdegi[0] <<" "<< rdegi[1] <<" "<< rdegi[2] <<" cos: "<< rdegi[3] <<" "<< rdegi[4] <<" "<< rdegi[5] <<endl;

      // // degrader output
      // ro[0] = x3;
      // ro[1] = y3;
      // ro[2] = z3;
      // // compute the direction cosines
      // ro[3+0] = x4 - x3;
      // ro[3+1] = y4 - y3;
      // ro[3+2] = z4 - z3;
      // Double_t ro_length = TMath::Sqrt(ro[3+0]*ro[3+0] + ro[3+1]*ro[3+1] + ro[3+2]*ro[3+2]);
      // for (int i=0; i<3; ++i) ro[3+i] /= ro_length;
      // if (debug1) cout<< "ro: " << ro[0] <<" "<< ro[1] <<" "<< ro[2] <<" cos: "<< ro[3] <<" "<< ro[4] <<" "<< ro[5] <<endl;

      // // coordinates of the degrader exit
      // t = TMath::Abs(ro[3+2]) > eps? (0.5*d - z3) / ro[3+2]: 0;   // parameter. NB: move backward
      // for (int i=0; i<3; ++i) {
      //         rdego[i] = ro[i] + ro[3+i]*t;
      //         rdego[3+i] = ro[3+i];
      // }
      // if (debug1) cout<< "rdego: t " << t <<" "<< rdego[0] <<" "<< rdego[1] <<" "<< rdego[2] <<" cos: "<< rdego[3] <<" "<< rdego[4] <<" "<< rdego[5] <<endl;

      // // multiple scattering
      // dr = TMath::Sqrt((rdego[0]-rdegi[0])*(rdego[0]-rdegi[0]) + (rdego[1]-rdegi[1])*(rdego[1]-rdegi[1]));
      // dtheta = TMath::ACos(rdegi[3+0]*rdego[3+0] + rdegi[3+1]*rdego[3+1] + rdegi[3+2]*rdego[3+2]);
      // if (debug1) cout<< "dr = " << dr << " dtheta = " << dtheta <<endl;

      // // if (dr > 70.) {
      // //      cout<<endl;
      // //    cout<< "x1 " << x1 << " y1 " << y1 << " z1 " << z1 << " x2 " << x2 << " y2 " << y2 << " z2 " << z2 <<endl;
      // //    cout<< "x3 " << x3 << " y3 " << y3 << " z3 " << z3 << " x4 " << x4 << " y4 " << y4 << " z4 " << z4 <<endl;
      // //    cout<< "ri: " << ri[0] <<" "<< ri[1] <<" "<< ri[2] <<" cos: "<< ri[3] <<" "<< ri[4] <<" "<< ri[5] <<endl;
      // //    cout<< "rdegi: t " << t <<" "<< rdegi[0] <<" "<< rdegi[1] <<" "<< rdegi[2] <<" cos: "<< rdegi[3] <<" "<< rdegi[4] <<" "<< rdegi[5] <<endl;
      // //    cout<< "ro: " << ro[0] <<" "<< ro[1] <<" "<< ro[2] <<" cos: "<< ro[3] <<" "<< ro[4] <<" "<< ro[5] <<endl;
      // //    cout<< "rdego: t " << t <<" "<< rdego[0] <<" "<< rdego[1] <<" "<< rdego[2] <<" cos: "<< rdego[3] <<" "<< rdego[4] <<" "<< rdego[5] <<endl;
      // //    cout<< "dr = " << dr << " dtheta = " << dtheta <<endl;
      // // }
   }
   void Init(std::string& line, bool debug1=false)
   {
      std::stringstream ss(line);
      Int_t number;                    // number of the sensitive element in the geant4 txt data file: will be ignored
      while (ss
            >> x1
            >> y1
            >> x2
            >> y2
            >> x3
            >> y3
            >> x4
            >> y4
            >> etot        // NB: degrader + bulky
            >> last
            >> number
            >> edeg
            >> number
            >> ebulky
            )
      {
         // read the rest of the line
         for (int itile=0; itile<ntiles; ++itile) {
            ss >> number >> e[itile];
         }
      }

      const Double_t z1 = -(15.225 + 10.3);  // cm
      const Double_t z2 = -15.225;           // cm
      const Double_t z3 = 15.225;            // cm
      const Double_t z4 = 15.225 + 10.3;     // cm
      Double_t eps = 1e-7;
      Double_t t;

      if (debug1) cout<< "x1 " << x1 << " y1 " << y1 << " z1 " << z1 << " x2 " << x2 << " y2 " << y2 << " z2 " << z2 <<endl;
      if (debug1) cout<< "x3 " << x3 << " y3 " << y3 << " z3 " << z3 << " x4 " << x4 << " y4 " << y4 << " z4 " << z4 <<endl;
      if (debug1) cout<< "etot " << etot << " last " << last << " ebulky " << ebulky <<endl;

      // degrader input
      ri[0] = x1;
      ri[1] = y1;
      ri[2] = z1;
      // compute the direction cosines
      ri[3+0] = x2 - x1;
      ri[3+1] = y2 - y1;
      ri[3+2] = z2 - z1;
      Double_t ri_length = TMath::Sqrt(ri[3+0]*ri[3+0] + ri[3+1]*ri[3+1] + ri[3+2]*ri[3+2]);
      for (int i=0; i<3; ++i) ri[3+i] /= ri_length;
      if (debug1) cout<< "ri: " << ri[0] <<" "<< ri[1] <<" "<< ri[2] <<" cos: "<< ri[3] <<" "<< ri[4] <<" "<< ri[5] <<endl;

      // coordinates of the degrader entrance
      t = TMath::Abs(ri[3+2]) > eps? (-0.5*d - z1) / ri[3+2]: 0;   // parameter
      for (int i=0; i<3; ++i) {
         rdegi[i] = ri[i] + ri[3+i]*t;
         rdegi[3+i] = ri[3+i];
      }
      if (debug1) cout<< "rdegi: t " << t <<" "<< rdegi[0] <<" "<< rdegi[1] <<" "<< rdegi[2] <<" cos: "<< rdegi[3] <<" "<< rdegi[4] <<" "<< rdegi[5] <<endl;

      // degrader output
      ro[0] = x3;
      ro[1] = y3;
      ro[2] = z3;
      // compute the direction cosines
      ro[3+0] = x4 - x3;
      ro[3+1] = y4 - y3;
      ro[3+2] = z4 - z3;
      Double_t ro_length = TMath::Sqrt(ro[3+0]*ro[3+0] + ro[3+1]*ro[3+1] + ro[3+2]*ro[3+2]);
      for (int i=0; i<3; ++i) ro[3+i] /= ro_length;
      if (debug1) cout<< "ro: " << ro[0] <<" "<< ro[1] <<" "<< ro[2] <<" cos: "<< ro[3] <<" "<< ro[4] <<" "<< ro[5] <<endl;

      // coordinates of the degrader exit
      t = TMath::Abs(ro[3+2]) > eps? (0.5*d - z3) / ro[3+2]: 0;   // parameter. NB: move backward
      for (int i=0; i<3; ++i) {
         rdego[i] = ro[i] + ro[3+i]*t;
         rdego[3+i] = ro[3+i];
      }
      if (debug1) cout<< "rdego: t " << t <<" "<< rdego[0] <<" "<< rdego[1] <<" "<< rdego[2] <<" cos: "<< rdego[3] <<" "<< rdego[4] <<" "<< rdego[5] <<endl;

      // multiple scattering
      dr = TMath::Sqrt((rdego[0]-rdegi[0])*(rdego[0]-rdegi[0]) + (rdego[1]-rdegi[1])*(rdego[1]-rdegi[1]));
      dtheta = TMath::ACos(rdegi[3+0]*rdego[3+0] + rdegi[3+1]*rdego[3+1] + rdegi[3+2]*rdego[3+2]);
      if (debug1) cout<< "dr = " << dr << " dtheta = " << dtheta <<endl;

      // if (dr > 70.) {
      //         cout<<endl;
      //    cout<< "x1 " << x1 << " y1 " << y1 << " z1 " << z1 << " x2 " << x2 << " y2 " << y2 << " z2 " << z2 <<endl;
      //    cout<< "x3 " << x3 << " y3 " << y3 << " z3 " << z3 << " x4 " << x4 << " y4 " << y4 << " z4 " << z4 <<endl;
      //    cout<< "ri: " << ri[0] <<" "<< ri[1] <<" "<< ri[2] <<" cos: "<< ri[3] <<" "<< ri[4] <<" "<< ri[5] <<endl;
      //    cout<< "rdegi: t " << t <<" "<< rdegi[0] <<" "<< rdegi[1] <<" "<< rdegi[2] <<" cos: "<< rdegi[3] <<" "<< rdegi[4] <<" "<< rdegi[5] <<endl;
      //    cout<< "ro: " << ro[0] <<" "<< ro[1] <<" "<< ro[2] <<" cos: "<< ro[3] <<" "<< ro[4] <<" "<< ro[5] <<endl;
      //    cout<< "rdego: t " << t <<" "<< rdego[0] <<" "<< rdego[1] <<" "<< rdego[2] <<" cos: "<< rdego[3] <<" "<< rdego[4] <<" "<< rdego[5] <<endl;
      //    cout<< "dr = " << dr << " dtheta = " << dtheta <<endl;
      // }
   }

   friend bool process_file(Int_t,Int_t);

   ClassDef(BulkyRangeDetector, 5);
};

ClassImp(Ray);
ClassImp(BulkyRangeDetector);

bool process_file(Int_t run, Int_t idegrader, bool debug1=false)
{
   TTree* tree = new TTree(Form("t%d",idegrader), Form("degrader thickness %d mm WET",idegrader));
   BulkyRangeDetector bulkyRangeDetector;
   tree->Branch("BulkyRangeDetector", "BulkyRangeDetector", &bulkyRangeDetector);
   tree->SetFillStyle(3001);
   tree->SetMarkerStyle(7);
   tree->SetMarkerColor(2);

   TDirectory* current_dir = gDirectory;
   TFile* ifile = TFile::Open(Form("run_%d.dat.root-ft_0.004_0.005_0.004_0.020.root",run));
   TTree* ft = (TTree*) ifile->Get("ft");

   OscFit* oscFit = 0;
   if (ft) ft->SetBranchAddress("oscFit", &oscFit);

   bulkyRangeDetector.SetDegrader(0.1*idegrader);     // convert from mm to cm
   for (int jentry=0; jentry<ft->GetEntries(); ++jentry)
   {
      ft->LoadTree(jentry);
      ft->GetEntry(jentry);
      bulkyRangeDetector.clear();
      bulkyRangeDetector.Init(oscFit, debug1);
      tree->Fill();
   }

   ifile->Close();
   current_dir->cd();

   return true;
}

void geant4bulkyrange_data()
{
   bool debug1 = false;
   // bool debug1 = true;

   const char* ofname = "deg_bulkyrange_data.root";
   TFile* ofile = TFile::Open(ofname, "recreate");

   // // for (int idegrader=100; idegrader<=100; ++idegrader)
   // // for (int idegrader=100; idegrader<=101; ++idegrader)
   // for (int idegrader=0; idegrader<=245; ++idegrader)
   // {
   //    if (process_file(idegrader, debug1)) cout<< "done idegrader = " << idegrader <<endl;
   //    else return;
   // }

   // if (process_file(40, 0, debug1)) cout<< "done idegrader = " << 0 <<endl;
   // ofile->cd();
   // if (process_file(42, 137, debug1)) cout<< "done idegrader = " << 137 <<endl;
   // ofile->cd();
   // if (process_file(43, 142, debug1)) cout<< "done idegrader = " << 142 <<endl;
   // ofile->cd();
   // if (process_file(44, 145, debug1)) cout<< "done idegrader = " << 145 <<endl;
   // ofile->cd();
   // if (process_file(45, 149, debug1)) cout<< "done idegrader = " << 149 <<endl;
   // ofile->cd();
   // if (process_file(46, 153, debug1)) cout<< "done idegrader = " << 153 <<endl;
   // ofile->cd();
   // if (process_file(47, 154, debug1)) cout<< "done idegrader = " << 154 <<endl;
   // ofile->cd();
   // if (process_file(48, 155, debug1)) cout<< "done idegrader = " << 155 <<endl;
   // ofile->cd();
   // if (process_file(49, 158, debug1)) cout<< "done idegrader = " << 158 <<endl;
   // ofile->cd();
   // if (process_file(50, 162, debug1)) cout<< "done idegrader = " << 162 <<endl;
   // ofile->cd();

   if (process_file(40, 0, debug1)) cout<< "done idegrader = " << 0 <<endl;
   if (process_file(42, 137, debug1)) cout<< "done idegrader = " << 137 <<endl;
   if (process_file(43, 142, debug1)) cout<< "done idegrader = " << 142 <<endl;
   if (process_file(44, 145, debug1)) cout<< "done idegrader = " << 145 <<endl;
   if (process_file(45, 149, debug1)) cout<< "done idegrader = " << 149 <<endl;
   if (process_file(46, 153, debug1)) cout<< "done idegrader = " << 153 <<endl;
   if (process_file(47, 154, debug1)) cout<< "done idegrader = " << 154 <<endl;
   if (process_file(48, 155, debug1)) cout<< "done idegrader = " << 155 <<endl;
   if (process_file(49, 158, debug1)) cout<< "done idegrader = " << 158 <<endl;
   if (process_file(50, 162, debug1)) cout<< "done idegrader = " << 162 <<endl;

   ofile->Write();
}

TGraphAsymmErrors* counting_rate(Int_t itile, Float_t thres, Int_t ideg1, Int_t ideg2)
{
   //TH1F* hcount = new TH1F(Form("hcount_tile_%d_thres_%0.1f",itile,thres), Form("hcount_tile_%d_thres_%0.1f",itile,thres), (ideg2+1)-ideg1+1,ideg1,(ideg2+1));

   Double_t d[1000];
   Double_t ratio[1000];
   Double_t eLratio[1000];
   Double_t eHratio[1000];
   Int_t np = 0;

   TEfficiency tEfficiency;
   for (int ideg=ideg1; ideg<=ideg2; ++ideg)
   {
      TTree* tree = (TTree*) gDirectory->Get(Form("t%d",ideg));
      if (!tree) {
         cout<< "Tree not found: t" << ideg <<endl;
         return 0;
      }
      Int_t nbulky = tree->Draw("ebulky","ebulky>0","goff");
      Int_t ntile = tree->Draw(Form("e[%d]",itile),Form("e[%d]>%f",itile,thres),"goff");
      nbulky > 0? ratio[np] = Double_t(ntile)/Double_t(nbulky): 0;
      eLratio[np] = ratio[np] - tEfficiency.ClopperPearson(nbulky, ntile, 0.683, kFALSE);
      eHratio[np] = tEfficiency.ClopperPearson(nbulky, ntile, 0.683, kTRUE) - ratio[np];
      d[np] = ideg;
      //Int_t ibin = ideg - ideg1 + 1;
      //hcount->SetBinContent(ibin, ratio);
      ++np;
   }

   TGraphAsymmErrors* gr_tile = new TGraphAsymmErrors(np, d, ratio, 0,0, eLratio,eHratio);
   gr_tile->SetNameTitle("gr_tile", Form("Efficiency of tile No. %d vs degrader, threshold %0.1f MeV",itile,thres));
   gr_tile->SetMarkerStyle(20);
   gr_tile->SetMarkerColor(1);
   gr_tile->SetLineColor(1);

   new TCanvas;
   //hcount->Draw();
   gr_tile->Draw("ap");

   // find position of the halfmax
   //Double_t halfmax = hcount->GetBinContent(1) / 2.;
   //Double_t halfmax_x = hcount->GetBinCenter(1);
   //for (int i=1; i<=hcount->GetNbinsX(); ++i) {
   //   halfmax_x = hcount->GetBinCenter(i);
   //   if (hcount->GetBinContent(i) < halfmax) break;
   //}
   Double_t halfmax = gr_tile->GetY()[0] / 2;
   Double_t halfmax_x = gr_tile->GetX()[0];
   for (int i=1; i<=gr_tile->GetN(); ++i) {
      halfmax_x = gr_tile->GetX()[i];
      if (gr_tile->GetY()[i] < halfmax) break;
   }

   Double_t eps = 1e-7;
   TF1* feff = new TF1("feff","0.5*[0]*(1 - TMath::Erf((x-[1])/(TMath::Sqrt(2)*[2])))", ideg1-eps,ideg2+eps);
   feff->SetParName(0, "A");
   feff->SetParName(1, "mean");
   feff->SetParName(2, "#sigma");
   feff->SetParameter(0, 1.);
   feff->SetParameter(1, halfmax_x);
   feff->SetParameter(2, 4.);
   
   //hcount->Fit(feff);
   gr_tile->Fit(feff);

   return gr_tile;
}
