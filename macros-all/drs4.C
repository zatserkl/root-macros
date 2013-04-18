#include <TROOT.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <TMath.h>

#include <iostream>
#include <string>
#include <cstring>
#include <ctime>
#include <cstdio>

using std::cout;     using std::endl;

// ~/macros/utils.C
TGraph* gtemp(TCanvas* can=0);
TH1* htemp(TCanvas* can=0);
const char* nextname(const char* base="h");

Float_t ipulse(Float_t xmin, Float_t xmax, Float_t bkg=0, Float_t sign=1, TGraph* gr=0)
{
   if (!gr) gr = gtemp();
   if (!gr) {
      cout<< "Could not found TGraph in the current canvas " << gPad->GetName() <<endl;
      return 0;
   }

   Float_t sum=0;
   for (int i=0; i<gr->GetN()-1; ++i) {
      if (gr->GetX()[i] >= xmin && gr->GetX()[i] <= xmax) {
         Float_t val = sign*(gr->GetY()[i] - bkg);
         Float_t dx = gr->GetX()[i+1] - gr->GetX()[i];
         sum += val*dx;
      }
   }

   return sum;
}

TGraph* gdifferentiate(Double_t xmin=0, Double_t xmax=0, Double_t sign=1, TGraph* gr=0)
{
   const char* title = 0;

   if (gr) title = gr->GetTitle();
   else {
      gr = gtemp();
      title = gr->GetTitle();
      if (std::strcmp(gr->GetTitle(), "Graph") == 0) {
         // try to use title of the histogram embedded into the current canvas
         TH1* h = htemp();
         if (h) title = h->GetTitle();
      }
   }

   if (!gr) {
      cout<< "Could not found TGraph in the current canvas " << gPad->GetName() <<endl;
      return 0;
   }

   if (xmin == 0 and xmax == 0) {
      xmin = gr->GetX()[0];
      xmax = gr->GetX()[gr->GetN()-1];
   }
   if (xmax < xmin) {
      xmax = gr->GetX()[gr->GetN()-1];
   }
   if (xmin >= xmax) {
      cout<< "Limit(s) are out of range (" << gr->GetX()[0] << ", " << gr->GetX()[gr->GetN()-1] << ")" <<endl;
      return 0;
   }
   // cout<< "check imin imax" <<endl;

   Double_t* x = new Double_t[gr->GetN()];
   Double_t* y = new Double_t[gr->GetN()];
   Int_t np = 0;

   for (int i=0; i<gr->GetN(); ++i) {
      if (gr->GetX()[i] < xmin) continue;
      if (gr->GetX()[i] > xmax) break;
      x[np] = gr->GetX()[i];
      Double_t dx = i > 0? gr->GetX()[i] - gr->GetX()[i-1]: gr->GetX()[i+1] - gr->GetX()[i];
      Double_t dy = i > 0? sign*(gr->GetY()[i] - gr->GetY()[i-1]): 0;
      y[np] = dy / dx;
      //cout<< i << "\t" << " x[" << np << "] " << x[np] << " y " << y[np] <<endl;
      ++np;
   }

   TGraph* grdiff = new TGraph(np, x, y);
   grdiff->SetNameTitle(Form("i_%s", gr->GetName()), Form("diff of %s", title));
   grdiff->SetMarkerStyle(gr->GetMarkerStyle());
   grdiff->SetMarkerColor(gr->GetMarkerColor());
   grdiff->SetLineColor(gr->GetLineColor());

   new TCanvas;
   grdiff->Draw("apl");

   delete[] x;
   delete[] y;
   return grdiff;
}

void gintegrate(Double_t xmin=0, Double_t xmax=0, Double_t bkg=0, Double_t sign=1, TGraph* gr=0)
{
   const char* title = 0;

   if (gr) title = gr->GetTitle();
   else {
      gr = gtemp();
      title = gr->GetTitle();
      if (std::strcmp(gr->GetTitle(), "Graph") == 0) {
         // try to use title of the histogram embedded into the current canvas
         TH1* h = htemp();
         if (h) title = h->GetTitle();
      }
   }

   if (!gr) {
      cout<< "Could not found TGraph in the current canvas " << gPad->GetName() <<endl;
      return;
   }

   if (xmin >= xmax and xmax == 0) {
      xmin = gr->GetX()[0];
      xmax = gr->GetX()[gr->GetN()-1];
   }
   //cout<< "xmin = " << xmin << " xmax = " << xmax <<endl;

   Double_t x[1024], y[1024];
   Int_t np = 0;

   Double_t integral = 0;
   for (int i=0; i<gr->GetN()-1; ++i) {
      if (gr->GetX()[i] < xmin) continue;
      if (gr->GetX()[i] > xmax) break;
      x[np] = gr->GetX()[i];
      Double_t dx = i < gr->GetN()-1? gr->GetX()[i+1] - gr->GetX()[i]: gr->GetX()[i] - gr->GetX()[i-1];
      integral += (sign*gr->GetY()[i] - bkg) * dx;
      y[np] = integral;
      //cout<< i << "\t" << " x[" << np << "] " << x[np] << " y " << y[np] <<endl;
      ++np;
   }

   TGraph* gri = new TGraph(np, x, y);
   gri->SetNameTitle(Form("i_%s", gr->GetName()), Form("integral of %s", title));
   gri->SetMarkerStyle(gr->GetMarkerStyle());
   gri->SetMarkerColor(gr->GetMarkerColor());
   gri->SetLineColor(gr->GetLineColor());

   new TCanvas;
   gri->Draw("apl");
}

/*
root -l mppc2983-71.8V_trig-71.2V_dark_1.dat.root
p->Draw("-c1-4.02190e-03:t1","Entry$==177","lp")
accumulate()
*/
void accumulate(Float_t bkg=0, Float_t sign=1, TGraph* gr=0)
{
   if (!gr) gr = gtemp();
   if (!gr) {
      cout<< "Could not found TGraph in the current canvas " << gPad->GetName() <<endl;
      return;
   }

   Float_t x[1024], y[1024];
   Int_t np = 0;

   Float_t integral=0;
   for (int i=0; i<gr->GetN(); ++i) {
      Float_t val = sign*(gr->GetY()[i] - bkg);
      Float_t dx = (i < gr->GetN()-1)? gr->GetX()[i+1] - gr->GetX()[i]: gr->GetX()[i] - gr->GetX()[i-1];
      integral += val*dx;

      x[np] = gr->GetX()[i];
      y[np] = integral;
      ++np;
   }

   // double differential to kill the background y = kx + b
   //
   // Float_t xd2[1024], yd2[1024];
   // Int_t npd2 = 0;

   // for (int i=1; i<gr->GetN()-1; ++i) {
   //    // NB: assumes the same step in x
   //    Float_t Nim1 = gr->GetY()[i-1];
   //    Float_t Ni = gr->GetY()[i];
   //    Float_t Nip1 = gr->GetY()[i+1];
   //    xd2[npd2] = gr->GetX()[i];
   //    yd2[npd2] = Nip1 - 2*Ni + Nim1;
   //    ++npd2;
   // }

   // TGraph* gr_d2 = new TGraph(npd2, xd2, yd2);
   // gr_d2->SetNameTitle("gr_d2", "gr_d2");
   // gr_d2->SetMarkerStyle(7);
   // gr_d2->SetMarkerColor(2);

   // new TCanvas;
   // gr_d2->Draw("ap");

   Float_t x32[32], y32[32];
   Float_t x16[64], y16[64];
   Float_t x8[128], y8[128];
   for (int i=0; i<32; ++i) {
      x32[i] = 0;
      y32[i] = 0;
   }
   for (int i=0; i<64; ++i) {
      x16[i] = 0;
      y16[i] = 0;
   }
   for (int i=0; i<128; ++i) {
      x8[i] = 0;
      y8[i] = 0;
   }

   for (int i=0; i<gr->GetN(); ++i) {     // worked for 1024
      x32[i/32] += gr->GetX()[i];
      y32[i/32] += gr->GetY()[i];
      x16[i/16] += gr->GetX()[i];
      y16[i/16] += gr->GetY()[i];
      x8[i/8] += gr->GetX()[i];
      y8[i/8] += gr->GetY()[i];
   }
   Int_t np32 = gr->GetN() / 32;        // the number of cells filled with 32 values summed
   Int_t np32last = gr->GetN() % 32;    // the number of values in the last cell
   for (int i=0; i<np32; ++i) {
      x32[i] /= 32;
      y32[i] /= 32;
   }
   if (np32last > 0) {                  // process the last cell with sum of less than 32 values
      x32[np32] /= np32last;
      y32[np32] /= np32last;
      ++np32;                           // add this last cell
   }
   Int_t np16 = gr->GetN() / 16;        // the number of cells filled with 16 values summed
   Int_t np16last = gr->GetN() % 16;    // the number of values in the last cell
   for (int i=0; i<np16; ++i) {
      x16[i] /= 16;
      y16[i] /= 16;
   }
   if (np16last > 0) {                  // process the last cell with sum of less than 32 values
      x16[np16] /= np16last;
      y16[np16] /= np16last;
      ++np16;                           // add this last cell
   }
   Int_t np8 = gr->GetN() / 8;        // the number of cells filled with 8 values summed
   Int_t np8last = gr->GetN() % 8;    // the number of values in the last cell
   for (int i=0; i<np8; ++i) {
      x8[i] /= 8;
      y8[i] /= 8;
   }
   if (np8last > 0) {                  // process the last cell with sum of less than 8 values
      x8[np8] /= np8last;
      y8[np8] /= np8last;
      ++np8;                           // add this last cell
   }

   cout<< "np32 = " << np32 << " np32last = " << np32last << " np16 = " << np16 << " np16last = " << np16last << " np8 = " << np8 << " np8last = " << np8last <<endl;

   TGraph* gr_int = new TGraph(np, x, y);
   gr_int->SetNameTitle("gr_int", "accumulated  charge (integral)");
   gr_int->SetMarkerStyle(7);

   new TCanvas;
   gr_int->Draw("ap");

   TGraph* gr32 = new TGraph(np32, x32, y32);
   gr32->SetNameTitle("gr32", "gr32");
   gr32->SetMarkerStyle(7);

   new TCanvas;
   gr32->Draw("ap");

   TGraph* gr16 = new TGraph(np16, x16, y16);
   gr16->SetNameTitle("gr16", "gr16");
   gr16->SetMarkerStyle(7);

   new TCanvas;
   gr16->Draw("ap");

   TGraph* gr8 = new TGraph(np8, x8, y8);
   gr8->SetNameTitle("gr8", "gr8");
   gr8->SetMarkerStyle(7);

   new TCanvas;
   gr8->Draw("ap");

   Float_t xs[1024], ys[1024];
   Int_t nps = 0;
   for (int i=0+2; i<gr->GetN()-2; ++i) {
      Float_t ysmooth = (
            1.*gr->GetY()[i-2] + 
            2.*gr->GetY()[i-1] +
            3.*gr->GetY()[i] +
            2.*gr->GetY()[i+1] +
            1.*gr->GetY()[i+2]
      ) / 9;
      ys[nps] = ysmooth;
      xs[nps] = gr->GetX()[i];
      ++nps;
   }

   TGraph* gr_s = new TGraph(nps, xs, ys);
   gr_s->SetNameTitle("gr_s", "smoothed on 5 points");
   gr_s->SetMarkerStyle(7);
   gr_s->SetMarkerColor(2);

   new TCanvas;
   gr_s->Draw("apl");

   Float_t xsg[1024], ysg[1024];
   Int_t npsg = 0;
   for (int i=0+2; i<gr->GetN()-2; ++i) {
      Float_t ysmooth = (
            0.043*gr->GetY()[i-2] + 
            0.241*gr->GetY()[i-1] +
            0.428*gr->GetY()[i] +
            0.241*gr->GetY()[i+1] +
            0.043*gr->GetY()[i+2]
      ) / 0.936;
      ysg[npsg] = ysmooth;
      xsg[npsg] = gr->GetX()[i];
      ++npsg;
   }

   TGraph* gr_sg = new TGraph(npsg, xsg, ysg);
   gr_sg->SetNameTitle("gr_sg", "gaussian smoothed on 5 points");
   gr_sg->SetMarkerStyle(7);
   gr_sg->SetMarkerColor(8);

   new TCanvas;
   gr_sg->Draw("apl");

   Float_t xs2[1024], ys2[1024];
   Int_t nps2 = 0;
   for (int i=0+2; i<gr_s->GetN()-2; ++i) {
      Float_t ysmooth = (
            1.*gr_s->GetY()[i-2] + 
            2.*gr_s->GetY()[i-1] +
            3.*gr_s->GetY()[i] +
            2.*gr_s->GetY()[i+1] +
            1.*gr_s->GetY()[i+2]
      ) / 9;
      ys2[nps2] = ysmooth;
      xs2[nps2] = gr_s->GetX()[i];
      ++nps2;
   }

   TGraph* gr_s2 = new TGraph(nps2, xs2, ys2);
   gr_s2->SetNameTitle("gr_s2", "smoothed on 5 points");
   gr_s2->SetMarkerStyle(7);
   gr_s2->SetMarkerColor(4);

   new TCanvas;
   gr_s2->Draw("apl");
}

void smooth(TGraph* gr=0)
{
   if (!gr) gr = gtemp();
   if (!gr) {
      cout<< "Could not found TGraph in the current canvas " << gPad->GetName() <<endl;
      return;
   }

   Float_t xs[1024], ys[1024];
   Int_t nps = 0;
   ys[nps] = (3.*gr->GetY()[nps] + 2.*gr->GetY()[nps+1] + 1.*gr->GetY()[nps+2]) / 6.;
   xs[nps] = gr->GetX()[nps];
   nps++;
   ys[nps] = (2.*gr->GetY()[nps-1] + 3.*gr->GetY()[nps] + 2.*gr->GetY()[nps+1] + 1.*gr->GetY()[nps+2]) / 8.;
   xs[nps] = gr->GetX()[nps];
   nps++;
   for (int i=0+2; i<gr->GetN()-2; ++i) {
      Float_t ysmooth = (
            1*gr->GetY()[i-2] + 
            2*gr->GetY()[i-1] +
            3*gr->GetY()[i] +
            2*gr->GetY()[i+1] +
            1*gr->GetY()[i+2]
      ) / 9;
      ys[nps] = ysmooth;
      xs[nps] = gr->GetX()[i];
      ++nps;
   }
   ys[nps] = (2.*gr->GetY()[nps+1] + 3.*gr->GetY()[nps] + 2.*gr->GetY()[nps-1] + 1.*gr->GetY()[nps-2]) / 8.;
   xs[nps] = gr->GetX()[nps];
   nps++;
   ys[nps] = (3.*gr->GetY()[nps] + 2.*gr->GetY()[nps-1] + 1.*gr->GetY()[nps-2]) / 6.;
   xs[nps] = gr->GetX()[nps];
   nps++;

   TGraph* gr_smooth = new TGraph(nps, xs, ys);
   gr_smooth->SetNameTitle(Form("%s_smooth",gr->GetName()), Form("%s smoothed on 5 points",gr->GetTitle()));
   gr_smooth->SetMarkerStyle(gr->GetMarkerStyle());
   gr_smooth->SetMarkerColor(8);
   gr_smooth->SetLineColor(8);

   new TCanvas;
   gr_smooth->Draw("apl");
}

void hgr(Float_t xmin=0, Float_t xmax=0, TGraph* gr=0)
{
   if (!gr) gr = gtemp();
   if (!gr) {
      cout<< "Could not found TGraph in the current canvas " << gPad->GetName() <<endl;
      return;
   }

   if (xmin == 0 and xmax == 0) {
      xmin = gr->GetX()[0];
      xmax = gr->GetX()[gr->GetN()-1];
   }
   if (xmax < xmin) {
      xmax = gr->GetX()[gr->GetN()-1];
   }
   if (xmin >= xmax) {
      cout<< "Limit(s) are out of range (" << gr->GetX()[0] << ", " << gr->GetX()[gr->GetN()-1] << ")" <<endl;
      return;
   }

   Int_t imin = -1;
   Int_t imax = -1;
   for (int i=0; i<gr->GetN(); ++i) {
      if (imin < 0 and gr->GetX()[i] >= xmin) imin = i;
      if (gr->GetX()[i] > xmax) break;
      imax = i;
   }
   //cout<< "imin = " << imin << " imax = " << imax <<endl;

   Double_t ymin = gr->GetY()[imin];
   Double_t ymax = gr->GetY()[imin];
   for (int i=imin+1; i<=imax; ++i) {
      if (gr->GetY()[i] < ymin) ymin = gr->GetY()[i];
      if (gr->GetY()[i] > ymax) ymax = gr->GetY()[i];
   }

   std::string hname = nextname("hgr");
   TH1F* hgr = new TH1F(hname.c_str(), Form("h for %s", gr->GetTitle()), 400, ymin, ymax);

   for (int i=imin; i<=imax; ++i) hgr->Fill(gr->GetY()[i]);

   new TCanvas;
   hgr->Draw();
   hgr->Fit("gaus","L");
}

void ptime(Int_t entry1=0, Int_t entry2=-1, TTree* p=0)
{
   if (!p) p = (TTree*) gDirectory->Get("p");
   if (!p) {
      cout<< "Could not find tree p" <<endl;
      return;
   }

   if (entry2 < entry1) entry2 = p->GetEntries()-1;

   /*
      tm_sec	seconds after the minute	0-61*
      tm_min	minutes after the hour	0-59
      tm_hour	hours since midnight	0-23
      tm_mday	day of the month	1-31
      tm_mon	months since January	0-11
      tm_year	years since 1900	
      tm_wday	days since Sunday	0-6
      tm_yday	days since January 1	0-365
      tm_isdst	Daylight Saving Time flag	
      The Daylight Saving Time flag (tm_isdst) is greater than zero if Daylight Saving Time is in effect, zero if Daylight Saving Time is not in effect, and less than zero if the information is not available.
    * tm_sec is generally 0-59. Extra range to accommodate for leap seconds in certain systems.
    */

   struct std::tm time1;
   time1.tm_isdst = -1;       // set to "not available" otherwise mktime may corrupt the struct

   Int_t nbytes = 0;
   nbytes = p->GetEntry(entry1);
   if (nbytes <= 0) {
      cout<< "Could not load entry " << entry1 <<endl;
      return;
   }

   time1.tm_year = (Int_t) p->GetLeaf("year")->GetValue() - 1900;
   time1.tm_mon = (Int_t) p->GetLeaf("month")->GetValue() - 1;
   time1.tm_mday = (Int_t) p->GetLeaf("day")->GetValue() - 1;
   time1.tm_hour = (Int_t) p->GetLeaf("hour")->GetValue();
   time1.tm_min = (Int_t) p->GetLeaf("minute")->GetValue();
   time1.tm_sec = (Int_t) p->GetLeaf("second")->GetValue();
   Int_t millisecond1 = (Int_t) p->GetLeaf("millisecond")->GetValue();

   printf("Entry %8d date: %02d/%02d/%d %02d:%02d:%02d.%03d\n", entry1, time1.tm_mon,time1.tm_mday,time1.tm_year+1900, time1.tm_hour,time1.tm_min,time1.tm_sec,millisecond1);

   struct std::tm time2;
   time2.tm_isdst = -1;       // set to "not available" otherwise mktime may corrupt the struct

   nbytes = 0;
   nbytes = p->GetEntry(entry2);
   if (nbytes <= 0) {
      cout<< "Could not load entry " << entry2 <<endl;
      return;
   }

   time2.tm_year = (Int_t) p->GetLeaf("year")->GetValue() - 1900;
   time2.tm_mon = (Int_t) p->GetLeaf("month")->GetValue() - 1;
   time2.tm_mday = (Int_t) p->GetLeaf("day")->GetValue() - 1;
   time2.tm_hour = (Int_t) p->GetLeaf("hour")->GetValue();
   time2.tm_min = (Int_t) p->GetLeaf("minute")->GetValue();
   time2.tm_sec = (Int_t) p->GetLeaf("second")->GetValue();
   Int_t millisecond2 = (Int_t) p->GetLeaf("millisecond")->GetValue();

   printf("Entry %8d date: %02d/%02d/%d %02d:%02d:%02d.%03d\n", entry2, time2.tm_mon,time2.tm_mday,time2.tm_year+1900, time2.tm_hour,time2.tm_min,time2.tm_sec,millisecond2);

   // time difference

   time_t abs_time1 = std::mktime(&time1);   // NB: mktime may change the argument: make sure you set time1.tm_isdst = -1;
   time_t abs_time2 = std::mktime(&time2);
   Double_t dtime = (abs_time2 - abs_time1) + 0.001*(millisecond2 - millisecond1);
   cout<< "time difference is " << dtime << " s" <<endl;
}

TGraph* average(Int_t nchan, Double_t xmin, Double_t xmax, TGraph* gr);
TGraph* average(Int_t nchan) {return average(nchan,0,0,0);}
TGraph* average(Int_t nchan, TGraph* gr) {return average(nchan,0,0,gr);}
TGraph* average(Int_t nchan, Int_t xmin, TGraph* gr) {return average(nchan,xmin,0,gr);}
// TGraph* average(Int_t nchan, Double_t xmin=0, Double_t xmax=0, TGraph* gr=0)
TGraph* average(Int_t nchan, Double_t xmin, Double_t xmax, TGraph* gr)
{
   if (!gr) gr = gtemp();
   if (!gr) {
      cout<< "Could not find TGraph object" <<endl;
      return 0;
   }

   if (xmin == 0 and xmax == 0) {
      xmin = gr->GetX()[0];
      xmax = gr->GetX()[gr->GetN()-1];
   }
   if (xmax < xmin) {
      xmax = gr->GetX()[gr->GetN()-1];
   }
   if (xmin >= xmax) {
      cout<< "Limit(s) are out of range (" << gr->GetX()[0] << ", " << gr->GetX()[gr->GetN()-1] << ")" <<endl;
      return 0;
   }
   // cout<< "check imin imax" <<endl;

   Int_t imin = -1;
   Int_t imax = -1;
   for (int i=0; i<gr->GetN(); ++i) {
      if (imin < 0 and gr->GetX()[i] >= xmin) imin = i;
      if (gr->GetX()[i] <= xmax) imax = i;
      else break;
   }
   if (imin < 0 or imax < 0) {
      cout<< "Limit(s) are out of range (" << gr->GetX()[0] << ", " << gr->GetX()[gr->GetN()-1] << ")" <<endl;
      return 0;
   }
   cout<< "imin = " << imin << " imax = " << imax <<endl;

   Int_t np_orig = imax - imin + 1;
   TGraph* gra = new TGraph(np_orig % nchan? np_orig/nchan+1: np_orig/nchan);

   gra->SetName(Form("%s_%d", gr->GetName(),nchan));
   gra->SetTitle(Form("average %d chan for %s", nchan, gr->GetTitle()));
   gra->SetMarkerStyle(gr->GetMarkerStyle());
   gra->SetMarkerColor(gr->GetMarkerColor());
   gra->SetLineColor(gr->GetLineColor());

   Int_t np = 0;
   Double_t xsum = 0;
   Double_t ysum = 0;
   Int_t nsum = 0;
   for (int i=imin; i<=imax; ++i) {
      xsum += gr->GetX()[i];
      ysum += gr->GetY()[i];
      nsum++;
      //cout<< i << "\t nsum = " << nsum << " xsum = " << xsum << " ysum = " << ysum <<endl;
      if (nsum == nchan or i == imax) {
         xsum /= nsum;
         ysum /= nsum;
         gra->SetPoint(np, xsum, ysum);
         // cout<< np << "\t x = " << gra->GetX()[np] << " y = " << gra->GetY()[np] <<endl;
         ++np;
         xsum = 0;
         ysum = 0;
         nsum = 0;
      }
   }
   cout<< "np = " << np <<endl;

   new TCanvas;
   // gra->Draw("ap");
   gra->Draw("apl");

   return gra;
}

TGraph* baseline(Int_t nchan, Double_t xmin=0, Double_t xmax=0, TGraph* gr=0)
{
   if (!gr) gr = gtemp();
   if (!gr) {
      cout<< "Could not find TGraph object" <<endl;
      return 0;
   }

   if (xmin == 0 and xmax == 0) {
      xmin = gr->GetX()[0];
      xmax = gr->GetX()[gr->GetN()-1];
   }
   if (xmax < xmin) {
      xmax = gr->GetX()[gr->GetN()-1];
   }
   if (xmin >= xmax) {
      cout<< "Limit(s) are out of range (" << gr->GetX()[0] << ", " << gr->GetX()[gr->GetN()-1] << ")" <<endl;
      return 0;
   }
   // cout<< "check imin imax" <<endl;

   Int_t imin = -1;
   Int_t imax = -1;
   for (int i=0; i<gr->GetN(); ++i) {
      if (imin < 0 and gr->GetX()[i] >= xmin) imin = i;
      if (gr->GetX()[i] <= xmax) imax = i;
      else break;
   }
   if (imin < 0 or imax < 0) {
      cout<< "Limit(s) are out of range (" << gr->GetX()[0] << ", " << gr->GetX()[gr->GetN()-1] << ")" <<endl;
      return 0;
   }
   cout<< "imin = " << imin << " imax = " << imax <<endl;

   Int_t np_orig = imax - imin + 1;
   TGraph* gra = new TGraph(np_orig % nchan? np_orig/nchan+1: np_orig/nchan);

   gra->SetName(Form("%s_%d", gr->GetName(),nchan));
   gra->SetTitle(Form("average %d chan for %s", nchan, gr->GetTitle()));
   gra->SetMarkerStyle(gr->GetMarkerStyle());
   gra->SetMarkerColor(gr->GetMarkerColor());
   gra->SetLineColor(gr->GetLineColor());

   TGraph* grs = (TGraph*) gra->Clone(Form("sigma_%s_%d", gr->GetName(),nchan));
   grs->SetTitle(Form("sigma of average %d chan for %s", nchan, gr->GetTitle()));

   Int_t np = 0;
   Double_t xsum = 0;

   Double_t mean = 0;
   Double_t M2 = 0;
   Int_t nsum = 0;
   for (int i=imin; i<=imax; ++i) {
      xsum += gr->GetX()[i];
      Double_t value = gr->GetY()[i];
      nsum++;
      Double_t delta = value - mean;
      mean += delta / nsum;
      M2 += delta*(value - mean);
      if (nsum == nchan or i == imax) {
         xsum /= nsum;
         gra->SetPoint(np, xsum, mean);
         Double_t variance = nsum > 1? M2/(nsum - 1): 0;    // NB: variance can be 0
         Double_t sigma = TMath::Sqrt(variance);
         grs->SetPoint(np, xsum, sigma);
         ++np;
         xsum = 0;
         mean = 0;
         M2 = 0;
         nsum = 0;
      }
   }
   cout<< "np = " << np <<endl;

   new TCanvas;
   gra->Draw("apl");

   new TCanvas;
   grs->Draw("apl");

   return gra;
}
