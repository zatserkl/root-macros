#include <TROOT.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>

#include <iostream>

using std::cout;     using std::endl;

// demo version, see below the main one
Double_t pol1fast(const Float_t x[], const Float_t y[], Int_t ifirst, Int_t np, Double_t& am, Double_t& bm) 
{
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
   Double_t W = np;     // if weights are different w = 1/(dy*dy) and W = sum(w)

   for (int i=0; i<np; ++i)
   {
      Int_t icurr = ifirst+i;
      wx += x[icurr];
      wx2 += x[icurr]*x[icurr];
      wy += y[icurr];
      wxy += x[icurr]*y[icurr];
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
      chi2 += r*r;
   }

   Double_t chi2ndf = np > 2? chi2/(np - 2): 0;
   return chi2ndf;
}

// demo version, see below the main one
Double_t pol1fast(Int_t np, const Float_t x[], const Float_t y[], Double_t& am, Double_t& bm) 
{
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
   Double_t W = np;     // if weights are different w = 1/(dy*dy) and W = sum(w)

   for (int i=0; i<np; ++i)
   {
      wx += x[i];
      wx2 += x[i]*x[i];
      wy += y[i];
      wxy += x[i]*y[i];
   }

   Double_t Discriminant = W*wx2 - wx*wx;
   if (TMath::Abs(Discriminant) < eps) return huge;   // should not be happen, actually

   am = (wy*wx2 - wxy*wx) / Discriminant;    // y = am + bm*x
   bm = (W*wxy - wx*wy) / Discriminant;

   chi2 = 0;
   for (int i=0; i<np; ++i)
   {
      Double_t yfit = am + bm*x[i];
      Double_t r = y[i] - yfit;
      chi2 += r*r;
   }

   Double_t chi2ndf = np > 2? chi2/(np - 2): 0;
   return chi2ndf;
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

///////////////////////// linfit ///////////////////////////

void linfit()
{
   gStyle->SetOptFit();

   Double_t x2[1000], y2[1000], ex2[1000], ey2[1000];
   Int_t n2 = 10;

   TRandom3 random;
   Double_t x0, delta;
   Double_t dy = 10;
   Double_t sigma = 10;

   // parameters of line y = k*x + b
   Double_t k;
   Double_t b;
   // errors on parameters
   Double_t dk;
   Double_t db;
   // intersection with x-axis
   Double_t xint;
   Double_t dx;

   TF1 *fitfunction = 0;

   delta = 100./n2;
   x0 = 0 - delta;
   for (int i=0; i<n2; ++i) {
      x0 += delta;
      x2[i] = x0;
      y2[i] = random.Gaus(x2[i], sigma);
      ex2[i] = 0;
      ey2[i] = dy;         // all weights are equal
      ey2[i] = TMath::Sqrt(y2[i]);
   }

   TGraphErrors *gr2 = new TGraphErrors(n2, x2,y2, ex2, ey2);
   gr2->SetNameTitle("gr2", "gr2");
   gr2->SetMarkerStyle(20);
   gr2->SetMarkerColor(4);
   gr2->SetLineColor(4);

   new TCanvas;
   gr2->Draw("ap");
   gr2->Fit("pol1");

   fitfunction = gr2->GetFunction("pol1");
   b = fitfunction->GetParameter(0);
   db = fitfunction->GetParError(0);
   k = fitfunction->GetParameter(1);
   dk = fitfunction->GetParError(1);
   xint = -b / k;
   dx = TMath::Abs(b/k) * TMath::Sqrt( (db/b)*(db/b) + (dk/k)*(dk/k) );
   gr2->SetTitle(Form("intersection with x-axis: %0.3f #pm %0.3f", xint,dx));
   cout<< "--> intersection with x-axis: " << xint << " +- " << dx <<endl;

   cout<< "--> chi2ndf = " << fitfunction->GetChisquare() / fitfunction->GetNDF();
   cout<<endl<<endl;

   // fast straight line fit y = am + bm*x
   // Source: Kalashnikova, Kozodaev, Detektory Elementarnyx chastic.

   //-- Float_t version of the data
   Float_t xf[1000], yf[1000], eyf[1000];
   Int_t nf = n2;
   for (int i=0; i<nf; ++i) {
      xf[i] = x2[i];
      yf[i] = y2[i];
      eyf[i] = ey2[i];
   }

   Double_t wx = 0;
   Double_t wy = 0;
   Double_t wxy = 0;
   Double_t wx2 = 0;
   Double_t W = 0;

   wx = wy = wxy = wx2 = W = 0;
   for (int i=0; i<gr2->GetN(); ++i)
   {
      Double_t x = gr2->GetX()[i];
      Double_t y = gr2->GetY()[i];
      Double_t yerror = gr2->GetEY()[i];
      Double_t w = 1./(yerror*yerror);       // weights are 1/ey

      W += w;
      wx += w*x;
      wx2 += w*x*x;
      wy += w*y;
      wxy += w*x*y;
   }

   Double_t Discriminant = W*wx2 - wx*wx;
   Double_t am = (wy*wx2 - wxy*wx) / Discriminant;    // y = am + bm*x
   Double_t bm = (W*wxy - wx*wy) / Discriminant;

   Double_t chi2 = 0;
   Double_t chi2ndf = 0;
   for (int i=0; i<gr2->GetN(); ++i)
   {
      Double_t x = gr2->GetX()[i];
      Double_t y = gr2->GetY()[i];
      Double_t yerror = gr2->GetEY()[i];
      Double_t yfit = am + bm*x;
      Double_t r = (y - yfit) / yerror;
      chi2 += r*r;
   }
   chi2ndf = chi2 / (gr2->GetN() - 2);

   cout<< "--> fast straight line fit: am " << am << " bm " << bm << " chi2ndf = " << chi2ndf <<endl;
   //--
   //-- NB: operator
   //-- cout<< "pol1fast(xf,yf,eyf,0,nf, am,bm) = " << pol1fast(xf,yf,eyf,0,nf, am,bm) << " am " << am << " bm " << bm <<endl;
   //-- will be executed from the tail: 
   //-- print __current__ version (before call of pol1fast) of the am and bm first
   //-- and will call pol1fast (and will obtain updated versions of am and bm) after that!!!!!!!!!!!
   //--
   chi2ndf = pol1fast(xf,yf,eyf,0,nf, am,bm);
   cout<< "pol1fast(xf,yf,eyf,0,nf, am,bm) = " << chi2ndf << " am " << am << " bm " << bm <<endl;
   cout<< "for comparizon eyf=0 generates fit with equal weights:" <<endl;
   chi2ndf = pol1fast(xf,yf,0,0,nf, am,bm);
   cout<< "pol1fast(xf,yf,0,0,nf, am,bm) / dy = " << chi2ndf / dy << " am " << am << " bm " << bm <<endl;

   cout<< "all weights are equal to 1" <<endl;
   wx = wy = wxy = wx2 = W = 0;
   for (int i=0; i<gr2->GetN(); ++i)
   {
      Double_t x = gr2->GetX()[i];
      Double_t y = gr2->GetY()[i];
      // Double_t yerror = gr2->GetEY()[i];
      // Double_t w = 1./(yerror*yerror);       // weights are 1/ey
      Double_t w = 1.;

      W += w;
      wx += w*x;
      wx2 += w*x*x;
      wy += w*y;
      wxy += w*x*y;
   }

   Discriminant = W*wx2 - wx*wx;
   am = (wy*wx2 - wxy*wx) / Discriminant;
   bm = (W*wxy - wx*wy) / Discriminant;

   chi2 = 0;
   for (int i=0; i<gr2->GetN(); ++i)
   {
      Double_t x = gr2->GetX()[i];
      Double_t y = gr2->GetY()[i];
      Double_t yerror = gr2->GetEY()[i];
      Double_t yfit = am + bm*x;
      Double_t r = (y - yfit) / yerror;
      chi2 += r*r;
   }
   chi2ndf = chi2 / (gr2->GetN() - 2);

   cout<< "--> fast straight line fit: am " << am << " bm " << bm << " chi2ndf = " << chi2ndf / dy <<endl;

   // optimized

   cout<< "optimized fit with equal weights" <<endl;

   wx = wy = wxy = wx2 = W = 0;
   W = n2;              // the number of points
   for (int i=0; i<n2; ++i)
   {
      Double_t x = x2[i];
      Double_t y = y2[i];
      wx += x;
      wx2 += x*x;
      wy += y;
      wxy += x*y;
   }

   Discriminant = W*wx2 - wx*wx;
   am = (wy*wx2 - wxy*wx) / Discriminant;    // y = am + bm*x
   bm = (W*wxy - wx*wy) / Discriminant;

   chi2 = 0;
   for (int i=0; i<n2; ++i)
   {
      Double_t x = x2[i];
      Double_t y = y2[i];
      Double_t yfit = am + bm*x;
      Double_t r = y - yfit;
      chi2 += r*r;
   }
   chi2ndf = chi2 / (n2-2);
   cout<< "normalize chi2ndf to dy = 10. NB: in current dataset the weights are different!" <<endl;
   chi2ndf /= dy;

   cout<< "--> fast straight line fit: am " << am << " bm " << bm << " chi2ndf = " << chi2ndf <<endl;

   chi2ndf = pol1fast(nf,xf,yf, am,bm);
   cout<< "pol1fast(nf,xf,yf, am,bm) / dy = " << chi2ndf / dy << " am " << am << " bm " << bm <<endl;

   // version optimized for point.C
   chi2ndf = pol1fast(xf,yf,0,nf, am,bm);
   cout<< "pol1fast(xf,yf,0,nf, am,bm) / dy = " << chi2ndf / dy << " am " << am << " bm " << bm <<endl;
}
