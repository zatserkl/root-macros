#include <TROOT.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TRandom3.h>
#include <TMath.h>

#include <iostream>

using std::cout;     using std::endl;

/*
   http://root.cern.ch/root/html532/TH1.html

   Histograms with automatic bins

   When an histogram is created with an axis lower limit greater or equal to its upper 
   limit, the SetBuffer is automatically called with an argument fBufferSize equal to 
   fgBufferSize (default value=1000). fgBufferSize may be reset via the static function 
   TH1::SetDefaultBufferSize. The axis limits will be automatically computed when the 
   buffer will be full or when the function BufferEmpty is called
*/

TH1F* histo_autolimits()
{
   TH1::SetDefaultBufferSize(1000000);          // The parameter type is Int_t
   TH1F* h_autolimits = new TH1F("h_autolimits", "xlow >= xup: compute limits automatically", 100, 0, 0);

   // h_autolimits->FillRandom("gaus", 1000);   // Does not work with auto limits

   TRandom3 random(4357);                       // 4357 is the default seed

   // parameters of the Gauss distribution
   Double_t mean = 100.;
   Double_t sigma = TMath::Sqrt(mean);

   Double_t vmin = 0, vmax = 0;
   for (int i=0; i<100; ++i)
   {
      Double_t r = random.Gaus(mean, sigma);
      if (i == 0) {
         vmin = r;
         vmax = r;
      }
      else {
         if (r < vmin) vmin = r;
         if (r > vmax) vmax = r;
      }
      h_autolimits->Fill(r);
   }

   new TCanvas;
   h_autolimits->Draw();

   // The option "L" gives the correct sigma 
   // because of limited statistics
   h_autolimits->Fit("gaus", "L");

   cout<< "vmin = " << vmin << " vmax = " << vmax <<endl;

   //-- Side effect of using of the generator

   // NB: the option "L" gives better result for this histo too
   TH1F* h_bias = new TH1F("h_bias", "Side effect of the random generator", 100, -200, 200);
   h_bias->FillRandom("gaus", 100);

   new TCanvas;
   h_bias->Draw();

   return h_autolimits;
}
