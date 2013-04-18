#include <TROOT.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TROOT.h>

#include <iostream>

using std::cout;        using std::endl;

void shift()
{
   TRandom3 rand;

   Double_t x[2000], y[2000];

   Int_t np = 20;
   Double_t mean = 100;
   Double_t sigma = 10;
   Double_t ashift = 1.40;

   for (int i=0; i<np; ++i) x[i] = i;
   for (int i=0; i<np; i+=2) {
      cout<< "i = " << i <<endl;
      y[i] = rand.Gaus(mean, sigma);
      y[i+1] = ashift*y[i];
   }

   TGraph* gr1 = new TGraph(np, x, y);
   gr1->SetNameTitle("gr1", "splitted data");
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerColor(2);

   new TCanvas;
   gr1->Draw("ap");

   // calculate shift factor
   Double_t y1y2 = 0;
   Double_t y2y2 = 0;
   for (int i=0; i<np; i+=2) {
      Double_t y1 = y[i];
      Double_t y2 = y[i+1];
      y1y2 += y1*y2;
      y2y2 += y2*y2;
   }
   Double_t alpha = y2y2/y1y2;
   cout<< "alpha = " << alpha <<endl;
   for (int i=0; i<np; i+=2) {
      y[i] *= alpha;
   }

   TGraph* gr2 = new TGraph(np, x, y);
   gr2->SetNameTitle("gr2", "shifted data");
   gr2->SetMarkerStyle(20);
   gr2->SetMarkerColor(4);

   new TCanvas;
   gr2->Draw("ap");
}
