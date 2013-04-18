#include <TH1.h>
#include <TF1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TList.h>

void refit()
{
   TH1F* h = new TH1F("h","h", 100,-3,3);
   h->FillRandom("gaus", 100);
   gStyle->SetOptFit();

   h->Fit("gaus");
   h->Fit("gaus", "+", "", -1,1);

   // swap functions such that the stats box shows the last fitted function
   TF1* f2 = (TF1*) h->GetListOfFunctions()->Last();
   h->GetListOfFunctions()->Remove(f2);
   h->GetListOfFunctions()->AddFirst(f2);
}
