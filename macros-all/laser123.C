#include <TH1.h>
#include <TTree.h>
#include <TCanvas.h>

#include <iostream>

using std::cout;     using std::endl;

// declarations of the external functions to pass ACLiC
TTree* sp3(const char* ifname, const char* tree_name, Bool_t plot);
void rightgaus();

void laser123(const char* ifname="123-close-tuneV-30mV.out")
{
   Double_t ped[] = {50, 37, 3};

   TTree* tree = sp3(ifname, "tree", kFALSE);

   TCanvas* can = new TCanvas("can", "can", 1100, 700);
   can->Divide(3,2);

   for (int i=0; i<3; ++i)
   {
      TH1* h = 0;
      can->cd(1+i);
      tree->Draw(Form("t%d",i+1), Form("t%d>0",i+1));
      h = tree->GetHistogram();
      h->Fit("gaus");
      Double_t tsigma = h->GetFunction("gaus")->GetParameter(2);
      rightgaus();
      h->SetTitle(Form("t%d: #sigma = %0.2f channel;channels", i+1, tsigma));

      can->cd(4+i);
      tree->Draw(Form("a%d",i+1), Form("a%d>1",i+1));
      h = tree->GetHistogram();
      h->Fit("gaus");
      rightgaus();
      Double_t amean = h->GetFunction("gaus")->GetParameter(1);
      Double_t asigma = h->GetFunction("gaus")->GetParameter(2);
      Double_t Npe = ((amean - ped[i])/asigma)*((amean - ped[i])/asigma);
      cout<< "--> Fit parameters: amean = " << amean << " asigma = " << asigma <<endl;
      cout<< "Npe = " << Npe <<endl;
      h->SetTitle(Form("a%d: Npe = %0.2e", i+1, Npe));
   }
   can->cd();
}
