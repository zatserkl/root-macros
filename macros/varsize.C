#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TCanvas.h>

#include <iostream>
using std::cout;        using std::endl;

class Store: public TObject {
public:
   static const Int_t dim = 100000;
   Int_t n;
   // Double_t amember[dim];       //[n]     does not work with array declaration
   Double_t *amember;              //[n]     works with pointer declaration
   Store() {
      cout<< "Store::Store: create array dim = " << dim <<endl;
      amember = new Double_t[dim];
   }
   ~Store() {
      cout<< "Store::~Store" <<endl;
      delete[] amember;
   }

   ClassDef(Store, 3);
};

ClassImp(Store);

void varsize()
{
   TFile* ofile = TFile::Open("varsize.root", "recreate");

   Store* store = new Store;
   TTree* tree = new TTree("tree", "Store tree");              // branch with class
   tree->Branch("store", "Store", &store);

   Double_t asimple[Store::dim];                               // simple static array
   Int_t nelements = 0;
   tree->Branch("nelements", &nelements, "nelements/I");
   tree->Branch("asimple", &asimple, "asimple[nelements]/D");  // nelements in static array

   // Int_t nevents = 10000000;
   Int_t nevents = 100000;
   TRandom3 random;

   Double_t mean = 100;
   Double_t sigma = TMath::Sqrt(mean);

   for (int ievt=0; ievt<nevents; ++ievt) {
      if (ievt % (nevents/10) == 0) cout<< "processing event " << ievt << " out of " << nevents <<endl;
      nelements = random.Integer(store->dim/4 - 1);
      store->n = nelements;
      for (int i=0; i<store->n; ++i) {
         Double_t value = random.Gaus(mean, sigma);
         asimple[i] = value;
         store->amember[i] = value;
      }
      tree->Fill();
   }

   ofile->Write();

   cout<< "Finished filling of the tree and writing of the output file" <<endl;

   return;

   cout<< "--- NB: plotting of the tree takes extra memory!" <<endl;

   new TCanvas;
   tree->Draw("asimple", "");

   new TCanvas;
   tree->Draw("amember", "");
}
