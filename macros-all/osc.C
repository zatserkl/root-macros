#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <iostream>

using std::cout;     using std::endl;

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
      TClonesArray* oscCh;                         //-> array of OscChannel

   public:
      void clear() {
         evt = 0;
         oscCh->Clear();
      }
      OscEvent(): TObject() {
         oscCh = new TClonesArray("OscChannel");
         clear();
      }
      ~OscEvent() {
         delete oscCh;   oscCh = 0;
      }
      Int_t nCh() const {return oscCh->GetEntries();}
      const OscChannel* oscChannel(int ich) const {return (const OscChannel*) oscCh->At(ich);}
      ClassDef(OscEvent, 1)
};

#ifdef __MAKECINT__
#pragma link C++ class OscChannel+;
#pragma link C++ class OscEvent+;
#endif

ClassImp(OscChannel);
ClassImp(OscEvent);

/*
root -l osc.C+
t->Draw("y:x","ch==1")
*/

TTree* osc()
{
   TFile* ofile = TFile::Open("osc.root", "recreate");
   TTree* tree = new TTree("t", "DRS4 tree");

   cout<< "-- create an instance of OscEvent to be used for the Branch" <<endl;
   OscEvent* oscEvent = new OscEvent;
   cout<< "-- tree->Branch" <<endl;
   tree->Branch("oscEvent","OscEvent",&oscEvent);

   for (int ievent=0; ievent<3; ++ievent)
   {
      cout<< "\n-- event " << ievent <<endl;

      cout<< "-- oscEvent->clear()" <<endl;
      oscEvent->clear();
      oscEvent->evt = 1;

      for (int ichan=0; ichan<2; ++ichan)
      {
         cout<< "-- create an instance of OscChannel using operator new with placement" <<endl;

         // Int_t n = oscEvent->oscCh->GetEntries();
         // new ((*oscEvent->oscCh)[n]) OscChannel;
         // OscChannel* oscChannel = (OscChannel*) oscEvent->oscCh->At(n);
         OscChannel* oscChannel = new ((*oscEvent->oscCh)[oscEvent->oscCh->GetEntries()]) OscChannel;

         oscChannel->ch = ichan+1;
         for (int i=0; i<1024; ++i) {
            oscChannel->x[i] = i;
            oscChannel->y[i] = 100*ievent + (ichan+1)*10*i;
         }
      }

      // fill the event
      tree->Fill();
   }

   ofile->Write();

   cout<< "-- delete oscEvent" <<endl;
   tree->ResetBranchAddresses();       // drop branches current objects and allocate new ones
   delete oscEvent;

   return tree;
}
