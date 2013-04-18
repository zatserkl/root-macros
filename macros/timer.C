#include <TH1.h>
#include <TCanvas.h>
#include <TTimer.h>
#include <TROOT.h>
#include <TSystem.h>

#include <iostream>

using std::cout;     using std::endl;

void timer()
{
   TH1F* h = new TH1F("h","h",100,-3,3);
   h->FillRandom("gaus",1000);
   new TCanvas;
   h->Draw();

   TTimer timer("gSystem->ProcessEvents();",50,kFALSE);
   // timer.TurnOn();         // works here but does not work in evt.C
   timer.Start();

   cout<< "<CR> = Quit" <<endl;
   getchar();

   timer.TurnOff();
}
