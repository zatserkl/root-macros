/*
  NB: inverse quotes ``!
  g++ `$ROOTSYS/bin/root-config --cflags --glibs` root.cpp
  
To build dictionaries see, e.g. xmltree.cpp. I compiled it using
rootcint -f xmltree_dict.cxx -c xmltree.h xmltree_linkdef.h
g++ `$ROOTSYS/bin/root-config --cflags --glibs` -o xmltree xmltree.cpp xmltree_dict.cxx

See details about dictionaries in ROOT User Guide,
chapter "Writing a Graphical User Interface", part "A Standalone Example"

My comments
-----------
I found that 
TRint* theApp = new TRint("Rint", &argc, argv);
is much better than recommended there
TApplication* theApp = new TApplication("App", &argc,argv);
because TRint gives you usual environment and allows to usual work with you stuff.

After your private code run theApp->Run();
If you'll use something like getchar(), you will not see graphics.

That chapter says that if you want to make dictionary for you own class you should use something like

rootcint -f ex2aDict.cxx -c example2a.h ex2aLinkDef.h
g++ `root-config --cflags --glibs` -o example2a example2a.cxx ex2aDict.cxx
*/

#include <TRint.h>          // the only include you need
//#include <TROOT.h>        // to use gPad, etc.

#include <iostream>
#include <cstdarg>

#include <TH1.h>
#include <TCanvas.h>

using std::cout;      using std::endl;

main(int argc, char *argv[])
{
  TRint* theApp = new TRint("Rint", &argc, argv, 0, 0, 1);  // do not show splash screen

  // code starts here

  TH1F* h = new TH1F("h", "h", 100, -3, 3);

  h->FillRandom("gaus", 10000);
  
  cout<< "h->GetMean() = " << h->GetMean() << " h->GetRMS() = " << h->GetRMS() <<endl;

  TCanvas* can = new TCanvas();
  
  h->Draw();

  // work in command line mode
  theApp->Run();
  delete theApp;
}

