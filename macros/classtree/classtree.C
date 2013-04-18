/*
// NB: to write file ACLiC is mandatory:      .L classtree.C+
//     to read file is enough just load it:   .L classtree.C
.q
root -l
.L classtree.C+

classtree()

t->Draw("vccu")
t->Draw("disk.vpc", "d==1")
t->Draw("disk.rtd")
t->Draw("disk.rtd", "d==1")
t->Draw("disk.rtd", "d==2 && bp==12")
*/

//#ifndef classtree_h
//#define classtree_h

#include "classtree.h"

#include <TRint.h>        // to compile as standalone C++ program

#include <TROOT.h>        // to use gPad, etc.

#include <TFile.h>
#include <TTree.h>

#include <TH1.h>
#include <TCanvas.h>

#include <iostream>
#include <cstdarg>

using std::cout;        using std::endl;

// class Disk: public TObject {
// public:
//    Float_t  vpc;        // Vpc for this disk
//    Float_t  rtd[8];     // rtd
//    Int_t    bp[8];      // blade/panel for rtd
//    Disk(): TObject() {
//       // set blades and panels labels
//       bp[0]   = 11;
//       bp[1]   = 12;
//       bp[2]   = 61;
//       bp[3]   = 62;
//       bp[4]   = 71;
//       bp[5]   = 72;
//       bp[6]   = 121;
//       bp[7]   = 122;
//       clear();
//    }
//    void clear() {
//       vpc = 0;
//       for (int i=0; i<8; ++i) {
//          rtd[i] = 0;
//       }
//    }
//    ClassDef(Disk, 1)
// };
// 
// class HCyl: public TObject {
// public:
//    Float_t  vccu;       // Vccu for the whole detector
//    Disk     disk[2];    //
//    Int_t    d[2];       // disk #
//    HCyl(): TObject() {
//       // set disk labels
//       d[0] = 1;
//       d[1] = 2;
//       clear();
//    }
//    void clear() {
//       vccu = 0;
//       for (int i=0; i<2; ++i) {
//          disk[i].clear();
//       }
//    }
//    ClassDef(HCyl, 1)
// };

ClassImp(Disk);
ClassImp(HCyl);

void classtree()
{
  const char ofname[] = "classtree.root";
  
  TFile* ofile = TFile::Open(ofname, "recreate");
  TTree* tree = new TTree("t", "test tree");
  
  HCyl* hcyl = new HCyl();
  
  tree->Branch("hc", "HCyl", &hcyl);
  
  // event #0
  hcyl->clear();
    hcyl->vccu = 0.7;
    // disk 1
    hcyl->disk[0].vpc = 2.5;
    hcyl->disk[0].rtd[0]  = -9.95;
    hcyl->disk[0].rtd[1]  = -8.66;
    hcyl->disk[0].rtd[2]  = -7.64;
    hcyl->disk[0].rtd[3]  = -8.98;
    hcyl->disk[0].rtd[4]  = -8.30;
    hcyl->disk[0].rtd[5]  = -9.00;
    hcyl->disk[0].rtd[6]  = -9.27;
    hcyl->disk[0].rtd[7]  = -8.15;
    // disk 2
    hcyl->disk[1].vpc = 2.7;
    hcyl->disk[1].rtd[0]  = -10.32;
    hcyl->disk[1].rtd[1]  = -9.83;
    hcyl->disk[1].rtd[2]  = -8.71;
    hcyl->disk[1].rtd[3]  = -9.38;
    hcyl->disk[1].rtd[4]  = -8.88;
    hcyl->disk[1].rtd[5]  = -9.76;
    hcyl->disk[1].rtd[6]  = -9.86;
    hcyl->disk[1].rtd[7]  = -9.63;
  tree->Fill();             
                            
  // event #1                
  hcyl->clear();            
    hcyl->vccu = 0.8;
    // disk 1
    hcyl->disk[0].vpc = 2.5;
    hcyl->disk[0].rtd[0]  = -9.92;
    hcyl->disk[0].rtd[1]  = -8.66;
    hcyl->disk[0].rtd[2]  = -7.67;
    hcyl->disk[0].rtd[3]  = -8.98;
    hcyl->disk[0].rtd[4]  = -8.30;
    hcyl->disk[0].rtd[5]  = -9.00;
    hcyl->disk[0].rtd[6]  = -9.27;
    hcyl->disk[0].rtd[7]  = -8.17;
    // disk 2
    hcyl->disk[1].vpc = 2.6;
    hcyl->disk[1].rtd[0]  = -10.32;
    hcyl->disk[1].rtd[1]  = -9.86;
    hcyl->disk[1].rtd[2]  = -8.71;
    hcyl->disk[1].rtd[3]  = -9.38;
    hcyl->disk[1].rtd[4]  = -8.88;
    hcyl->disk[1].rtd[5]  = -9.76;
    hcyl->disk[1].rtd[6]  = -9.86;
    hcyl->disk[1].rtd[7]  = -9.63;
  tree->Fill();             
                            
  // event #2                
  hcyl->clear();            
    hcyl->vccu = 0.9;
    // disk 1
    hcyl->disk[0].vpc = 2.6;
    hcyl->disk[0].rtd[0]  = -9.92;
    hcyl->disk[0].rtd[1]  = -8.66;
    hcyl->disk[0].rtd[2]  = -7.64;
    hcyl->disk[0].rtd[3]  = -8.98;
    hcyl->disk[0].rtd[4]  = -8.30;
    hcyl->disk[0].rtd[5]  = -9.00;
    hcyl->disk[0].rtd[6]  = -9.27;
    hcyl->disk[0].rtd[7]  = -8.17;
    // disk 2
    hcyl->disk[1].vpc = 2.5;
    hcyl->disk[1].rtd[0]  = -10.32;
    hcyl->disk[1].rtd[1]  = -9.86;
    hcyl->disk[1].rtd[2]  = -8.71;
    hcyl->disk[1].rtd[3]  = -9.38;
    hcyl->disk[1].rtd[4]  = -8.88;
    hcyl->disk[1].rtd[5]  = -9.76;
    hcyl->disk[1].rtd[6]  = -9.86;
    hcyl->disk[1].rtd[7]  = -9.63;
  tree->Fill();             
                            
  ofile->Write();           

  new TCanvas;
  tree->Draw("vccu");

  new TCanvas;
  tree->Draw("disk.vpc", "d==1");

  new TCanvas;
  tree->Draw("disk.rtd");

  new TCanvas;
  tree->Draw("disk.rtd", "d==1");

  new TCanvas;
  tree->Draw("disk.rtd", "d==2 && bp==12");
}

/*
rootcint -f classtree_dict.cxx -c classtree.h classtree_linkdef.h
g++ `root-config --cflags --glibs` -o classtree classtree.C classtree_dict.cxx
*/

#if !defined(__CINT__) && !defined(__MAKECINT__)
int main(int argc, char* argv[])
{
  // argc = 1;  // if you need to use command line arguments for your needs
  TRint* theApp = new TRint("Rint", &argc, argv, 0, 0, 1);  // do not show splash screen and print ROOT version
  
  classtree();
  
  // work in command line mode
  theApp->Run();
  delete theApp;
}
#endif // #ifndef __CINT__

