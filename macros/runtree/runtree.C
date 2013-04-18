// zatserkl@fnal.gov

/*
  NB: inverse quotes ``!
  g++ `$ROOTSYS/bin/root-config --cflags --glibs` -o runtree runtree.C
*/
#include <TRint.h>          // the only include you need

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

//
// namespace Tree
//
namespace Tree
{
   Int_t run;
   Int_t evt;
   Int_t t1;
   Int_t t2;
   Int_t t3;
   Int_t a1;
   Int_t a2;
   Int_t a3;
   
   void clear()
   {
      run      = 0;
      evt      = 0;
      t1       = 0;
      t2       = 0;
      t3       = 0;
      a1       = 0;
      a2       = 0;
      a3       = 0;
   }
   
   void book(TTree* tree) {
      tree->Branch("run",         &run,        "run/I");
      tree->Branch("evt",         &evt,        "evt/I");
      tree->Branch("t1",          &t1,         "t1/I");
      tree->Branch("t2",          &t2,         "t2/I");
      tree->Branch("t3",          &t3,         "t3/I");
      tree->Branch("a1",          &a1,         "a1/I");
      tree->Branch("a2",          &a2,         "a2/I");
      tree->Branch("a3",          &a3,         "a3/I");
   }
   void connect(TTree* tree)                                // need for event-by-event analysis
   {   
      // connects tree buffers with variables to use for event-by-event analysis
      tree->SetBranchAddress("run",           &run);
      tree->SetBranchAddress("evt",           &evt);
      tree->SetBranchAddress("t1",            &t1);
      tree->SetBranchAddress("t2",            &t2);
      tree->SetBranchAddress("t3",            &t3);
      tree->SetBranchAddress("a1",            &a1);
      tree->SetBranchAddress("a2",            &a2);
      tree->SetBranchAddress("a3",            &a3);
   }
}  // namespace Tree

int runNum(std::string fname)
{
  // const std::string ifname = "/home/zatserkl/work/time/2009_run/run172.out";
  // const std::string run = "run";
  // unsigned pos = ifname.rfind(run);
  // cout<< "runnum(ifname.substr(pos+run.size())) = " << runnum(ifname.substr(pos+run.size())) <<endl;

  std::stringstream ss(fname);
  int n = 0;
  ss >> n;
  return n;
}

void addrun(TTree* tree, const std::string& ifname)
{
   const std::string run = "run";
   unsigned pos = ifname.rfind(run);
   int runnum = runNum(ifname.substr(pos+run.size()));
   //cout<< "runnum(ifname.substr(pos+run.size())) = " << runnum <<endl;
   cout<< "runnum = " << runnum <<endl;

   std::ifstream ifile(ifname.c_str());
   if (!ifile) {
      cout<< "***ERROR addrun: file not found: " << ifname <<endl;
      return;
   }

   Tree::clear();
   Tree::run = runnum;
   --Tree::evt;
   while (ifile
         >> Tree::t1
         >> Tree::t2
         >> Tree::t3
         >> Tree::a1
         >> Tree::a2
         >> Tree::a3
         )
   {
      ++Tree::evt;
      tree->Fill();
   }
   ifile.close();
}

TTree* tree_write(const char* runlist_fname, const char* ofname)
{
   //bool debug = true;
   bool debug = false;
   if (debug) cout<< "debug is on" <<endl;  // to use identificator debug
   
   // create tree
   TFile* ofile = TFile::Open(ofname,"recreate");
   TTree* tree = new TTree("t","tree example");
   Tree::book(tree);
   
   // loop over list of run files
   std::ifstream runlist(runlist_fname);
   if (!runlist) {
      cout<< "***ERROR tree_write: file not found: " << runlist_fname <<endl;
      return 0;
   }

   std::string ifname;
   while (runlist >> ifname) {
      addrun(tree, ifname);
   }
   runlist.close();
   
   ofile = tree->GetCurrentFile();  // update current file pointer
   ofile->Write();
   cout<< "Wrote tree " << tree->GetName() << " in file " << ofile->GetName() <<endl;
   return tree;
}

// void tree_draw(const char* ifname="tree_example.root")
// {
//    // uses TTree::Draw
//    
//    //bool debug = true;
//    bool debug = false;
//    if (debug) cout<< "debug is on" <<endl;  // to use identificator debug
//     
//    // arrange tree
//    TFile* ifile = TFile::Open(ifname);
//    TTree* tree  = (TTree*) ifile->Get("t");
// 
//    TH1F* h = new TH1F("h", "h", 100, 0, 100);
//    if (debug) cout<< "Fill histogram " << h->GetName() <<endl;  // to use h
//    
//    new TCanvas;
//    tree->Draw("a>>h", "");
// }
// 
// void tree_read(const char* ifname="tree_example.root")
// {
//    // calls Tree::connect(tree) to connect tree buffers with variables for event-by-event analysis
//    
//    //bool debug = true;
//    bool debug = false;
//    if (debug) cout<< "debug is on" <<endl;  // to use identificator debug
//     
//    // arrange tree
//    TFile* ifile = TFile::Open(ifname);
//    if (!ifile) {
//       cout<< "File not found: " << ifname <<endl;
//       return;
//    }
//    TTree* tree  = (TTree*) ifile->Get("t");
//    Tree::connect(tree);
// 
//    TH1F* h = new TH1F("h", "h", 100, 0, 100);
//    if (debug) cout<< "h name = " << h->GetName() <<endl;    // to use debug and h
//  
//    for (int ievent=0; ievent<tree->GetEntries(); ++ievent)  // loop over events
//    {
//       tree->GetEntry(ievent);
//       h->Fill(Tree::a);
//    }
//    
//    new TCanvas;
//    h->Draw();
// }

#if !defined(__CINT__) && !defined(__MAKECINT__)
int main(int argc, char *argv[])                   // w/o int will be warning ~ "no type"
{
  TRint* theApp = new TRint("Rint", &argc, argv, 0, 0, 1);  // do not show splash screen

  // code starts here

  const char runlist_fname[] = "runtree.lst";
  const char ofname[] = "runtree.root";

  tree_write(runlist_fname, ofname);

  // TH1F* h = new TH1F("h", "h", 100, -3, 3);

  // h->FillRandom("gaus", 10000);
  // 
  // cout<< "h->GetMean() = " << h->GetMean() << " h->GetRMS() = " << h->GetRMS() <<endl;

  // TCanvas* can = new TCanvas();
  // 
  // h->Draw();

  // work in command line mode
  theApp->Run();
  delete theApp;
}
#endif
