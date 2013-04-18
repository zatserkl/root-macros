// zatserkl@fnal.gov

#include <TTree.h>
//
// namespace Tree
//
namespace Tree
{
   Int_t    n;
   Float_t  a;
   
   void clear()
   {
      n        = 0;
      a        = 0;
   }
   
   void book(TTree* tree) {
      tree->Branch("n",          &n,         "n/I");
      tree->Branch("a",          &a,         "a");          // By default equivalent to "a/F"
   }
   void connect(TTree* tree)                                // need for event-by-event analysis
   {   
      // connects tree buffers with variables to use for event-by-event analysis
      tree->SetBranchAddress("n",            &n);
      tree->SetBranchAddress("a",            &a);
   }
}  // namespace Tree

/*
include guards to include above definition of the tree
Code example:
#define tree_example_C
#include "tree_example.C"
*/
#ifndef tree_example_C
#define tree_example_C

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

#include <iostream>

using std::cout;     using std::endl;

TTree* tree_write(const char* ofname="tree_example.root")
{
   //bool debug = true;
   bool debug = false;
   if (debug) cout<< "debug is on" <<endl;  // to use identificator debug
   
   // create tree
   TFile* ofile = TFile::Open(ofname,"recreate");
   TTree* tree = new TTree("t","tree example");
   Tree::book(tree);
   
   // fill tree
   for (int ievent=0; ievent<10; ++ievent) {
      Tree::clear();          // optional
      Tree::n = ievent;       // do not use(?) Tree::n as a loop variable: for(Tree::n=0; Tree::n<10; ++Tree::n)
      Tree::a = 10.*ievent;
      tree->Fill();
   }
   
   ofile = tree->GetCurrentFile();  // update current file pointer
   ofile->Write();
   cout<< "Wrote tree " << tree->GetName() << " in file " << ofile->GetName() <<endl;
   return tree;
}

void tree_draw(const char* ifname="tree_example.root")
{
   // uses TTree::Draw
   
   //bool debug = true;
   bool debug = false;
   if (debug) cout<< "debug is on" <<endl;  // to use identificator debug
    
   // arrange tree
   TFile* ifile = TFile::Open(ifname);
   TTree* tree  = (TTree*) ifile->Get("t");

   TH1F* h = new TH1F("h", "h", 100, 0, 100);
   if (debug) cout<< "Fill histogram " << h->GetName() <<endl;  // to use h
   
   new TCanvas;
   tree->Draw("a>>h", "");
}

void tree_read(const char* ifname="tree_example.root")
{
   // calls Tree::connect(tree) to connect tree buffers with variables for event-by-event analysis
   
   //bool debug = true;
   bool debug = false;
   if (debug) cout<< "debug is on" <<endl;  // to use identificator debug
    
   // arrange tree
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return;
   }
   TTree* tree  = (TTree*) ifile->Get("t");
   Tree::connect(tree);

   TH1F* h = new TH1F("h", "h", 100, 0, 100);
   if (debug) cout<< "h name = " << h->GetName() <<endl;    // to use debug and h
 
   for (int ievent=0; ievent<tree->GetEntries(); ++ievent)  // loop over events
   {
      tree->GetEntry(ievent);
      h->Fill(Tree::a);
   }
   
   new TCanvas;
   h->Draw();
}
#endif
