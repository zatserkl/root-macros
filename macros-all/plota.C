// Andriy Zatserklyaniy <zatserkl@fnal.gov>

#include <TROOT.h>
#include <TRint.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TEnv.h>
#include <TMath.h>
#include <TF1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TEventList.h>
#include <TBrowser.h>
#include <TPaveStats.h>
#include <TCanvas.h>
#include <TList.h>
#include <TColor.h>
#include <TTimer.h>
#include <TLatex.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <cstdio>

using std::cout;                using std::endl;

Long64_t plot(TTree* tree, const char* varexp, const char* selection="", Option_t* option="", Long64_t nentries=1000000000, Long64_t firstentry=0);
Long64_t plot(TTree* tree, const char* varexp, const TCut& selection, Option_t* option="", Long64_t nentries=1000000000, Long64_t firstentry=0);
Long64_t project(TTree* tree, const char* hname, const char* varexp, const char* selection="", Option_t* option="", Long64_t nentries=1000000000, Long64_t firstentry=0);
Long64_t plota(TTree* tree, const char* varexp, const char* selection="", Option_t* option="", Long64_t nentries=1000000000, Long64_t firstentry=0);
Long64_t plota(TTree* tree, const char* varexp, const TCut& selection, Option_t* option="", Long64_t nentries=1000000000, Long64_t firstentry=0);
Long64_t projecta(TTree* tree, const char* hname, const char* varexp, const char* selection="", Option_t* option="", Long64_t nentries=1000000000, Long64_t firstentry=0);


// 1) TTree::Draw gives wrong result for plotting of single varexp when selection includes arrays
// 2) plot works correctly in this case (single varexp when selection includes arrays)
// 3) plot gives wrong result for plotting of array variable with selection with arrays
// 4) TTree::Draw correctly plots array varexp with selection with arrays (the same true for plota)
//
// functions plot and project: for single (not array) varexp and arrays in selection
// TTree::Draw and TTree::Project give wrong result in this case
// NB: do not apply for array varexp: gives wrong result. Use TTree::Draw and TTree::Project or plota and projecta.
//
Long64_t plot(TTree* tree, const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   TEventList* elist_curr = tree->GetEventList();
   // create list of events for given selection
   tree->Draw(">>elist_temp", selection, option, nentries, firstentry);
   TEventList* elist_temp = (TEventList*) gDirectory->Get("elist_temp");
   tree->SetEventList(elist_temp);
   // plot varexp without selection to provide exactly one entry per event
   Long64_t nselected = tree->Draw(varexp, "", option, nentries, firstentry);
   tree->SetEventList(elist_curr);
   elist_temp->SetDirectory(0);
   delete elist_temp;
   return nselected;
}
Long64_t plot(TTree* tree, const char* varexp, const TCut& selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   TEventList* elist_curr = tree->GetEventList();
   tree->Draw(">>elist_temp", selection, option, nentries, firstentry);
   TEventList* elist_temp = (TEventList*) gDirectory->Get("elist_temp");
   tree->SetEventList(elist_temp);
   Long64_t nselected = tree->Draw(varexp, "", option, nentries, firstentry);
   tree->SetEventList(elist_curr);
   elist_temp->SetDirectory(0);
   delete elist_temp;
   return nselected;
}
Long64_t project(TTree* tree, const char* hname, const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   TEventList* elist_curr = tree->GetEventList();
   tree->Draw(">>elist_temp", selection, option, nentries, firstentry);
   TEventList* elist_temp = (TEventList*) gDirectory->Get("elist_temp");
   tree->SetEventList(elist_temp);
   Long64_t nselected = tree->Project(hname, varexp, "", option, nentries, firstentry);
   tree->SetEventList(elist_curr);
   elist_temp->SetDirectory(0);
   delete elist_temp;
   return nselected;
}
//-- NB: there is no TTree::Project with TCut selection
//
// functions plota and projecta: for array varexp and arrays in selection. Give the same result as TTree::Draw and TTree::Project
// NB: do not apply for single (not array) varexp: gives wrong result. Use plot and project.
//
Long64_t plota(TTree* tree, const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   // array version of plot: apply for array varexp. Gives the same result as just TTree::Draw
   TEventList* elist_curr = tree->GetEventList();
   // create list of events for given selection
   tree->Draw(">>elist_temp", selection, option, nentries, firstentry);
   TEventList* elist_temp = (TEventList*) gDirectory->Get("elist_temp");
   tree->SetEventList(elist_temp);
   //-- substitute "" by selection in 'Long64_t nselected = tree->Draw(varexp, "", option, nentries, firstentry);'
   Long64_t nselected = tree->Draw(varexp, selection, option, nentries, firstentry);
   tree->SetEventList(elist_curr);
   elist_temp->SetDirectory(0);
   delete elist_temp;
   return nselected;
}
Long64_t plota(TTree* tree, const char* varexp, const TCut& selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   // array version of plot: apply for array varexp. Gives the same result as just TTree::Draw
   TEventList* elist_curr = tree->GetEventList();
   tree->Draw(">>elist_temp", selection, option, nentries, firstentry);
   TEventList* elist_temp = (TEventList*) gDirectory->Get("elist_temp");
   tree->SetEventList(elist_temp);
   //-- substitute "" by selection in 'Long64_t nselected = tree->Draw(varexp, "", option, nentries, firstentry);'
   Long64_t nselected = tree->Draw(varexp, selection, option, nentries, firstentry);
   tree->SetEventList(elist_curr);
   elist_temp->SetDirectory(0);
   delete elist_temp;
   return nselected;
}
Long64_t projecta(TTree* tree, const char* hname, const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   // array version of project: apply for array varexp. Gives the same result as just TTree::Project
   TEventList* elist_curr = tree->GetEventList();
   tree->Draw(">>elist_temp", selection, option, nentries, firstentry);
   TEventList* elist_temp = (TEventList*) gDirectory->Get("elist_temp");
   tree->SetEventList(elist_temp);
   //-- substitute "" by selection in 'Long64_t nselected = tree->Project(hname, varexp, "", option, nentries, firstentry);'
   Long64_t nselected = tree->Project(hname, varexp, selection, option, nentries, firstentry);
   tree->SetEventList(elist_curr);
   elist_temp->SetDirectory(0);
   delete elist_temp;
   return nselected;
}
