#define MkProPlain_C
#include "MkProPlain.C"

#include <TFile.h>
#include <TCut.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TEventList.h>
#include <TStyle.h>

#include <iostream>
#include <string>

using std::cout;     using std::endl;

//
//    functions plot and project: for single (not array) varexp and arrays in selection
//    TTree::Draw and TTree::Project give wrong result in this case
//    NB: do not apply for array varexp: gives wrong result. Use TTree::Draw and TTree::Project or plota and projecta.
//
Long64_t plot(TTree* tree, const char* varexp, const char* selection="", Option_t* option="", Long64_t nentries = 1000000000, Long64_t firstentry = 0);
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
Long64_t plot(TTree* tree, const char* varexp, const TCut& selection, Option_t* option="", Long64_t nentries = 1000000000, Long64_t firstentry = 0); // *MENU*
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
Long64_t project(TTree* tree, const char* hname, const char* varexp, const char* selection = "", Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0); // *MENU*
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
//
//    functions plota and projecta: for array varexp and arrays in selection. Give the same result as TTree::Draw and TTree::Project
//    NB: do not apply for single (not array) varexp: gives wrong result. Use plot and project.
//
Long64_t plota(TTree* tree, const char* varexp, const char* selection="", Option_t* option="", Long64_t nentries = 1000000000, Long64_t firstentry = 0);
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
Long64_t plota(TTree* tree, const char* varexp, const TCut& selection, Option_t* option="", Long64_t nentries = 1000000000, Long64_t firstentry = 0); // *MENU*
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
Long64_t projecta(TTree* tree, const char* hname, const char* varexp, const char* selection = "", Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0); // *MENU*
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

/*
   Results on TTree::Draw, plot and plota
   --------------------------------------

1) TTree::Draw gives wrong result for plotting of single varexp when selection includes arrays
2) plot works correctly plots this case (single varexp when selection includes arrays)
3) plot gives wrong result for plotting of array variable with selection with arrays
4) TTree::Draw correctly plots array varexp with selection with arrays(the same true for plota)

bt$ root -l mkplain.C+
*-- Local rootlogon
*-- Execute default rootlogon
*-- Default rootlogon
Local rootlogon.C: Defined const char s[] = "same"
Local rootlogon.C: Defined const char sh[] = "samehist"
Local rootlogon.C: Defined const char sa[] = "sameaxis"
// 
Processing mkplain.C+...
Info in <TUnixSystem::ACLiC>: creating shared library /home/zatserkl/work/bt/./mkplain_C.so

Plot array variable

varexp = mu_mt
TCut selection = mu_nseg==3 && mu_dR>0.5 && mu_pt>20&&mu_pt<200 && mu_mt>0&&mu_mt<200
Loop: N entries = 6447 out of N events passed = 6162. The total number of muons in passed events is 7398 # nmuons_which_passed_selection_in_passed_events = 0
Draw: N entries = 6447
plot: N entries = 7398

Plot single variable (good example: missing ET)

varexp = met_pt
TCut selection = met_pt<200 && mu_nseg==3 && mu_dR>0.5 && mu_pt>20&&mu_pt<200
Loop: N entries = 6142 out of N events passed = 6142. The total number of muons in passed events is 7382 # nmuons_which_passed_selection_in_passed_events = 6430
Draw: N entries = 6430
plot: N entries = 6142

Plot array variable. Compare functions plot and plota

varexp = mu_mt
TCut selection = mu_nseg==3 && mu_dR>0.5 && mu_pt>20&&mu_pt<200 && mu_mt>0&&mu_mt<200
Loop: N entries = 6447 out of N events passed = 6162. The total number of muons in passed events is 7398 # nmuons_which_passed_selection_in_passed_events = 6447
Draw: N entries = 6447
plota: N entries = 6447

Plot single variable (good example: missing ET). Compare functions plot and plota

varexp = met_pt
TCut selection = met_pt<200 && mu_nseg==3 && mu_dR>0.5 && mu_pt>20&&mu_pt<200
Loop: N entries = 6142 out of N events passed = 6142. The total number of muons in passed events is 7382 # nmuons_which_passed_selection_in_passed_events = 6430
Draw: N entries = 6430
plota: N entries = 6430
(class TTree*)0x9dc8860
// .q
*/

/*
root -l mkplain.C+
*/
TTree* mkplain(TTree* tree);
TTree* mkplain(const char* ifname="tmb_plain_tree/data-out.root")
{
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return 0;
   }

   TTree* tree = (TTree*) ifile->Get("t");
   Tree::connect(tree);

   return mkplain(tree);
}

TTree* mkplain(TTree* tree)
{
   std::string varexp;
   TCut selection;
   bool event_passed = false;

   Int_t nmuons_in_passed_events = 0;
   Int_t nmuons_which_passed_selection_in_passed_events = 0;
   Long64_t nEventsLoop = 0;
   Long64_t nEntriesLoop = 0;
   Long64_t nEntriesDraw = 0;
   Long64_t nEntriesPlot = 0;

   gStyle->SetOptStat(1101111);

   ////////////// array variable ///////////////////

   cout<< "\nPlot array variable\n" <<endl;

   varexp = "mu_mt";
   selection = "mu_nseg==3 && mu_dR>0.5 && mu_pt>20&&mu_pt<200 && mu_mt>0&&mu_mt<200";

   cout<< "varexp = " << varexp <<endl;
   cout<< "TCut selection = " << selection.GetTitle() <<endl;

   TH1F* h_mu_mt_loop = new TH1F("h_mu_mt_loop", Form("loop for %s", varexp.c_str()), 200, 0, 2200);
   TH1F* h_mu_mt_Draw = new TH1F("h_mu_mt_Draw", Form("Draw for %s", varexp.c_str()), 200, 0, 2200);
   TH1F* h_mu_mt_plot = new TH1F("h_mu_mt_plot", Form("plot for %s", varexp.c_str()), 200, 0, 2200);
 
   nEventsLoop = nEntriesLoop = 0;
   nmuons_in_passed_events = 0;
   nmuons_which_passed_selection_in_passed_events = 0;
   for (int ievent=0; ievent<tree->GetEntries(); ++ievent)  // loop over events
   {
      tree->GetEntry(ievent);
      event_passed = false;

      for (int imu=0; imu<Tree::nmu; ++imu) {
         if (true
               // TCut selection = "mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200 &&mu_mt<200";
               && Tree::mu_nseg[imu]   == 3
               && Tree::mu_dR[imu]     > 0.5
               && Tree::mu_pt[imu]     > 20
               && Tree::mu_pt[imu]     < 200
               && Tree::mu_mt[imu]     > 0
               && Tree::mu_mt[imu]     < 200
         ) {
            event_passed = true;
            ++nEntriesLoop;
            h_mu_mt_loop->Fill(Tree::mu_mt[imu]);
         }
      }
      if (event_passed) {
         ++nEventsLoop;
         nmuons_in_passed_events += Tree::nmu;
      }
   }
   cout<< "Loop: N entries = " << nEntriesLoop << " out of N events passed = " << nEventsLoop << ". The total number of muons in passed events is " << nmuons_in_passed_events << " # nmuons_which_passed_selection_in_passed_events = " << nmuons_which_passed_selection_in_passed_events <<endl;
   
   new TCanvas;
   gPad->SetLogy();
   h_mu_mt_loop->Draw();

   new TCanvas;
   gPad->SetLogy();
   nEntriesDraw = tree->Draw(Form("%s>>%s",varexp.c_str(),h_mu_mt_Draw->GetName()), selection);
   // Long64_t nEntriesDraw = tree->Draw(varexp.c_str(), selection);
   cout<< "Draw: N entries = " << nEntriesDraw <<endl;

   new TCanvas;
   gPad->SetLogy();
   nEntriesPlot = plot(tree, Form("%s>>%s",varexp.c_str(),h_mu_mt_plot->GetName()), selection);
   cout<< "plot: N entries = " << nEntriesPlot <<endl;

   ////////// single variable (good example: missing ET) ///////////////////

   cout<< "\nPlot single variable (good example: missing ET)\n" <<endl;

   varexp = "met_pt";
   selection = "met_pt<200 && mu_nseg==3 && mu_dR>0.5 && mu_pt>20&&mu_pt<200";

   cout<< "varexp = " << varexp <<endl;
   cout<< "TCut selection = " << selection.GetTitle() <<endl;

   TH1F* h_met_pt_loop = new TH1F("h_met_pt_loop", Form("loop for %s", varexp.c_str()), 200, 0, 2200);
   TH1F* h_met_pt_Draw = new TH1F("h_met_pt_Draw", Form("Draw for %s", varexp.c_str()), 200, 0, 2200);
   TH1F* h_met_pt_plot = new TH1F("h_met_pt_plot", Form("plot for %s", varexp.c_str()), 200, 0, 2200);
 
   nEventsLoop = nEntriesLoop = 0;
   nmuons_in_passed_events = 0;
   nmuons_which_passed_selection_in_passed_events = 0;
   for (int ievent=0; ievent<tree->GetEntries(); ++ievent)  // loop over events
   {
      tree->GetEntry(ievent);
      event_passed = false;

      // selection = "met_pt<200 && mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200";

      //-- cut on single variable
      event_passed = true
         && Tree::met_pt < 200
      ;
      //-- cut on muon array: need at least one muon
      Int_t nmu_passed = 0;
      if (event_passed) {
         event_passed = false;
         for (int imu=0; imu<Tree::nmu; ++imu) {
            if (true
                  && Tree::mu_nseg[imu]   == 3
                  && Tree::mu_dR[imu]     > 0.5
                  && Tree::mu_pt[imu]     > 20
                  && Tree::mu_pt[imu]     < 200
            ) {
               event_passed = true;
               ++nmu_passed;
               // break;      <---- do not break, we are counting muons
            }
         }
      }
      //-- fill histo
      if (event_passed) {
         ++nEventsLoop;
         ++nEntriesLoop;
         nmuons_in_passed_events += Tree::nmu;
         nmuons_which_passed_selection_in_passed_events += nmu_passed;
         h_met_pt_loop->Fill(Tree::met_pt);
      }
   }
   cout<< "Loop: N entries = " << nEntriesLoop << " out of N events passed = " << nEventsLoop << ". The total number of muons in passed events is " << nmuons_in_passed_events << " # nmuons_which_passed_selection_in_passed_events = " << nmuons_which_passed_selection_in_passed_events <<endl;
   
   new TCanvas;
   gPad->SetLogy();
   h_met_pt_loop->Draw();

   new TCanvas;
   gPad->SetLogy();
   nEntriesDraw = tree->Draw(Form("%s>>%s",varexp.c_str(),h_met_pt_Draw->GetName()), selection);
   // Long64_t nEntriesDraw = tree->Draw(varexp.c_str(), selection);
   cout<< "Draw: N entries = " << nEntriesDraw <<endl;

   new TCanvas;
   gPad->SetLogy();
   nEntriesPlot = plot(tree, Form("%s>>%s",varexp.c_str(),h_met_pt_plot->GetName()), selection);
   cout<< "plot: N entries = " << nEntriesPlot <<endl;

   //---------------------------------------------
   //--
   //------------- plota
   //--
   //---------------------------------------------

   ////////////// array variable ///////////////////

   cout<< "\nPlot array variable. Compare functions plot and plota\n" <<endl;

   varexp = "mu_mt";
   selection = "mu_nseg==3 && mu_dR>0.5 && mu_pt>20&&mu_pt<200 && mu_mt>0&&mu_mt<200";

   cout<< "varexp = " << varexp <<endl;
   cout<< "TCut selection = " << selection.GetTitle() <<endl;

   TH1F* h_mu_mt_loop1 = new TH1F("h_mu_mt_loop1", Form("loop for %s", varexp.c_str()), 200, 0, 2200);
   TH1F* h_mu_mt_Draw1 = new TH1F("h_mu_mt_Draw1", Form("Draw for %s", varexp.c_str()), 200, 0, 2200);
   TH1F* h_mu_mt_plota = new TH1F("h_mu_mt_plota", Form("plot for %s", varexp.c_str()), 200, 0, 2200);
 
   nEventsLoop = nEntriesLoop = 0;
   nmuons_in_passed_events = 0;
   nmuons_which_passed_selection_in_passed_events = 0;
   for (int ievent=0; ievent<tree->GetEntries(); ++ievent)  // loop over events
   {
      tree->GetEntry(ievent);
      event_passed = false;

      for (int imu=0; imu<Tree::nmu; ++imu) {
         if (true
               // TCut selection = "mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200 &&mu_mt<200";
               && Tree::mu_nseg[imu]   == 3
               && Tree::mu_dR[imu]     > 0.5
               && Tree::mu_pt[imu]     > 20
               && Tree::mu_pt[imu]     < 200
               && Tree::mu_mt[imu]     > 0
               && Tree::mu_mt[imu]     < 200
         ) {
            event_passed = true;
            ++nEntriesLoop;
            ++nmuons_which_passed_selection_in_passed_events;
            h_mu_mt_loop1->Fill(Tree::mu_mt[imu]);
         }
      }
      if (event_passed) {
         ++nEventsLoop;
         nmuons_in_passed_events += Tree::nmu;
      }
   }
   cout<< "Loop: N entries = " << nEntriesLoop << " out of N events passed = " << nEventsLoop << ". The total number of muons in passed events is " << nmuons_in_passed_events << " # nmuons_which_passed_selection_in_passed_events = " << nmuons_which_passed_selection_in_passed_events <<endl;
   
   new TCanvas;
   gPad->SetLogy();
   h_mu_mt_loop1->Draw();

   new TCanvas;
   gPad->SetLogy();
   nEntriesDraw = tree->Draw(Form("%s>>%s",varexp.c_str(),h_mu_mt_Draw1->GetName()), selection);
   // Long64_t nEntriesDraw = tree->Draw(varexp.c_str(), selection);
   cout<< "Draw: N entries = " << nEntriesDraw <<endl;

   new TCanvas;
   gPad->SetLogy();
   //nEntriesPlot = plot(tree, Form("%s>>%s",varexp.c_str(),h_mu_mt_plota->GetName()), selection);
   nEntriesPlot = plota(tree, Form("%s>>%s",varexp.c_str(),h_mu_mt_plota->GetName()), selection);
   cout<< "plota: N entries = " << nEntriesPlot <<endl;

   ////////////// single variable (good example: missing ET) ///////////////////

   cout<< "\nPlot single variable (good example: missing ET). Compare functions plot and plota\n" <<endl;

   varexp = "met_pt";
   selection = "met_pt<200 && mu_nseg==3 && mu_dR>0.5 && mu_pt>20&&mu_pt<200";

   cout<< "varexp = " << varexp <<endl;
   cout<< "TCut selection = " << selection.GetTitle() <<endl;

   TH1F* h_met_pt_loop1 = new TH1F("h_met_pt_loop1", Form("loop for %s", varexp.c_str()), 200, 0, 2200);
   TH1F* h_met_pt_Draw1 = new TH1F("h_met_pt_Draw1", Form("Draw for %s", varexp.c_str()), 200, 0, 2200);
   TH1F* h_met_pt_plota = new TH1F("h_met_pt_plota", Form("plota for %s", varexp.c_str()), 200, 0, 2200);
 
   nEventsLoop = nEntriesLoop = 0;
   nmuons_in_passed_events = 0;
   nmuons_which_passed_selection_in_passed_events = 0;
   for (int ievent=0; ievent<tree->GetEntries(); ++ievent)  // loop over events
   {
      tree->GetEntry(ievent);
      event_passed = false;

      // selection = "met_pt<200 && mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200";

      //-- cut on single variable
      event_passed = true
         && Tree::met_pt < 200
      ;
      //-- cut on muon array: need at least one muon
      Int_t nmu_passed = 0;
      if (event_passed) {
         event_passed = false;
         for (int imu=0; imu<Tree::nmu; ++imu) {
            if (true
                  && Tree::mu_nseg[imu]   == 3
                  && Tree::mu_dR[imu]     > 0.5
                  && Tree::mu_pt[imu]     > 20
                  && Tree::mu_pt[imu]     < 200
            ) {
               event_passed = true;
               ++nmu_passed;
               // break;      <---- do not break, we are counting muons
            }
         }
      }
      //-- fill histo
      if (event_passed) {
         ++nEventsLoop;
         ++nEntriesLoop;
         nmuons_in_passed_events += Tree::nmu;
         nmuons_which_passed_selection_in_passed_events += nmu_passed;
         h_met_pt_loop1->Fill(Tree::met_pt);
      }
   }
   cout<< "Loop: N entries = " << nEntriesLoop << " out of N events passed = " << nEventsLoop << ". The total number of muons in passed events is " << nmuons_in_passed_events << " # nmuons_which_passed_selection_in_passed_events = " << nmuons_which_passed_selection_in_passed_events <<endl;
   
   new TCanvas;
   gPad->SetLogy();
   h_met_pt_loop1->Draw();

   new TCanvas;
   gPad->SetLogy();
   nEntriesDraw = tree->Draw(Form("%s>>%s",varexp.c_str(),h_met_pt_Draw1->GetName()), selection);
   // Long64_t nEntriesDraw = tree->Draw(varexp.c_str(), selection);
   cout<< "Draw: N entries = " << nEntriesDraw <<endl;

   new TCanvas;
   gPad->SetLogy();
   //nEntriesPlot = plot(tree, Form("%s>>%s",varexp.c_str(),h_met_pt_plota->GetName()), selection);
   nEntriesPlot = plota(tree, Form("%s>>%s",varexp.c_str(),h_met_pt_plota->GetName()), selection);
   cout<< "plota: N entries = " << nEntriesPlot <<endl;

   return tree;

   //-------------------------- BACKUP examples NB: (slightly) obsolete ------------------------------------//

   // ////////////// single variable ///////////////////

   // cout<< "\nPlot single variable\n" <<endl;

   // varexp = "wmunu_mt";
   // selection = "wmunu_mt>20 && mu_pt>20&&mu_pt<200";
   // // selection = "mu_pt>20&&mu_pt<200";

   // cout<< "varexp = " << varexp <<endl;
   // cout<< "TCut selection = " << selection.GetTitle() <<endl;

   // gStyle->SetOptStat(1101111);

   // TH1F* h_wmunu_mt_loop = new TH1F("h_wmunu_mt_loop", Form("loop for %s", varexp.c_str()), 200, 0, 2200);
   // TH1F* h_wmunu_mt_Draw = new TH1F("h_wmunu_mt_Draw", Form("Draw for %s", varexp.c_str()), 200, 0, 2200);
   // TH1F* h_wmunu_mt_plot = new TH1F("h_wmunu_mt_plot", Form("plot for %s", varexp.c_str()), 200, 0, 2200);
 
   // nEventsLoop = nEntriesLoop = 0;
   // for (int ievent=0; ievent<tree->GetEntries(); ++ievent)  // loop over events
   // {
   //    tree->GetEntry(ievent);
   //    event_passed = false;

   //    // TCut selection = "wmunu_mt>20 && mu_pt>20&&mu_pt<200";
   //    event_passed = true
   //       // cut on single variable
   //       && Tree::wmunu_mt       > 20
   //    ;
   //    // cut on muon array
   //    if (event_passed)
   //    {
   //       event_passed = false;
   //       for (int imu=0; imu<Tree::nmu; ++imu) {
   //          if (true
   //                && Tree::mu_pt[imu]     > 20
   //                && Tree::mu_pt[imu]     < 200
   //          ) {
   //             event_passed = true;
   //             break;
   //          }
   //       }
   //    }
   //    if (event_passed) h_wmunu_mt_loop->Fill(Tree::wmunu_mt);
   // }
   // cout<< "Loop: N entries = " << h_wmunu_mt_loop->GetEntries() <<endl;
   // 
   // new TCanvas;
   // gPad->SetLogy();
   // h_wmunu_mt_loop->Draw();

   // new TCanvas;
   // gPad->SetLogy();
   // nEntriesDraw = tree->Draw(Form("%s>>%s",varexp.c_str(),h_wmunu_mt_Draw->GetName()), selection);
   // // Long64_t nEntriesDraw = tree->Draw(varexp.c_str(), selection);
   // cout<< "Draw: N entries = " << nEntriesDraw <<endl;

   // new TCanvas;
   // gPad->SetLogy();
   // nEntriesPlot = plot(tree, Form("%s>>%s",varexp.c_str(),h_wmunu_mt_plot->GetName()), selection);
   // cout<< "plot: N entries = " << nEntriesPlot <<endl;

   // ////////////// single variable 2 (confusing case) ///////////////////

   // cout<< "\nPlot single variable (some confusing case)\n" <<endl;

   // varexp = "zmumu_m";
   // selection = "mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200 &&zmumu_m>0&&zmumu_m<200";

   // cout<< "varexp = " << varexp <<endl;
   // cout<< "TCut selection = " << selection.GetTitle() <<endl;

   // gStyle->SetOptStat(1101111);

   // TH1F* h_zmumu_m_loop = new TH1F("h_zmumu_m_loop", Form("loop for %s", varexp.c_str()), 200, 0, 2200);
   // TH1F* h_zmumu_m_Draw = new TH1F("h_zmumu_m_Draw", Form("Draw for %s", varexp.c_str()), 200, 0, 2200);
   // TH1F* h_zmumu_m_plot = new TH1F("h_zmumu_m_plot", Form("plot for %s", varexp.c_str()), 200, 0, 2200);
 
   // nEventsLoop = nEntriesLoop = 0;
   // for (int ievent=0; ievent<tree->GetEntries(); ++ievent)  // loop over events
   // {
   //    tree->GetEntry(ievent);
   //    event_passed = false;

   //    // selection = "mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200 &&zmumu_m>0&&zmumu_m<200";

   //    //-- cut on single variable
   //    event_passed = true
   //       && Tree::zmumu_m > 0
   //       && Tree::zmumu_m < 200
   //    ;
   //    //-- cut on muon array: need at least one muon
   //    if (event_passed) {
   //       event_passed = false;
   //       for (int imu=0; imu<Tree::nmu; ++imu) {
   //          if (true
   //                && Tree::mu_nseg[imu]   == 3
   //                && Tree::mu_dR[imu]     > 0.5
   //                && Tree::mu_pt[imu]     > 20
   //                && Tree::mu_pt[imu]     < 200
   //          ) {
   //             event_passed = true;
   //             break;
   //          }
   //       }
   //    }
   //    //-- fill histo
   //    if (event_passed) h_zmumu_m_loop->Fill(Tree::zmumu_m);
   // }
   // cout<< "Loop: N entries = " << h_zmumu_m_loop->GetEntries() <<endl;
   // 
   // new TCanvas;
   // gPad->SetLogy();
   // h_zmumu_m_loop->Draw();

   // new TCanvas;
   // gPad->SetLogy();
   // nEntriesDraw = tree->Draw(Form("%s>>%s",varexp.c_str(),h_zmumu_m_Draw->GetName()), selection);
   // // Long64_t nEntriesDraw = tree->Draw(varexp.c_str(), selection);
   // cout<< "Draw: N entries = " << nEntriesDraw <<endl;

   // new TCanvas;
   // gPad->SetLogy();
   // nEntriesPlot = plot(tree, Form("%s>>%s",varexp.c_str(),h_zmumu_m_plot->GetName()), selection);
   // cout<< "plot: N entries = " << nEntriesPlot <<endl;

   //---------------- BACKUP 2: used to catch a bug -------------------//

   // //////////////////////////////////////////////////////////
   // ////////////// single variable 2 (confusing case) ///////////////////

   // cout<< "\nPlot single variable -- code from backup\n" <<endl;

   // varexp = "met_pt";
   // selection = "met_pt<200 && mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200";

   // cout<< "varexp = " << varexp <<endl;
   // cout<< "TCut selection = " << selection.GetTitle() <<endl;

   // TH1F* h_met_pt_1_loop = new TH1F("h_met_pt_1_loop", Form("loop for %s", varexp.c_str()), 200, 0, 2200);
   // TH1F* h_met_pt_1_Draw = new TH1F("h_met_pt_1_Draw", Form("Draw for %s", varexp.c_str()), 200, 0, 2200);
   // TH1F* h_met_pt_1_plot = new TH1F("h_met_pt_1_plot", Form("plot for %s", varexp.c_str()), 200, 0, 2200);
 
   // nEventsLoop = nEntriesLoop = 0;
   // for (int ievent=0; ievent<tree->GetEntries(); ++ievent)  // loop over events
   // {
   //    tree->GetEntry(ievent);
   //    event_passed = false;

   //    // selection = "met_pt<200 && mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200";

   //    //-- cut on single variable
   //    event_passed = true
   //       && Tree::met_pt < 200
   //    ;
   //    //-- cut on muon array: need at least one muon
   //    if (event_passed) {
   //       event_passed = false;
   //       for (int imu=0; imu<Tree::nmu; ++imu) {
   //          if (true
   //                && Tree::mu_nseg[imu]   == 3
   //                && Tree::mu_dR[imu]     > 0.5
   //                && Tree::mu_pt[imu]     > 20
   //                && Tree::mu_pt[imu]     < 200
   //          ) {
   //             event_passed = true;
   //             break;
   //          }
   //       }
   //    }
   //    //-- fill histo
   //    if (event_passed) h_met_pt_1_loop->Fill(Tree::zmumu_m);
   // }
   // cout<< "Loop: N entries = " << h_met_pt_1_loop->GetEntries() <<endl;
   // 
   // new TCanvas;
   // gPad->SetLogy();
   // h_met_pt_1_loop->Draw();

   // new TCanvas;
   // gPad->SetLogy();
   // nEntriesDraw = tree->Draw(Form("%s>>%s",varexp.c_str(),h_met_pt_1_Draw->GetName()), selection);
   // // Long64_t nEntriesDraw = tree->Draw(varexp.c_str(), selection);
   // cout<< "Draw: N entries = " << nEntriesDraw <<endl;

   // new TCanvas;
   // gPad->SetLogy();
   // nEntriesPlot = plot(tree, Form("%s>>%s",varexp.c_str(),h_met_pt_1_plot->GetName()), selection);
   // cout<< "plot: N entries = " << nEntriesPlot <<endl;
   // //////////////////////////////////////////////////////////

}
