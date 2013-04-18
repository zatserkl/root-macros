#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TTimer.h>

#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <vector>

using std::cout;     using std::endl;

// utils.C stuff
TGraph* gtemp(TCanvas* can=0);

TTree* getTree()
{
   TTree* t = (TTree*) gDirectory->Get("p");
   if (!t) t = (TTree*) gDirectory->Get("pulse");
   if (!t) {
      // the pulse tree may be a friend of the ft tree
      TTree* ft = (TTree*) gDirectory->Get("ft");
      if (ft) {
         t = ft->GetFriend("p");
         if (!t) t = ft->GetFriend("pulse");
      }
   }
   return t;
}

void compare(Int_t evt, Int_t chan1=1, Int_t chan2=2, Int_t chan3=0, Int_t chan4=0, Int_t chan5=0, Int_t chan6=0, Int_t chan7=0, Int_t chan8=0, TTree* tree=0, const char* wname="compare")
{
   if (!tree) tree = getTree();
   if (!tree) {
      cout<< "Could not find pulse tree" <<endl;
      return;
   }

   bool newformat = tree->GetLeaf(Form("t%d",(chan1-1)/4 + 1));
   if (!newformat) {
      if (tree->GetLeaf(Form("b%d_t",(chan1-1)/4 + 1))) newformat = false;
      else {
         cout<< "Could not find time branch" <<endl;
         return;
      }
   }

   Int_t channels[8];
   Int_t nchannels = 0;
   if (newformat) {
      if (chan1 > 0 and tree->GetLeaf("c1")) channels[nchannels++] = chan1;
      if (chan2 > 0 and tree->GetLeaf("c2")) channels[nchannels++] = chan2;
      if (chan3 > 0 and tree->GetLeaf("c3")) channels[nchannels++] = chan3;
      if (chan4 > 0 and tree->GetLeaf("c4")) channels[nchannels++] = chan4;
      if (chan5 > 0 and tree->GetLeaf("c5")) channels[nchannels++] = chan5;
      if (chan6 > 0 and tree->GetLeaf("c6")) channels[nchannels++] = chan6;
      if (chan7 > 0 and tree->GetLeaf("c7")) channels[nchannels++] = chan7;
      if (chan8 > 0 and tree->GetLeaf("c8")) channels[nchannels++] = chan8;
   }
   else {
      if (chan1 > 0 and tree->GetLeaf("b1_c1")) channels[nchannels++] = chan1;
      if (chan2 > 0 and tree->GetLeaf("b1_c2")) channels[nchannels++] = chan2;
      if (chan3 > 0 and tree->GetLeaf("b1_c3")) channels[nchannels++] = chan3;
      if (chan4 > 0 and tree->GetLeaf("b1_c4")) channels[nchannels++] = chan4;
      if (chan5 > 0 and tree->GetLeaf("b2_c1")) channels[nchannels++] = chan5;
      if (chan6 > 0 and tree->GetLeaf("b2_c2")) channels[nchannels++] = chan6;
      if (chan7 > 0 and tree->GetLeaf("b2_c3")) channels[nchannels++] = chan7;
      if (chan8 > 0 and tree->GetLeaf("b2_c4")) channels[nchannels++] = chan8;
   }

   Int_t ndivy = 2;                          // 2 divisions along the y-axis
   Int_t ndivx = (nchannels + 1) / ndivy;
   if (nchannels == 1) {
      ndivx = 1;
      ndivy = 1;
   }
   TCanvas* can = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(wname);
   if (!can) {
      can = new TCanvas(wname,wname, 0,0, ndivx*gStyle->GetCanvasDefW(), ndivy*gStyle->GetCanvasDefH());
      if (nchannels > 1) can->Divide(ndivx,ndivy);
   }
   else gPad = can;

   for (int ix=0; ix<ndivx; ++ix) {
      for (int iy=0; iy<ndivy; ++iy)
      {
         can->cd(ix + ndivx*iy + 1);

         Int_t ich = ix*ndivy + iy;
         if (ich < nchannels) {
            Int_t board = (channels[ich]-1)/4 + 1;
            if (newformat) tree->Draw(Form("c%d:t%d",channels[ich],board), Form("Entry$==%d",evt), "lp");
            else {
               Int_t channel = (channels[ich]-1) % 4 + 1;      // old channels from 1 to 4
               tree->Draw(Form("b%d_c%d:b%d_t", board,channel,board), Form("Entry$==%d",evt), "lp");
            }
            gtemp()->SetTitle(tree->GetHistogram()->GetTitle());
            // cout<< "ix = " << ix << " iy = " << iy << " ich = " << ich << " old channel = " << (channels[ich]-1) % 4 + 1 << " board = " << board <<endl;
         }
      }
   }
   can->cd();
   can->Update();
}

void eloop(Int_t evtNo=0, Int_t chan1=1, Int_t chan2=2, Int_t chan3=0, Int_t chan4=0, Int_t chan5=0, Int_t chan6=0, Int_t chan7=0, Int_t chan8=0, TTree* tree=0, const char* wname="compare")
{
   TTimer timer("gSystem->ProcessEvents();",50,kFALSE);  //-- process mouse events every 50 ms

   std::vector<std::string> commands;

   bool first_event = true;

   while (true)
   {
      timer.Start();                                     //-- start processing of mouse events

      std::string line;
      if (first_event) {
         --evtNo;
         first_event = false;
      }
      else {
         cout<< "<CR>: " << evtNo << ", Clone, Save, Quit, command: ";
         std::getline(cin, line);
      }

      timer.Stop();                                      //-- disable processing of mouse events

      // Possible inputs
      //
      // 0) <CR>: line.size() == 0:
      //    show next event
      // 1) number:
      //    interpret as event number: show this event
      // 2) not a number:
      //    can be one character or more than one character
      //       2.1) one character, line.size() == 1:
      //            interpret as menu command: C or S or Q
      //       2.2) more than one character, line.size() > 1:
      //            interpret as a ROOT command, try to execute

      // 0) <CR>
      if (line.size() == 0)
      {
         ++evtNo;
         compare(evtNo,chan1,chan2,chan3,chan4,chan5,chan6,chan7,chan8,tree,wname);
         continue;
      }

      // 1) number
      std::stringstream ss(line);
      Int_t number = -1;
      if (ss >> number and ss.eof()) {
         evtNo = number;
         compare(evtNo,chan1,chan2,chan3,chan4,chan5,chan6,chan7,chan8,tree,wname);
         continue;
      }

      // 2) not a number
      if (line.size() == 1)
      {
         // 2.1) input was just one character: interpret as a menu command
         switch (toupper(line[0])) {
            case 'C':
               gPad->DrawClone();
               break;
            case 'S':
               if (getTree()) gPad->SaveAs(Form("%s_evt_%d.png", getTree()->GetCurrentFile()->GetName(), evtNo));
               else gPad->SaveAs(Form("evt_%d.png", evtNo));
               break;
            case 'Q':
               return;
            default: cout<< "No such key" <<endl;
         }
         continue;
      }
      else {
         // 2.2) input was more than one character: interpret as a ROOT command
         if (line == std::string(".q")) {
            cout<< "To terminate the ROOT session exit the macro first" <<endl;
            break;
         }
         else {
            if (unsigned(line[0]) == 27 and unsigned(line[1]) == 91) {
               // input was some arrow: ESC + '[' + letter A or B or C or D
               for (unsigned i=0; i<commands.size(); ++i) cout<< commands[i] <<endl;
               continue;
            }

            commands.push_back(line);
            gROOT->ProcessLine(line.c_str());
            continue;
         }
      }
   }
}
