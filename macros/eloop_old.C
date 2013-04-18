#include "drs.C"

#include <TROOT.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TCut.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TTimer.h>
#include <TAxis.h>

#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <vector>

using std::cout;     using std::endl;

// utils.C stuff
TGraph* gtempClone(const char* name=0, const char* title=0);

void evt(TCut cut="", Int_t evtNo=0, Int_t ch1=-1, Int_t ch2=-1, Int_t ch3=-1, Int_t ch4=-1, Int_t ch5=-1, Int_t ch6=-1, Int_t ch7=-1, Int_t ch8=-1)
{
   if (ch1==-1 && ch2==-1 && ch3==-1 && ch4==-1 && ch5==-1 && ch6==-1 && ch7==-1 && ch8==-1) {
      cout<< "Draws:    y vs x for tree with name \"t\"" <<endl;
      cout<< "Usage:    evt(TCut cut=\"\", Int_t evtNo=0, ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8)" <<endl;
      cout<< "Example:  evt(\"x>0&&x<20\", 9, 1,2,3)" <<endl;
      return;
   }

   TTree *t = (TTree*) gDirectory->Get("t");
   if (!t) {
      cout<< "No tree \"t\" was found in the current directory " << gDirectory->GetName() <<endl;
      return;
   }

   std::vector<Int_t> channels;
   if (ch1 != -1) channels.push_back(ch1);
   if (ch2 != -1) channels.push_back(ch2);
   if (ch3 != -1) channels.push_back(ch3);
   if (ch4 != -1) channels.push_back(ch4);
   if (ch5 != -1) channels.push_back(ch5);
   if (ch6 != -1) channels.push_back(ch6);
   if (ch7 != -1) channels.push_back(ch7);
   if (ch8 != -1) channels.push_back(ch8);

   const Int_t colors[] = {2, 4, 6, 8, 5, 7, 46, 1};

   std::vector<TGraph*> graphs;

   Double_t xmin = 0;
   Double_t xmax = 0;
   Double_t ymin = 0;
   Double_t ymax = 0;
   for (unsigned ichannel=0; ichannel<channels.size(); ++ichannel) {
      Int_t ch = channels[ichannel];
      TCut cutevt = TCut(Form("Entry$==%d&&ch==%d",evtNo,ch)) + cut;
      t->Draw("-y:x", cutevt, "goff");

      TGraph* gr = new TGraph(t->GetSelectedRows(), t->GetV2(), t->GetV1());
      gr->SetNameTitle(Form("gr_evt_%d_ch_%d",evtNo,ch), Form("evt %d ch %d",evtNo,ch));
      gr->SetMarkerStyle(7);
      gr->SetMarkerColor(colors[(ch - 1) % sizeof(colors)]);
      graphs.push_back(gr);

      if (xmin == 0 && xmax == 0) {
         xmin = gr->GetXaxis()->GetXmin();
         xmax = gr->GetXaxis()->GetXmax();
         ymin = gr->GetYaxis()->GetXmin();
         ymax = gr->GetYaxis()->GetXmax();
      }
      if (ymin > gr->GetYaxis()->GetXmin()) ymin = gr->GetYaxis()->GetXmin();
      if (ymax < gr->GetYaxis()->GetXmax()) ymax = gr->GetYaxis()->GetXmax();
   }

   if (gROOT->GetListOfCanvases()->GetEntries() == 0) new TCanvas;
   std::stringstream ss_title;
   ss_title << "event " << evtNo << " ch ";
   for (unsigned ichannel=0; ichannel<channels.size(); ++ichannel) ss_title << channels[ichannel] << " ";
   gPad->DrawFrame(xmin,ymin, xmax,ymax, ss_title.str().c_str());
   for (unsigned ichannel=0; ichannel<graphs.size(); ++ichannel) graphs[ichannel]->Draw("p");

   gPad->Modified();
   gPad->Update();
}

void evt4(Int_t evtNo=0, TCut cut="", TCanvas* can=0)
{
   // if (ch1==-1 && ch2==-1 && ch3==-1 && ch4==-1) {
   //    cout<< "Draws:    y vs x for tree with name \"t\"" <<endl;
   //    cout<< "Usage:    evt(TCut cut=\"\", Int_t evtNo=0, ch1, ch2, ch3, ch4)" <<endl;
   //    cout<< "Example:  evt(\"x>0&&x<20\", 9, 1,2,3)" <<endl;
   //    return;
   // }

   TTree *t = (TTree*) gDirectory->Get("t");
   if (!t) {
      cout<< "No tree \"t\" was found in the current directory " << gDirectory->GetName() <<endl;
      return;
   }

   if (!can) {
      TCanvas* can_old = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("can4");
      if (can_old) {
         can_old->cd();
         can_old->Clear();
         can = can_old;
      }
      else can = new TCanvas("can4","can4",800,800);     // width, height
   }
   can->Divide(2,2);

   Int_t marker0 = t->GetMarkerStyle();
   t->SetMarkerStyle(7);
   Int_t color0 = t->GetMarkerColor();
   t->SetMarkerColor(2);

   for (int i=0; i<4; ++i) {
      Int_t ch = i+1;
      can->cd(ch);
      TCut cutevt = TCut(Form("Entry$==%d&&ch==%d",evtNo,ch)) + cut;
      t->Draw("-y:x", cutevt, "");
      t->GetHistogram()->SetTitle(Form("evt_%d_ch_%d",evtNo,ch));
      gPad->Update();
      gPad->Modified();
   }
   can->cd();
   t->SetMarkerStyle(marker0);
   t->SetMarkerColor(color0);

   gPad->Modified();
   gPad->Update();
}

void evt8(Int_t evtNo=0, TCut cut="", TCanvas* can=0)
{
   // if (ch1==-1 && ch2==-1 && ch3==-1 && ch4==-1) {
   //    cout<< "Draws:    y vs x for tree with name \"t\"" <<endl;
   //    cout<< "Usage:    evt(TCut cut=\"\", Int_t evtNo=0, ch1, ch2, ch3, ch4)" <<endl;
   //    cout<< "Example:  evt(\"x>0&&x<20\", 9, 1,2,3)" <<endl;
   //    return;
   // }

   TTree *t = (TTree*) gDirectory->Get("t");
   if (!t) {
      cout<< "No tree \"t\" was found in the current directory " << gDirectory->GetName() <<endl;
      return;
   }

   if (!can) {
      TCanvas* can_old = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("can8");
      if (can_old) {
         can_old->cd();
         can_old->Clear();
         can = can_old;
      }
      else can = new TCanvas("can8","can8",800,800);     // width, height
   }
   can->Divide(2,2);

   Int_t marker0 = t->GetMarkerStyle();
   t->SetMarkerStyle(7);
   Int_t color0 = t->GetMarkerColor();
   t->SetMarkerColor(4);

   for (int i=0; i<4; ++i) {
      Int_t ch = i+5;
      can->cd(i+1);
      TCut cutevt = TCut(Form("Entry$==%d&&ch==%d",evtNo,ch)) + cut;
      t->Draw("-y:x", cutevt, "");
      t->GetHistogram()->SetTitle(Form("evt_%d_ch_%d",evtNo,ch));
      //gPad->Update();
      //gPad->Modified();
   }
   can->cd();
   t->SetMarkerStyle(marker0);
   t->SetMarkerColor(color0);

   gPad->Modified();
   gPad->Update();
}

void eloop(TCut cut="", Int_t evtNo=0, Int_t ch1=-1, Int_t ch2=-1, Int_t ch3=-1, Int_t ch4=-1, Int_t ch5=-1, Int_t ch6=-1, Int_t ch7=-1, Int_t ch8=-1)
{
   if (ch1==-1 || ch2==-1) {
      cout<< "Draws:    (in loop) y vs x for tree with name \"t\"" <<endl;
      cout<< "Usage:    evtloop(TCut cut=\"\", Int_t evtNo=0, ch1,ch2,...)" <<endl;
      cout<< "Example:  evtloop(\"x>0&&x<20\", 9, 1,2,3)" <<endl;
      return;
   }

   TTimer timer("gSystem->ProcessEvents();",50,kFALSE);  // process mouse events every 50 ms

   Int_t ievtNo = evtNo;
   while (ievtNo >= 0)
   {
      evt(cut, ievtNo, ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8);

      timer.Start();    // start processing of mouse events on TCanvas

      ++ievtNo;
      cout<< "entry (<CR>: next entry " << ievtNo << ", N: new TCanvas, S: save event  Q: stop, !command: execute command) ";
      std::string line;
      std::getline(cin, line);

      timer.Stop();     // disable timer after <CR>

      if (line.size() == 0) continue;                                                           // <CR>
      if (line.find("!") != std::string::npos) {                                                // !command
         line.erase(line.find("!"), 1);
         gROOT->ProcessLine(line.c_str());
         --ievtNo;
         continue;
      }
      if (line.find("N") != std::string::npos || line.find("n") != std::string::npos) {         // N or n
         new TCanvas;
         --ievtNo;
         continue;
      }
      if (line.find("S") != std::string::npos || line.find("s") != std::string::npos) {         // N or n
         --ievtNo;
         TTree* t = (TTree*) gDirectory->Get("t");
         gPad->SaveAs(Form("%s_event_%d_ch_%d_%d.png", t->GetCurrentFile()->GetName(),ievtNo,ch1,ch2));
         continue;
      }
      if (line.find("Q") != std::string::npos || line.find("q") != std::string::npos) break;    // Q or q
      std::stringstream ss(line);                                                               // number
      Int_t number = -1;
      if (ss >> number) ievtNo = number;
      else break;
   }
}

void eloop4(TCut cut="", Int_t evtNo=0)
{
   // if (ch1==-1 || ch2==-1) {
   //    cout<< "Draws:    (in loop) y vs x for tree with name \"t\"" <<endl;
   //    cout<< "Usage:    evtloop(TCut cut=\"\", Int_t evtNo=0, ch1,ch2,...)" <<endl;
   //    cout<< "Example:  evtloop(\"x>0&&x<20\", 9, 1,2,3)" <<endl;
   //    return;
   // }

   TTimer timer("gSystem->ProcessEvents();",50,kFALSE);  // process mouse events every 50 ms

   TCanvas* can = new TCanvas("can4", "can4", 800,800);

   Int_t ievtNo = evtNo;
   while (ievtNo >= 0)
   {
      evt4(ievtNo, cut, can);

      timer.Start();    // start processing of mouse events on TCanvas

      ++ievtNo;
      cout<< "entry (<CR>: next entry " << ievtNo << ", N: new TCanvas, S: save event  Q: stop, !command: execute command) ";
      std::string line;
      std::getline(cin, line);

      timer.Stop();     // disable timer after <CR>

      if (line.size() == 0) continue;                                                           // <CR>
      if (line.find("!") != std::string::npos) {                                                // !command
         line.erase(line.find("!"), 1);
         gROOT->ProcessLine(line.c_str());
         --ievtNo;
         continue;
      }
      if (line.find("N") != std::string::npos || line.find("n") != std::string::npos) {         // N or n
         new TCanvas;
         --ievtNo;
         continue;
      }
      // if (line.find("S") != std::string::npos || line.find("s") != std::string::npos) {         // N or n
      //    --ievtNo;
      //    TTree* t = (TTree*) gDirectory->Get("t");
      //    gPad->SaveAs(Form("%s_event_%d_ch_%d_%d.png", t->GetCurrentFile()->GetName(),ievtNo,ch1,ch2));
      //    continue;
      // }
      if (line.find("Q") != std::string::npos || line.find("q") != std::string::npos) break;    // Q or q
      std::stringstream ss(line);                                                               // number
      Int_t number = -1;
      if (ss >> number) ievtNo = number;
      else break;
   }
}
