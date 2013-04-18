//
// File: analysis.C (Selcuk & Zongru, June 8, 2010)
// Modified by Andriy Zatserklyaniy, May 9, 2011
//
// (1) Change the inputs
//
//         fileName = Raw_Data_July1_Ped_300V
//         nAPVs    = 4
//         nEvents  = -1
//
// (2) Run this script by the following command
//
//         root analysis.C
//
// (3) Output
//
//         Raw_Data_July1_Ped_300V.eps
//         Raw_Data_July1_Ped_300V.root
//
// (4) To play with the ntuple in analysis.root with script, for example,
//
//         root
//         root [0] f = new TFile("Raw_Data_July1_Ped_300V.root")
//         root [1] f->ls()
//         root [2] h2d->Draw()
//         root [3] tree->Print()
//         root [4] tree->Draw("adc")
//         root [5] tree->Draw("adc-90","apv==1&&adc>90")
//         root [6] tree->Draw("adc:channel","apv==1&&event<=10")
//         root [7] .q
//

#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH2.h>
#include <TGraph.h>
#include <TProfile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>

using std::cout;        using std::endl;

// stuff from ~/macros/utils.C
std::ostream& exitl(std::ostream&);
TH1* htemp(const char* name=0, const char* title=0);
void pic(TString pathname="", TString suffix="");
void png(TString pathname="", TString suffix="");
void fitgl(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void fitgr(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void pargaus(Double_t& a, Double_t& mean, Double_t& sigma, const char* hname=0);

void striptree(TString fileName="", Int_t nAPVs=4)
{
  //===========================================================================
  // Change the inputs
  //===========================================================================

  //TString fileName = "Raw_Data_July1_Ped_300V";
  //TString fileName = "Raw_Data_July1_Source_300V";
  //TString fileName = "Raw_Data_July1_Ped_200V";
  //TString fileName = "Raw_Data_Pedestal_300V";
  // TString fileName = "Raw_Data_Source_300V";

  // TString fileName = "Raw_Data_Sept_3_1";
  // TString fileName = "Raw_Data_test_3_7";
  // TString fileName = "Raw_Data_Albox_2011_03_31_1";
  //-- TString fileName = "Raw_Data_04_01_AL_center_1";
  // TString fileName = "Raw_Data_04_01_AL_7820_1";
  //-- TString fileName = "Raw_Data_7817_04_01_1";
  // TString fileName = "Raw_Data_7817_04_01_bkg_1";

  //int nAPVs = 4;

  //TString fileName = "Raw_Data_July1_Source_200V";
  //int nAPVs = 1;

  int nEvents = -1; // -1 means all

  TString dat = fileName;
  TString eps = fileName + "-events.eps";
  TString root = fileName + "-events.root";

  TFile* file = new TFile(root, "RECREATE");

  //===========================================================================
  // Read the data and make a 2D scatter plot
  //===========================================================================

  ifstream in(dat);
  if (!in) {
     cout<< "File " << dat << " not found" <<endl;
     return;
  }

  TString comments = "vertical --> APV0, APV1, ... of event 1, ...";
  TString s;
  while (1) {
    s.ReadLine(in);
    if (s.BeginsWith(comments))
      break;
  }
  gStyle->SetOptStat(10);
  TString title = fileName+"   ADC Counts vs Channel";
  TString title32 = fileName+"   ADC Counts vs Channel for 32 ch";
  TH2D *h2d = new TH2D("h2d", title.Data(), nAPVs*128, 0, nAPVs*128, 200, 0, 200);
  TH2D *h2d32 = new TH2D("h2d32", title.Data(), 16, 0, 16, 200, 0, 200);
  TProfile *profile = new TProfile("profile", title32.Data(), nAPVs*128, 0, nAPVs*128);

  int event, apv, channel, adc;
  // TTree* tree = new TTree("tree", "tree");
  // tree->Branch("event",   &event,   "event/I");
  // tree->Branch("apv",     &apv,     "apv/I");
  // tree->Branch("channel", &channel, "channel/I");
  // tree->Branch("adc",     &adc,     "adc/I");

  Int_t raw[512];
  TTree* etree = new TTree("etree", "etree");
  etree->Branch("raw",   &raw,   "raw[512]/I");    // etree->Draw("raw:Iteration$","Entry$==0")
  etree->SetMarkerStyle(6);
  etree->SetMarkerColor(2);

  int n = 0;
  while (1) {
    in >> s;
    if (!in.good())
      break;
    if (s.BeginsWith("]DATA]")) // last line
      break;
    n++;
    if ((nEvents != -1) && (n > nEvents))
      break;
    event = atoi(s.Data());
    if (event%100 == 0)
      cout << "Processing event " << event << endl;
    for (int i=0; i<2; i++)
      in >> s; // time

    Int_t ichannel = 0;
    for (apv=1; apv<=nAPVs; apv++) {
      in >> s; // header + analog APV data + one tick mark
      int sequenceNumber = 0;
      TString subString = "";
      int length = s.Length();
      for (int j=0; j<length; j++) {
        char c = s[j];
        if (c != ',') {
          subString += c;
        } else {
          if (sequenceNumber >= 12) {
            channel = sequenceNumber - 12;
            adc = atoi(subString.Data());
            h2d->Fill((apv - 1) * 128 + channel + 0.5, adc);
            h2d32->Fill(((apv - 1) * 128 + channel)/32 + 0.5, adc);
            profile->Fill((apv - 1) * 128 + channel + 0.5, adc);
            //tree->Fill();
            raw[ichannel] = adc;
            ++ichannel;
          }
          subString = "";
          sequenceNumber++;
        }
      }
    }
    etree->Fill();
  }
  in.close();

  TCanvas *canvas = new TCanvas;
  canvas->Divide(1,2);
  canvas->cd();
  canvas->cd(1); h2d->Draw();
  canvas->cd(2); profile->Draw();
  canvas->SaveAs(eps);
  // //-- TFile file(root, "RECREATE");
  // h2d->Write();
  // profile->Write();
  // tree->Write();
  // //file.Close();
  cout<< "\nWrite " << etree->GetEntries() << " of tree " << etree->GetName() << " into file " << file->GetName() <<endl;
  cout<< "To plot strip content use e.g. etree->Draw(\"raw:Iteration$\",\"Entry$==0\")" <<endl;
  file->Write();
}

void convert()
{
   striptree("Raw_Data_FZ320P_05_MSSD_2_250V_K237_Pedestal.dat");
   striptree("Raw_Data_FZ320P_05_MSSD_250V_K237_Position_1.dat");
   striptree("Raw_Data_FZ320P_05_MSSD_250V_K237_Position_2.dat");
   striptree("Raw_Data_FZ320P_05_MSSD_2_250V_K237_Position_3.dat");
   striptree("Raw_Data_FZ320P_05_MSSD_2_250V_K237_Position_4.dat");
}

TChain* init_chain()
{
   TChain* chain = new TChain("etree");
   chain->Add("Raw_Data_FZ320P_05_MSSD_250V_K237_Position_1.dat-events.root");
   chain->Add("Raw_Data_FZ320P_05_MSSD_250V_K237_Position_2.dat-events.root");
   chain->Add("Raw_Data_FZ320P_05_MSSD_2_250V_K237_Position_3.dat-events.root");
   chain->Add("Raw_Data_FZ320P_05_MSSD_2_250V_K237_Position_4.dat-events.root");
   return chain;
}

void allstrips()
{
   TChain* chain = new TChain("etree");
   chain->Add("Raw_Data_FZ320P_05_MSSD_250V_K237_Position_1.dat-events.root");
   chain->Add("Raw_Data_FZ320P_05_MSSD_250V_K237_Position_2.dat-events.root");
   chain->Add("Raw_Data_FZ320P_05_MSSD_2_250V_K237_Position_3.dat-events.root");
   chain->Add("Raw_Data_FZ320P_05_MSSD_2_250V_K237_Position_4.dat-events.root");

   new TCanvas;
   chain->Draw("raw:Iteration$","");
   // change FZ320P to FZ320N in the title
   htemp()->SetTitle("Raw Data FZ320N_05_MSSD_2 250V K237 all events;strip No.;ADC");
   //-- png("Raw_Data_FZ320P_05_MSSD_2_250V_K237-allstrips");
}

void process()
{
   Int_t raw[512];   // buffer for input signal and bkg trees
   // buffers for output trees
   Int_t sig[512];
   Int_t cmsig[512];
   Int_t cm[16];

   // pedestal
   const char* fbkg_name = "Raw_Data_FZ320P_05_MSSD_2_250V_K237_Pedestal.dat-events.root";
   TFile* fbkg = TFile::Open(fbkg_name);
   if (!fbkg) cout<< "File not found: " << fbkg <<endl<<exitl;

   TTree* tree = (TTree*) fbkg->Get("etree");
   tree->SetBranchAddress("raw", &raw);

   TH2* h2d = (TH2*) fbkg->Get("h2d");
   new TCanvas;
   h2d->Draw();

   TProfile* profile = (TProfile*) fbkg->Get("profile"); 
   //new TCanvas;
   //profile->Draw();

   Int_t pedestal[512];
   for (int i=0; i<512; i++) {
      // pedestal[i] = profile->GetBinContent(i+1) - 0.5;
      // pedestal[i] = profile->GetBinContent(i+1) + 0.5;
      pedestal[i] = profile->GetBinContent(i+1);
   }

   // cout << "\nPedestals for every channel\n" << endl;
   // for (int i=0; i<512; i++) {
   //    cout << pedestal[i] << " ";
   //    if (i>0 && (i+1)%128==0)
   //       cout << endl;
   // }

   cout<< "processing bkg" <<endl;

   // output file with tree
   const char* obfname = "FZ320P_05_MSSD_2-bkg.root";
   TFile* obfile = TFile::Open(obfname, "recreate");

   TTree* btree = new TTree("btree", "btree");
   btree->Branch("sig", &sig, "sig[512]/I");
   btree->Branch("cmsig", &cmsig, "cmsig[512]/I");
   btree->Branch("cm", &cm, "cm[16]/I");
   btree->SetMarkerStyle(6);
   btree->SetMarkerColor(2);

   for (int jentry=0; jentry<tree->GetEntries(); ++jentry)
   {
      tree->GetEvent(jentry);
      // sig
      for (int i=0; i<512; ++i) {
         sig[i] = raw[i] - pedestal[i];
      }
      // calc common mode
      Int_t group32[32];
      Int_t index32[32];
      for (int igroup=0; igroup<16; ++igroup) {
         for (int istrip=0; istrip<32; ++istrip)   // istrip is number inside group of 32
         {
            group32[istrip] = sig[igroup*32 + istrip];
         }
         // sort array group32 in ascending order
         TMath::Sort(32, group32, index32, kFALSE);
         Int_t median = group32[index32[14]];
         cm[igroup] = median;
      }
      // subtract common mode
      for (int istrip=0; istrip<512; ++istrip) {
         Int_t igroup = istrip/32;
         cmsig[istrip] = sig[istrip] - cm[igroup];
      }
      // Fill sig, cmsig, cm
      btree->Fill();
   }
   obfile->Write();

   /////////////////////////////////////////////////
   //
   // signal tree
   //
   /////////////////////////////////////////////////

   cout<< "processing signal" <<endl;
   
   TChain* chain = new TChain("etree");
   chain->Add("Raw_Data_FZ320P_05_MSSD_250V_K237_Position_1.dat-events.root");
   chain->Add("Raw_Data_FZ320P_05_MSSD_250V_K237_Position_2.dat-events.root");
   chain->Add("Raw_Data_FZ320P_05_MSSD_2_250V_K237_Position_3.dat-events.root");
   chain->Add("Raw_Data_FZ320P_05_MSSD_2_250V_K237_Position_4.dat-events.root");
   chain->SetBranchAddress("raw", &raw);

   // output file with tree
   const char* osfname = "FZ320P_05_MSSD_2-signal.root";
   TFile* osfile = TFile::Open(osfname, "recreate");

   TTree* stree = new TTree("stree", "stree");
   stree->Branch("sig", &sig, "sig[512]/I");
   stree->Branch("cmsig", &cmsig, "cmsig[512]/I");
   stree->Branch("cm", &cm, "cm[16]/I");
   stree->SetMarkerStyle(6);
   stree->SetMarkerColor(2);

   for (int jentry=0; jentry<chain->GetEntries(); ++jentry)
   {
      chain->GetEvent(jentry);
      // sig
      for (int i=0; i<512; ++i) {
         sig[i] = raw[i] - pedestal[i];
      }
      // calc common mode
      Int_t group32[32];
      Int_t index32[32];
      for (int igroup=0; igroup<16; ++igroup) {
         for (int istrip=0; istrip<32; ++istrip)   // istrip is number inside group of 32
         {
            group32[istrip] = sig[igroup*32 + istrip];
         }
         // sort array group32 in ascending order
         TMath::Sort(32, group32, index32, kFALSE);
         Int_t median = group32[index32[14]];
         cm[igroup] = median;
      }
      // subtract common mode
      for (int istrip=0; istrip<512; ++istrip) {
         Int_t igroup = istrip/32;
         cmsig[istrip] = sig[istrip] - cm[igroup];
      }
      // Fill sig, cmsig, cm
      stree->Fill();
   }
   osfile->Write();

   ////////////////////////////////////////////////////////
   //
   // process trees
   //
   ///////////////////////////////////////////////////////

   cout<< "results" <<endl;

   Double_t a, mean, sigma;

   TH1F* h_sigma_bkg = new TH1F("h_sigma_bkg","CM subtr. noise for groups", 16,0,16);
   TH1F* h_mean_sig = new TH1F("h_mean_sig","CM subtr. signal for groups", 16,0,16);
   TH1F* h_SN = new TH1F("h_SN","Signal to Noise Ratio for groups", 16,0,16);

   new TCanvas;
   btree->Draw("cmsig","Iteration$>=0&&Iteration$<32");
   fitgr(0,0, "Q", "goff", btree->GetHistogram());
   pargaus(a,mean,sigma,"htemp");
   //cout<< "mean = " << mean << " sigma = " << sigma <<endl;
   //-- png("FZ320P_05_MSSD_2-bkg-ex");

   new TCanvas;
   // for (int igroup=0; igroup<16; ++igroup) {
   for (int igroup=0; igroup<15; ++igroup) {
      Int_t ch1 = igroup*32;
      Int_t ch2 = (igroup+1)*32;
      btree->Draw("cmsig",Form("Iteration$>=%d&&Iteration$<%d",ch1,ch2),"");
      fitgr(0,0, "", "", btree->GetHistogram());
      pargaus(a,mean,sigma,"htemp");
      // h_sigma_bkg->Fill(igroup, sigma);
      h_sigma_bkg->SetBinContent(igroup+1, sigma);
   }
   new TCanvas;
   h_sigma_bkg->Draw();
   //-- png("FZ320P_05_MSSD_2-bkg-allgroups");

   // signal
   new TCanvas;
   stree->Draw("cmsig","cmsig>8 &&Iteration$>=0&&Iteration$<32");
   //-- png("FZ320P_05_MSSD_2-sig-ex");

   new TCanvas;
   // for (int igroup=0; igroup<16; ++igroup) {
   for (int igroup=0; igroup<15; ++igroup) {
      Int_t ch1 = igroup*32;
      Int_t ch2 = (igroup+1)*32;
      stree->Draw("cmsig",Form("cmsig>8 &&Iteration$>=%d&&Iteration$<%d",ch1,ch2),"");
      gPad->Update();
      gPad->Modified();
      mean = stree->GetHistogram()->GetMean();
      h_mean_sig->SetBinContent(igroup+1, mean);
   }
   new TCanvas;
   h_mean_sig->Draw();
   //-- png("FZ320P_05_MSSD_2-signal-allgroups");

   for (int igroup=0; igroup<16; ++igroup) {
      Double_t signal32 =  h_mean_sig->GetBinContent(igroup+1);
      Double_t noise32 = h_sigma_bkg->GetBinContent(igroup+1);
      Double_t snr = 0;
      if (noise32 > 0) snr = signal32 / noise32;
      h_SN->SetBinContent(igroup+1, snr);
   }
   new TCanvas;
   h_SN->Draw();
   //-- png("FZ320P_05_MSSD_2-SN-allgroups");
}
