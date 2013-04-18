#include <TTree.h>
//
// namespace Tree
//
namespace Tree
{
   Float_t  met_pt;              // missingET
   Float_t  met_phi;
   Float_t  mht_pt;              // missingHT
   Float_t  mht_phi;
   Float_t  A;                   // Asymmetry = (met - mht)/(met + mht)

   Float_t  zmumu_m;
   Float_t  zee_m;

   Float_t  wmunu_mt;
   Float_t  wenu_mt;

   Float_t  dphi_jet01;          // azimuthal angle between two leading jets, cut < 2.88
   Float_t  m_jet01;             // invariant mass of two leading jets

   Int_t    ntag;                // the number of b-tagged jets

   Int_t    njet;
   Float_t  jet_pt[100];
   Float_t  jet_phi[100];
   Float_t  jet_eta[100];
   Float_t  jet_dphi_met[100];
   Float_t  jet_mt[100];
   Int_t    nmu;
   Float_t  mu_pt[100];
   Float_t  mu_phi[100];
   Float_t  mu_eta[100];
   Float_t  mu_dphi_met[100];
   Float_t  mu_dR[100];
   Float_t  mu_dRpt[100];
   Int_t    mu_nseg[100];
   Float_t  mu_mt[100];
   Float_t  mu_isocal[100];      // calorimeter isolation: difference of calorimeter energy in cones 0.6 and 0.4 around the muon, for isolated muons < 3.5 GeV
   Float_t  mu_etHalo[100];      // total transverse energy measured in the calorimeter in the anulus between 0.1 and 0.4 centered around the muon, for isolated muons cut < 2.5 GeV
   Float_t  mu_ettrack[100];     // track isolation: sum of tracks pT in cone 0.5, for isolated muon < 2.5 GeV
   Int_t    nem;
   Float_t  em_pt[100];
   Float_t  em_phi[100];
   Float_t  em_eta[100];
   Float_t  em_dphi_met[100];
   Float_t  em_dR[100];
   Float_t  em_mt[100];
   
   void clear()
   {
      met_pt         = 0;
      met_phi        = 0;
      mht_pt         = 0;
      mht_phi        = 0;
      A              = 0;
      zmumu_m     = 0;
      zee_m     = 0;
      wmunu_mt    = 0;
      wenu_mt    = 0;
      dphi_jet01  = 0;
      m_jet01  = 0;

      ntag        = 0;

      njet        = 0;
      nmu        = 0;
      nem        = 0;
   }
   
   void book(TTree* tree) {
      tree->Branch("met_pt",          &met_pt,         "met_pt/F");
      tree->Branch("met_phi",          &met_phi,         "met_phi/F");
      tree->Branch("mht_pt",          &mht_pt,         "mht_pt/F");
      tree->Branch("mht_phi",          &mht_phi,         "mht_phi/F");
      tree->Branch("A",          &A,         "A/F");
      tree->Branch("zmumu_m",          &zmumu_m,         "zmumu_m/F");
      tree->Branch("zee_m",          &zee_m,         "zee_m/F");
      tree->Branch("wmunu_mt",          &wmunu_mt,         "wmunu_mt/F");
      tree->Branch("wenu_mt",          &wenu_mt,         "wenu_mt/F");
      tree->Branch("dphi_jet01",          &dphi_jet01,         "dphi_jet01/F");
      tree->Branch("m_jet01",          &m_jet01,         "m_jet01/F");
      tree->Branch("ntag",          &ntag,         "ntag/I");
      tree->Branch("njet",          &njet,         "njet/I");
      tree->Branch("jet_pt",          &jet_pt,         "jet_pt[njet]/F");
      tree->Branch("jet_phi",          &jet_phi,         "jet_phi[njet]/F");
      tree->Branch("jet_eta",          &jet_eta,         "jet_eta[njet]/F");
      tree->Branch("jet_dphi_met",          &jet_dphi_met,         "jet_dphi_met[njet]/F");
      tree->Branch("jet_mt",          &jet_mt,         "jet_mt[njet]/F");
      tree->Branch("nmu",          &nmu,         "nmu/I");
      tree->Branch("mu_pt",          &mu_pt,         "mu_pt[nmu]/F");
      tree->Branch("mu_phi",          &mu_phi,         "mu_phi[nmu]/F");
      tree->Branch("mu_eta",          &mu_eta,         "mu_eta[nmu]/F");
      tree->Branch("mu_dphi_met",          &mu_dphi_met,         "mu_dphi_met[nmu]/F");
      tree->Branch("mu_dR",          &mu_dR,         "mu_dR[nmu]/F");
      tree->Branch("mu_dRpt",          &mu_dRpt,         "mu_dRpt[nmu]/F");
      tree->Branch("mu_nseg",          &mu_nseg,         "mu_nseg[nmu]/I");
      tree->Branch("mu_mt",          &mu_mt,         "mu_mt[nmu]/F");
      tree->Branch("mu_isocal",          &mu_isocal,         "mu_isocal[nmu]/F");
      tree->Branch("mu_etHalo",          &mu_etHalo,         "mu_etHalo[nmu]/F");
      tree->Branch("mu_ettrack",          &mu_ettrack,         "mu_ettrack[nmu]/F");
      tree->Branch("nem",          &nem,         "nem/I");
      tree->Branch("em_pt",          &em_pt,         "em_pt[nem]/F");
      tree->Branch("em_phi",          &em_phi,         "em_phi[nem]/F");
      tree->Branch("em_eta",          &em_eta,         "em_eta[nem]/F");
      tree->Branch("em_dphi_met",          &em_dphi_met,         "em_dphi_met[nem]/F");
      tree->Branch("em_dR",          &em_dR,         "em_dR[nem]/F");
      tree->Branch("em_mt",          &em_mt,         "em_mt[nem]/F");
   }
   void connect(TTree* tree)                                // need for event-by-event analysis
   {   
      // connects tree buffers with variables to use for event-by-event analysis
      tree->SetBranchAddress("met_pt",          &met_pt);
      tree->SetBranchAddress("met_phi",          &met_phi);
      tree->SetBranchAddress("mht_pt",          &mht_pt);
      tree->SetBranchAddress("mht_phi",          &mht_phi);
      tree->SetBranchAddress("A",          &A);
      tree->SetBranchAddress("zmumu_m",          &zmumu_m);
      tree->SetBranchAddress("zee_m",          &zee_m);
      tree->SetBranchAddress("wmunu_mt",          &wmunu_mt);
      tree->SetBranchAddress("wenu_mt",          &wenu_mt);
      tree->SetBranchAddress("dphi_jet01",          &dphi_jet01);
      tree->SetBranchAddress("m_jet01",          &m_jet01);
      tree->SetBranchAddress("ntag",          &ntag);
      tree->SetBranchAddress("njet",          &njet);
      tree->SetBranchAddress("jet_pt",          &jet_pt);
      tree->SetBranchAddress("jet_phi",          &jet_phi);
      tree->SetBranchAddress("jet_eta",          &jet_eta);
      tree->SetBranchAddress("jet_dphi_met",          &jet_dphi_met);
      tree->SetBranchAddress("jet_mt",          &jet_mt);
      tree->SetBranchAddress("nmu",          &nmu);
      tree->SetBranchAddress("mu_pt",          &mu_pt);
      tree->SetBranchAddress("mu_phi",          &mu_phi);
      tree->SetBranchAddress("mu_eta",          &mu_eta);
      tree->SetBranchAddress("mu_dphi_met",          &mu_dphi_met);
      tree->SetBranchAddress("mu_dR",          &mu_dR);
      tree->SetBranchAddress("mu_dRpt",          &mu_dRpt);
      tree->SetBranchAddress("mu_nseg",          &mu_nseg);
      tree->SetBranchAddress("mu_mt",          &mu_mt);
      tree->SetBranchAddress("mu_isocal",          &mu_isocal);
      tree->SetBranchAddress("mu_etHalo",          &mu_etHalo);
      tree->SetBranchAddress("mu_ettrack",          &mu_ettrack);
      tree->SetBranchAddress("nem",          &nem);
      tree->SetBranchAddress("em_pt",          &em_pt);
      tree->SetBranchAddress("em_phi",          &em_phi);
      tree->SetBranchAddress("em_eta",          &em_eta);
      tree->SetBranchAddress("em_dphi_met",          &em_dphi_met);
      tree->SetBranchAddress("em_dR",          &em_dR);
      tree->SetBranchAddress("em_mt",          &em_mt);
   }
}  // namespace Tree

#ifndef MkProPlain_C
#define MkProPlain_C

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "tmbpro/MkPro.h"
#else
#include "tmbpro/tmbpro/MkPro.h"
#endif

#include <TH2.h>
#include <TCanvas.h>

#include <TROOT.h>         // for ct()
#include <TIterator.h>     // for ct()

#include <TText.h>

#include <iostream>
#include <sstream>         // for updating histogram title

using std::cout;     using std::endl;

////////////////////////////////////
//
//   global pointer to current tree
//
TTree* t;
////////////////////////////////////

/*
root -l tmbpro/macros/Make_so.C
.L MkProPlain.C+

MkPro t_data("wbase/bt-np.root")
t_data.Select("data-out.root")

MkPro t_zh("wbase/bt-zh125nunubb-wbase_mkTree.root")
t_zh.Select("zh125nunubb-out.root")

MkPro t_znunujj("wbase/bt-znunujj-wbase_mkTree.root")
t_znunujj.Select("znunujj-out.root")

MkPro t_wtaunujj("wbase/bt-wtaunujj-wbase_mkTree.root")
t_wtaunujj.Select("wtaunujj-out.root")

MkPro t_zz("wbase/bt-zznunubb-wbase_mkTree.root")
t_zz.Select("zznunubb-out.root")

MkPro t_zmumubb("wbase/bt-zmumubb-wbase_mkTree.root")
t_zmumubb.Select("zmumubb-out.root")

MkPro t_zmumujj("wbase/bt-zmumujj-wbase_mkTree.root")
t_zmumujj.Select("zmumujj-out.root")

MkPro t_ttbbl2j("wbase/bt-ttbbl2j-wbase_mkTree.root")
t_ttbbl2j.Select("ttbbl2j-out.root")

MkPro t_wh("wbase/bt-wh125munubb-wbase_mkTree.root")
t_wh.Select("wh125munubb-out.root")

MkPro t_wzmunubb("wbase/bt-wzmunubb-wbase_mkTree.root")
t_wzmunubb.Select("wzmunubb-out.root")

.q
*/

// need to build shared libraries
void  TMBRoot::Loop(Int_t num, Int_t first) {cout<< "num=" << num << " first=" << first <<endl;}
Int_t TMBRoot::Select(const char* outfile, Int_t num, Int_t first) {cout<< "outfile=" << outfile << " num=" << num << " first=" << first <<endl; return 0;}
Int_t MkPro::Select(const char* outfile, Int_t num, Int_t first)
{
   cout<< "outfile=" << outfile << " num=" << num << " first=" << first <<endl;

   TFile* ofile = TFile::Open(outfile, "recreate");
   TTree* tree = new TTree("t", "subtree");
   Tree::book(tree);

   t = tree;                                          // assign the global pointer to current tree

   Int_t nentries = Int_t(_mktree->GetEntries());
   cout<< "Loop: nentries = " << nentries <<endl;
   if (num == 0) num = nentries;
   
   Int_t jentry=first;
   while ((jentry<num) && GetEntry(jentry++))
   {
      if (jentry == 10)       cout << "Loop: jentry: " << jentry << endl;
      if (jentry%1000 == 0)   cout << "Loop: jentry: " << jentry << endl;

      // insert code here

      // bool zmumu  = _mk->zmu.Mu(0) && (_mk->zmu.Mu(0)->pT() > 15. || _mk->zmu.Mu(1)->pT() > 15.);
      // Z bold
      bool zmumu_Z_iso = true
         && _mk->zmu.Mu(0)
         && _mk->zmu.Mu(0)->pT() > 15.
         && _mk->zmu.Mu(1)->pT() > 15.
         && ((_mk->zmu.Mu(0)->isocal < 4.5 && _mk->zmu.Mu(0)->ettrack < 4.0) || (_mk->zmu.Mu(1)->isocal < 4.5 && _mk->zmu.Mu(1)->ettrack < 4.0))
      ;

      // bool wmu    = _mk->wmu.Mu()  && _mk->met.pt > 15. && _mk->wmu.Mu()->pT()  > 15.;
      // W bold
      bool wmu_W_iso = true
         && _mk->wmu.Mu()
         && _mk->met.pt            > 20.
         && _mk->wmu.Mu()->pT()    > 20.
         && _mk->wmu.Mu()->isocal  < 4.5
         && _mk->wmu.Mu()->ettrack < 4.0
      ;

      // bool zee    = _mk->ze.EM(0)  && (_mk->ze.EM(0)->pT()  > 15. || _mk->ze.EM(1)->pT()  > 15.);
      // Zee bold
      bool Zee    = _mk->ze.EM(0)  && (_mk->ze.EM(0)->pT()  > 15. || _mk->ze.EM(1)->pT()  > 15.) && _mk->ze.m > 60. && _mk->ze.m < 200.;

      // We bold
      bool We     = _mk->we.EM()   && _mk->met.pt > 25. && _mk->we.EM()->pT()   > 25. && _mk->we.EM()->prob > .01 && _mk->we.EM()->dphi_met > .4; // tight quality but no cleaning

      // regular selection
      bool base = true
         // TCut clean     ("clean",      "cpf.clean15>0 && badu.ptmax<15 && muok");
         // TCut tmin      ("tmin",       "jet1pt>15&&met.pt>40&&mht.jet>40&&jet[0].pt>40");
         // TCut noiso     ("noiso",      "t.iemmax<5&&t.imumax<5&&t.imulmax<10");
         // TCut base      ("base",       clean+tmin20+noiso);

         && _mk->mumax < 200
         && _mk->cpf.clean15 > 0
         && _mk->badu.ptmax < 15
         && _mk->njet > 0
         && _mk->Jet(0)->pt > 40
         && _mk->met.pt > 40
      ;

      if (wmu_W_iso || zmumu_Z_iso || We || Zee || base)
      {
         Tree::clear();

         Tree::met_pt = _mk->met.pt;
         Tree::met_phi = _mk->met.phi;
         Tree::mht_pt = _mk->mht.pt;
         Tree::mht_phi = _mk->mht.phi;
         Tree::A = (Tree::met_pt - Tree::mht_pt) / (Tree::met_pt + Tree::mht_pt);

         Tree::zmumu_m = _mk->zmu.m;
         Tree::zee_m = _mk->ze.m;

         Tree::wmunu_mt = _mk->wmu.mt;
         Tree::wenu_mt = _mk->we.mt;

         Tree::dphi_jet01 = _mk->dphi.j01;
         Tree::m_jet01 = _mk->Mjet;
         Tree::ntag = _mk->jlip.npos;

         Tree::njet = _mk->Njet();
         for (int ijet=0; ijet<_mk->Njet(); ++ijet) {
            Tree::jet_pt[ijet] = _mk->Jet(ijet)->pt;
            Tree::jet_phi[ijet] = _mk->Jet(ijet)->phi;
            Tree::jet_eta[ijet] = _mk->Jet(ijet)->eta;
            Tree::jet_mt[ijet] = _mk->Jet(ijet)->mt;
         }

         Tree::nmu = _mk->Nmu();
         for (int imu=0; imu<_mk->Nmu(); ++imu) {
            Tree::mu_pt[imu] = _mk->Mu(imu)->pt;
            Tree::mu_phi[imu] = _mk->Mu(imu)->phi;
            Tree::mu_eta[imu] = _mk->Mu(imu)->eta;
            Tree::mu_dphi_met[imu] = _mk->Mu(imu)->dphi_met;
            Tree::mu_dR[imu] = _mk->Mu(imu)->dR;
            Tree::mu_dRpt[imu] = _mk->Mu(imu)->dR*_mk->Mu(imu)->pt;
            Tree::mu_nseg[imu] = _mk->Mu(imu)->nseg;
            Tree::mu_mt[imu] = _mk->Mu(imu)->mt;
            Tree::mu_isocal[imu] = _mk->Mu(imu)->isocal;
            Tree::mu_etHalo[imu] = _mk->Mu(imu)->etHalo;
            Tree::mu_ettrack[imu] = _mk->Mu(imu)->ettrack;
         }

         Tree::nem = _mk->Nem();
         for (int iem=0; iem<_mk->Nem(); ++iem) {
            Tree::em_pt[iem] = _mk->EM(iem)->pt;
            Tree::em_phi[iem] = _mk->EM(iem)->phi;
            Tree::em_eta[iem] = _mk->EM(iem)->eta;
            Tree::em_dphi_met[iem] = _mk->EM(iem)->dphi_met;
            Tree::em_dR[iem] = _mk->EM(iem)->dR;
            Tree::em_mt[iem] = _mk->EM(iem)->mt;
         }
         
         tree->Fill();
      }
   }

   ofile->Write();
   return tree->GetEntries();
}

void MkPro::Loop(Int_t num, Int_t first)
{
   TH1F* h_mu1_pt = new TH1F("h_mu1_pt", "h_mu1_pt", 100, 0, 100);
   TH1F* h_mu1_pt1 = new TH1F("h_mu1_pt1", "h_mu1_pt1", 100, 0, 100);
   TH1F* h_mu1_pt2 = new TH1F("h_mu1_pt2", "h_mu1_pt2", 100, 0, 100);

   _mktree->Project(h_mu1_pt2->GetName(), "sep.mu1.pt", "sep.mu1_pt>0");

   Int_t nentries = Int_t(_mktree->GetEntries());
   cout<< "Loop: nentries = " << nentries <<endl;
   if (num == 0) num = nentries;
   
   Int_t jentry=first;
   while ((jentry<num) && GetEntry(jentry++))
   {
      if (jentry == 10)       cout << "Loop: jentry: " << jentry << endl;
      if (jentry%1000 == 0)   cout << "Loop: jentry: " << jentry << endl;

      // insert code here

      if (true
            //&& _mk->Njet() == 2
            //&& _mk->Jet(0)->phi == _mk->Jet(1)->phi
         )
      {
         //if (_mk->sep.nmu > 0) {
         if (_mk->sep.mu1_pt > 0) {
            // cout<< "_mk->sep.mu1_pt = " << _mk->sep.mu1_pt <<endl;
            h_mu1_pt->Fill(_mk->sep.mu1_pt);
            h_mu1_pt1->Fill(_mk->sep.get_mu1()->pt);
         }
      }
   }

   new TCanvas;
   h_mu1_pt->Draw();
   new TCanvas;
   h_mu1_pt1->Draw();
   new TCanvas;
   h_mu1_pt2->Draw();
}
#endif
