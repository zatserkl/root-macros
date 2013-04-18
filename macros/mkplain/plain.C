{
/*
root -l plain.C
*/
   TFile* file_data = TFile::Open("tmb_plain_tree/data-out.root");
   TFile* file_zh = TFile::Open("tmb_plain_tree/zh125nunubb-out.root");
   TFile* file_znunujj = TFile::Open("tmb_plain_tree/znunujj-out.root");
   TFile* file_wtaunujj = TFile::Open("tmb_plain_tree/wtaunujj-out.root");
   TFile* file_zz = TFile::Open("tmb_plain_tree/zznunubb-out.root");
   TFile* file_zmumubb = TFile::Open("tmb_plain_tree/zmumubb-out.root");
   TFile* file_zmumujj = TFile::Open("tmb_plain_tree/zmumujj-out.root");
   TFile* file_ttbbl2j = TFile::Open("tmb_plain_tree/ttbbl2j-out.root");
   TFile* file_wh = TFile::Open("tmb_plain_tree/wh125munubb-out.root");
   TFile* file_wzmunubb = TFile::Open("tmb_plain_tree/wzmunubb-out.root");
   gROOT->cd();

   TTree* data = (TTree*) file_data->Get("t");
   TTree* zh = (TTree*) file_zh->Get("t");
   TTree* znunujj = (TTree*) file_znunujj->Get("t");
   TTree* wtaunujj = (TTree*) file_wtaunujj->Get("t");
   TTree* zz = (TTree*) file_zz->Get("t");
   TTree* zmumubb = (TTree*) file_zmumubb->Get("t");
   TTree* zmumujj = (TTree*) file_zmumujj->Get("t");
   TTree* ttbbl2j = (TTree*) file_ttbbl2j->Get("t");
   TTree* wh = (TTree*) file_wh->Get("t");
   TTree* wzmunubb = (TTree*) file_wzmunubb->Get("t");

   data->SetMarkerStyle(6);
   data->SetMarkerColor(1);
   zh->SetMarkerStyle(6);
   zh->SetMarkerColor(2);
   znunujj->SetMarkerStyle(6);
   znunujj->SetMarkerColor(4);
   wtaunujj->SetMarkerStyle(6);
   wtaunujj->SetMarkerColor(5);
   zz->SetMarkerStyle(6);
   zz->SetMarkerColor(6);
   zmumubb->SetMarkerStyle(6);
   zmumubb->SetMarkerColor(8);
   zmumujj->SetMarkerStyle(6);
   zmumujj->SetMarkerColor(9);
   ttbbl2j->SetMarkerStyle(6);
   ttbbl2j->SetMarkerColor(3);
   wh->SetMarkerStyle(6);
   wh->SetMarkerColor(2);
   wzmunubb->SetMarkerStyle(6);
   wzmunubb->SetMarkerColor(4);

   //-- Cuts --//
/*
data->Draw("mu_dphi_met[0]","mu_nseg[0]==3 &&mu_dRpt[0]<3.5")
zh->Draw("mu_dphi_met[0]","mu_nseg[0]==3 &&mu_dRpt[0]<3.5")

data->Draw("mu_dphi_met[0]:met_pt","mu_nseg[0]==3 &&mu_dRpt[0]<3.5 &&met_pt<300")
zh->Draw("mu_dphi_met[0]:met_pt","mu_nseg[0]==3 &&mu_dRpt[0]<3.5 &&met_pt<300")

// W --> munu
data->Draw("mu_mt","mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200 &&mu_mt>0&&mu_mt<200 &&met_pt>20")
// Z --> mumu
data->Draw("zmumu_m","mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200 &&zmumu_m>0&&zmumu_m<200")

   ROOT's double counting

// t->Draw("zmumu_m","mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200 &&zmumu_m>0&&zmumu_m<200")
<TCanvas::MakeDefCanvas>: created default TCanvas with name c1                               
(Long64_t)1388                                                                               
// t->Draw(">>elist_zmumu_m","mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200 &&zmumu_m>0&&zmumu_m<200")
(Long64_t)1388                                                                                       
// elist_zmumu_m=elist_zmumu_m                                                                       
(class TEventList*)0xa70a750                                                                         
// elist_zmumu_m->GetN()                                                                             
(const Int_t)1102                                                                                    
//                                                                                                   
// t->SetEventList(elist_zmumu_m)
// <TCanvas::MakeDefCanvas>: created default TCanvas with name c1_n2
// t->Draw("zmumu_m","mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200 &&zmumu_m>0&&zmumu_m<200")
(Long64_t)1388                                                                               
// t->Draw("zmumu_m","")
(Long64_t)1102          
//  FCN=378.554 FROM MIGRAD    STATUS=CONVERGED     113 CALLS         114 TOTAL
                     EDM=1.61628e-10    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST           
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE         
   1  Constant     1.57604e+01   6.67771e-01   5.57406e-03  -2.75398e-05       
   2  Mean         4.08643e+01   4.79062e+00   2.33603e-02  -2.14307e-06       
   3  Sigma        5.11155e+01   3.06782e+00   1.12915e-04  -8.75203e-04       
 FCN=23.8001 FROM MIGRAD    STATUS=CONVERGED     164 CALLS         165 TOTAL   
                     EDM=3.89382e-12    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST           
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE         
   1  Constant     3.70882e+01   2.24249e+00   4.00050e-03   2.06889e-07       
   2  Mean         8.94410e+01   6.41268e-01   1.52177e-03   7.03735e-07       
   3  Sigma        1.24307e+01   6.67686e-01   7.24790e-05  -2.34021e-04       
right()                                                                        
//  FCN=43.4245 FROM MIGRAD    STATUS=CONVERGED      68 CALLS          69 TOTAL
                     EDM=5.33869e-09    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST           
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE         
   1  Constant     5.46732e+01   2.81165e+00   6.62445e-03  -2.41021e-05       
   2  Mean         9.00220e+01   4.89370e-01   1.58588e-03  -1.92579e-04       
   3  Sigma        1.18025e+01   5.15554e-01   1.86111e-05  -7.70395e-03       
// right()
//                                   
// htemp()->GetBinWidth(1)
(const Double_t)2.14999999999999991e+00
// 37.09*12.43 / 2.15
(const double)2.14431953488372102e+02
//                                   
// 54.67*11.8 / 2.15                
(const double)3.00049302325581436e+02
// htemp()->GetBinWidth(1)
(const Double_t)2.14999999999999991e+00
// <TCanvas::MakeDefCanvas>: created default TCanvas with name c1_n3
// t->Draw("zmumu_m","mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200 &&zmumu_m>0&&zmumu_m<200")
(Long64_t)1388
// t->Draw("zmumu_m","mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200")
(Long64_t)1388
// t->Draw("zmumu_m","mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20")
(Long64_t)1390
// t->Draw("zmumu_m","")
(Long64_t)1102
//
//
// t->Draw("zmumu_m","mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20&&mu_pt<200")
(Long64_t)1388
// t->Draw("zmumu_m","mu_nseg==3 &&mu_dR>0.5 &&mu_pt>20")
(Long64_t)1390
// <TCanvas::MakeDefCanvas>: created default TCanvas with name c1_n4
// t->Draw("zmumu_m","mu_nseg==3 &&mu_dR>0.5")
(Long64_t)1499
// <TCanvas::MakeDefCanvas>: created default TCanvas with name c1_n5
// t->Draw("zmumu_m","mu_nseg==3")
(Long64_t)1631
// <TCanvas::MakeDefCanvas>: created default TCanvas with name c1_n6
// t->Draw("zmumu_m","")
(Long64_t)1102
// <TCanvas::MakeDefCanvas>: created default TCanvas with name c1_n7
// t->Draw("zmumu_m","(mu_nseg==3)")
(Long64_t)1631


--------------------> the same with mk tree:

bt$ root -l tmbpro/macros/Make_so.C
...
...
// TFile::Open("wbase/bt-np.root")
(class TFile*)0x9b107f0
// mk->Draw("zmu.m", "mu.nseg==3 &&mu.dR>0.5 &&mu.pt>20&&mu.pt<200 &&zmu.m>0&&zmu.m<200")
<TCanvas::MakeDefCanvas>: created default TCanvas with name c1
(Long64_t)1388
// mk->Draw(">>zmu_elist", "mu.nseg==3 &&mu.dR>0.5 &&mu.pt>20&&mu.pt<200 &&zmu.m>0&&zmu.m<200")
(Long64_t)1388
// zmu_elist=zmu_elist
(class TEventList*)0x99d1480
// mk->SetEventList(zmu_elist)
// <TCanvas::MakeDefCanvas>: created default TCanvas with name c1_n2
// zmu_elist->GetN()
(const Int_t)1102
// mk->Draw("zmu.m", "")
(Long64_t)1102
// <TCanvas::MakeDefCanvas>: created default TCanvas with name c1_n3
// mk->Draw("zmu.m", "mu.nseg==3 &&mu.dR>0.5 &&mu.pt>20&&mu.pt<200 &&zmu.m>0&&zmu.m<200")
(Long64_t)1388
// 


   */
}
