void eventsum()
{
   TFile* f = 0;
   TH1* h = 0;
   TLegend* legend = 0;

   gStyle->SetPadGridY(0);
   gStyle->SetPadGridX(0);
   gStyle->SetOptStat(0);

   Double_t lum_sig_900 = 384;
   Double_t lum_sig_2360 = 261;
   Double_t lum_bkg_900 = 6.87;
   Double_t lum_bkg_2360 = 5.46;

   //--
   //-- Opposite Sign
   //--

   // signal

   f = TFile::Open("histos_jPsi_900GeV_STARTUP.root");

   h = (TH1*) f->Get("hReco_mass_OS_glgl_Cut4");
   Double_t n_OS_glgl_sig_900 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_OS_glgl_sig_norm_900 = n_OS_glgl_sig_900/lum_sig_900;

   h = (TH1*) f->Get("hReco_mass_OS_gltr_Cut4");
   Double_t n_OS_gltr_sig_900 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_OS_gltr_sig_norm_900 = n_OS_gltr_sig_900/lum_sig_900;

   h = (TH1*) f->Get("hReco_mass_OS_trtr_Cut4");
   Double_t n_OS_trtr_sig_900 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_OS_trtr_sig_norm_900 = n_OS_trtr_sig_900/lum_sig_900;

   cout<< "n_OS_glgl_sig_900 = " << n_OS_glgl_sig_900 << " n_OS_glgl_sig_norm_900 = " << n_OS_glgl_sig_norm_900 <<endl;
   cout<< "n_OS_gltr_sig_900 = " << n_OS_gltr_sig_900 << " n_OS_gltr_sig_norm_900 = " << n_OS_gltr_sig_norm_900 <<endl;
   cout<< "n_OS_trtr_sig_900 = " << n_OS_trtr_sig_900 << " n_OS_trtr_sig_norm_900 = " << n_OS_trtr_sig_norm_900 <<endl;

   f = TFile::Open("histos_jPsi_2360GeV_STARTUP.root");

   h = (TH1*) f->Get("hReco_mass_OS_glgl_Cut4");
   Double_t n_OS_glgl_sig_2360 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_OS_glgl_sig_norm_2360 = n_OS_glgl_sig_2360/lum_sig_2360;

   h = (TH1*) f->Get("hReco_mass_OS_gltr_Cut4");
   Double_t n_OS_gltr_sig_2360 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_OS_gltr_sig_norm_2360 = n_OS_gltr_sig_2360/lum_sig_2360;

   h = (TH1*) f->Get("hReco_mass_OS_trtr_Cut4");
   Double_t n_OS_trtr_sig_2360 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_OS_trtr_sig_norm_2360 = n_OS_trtr_sig_2360/lum_sig_2360;

   cout<< "n_OS_glgl_sig_2360 = " << n_OS_glgl_sig_2360 << " n_OS_glgl_sig_norm_2360 = " << n_OS_glgl_sig_norm_2360 <<endl;
   cout<< "n_OS_gltr_sig_2360 = " << n_OS_gltr_sig_2360 << " n_OS_gltr_sig_norm_2360 = " << n_OS_gltr_sig_norm_2360 <<endl;
   cout<< "n_OS_trtr_sig_2360 = " << n_OS_trtr_sig_2360 << " n_OS_trtr_sig_norm_2360 = " << n_OS_trtr_sig_norm_2360 <<endl;

   // background

   f = TFile::Open("histos_ppMuXLoose_900GeV_STARTUP.root");

   h = (TH1*) f->Get("hReco_mass_OS_glgl_Cut4");
   Double_t n_OS_glgl_bkg_900 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_OS_glgl_bkg_norm_900 = n_OS_glgl_bkg_900/lum_bkg_900;

   h = (TH1*) f->Get("hReco_mass_OS_gltr_Cut4");
   Double_t n_OS_gltr_bkg_900 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_OS_gltr_bkg_norm_900 = n_OS_gltr_bkg_900/lum_bkg_900;

   h = (TH1*) f->Get("hReco_mass_OS_trtr_Cut4");
   Double_t n_OS_trtr_bkg_900 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_OS_trtr_bkg_norm_900 = n_OS_trtr_bkg_900/lum_bkg_900;

   cout<< "n_OS_glgl_bkg_900 = " << n_OS_glgl_bkg_900 << " n_OS_glgl_bkg_norm_900 = " << n_OS_glgl_bkg_norm_900 <<endl;
   cout<< "n_OS_gltr_bkg_900 = " << n_OS_gltr_bkg_900 << " n_OS_gltr_bkg_norm_900 = " << n_OS_gltr_bkg_norm_900 <<endl;
   cout<< "n_OS_trtr_bkg_900 = " << n_OS_trtr_bkg_900 << " n_OS_trtr_bkg_norm_900 = " << n_OS_trtr_bkg_norm_900 <<endl;

   f = TFile::Open("histos_ppMuXLoose_2360GeV_STARTUP.root");

   h = (TH1*) f->Get("hReco_mass_OS_glgl_Cut4");
   Double_t n_OS_glgl_bkg_2360 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_OS_glgl_bkg_norm_2360 = n_OS_glgl_bkg_2360/lum_bkg_2360;

   h = (TH1*) f->Get("hReco_mass_OS_gltr_Cut4");
   Double_t n_OS_gltr_bkg_2360 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_OS_gltr_bkg_norm_2360 = n_OS_gltr_bkg_2360/lum_bkg_2360;

   h = (TH1*) f->Get("hReco_mass_OS_trtr_Cut4");
   Double_t n_OS_trtr_bkg_2360 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_OS_trtr_bkg_norm_2360 = n_OS_trtr_bkg_2360/lum_bkg_2360;

   cout<< "n_OS_glgl_bkg_2360 = " << n_OS_glgl_bkg_2360 << " n_OS_glgl_bkg_norm_2360 = " << n_OS_glgl_bkg_norm_2360 <<endl;
   cout<< "n_OS_gltr_bkg_2360 = " << n_OS_gltr_bkg_2360 << " n_OS_gltr_bkg_norm_2360 = " << n_OS_gltr_bkg_norm_2360 <<endl;
   cout<< "n_OS_trtr_bkg_2360 = " << n_OS_trtr_bkg_2360 << " n_OS_trtr_bkg_norm_2360 = " << n_OS_trtr_bkg_norm_2360 <<endl;

   //--
   //-- Same Sign
   //--
   cout<<endl<< "Same Sign" <<endl;

   // signal

   f = TFile::Open("histos_jPsi_900GeV_STARTUP.root");

   h = (TH1*) f->Get("hReco_mass_LS_glgl_Cut4");
   Double_t n_SS_glgl_sig_900 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_SS_glgl_sig_norm_900 = n_SS_glgl_sig_900/lum_sig_900;

   h = (TH1*) f->Get("hReco_mass_LS_gltr_Cut4");
   Double_t n_SS_gltr_sig_900 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_SS_gltr_sig_norm_900 = n_SS_gltr_sig_900/lum_sig_900;

   h = (TH1*) f->Get("hReco_mass_LS_trtr_Cut4");
   Double_t n_SS_trtr_sig_900 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_SS_trtr_sig_norm_900 = n_SS_trtr_sig_900/lum_sig_900;

   cout<< "n_SS_glgl_sig_900 = " << n_SS_glgl_sig_900 << " n_SS_glgl_sig_norm_900 = " << n_SS_glgl_sig_norm_900 <<endl;
   cout<< "n_SS_gltr_sig_900 = " << n_SS_gltr_sig_900 << " n_SS_gltr_sig_norm_900 = " << n_SS_gltr_sig_norm_900 <<endl;
   cout<< "n_SS_trtr_sig_900 = " << n_SS_trtr_sig_900 << " n_SS_trtr_sig_norm_900 = " << n_SS_trtr_sig_norm_900 <<endl;

   f = TFile::Open("histos_jPsi_2360GeV_STARTUP.root");

   h = (TH1*) f->Get("hReco_mass_LS_glgl_Cut4");
   Double_t n_SS_glgl_sig_2360 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_SS_glgl_sig_norm_2360 = n_SS_glgl_sig_2360/lum_sig_2360;

   h = (TH1*) f->Get("hReco_mass_LS_gltr_Cut4");
   Double_t n_SS_gltr_sig_2360 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_SS_gltr_sig_norm_2360 = n_SS_gltr_sig_2360/lum_sig_2360;

   h = (TH1*) f->Get("hReco_mass_LS_trtr_Cut4");
   Double_t n_SS_trtr_sig_2360 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_SS_trtr_sig_norm_2360 = n_SS_trtr_sig_2360/lum_sig_2360;

   cout<< "n_SS_glgl_sig_2360 = " << n_SS_glgl_sig_2360 << " n_SS_glgl_sig_norm_2360 = " << n_SS_glgl_sig_norm_2360 <<endl;
   cout<< "n_SS_gltr_sig_2360 = " << n_SS_gltr_sig_2360 << " n_SS_gltr_sig_norm_2360 = " << n_SS_gltr_sig_norm_2360 <<endl;
   cout<< "n_SS_trtr_sig_2360 = " << n_SS_trtr_sig_2360 << " n_SS_trtr_sig_norm_2360 = " << n_SS_trtr_sig_norm_2360 <<endl;

   // background

   f = TFile::Open("histos_ppMuXLoose_900GeV_STARTUP.root");

   h = (TH1*) f->Get("hReco_mass_LS_glgl_Cut4");
   Double_t n_SS_glgl_bkg_900 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_SS_glgl_bkg_norm_900 = n_SS_glgl_bkg_900/lum_bkg_900;

   h = (TH1*) f->Get("hReco_mass_LS_gltr_Cut4");
   Double_t n_SS_gltr_bkg_900 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_SS_gltr_bkg_norm_900 = n_SS_gltr_bkg_900/lum_bkg_900;

   h = (TH1*) f->Get("hReco_mass_LS_trtr_Cut4");
   Double_t n_SS_trtr_bkg_900 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_SS_trtr_bkg_norm_900 = n_SS_trtr_bkg_900/lum_bkg_900;

   cout<< "n_SS_glgl_bkg_900 = " << n_SS_glgl_bkg_900 << " n_SS_glgl_bkg_norm_900 = " << n_SS_glgl_bkg_norm_900 <<endl;
   cout<< "n_SS_gltr_bkg_900 = " << n_SS_gltr_bkg_900 << " n_SS_gltr_bkg_norm_900 = " << n_SS_gltr_bkg_norm_900 <<endl;
   cout<< "n_SS_trtr_bkg_900 = " << n_SS_trtr_bkg_900 << " n_SS_trtr_bkg_norm_900 = " << n_SS_trtr_bkg_norm_900 <<endl;

   f = TFile::Open("histos_ppMuXLoose_2360GeV_STARTUP.root");

   h = (TH1*) f->Get("hReco_mass_LS_glgl_Cut4");
   Double_t n_SS_glgl_bkg_2360 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_SS_glgl_bkg_norm_2360 = n_SS_glgl_bkg_2360/lum_bkg_2360;

   h = (TH1*) f->Get("hReco_mass_LS_gltr_Cut4");
   Double_t n_SS_gltr_bkg_2360 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_SS_gltr_bkg_norm_2360 = n_SS_gltr_bkg_2360/lum_bkg_2360;

   h = (TH1*) f->Get("hReco_mass_LS_trtr_Cut4");
   Double_t n_SS_trtr_bkg_2360 = h->Integral(h->GetXaxis()->FindBin(3.0), h->GetXaxis()->FindBin(3.2));
   Double_t n_SS_trtr_bkg_norm_2360 = n_SS_trtr_bkg_2360/lum_bkg_2360;

   cout<< "n_SS_glgl_bkg_2360 = " << n_SS_glgl_bkg_2360 << " n_SS_glgl_bkg_norm_2360 = " << n_SS_glgl_bkg_norm_2360 <<endl;
   cout<< "n_SS_gltr_bkg_2360 = " << n_SS_gltr_bkg_2360 << " n_SS_gltr_bkg_norm_2360 = " << n_SS_gltr_bkg_norm_2360 <<endl;
   cout<< "n_SS_trtr_bkg_2360 = " << n_SS_trtr_bkg_2360 << " n_SS_trtr_bkg_norm_2360 = " << n_SS_trtr_bkg_norm_2360 <<endl;

   //--
   //-- OS - SS
   //--
   cout<<endl<< "Same Sign - Opposite Sign" <<endl;

   // signal
   Double_t n_OSSS_glgl_sig_norm_900 = n_OS_glgl_sig_norm_900 - n_SS_glgl_sig_norm_900;
   Double_t n_OSSS_gltr_sig_norm_900 = n_OS_gltr_sig_norm_900 - n_SS_gltr_sig_norm_900;
   Double_t n_OSSS_trtr_sig_norm_900 = n_OS_trtr_sig_norm_900 - n_SS_trtr_sig_norm_900;
   Double_t n_OSSS_glgl_sig_norm_2360 = n_OS_glgl_sig_norm_2360 - n_SS_glgl_sig_norm_2360;
   Double_t n_OSSS_gltr_sig_norm_2360 = n_OS_gltr_sig_norm_2360 - n_SS_gltr_sig_norm_2360;
   Double_t n_OSSS_trtr_sig_norm_2360 = n_OS_trtr_sig_norm_2360 - n_SS_trtr_sig_norm_2360;
   cout<< "OS-SS sig 900 GeV:  glgl = " << n_OSSS_glgl_sig_norm_900 << " gltr = " << n_OSSS_gltr_sig_norm_900 << " trtr = " << n_OSSS_trtr_sig_norm_900 <<endl;
   cout<< "OS-SS sig 2360 GeV: glgl = " << n_OSSS_glgl_sig_norm_2360 << " gltr = " << n_OSSS_gltr_sig_norm_2360 << " trtr = " << n_OSSS_trtr_sig_norm_2360 <<endl;
   // bkg
   Double_t n_OSSS_glgl_bkg_norm_900 = n_OS_glgl_bkg_norm_900 - n_SS_glgl_bkg_norm_900;
   Double_t n_OSSS_gltr_bkg_norm_900 = n_OS_gltr_bkg_norm_900 - n_SS_gltr_bkg_norm_900;
   Double_t n_OSSS_trtr_bkg_norm_900 = n_OS_trtr_bkg_norm_900 - n_SS_trtr_bkg_norm_900;
   Double_t n_OSSS_glgl_bkg_norm_2360 = n_OS_glgl_bkg_norm_2360 - n_SS_glgl_bkg_norm_2360;
   Double_t n_OSSS_gltr_bkg_norm_2360 = n_OS_gltr_bkg_norm_2360 - n_SS_gltr_bkg_norm_2360;
   Double_t n_OSSS_trtr_bkg_norm_2360 = n_OS_trtr_bkg_norm_2360 - n_SS_trtr_bkg_norm_2360;
   cout<< "OS-SS bkg 900 GeV:  glgl = " << n_OSSS_glgl_bkg_norm_900 << " gltr = " << n_OSSS_gltr_bkg_norm_900 << " trtr = " << n_OSSS_trtr_bkg_norm_900 <<endl;
   cout<< "OS-SS bkg 2360 GeV: glgl = " << n_OSSS_glgl_bkg_norm_2360 << " gltr = " << n_OSSS_gltr_bkg_norm_2360 << " trtr = " << n_OSSS_trtr_bkg_norm_2360 <<endl;

   //--
   //-- Plot
   //--

   TH1F* h_OS_sig_900  = new TH1F("h_OS_sig_900", "OS J/#psi 900 GeV", 3, 0, 3);
   h_OS_sig_900->SetMarkerStyle(24);
   h_OS_sig_900->SetMarkerColor(2);
   h_OS_sig_900->GetXaxis()->SetBinLabel(1, "glgl");
   h_OS_sig_900->GetXaxis()->SetBinLabel(2, "gltr");
   h_OS_sig_900->GetXaxis()->SetBinLabel(3, "trtr");
   h_OS_sig_900->SetBinContent(1, n_OS_glgl_sig_norm_900);
   h_OS_sig_900->SetBinContent(2, n_OS_gltr_sig_norm_900);
   h_OS_sig_900->SetBinContent(3, n_OS_trtr_sig_norm_900);
   //new TCanvas;
   //h_OS_sig_900->Draw("p0");

   TH1F* h_SS_sig_900  = new TH1F("h_SS_sig_900", "SS J/#psi 900 GeV", 3, 0, 3);
   h_SS_sig_900->SetMarkerStyle(24);
   h_SS_sig_900->SetMarkerColor(4);
   h_SS_sig_900->GetXaxis()->SetBinLabel(1, "glgl");
   h_SS_sig_900->GetXaxis()->SetBinLabel(2, "gltr");
   h_SS_sig_900->GetXaxis()->SetBinLabel(3, "trtr");
   h_SS_sig_900->SetBinContent(1, n_SS_glgl_sig_norm_900);
   h_SS_sig_900->SetBinContent(2, n_SS_gltr_sig_norm_900);
   h_SS_sig_900->SetBinContent(3, n_SS_trtr_sig_norm_900);
   //new TCanvas;
   //h_SS_sig_900->Draw("p0");

   TH1F* h_OSSS_sig_900  = new TH1F("h_OSSS_sig_900", "OS - SS J/#psi 900 GeV", 3, 0, 3);
   h_OSSS_sig_900->SetMarkerStyle(20);
   h_OSSS_sig_900->SetMarkerColor(2);
   h_OSSS_sig_900->GetXaxis()->SetBinLabel(1, "glgl");
   h_OSSS_sig_900->GetXaxis()->SetBinLabel(2, "gltr");
   h_OSSS_sig_900->GetXaxis()->SetBinLabel(3, "trtr");
   h_OSSS_sig_900->SetBinContent(1, n_OSSS_glgl_sig_norm_900);
   h_OSSS_sig_900->SetBinContent(2, n_OSSS_gltr_sig_norm_900);
   h_OSSS_sig_900->SetBinContent(3, n_OSSS_trtr_sig_norm_900);

   new TCanvas;
   h_OS_sig_900->Draw("p0");
   h_SS_sig_900->Draw("p0same");
   h_OSSS_sig_900->Draw("p0same");
   legend = new TLegend(0.20,0.50, 0.40,0.75);
   legend->SetFillColor(0);
   legend->SetBorderSize(0);
   legend->AddEntry(h_OS_sig_900, "OS");
   legend->AddEntry(h_SS_sig_900, "SS");
   legend->AddEntry(h_OSSS_sig_900, "OS-SS");
   legend->Draw();
   gPad->SaveAs("h_sig_900.png");

   TH1F* h_OS_bkg_900  = new TH1F("h_OS_bkg_900", "OS bkg 900 GeV", 3, 0, 3);
   h_OS_bkg_900->SetMarkerStyle(24);
   h_OS_bkg_900->SetMarkerColor(2);
   h_OS_bkg_900->GetXaxis()->SetBinLabel(1, "glgl");
   h_OS_bkg_900->GetXaxis()->SetBinLabel(2, "gltr");
   h_OS_bkg_900->GetXaxis()->SetBinLabel(3, "trtr");
   h_OS_bkg_900->SetBinContent(1, n_OS_glgl_bkg_norm_900);
   h_OS_bkg_900->SetBinContent(2, n_OS_gltr_bkg_norm_900);
   h_OS_bkg_900->SetBinContent(3, n_OS_trtr_bkg_norm_900);
   //new TCanvas;
   //h_OS_bkg_900->Draw("p0");

   TH1F* h_SS_bkg_900  = new TH1F("h_SS_bkg_900", "SS bkg 900 GeV", 3, 0, 3);
   h_SS_bkg_900->SetMarkerStyle(24);
   h_SS_bkg_900->SetMarkerColor(4);
   h_SS_bkg_900->GetXaxis()->SetBinLabel(1, "glgl");
   h_SS_bkg_900->GetXaxis()->SetBinLabel(2, "gltr");
   h_SS_bkg_900->GetXaxis()->SetBinLabel(3, "trtr");
   h_SS_bkg_900->SetBinContent(1, n_SS_glgl_bkg_norm_900);
   h_SS_bkg_900->SetBinContent(2, n_SS_gltr_bkg_norm_900);
   h_SS_bkg_900->SetBinContent(3, n_SS_trtr_bkg_norm_900);
   //new TCanvas;
   //h_SS_bkg_900->Draw("p0");

   TH1F* h_OSSS_bkg_900  = new TH1F("h_OSSS_bkg_900", "OS - SS bkg 900 GeV", 3, 0, 3);
   h_OSSS_bkg_900->SetMarkerStyle(20);
   h_OSSS_bkg_900->SetMarkerColor(2);
   h_OSSS_bkg_900->GetXaxis()->SetBinLabel(1, "glgl");
   h_OSSS_bkg_900->GetXaxis()->SetBinLabel(2, "gltr");
   h_OSSS_bkg_900->GetXaxis()->SetBinLabel(3, "trtr");
   h_OSSS_bkg_900->SetBinContent(1, n_OSSS_glgl_bkg_norm_900);
   h_OSSS_bkg_900->SetBinContent(2, n_OSSS_gltr_bkg_norm_900);
   h_OSSS_bkg_900->SetBinContent(3, n_OSSS_trtr_bkg_norm_900);

   new TCanvas;
   h_OS_bkg_900->Draw("p0");
   h_SS_bkg_900->Draw("p0same");
   h_OSSS_bkg_900->Draw("p0same");
   legend = new TLegend(0.20,0.50, 0.40,0.75);
   legend->SetFillColor(0);
   legend->SetBorderSize(0);
   legend->AddEntry(h_OS_bkg_900, "OS");
   legend->AddEntry(h_SS_bkg_900, "SS");
   legend->AddEntry(h_OSSS_bkg_900, "OS-SS");
   legend->Draw();
   gPad->SaveAs("h_bkg_900.png");

   //-----------------
   //-- 2360 GeV
   //-----------------

   TH1F* h_OS_sig_2360  = new TH1F("h_OS_sig_2360", "OS J/#psi 2360 GeV", 3, 0, 3);
   h_OS_sig_2360->SetMarkerStyle(24);
   h_OS_sig_2360->SetMarkerColor(2);
   h_OS_sig_2360->GetXaxis()->SetBinLabel(1, "glgl");
   h_OS_sig_2360->GetXaxis()->SetBinLabel(2, "gltr");
   h_OS_sig_2360->GetXaxis()->SetBinLabel(3, "trtr");
   h_OS_sig_2360->SetBinContent(1, n_OS_glgl_sig_norm_2360);
   h_OS_sig_2360->SetBinContent(2, n_OS_gltr_sig_norm_2360);
   h_OS_sig_2360->SetBinContent(3, n_OS_trtr_sig_norm_2360);
   //new TCanvas;
   //h_OS_sig_2360->Draw("p0");

   TH1F* h_SS_sig_2360  = new TH1F("h_SS_sig_2360", "SS J/#psi 2360 GeV", 3, 0, 3);
   h_SS_sig_2360->SetMarkerStyle(24);
   h_SS_sig_2360->SetMarkerColor(4);
   h_SS_sig_2360->GetXaxis()->SetBinLabel(1, "glgl");
   h_SS_sig_2360->GetXaxis()->SetBinLabel(2, "gltr");
   h_SS_sig_2360->GetXaxis()->SetBinLabel(3, "trtr");
   h_SS_sig_2360->SetBinContent(1, n_SS_glgl_sig_norm_2360);
   h_SS_sig_2360->SetBinContent(2, n_SS_gltr_sig_norm_2360);
   h_SS_sig_2360->SetBinContent(3, n_SS_trtr_sig_norm_2360);
   //new TCanvas;
   //h_SS_sig_2360->Draw("p0");

   TH1F* h_OSSS_sig_2360  = new TH1F("h_OSSS_sig_2360", "OS - SS J/#psi 2360 GeV", 3, 0, 3);
   h_OSSS_sig_2360->SetMarkerStyle(20);
   h_OSSS_sig_2360->SetMarkerColor(2);
   h_OSSS_sig_2360->GetXaxis()->SetBinLabel(1, "glgl");
   h_OSSS_sig_2360->GetXaxis()->SetBinLabel(2, "gltr");
   h_OSSS_sig_2360->GetXaxis()->SetBinLabel(3, "trtr");
   h_OSSS_sig_2360->SetBinContent(1, n_OSSS_glgl_sig_norm_2360);
   h_OSSS_sig_2360->SetBinContent(2, n_OSSS_gltr_sig_norm_2360);
   h_OSSS_sig_2360->SetBinContent(3, n_OSSS_trtr_sig_norm_2360);

   new TCanvas;
   h_OS_sig_2360->Draw("p0");
   h_SS_sig_2360->Draw("p0same");
   h_OSSS_sig_2360->Draw("p0same");
   legend = new TLegend(0.20,0.50, 0.40,0.75);
   legend->SetFillColor(0);
   legend->SetBorderSize(0);
   legend->AddEntry(h_OS_sig_2360, "OS");
   legend->AddEntry(h_SS_sig_2360, "SS");
   legend->AddEntry(h_OSSS_sig_2360, "OS-SS");
   legend->Draw();
   gPad->SaveAs("h_sig_2360.png");

   TH1F* h_OS_bkg_2360  = new TH1F("h_OS_bkg_2360", "OS bkg 2360 GeV", 3, 0, 3);
   h_OS_bkg_2360->SetMarkerStyle(24);
   h_OS_bkg_2360->SetMarkerColor(2);
   h_OS_bkg_2360->GetXaxis()->SetBinLabel(1, "glgl");
   h_OS_bkg_2360->GetXaxis()->SetBinLabel(2, "gltr");
   h_OS_bkg_2360->GetXaxis()->SetBinLabel(3, "trtr");
   h_OS_bkg_2360->SetBinContent(1, n_OS_glgl_bkg_norm_2360);
   h_OS_bkg_2360->SetBinContent(2, n_OS_gltr_bkg_norm_2360);
   h_OS_bkg_2360->SetBinContent(3, n_OS_trtr_bkg_norm_2360);
   //new TCanvas;
   //h_OS_bkg_2360->Draw("p0");

   TH1F* h_SS_bkg_2360  = new TH1F("h_SS_bkg_2360", "SS bkg 2360 GeV", 3, 0, 3);
   h_SS_bkg_2360->SetMarkerStyle(24);
   h_SS_bkg_2360->SetMarkerColor(4);
   h_SS_bkg_2360->GetXaxis()->SetBinLabel(1, "glgl");
   h_SS_bkg_2360->GetXaxis()->SetBinLabel(2, "gltr");
   h_SS_bkg_2360->GetXaxis()->SetBinLabel(3, "trtr");
   h_SS_bkg_2360->SetBinContent(1, n_SS_glgl_bkg_norm_2360);
   h_SS_bkg_2360->SetBinContent(2, n_SS_gltr_bkg_norm_2360);
   h_SS_bkg_2360->SetBinContent(3, n_SS_trtr_bkg_norm_2360);
   //new TCanvas;
   //h_SS_bkg_2360->Draw("p0");

   TH1F* h_OSSS_bkg_2360  = new TH1F("h_OSSS_bkg_2360", "OS - SS bkg 2360 GeV", 3, 0, 3);
   h_OSSS_bkg_2360->SetMarkerStyle(20);
   h_OSSS_bkg_2360->SetMarkerColor(2);
   h_OSSS_bkg_2360->GetXaxis()->SetBinLabel(1, "glgl");
   h_OSSS_bkg_2360->GetXaxis()->SetBinLabel(2, "gltr");
   h_OSSS_bkg_2360->GetXaxis()->SetBinLabel(3, "trtr");
   h_OSSS_bkg_2360->SetBinContent(1, n_OSSS_glgl_bkg_norm_2360);
   h_OSSS_bkg_2360->SetBinContent(2, n_OSSS_gltr_bkg_norm_2360);
   h_OSSS_bkg_2360->SetBinContent(3, n_OSSS_trtr_bkg_norm_2360);

   new TCanvas;
   h_OS_bkg_2360->Draw("p0");
   h_SS_bkg_2360->Draw("p0same");
   h_OSSS_bkg_2360->Draw("p0same");
   legend = new TLegend(0.20,0.50, 0.40,0.75);
   legend->SetFillColor(0);
   legend->SetBorderSize(0);
   legend->AddEntry(h_OS_bkg_2360, "OS");
   legend->AddEntry(h_SS_bkg_2360, "SS");
   legend->AddEntry(h_OSSS_bkg_2360, "OS-SS");
   legend->Draw();
   gPad->SaveAs("h_bkg_2360.png");
}
