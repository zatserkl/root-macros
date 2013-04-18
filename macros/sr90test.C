void sr90test()
{
   // Double_t hv[]= {1000,   600,  700,  800,     900,     1000,    1100,    1200,    1000};
   // Double_t n1[]= {35436,  0,    485,  17206,   31313,   35046,   37844,   59100,   34364};
   // Double_t n2[]= {38696,  5084, 30017,33761,   36086,   38201,   47269,   88931,   37460};
   // Double_t cc[]= {35119,  0,    485,  17204,   31307,   34752,   35537,   41405,   34072};
   // const Int_t np = 9;
   Double_t hv[]= {1000,   700,  800,     900,     1000,    1100,    1200,    1000};
   Double_t n1[]= {35436,  485,  17206,   31313,   35046,   37844,   59100,   34364};
   Double_t n2[]= {38696,  30017,33761,   36086,   38201,   47269,   88931,   37460};
   Double_t cc[]= {35119,  485,  17204,   31307,   34752,   35537,   41405,   34072};
   const Int_t np = 8;

   Double_t ehv[np];
   Double_t en1[np];
   Double_t en2[np];
   Double_t ecc[np];

   for (int i=0; i<np; ++i) {
      ehv[i] = 0;
      en1[i] = TMath::Sqrt(n1[i]);
      en2[i] = TMath::Sqrt(n2[i]);
      ecc[i] = TMath::Sqrt(cc[i]);
   }

   TGraphErrors* gr_n1 = new TGraphErrors(np, hv,n1, ehv,en1);
   TGraphErrors* gr_n2 = new TGraphErrors(np, hv,n2, ehv,en2);
   TGraphErrors* gr_cc = new TGraphErrors(np, hv,cc, ehv,ecc);

   gr_n1->SetMarkerStyle(20);
   gr_n2->SetMarkerStyle(20);
   // gr_cc->SetMarkerStyle(24);
   gr_cc->SetMarkerStyle(5);
   // gr_cc->SetMarkerSize(1.50);
   gr_cc->SetMarkerSize(2.00);

   gr_n1->SetMarkerColor(4);
   gr_n2->SetMarkerColor(8);
   gr_cc->SetMarkerColor(2);

   gr_n1->SetFillColor(0);
   gr_n2->SetFillColor(0);
   gr_cc->SetFillColor(0);

   gr_n1->SetLineColor(0);
   gr_n2->SetLineColor(0);
   gr_cc->SetLineColor(0);

   TMultiGraph *mg = new TMultiGraph();
   mg->SetNameTitle("mg", "HV dependence with source ^{90}Sr;HV, V;counts per 10 s   ");
   mg->Add(gr_n1);
   mg->Add(gr_n2);
   mg->Add(gr_cc);

   new TCanvas;
   gStyle->SetTitleBorderSize(1);   
   gPad->SetLogy();
   gPad->SetGrid(0,0);
   mg->Draw("ap");

   TLegend* legend = new TLegend(0.65,0.30, 0.85,0.60);
   legend->SetFillColor(0);
   legend->AddEntry(gr_n1, "PMT #1");
   legend->AddEntry(gr_n2, "PMT #2");
   legend->AddEntry(gr_cc, "CC");
   legend->Draw();

   TArrow* arrow = new TArrow(1000,2000, 1000,400);
   arrow->SetLineWidth(2);
   arrow->SetLineColor(2);
   arrow->SetAngle(45);
   arrow->Draw();
}
