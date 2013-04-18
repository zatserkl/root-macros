void axis_invert()
{
   // Rene Brun
   TF1 *f1 = new TF1("f1","-x",-3,3);
   TH1F *h = new TH1F("h","test",100,-3,3);
   h->FillRandom("gaus",5000);
   TCanvas *c1 = new TCanvas("c1");
   h->Draw();
   c1->Update();
   h->GetXaxis()->SetLabelOffset(99);
   Double_t xmin = c1->GetUxmin();
   Double_t xmax = c1->GetUxmax();
   Double_t ymin = c1->GetUymin();
   TGaxis *axis = new TGaxis(xmin,ymin,xmax,ymin,"f1");
   axis->Draw();
}
