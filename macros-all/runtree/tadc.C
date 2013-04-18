void tadc(Int_t run)
{
   TCut runno = Form("run==%d", run);
   std::string cname = nextcan(Form("run%d",run));
   TCanvas* can = new TCanvas(cname.c_str(), cname.c_str(), 1000,500);
   can->Divide(3,2);
   // can->cd(1); t->Draw("t1>>ht1", runno+"t1>0&&t1<6000");
   // can->cd(2); t->Draw("t2>>ht2", runno+"t2>0&&t2<6000");
   // can->cd(3); t->Draw("t3>>ht3", runno+"t3>0&&t3<6000");
   // can->cd(4); t->Draw("a1>>ha1", runno+"a1>1");
   // can->cd(5); t->Draw("a2>>ha2", runno+"a2>1");
   // can->cd(6); t->Draw("a3>>ha3", runno+"a3>1");
   can->cd(1); t->Draw("t1", runno+"t1>0&&t1<6000");
   can->cd(2); t->Draw("t2", runno+"t2>0&&t2<6000");
   can->cd(3); t->Draw("t3", runno+"t3>0&&t3<6000");
   can->cd(4); t->Draw("a1", runno+"a1>1");
   can->cd(5); t->Draw("a2", runno+"a2>1");
   can->cd(6); t->Draw("a3", runno+"a3>1");
   can->cd();
}

void run(Int_t run) {tadc(run);}

// class TTreeSettins {
// public:
//    TreeSettings() {
//       cout<< "Set styles and colors for tree t" <<endl;
//       t->SetMarkerStyle(6);
//       t->SetMarkerColor(2);
//       t->SetFillStyle(3001);
//       t->SetFillColor(38);
//    }
// } tTreeSettings;
//
void ttree() {
   cout<< "Set styles and colors for tree t" <<endl;
   t->SetMarkerStyle(6);
   t->SetMarkerColor(2);
   t->SetFillStyle(3001);
   t->SetFillColor(38);
}
