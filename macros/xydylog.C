void xydylog(char* fname)
{
   //gROOT->Reset();

   //char* fname = "/tmp_root/759/niu_1/zatserkl/L1MET/L1MET.dat";
   //char* fname = "L1RecoMET.dat";
   ifstream fin(fname, ios_base::in);   // declare and initilize at once
   
   if ( (fin.rdstate() & ifstream::failbit ) != 0 )
   {
      cerr << "Error opening " << fname << endl;
      exit(1);
   }

   const Int_t dim = 1000;
   Float_t x[dim], y[dim],dy[dim];
   Int_t np = 0;
   
   char line[256];
   Int_t len = sizeof(line);

   stringstream strline;
   strline << ends;     // to avoid storing all input in strigstream

   while (fin.getline(line,len))
   {
      //cout << line << endl;
      
      strline << line;
      strline >> x[np] >> y[np] >> dy[np];
      
      if (strline.fail()) {
         //   Note: in case of fail no symbols are extracted from stream!
         strline.clear();       // clear failbit to work with stream
         //cout << "Ignored line:\n" << strline.str().c_str() << endl;
         strline.ignore(len);   // get rid failed stream
         strline.clear();       // clear eofbit after ignore() or getline()
         continue;
      }
      
      if (!strline.eof()) {
         cout << "there is some extra data in line" << endl;
         strline.ignore(len);   // get rid failed stream
         strline.clear();       // clear eofbit after ignore() or getline()
      }
      else strline.clear();     // clear eofbit
           
      cout << "   " << np << ": " << x[np] <<" "<< y[np]   <<" "<< dy[np] << endl;
      np++;
   }
   cout << "Read " << np << " lines" << endl;
   
   fin.close();
   
   //cout << "Contents of strline:\n" << strline.str().c_str() << endl;

   // Plot
   if (np == 0) {
      cout << "No points were read. Stop" << endl;
      return;
   }

   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);

   c1->SetFillColor(10);
   c1->SetGridx();
   c1->SetGridy();

   c1->SetLogy();
      
   TGraph *gr1 = new TGraph(np,x,y);
   gr1->SetMarkerColor(4);
   gr1->SetMarkerStyle(20);
   gr1->SetTitle("y vs x");
   gr1->Draw("AP");

   //gr1->GetHistogram()->SetXTitle("MET cut, GeV");
   //gr1->GetHistogram()->SetYTitle("Efficiency");
   c1->Modified();   
   
   //leg = new TLegend(0.2,0.2, 0.4,0.3);
   //leg->AddEntry(gr1,"L1","p");
   //leg->AddEntry(gr2,"Reco","p");
   //leg->Draw();
}
