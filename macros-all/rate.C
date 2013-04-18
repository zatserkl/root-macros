#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include <iostream>
#include <string>
#include <map>

using std::cout;     using std::endl;

void fitg(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);

namespace Tree
{
   Int_t    n;
   Float_t  thres;
   
   void clear()
   {
      n        = 0;
      thres        = 0;
   }
   
   void book(TTree* tree) {
      tree->Branch("n",          &n,         "n/I");
      tree->Branch("thres",          &thres,         "thres");          // By default equivalent to "a/F"
   }
   void connect(TTree* tree)                                // need for event-by-event analysis
   {   
      // connects tree buffers with variables to use for event-by-event analysis
      tree->SetBranchAddress("n",            &n);
      tree->SetBranchAddress("thres",            &thres);
   }
}  // namespace Tree

bool peak(Float_t x[], Float_t y[], Float_t thres, Float_t xmin, Float_t xmax, Float_t window
      , Int_t& isigmin, Int_t& isigmax
      )
{
   Int_t trig_i = -1;

   for (int i=isigmin; i<1024; ++i)
   {
      if (x[i] < xmin) continue;
      if (x[i] > xmax) break;

      if (y[i] > thres) {
	 trig_i = i;
	 break;
      }
      else continue;
   }
   if (trig_i < 0) return false;

   // the threshold was exceeded. Look for the maximum in nearest 5 ns

   isigmin = trig_i;

   // find the maximum

   Float_t xmax_local = x[trig_i] + window;
   Int_t maximum_i = trig_i;
   for (int i=maximum_i+1; i<1024; ++i)
   {
      if (x[i] > xmax) break;	    // gloabal range
      if (x[i] > xmax_local) break; // 5 ns window to find the new maximum
      if (y[i] > y[maximum_i]) {
      	 maximum_i = i;
      	 xmax_local = x[i] + window;	  // start 5 ns window to find new maximum
      	 //cout<< "new maximum x = " << x[maximum_i] << " xmax_local = " << xmax_local <<endl;
      }
   }

   // find the tail

   Float_t ymin = 0.1 * y[maximum_i];
   if (thres < ymin) ymin = thres;

   isigmax = maximum_i;
   for (int i=maximum_i; i<1024; ++i)
   {
      if (x[i] > xmax) break;		  // gloabal range
      if (x[i] > xmax_local) break;	  // 5 ns window to find the new maximum
      if (y[i] < 0 or y[i] < ymin) break; // pulse tail
      isigmax = i;
   }

   // find the leading edge

   for (int i=trig_i; i>=0; --i) {
      if (x[i] < xmin) break;		  // gloabal range
      if (y[i] < 0 or y[i] < ymin) break; // pulse tail
      isigmin = i;
   }

   // // fit (for debug)
   // cout<< "fit from " << x[isigmin] << " to " << x[isigmax] <<endl;
   // fitg(x[isigmin], x[isigmax], "+");

   // cout<< "peak: thres " << thres << " xmin " << xmin << " xmax " << xmax << " isigmin " << isigmin << " isigmax " << isigmax <<endl;

   return true;
}

Int_t event_peaks(Float_t x[], Float_t y[], Float_t thres, Float_t trig, Float_t gate, Float_t window=5., Float_t pulse_width_min=5.)
{
   Int_t npeaks = 0;

   Int_t isigmin = 0;
   Int_t isigmax = 0;

   Float_t xmin = trig;
   Float_t xmax = trig + gate;
   while (peak(x, y, thres, xmin, xmax, window, isigmin, isigmax))
   {
      // cout<< "event_peaks: thres " << thres << " xmin " << xmin << " xmax " << xmax << " isigmin " << isigmin << " isigmax " << isigmax <<endl;
      if (pulse_width_min == 0) npeaks++;
      else if (x[isigmax] - x[isigmin] > pulse_width_min) npeaks++;   // ignore spikes if pulse_width_min has been specified

      if (x[isigmax+1] > x[isigmin]+2.*window) xmin = x[isigmax+1];
      else xmin = x[isigmin+1] + 2.*window;
      isigmin = isigmax+1;
   }

   return npeaks;
}

Float_t rate_thres(TTree* pulse, Int_t chan, Float_t thres, Float_t bkg, Float_t trig, Float_t gate, Int_t evt1=0, Int_t evt2=-1)
{
   // NB: bkg is background for the inverted pulse

   const Int_t nsamples = 1024;

   Float_t window = 5.;		       // ns, window to search for the local maximum
   //-- Float_t pulse_width_min = 0;   // min pulse width = 5 ns
   // Float_t pulse_width_min = 5.;       // min pulse width = 5 ns
   Float_t pulse_width_min = 3.;       // min pulse width = 5 ns

   Float_t x[nsamples];
   Float_t y[nsamples];
   Int_t np = 0;

   Int_t board = (chan-1)/4 + 1;
   //Int_t channel = (chan-1) % 4 + 1;

   std::string tleaf = Form("t%d",board);
   std::string vleaf = Form("c%d",chan);
   if (pulse->GetLeaf(tleaf.c_str()))
   {
      // found time brach (new name), look for the voltage branch
      if (!pulse->GetLeaf(vleaf.c_str())) {
	 cout<< "Could not find voltage branch" <<endl;
	 return 0;
      }
   }
   else {
      // look for the old branch names
      Int_t chan_1_4 = (chan-1) % 4 + 1;
      tleaf = Form("b%d_t", board);
      vleaf = Form("b%d_c%d",board,chan_1_4);
      if (!pulse->GetLeaf(tleaf.c_str())) {
      	 cout<< "Could not find time branch" <<endl;
      	 return 0;
      }
      if (!pulse->GetLeaf(tleaf.c_str())) {
      	 cout<< "Could not find voltage branch" <<endl;
      	 return 0;
      }
   }

   Float_t* v = 0;
   Float_t* t = 0;

   pulse->SetBranchStatus("*", 0);
   pulse->SetBranchStatus(tleaf.c_str(), 1);
   pulse->SetBranchStatus(vleaf.c_str(), 1);

   static TH1F* h_npeaks;
   if (!h_npeaks) {
      h_npeaks = new TH1F("h_npeaks", "The No. of pulses per event", 500, 0, 50);
      //h_npeaks->SetDirectory(0);
   }
   else h_npeaks->Reset();

   if (evt2 < evt1) evt2 = pulse->GetEntries()-1;

   for (Int_t treeNumber=-1, evt=evt1; evt<=evt2; ++evt)    // NB init of the treeNumber
   {
      //-- if (evt and evt % 1000 == 0) cout<< "--- evt = " << evt <<endl;

      if (pulse->LoadTree(evt) < 0) break;
      if (pulse->GetTreeNumber() != treeNumber)
      {
      	 // update addresses after loading new file in the chain
	 t = (Float_t*) pulse->GetLeaf(tleaf.c_str())->GetValuePointer();
	 v = (Float_t*) pulse->GetLeaf(vleaf.c_str())->GetValuePointer();
      }
      pulse->GetEntry(evt);

      //cout<< "t[0] " << t[0] << " v[0] " << v[0] <<endl;

      np = 0;
      for (int i=0; i<nsamples; ++i)
      {
	 x[np] = t[i];
	 y[np] = -1.*v[i] - bkg;
	 ++np;
      }

      Int_t npeaks = event_peaks(x, y, thres, trig, gate, window, pulse_width_min);
      h_npeaks->Fill(npeaks);
   }
   pulse->SetBranchStatus("*", 1);

   cout<< "thres = " << thres << "\t h_npeaks->GetMean() = " << h_npeaks->GetMean() <<endl;

   return h_npeaks->GetMean();
}

TGraph* rate(TTree* pulse, Int_t chan, Float_t bkg, Float_t trig, Float_t gate, Int_t evt1=0, Int_t evt2=-1, const char* title="rate")
{
   Float_t thres1 = 0.001;
   //-- Float_t thres2 = 0.100;
   Float_t thres2 = 0.200;
   //-- Float_t dthres = 0.001;
   Float_t dthres = 0.0005;

   Float_t threshold[1000];
   Float_t npeaks[1000];
   Int_t np = 0;

   Float_t thres = thres1;
   while (thres < thres2) {
      threshold[np] = thres;
      npeaks[np] = rate_thres(pulse,chan, thres,bkg, trig,gate, evt1,evt2);
      //-- if (npeaks[np] < 0.01) break;	     // stop when rate dropped below 0.01 for the gate (typical gate ~200)
      if (npeaks[np] < 0.001) break;	     // stop when rate dropped below 0.01 for the gate (typical gate ~200)
      ++np;
      thres += dthres;
   }

   TGraph* gr_thres = new TGraph(np, threshold, npeaks);
   gr_thres->SetNameTitle("gr_thres", title);
   gr_thres->SetMarkerStyle(20);
   gr_thres->SetMarkerSize(0.5);
   //gr_thres->SetMarkerStyle(7);
   gr_thres->SetMarkerColor(2);
   gr_thres->GetXaxis()->SetTitle("threshold, V");
   gr_thres->GetYaxis()->SetTitle(Form("pulses in %0.0f ns",gate));

   new TCanvas;
   gPad->SetLogy();
   gr_thres->Draw("ap");

   return gr_thres;
}

TGraph* rate(const char* ifname, Int_t chan, Float_t bkg, Float_t trig, Float_t gate, Int_t evt1=0, Int_t evt2=-1)
{
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return 0;
   }

   TTree* pulse = (TTree*) ifile->Get("p");
   if (!pulse) pulse = (TTree*) ifile->Get("pulse");
   if (!pulse) {
      cout<< "Could not find pulse tree" <<endl;
      return 0;
   }

   return rate(pulse, chan, bkg, trig, gate, evt1, evt2, Form("rate for %s", ifname));
}

// TGraph* differentiate(const TGraph* gr)	//-- mistake: forgot to divide by dx
// {
//    TGraph* gr_diff = (TGraph*) gr->Clone(Form("diff_%s",gr->GetName()));
//    gr_diff->SetTitle(Form("diff %s",gr->GetName()));
//    for (int i=0; i<gr->GetN()-1; ++i) {
//       gr_diff->SetPoint(i, gr->GetX()[i], gr->GetY()[i] - gr->GetY()[i+1]);
//    }
//    gr_diff->Set(gr->GetN()-1);
//    return gr_diff;
// }
TGraph* differentiate_down(const TGraph* gr)
{
   TGraph* gr_diff = (TGraph*) gr->Clone(Form("diff_%s",gr->GetName()));
   gr_diff->SetTitle(Form("diff %s",gr->GetName()));
   for (int i=0; i<gr->GetN()-1; ++i) {
      Double_t dx = gr->GetX()[i+1] - gr->GetX()[i];
      Double_t dy = -1. * (gr->GetY()[i+1] - gr->GetY()[i]);   // NB: *(-1)
      const Double_t eps = 1e-7;
      Double_t dydx = TMath::Abs(dx) > eps? dy/dx: dy/eps;
      gr_diff->SetPoint(i, gr->GetX()[i], dydx);
   }
   gr_diff->Set(gr->GetN()-1);
   return gr_diff;
}

// void rate_sig_bkg(const char* ifname)
void rate_sig_bkg(const char* ifname, Int_t chan, Float_t bkg, Float_t trig_sig, Float_t trig_bkg, Float_t gate, Int_t evt1=0, Int_t evt2=-1)
{
   // TGraph* gr_sig = rate(ifname, 2, 0.007, 150,200);
   // TGraph* gr_bkg = rate(ifname, 2, 0.007, 50,200);
   TGraph* gr_sig = rate(ifname, chan, bkg, trig_sig,gate, evt1,evt2);
   TGraph* gr_bkg = rate(ifname, chan, bkg, trig_bkg,gate, evt1,evt2);
   gr_bkg->SetMarkerColor(9);

   TGraph* gr_sig_diff = differentiate_down(gr_sig);
   TGraph* gr_bkg_diff = differentiate_down(gr_bkg);
   new TCanvas;
   gr_sig_diff->Draw("ap");
   new TCanvas;
   gr_bkg_diff->Draw("ap");

   Float_t x[1000], ydiff[1000];
   for (int i=0; i<gr_sig->GetN(); ++i) {
      x[i] = gr_sig->GetX()[i];
      ydiff[i] = i < gr_bkg->GetN()? gr_sig->GetY()[i] - gr_bkg->GetY()[i]: gr_sig->GetY()[i];
      if (ydiff[i] < 0) ydiff[i] = 0;
   }

   TGraph* gr_diff = new TGraph(gr_sig->GetN(), x, ydiff);
   gr_diff->SetNameTitle("gr_diff", Form("signal - bkg for %s", ifname));
   gr_diff->SetMarkerStyle(24);
   gr_diff->SetMarkerSize(0.5);
   gr_diff->SetMarkerColor(2);

   new TCanvas;
   gr_diff->Draw("ap");
}

//------------------------

#include <TROOT.h>
#include <TGraph.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH1.h>
#include <TText.h>

#include <iostream>

using std::cout;     using std::endl;

void show_peaks()
{
   const char* ifname = "direct_25x15_1x1_69.5V_3x3_72.0V_23cm_trig_110_1.dat.root-rate_graphs.root";
   TFile* ifile = TFile::Open(ifname);

   TGraph* gr = (TGraph*) ifile->Get("gr_thres");
   new TCanvas;
   gr->Draw("ap");

   TGraph* gr_diff = (TGraph*) ifile->Get("diff_gr_thres");
   new TCanvas;
   gr_diff->Draw("ap");

   TClonesArray* graphs = new TClonesArray("TGraph");
   TGraph* graph = 0;

   Double_t* x = gr_diff->GetX();
   Double_t* y = gr_diff->GetY();
   Double_t xmin = 0;
   Double_t xmax = 0;

   // while (xmax < gr_diff->GetX()[gr_diff->GetN()-1])
   // {
   //    Int_t ifirst = 0;
   //    Int_t npoints = 0;
   //    for (int i=0; i<gr_diff->GetN(); ++i)
   //    {
   //       if (x[i] < xmin) continue;
   //       if (x[i] > xmax) break;
   //       if (npoints == 0) ifirst = i;
   //       ++npoints;
   //    }
   //    graph = new((*graphs)[graphs->GetLast()+1]) TGraph(npoints, &x[ifirst], &y[ifirst]); 
   //    graph->SetNameTitle(Form("graph%d",graphs->GetLast()+1), Form("xmin = %f xmax = %f", xmin,xmax));
   //    graph->SetMarkerStyle(20);
   //    graph->SetMarkerColor(2);

   //    Double_t dx = xmax - xmin;
   //    xmin = xmax;
   //    xmax += dx;
   // }

   Int_t np = 0;
   Int_t ifirst = 0;
   
   xmin = 0.0;
   xmax = 0.005;
   np = 0;
   for (int i=0; i<gr_diff->GetN(); ++i)
   {
      if (x[i] < xmin) continue;
      if (x[i] > xmax) break;
      if (np == 0) ifirst = i;
      ++np;
   }
   graph = new((*graphs)[graphs->GetLast()+1]) TGraph(np, &x[ifirst], &y[ifirst]); 
   graph->SetNameTitle(Form("graph%d",graphs->GetLast()+1), Form("xmin = %f xmax = %f", xmin,xmax));
   graph->SetMarkerStyle(20);
   graph->SetMarkerColor(graphs->GetLast()+1);
   
   xmin = 0.005;
   xmax = 0.010;
   np = 0;
   for (int i=0; i<gr_diff->GetN(); ++i)
   {
      if (x[i] < xmin) continue;
      if (x[i] > xmax) break;
      if (np == 0) ifirst = i;
      ++np;
   }
   graph = new((*graphs)[graphs->GetLast()+1]) TGraph(np, &x[ifirst], &y[ifirst]); 
   graph->SetNameTitle(Form("graph%d",graphs->GetLast()+1), Form("xmin = %f xmax = %f", xmin,xmax));
   graph->SetMarkerStyle(20);
   graph->SetMarkerColor(graphs->GetLast()+1);

   xmin = 0.010;
   xmax = 0.015;
   np = 0;
   for (int i=0; i<gr_diff->GetN(); ++i)
   {
      if (x[i] < xmin) continue;
      if (x[i] > xmax) break;
      if (np == 0) ifirst = i;
      ++np;
   }
   graph = new((*graphs)[graphs->GetLast()+1]) TGraph(np, &x[ifirst], &y[ifirst]); 
   graph->SetNameTitle(Form("graph%d",graphs->GetLast()+1), Form("xmin = %f xmax = %f", xmin,xmax));
   graph->SetMarkerStyle(20);
   graph->SetMarkerColor(graphs->GetLast()+1);

   xmin = 0.015;
   xmax = 0.0186;
   np = 0;
   for (int i=0; i<gr_diff->GetN(); ++i)
   {
      if (x[i] < xmin) continue;
      if (x[i] > xmax) break;
      if (np == 0) ifirst = i;
      ++np;
   }
   graph = new((*graphs)[graphs->GetLast()+1]) TGraph(np, &x[ifirst], &y[ifirst]); 
   graph->SetNameTitle(Form("graph%d",graphs->GetLast()+1), Form("xmin = %f xmax = %f", xmin,xmax));
   graph->SetMarkerStyle(20);
   graph->SetMarkerColor(graphs->GetLast()+1);

   xmin = 0.0186;
   xmax = 0.0238;
   np = 0;
   for (int i=0; i<gr_diff->GetN(); ++i)
   {
      if (x[i] < xmin) continue;
      if (x[i] > xmax) break;
      if (np == 0) ifirst = i;
      ++np;
   }
   graph = new((*graphs)[graphs->GetLast()+1]) TGraph(np, &x[ifirst], &y[ifirst]); 
   graph->SetNameTitle(Form("graph%d",graphs->GetLast()+1), Form("xmin = %f xmax = %f", xmin,xmax));
   graph->SetMarkerStyle(20);
   graph->SetMarkerColor(graphs->GetLast()+1);

   xmin = 0.0238;
   xmax = 0.0278;
   np = 0;
   for (int i=0; i<gr_diff->GetN(); ++i)
   {
      if (x[i] < xmin) continue;
      if (x[i] > xmax) break;
      if (np == 0) ifirst = i;
      ++np;
   }
   graph = new((*graphs)[graphs->GetLast()+1]) TGraph(np, &x[ifirst], &y[ifirst]); 
   graph->SetNameTitle(Form("graph%d",graphs->GetLast()+1), Form("xmin = %f xmax = %f", xmin,xmax));
   graph->SetMarkerStyle(20);
   graph->SetMarkerColor(graphs->GetLast()+1);

   xmin = 0.0278;
   xmax = 0.0328;
   np = 0;
   for (int i=0; i<gr_diff->GetN(); ++i)
   {
      if (x[i] < xmin) continue;
      if (x[i] > xmax) break;
      if (np == 0) ifirst = i;
      ++np;
   }
   graph = new((*graphs)[graphs->GetLast()+1]) TGraph(np, &x[ifirst], &y[ifirst]); 
   graph->SetNameTitle(Form("graph%d",graphs->GetLast()+1), Form("xmin = %f xmax = %f", xmin,xmax));
   graph->SetMarkerStyle(20);
   graph->SetMarkerColor(graphs->GetLast()+1);

   cout<< "plot graphs. graphs->GetLast()+1 = " << graphs->GetLast()+1 <<endl;

   for (int i=0; i<graphs->GetLast()+1; ++i)
   {
      graph = (TGraph*) graphs->At(i);
      // cout<< graph->GetTitle() <<endl;
      // for (int ipoint=0; ipoint<graph->GetN(); ++ipoint) {
      // 	 cout<< ipoint <<"\t x " << graph->GetX()[ipoint] << "\t y " << graph->GetY()[ipoint] <<endl;
      // }
      new TCanvas;
      graph->Draw("ap");
   }

   new TCanvas;
   TH1F* hframe = gPad->DrawFrame(0,0, gr_diff->GetXaxis()->GetXmax(), gr_diff->GetYaxis()->GetXmax());
   hframe->SetNameTitle("hframe", "Partial photoelectron peaks (scaled)");

   TText* text = new TText(0,0,"");

   Double_t ymax = ((TGraph*) graphs->At(0))->GetYaxis()->GetXmax();
   cout<< "ymax = " << ymax <<endl;

   for (int igraph=0; igraph<graphs->GetLast()+1; ++igraph)
   {
      graph = (TGraph*) graphs->At(igraph);
      Double_t ymax_local = ((TGraph*) graphs->At(igraph))->GetYaxis()->GetXmax();
      // Double_t scal = ymax / ymax_local;
      Int_t scale = ymax / ymax_local + 0.5;
      cout<< igraph <<"\t scale = " << scale <<endl;
      for (int ipoint=0; ipoint<graph->GetN(); ++ipoint) {
      	 graph->SetPoint(ipoint, graph->GetX()[ipoint], scale*graph->GetY()[ipoint]);
      }
      graph->Draw("pl");
      text->DrawText(graph->GetX()[0], 0.95*ymax, Form("x%d", scale));
   }
}
