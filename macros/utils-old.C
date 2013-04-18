// Andriy Zatserklyaniy <zatserkl@fnal.gov>

#include <TROOT.h>
#include <TRint.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TEnv.h>
#include <TMath.h>
#include <TF1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TEventList.h>
#include <TBrowser.h>
#include <TPaveStats.h>
#include <TCanvas.h>
#include <TList.h>
#include <TColor.h>
#include <TTimer.h>
#include <TLatex.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <cstdio>

using std::cout;                using std::endl;

#define rootalias

#ifndef utils_h
#define utils_h

// 
// to use utils.C function in the root macro declare them like
// TH1* htemp(const char* name=0, const char* title=0);
// 

void ali()
{
   TTree* pulse = (TTree*) gDirectory->Get("pulse");
   if (pulse) {
      // cout<< "Set aliases for the tree pulse" <<endl;

      pulse->SetAlias("t1", "b1_t");
      pulse->SetAlias("t2", "b2_t");

      pulse->SetAlias("c1", "b1_c1");
      pulse->SetAlias("c2", "b1_c2");
      pulse->SetAlias("c3", "b1_c3");
      pulse->SetAlias("c4", "b1_c4");
      pulse->SetAlias("c5", "b2_c1");
      pulse->SetAlias("c6", "b2_c2");
      pulse->SetAlias("c7", "b2_c3");
      pulse->SetAlias("c8", "b2_c4");
   }
}

TGraph* differentiate(const TGraph* gr);
TGraph* differentiate(const TGraph* gr)
{
   TGraph* gr_diff = (TGraph*) gr->Clone(Form("diff_%s",gr->GetName()));
   gr_diff->SetTitle(Form("diff %s",gr->GetName()));
   for (int i=0; i<gr->GetN()-1; ++i) {
      Double_t dx = i > 0? gr->GetX()[i] - gr->GetX()[i-1]: gr->GetX()[i+1] - gr->GetX()[i];
      Double_t dy = i > 0? gr->GetY()[i] - gr->GetY()[i-1]: 0;
      const Double_t eps = 1e-7;
      Double_t dydx = TMath::Abs(dx) > eps? dy/dx: dy/eps;
      gr_diff->SetPoint(i, gr->GetX()[i], dydx);
   }
   gr_diff->Set(gr->GetN()-1);
   return gr_diff;
}

// 
// modified functor from /srv/opt/root/root/tutorials/math/exampleFunctor.C
// Comments from here:
//
// in order to work with interpreter the function object must be created and lived all time for all time 
// of the TF1. In compiled mode, the function object can be passed by value (reccomended) and there 
// is also no need to specify the type of the function class. Example is as follow: 
// TF1 * f2 = new TF1("f2",MyDerivFunc(f1), xmin, xmax,0); // only for C++ compiled mode

/*
DerivativeFunctor* dfunctor = new DerivativeFunctor(fun);
TF1* dfun = new TF1("dfun", dfunctor, fun->GetXmin(), fun->GetXmax(), 0, "DerivativeFunctor");
*/
class DerivativeFunctor { 
   TF1* tf1_; 
public:
   DerivativeFunctor(TF1* tf1): tf1_(tf1) {}
   double operator() (double* x, double*) const { 
      return tf1_->Derivative(*x);
   }
};
/*
IntegralFunctor* ifunctor = new IntegralFunctor(fun);
TF1* ifun = new TF1("ifun", ifunctor, fun->GetXmin(), fun->GetXmax(), 0, "IntegralFunctor");
*/
class IntegralFunctor { 
   TF1* tf1_; 
public:
   IntegralFunctor(TF1* tf1): tf1_(tf1) {}
   double operator() (double* x, double*) const { 
      double a = tf1_->GetXmin();
      return tf1_->Integral(a, *x);
   }
};
TF1* derfun(TF1* fun);
TF1* derfun(TF1* fun) {
   DerivativeFunctor* dfunctor = new DerivativeFunctor(fun);
   TF1* dfun = new TF1(Form("d%s",fun->GetName()), dfunctor, fun->GetXmin(), fun->GetXmax(), 0, "DerivativeFunctor");
   return dfun;
}
TF1* intfun(TF1* fun);
TF1* intfun(TF1* fun) {
   IntegralFunctor* ifunctor = new IntegralFunctor(fun);
   TF1* ifun = new TF1(Form("i%s",fun->GetName()), ifunctor, fun->GetXmin(), fun->GetXmax(), 0, "IntegralFunctor");
   return ifun;
}

std::ostream& endn(std::ostream& os)
{ 
   // example of usage
   //
   // bool debug = true;
   // std::ostream& (*end)(std::ostream&);      // end is a pointer to function
   // if (debug) end = std::endl;
   // else end = endn;

   os << "\n";
   return os;
}

std::ostream& exitl(std::ostream&);
std::ostream& exitl(std::ostream&) {exit(0);}
// std::ostream& exitl(std::ostream& os) {
//    os << endl;
//    exit(0);
// }

Long64_t plot(TTree* tree, const char* varexp, const char* selection="", Option_t* option="", Long64_t nentries=1000000000, Long64_t firstentry=0);
Long64_t plot(TTree* tree, const char* varexp, const TCut& selection, Option_t* option="", Long64_t nentries=1000000000, Long64_t firstentry=0);
Long64_t project(TTree* tree, const char* hname, const char* varexp, const char* selection="", Option_t* option="", Long64_t nentries=1000000000, Long64_t firstentry=0);
Long64_t plota(TTree* tree, const char* varexp, const char* selection="", Option_t* option="", Long64_t nentries=1000000000, Long64_t firstentry=0);
Long64_t plota(TTree* tree, const char* varexp, const TCut& selection, Option_t* option="", Long64_t nentries=1000000000, Long64_t firstentry=0);
Long64_t projecta(TTree* tree, const char* hname, const char* varexp, const char* selection="", Option_t* option="", Long64_t nentries=1000000000, Long64_t firstentry=0);
bool fexists(const char* ifname);
void showbra(TTree* tree);
void history(Int_t nlines=10);
void history(const char* pattern, Int_t nlines=10);
void hi(Int_t nlines=10);
void hi(const char* pattern, Int_t nlines=10);
void historyDraw(Int_t nlines=10);
void historyDraw(const char* pattern, Int_t nlines=10);
void hiDraw(Int_t nlines=10);
void hiDraw(const char* pattern, Int_t nlines=10);
Double_t tdcu(Int_t rtd, Int_t Vpc=0);
void ref();
void log();
void logy();
void logz();
void lin();
void zoomx(Axis_t xmin=0, Axis_t xmax=0);
void zoomy(Axis_t ymin=0, Axis_t ymax=0);
void zoom(Axis_t xmin=0, Axis_t xmax=0, Axis_t ymin=0, Axis_t ymax=0);
void unzoom();          // same as unzoomx()
void unzoomx();
void unzoomy();
void unzoomxy();
void minmax(Axis_t ymin, Axis_t ymax);
TCanvas* fonts();
void colors();
Int_t Color(Int_t number);
void hcolor(TH1* h, Int_t color, Int_t fillstyle=1001);
void tcolor(TTree* tree, Int_t color, Int_t hatch=3001);
void tmarker(TTree* tree, Int_t marker);
void bra(TTree* tree, const char* pattern="");
void fwlite();
void tb(const char* fname=0);
const char* substitute(const std::string& str0, const std::string what, const std::string by);
int nlen(const std::string str);
Int_t canvas_max_number(const char* pattern="c1_n");
const char* nextname(const char* base="h");
const char* nextcan(const char* base="can");
const char* defcan();
// TCanvas* newcan(Int_t width=525, Int_t height=375, const char* name="", const char* title="");
// TCanvas* newcan(Int_t width=522, Int_t height=354, const char* name="", const char* title="");
TCanvas* newcan(Int_t width=0, Int_t height=0, const char* name="", const char* title="");
TString NextName(TString base="h");
Float_t qsum(Float_t e1, Float_t e2);
void nostatover();
void gri();
void nogri();
void nbi(Int_t nx=0, Int_t ny=0);     // gEnv->SetValue("Hist.Binning.1D.x", n)
void pro(const char* prompt="// ");
void rex();
TPaveStats* movestat(Double_t x1NDC=0, Double_t y1NDC=0, Double_t xlen=0, Double_t ylen=0);
void left();
void right();
void leftgaus();
void rightgaus();
void lg();
void rg();
Double_t ngaus(TH1* h=0);
const char* pwd();
TFile* cd(const char* filename="");
TFile* cd(Int_t fnumber);                                       // NB: no default value, it's in string version
TGraphErrors* getemp(const char* name=0, const char* title=0);
TGraphErrors* getempClone(const char* name=0, const char* title=0);
TCanvas* gcan();
void resize(Int_t width=0, Int_t height=0, Bool_t internal=kFALSE);
TGraph* gtemp(const char* name=0, const char* title=0);
TGraph* gtempClone(const char* name=0, const char* title=0);
TGraph* grtemp(const char* name=0, const char* title=0);
TGraph* grtempClone(const char* name=0, const char* title=0);
TMultiGraph* mgtemp(const char* name=0, const char* title=0);
TPaveStats* stemp(TH1* h=0);     // return pointer to the stat box
TGraphErrors* gsubf(const TGraph* gr=0, const TF1* f=0, const char* name="", const char* title="");
Int_t axisNdivisions(Int_t Ndivisions=0);
Int_t ndi(Int_t Ndivisions=0);
void clone();
TPaveText* statstemp(TVirtualPad* can=0);
// TPaveText* titletemp(TVirtualPad* can=0);
// TPaveText* titletemp(const char* title, TVirtualPad* can=0);
TH1* hlast(const char* name=0, const char* title=0);
TH1* hfirst(const char* name=0, const char* title=0);
TH1* htemp(const char* name=0, const char* title=0);
TH1* htempClone(const char* name=0, const char* title=0);
TH2* htemp2(const char* name=0, const char* title=0);
TH2* htemp2Clone(const char* name=0, const char* title=0);
TF1* ftemp(const char* name=0, const char* title=0);
TF1* ftempClone(const char* name=0, const char* title=0);
TGraph* tgraph(TTree* tree, const char* name="", const char* title="");
void drawpan(TH1* h=0);
void fitpan(TH1* h=0);
void InitPredefinedFunctions();
void fitgaus(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void fitg(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void fitgl(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void fitgr(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void fitexpo(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void fite(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void fitel(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void fiter(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void fitpol(Double_t xmin=0, Double_t xmax=0, Int_t power=1, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void fitp(Double_t xmin=0, Double_t xmax=0, Int_t power=1, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
Double_t gaus_ampl();
Double_t gaus_mean();
Double_t gaus_sigma();
Double_t gaus_eampl();
Double_t gaus_emean();
Double_t gaus_esigma();
void pargaus(Double_t& a, Double_t& mean, Double_t& sigma, const char* hname=0);
void npe(TH1* h, Double_t* p_Npe=0, Double_t* p_eNpe=0);
void npe(const char* hname=0, Double_t* p_Npe=0, Double_t* p_eNpe=0);
void pic(TString pathname="", TString suffix="");
void eps(TString pathname="", TString suffix="");
void png(TString pathname="", TString suffix="", bool check_fexists=true);
void pngall();
Double_t intgaus(TH1* h=0);
TH1* hname(const char* name=0);
void print_htemp(TH1* h=0);
void hprint(TH1* h=0);
void tax(const char* format="");
void nul(Axis_t min=0);
void ct();
void ctadd();
void ctlist();
const char* addtit(const char* title, TCanvas* can=0);
TChain* chlist(const char* tree_name, const char* filelist);
TTree* mkread(const char* fname, bool addFriend=false);
Stat_t hint(char* hname, Float_t x1=0, Float_t x2=0);
Stat_t hint(TH1F* h, Float_t x1=0, Float_t x2=0);
void cloop(const char* command);

#endif  // define utils_h

#ifndef utils_cxx
#define utils_cxx

// 1) TTree::Draw gives wrong result for plotting of single varexp when selection includes arrays
// 2) plot works correctly in this case (single varexp when selection includes arrays)
// 3) plot gives wrong result for plotting of array variable with selection with arrays
// 4) TTree::Draw correctly plots array varexp with selection with arrays (the same true for plota)
//
// functions plot and project: for single (not array) varexp and arrays in selection
// TTree::Draw and TTree::Project give wrong result in this case
// NB: do not apply for array varexp: gives wrong result. Use TTree::Draw and TTree::Project or plota and projecta.
//
Long64_t plot(TTree* tree, const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   TEventList* elist_curr = tree->GetEventList();
   // create list of events for given selection
   tree->Draw(">>elist_temp", selection, option, nentries, firstentry);
   TEventList* elist_temp = (TEventList*) gDirectory->Get("elist_temp");
   tree->SetEventList(elist_temp);
   // plot varexp without selection to provide exactly one entry per event
   Long64_t nselected = tree->Draw(varexp, "", option, nentries, firstentry);
   tree->SetEventList(elist_curr);
   elist_temp->SetDirectory(0);
   delete elist_temp;
   return nselected;
}
Long64_t plot(TTree* tree, const char* varexp, const TCut& selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   TEventList* elist_curr = tree->GetEventList();
   tree->Draw(">>elist_temp", selection, option, nentries, firstentry);
   TEventList* elist_temp = (TEventList*) gDirectory->Get("elist_temp");
   tree->SetEventList(elist_temp);
   Long64_t nselected = tree->Draw(varexp, "", option, nentries, firstentry);
   tree->SetEventList(elist_curr);
   elist_temp->SetDirectory(0);
   delete elist_temp;
   return nselected;
}
Long64_t project(TTree* tree, const char* hname, const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   TEventList* elist_curr = tree->GetEventList();
   tree->Draw(">>elist_temp", selection, option, nentries, firstentry);
   TEventList* elist_temp = (TEventList*) gDirectory->Get("elist_temp");
   tree->SetEventList(elist_temp);
   Long64_t nselected = tree->Project(hname, varexp, "", option, nentries, firstentry);
   tree->SetEventList(elist_curr);
   elist_temp->SetDirectory(0);
   delete elist_temp;
   return nselected;
}
//-- NB: there is no TTree::Project with TCut selection
//
// functions plota and projecta: for array varexp and arrays in selection. Give the same result as TTree::Draw and TTree::Project
// NB: do not apply for single (not array) varexp: gives wrong result. Use plot and project.
//
Long64_t plota(TTree* tree, const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   // array version of plot: apply for array varexp. Gives the same result as just TTree::Draw
   TEventList* elist_curr = tree->GetEventList();
   // create list of events for given selection
   tree->Draw(">>elist_temp", selection, option, nentries, firstentry);
   TEventList* elist_temp = (TEventList*) gDirectory->Get("elist_temp");
   tree->SetEventList(elist_temp);
   //-- substitute "" by selection in 'Long64_t nselected = tree->Draw(varexp, "", option, nentries, firstentry);'
   Long64_t nselected = tree->Draw(varexp, selection, option, nentries, firstentry);
   tree->SetEventList(elist_curr);
   elist_temp->SetDirectory(0);
   delete elist_temp;
   return nselected;
}
Long64_t plota(TTree* tree, const char* varexp, const TCut& selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   // array version of plot: apply for array varexp. Gives the same result as just TTree::Draw
   TEventList* elist_curr = tree->GetEventList();
   tree->Draw(">>elist_temp", selection, option, nentries, firstentry);
   TEventList* elist_temp = (TEventList*) gDirectory->Get("elist_temp");
   tree->SetEventList(elist_temp);
   //-- substitute "" by selection in 'Long64_t nselected = tree->Draw(varexp, "", option, nentries, firstentry);'
   Long64_t nselected = tree->Draw(varexp, selection, option, nentries, firstentry);
   tree->SetEventList(elist_curr);
   elist_temp->SetDirectory(0);
   delete elist_temp;
   return nselected;
}
Long64_t projecta(TTree* tree, const char* hname, const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   // array version of project: apply for array varexp. Gives the same result as just TTree::Project
   TEventList* elist_curr = tree->GetEventList();
   tree->Draw(">>elist_temp", selection, option, nentries, firstentry);
   TEventList* elist_temp = (TEventList*) gDirectory->Get("elist_temp");
   tree->SetEventList(elist_temp);
   //-- substitute "" by selection in 'Long64_t nselected = tree->Project(hname, varexp, "", option, nentries, firstentry);'
   Long64_t nselected = tree->Project(hname, varexp, selection, option, nentries, firstentry);
   tree->SetEventList(elist_curr);
   elist_temp->SetDirectory(0);
   delete elist_temp;
   return nselected;
}

bool fexists(const char* ifname) {
  bool fexists = false;
  std::ifstream ifile(ifname, std::ifstream::in|std::ifstream::binary);
  if (ifile.is_open()) {
     ifile.close();
     fexists = true;
  }
  return fexists;
}

void showbra(TTree* tree) {
   for (int i=0; i<tree->GetListOfBranches()->GetEntries(); ++i) {
      TBranch* branch = (TBranch*) tree->GetListOfBranches()->At(i);
      void* address = branch->GetAddress();
      cout<< "branch " << branch->GetName() << " address: " << address <<endl;
   }
}

void history(Int_t nlines)
{
   gSystem->Exec(Form("tail -%d root_hist",nlines));
}
void history(const char* pattern, Int_t nlines)
{
   if (nlines > 0) gSystem->Exec(Form("fgrep %s root_hist | tail -%d",pattern,nlines));
   else gSystem->Exec(Form("fgrep %s root_hist",pattern));
}
void hi(Int_t nlines) {history(nlines);}
void hi(const char* pattern, Int_t nlines) {history(pattern,nlines);}

void historyDraw(Int_t nlines)
{
   gSystem->Exec(Form("fgrep Draw root_hist | tail -%d",nlines));
}
void historyDraw(const char* pattern, Int_t nlines)
{
   if (nlines > 0) gSystem->Exec(Form("fgrep Draw root_hist | fgrep %s | tail -%d",pattern,nlines));
   else gSystem->Exec(Form("fgrep %s root_hist",pattern));
}
void hiDraw(Int_t nlines) {historyDraw(nlines);}
void hiDraw(const char* pattern, Int_t nlines) {historyDraw(pattern,nlines);}

Double_t tdcu(Int_t rtd, Int_t Vpc) {
   static Int_t Vportcard;
   if (Vpc > 0) Vportcard = Vpc;
   if (Vportcard == 0) {
      cout<< "Usage: tdcu(Int_t rtd, Int_t Vpc)" <<endl;
      return 0;
   }
   Double_t T = (2/3850.e-6) * (1. - Double_t(rtd)/Double_t(Vportcard));
   return T;
}

void ref() {
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) return;
   gPad->Modified();
   gPad->Update();
}

void log() {
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) return;
   gPad->SetLogy(!gPad->GetLogy());
   gPad->Update();
}

void logy() {log();}

void logz() {
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) return;
   gPad->SetLogz(!gPad->GetLogz());
   gPad->Update();
}

void logx() {
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) return;
   gPad->SetLogx(!gPad->GetLogx());
   gPad->Update();
}

void lin() {
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) return;
   gPad->SetLogy(0);
   gPad->Update();
}

void zoomx(Axis_t xmin, Axis_t xmax) {
   TAxis *xaxis = 0;
   if (htemp()) xaxis = htemp()->GetXaxis();
   else if (gtemp()) xaxis = gtemp()->GetXaxis();
   else if (ftemp()) xaxis = ftemp()->GetXaxis();
   if (xaxis) {
      if (xmax == 0 and xmax < xmin) xmax = xaxis->GetXmax();
      xaxis->SetRangeUser(xmin,xmax);
      gPad->Modified();
      gPad->Update();
   }
}

void zoomy(Axis_t ymin, Axis_t ymax) {
   TAxis *yaxis = 0;
   if (htemp()) yaxis = htemp()->GetYaxis();
   else if (gtemp()) yaxis = gtemp()->GetYaxis();
   else if (ftemp()) yaxis = ftemp()->GetYaxis();
   if (yaxis) {
      if (ymax == 0 and ymax < ymin) ymax = yaxis->GetXmax();
      yaxis->SetRangeUser(ymin,ymax);
      gPad->Modified();
      gPad->Update();
   }
}

void zoom(Axis_t xmin, Axis_t xmax, Axis_t ymin, Axis_t ymax) {
   TAxis *xaxis = 0;
   TAxis *yaxis = 0;
   if (htemp()) {
      xaxis = htemp()->GetXaxis();
      yaxis = htemp()->GetYaxis();
   }
   else if (gtemp()) {
      xaxis = gtemp()->GetXaxis();
      yaxis = gtemp()->GetYaxis();
   }
   else if (ftemp()) {
      xaxis = ftemp()->GetXaxis();
      yaxis = ftemp()->GetYaxis();
   }
   if (xmax > xmin and xaxis) {
      if (xmax == 0 and xmax < xmin) xmax = xaxis->GetXmax();
      xaxis->SetRangeUser(xmin,xmax);
      gPad->Modified();
      gPad->Update();
   }
   if (ymax > ymin and yaxis) {
      if (ymax == 0 and ymax < ymin) ymax = yaxis->GetXmax();
      yaxis->SetRangeUser(ymin,ymax);
      gPad->Modified();
      gPad->Update();
   }
}

void unzoom() {                                 // same as unzoomx()
   TAxis *xaxis = 0;
   if (htemp()) xaxis = htemp()->GetXaxis();
   else if (gtemp()) xaxis = gtemp()->GetXaxis();
   else if (ftemp()) xaxis = ftemp()->GetXaxis();
   if (xaxis) {
      xaxis->UnZoom();
      gPad->Modified();
      gPad->Update();
   }
}
void unzoomx() {
   TAxis *xaxis = 0;
   if (htemp()) xaxis = htemp()->GetXaxis();
   else if (gtemp()) xaxis = gtemp()->GetXaxis();
   else if (ftemp()) xaxis = ftemp()->GetXaxis();
   if (xaxis) {
      xaxis->UnZoom();
      gPad->Modified();
      gPad->Update();
   }
}
void unzoomy() {
   TAxis *yaxis = 0;
   if (htemp()) yaxis = htemp()->GetYaxis();
   else if (gtemp()) yaxis = gtemp()->GetYaxis();
   else if (ftemp()) yaxis = ftemp()->GetYaxis();
   if (yaxis) {
      yaxis->UnZoom();
      gPad->Modified();
      gPad->Update();
   }
}
void unzoomxy() {
   TAxis *xaxis = 0;
   TAxis *yaxis = 0;
   if (htemp()) {
      xaxis = htemp()->GetXaxis();
      yaxis = htemp()->GetYaxis();
   }
   else if (gtemp()) {
      xaxis = gtemp()->GetXaxis();
      yaxis = gtemp()->GetYaxis();
   }
   else if (ftemp()) {
      xaxis = ftemp()->GetXaxis();
      yaxis = ftemp()->GetYaxis();
   }
   if (xaxis && yaxis) {
      xaxis->UnZoom();
      yaxis->UnZoom();
      gPad->Modified();
      gPad->Update();
   }
}

void minmax(Axis_t ymin, Axis_t ymax) {
   if (htemp()) {
      htemp()->SetMinimum(ymin);
      htemp()->SetMaximum(ymax);
   }
   else if (gtemp()) {
      gtemp()->SetMinimum(ymin);
      gtemp()->SetMaximum(ymax);
   }
   else if (ftemp()) {
      ftemp()->SetMinimum(ymin);
      ftemp()->SetMaximum(ymax);
   }
   gPad->Modified();
   gPad->Update();
}

void fontsdrawtext(double x, double y, int f, const char *s)
{
   TLatex* t = new TLatex(x,y,Form("#font[41]{%d :} %s",f,s));
   t->SetTextFont(f);
   t->SetTextAlign(12);
   t->SetTextSize(0.048);
   t->Draw();
}
TCanvas* fonts()
{
   TCanvas* Tf = new TCanvas("Tf", "Tf",0,0,800,400);
   Tf->Range(0,0,1,1);
   Tf->SetBorderSize(2);
   Tf->SetFrameFillColor(0);
                                                                                
   double y = 0.95;
   for (int f = 12; f<=152; f+=10) {
      if (f!=142) fontsdrawtext(0.02,y, f,"ABCDEFGHIJKLMNOPQRSTUVWXYZ abcdefghijklmnopqrstuvwxyz 0123456789 @#$");
      else fontsdrawtext(0.02,y, f,"ABCD efgh 01234 @#$");
      y -= 0.065;
   }
   return Tf;
}

void colors()
{
   TCanvas *colors = new TCanvas("colors","Fill Area colors",0,0,500,200);
   colors->DrawColorTable();
}

Int_t Color(Int_t number)
{
   if (number == 0) return 0;
   if (number == 1) return 602;
   if (number == 2) return 2;
   if (number == 3) return 30;
   if (number == 4) return 4;
   if (number == 5) return 41;
   if (number == 6) return 6;
   if (number == 7) return 38;
   if (number == 8) return 8;
   if (number == 9) return 9;
   if (number == 10) return 29;
   if (number == 11) return 12;
   if (number == 12) return 46;
   if (number == 13) return 3;
   if (number == 14) return 39;
   if (number == 15) return 5;
   if (number == 16) return 48;
   if (number == 17) return 36;
   if (number == 18) return 29;
   if (number == 19) return 33;
   if (number == 20) return 1;
   if (number == 21) return 37;
   if (number == 22) return 28;
   if (number == 23) return 32;
   if (number == 24) return 40;
   if (number == 25) return 44;
   if (number == 26) return 47;
   if (number == 27) return 7;
   if (number == 28) return 34;
   if (number == 29) return 49;
   if (number == 30) return 11;
   return number;
}

void hcolor(TH1* h, Int_t color, Int_t fillstyle) {
   h->SetFillColor(color);
   h->SetLineColor(color);
   h->SetMarkerColor(color);
   h->SetFillStyle(fillstyle);
}

void hline(Int_t color) {
   TH1* h = hlast();
   if (h) {
      h->SetLineColor(color);
      h->SetMarkerColor(color);
      gPad->Modified();
      gPad->Update();
   }
}

void hfill(Int_t color) {
   TH1* h = hlast();
   if (h) {
      h->SetFillColor(color);
      gPad->Modified();
      gPad->Update();
   }
}

void tcolor(TTree* tree, Int_t color, Int_t hatch)
{
   tree->SetMarkerColor(color);
   tree->SetLineColor(color);
   tree->SetFillColor(color);
   tree->SetFillStyle(hatch);
   if (color == 0 || color == 1) {
      tree->SetMarkerColor(1);
      tree->SetLineColor(1);
      tree->SetFillColor(0);
   }
}

void tmarker(TTree* tree, Int_t marker)
{
   tree->SetMarkerStyle(marker);
}

// prints tree's list of branches
void bra(TTree* tree, const char* pattern)
{
   TIter next(tree->GetListOfBranches());
   TBranch* b;
   cout<< "----- TTree " << tree->GetName() << " -----" <<endl;
   while ((b = (TBranch*) next())) {
      if (pattern && !strstr(b->GetName(), pattern)) continue;
      cout<< b->GetName() <<endl;
   }
}

void fwlite() {
   gSystem->Load("libFWCoreFWLite.so");
   //-- AutoLibraryLoader::enable();
}

// TF1* erf = new TF1("erf", "TMath::Erf(x)", -10, 10);

// void hmove(TDirectory* from, TDirectory* to=0, bool verb=true)
// {
//    if (!to) {
//       TDirectory* dir = gDirectory;
//       gROOT->cd();
//       to = gDirectory;
//       dir->cd();
//    }
//    from->cd();
//    TList* lkey = from->GetListOfKeys();
//    TIter next_key(lkey);
//    TKey* key = 0;
//    Int_t nmoved = 0;
//    while ((key = (TKey*) next_key())) {
//       TObject* obj = from->Get(key->GetName());
//       if (obj->InheritsFrom("TH1")) {
//       //TH1* h = (TH1*) obj->Clone();
//       TH1* h = (TH1*) obj;
//       h->SetDirectory(to);
//       ++nmoved;
//       }
//    }
//    if (verb) cout<< "The number of histos moved is " << nmoved <<endl;
// }

void tb(const char* fname) {
   if (fname) TFile::Open(fname);
   new TBrowser("TBrowser", "TBrowser", 600, 0, 670, 500);
}

// void substitute(string& str, const string what, const string by)
// {
//    for (string::size_type loc = str.find(what,0); loc != string::npos; ) {
//       str.replace(loc, what.length(), by);
//       loc = str.find(what,0);
//    }
// }
// std::string substitute(const std::string& str0, const std::string what, const std::string by)
// {
//    std::string str = str0;
//    for (std::string::size_type loc = str.find(what,0); loc != std::string::npos; ) {
//       str.replace(loc, what.length(), by);
//       loc = str.find(what,0);
//    }
//    return str;
// }
const char* substitute(const std::string& str0, const std::string what, const std::string by)
{
   if (what == by) return str0.c_str();
   std::string str = str0;
   for (std::string::size_type loc = str.find(what,0); loc != std::string::npos; ) {
      str.replace(loc, what.length(), by);
      loc = str.find(what,0);
   }
   //cout<< "str0 = " << str0 << " str = " << str <<endl;
   return Form("%s",str.c_str());       // put in Form circular buffer
}
int nlen(const std::string str)
{       // string length without leading spaces
   unsigned i = 0;
   for (; i<str.length(); ++i) {
      if (isspace(str[i])) continue;
      else break;
   }
   return strlen(str.c_str()) - i;
}

Int_t canvas_max_number(const char* pattern)
{
   TSeqCollection* list = gROOT->GetListOfCanvases();

   Int_t number_max = 0;
   for (int i=0; i<list->GetEntries(); ++i)
   {
      Int_t number = 0;
      std::string format_string = std::string(pattern) + std::string("%d");
      std::sscanf(list->At(i)->GetName(), format_string.c_str(), &number);
      if (number > number_max) number_max = number;
   }
   return number_max;
}

const char* nextname(const char* base) {
   // if name 'base' already exists increaments it like h_1, h_2, etc.
   // returns pointer to circular buffer owned by ROOT function Form
   bool found = gDirectory->Get(base);
   if (!found) return base;
   Int_t icycle = 0;
   while (found) {
      std::stringstream ss;
      ss.str("");
      ss << base << "_" << ++icycle;
      //cout<< "         current hname: " << ss.str() <<endl;
      found = gDirectory->Get(ss.str().c_str());
   }
   //cout<< "new hname: " << Form("%s_%d",base,icycle) <<endl;
   return Form("%s_%d",base,icycle);
}

const char* nextcan(const char* base) {
   const char* cname = base;
   // if canvas 'cname' already exists increaments it like can_1, can_2, etc.
   // returns pointer to circular buffer owned by ROOT function Form
   int icycle = 0;
   while (gROOT->GetListOfCanvases()->FindObject(cname)) {
      cname = Form("%s_%d", base, ++icycle);
   }
   //cout<< "cname = " << cname <<endl;
   return cname;
}

const char* defcan() {
   TList *lc = (TList*)gROOT->GetListOfCanvases();
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) return Form("c1");

   const char *defcanvas = gROOT->GetDefCanvasName();
   char *cdef = 0;
   if (lc->FindObject(defcanvas)) {
      // taken from my correction to TCanvas::TCanvas(Bool_t build)
      //-- original: cdef = StrDup(Form("%s_n%d",defcanvas,lc->GetSize()+1));
      //-- Andriy Zatserklyaniy, Jan 26, 2012
      Int_t n = lc->GetSize()+1;
      while (lc->FindObject(Form("%s_n%d",defcanvas,n))) n++;
      cdef = Form("%s_n%d",defcanvas,n);
   }
   return cdef;
}

TCanvas* newcan(Int_t width, Int_t height, const char* name, const char* title)
{
   TString name_str = name;
   // if (!name || strlen(name)==0) name_str = Form("c_%dx%d",width,height);
   //-- if (!name || !*name) name_str = Form("c_1");
   if (!name || !*name) name_str = Form("c");
   name_str = nextcan(name_str.Data());
   TString title_str = title;
   if (!title || !*title) title_str = name_str;
   if (width == 0) width = gStyle->GetCanvasDefW();
   if (height == 0) height = gStyle->GetCanvasDefH();
   TCanvas * can = new TCanvas(name_str.Data(), title_str.Data(), width, height);
   //-- can->SetWindowSize(width + (width - can->GetWw()), height + (height - can->GetWh()));
   return can;
}

TString NextName(TString base)
{
   // if name 'base' already exists increaments it like h_1, h_2, etc.
   // returns TString which holds buffer with name
   TString out = base;
   Int_t i = 0;
   bool found = gDirectory->Get(out.Data());
   while (found) {
      std::stringstream ss;
      ss.str("");
      ss << base << "_" << ++i;
      //cout<< "         current hname: " << out <<endl;
      out = ss.str().c_str();
      TObject* obj = gDirectory->Get(out.Data());
      //cout<< "obj = " << obj <<endl;
      //found = gDirectory->Get(out.Data());
      found = (obj != 0);
   }
   //cout<< "new hname: " << out <<endl;
   return out;
}

// TObjArray* ahist(const char* ifname, const char* startpat=0, const char* pattern=0, bool verb=true);
// TObjArray* ahist(const char* ifname, const char* startpat, const char* pattern, bool verb)
// {
//    TDirectory* dir = gDirectory;
//    TFile* ifile = TFile::Open(ifname);
//    if (!ifile) {
//       cout<< "*ERROR ahist: file not found: " << ifname <<endl;
//       return 0;
//    }
//    TObjArray* ar = ahist(ifile,startpat,pattern,verb);
//    ifile->Close();
//    dir->cd();
//    return ar;
// }
// TObjArray* ahist(TDirectory* ifile, const char* startpat=0, const char* pattern=0, bool verb=true);
// TObjArray* ahist(TDirectory* ifile, const char* startpat, const char* pattern, bool verb)
// {
//    TObjArray* ahist = new TObjArray();
//    TList* lkey = ifile->GetListOfKeys();
//    TIter next_key(lkey);
//    TKey* key = 0;
//    while (key = (TKey*) next_key()) {
//       TObject* obj = ifile->Get(key->GetName());
//       if (obj->InheritsFrom("TH1")) {
//       if (startpat && (strstr(obj->GetName(),startpat) != obj->GetName())) continue;
//       if (pattern && !strstr(obj->GetName(),pattern)) continue;
//       // cout<< ahist.GetEntries() << "\t " << obj->GetName() <<endl;
//       TH1* h = (TH1*) obj->Clone();
//       h->SetDirectory(0);
//       ahist->Add(h);
//       }
//    }
//    if (verb) cout<< "ahist: The number of histos in the array is " << ahist->GetEntries() <<endl;
//    return ahist;
// }

Float_t qsum(Float_t e1, Float_t e2) {return TMath::Sqrt(e1*e1 + e2*e2);}

void statover();
void statover() {
   //-- used by lq
   // underflows and overflows are used by the Fill functions in the computation of statistics (mean value, RMS). 
   // By default, underflows or overflows are not used.
   TH1::StatOverflows(kTRUE);
}
void nostatover() {
   //-- default in ROOT
   // underflows and overflows are used by the Fill functions in the computation of statistics (mean value, RMS). 
   // By default, underflows or overflows are not used.
   TH1::StatOverflows(kFALSE);
}
void gri() {
   gStyle->SetPadGridY(1);
   gStyle->SetPadGridX(1);
}
void nogri() {
   gStyle->SetPadGridY(0);
   gStyle->SetPadGridX(0);
}
// Int_t nbi(Int_t nx, Int_t ny) {
//    if (nx == 0) return gEnv->GetValue("Hist.Binning.1D.x", 0);
//    gEnv->SetValue("Hist.Binning.1D.x", nx);
//    if (ny > 0) gEnv->SetValue("Hist.Binning.2D.y", ny);
//    return nx;
// }
void nbi(Int_t nx, Int_t ny) {
   if (nx == 0) cout<< "nx = " << gEnv->GetValue("Hist.Binning.1D.x", 0) << " ny = " << gEnv->GetValue("Hist.Binning.2D.y", 0) <<endl;
   if (nx > 0) gEnv->SetValue("Hist.Binning.1D.x", nx);
   if (ny > 0) gEnv->SetValue("Hist.Binning.2D.y", ny);
}
void pro(const char* prompt) {((TRint*) gROOT->GetApplication())->SetPrompt(prompt);}
void rex() {gPad->RedrawAxis("g");}

//
// movestat with leftgaus and rightgaus
//
TPaveStats* movestat(Double_t x1NDC, Double_t y1NDC, Double_t xlen, Double_t ylen)
{
   // Good settings to move box at the left upper corner is movestat(0.15,0.5, 0.3,0.4)
   // i.e. y1NDC + ylen = 0.9

   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   TPaveStats* stat = 0;
   TIter next(gPad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) {
      if (obj->IsA()->InheritsFrom("TH1")) {
         TH1* ph = (TH1*) obj;
         stat = (TPaveStats*) ph->GetListOfFunctions()->FindObject("stats");
         break;
      }
      if (obj->IsA()->InheritsFrom("TGraph")) {
         TGraph* pg = (TGraph*) obj;
         stat = (TPaveStats*) pg->GetListOfFunctions()->FindObject("stats");
         break;
      }
      obj=0;
   }
   if (!stat) {
      cout<< "movestat: No stat box found. Try \ngPad->Modified();\ngPad->Update();" <<endl;
      return 0;
   }

   if (x1NDC == 0 && y1NDC == 0 && xlen == 0 && ylen == 0) {
      cout<< "Current: stat->GetX1NDC() = " << stat->GetX1NDC() << " stat->GetY1NDC() = " << stat->GetY1NDC() << " xlen " << stat->GetX2NDC() - stat->GetX1NDC() << " ylen " << stat->GetY2NDC() - stat->GetY1NDC() <<endl;
      return stat;
   }

   if (x1NDC == 0) x1NDC = stat->GetX1NDC();
   if (y1NDC == 0) y1NDC = stat->GetY1NDC();
   if (xlen == 0) xlen = stat->GetX2NDC() - stat->GetX1NDC();
   if (ylen == 0) ylen = stat->GetY2NDC() - stat->GetY1NDC();

   stat->SetX1NDC(x1NDC);
   stat->SetX2NDC(x1NDC + xlen);
   stat->SetY1NDC(y1NDC);
   stat->SetY2NDC(y1NDC + ylen);

   gPad->Modified();
   gPad->GetMother()->Update();

   return stat;
}

void left()               { movestat(0.15,0.5, 0.35,0.4); }
void right()      { movestat(0.64,0.6, 0.35,0.4); }
void leftgaus()  { movestat(0.15,0.5, 0.35,0.4); }
void rightgaus() { movestat(0.64,0.5, 0.35,0.4); }

// short names for rightgaus and leftgaus
void lg() {leftgaus();}
void rg() {rightgaus();}

Double_t ngaus(TH1* h)
{
   if (!h) h = htemp();
   if (!h) return 0;

   TF1* f = h->GetFunction("gaus");
   if (!f) return 0;

   // access histogram to obtain a bin width (assumed const)
   Double_t w = h->GetBinWidth(1);

   Double_t A = f->GetParameter(0);
   //Double_t mean = f->GetParameter(1);
   Double_t sigma = f->GetParameter(2);

   Double_t sum = TMath::Sqrt(2.*TMath::Pi()) * A * sigma / w;

   cout<< "The number of events in the gaussian fit is " << sum <<endl;

   return sum;
}

const char* pwd()
{
   return Form(gDirectory->GetName());
}

TFile* cd(const char* filename)
{
   if (gROOT->GetListOfFiles()->GetEntries() == 0) {
      cout<< "No file open" <<endl;
      return 0;
   }

   TFile* file = (TFile*) gDirectory;

   if (filename and *filename)
   {
      file = (TFile*) gROOT->FindObject(filename);
      if (file) {
         file->cd();
         cout<< "--> " << "     " << gDirectory->GetName() <<endl;
         return file;
      }
      else {
         cout<< "File not found: " << filename <<endl;
         cout<< "Current dir is still " << gDirectory->GetName() <<endl;
         return 0;
      }
   }
   else {
      // print names of all opened files
      if (gROOT->GetListOfFiles()->GetEntries() == 0) return file;
      TIter next(gROOT->GetListOfFiles());
      int i = 0;
      for (TObject* obj=0; (obj = next());) {
         //-- if (strcmp(obj->GetName(), gDirectory->GetName())==0) cout<< "--> " << std::setw(2) << i << " " << obj->GetName() <<endl;
         if ((TFile*) obj == file) cout<< "--> " << std::setw(2) << i << " " << obj->GetName() <<endl;
         else cout<< std::setw(6) << i << " " << obj->GetName() <<endl;
         ++i;
      }

      cout<< "Enter file/dir number (<CR>=Quit): ";
      std::string snumber;
      std::getline(cin,snumber);

      // read to snumber the first non-blank string (if any)
      std::istringstream ss(snumber);
      snumber.clear();                                                  // clear to reuse
      ss >> snumber;
      if (snumber.size() == 0) return file;

      ss.clear();                                                                       // clear EOF
      ss.seekg(0);                                                              // back to the beginning of the stream

      int number = -1;
      ss >> number;
      if (number >= 0 && number < gROOT->GetListOfFiles()->GetEntries()) {
         file = (TFile*) gROOT->GetListOfFiles()->At(number);
         file->cd();
         cout<< "--> " << setw(3) << " " << gDirectory->GetName() <<endl;
      }
      else {
         cout<< "Out of range. Current dir is still " << gDirectory->GetName() <<endl;
      }
      return file;
   }

   return (TFile*) gDirectory;
}

TFile* cd(Int_t fnumber)
{
   if (gROOT->GetListOfFiles()->GetEntries() == 0) {
      cout<< "No file open" <<endl;
      return 0;
   }

   TFile* file = (TFile*) gDirectory;

   // print names of all opened files
   if (gROOT->GetListOfFiles()->GetEntries() == 0) return file;
   TIter next(gROOT->GetListOfFiles());
   int i = 0;
   for (TObject* obj=0; (obj = next());) {
      //-- if (strcmp(obj->GetName(), gDirectory->GetName())==0) cout<< "--> " << std::setw(2) << i << " " << obj->GetName() <<endl;
      if ((TFile*) obj == file) cout<< "--> " << std::setw(2) << i << " " << obj->GetName() <<endl;
      else cout<< std::setw(6) << i << " " << obj->GetName() <<endl;
      ++i;
   }

   if (fnumber >= 0 && fnumber < gROOT->GetListOfFiles()->GetEntries()) {
      file = (TFile*) gROOT->GetListOfFiles()->At(fnumber);
      file->cd();
      cout<< "--> " << setw(3) << " " << gDirectory->GetName() <<endl;
   }
   else {
      cout<< "Out of range. Current dir is still " << gDirectory->GetName() <<endl;
   }
   return file;
}

// gtemp(), htemp() and ftemp() functions

TCanvas* gcan()
{
   // can be used for e.g. gcan()->SetWindowSize(300, 300);      // NB: does not work: gPad->SetWindowSize(300, 300);
   if (gROOT->GetListOfCanvases()->GetEntries() > 0) return (TCanvas*) gROOT->GetListOfCanvases()->FindObject(gPad);
   else return 0;
}

void resize(Int_t width, Int_t height, Bool_t internal)
{
   if (width == 0) width = gStyle->GetCanvasDefW();
   if (height == 0) height = gStyle->GetCanvasDefH();
   // if (gROOT->GetListOfCanvases()->GetEntries() > 0) static_cast<TCanvas*>(gROOT->GetListOfCanvases()->FindObject(gPad))->SetWindowSize(width,height);
   if (gROOT->GetListOfCanvases()->GetEntries() > 0) {
      TCanvas* can = static_cast<TCanvas*>(gROOT->GetListOfCanvases()->FindObject(gPad));
      if (internal) can->SetWindowSize(width + (width - can->GetWw()), height + (height - can->GetWh()));
      else can->SetWindowSize(width,height);
   }
}

TGraphErrors* getemp(const char* name, const char* title)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   TIter next(gPad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) if (obj->IsA()->InheritsFrom("TGraphErrors")) break; else obj=0;
   if (!obj) return 0;
   TGraphErrors* gr = (TGraphErrors*) obj;
   if (gr && name && *name) gr->SetName(name);
   if (gr && title && *title) gr->SetTitle(title);
   return gr;
}
TGraphErrors* getempClone(const char* name, const char* title)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   TIter next(gPad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) if (obj->IsA()->InheritsFrom("TGraphErrors")) break; else obj=0;
   if (!obj) return 0;
   TGraphErrors* gr = (TGraphErrors*) obj->Clone();
   if (gr && name && *name) gr->SetName(name);
   if (gr && title && *title) gr->SetTitle(title);
   return gr;
}

TGraph* gtemp(const char* name, const char* title)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   TIter next(gPad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) if (obj->IsA()->InheritsFrom("TGraph")) break; else obj=0;
   if (!obj) return 0;
   TGraph* gr = (TGraph*) obj;
   if (gr && name && *name) gr->SetName(name);
   if (gr && title && *title) gr->SetTitle(title);
   return gr;
}
TGraph* gtempClone(const char* name, const char* title)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   TIter next(gPad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) if (obj->IsA()->InheritsFrom("TGraph")) break; else obj=0;
   if (!obj) return 0;
   TGraph* gr = (TGraph*) obj->Clone();
   if (gr && name && *name) gr->SetName(name);
   if (gr && title && *title) gr->SetTitle(title);
   return gr;
}
// grtemp, grtempClone: the same as gtemp, gtempClone
TGraph* grtemp(const char* name, const char* title) {return gtemp(name,title);}
TGraph* grtempClone(const char* name, const char* title) {return gtempClone(name,title);}

TH1* hlast(const char* name, const char* title)
{
   TObjArray histos;
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   TIter next(gPad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) if (obj->IsA()->InheritsFrom("TH1")) histos.Add(obj); else obj=0;
   // pick up the last object from the array
   obj = histos.Last();
   if (!obj) return 0;
   TH1* h = (TH1*) obj;
   if (h && name && *name) h->SetName(name);
   if (h && title && *title) h->SetTitle(title);
   return h;
}
TH1* hfirst(const char* name, const char* title)   // was htemp originally
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   TIter next(gPad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) if (obj->IsA()->InheritsFrom("TH1")) break; else obj=0;
   if (!obj) return 0;
   TH1* h = (TH1*) obj;
   if (h && name && *name) h->SetName(name);
   if (h && title && *title) h->SetTitle(title);
   return h;
}
TH1* htemp(const char* name, const char* title)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   TIter next(gPad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) if (obj->IsA()->InheritsFrom("TH1")) break; else obj=0;
   if (!obj) return 0;
   TH1* h = (TH1*) obj;
   if (h && name && *name) h->SetName(name);
   if (h && title && *title) h->SetTitle(title);
   return h;
}
TH1* htempClone(const char* name, const char* title)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }

   TIter next(gPad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) if (obj->IsA()->InheritsFrom("TH1")) break; else obj=0;
   if (!obj) return 0;
   TH1* h = (TH1*) obj->Clone();
   if (h && name && *name) h->SetName(name);
   if (h && title && *title) h->SetTitle(title);
   return h;
}
TH2* htemp2(const char* name, const char* title)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   // there are 2 possible pad which appears as 2-dim histogram:
   // 1) in case of TH2::Draw the pad contains 2-dim histogram indeed
   // 2) in case of TTree::Draw the pad contains empty 2-dim histogram and TGraph
   TH2* h2 = 0;
   TGraph* gr = 0;
   //cout<< "gPad->GetName() = " << gPad->GetName() << endl;
   TIter next(gPad->GetListOfPrimitives());
   TObject* obj;
   while ((obj = next())) {
      //cout<< "obj->GetName() = " << obj->GetName() <<endl;
      if (obj->IsA()->InheritsFrom("TH2") && !h2) {
         h2 = (TH2*) obj;
      }
      if (obj->IsA()->InheritsFrom("TGraph") && !gr) {
         gr = (TGraph*) obj;
      }
   }
   if (!h2) return 0;
   if (h2 && !gr) {
      // looks like pure 2-dim histogram
      if (name && *name) h2->SetName(name);
      if (title && *title) h2->SetTitle(title);
      return h2;
   }
   if (h2 && gr) {
      // case of TTree::Draw
      // constract "normal" histogram filling up the h2
      // Note that we do not have information about the number of entries!
      for (int i=0; i<gr->GetN(); ++i) {
         h2->Fill(*(gr->GetX()+i), *(gr->GetY()+i));
      }
   }
   if (h2 && name && *name) h2->SetName(name);
   if (h2 && title && *title) h2->SetTitle(title);
   return h2;
}
TH2* htemp2Clone(const char* name, const char* title)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   // there are 2 possible pad which appears as 2-dim histogram:
   // 1) in case of TH2::Draw the pad contains 2-dim histogram indeed
   // 2) in case of TTree::Draw the pad contains empty 2-dim histogram and TGraph
   TH2* h2 = 0;
   TGraph* gr = 0;
   //cout<< "gPad->GetName() = " << gPad->GetName() << endl;
   TIter next(gPad->GetListOfPrimitives());
   TObject* obj;
   while ((obj = next())) {
      //cout<< "obj->GetName() = " << obj->GetName() <<endl;
      if (obj->IsA()->InheritsFrom("TH2") && !h2) {
         h2 = (TH2*) obj->Clone();      // clone the histo
      }
      if (obj->IsA()->InheritsFrom("TGraph") && !gr) {
         gr = (TGraph*) obj;
      }
   }
   if (!h2) return 0;
   if (h2 && !gr) {
      // looks like pure 2-dim histogram
      return h2;
   }
   if (h2 && gr) {
      // case of TTree::Draw
      // constract "normal" histogram filling up the h2
      // Note that we do not have information about the number of entries!
      for (int i=0; i<gr->GetN(); ++i) {
         h2->Fill(*(gr->GetX()+i), *(gr->GetY()+i));
      }
   }
   if (h2 && name && *name) h2->SetName(name);
   if (h2 && title && *title) h2->SetTitle(title);
   return h2;
}

TF1* ftemp(const char* name, const char* title)
{
   //cout<< "ftemp" <<endl;
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   //
   // Check TCanvas in this order:
   // 1) TF1: check the name if specified
   // 2) TGraph + empty TH2: check TGraph for TF1 (check the name if specified)
   // 3) TH1: check it for TF1 (check the name if specified)
   //

   // looking for TF1
   TIter next(gPad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) {
      if (obj->IsA()->InheritsFrom("TF1")) {
         if (name && *name) {
            if (strcmp(obj->GetName(),name)==0) {
               TF1* fun = (TF1*) obj;
               if (title) fun->SetTitle(title);
               return fun;              //-- success
            }
            else continue;              // wrong name
         }
         else {
            TF1* fun = (TF1*) obj;
            if (strcmp(fun->GetName(),"stats")==0) continue;
            if (title) fun->SetTitle(title);
            return fun;                 //-- success: name was not specified: take just first TF1
         }
      }
      else obj = 0;     // this is not a TF1
   }

   // looking for TGraph
   //cout<< "looking for TGraph" <<endl;
   next.Reset();
   while ((obj = next())) {
      if (obj->IsA()->InheritsFrom("TGraph")) {
         TGraph* gr = (TGraph*) obj;
         TIter next_fun(gr->GetListOfFunctions());
         TF1* fun;
         while ((fun = (TF1*) next_fun())) {
            if (name && *name) {
               if (strcmp(fun->GetName(), name)==0) {
                  if (title) fun->SetTitle(title);
                  return fun;   //-- success
               }
               else continue;           // wrong name
            }
            else {
               if (strcmp(fun->GetName(),"stats")==0) continue;
               if (title) fun->SetTitle(title);
               return fun;  //-- success: name was not specified: take just first TF1
            }
         }
      }
      else obj=0;
   }

   // looking for TH1
   //cout<< "looking for TH1" <<endl;
   next.Reset();
   while ((obj = next())) {
      if (obj->IsA()->InheritsFrom("TH1")) {
         TH1* h = (TH1*) obj;
         TIter next_fun(h->GetListOfFunctions());
         TF1* fun;
         while ((fun = (TF1*) next_fun())) {
            if (name && *name) {
               if (strcmp(fun->GetName(), name)==0) {
                  if (title) fun->SetTitle(title);
                  return fun;   //-- success
               }
               else continue;           // wrong name
            }
            else {
               if (strcmp(fun->GetName(),"stats")==0) continue;
               if (title) fun->SetTitle(title);
               return fun;  //-- success: name was not specified: take just first TF1
            }
         }
      }
      else obj=0;
   }
   return 0;
}
TF1* ftempClone(const char* name, const char* title)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   TF1* fun = ftemp(name, title);
   if (!fun) return 0;

   // try to clone this function.
   // NB: one can clone only these TF1 function which were created through 
   // constructor with formula. E.g. one can clone all predefined functions
   // See documentation on TF1 constructors.

   TF1* funClone = (TF1*) fun->Clone();
   if (funClone == 0) cout<< "*** ftempClone WARNING: this function cannot be cloned!" <<endl;
   return funClone;
}

TGraph* tgraph(TTree* tree, const char* name, const char* title) {
   // tree->Draw("y:x") creates two objects on the canvas:
   // 1) empty TH2F as just a graphical frame
   // 2) TGraph object which keeps data
   TGraph* gr = 0;
   TH1* h = tree->GetHistogram();
   if (h && h->IsA()->InheritsFrom("TH2"))      // sanity check
   {
      gr = new TGraph(tree->GetSelectedRows(), tree->GetV2(), tree->GetV1());
      gr->SetMarkerStyle(tree->GetMarkerStyle());
      gr->SetMarkerColor(tree->GetMarkerColor());
      if (name && *name) gr->SetTitle(name);
      else gr->SetName("grtemp");
      if (title && *title) gr->SetTitle(title);
      else gr->SetTitle(h->GetTitle());
   }
   return gr;
}

TMultiGraph* mgtemp(const char* name, const char* title)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   TIter next(gPad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) if (obj->IsA()->InheritsFrom("TMultiGraph")) break; else obj=0;
   if (!obj) return 0;
   TMultiGraph* mg = (TMultiGraph*) obj;
   if (mg && name && *name) mg->SetName(name);
   if (mg && title && *title) mg->SetTitle(title);
   return mg;
}

TPaveStats* stemp(TH1* h)     // return pointer to the stat box
{
   TPaveStats* stats = 0;

   if (!h) h = htemp();
   if (h) {
      stats = (TPaveStats*) h->GetFunction("stats");
   }
   else {
      // try to access the whole canvas
      if (gROOT->GetListOfCanvases()->GetEntries()) {
         stats = (TPaveStats*) gPad->GetPrimitive("stats");
      }
   }

   return stats;
}

TGraphErrors* gsubf(const TGraph* gr, const TF1* f, const char* name, const char* title)
{
   if (!gr) gr = gtemp();
   if (!gr) {
      cout<< "No TGraph found" <<endl;
      return 0;
   }
   if (!f) f = (TF1*) gr->GetListOfFunctions()->First();
   if (!f) {
      cout<< "No TF1 found" <<endl;
      return 0;
   }

   Double_t* x = gr->GetX();
   Double_t* y = gr->GetY();
   Double_t* ex = 0;
   Double_t* ey = 0;
   if (gr->InheritsFrom("TGraphErrors")) {
      ex = gr->GetEX();
      ey = gr->GetEY();
   }
   TGraphErrors* grsub = new TGraphErrors(gr->GetN(),x,y,ex,ey);
   if (name && *name) grsub->SetName(name);
   else grsub->SetName(Form("%s-%s",gr->GetName(),f->GetName()));
   if (title && *title) grsub->SetTitle(title);
   else grsub->SetTitle(Form("%s - %s",gr->GetTitle(),f->GetName()));
   grsub->SetMarkerStyle(gr->GetMarkerStyle());
   grsub->SetMarkerColor(gr->GetMarkerColor());
   grsub->SetLineStyle(gr->GetLineStyle());
   grsub->SetLineColor(gr->GetLineColor());

   for (int i=0; i<grsub->GetN(); ++i) {
      Double_t val_old = grsub->GetY()[i];
      Double_t val_fun = f->Eval(grsub->GetX()[i]);
      Double_t val_new = val_old - val_fun;
      grsub->SetPoint(i, grsub->GetX()[i], val_new);
   }
   return grsub;
}

Int_t axisNdivisions(Int_t Ndivisions) {
   TH1* h = htemp();
   if (h) {
      //if (Ndivisions == 0) cout<< "Ndivisions = " << h->GetXaxis()->GetNdivisions() <<endl;
      if (Ndivisions == 0) Ndivisions = h->GetXaxis()->GetNdivisions();
      else h->GetXaxis()->SetNdivisions(Ndivisions);
      return Ndivisions;
   }
   TGraph* gr = grtemp();
   if (gr) {
      //if (Ndivisions == 0) cout<< "Ndivisions = " << gr->GetXaxis()->GetNdivisions() <<endl;
      if (Ndivisions == 0) Ndivisions = gr->GetXaxis()->GetNdivisions();
      else gr->GetXaxis()->SetNdivisions(Ndivisions);
      return Ndivisions;
   }
   TF1* f = ftemp();
   if (f) {
      //if (Ndivisions == 0) cout<< "Ndivisions = " << f->GetXaxis()->GetNdivisions() <<endl;
      if (Ndivisions == 0) Ndivisions = f->GetXaxis()->GetNdivisions();
      else f->GetXaxis()->SetNdivisions(Ndivisions);
      return Ndivisions;
   }
   return 0;
}
Int_t ndi(Int_t Ndivisions) {
   if (Ndivisions > 0) axisNdivisions(Ndivisions);
   else {
      Ndivisions = axisNdivisions(0);
      axisNdivisions(--Ndivisions);
   }
   return Ndivisions;
}
void clone() {gPad->DrawClone();}

TPaveText* statstemp(TVirtualPad* can)
{
   if (can == 0) {
      if (gROOT->GetListOfCanvases()->GetEntries() > 0) can = gPad;
      else {
         cout<< "No canvas found" <<endl;
         return 0;
      }
   }
   TPaveText* stats = (TPaveText*) can->GetPrimitive("stats");
   return stats;
}

// TPaveText* titletemp(TVirtualPad* can)
// {
//    if (can == 0) {
//       if (gROOT->GetListOfCanvases()->GetEntries() > 0) can = gPad;
//       else {
//       cout<< "No canvas found" <<endl;
//       return 0;
//       }
//    }
//    TPaveText* title = (TPaveText*) can->GetPrimitive("title");
//    return title;
// }
// TPaveText* titletemp(const char* title, TVirtualPad* can)
// {
//    if (can == 0) {
//       if (gROOT->GetListOfCanvases()->GetEntries() > 0) can = gPad;
//       else {
//       cout<< "No canvas found" <<endl;
//       return 0;
//       }
//    }
//    TPaveText* tit = (TPaveText*) can->GetPrimitive("title");
//    tit->Clear();
//    tit->AddText(title);
//    return tit;
// }
// ----------- end of gtemp, htemp and ftemp ---------------

void drawpan(TH1* h) {
   if (h == 0) h = htemp();
   if (h) {
      h->DrawPanel();
      return;
   }
   TGraph* gr = gtemp();
   if (gr) {
      gr->DrawPanel();
      return;
   }
}

void fitpan(TH1* h) {
   if (h == 0) h = htemp();
   if (h) {
      h->FitPanel();
      return;
   }
   TGraph* gr = gtemp();
   if (gr) {
      gr->FitPanel();
      return;
   }
}

void InitPredefinedFunctions()
{
   TF1 *f1;
   if (!gROOT->GetListOfFunctions()->FindObject("gaus"))    {f1 = new TF1("gaus","gaus",-1,1);       f1->SetParameters(1,0,1);}
   if (!gROOT->GetListOfFunctions()->FindObject("gausn"))   {f1 = new TF1("gausn","gausn",-1,1);     f1->SetParameters(1,0,1);}
   if (!gROOT->GetListOfFunctions()->FindObject("landau"))  {f1 = new TF1("landau","landau",-1,1);   f1->SetParameters(1,0,1);}
   if (!gROOT->GetListOfFunctions()->FindObject("landaun")) {f1 = new TF1("landaun","landaun",-1,1); f1->SetParameters(1,0,1);}
   if (!gROOT->GetListOfFunctions()->FindObject("expo"))    {f1 = new TF1("expo","expo",-1,1);       f1->SetParameters(1,1);}
   if (!gROOT->GetListOfFunctions()->FindObject("pol0")) {
      for (Int_t i=0;i<10;i++) {
         f1 = new TF1(Form("pol%d",i),Form("pol%d",i),-1,1);
         f1->SetParameters(1,1,1,1,1,1,1,1,1,1);
      }
   }
}

void fitgaus(Double_t xmin, Double_t xmax, const char* opt, const char* gopt, TH1* h, TGraph* gr)
{
   // histogram part
   if (h == 0) h = htemp();
   if (h && h->GetEntries() == 0) h = 0;        // seems, this is just 2-dim frame for the plot
   if (h)
   {
      if (xmax == 0 and xmax < xmin) xmax = h->GetXaxis()->GetXmax();
      if (xmin == 0 and xmax == 0) h->Fit("gaus", opt, gopt);
      else {
         TF1* gaus = new TF1("gaus", "gaus", xmin, xmax);
         h->Fit(gaus, Form("%sR",opt), gopt);
      }
      gPad->Modified();
      gPad->Update();
      return;
   }

   // graph part
   if (gr == 0) gr = gtemp();
   if (gr)
   {
      if (xmax == 0 and xmax < xmin) xmax = gr->GetXaxis()->GetXmax();
      if (xmin == 0 and xmax == 0) gr->Fit("gaus", opt, gopt);
      else {
         TF1* gaus = new TF1("gaus", "gaus", xmin, xmax);
         gr->Fit(gaus, Form("%sR",opt), gopt);
      }
      gPad->Modified();
      gPad->Update();
      return;
   }
} 
void fitg(Double_t xmin, Double_t xmax, const char* opt, const char* gopt, TH1* h, TGraph* gr) {
   fitgaus(xmin,xmax,opt,gopt,h,gr);
}
void fitgl(Double_t xmin, Double_t xmax, const char* opt, const char* gopt, TH1* h, TGraph* gr) {
   fitgaus(xmin,xmax,opt,gopt,h,gr);
   leftgaus();
}
void fitgr(Double_t xmin, Double_t xmax, const char* opt, const char* gopt, TH1* h, TGraph* gr) {
   fitgaus(xmin,xmax,opt,gopt,h,gr);
   rightgaus();
}

void fitexpo(Double_t xmin, Double_t xmax, const char* opt, const char* gopt, TH1* h, TGraph* gr)
{
   // histogram part
   if (h == 0) h = htemp();
   if (h && h->GetEntries() == 0) h = 0;        // seems, this is just 2-dim frame for the plot
   if (h)
   {
      if (xmax == 0 and xmax < xmin) xmax = h->GetXaxis()->GetXmax();
      if (xmin == 0 and xmax == 0) h->Fit("expo", opt, gopt);
      else {
         TF1* expo = new TF1("expo", "expo", xmin, xmax);
         h->Fit(expo, Form("%sR",opt), gopt);
      }
      gPad->Modified();
      gPad->Update();
      return;
   }

   // graph part
   if (gr == 0) gr = gtemp();
   if (gr)
   {
      if (xmax == 0 and xmax < xmin) xmax = gr->GetXaxis()->GetXmax();
      if (xmin == 0 and xmax == 0) gr->Fit("expo", opt, gopt);
      else {
         TF1* expo = new TF1("expo", "expo", xmin, xmax);
         gr->Fit(expo, Form("%sR",opt), gopt);
      }
      gPad->Modified();
      gPad->Update();
      return;
   }
} 
void fite(Double_t xmin, Double_t xmax, const char* opt, const char* gopt, TH1* h, TGraph* gr) {
   fitexpo(xmin,xmax,opt,gopt,h,gr);
}
void fitel(Double_t xmin, Double_t xmax, const char* opt, const char* gopt, TH1* h, TGraph* gr) {
   fitexpo(xmin,xmax,opt,gopt,h,gr);
   leftgaus();
}
void fiter(Double_t xmin, Double_t xmax, const char* opt, const char* gopt, TH1* h, TGraph* gr) {
   fitexpo(xmin,xmax,opt,gopt,h,gr);
   rightgaus();
}

// void fitpol(Double_t xmin, Double_t xmax, Int_t power, const char* opt, const char* gopt, TH1* h, TGraph* gr)
// {
//    if (h == 0 && gr == 0) h = htemp();
//    if (h && h->GetEntries() == 0) h = 0;             // case of Tree::Draw for 2-dim histo
//    if (h) {
//       if (xmax == 0 and xmax < xmin) xmax = h->GetXaxis()->GetXmax();
//       if (xmin == 0 and xmax == 0) h->Fit(Form("pol%d",power), opt, gopt);
//       else {
//       TF1* pol = new TF1(Form("pol%d",power),Form("pol%d",power), xmin, xmax);
//       h->Fit(pol, Form("%sR",opt), gopt);
//       }
//    }
//    else {
//       if (gr == 0) gr = gtemp();
//       if (gr) {
//       if (xmax == 0 and xmax < xmin) xmax = gr->GetXaxis()->GetXmax();
//       if (xmin == 0 and xmax == 0) gr->Fit(Form("pol%d",power), opt, gopt);
//       else {
//          TF1* pol = new TF1(Form("pol%d",power),Form("pol%d",power), xmin, xmax);
//          gr->Fit(pol, Form("%sR",opt), gopt);
//       }
//       }
//    }
//    gPad->Modified();
//    gPad->Update();
// }
void fitpol(Double_t xmin, Double_t xmax, Int_t power, const char* opt, const char* gopt, TH1* h, TGraph* gr)
{
   if (!gROOT->GetListOfFunctions()->FindObject("pol0")) InitPredefinedFunctions();

   if (h == 0 && gr == 0) h = htemp();
   if (h && h->GetEntries() == 0) h = 0;                // case of Tree::Draw for 2-dim histo
   if (h) {
      if (xmax == 0 and xmax < xmin) xmax = h->GetXaxis()->GetXmax();
      if (xmin == 0 and xmax == 0) h->Fit(Form("pol%d",power), opt, gopt);
      else {
         // TF1* pol = new TF1(Form("pol%d",power),Form("pol%d",power), xmin, xmax);
         // h->Fit(pol, Form("%sR",opt), gopt);
         h->Fit(Form("pol%d",power), opt, gopt, xmin, xmax);
      }
   }
   else {
      if (gr == 0) gr = gtemp();
      if (gr) {
         if (xmax == 0 and xmax < xmin) xmax = gr->GetXaxis()->GetXmax();
         if (xmin == 0 and xmax == 0) gr->Fit(Form("pol%d",power), opt, gopt);
         else {
            // TF1* pol = new TF1(Form("pol%d",power),Form("pol%d",power), xmin, xmax);
            // gr->Fit(pol, Form("%sR",opt), gopt);
            gr->Fit(Form("pol%d",power), opt, gopt, xmin, xmax);
         }
      }
   }
   gPad->Modified();
   gPad->Update();
}
void fitp(Double_t xmin, Double_t xmax, Int_t power, const char* opt, const char* gopt, TH1* h, TGraph* gr) {
   fitpol(xmin,xmax,power,opt,gopt,h,gr);
}

Double_t fit_parameter(const char* fun_name, const char* par_name)
{
   // try TGraph
   TH1* h = htemp();
   if (h) {
      TF1* f = h->GetFunction(fun_name);
      if (f) return f->GetParameter(par_name);
   }
   // try TGraph
   TGraph* gr = gtemp();
   if (gr) {
      TF1* f = gr->GetFunction(fun_name);
      if (f) return f->GetParameter(par_name);
   }
   return 0;
}
Double_t fit_parameter(const char* fun_name, Int_t number)
{
   // try TGraph
   TH1* h = htemp();
   if (h) {
      TF1* f = h->GetFunction(fun_name);
      if (f) return f->GetParameter(number);
   }
   // try TGraph
   TGraph* gr = gtemp();
   if (gr) {
      TF1* f = gr->GetFunction(fun_name);
      if (f) return f->GetParameter(number);
   }
   return 0;
}
Double_t fit_parameter_error(const char* fun_name, Int_t number)
{
   // try TGraph
   TH1* h = htemp();
   if (h) {
      TF1* f = h->GetFunction(fun_name);
      if (f) return f->GetParError(number);
   }
   // try TGraph
   TGraph* gr = gtemp();
   if (gr) {
      TF1* f = gr->GetFunction(fun_name);
      if (f) return f->GetParError(number);
   }
   return 0;
}
// gaus parameters
Double_t gaus_ampl() {return fit_parameter("gaus", "Constant");}
Double_t gaus_mean() {return fit_parameter("gaus", "Mean");}
Double_t gaus_sigma() {return fit_parameter("gaus", "Sigma");}
// gaus parameter errors
Double_t gaus_eampl() {return fit_parameter_error("gaus", 0);}
Double_t gaus_emean() {return fit_parameter_error("gaus", 1);}
Double_t gaus_esigma() {return fit_parameter_error("gaus", 2);}

void pargaus(Double_t& a, Double_t& mean, Double_t& sigma, const char* hname)
{
   a = mean = sigma = 0;

   TH1* h = htemp(hname);
   if (h) {
      TF1* gaus = h->GetFunction("gaus");
      if (gaus) {
         a = gaus->GetParameter(1);
         mean = gaus->GetParameter(1);
         sigma = gaus->GetParameter(2);
      }
   }
}

void npe(TH1* h, Double_t* p_Npe, Double_t* p_eNpe)
{
   if (p_Npe) *p_Npe = 0;
   if (p_eNpe) *p_eNpe = 0;
   if (h) {
      TF1* gaus = h->GetFunction("gaus");
      if (gaus) {
         Double_t mean = gaus->GetParameter(1);
         Double_t sigma = gaus->GetParameter(2);
         Double_t emean = gaus->GetParError(1);
         Double_t esigma = gaus->GetParError(2);
         Double_t ratio = mean/sigma;
         Double_t Npe = ratio*ratio;
         Double_t ra = emean/mean;
         Double_t rb = esigma/sigma;
         Double_t eNpe = 2.*Npe*TMath::Sqrt(ra*ra + rb*rb);
         // cout<< "Npe = " << Npe << " +- " << eNpe << "\t mean = " << mean << " sigma = " << sigma <<endl;
         // cout<< setprecision(3) << "Npe = " << Npe << " +- " << eNpe << "\t mean = " << mean << " sigma = " << sigma <<endl;
         printf("Npe = %0.3f +- %0.3f \t mean = %0.3f sigma = %0.3f\n", Npe,eNpe, mean,sigma);
         if (p_Npe) *p_Npe = Npe;
         if (p_eNpe) *p_eNpe = eNpe;
      }
   }
}

void npe(const char* hname, Double_t* p_Npe, Double_t* p_eNpe)
{
   if (p_Npe) *p_Npe = 0;
   if (p_eNpe) *p_eNpe = 0;
   TH1* h = htemp(hname);
   if (h) {
      TF1* gaus = h->GetFunction("gaus");
      if (gaus) {
         Double_t mean = gaus->GetParameter(1);
         Double_t sigma = gaus->GetParameter(2);
         Double_t emean = gaus->GetParError(1);
         Double_t esigma = gaus->GetParError(2);
         Double_t ratio = mean/sigma;
         Double_t Npe = ratio*ratio;
         Double_t ra = emean/mean;
         Double_t rb = esigma/sigma;
         Double_t eNpe = 2.*Npe*TMath::Sqrt(ra*ra + rb*rb);
         // cout<< "Npe = " << Npe << " +- " << eNpe << "\t mean = " << mean << " sigma = " << sigma <<endl;
         // cout<< setprecision(3) << "Npe = " << Npe << " +- " << eNpe << "\t mean = " << mean << " sigma = " << sigma <<endl;
         printf("Npe = %0.3f +- %0.3f \t mean = %0.3f sigma = %0.3f\n", Npe,eNpe, mean,sigma);
         if (p_Npe) *p_Npe = Npe;
         if (p_eNpe) *p_eNpe = eNpe;
      }
   }
}

void pic(TString pathname, TString suffix)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return;
   }
   if (pathname == "") {
      if (htemp()) pathname = htemp()->GetName();
      else if (grtemp()) pathname = grtemp()->GetName();
      else if (ftemp()) pathname = ftemp()->GetName();
      else pathname = gPad->GetName();
   }
   gPad->SaveAs(pathname + suffix + ".eps");
   gPad->SaveAs(pathname + suffix + ".png");
}
// just eps
void eps(TString pathname, TString suffix)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return;
   }
   if (pathname == "") {
      if (htemp()) pathname = htemp()->GetName();
      else if (grtemp()) pathname = grtemp()->GetName();
      else if (ftemp()) pathname = ftemp()->GetName();
      else pathname = gPad->GetName();
   }
   gPad->SaveAs(pathname + suffix + ".eps");
}
// just png
void png(TString pathname, TString suffix, bool check_fexists)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return;
   }
   if (pathname == "") {
      if (htemp()) pathname = htemp()->GetName();
      else if (grtemp()) pathname = grtemp()->GetName();
      else if (ftemp()) pathname = ftemp()->GetName();
      else pathname = gPad->GetName();
   }
   TString ofname = pathname + suffix + ".png";
   if (check_fexists and fexists(ofname.Data())) {
      cout<< "File exists: \"" << ofname << "\"" <<endl;
      cout<< "Enter <CR> to overwrite, any other key to quit: ";
      std::string str;
      std::getline(cin, str);    // read complete line to avoid problems with the next call
      if (str.size() > 0) {
         cout<< "*** quit without saving" <<endl;
         return;
      }
   }
   gPad->SaveAs(pathname + suffix + ".png");
}
void pngall() {
   for (int i=0; i<gROOT->GetListOfCanvases()->GetEntries(); i++) {
      TCanvas* can=(TCanvas*)gROOT->GetListOfCanvases()->At(i);
      can->SaveAs(Form("%s.png",can->GetName()));
   }
}

Double_t intgaus(TH1* h)
{
   if (h == 0) h = htemp();

   Double_t gampl = h->GetFunction("gaus")->GetParameter(0);
   // Double_t gmean = h->GetFunction("gaus")->GetParameter(1);
   Double_t gsigma = h->GetFunction("gaus")->GetParameter(2);
   Double_t nevt_w = TMath::Sqrt(2.*TMath::Pi()) * gampl*gsigma / h->GetBinWidth(1);
   cout<< "TMath::Sqrt(2.*TMath::Pi()) * gampl*gsigma / h->GetBinWidth(1) = " << nevt_w <<endl;

   return nevt_w;
}

TH1* hname(const char* name)
{
   if (!name || !*name) {
      gDirectory->ls();
      return 0;
   }
   TH1* h = (TH1*) gDirectory->Get(name);
   if (!h) cout<< "--> Histogram not found: " << name <<endl;
   return h;
}

void print_htemp(TH1* h)
{
   if (!h) h = htemp();
   if (!h) return;

   cout<< "h->GetEntries() = " << h->GetEntries() <<endl;
   cout<< "sum of bin contents: h->Integral() = " << h->Integral() <<endl;

   Float_t integral = 0;
   Float_t integral_count = 0;
   for (int i=1; i<=h->GetNbinsX(); ++i) {
      Float_t con = h->GetBinContent(i);
      integral += con*h->GetBinCenter(i);
      integral_count += con*(i-1);
      if (con) {
         cout<< i << "\t " << h->GetBinCenter(i) << "\t " << con <<endl;
      }
   }
   cout<< "integer count integral con*(i-1) = " << integral_count <<endl;
   cout<< "integral con*h->GetBinCenter(i)  = " << integral <<endl;
}
// void hpri(TH1* h=0)
void hprint(TH1* h)
{
   if (!h) h = htemp();
   if (!h) return;

   //cout<< "h->GetEntries() = " << h->GetEntries() <<endl;
   //cout<< "sum of bin contents: h->Integral() = " << h->Integral() <<endl;

   Float_t integral = 0;
   //Float_t integral_count = 0;
   for (int i=1; i<=h->GetNbinsX(); ++i) {
      Float_t con = h->GetBinContent(i);
      integral += con*h->GetBinLowEdge(i);
      //integral_count += con*(i-1);
      if (con) {
         // cout<< i << "\t " << h->GetBinCenter(i) << "\t " << con <<endl;
         cout<< i << "\t " << h->GetBinLowEdge(i) << "\t " << con <<endl;
      }
   }
   //cout<< "integer count integral con*(i-1) = " << integral_count <<endl;
   cout<< "integral con*h->GetBinLowEdge(i)     = " << integral <<endl;
}

//void tax(const char* format="%d-%H")
void tax(const char* format)
{
   // special option format = 0 (NULL) to reset(remove) time axis
   /*
      Defines the format of the labels along the time axis. 
      It can be changed using the TAxis method SetTimeFormat. 
      The time format is the one used by the C function strftime(). 
      It's a string containing the following formatting characters :
    */
   if (strcmp(format, "help")==0)
   {
      cout<< "   * for date :                                                                                   " <<endl;
      cout<< "                    o %a abbreviated weekday name         " <<endl;
      cout<< "                    o %b abbreviated month name           " <<endl;
      cout<< "                    o %d day of the month (01-31)         " <<endl;
      cout<< "                    o %m month (01-12)                    " <<endl;
      cout<< "                    o %y year without century             " <<endl;
      cout<< "                    o %Y year with century                " <<endl;
      cout<< "   * for time :                                           " <<endl;
      cout<< "                    o %H hour (24-hour clock)             " <<endl;
      cout<< "                    o %I hour (12-hour clock)             " <<endl;
      cout<< "                    o %p local equivalent of AM or PM     " <<endl;
      cout<< "                    o %M minute (00-59)                   " <<endl;
      cout<< "                    o %S seconds (00-61)                  " <<endl;
      cout<< "                    o %% %                                " <<endl;
      return;
   }
   /*
      The other characters are output as is. For example to have a format like dd/mm/yyyy one should do:
      h->GetXaxis()->SetTimeFormat("%d\/%m\/%Y");
      If the time format is not defined, a default one will be computed automatically.
    */
   TAxis* axis = 0;
   if (htemp()) {
      axis = htemp()->GetXaxis();
   }
   if (!axis) {
      if (gtemp()) {
         axis = gtemp()->GetXaxis();
      }
   }
   if (axis) {
      if (format) {
         axis->SetTimeDisplay(1);
         axis->SetTimeFormat(format);
      }
      else {
         // special option format = 0 (NULL) to reset time axis
         axis->SetTimeDisplay(0);
      }
      gPad->Modified();
      gPad->Update();
   }
   else cout<< "***error: no hist or graph found" <<endl;
}

void nul(Axis_t min)
{
   TIter next(gPad->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) if (obj->IsA()->InheritsFrom("TH1")) break; else obj=0;
   if (obj) {
      TH1* h = (TH1*) obj;
      h->SetMinimum(min);
      h->Draw();
   }
}
// void nul(Axis_t min)
// {
//    TIter next(gPad->GetListOfPrimitives());
//    TObject* obj = 0;
//    while (obj = next()) if (obj->IsA()->InheritsFrom("TH1")) break; else obj=0;
//    if (obj) {
//       if (obj->IsA()->InheritsFrom("TH2")) {
//       TH2* h2 = (TH2*) obj;
//       Axis_t ymin = min;
//       Axis_t ymax = h2->GetYaxis()->GetXmax() * 1.10;
//       Axis_t xmin = h2->GetXaxis()->GetXmin();
//       Axis_t xmax = h2->GetXaxis()->GetXmax();
//       gPad->DrawFrame(xmin,ymin, xmax,ymax);
//       //h2->Draw("same hist");
//       h2->Draw("same");
//       }
//       else {
//       TH1* h1 = (TH1*) obj;
//       h1->SetMinimum(min);
//       h1->Draw();
//       }
//    }
// }
// void nul(Axis_t min)
// {
//    TIter next(gPad->GetListOfPrimitives());
//    TObject* obj = 0;
//    while (obj = next()) if (obj->IsA()->InheritsFrom("TH2")) break; else obj=0;
//    if (obj) {
//       if (obj->IsA()->InheritsFrom("TH2")) {
//       TH2* h2 = (TH2*) obj;
//       Axis_t ymin = min;
//       Axis_t ymax = h2->GetYaxis()->GetXmax() * 1.10;
//       Axis_t xmin = h2->GetXaxis()->GetXmin();
//       Axis_t xmax = h2->GetXaxis()->GetXmax();
//       gPad->DrawFrame(xmin,ymin, xmax,ymax);
//       //h2->Draw("same hist");
//       h2->Draw("same");
//       }
//       else {
//       while (obj = next()) if (obj->IsA()->InheritsFrom("TH1")) break; else obj=0;
//       if (obj) {
//          TH1* h1 = (TH1*) obj;
//          h1->SetMinimum(min);
//          h1->Draw();
//       }
//       }
//    }
// }

void ct()
{
   // substitute canvas window title by object title
   TIter next(gROOT->GetListOfCanvases());
   TCanvas* c;
   while ((c = (TCanvas*) next())) {
      //cout<< c->GetName() <<endl;
      TIter next_obj(c->GetListOfPrimitives());
      TObject* obj = 0;
      while ((obj = next_obj())) {
         if ( false
               || obj->IsA()->InheritsFrom("TH1")
               || obj->IsA()->InheritsFrom("TGraph")
               || obj->IsA()->InheritsFrom("TMultiGraph")
               || obj->IsA()->InheritsFrom("TF1")
            ) break;
         else obj=0;
      }
      if (obj) {
         const char* tit = obj->GetTitle();
         //cout<< "      " << tit <<endl;
         c->SetTitle(tit);
      }
   }
}
void ctadd()
{
   // add to canvas window title an object title
   TIter next(gROOT->GetListOfCanvases());
   TCanvas* c;
   while ((c = (TCanvas*) next())) {
      //cout<< c->GetName() <<endl;
      TIter next_obj(c->GetListOfPrimitives());
      TObject* obj = 0;
      while ((obj = next_obj())) {
         if ( false
               || obj->IsA()->InheritsFrom("TH1")
               || obj->IsA()->InheritsFrom("TGraph")
               || obj->IsA()->InheritsFrom("TMultiGraph")
               || obj->IsA()->InheritsFrom("TF1")
            ) break;
         else obj=0;
      }
      if (obj) {
         const char* tit = obj->GetTitle();
         //cout<< "      " << tit <<endl;
         c->SetTitle(Form("%s %s",c->GetName(),tit));
      }
   }
}
void ctlist()
{
   // list canvas window title and object title
   TIter next(gROOT->GetListOfCanvases());
   TCanvas* c;
   while ((c = (TCanvas*) next())) {
      //cout<< c->GetName() <<endl;
      TIter next_obj(c->GetListOfPrimitives());
      TObject* obj = 0;
      while ((obj = next_obj())) {
         if ( false
               || obj->IsA()->InheritsFrom("TH1")
               || obj->IsA()->InheritsFrom("TGraph")
               || obj->IsA()->InheritsFrom("TMultiGraph")
               || obj->IsA()->InheritsFrom("TF1")
            ) break;
         else obj=0;
      }
      if (obj) {
         //const char* tit = obj->GetTitle();
         //cout<< "   " << tit <<endl;
         //c->SetTitle(Form("%s %s",c->GetName(),tit));
         cout<< c->GetTitle() << "\t " << obj->GetTitle() <<endl;
      }
   }
}
const char* addtit(const char* title, TCanvas* can)
{
   // adds text to existing canvas title to have e.g. "c1_n3 before the fit"
   if (can == 0) can = (TCanvas*) gPad; // conversion from ‘TVirtualPad*’ to ‘TCanvas*
   if (can == 0) return 0;

   std::stringstream ss;
   ss.str("");
   ss << can->GetTitle() << " " << title;
   can->SetTitle(ss.str().c_str());
   return can->GetTitle();
}

TChain* chlist(const char* tree_name, const char* filelist) {
   TChain* ch = new TChain(tree_name);
   ifstream f(filelist);
   if (!f.is_open()) {
      cout<< "Filelist not found: " << filelist <<endl;
      return 0;
   }
   std::string str;
   while (f >> str) {
      // // check if file exist (TChain::Add do not check it)
      // ifstream ftmp(str.c_str());
      // if (!ftmp.is_open()) {
      //    cout<< "File not found: " << str <<endl;
      //    return 0;
      // }
      // else ftmp.close();
      //  add to chain
      ch->Add(str.c_str(),-1);
   }
   cout<< "Total events: " << ch->GetEntries() <<endl;
   return ch;
}

Int_t mkread_Ngen = 0;  // Global counter for Ngen
//
TTree* mkread(const char* fname, bool addFriend)
{
   mkread_Ngen = 0;     // clear global counter

   TTree* mk = 0;
   TDirectory* dir = gDirectory;                        // current dir
   TString ts(fname);
   if (ts.EndsWith(".root")) {
      // root file
      //mk = get_mk(fname);
      TFile* f = gROOT->GetFile(fname);
      if (!f) f = TFile::Open(fname);
      else        f->cd();
      if (!f) {
         cout<< "*** ERROR mkread: file not found: " << fname <<endl;
         return 0;
      }
      mk  = (TTree*) gDirectory->Get("mk");
      if (addFriend) {
         TTree* tmb = (TTree*) gDirectory->Get("TMBTree");
         if (tmb) {
            cout<< "   Adding friend TMBTree for " << fname <<endl;
            mk->AddFriend(tmb);
         }
      }
   }
   else {
      // filelist
      //mk = chlist("mk", fname);
      TChain* chain_mk  = 0;
      TChain* chain_tmb = 0;
      chain_mk = new TChain("mk");
      if (addFriend) chain_tmb = new TChain("TMBTree");
      std::ifstream f(fname);
      // if (!f.is_open()) {
      //    cout<< "chlist: File not found: " << fname <<endl;
      //    return 0;
      // }

      char line[255];
      Int_t len = sizeof(line);
      mkread_Ngen = 0;
      while (f.getline(line,len)) {
         if (strlen(line) == 0 || line[0] == '#') continue;
         std::stringstream ss(line);
         std::string fname_current;
         Int_t  ngen_current = 0;
         ss >> fname_current >> ngen_current;
         mkread_Ngen += ngen_current;
         //cout<< "fname_current = " << fname_current << " ngen_current = " << ngen_current << " mkread_Ngen = " << mkread_Ngen <<endl;
         Int_t n = chain_mk->Add(fname_current.c_str(), -1);  // second parameter to ensure file exists
         if (n == 0) {
            cout<< "*** ERROR chlist_lq while processing filelist " << fname << ". File not found: " << fname_current <<endl;
            return 0;
         }
         if (addFriend) n = chain_tmb->Add(fname_current.c_str(), -1);
         cout<< "Added " << fname_current << " Total entries: " << chain_mk->GetEntries();
         if (ngen_current) cout<< "   Ngen = " << mkread_Ngen;
         cout<<endl;
      }
      if (addFriend) chain_mk->AddFriend(chain_tmb,""); // ROOT 5.16 complains w/o ""
      mk = chain_mk;
   }

   dir->cd();           // back to starting dir
   Long64_t nentries = mk->GetEntries();
   cout<< "Total entries: " << nentries;
   if (mkread_Ngen) cout<< "     Ngen = " << mkread_Ngen;
   cout<<endl;
   //-- mk->SetEstimate(nentries);
   return mk;
}

Stat_t hint(char* hname, Float_t x1, Float_t x2)
{
   TObject* obj = (TObject*) gDirectory->Get(hname);
   if (obj && strcmp(obj->IsA()->GetName(),TH1F::Class()->GetName())==0)
   {
      TH1F* h = (TH1F*) obj;
      Int_t bin1 = h->GetXaxis()->FindBin(x1);
      if (x2 == 0) x2 = h->GetXaxis()->GetXmax();
      Int_t bin2 = h->GetXaxis()->FindBin(x2);
      Stat_t sum = h->Integral(bin1,bin2);
      Stat_t tot = h->Integral(h->GetXaxis()->GetFirst(),h->GetXaxis()->GetLast());
      //cout<< "Intergal is " << sum << " ==> " << 100.*Float_t(sum)/Float_t(tot) << " %" <<endl;
      cout<< "Intergal from " << x1 << " to " << x2 << " is " << sum << " ==> " << 100.*Float_t(sum)/Float_t(tot) << " %" <<endl;
      return sum;
   }
   else {
      cout<< "Could not found histogram " << hname <<endl;
      return 0;
   }
}

Stat_t hint(TH1F* h, Float_t x1, Float_t x2)
{
   Int_t bin1 = h->GetXaxis()->FindBin(x1);
   if (x2 == 0) x2 = h->GetXaxis()->GetXmax();
   Int_t bin2 = h->GetXaxis()->FindBin(x2);
   Stat_t sum = h->Integral(bin1,bin2);
   Stat_t tot = h->Integral(h->GetXaxis()->GetFirst(),h->GetXaxis()->GetLast());
   cout<< "Intergal is " << sum << " ==> " << 100.*Float_t(sum)/Float_t(tot) << " %" <<endl;
   return sum;
}

/*
t->Draw(">>elist","np==1&&Nt(0.1)>0&&Nt(0.1)<=3"); TEventList* elist=(TEventList*) gDirectory->Get("elist"); for(int i=0; i<10; ++i) {cout<< elist->GetEntry(i) <<endl;}
Int_t n = -1                           // to have n to be equal to the current event number
cloop("FillDrawCells(elist->GetEntry(++n))")

----> to see prev event issue command n-=2
*/
void cloop(const char* command)
{
   TTimer timer("gSystem->ProcessEvents();",50,kFALSE);  // process mouse events every 50 ms

   while (true)
   {
      gROOT->ProcessLine(command);

      timer.Start();    // start processing of mouse events on TCanvas
      cout<< "<CR>=next, Clone, Quit, command ";
      std::string line;
      std::getline(cin, line);
      timer.Stop();     // disable timer after <CR>

      if (line.size() == 0) continue;           // <CR>

      if (line == "Q" or line == "q") break;    // Q or q

      if (line == "C" or line == "c") {         // C or c
         gPad->DrawClone();
         continue;
      }

      if (line == ".q") {
         cout<< "--> Confirm command to quit ROOT: ";
         std::getline(cin, line);
      }
      gROOT->ProcessLine(line.c_str());         // interpret as a command
   }
}

///////////////////////////////////////////////////////////////////////////////////////

// Float_t dRphi_eta(Float_t phi1, Float_t eta1, Float_t phi2, Float_t eta2) {
//    Float_t dphi = TMath::Abs(phi1 - phi2);
//    if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
//    Float_t deta = eta1 - eta2;
//    return TMath::Sqrt(dphi*dphi + deta*deta);
// }
// 
// void wait(char* mess=0) {
//    if (mess) cout << mess <<endl;
//    cout << "<CR>=Continue, Q=Quit ";
//    char ch = 'Q';
//    ch=getchar();
//    if (ch == 'Q' || ch == 'q') exit(1);
// }
// 
// void hpan(Int_t x, Int_t y)
//    // Generates macro with control bar to draw all histos in memory
//    // zatserkl@fnal.gov Jun 8, 2002
//    // cf http://root.cern.ch/cgi-bin/print_hit_bold.pl/root/roottalk/roottalk99/1228.html
// {
//    TControlBar* hpan = new TControlBar("vertical", "hpan", x, y);
//    hpan->AddButton("grid", "gPad->SetGrid(); gPad->Draw();", "");
// 
//    TIter nextobj(gDirectory->GetList());
//    TObject *obj;
//    while ((obj = (TObject *) nextobj())) {
//       if (obj->IsA()->InheritsFrom("TH1")) {
//       const char* id = obj->GetName();
//       char action[256];
//       strcpy(action,id); strcat(action,"->Draw()");
//       hpan->AddButton(id, action, "");
//       }
//    }
//    hpan->Show();
//    gROOT->SaveContext();
// }
// 
// void hpan(char* file, Int_t x=0, Int_t y=0)
//    // Generates macro with control bar to draw all histos in memory
//    // zatserkl@fnal.gov Jun 8, 2002
//    // cf http://root.cern.ch/cgi-bin/print_hit_bold.pl/root/roottalk/roottalk99/1228.html
// {
//    ofstream fout(file, ios_base::out);
//    fout << "{" << endl;
//    fout << "histo = new TControlBar(\"vertical\", \"histo\","
//       << x << "," << y << ");" << endl;
//    fout << "histo->AddButton(\"grid\", \"gPad->SetGrid(); gPad->Draw();\", \"\");" << endl;
// 
//    TIter nextobj(gDirectory->GetList());
//    TObject *obj;
//    while ((obj = (TObject *) nextobj())) {
//       if (obj->IsA()->InheritsFrom("TH1")) {
//       const char* str = obj->GetName();
//       fout<< "histo->AddButton(" <<'\"'<< str
//          <<"\", \"" << str <<"->Draw();\", \"\");"
//          <<endl;
//       }
//    }
//    fout << "histo->Show();" << endl;
//    fout << "gROOT->SaveContext();" << endl;
//    fout << "}" << endl;
//    fout.close();
// }
// 
// void hsave(char* file)
//    // Saves all histos present in memory
//    // zatserkl@fnal.gov Jun 20, 2002
//    // cf http://root.cern.ch/cgi-bin/print_hit_bold.pl/root/roottalk/roottalk99/1228.html
// {
//    TObjArray hlist(0);
// 
//    TIter nextobj(gDirectory->GetList());
//    TObject *obj;
//    while ((obj = (TObject *) nextobj())) {
//       if (obj->IsA()->InheritsFrom("TH1")) {
//       hlist.Add(obj);
//       //const char* str = obj->GetName();
//       //cout<< "Adding to list " << str <<endl;
//       }
//    }
//    TFile f;
//    if (f.Open(file,"RECREATE")) {
//       hlist.Write();
//       f.Close();
//       cout << "Histos are saved in file " << file << endl;
//    }
//    else cout<< "Could not open file " << file << " to save histograms" <<endl;
//    gROOT->cd();
// }
// 
// void hps(char* fname, TCanvas* can=NULL)
//    // Saves in PostScript all histos present in memory
//    // zatserkl@fnal.gov Jun 23, 2002
//    // http://root.cern.ch/root/htmldoc/src/TPostScript.cxx.html
// {
//    char* file = new char[strlen(fname)+3+1]; strcpy(file,fname);
//    char* pext = fname + strlen(fname) - 3;
//    if (pext <= 0 || strcmp(pext, ".ps") != 0) strcat(file, ".ps");
// 
//    cout<< "hps: Creating PostScript file " << file <<endl;
//    // Create a canvas first!
//    if (can == NULL) {
//       cout<< "Creating canvas" <<endl;
//       can = new TCanvas("psmaker","psmaker");
//    }
//    else cout<< "hps: Using current canvas" <<endl;
// 
//    //Int_t type = 111;        //portrait  ps
//    Int_t type = 112;   //landscape ps
//    //Int_t type = 113;        //eps
// 
//    //-- NB: canvas should exist for successful TPostScript!
//    TPostScript ps(file,type);
//    //ps.Range(16,24);  //set x,y of printed page
// 
//    TIter nextobj(gDirectory->GetList());
//    TObject *obj;
//    while ((obj = (TObject *) nextobj())) {
//       if (obj->IsA()->InheritsFrom("TH1")) {
//       obj->Draw();
//       can->Update();
//       }
//    }
//    ps.Close();
//    delete file;
// }
//  
// /*
// class Histo {
//    char tag[20];
//    char ID[80];
// public:
//    Histo(char* atag=0) {strcpy(tag,atag);}
//    TH1F* h1d(char* id, Int_t nx, Float_t xmin, Float_t xmax, char* title, Int_t color=33)
//    {
//       strcpy(ID,id);
//       strcat(ID,tag);
//       if (gDirectory->Get(ID)) gDirectory->Delete(ID);
//       TH1F *h = new TH1F(ID, title, nx, xmin, xmax);
//       h->SetFillColor(color);
//       return h;
//    }
//    Bool_t Fill(char* id, Float_t val) {
//       strcpy(ID,id);
//       TH1F* h = (TH1F*) gDirectory->Get(strcat(ID,tag));
//       if (h) { h->Fill(val); return kTRUE;}
//       else return kFALSE;
//    }
//    Bool_t Fill(char* id, Int_t val) {
//       strcpy(ID,id);
//       TH1F* h = (TH1F*) gDirectory->Get(strcat(ID,tag));
//       if (h) { h->Fill(val); return kTRUE;}
//       else return kFALSE;
//    }
// 
//    TH2F* h2d(char* id, Int_t nx, Float_t xmin, Float_t xmax,
//       Int_t ny, Float_t ymin, Float_t ymax, char* title)
//    {
//       strcpy(ID,id);
//       strcat(ID,tag);
//       if (gDirectory->Get(ID)) gDirectory->Delete(ID);
//       TH2F *h = new TH2F(ID, title, nx, xmin, xmax, ny, ymin, ymax);
//       return h;
//    }
//    Bool_t Fill(char* id, Float_t x, Float_t y) {
//       strcpy(ID,id);
//       TH2F* h = (TH2F*) gDirectory->Get(strcat(ID,tag));
//       if (h) { h->Fill(x,y); return kTRUE;}
//       else return kFALSE;
//    }
// };
// 
//    template <class T>
// void triggerNames(T& t)
// {
//    //if (!t.fChain->GetBranch("EVENT")) t.fChain->SetBranchStatus("EVENT",1);
//    //cout << "EVENT branch address: "<< t.fChain->GetBranch("EVENT") <<endl;
//    t.fChain->SetBranchStatus("*",0);
//    t.fChain->SetBranchStatus("EVENT",1);
// 
//    if (t.fChain == 0) return;
//    Int_t nentries = Int_t(t.fChain->GetEntries());
//    Int_t nbytes = 0, nb = 0;
// 
//    struct TnameStat {
//       struct Tname {
//       //const Int_t dim;
//       char name[32];
//       Int_t count;
//       Int_t l3bit;
//       char l3bit_str[80];
//       Tname() {name[0]='\0'; count=0; l3bit=-1;}
//       //Tname():dim(80) {name[0]='\0'; count=0; l3bit=-1;}   // ROOT doesn't support
//       } tname[128];
//       Int_t n;
//       TnameStat() {n=0;}
//       Int_t GetCount() {return n;}
//       Int_t Add(char* name) {
//       for (Int_t i=0; i<n; i++) {
//          if (strcmp(tname[i].name, name) == 0) {
//             tname[i].count++;
//             return 1;
//          }
//       }
//       strcpy(tname[n].name, name);
//       n++;
//       return 0;
//       }
//       Int_t Add(char* name, Int_t l3bit) {
//       char buf[32];
//       for (Int_t i=0; i<n; i++) {
//          if (strcmp(tname[i].name, name) == 0) {
//             tname[i].count++;
//             // check l3bit
//             if (tname[i].l3bit != l3bit) {
//                //cout<<"l3bit for "<<tname[i].name<<" "<<tname[i].l3bit<<" --> "<<l3bit<<endl;
//                sprintf(buf,",%d",l3bit);
//                if (strlen(tname[i].l3bit_str)+strlen(buf) < 80-2) strcat(tname[i].l3bit_str, buf);
//                else tname[i].l3bit_str[strlen(tname[i].l3bit_str)-1] = '*';
//                tname[i].l3bit = l3bit;
//             }
//             return 1;
//          }
//       }
//       strcpy(tname[n].name, name);
//       tname[n].l3bit = l3bit;
//       sprintf(buf,"%d",l3bit);
//       strcpy(tname[n].l3bit_str, buf);
//       n++;
//       return 0;
//       }
//       Int_t Show() {
//       Int_t tot = 0;
//       if (n == 0) cout<<"No fired triggers" <<endl;
//       else cout<<"Fired triggers:" <<endl;
//       for (Int_t i=0; i<n; i++) {
//          tot += tname[i].count;
//          cout<<setw(3)<< i <<" "<<setw(20)<< tname[i].name
//             <<" fired "<<setw(6)<< tname[i].count <<" l3bit: "<<tname[i].l3bit_str<< endl;
//       }
//       cout << "Total " << tot << endl;
//       }
//    } tnameStat;
// 
//    char tname[32];
//    for (Int_t jentry=0; jentry<nentries; jentry++)
//    {
//       Int_t ientry = t.LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
//       nb = t.fChain->GetEntry(jentry);   nbytes += nb;
// 
//       Int_t npro = jentry+1;
//       if (npro%100 == 0) cout<<"\rEvent" <<setw(6)<< npro <<flush;
// 
//       for (Int_t l3bit=0; l3bit<128; l3bit++) {
//       if (t.EVENT_L3bits[l3bit]==1 || t.EVENT_L3MarkPass[l3bit]==1)
//       {
//          Int_t ichar = 0;
//          Int_t iint = l3bit*32;
//          while (tname[ichar++] = (char) t.EVENT_L3names[iint++]);
// 
//          tnameStat.Add(tname,l3bit);
//       }
//       }
//       //wait("next event");
//    }
//    cout << endl;     // CR after counter
//    tnameStat.Show();
// }
// */

#endif  // ifndef utils_cxx
