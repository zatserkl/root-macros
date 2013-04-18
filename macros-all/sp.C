// zatserkl@fnal.gov

#include <TTree.h>
#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TROOT.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iterator>

using std::cout;     using std::endl;

#ifndef sp_namespace_Tree
#define sp_namespace_Tree
//
// namespace Tree
//
namespace Tree
{
   Int_t t1;
   Int_t t2;
   Int_t t3;
   Int_t a1;
   Int_t a2;
   Int_t a3;
   
   void clear()
   {
      t1       = 0;
      t2       = 0;
      t3       = 0;
      a1       = 0;
      a2       = 0;
      a3       = 0;
   }
   
   void book(TTree* tree, Int_t nchan=0) {
      if (nchan > 2) {
         tree->Branch("t3",          &t3,         "t3/I");
         tree->Branch("a3",          &a3,         "a3/I");
      }
      if (nchan > 1) {
         tree->Branch("t2",          &t2,         "t2/I");
         tree->Branch("a2",          &a2,         "a2/I");
      }
      tree->Branch("t1",          &t1,         "t1/I");
      tree->Branch("a1",          &a1,         "a1/I");
   }
   void connect(TTree* tree, Int_t nchan=0)
   {   
      // connects tree buffers with variables to use for event-by-event analysis
      if (nchan > 2) {
         tree->SetBranchAddress("t3",            &t3);
         tree->SetBranchAddress("a3",            &a3);
      }
      if (nchan > 1) {
         tree->SetBranchAddress("t2",            &t2);
         tree->SetBranchAddress("a2",            &a2);
      }
      tree->SetBranchAddress("t1",            &t1);
      tree->SetBranchAddress("a1",            &a1);
   }
}  // namespace Tree
#endif   // #define sp_namespace_Tree

#ifndef sp_C_header
#define sp_C_header

// utils.C stuff
TH1* htemp(const char* name=0, const char* title=0);
const char* nextcan(const char* base="can");
const char* nextname(const char* base="h");

TTree* sp(const char* ifname="data.out", const char* tree_name="t", Bool_t plot=kTRUE, Int_t color=2);
void ampl(TTree* tree=0, bool fit=true);
void time(TTree* tree=0, bool fit=true);
TTree* findtree(const char* name="", const char* title="");
TTree* ampl_tree(TTree* tree, const char* a1="a1", const char* a2="a2", TCut cut_a1="", TCut cut_a2="");

#endif   // #ifndef sp_C_header

#ifndef sp_C_implementation
#define sp_C_implementation

//////////////////////////////////////
//
//    GLOBAL VARIABLES
//
// global pointer to the tree
TTree* t = 0;
//-- limits
// TCut sp_t1("sp_t1","t1>0&&t1<9000");
// TCut sp_t2("sp_t2","t2>0&&t2<9000");
// TCut sp_t3("sp_t3","t3>0&&t3<9000");
// TCut sp_a1("sp_a1","a1>0&&a1<1100");
// TCut sp_a2("sp_a2","a2>0&&a2<1100");
// TCut sp_a3("sp_a3","a3>0&&a3<1100");

//-- run0033 settings
TCut sp_t1("sp_t1","t1>0&&t1<9000");
TCut sp_t2("sp_t2","t2>6000&&t2<9000");
TCut sp_t3("sp_t3","t3>6000&&t3<9000");
TCut sp_a1("sp_a1","a1>0&&a1<1100");
TCut sp_a2("sp_a2","a2>0&&a2<1100");
TCut sp_a3("sp_a3","a3>0&&a3<1100");

//-- run0036 settings
// TCut sp_t1("sp_t1","t1>0&&t1<9000");
// TCut sp_t2("sp_t2","t2>6000&&t2<9000");
// TCut sp_t3("sp_t3","t3>6000&&t3<9000");
// TCut sp_a1("sp_a1","a1>0&&a1<1100");
// TCut sp_a2("sp_a2","a2>0&&a2<1100");
// TCut sp_a3("sp_a3","a3>0&&a3<1100");

// TCut sp_t1("sp_t1","t1>0&&t1<9000");
// TCut sp_t2("sp_t2","t2>4000&&t2<9000");
// TCut sp_t3("sp_t3","t3>0&&t3<9000");
// TCut sp_a1("sp_a1","a1>0&&a1<1100");
// TCut sp_a2("sp_a2","a2>0&&a2<1100");
// TCut sp_a3("sp_a3","a3>0&&a3<1100");

//-- SiPM
// TCut sp_t1("sp_t1","t1>0&&t1<9000");
// TCut sp_t2("sp_t2","t2>4000&&t2<9000");
// TCut sp_t3("sp_t3","t3>6000&&t3<9000");
// TCut sp_a1("sp_a1","a1>0&&a1<1100");
// TCut sp_a2("sp_a2","a2>0&&a2<1100");
// TCut sp_a3("sp_a3","a3>0&&a3<1100");

//-- qbars
// TCut sp_t1("sp_t1","t1>4000&&t1<9000");
// TCut sp_t2("sp_t2","t2>6000&&t2<9000");
// TCut sp_t3("sp_t3","t3>2000&&t3<7000");
// TCut sp_a1("sp_a1","a1>0&&a1<1100");
// TCut sp_a2("sp_a2","a2>0&&a2<1100");
// TCut sp_a3("sp_a3","a3>0&&a3<1100");

//-- TOF
// TCut sp_t1("sp_t1","t1>2500&&t1<3500");
// TCut sp_t2("sp_t2","t2>0&&t2<3800");
// TCut sp_t3("sp_t3","t3>7000&&t3<9000");
// TCut sp_a1("sp_a1","a1>60&&a1<1100");
// TCut sp_a2("sp_a2","a2>60&&a2<1100");
// TCut sp_a3("sp_a3","a3>60&&a3<1100");

//-- SiPM again, MPPC matrices
// TCut sp_t1("sp_t1","t1>0&&t1<6000");
// TCut sp_t2("sp_t2","t2>0&&t2<6000");
// TCut sp_t3("sp_t3","t3>5000&&t3<7000");
// TCut sp_a1("sp_a1","a1>0&&a1<1100");
// TCut sp_a2("sp_a2","a2>0&&a2<1100");
// TCut sp_a3("sp_a3","a3>0&&a3<1100");

// pedestal
//Double_t sp_ped[3] = {38.5, 31.4, 26.5};     // ped.out, channels 0, 2, 4
// Double_t sp_ped[3] = {37.6, 29.8, 26.0};     // ped.out, channels 0, 2, 4. run037.out
Double_t sp_ped[3] = {37.7, 29.8, 26.0};     // ped.out, channels 0, 2, 4. run078.out
//////////////////////////////////////

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
  if (name && strlen(name) > 0) h->SetName(name);
  if (title && strlen(title) > 0) h->SetTitle(title);
  return h;
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
      //cout<< "   current hname: " << ss.str() <<endl;
      found = gDirectory->Get(ss.str().c_str());
   }
   //cout<< "new hname: " << Form("%s_%d",base,icycle) <<endl;
   return Form("%s_%d",base,icycle);
}

TTree* sp(const char* ifname, const char* tree_name, Bool_t plot, Int_t color)
{
   std::ifstream ifile(ifname);
   if (!ifile) {
      cout<< "***ERROR : file not found: " << ifname <<endl;
      return 0;
   }

   std::string line;

   // How many numbers do we have in data line?
   // Look in the first line which starts from at least two numbers.
   Int_t nchan = 0;
   while (std::getline(ifile,line))
   {
      std::stringstream ss;   ss.str("");    ss << line;

      // try to read 6 numbers (3 channels)
      ss >> Tree::t1 >> Tree::t2 >> Tree::t3 >> Tree::a1 >> Tree::a2 >> Tree::a3;
      if (!ss.fail()) {
         nchan = 3;
         break;
      }
      ss.clear();                   // clear stream state after previous unsuccessful attempt
      ss.seekg(ios_base::beg);      // back to the beginning of the stream
      ss >> Tree::t1 >> Tree::t2 >> Tree::a1 >> Tree::a2;
      if (!ss.fail()) {
         nchan = 2;
         break;
      }
      ss.clear();                   // clear stream state after previous unsuccessful attempt
      ss.seekg(ios_base::beg);      // back to the beginning of the stream
      ss >> Tree::t1 >> Tree::a1;
      if (!ss.fail()) {
         nchan = 1;
         break;
      }
   }
   ifile.seekg(ios_base::beg);      // back to the beginning of the file

   if (nchan == 0 || nchan > 3) {
      cout<< "no valid data channels found" <<endl;
      return 0;
   }
   else cout<< "found " << nchan << " channel(s) in the data" <<endl;

   TTree* tree_old = (TTree*) gDirectory->Get(tree_name);
   if (tree_old && strcmp(tree_old->GetTitle(), ifname)==0) delete tree_old;

   std::string tree_title = ifname;
   tree_title.resize(tree_title.find("."));                 // title: ifname w/o ".out"
   TTree* tree = new TTree(tree_name, tree_title.c_str());
   Tree::book(tree, nchan);
   tree->SetMarkerStyle(6);
   tree->SetMarkerColor(color);
   tree->SetLineColor(color);

   //-- assign to the global pointer
   if (std::strcmp(tree_name,"t") == 0) t = tree;

   Tree::clear();

   if (nchan == 1) {
      while (std::getline(ifile,line)) {
         std::stringstream ss;   ss.str("");    ss << line;
         ss >> Tree::t1 >> Tree::a1;
         if (!ss.fail()) tree->Fill();
      }
   }
   if (nchan == 2) {
      while (std::getline(ifile,line)) {
         std::stringstream ss;   ss.str("");    ss << line;
         ss >> Tree::t1 >> Tree::t2 >> Tree::a1 >> Tree::a2;
         if (!ss.fail()) tree->Fill();
      }
   }
   if (nchan == 3) {
      while (std::getline(ifile,line)) {
         std::stringstream ss;   ss.str("");    ss << line;
         ss >> Tree::t1 >> Tree::t2 >> Tree::t3 >> Tree::a1 >> Tree::a2 >> Tree::a3;
         if (!ss.fail()) tree->Fill();
      }
   }

   ifile.close();
   cout<< "Read " << tree->GetEntries() << " entries" <<endl;

   if (plot)
   {
      Int_t width = nchan*320+2*(nchan-1);
      if (nchan == 1) width = 450;
      TCanvas* can = new TCanvas(nextcan(ifname), nextcan(ifname), width, 600);

      TCut sp_t[3];
      TCut sp_a[3];
      sp_t[0] = sp_t1;
      sp_t[1] = sp_t2;
      sp_t[2] = sp_t3;
      sp_a[0] = sp_a1;
      sp_a[1] = sp_a2;
      sp_a[2] = sp_a3;
   
      can->Divide(nchan,2);
      for (int i=0; i<nchan; ++i) {
         can->cd(i+1); tree->Draw(Form("t%d",i+1), sp_t[i]);
         can->cd(nchan+i+1); tree->Draw(Form("a%d",i+1), sp_a[i]);
      }
      can->cd();
   }

   return tree;
}

void ampl(TTree* tree, bool fit)
{
   if (tree == 0) tree = t;
   if (tree == 0) {
      cout<< "Default tree not found: t = " << t <<endl;
      return;
   }
   Int_t nchan = 0;
   if (tree->GetBranch("t3")) nchan = 3;
   else if (tree->GetBranch("t2")) nchan = 2;
   else if (tree->GetBranch("t1")) nchan = 1;

   std::string cname = nextcan(Form("%s_ampl", tree->GetTitle()));
   Int_t width = nchan*400+2*(nchan-1);
   if (nchan == 1) width = 560;
   TCanvas* can = new TCanvas(cname.c_str(), cname.c_str(), width,400);
   can->Divide(nchan,1);
   TCut ampl_cut[3];
   ampl_cut[0] = sp_a1;
   ampl_cut[1] = sp_a2;
   ampl_cut[2] = sp_a3;

   for (int i=0; i<nchan; ++i) {
      can->cd(i+1);
      tree->Draw(Form("a%d",i+1), ampl_cut[i]);
      if (fit) {
         // tree->GetHistogram()->Fit("gaus","L");
         tree->GetHistogram()->Fit("gaus");
         Double_t gmean = tree->GetHistogram()->GetFunction("gaus")->GetParameter(1);
         Double_t gsigma = tree->GetHistogram()->GetFunction("gaus")->GetParameter(2);
         Double_t Npe = ((gmean - sp_ped[i])/gsigma)*((gmean - sp_ped[i])/gsigma);
         cout<< "Npe = " << Npe <<endl;
      }
   }
   can->cd(0);
}

void time(TTree* tree, bool fit)
{
   if (tree == 0) tree = t;
   if (tree == 0) {
      cout<< "Default tree not found: t = " << t <<endl;
      return;
   }
   Int_t nchan = 0;
   if (tree->GetBranch("t3")) nchan = 3;
   else if (tree->GetBranch("t2")) nchan = 2;
   else if (tree->GetBranch("t1")) nchan = 1;

   std::string cname = nextcan(Form("%s_time", tree->GetTitle()));
   Int_t width = nchan*400+2*(nchan-1);
   if (nchan == 1) width = 560;
   TCanvas* can = new TCanvas(cname.c_str(), cname.c_str(), width,400);
   can->Divide(nchan,1);
   TCut time_cut[3];
   time_cut[0] = sp_t1;
   time_cut[1] = sp_t2;
   time_cut[2] = sp_t3;

   for (int i=0; i<3; ++i) {
      can->cd(i+1);
      tree->Draw(Form("t%d",i+1), time_cut[i]);
      if (fit) {
         tree->GetHistogram()->Fit("gaus");
      }
   }
   can->cd(0);
}

/*
.q
root -l
//sp("data-4.9_4.8_4.7-32.6-oldVT120-028.out")
//sp("data-4.9_4.8_4.7-32.6-oldVT120-swap2249_a2_a3.out")
.L findtree.C
findtree()
*/
TTree* findtree(const char* name, const char* title)
{
   TIter next(gDirectory->GetList());
   TTree* tree = 0;

   TList list_cand;

   // fill all trees
   //while ((tree = (TTree*) next())) list_cand.Add(tree);
   TObject* obj;
   while ((obj = next())) {
      if (obj->InheritsFrom("TTree")) list_cand.Add(tree);
   }

   if (strlen(name) > 0) {
      // remove trees with other names
      TIter next_sel(&list_cand);
      while ((tree = (TTree*) next_sel())) {
         if (strcmp(tree->GetName(), name) != 0) list_cand.Remove(tree);
      }
   }

   if (strlen(title) > 0) {
      // remove trees with other titles
      TIter next_sel(&list_cand);
      while ((tree = (TTree*) next_sel())) {
         if (strcmp(tree->GetTitle(), title) != 0) list_cand.Remove(tree);
      }
   }

   if (list_cand.GetEntries() > 1) cout<< "selected last of:" <<endl;
   else                            cout<< "selected" <<endl;
   next.Reset();
   while ((tree = (TTree*) next())) {
      if (tree) cout<< tree->GetName() <<"\t\t "<< tree->GetTitle() <<endl;
   }

   tree = (TTree*) list_cand.Last();
   return tree;
}

//////////////////////////////////////////
//
//    subtree
//
//////////////////////////////////////////

/*
t=sp("data/run033.out","t",true,2)
tt=ampl_tree(t,"a1","a2","a1>500")
*/
TTree* ampl_tree(TTree* tree, const char* a1, const char* a2, TCut cut_a1, TCut cut_a2)
{
   Int_t* a1ptr;
   if (strcmp(a1,"a1")==0) a1ptr = &Tree::a1;
   else if (strcmp(a1,"a2")==0) a1ptr = &Tree::a2;
   else if (strcmp(a1,"a3")==0) a1ptr = &Tree::a3;
   else {
      cout<< "ampl_tree: variable not found: " << a1 <<endl;
      return 0;
   }
   Int_t* a2ptr;
   if (strcmp(a2,"a1")==0) a2ptr = &Tree::a1;
   else if (strcmp(a2,"a2")==0) a2ptr = &Tree::a2;
   else if (strcmp(a2,"a3")==0) a2ptr = &Tree::a3;
   else {
      cout<< "ampl_tree: variable not found: " << a2 <<endl;
      return 0;
   }

   const Int_t nchan = 3;
   Tree::connect(tree, nchan);

   // const Int_t nsum = 10;
   // const Int_t nsum = 20;
   const Int_t nsum = 9;
   // Double_t thres = nsum;
   Double_t thres = 6;

   Double_t xmin,xmax, ymin,ymax;
   xmin = ymin = 0;
   xmax = ymax = 1100;
   Int_t ncx = (Int_t) xmax/nsum;
   Int_t ncy = (Int_t) ymax/nsum;

   std::string hname = nextname(Form("h2_%s%s_%s", a1,a2,tree->GetTitle()));
   TH2F* h2 = new TH2F(hname.c_str(), hname.c_str(), ncx,xmin,xmax, ncy,ymin,ymax);

   //new TCanvas;
   //tree->Draw(Form("%s:%s>>%s", a2,a1, hname.c_str()), cut_a1+cut_a2);
   tree->Draw(Form("%s:%s>>%s", a2,a1, hname.c_str()), cut_a1+cut_a2, "goff");

   // subtract nsum counts from every cell
   for (int ix=1; ix<=h2->GetNbinsX(); ++ix)
      for (int iy=1; iy<=h2->GetNbinsY(); ++iy) {
         Double_t content = h2->GetBinContent(ix,iy);
         // subtract thres
         content = content > thres? content-thres: 0;
         h2->SetBinContent(ix,iy, content);
   }
   h2->SetEntries(h2->Integral());

   // find contour

   TTree* tree_contour = new TTree(Form("%s_%s%s",tree->GetName(),a1,a2), Form("%s_%s%s",tree->GetTitle(),a1,a2));
   Tree::book(tree_contour, nchan);
   tree_contour->SetMarkerStyle(tree->GetMarkerStyle());
   tree_contour->SetMarkerColor(tree->GetMarkerColor());
   tree_contour->SetLineColor(tree->GetLineColor());

   for (Int_t ientry=0; ientry<tree->GetEntries(); ++ientry)
   {
      tree->GetEntry(ientry);
      Int_t global_bin = h2->FindBin(*a1ptr,*a2ptr);
      if (h2->IsBinOverflow(global_bin)) continue;
      if (h2->GetBinContent(global_bin) > 0) tree_contour->Fill();
   }
   cout<< "Event loop done, nsum = " << nsum << " thres = " << thres << " tree_contour->GetEntries() = " << tree_contour->GetEntries() <<endl;

   new TCanvas;
   h2->Draw("lego2");

   return tree_contour;
}

// /*
//    NB: inverse quotes ``!
//    g++ `$ROOTSYS/bin/root-config --cflags --glibs` -o sp sp.C
// 
//    NB: to use ./sp do not load file sp.C in your rootlogon.C
// */
// // to define your main apply #define sp_main
// #if !defined sp_main
// #define sp_main
// 
// #if !defined(__CINT__)
// #include <TRint.h>            // include TRint.h to build stand-along ROOT program
// #include <TROOT.h>
// #include <cstdarg>
// 
// int main(int argc, char *argv[])                      // w/o int ACLiC warning ~ "no type"
// {
//   TRint* theApp = new TRint("Rint", 0, 0, 0, 0, 1);   // do not show splash screen
// 
//   const char* ifname = "";
//   if (argc > 1) ifname = argv[1];
//   const char* tree_name = "t";
//   if (argc > 2) tree_name = argv[2];
//   Bool_t plot = kTRUE;
//   if (argc > 3) plot = std::atoi(argv[3]);
//   Int_t color = 1;
//   if (argc > 4) color = std::atoi(argv[4]);
// 
//   // sp(ifname, tree_name, plot);  //-- NB: global pointer t will NOT be seen by CINT!
// 
//   // work in command line mode
//   cout<< argv[0] << ": load macro sp3.C" <<endl;
//   gROOT->LoadMacro("sp.C");                       // add "+" to compile with ACLiC
//   if (*ifname != '\0')  gROOT->ProcessLine(Form("sp(\"%s\",\"%s\",%d,%d)", ifname,tree_name,plot,color));
//   else                  gROOT->ProcessLine("sp()");
//   theApp->Run();
//   delete theApp;
// }
// #endif   // #if !defined(__CINT__)
// #endif   // #if !defined sp_main

#endif   // sp_C_implementation
