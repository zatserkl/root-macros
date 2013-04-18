#include <TROOT.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TList.h>

#include <iostream>

using std::cout;     using std::endl;

void lstemp(TCanvas* can=0);
void lstemp(TCanvas* can)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return;
   }
   if (can == 0) can = (TCanvas*) gPad;

   TObjArray pads;
   TObjArray histos;
   TObjArray histos2d;
   TObjArray graphs;
   TObjArray functions;

   TIter next(can->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) {
      if (obj->IsA()->InheritsFrom("TPad")) {
         pads.Add(obj);
         continue;
      }
      if (obj->IsA()->InheritsFrom("TH2")) {
         histos2d.Add(obj);
         continue;
      }
      if (obj->IsA()->InheritsFrom("TH1")) {
         histos.Add(obj);
         continue;
      }
      if (obj->IsA()->InheritsFrom("TGraph")) {
         graphs.Add(obj);
         continue;
      }
      if (obj->IsA()->InheritsFrom("TF1")) {
         functions.Add(obj);
         continue;
      }
   }

   cout<< "Content of the canvas \"" << can->GetName() << "\" \"" << can->GetTitle() << "\"" <<endl;
   for (int i=0; i<pads.GetEntries(); ++i) cout<< "pad\t " << i << " " << pads[i]->GetName() << " " << pads[i]->GetTitle() <<endl;
   for (int i=0; i<histos.GetEntries(); ++i) cout<< "histo\t " << i << " " << histos[i]->GetName() << " " << histos[i]->GetTitle() <<endl;
   for (int i=0; i<histos2d.GetEntries(); ++i) cout<< "histo2d\t " << i << " " << histos2d[i]->GetName() << " " << histos2d[i]->GetTitle() <<endl;
   for (int i=0; i<graphs.GetEntries(); ++i) cout<< "graph\t " << i << " " << graphs[i]->GetName() << " " << graphs[i]->GetTitle() <<endl;
   for (int i=0; i<functions.GetEntries(); ++i) cout<< "function\t " << i << " " << functions[i]->GetName() << " " << functions[i]->GetTitle() <<endl;
}

// --------------------------------------------------------------------
//    gtemp. NB: Declaration of the gtemp should precede htemp2Clone
// --------------------------------------------------------------------

TGraph* gtemp(TCanvas* can, Int_t index);
TGraph* gtemp(TCanvas* can, Int_t index)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   if (can == 0) can = (TCanvas*) gPad;

   TObjArray graphs;
   TIter next(can->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) {
      if (obj->IsA()->InheritsFrom("TGraph")) graphs.Add(obj);
   }
   if (graphs.GetEntries() == 0) {
      cout<< "Could not find TGraph object. ";
      lstemp(can);
      return 0;
   }

   if (index < 0) return (TGraph*) graphs.Last();

   if (index < graphs.GetEntries()) return (TGraph*) graphs[index];
   else return 0;
}

TGraph* gtemp(TCanvas* can=0);
TGraph* gtemp(TCanvas* can) {return gtemp(can, -1);}

TGraph* gtemp(Int_t index);
TGraph* gtemp(Int_t index) {return gtemp(0, index);}

TGraph* gtempClone(TCanvas* can, Int_t index, const char* name="", const char* title="");
TGraph* gtempClone(TCanvas* can, Int_t index, const char* name, const char* title)
{
   TGraph* gr = gtemp(can, index);
   if (gr) {
      gr = (TGraph*) gr->Clone();
      if (gr && name && *name) gr->SetName(name);
      if (gr && title && *title) gr->SetTitle(title);
   }
   return gr;
}

TGraph* gtempClone(TCanvas* can=0, const char* name="", const char* title="");
TGraph* gtempClone(TCanvas* can, const char* name, const char* title) {
   return gtempClone(can, -1, name, title);
}

TGraph* gtempClone(Int_t index, const char* name="", const char* title="");
TGraph* gtempClone(Int_t index, const char* name, const char* title) {
   return gtempClone(0, index, name, title);
}

// --------------------------------------------
//    htemp. NB: htemp2Clone depends on gtemp
// --------------------------------------------

TH1* htemp(TCanvas* can, Int_t index);
TH1* htemp(TCanvas* can, Int_t index)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   if (can == 0) can = (TCanvas*) gPad;

   TObjArray histos;
   TIter next(can->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) {
      if (obj->IsA()->InheritsFrom("TH1")) histos.Add(obj);
   }
   if (histos.GetEntries() == 0) {
      cout<< "Could not find TH1 object. ";
      lstemp(can);
      return 0;
   }

   if (index < 0) return (TH1*) histos.Last();

   if (index < histos.GetEntries()) return (TH1*) histos[index];
   else return 0;
}

TH1* htemp(TCanvas* can=0);
TH1* htemp(TCanvas* can) {return htemp(can, -1);}

TH1* htemp(Int_t index);
TH1* htemp(Int_t index) {return htemp(0, index);}

TH1* htempClone(const char* name, const char* title, TCanvas* can=0);
TH1* htempClone(const char* name, const char* title, TCanvas* can) {
   TH1* h = htemp(can);
   if (!h) return 0;

   if (h && name && *name) h->SetName(name);
   if (h && title && *title) h->SetTitle(title);
   return h;
}

TH1* htempClone(const char* name, const char* title, Int_t index);
TH1* htempClone(const char* name, const char* title, Int_t index) {
   TH1* h = htemp(index);
   if (!h) return 0;

   if (h && name && *name) h->SetName(name);
   if (h && title && *title) h->SetTitle(title);
   return h;
}

TH1* htempClone(const char* name, const char* title, TCanvas* can, Int_t index);
TH1* htempClone(const char* name, const char* title, TCanvas* can, Int_t index) {
   TH1* h = htemp(can,index);
   if (!h) return 0;

   if (h && name && *name) h->SetName(name);
   if (h && title && *title) h->SetTitle(title);
   return h;
}

TH2* htemp2(TCanvas* can, Int_t index);
TH2* htemp2(TCanvas* can, Int_t index)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   if (can == 0) can = (TCanvas*) gPad;

   TObjArray histos2d;
   TIter next(can->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) {
      if (obj->IsA()->InheritsFrom("TH2")) histos2d.Add(obj);
   }
   if (histos2d.GetEntries() == 0) return 0;

   if (index < 0) return (TH2*) histos2d.Last();

   if (index < histos2d.GetEntries()) return (TH2*) histos2d[index];
   else return 0;
}

TH2* htemp2(TCanvas* can=0);
TH2* htemp2(TCanvas* can) {return htemp2(can, -1);}

TH2* htemp2(Int_t index);
TH2* htemp2(Int_t index) {return htemp2(0, index);}

TH1* htemp2Clone(const char* name, const char* title, TCanvas* can=0);
TH1* htemp2Clone(const char* name, const char* title, TCanvas* can) {
   TH1* h = htemp2(can);
   if (!h) return 0;

   if (h && name && *name) h->SetName(name);
   if (h && title && *title) h->SetTitle(title);
   return h;
}

TH1* htemp2Clone(const char* name, const char* title, Int_t index);
TH1* htemp2Clone(const char* name, const char* title, Int_t index) {
   TH1* h = htemp2(index);
   if (!h) return 0;

   if (h && name && *name) h->SetName(name);
   if (h && title && *title) h->SetTitle(title);
   return h;
}

TH1* htemp2Clone(const char* name, const char* title, TCanvas* can, Int_t index);
TH1* htemp2Clone(const char* name, const char* title, TCanvas* can, Int_t index) {
   TH1* h = htemp2(can,index);
   if (!h) return 0;

   if (h && name && *name) h->SetName(name);
   if (h && title && *title) h->SetTitle(title);
   return h;
}

// --------------------------------------------
//    ftemp
// --------------------------------------------

TF1* ftemp(TCanvas* can, Int_t index);
TF1* ftemp(TCanvas* can, Int_t index)
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "No canvas found" <<endl;
      return 0;
   }
   if (can == 0) can = (TCanvas*) gPad;

   TObjArray functions;
   TIter next(can->GetListOfPrimitives());
   TObject* obj = 0;
   while ((obj = next())) {
      if (obj->IsA()->InheritsFrom("TF1")) functions.Add(obj);
   }
   if (functions.GetEntries() == 0) {
      cout<< "Could not find TF1 object. ";
      lstemp(can);
      return 0;
   }

   if (index < 0) return (TF1*) functions.Last();

   if (index < functions.GetEntries()) return (TF1*) functions[index];
   else return 0;
}

TF1* ftemp(TCanvas* can=0);
TF1* ftemp(TCanvas* can) {return ftemp(can, -1);}

TF1* ftemp(Int_t index);
TF1* ftemp(Int_t index) {return ftemp(0, index);}
