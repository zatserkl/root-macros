
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
