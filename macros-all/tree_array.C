/*
Hi Andriy,

This is the expected behavior (See the TTree::Draw documentation for all
the gory details).

In the case where you introduce the TClonesArray, you now have 2 nested
containers/arrays.

Namely the request:


   to->Draw("x[0]","");

draws the content (all the values) of the first element of the outer array.

To get the same results has with the 'ti' tree, you need to request
the first element of the inner array of all the elements of the outer array:


   to->Draw("x[][0]","");

In less compact form (and possibly clear), the 2 statements above are
respectively equivalent to:

   to->Draw("arr[0].x","");
   to->Draw("arr[].x[0]","");

Cheers,
Philippe.

Hi Philippe,

Thank you for clarification.
I prefer your suggestion

to->Draw("arr[].x[0]","")

Best,
Andriy.
*/
#include <TROOT.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TCanvas.h>

#include <iostream>

using std::cout;     using std::endl;

class ArrClass: public TObject {
public:
   Float_t x[5];
public:
   ArrClass(): TObject() {
      for (int i=0; i<5; ++i) x[i] = 0;
   }

   ClassDef(ArrClass, 1)
};

class OuterClass: public TObject {
public:
   TClonesArray* arr;
public:
   OuterClass(): TObject() {
      arr = new TClonesArray("ArrClass");
   }
   ~OuterClass() {
      delete arr;
   }

   ClassDef(OuterClass, 1)
};

#ifdef __MAKECINT__
#pragma link C++ class ArrClass;
#pragma link C++ class OuterClass;
#endif

ClassImp(ArrClass);
ClassImp(OuterClass);

TTree* tree_inner()
{
   TTree* tree = new TTree("ti", "ArrClass tree");

   ArrClass* arrClass = new ArrClass;
   tree->Branch("arrClass", arrClass);

   for (int ievent=0; ievent<3; ++ievent)
   {
      for (int i=0; i<5; ++i) {
         Float_t val = (ievent+1)*(i+1);
         arrClass->x[i] = val;
         cout<< "arrClass->x[" << i << "] = " << arrClass->x[i] << "   ";
      }
      cout<<endl;
      tree->Fill();
   }

   return tree;
}

TTree* tree_outer()
{
   TTree* tree = new TTree("to", "OuterClass tree");

   OuterClass* outerClass = new OuterClass;
   tree->Branch("outerClass", outerClass);

   for (int ievent=0; ievent<3; ++ievent)
   {
      ArrClass* arrClass = new ((*outerClass->arr)[0]) ArrClass;
      for (int i=0; i<5; ++i) {
         Float_t val = (ievent+1)*(i+1);
         arrClass->x[i] = val;
         cout<< "arrClass->x[" << i << "] = " << arrClass->x[i] << "   ";
      }
      cout<<endl;
      tree->Fill();
   }

   return tree;
}

void tree_inner_read(TTree* ti)
{
   ArrClass* arrClass = 0;
   ti->SetBranchAddress("arrClass", &arrClass);

   for (int ientry=0; ientry<ti->GetEntries(); ++ientry) {
      if (ti->LoadTree(ientry) < 0) break;
      ti->GetEntry(ientry);
      cout<< "ientry = " << ientry << "\t ";
      for (int i=0; i<5; ++i) cout<< "   " << arrClass->x[i];
      cout<<endl;
   }
}

void tree_outer_read(TTree* to)
{
   OuterClass* outerClass = 0;
   to->SetBranchAddress("outerClass", &outerClass);

   for (int ientry=0; ientry<to->GetEntries(); ++ientry) {
      if (to->LoadTree(ientry) < 0) break;
      to->GetEntry(ientry);
      cout<< "ientry = " << ientry << "\t ";
      ArrClass* arrClass = (ArrClass*) outerClass->arr->At(0);
      for (int i=0; i<5; ++i) cout<< "   " << arrClass->x[i];
      cout<<endl;
   }
}

void tree_array()
{
   cout<< "-- create inner" <<endl;
   TTree* ti = tree_inner();
   cout<< "-- create outer" <<endl;
   TTree* to = tree_outer();

   cout<<endl<< "-- loop over inner" <<endl;
   tree_inner_read(ti);

   cout<<endl<< "-- loop over outer" <<endl;
   tree_outer_read(to);

   new TCanvas;
   ti->Draw("x[0]","");
   new TCanvas;
   ti->Draw("x[1]","");

   new TCanvas;
   to->Draw("x[0]","");
   new TCanvas;
   to->Draw("x[1]","");
}
