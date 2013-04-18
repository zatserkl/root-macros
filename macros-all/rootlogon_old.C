#include <TROOT.h>
#include <TRint.h>
#include <TSystem.h>
#include <TStyle.h>

#include <iostream>

using std::cout;     using std::endl;

void rootlogon()
{
   cout<< "*-- Default rootlogon" <<endl;

   gSystem->Exec("echo //-- `date` >> root_hist");

   gROOT->LoadMacro("utils.C+");

   /* Add to include path directory with utils.C to use it for ACLiC like
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <utils.C>
#endif
   */
   // gSystem->AddIncludePath(gSystem->ExpandPathName("-I$(HOME)/macros"));
   //cout<< "gSystem->GetIncludePath() = " << gSystem->GetIncludePath() <<endl;

   gROOT->SetStyle("Plain");
   gStyle->SetCanvasColor(0);

   ((TRint*)gROOT->GetApplication())->SetPrompt("// ");

   gStyle->SetOptFit(1);
      
   gStyle->SetMarkerColor(2);
   gStyle->SetFillStyle(3001);

   //-- my long used sizes, 80% of default size
   //gStyle->SetCanvasDefW(540); gStyle->SetCanvasDefH(400);
   //-- 3/4 of internal area of the default: 526x382 --> 522x354 as newcan in utils.C
   gStyle->SetCanvasDefW(526); gStyle->SetCanvasDefH(382);

   gStyle->SetCanvasDefX(0);
   gStyle->SetCanvasDefY(0);
   
   gStyle->SetPadGridX(kTRUE);
   gStyle->SetPadGridY(kTRUE);

   // no box around title
   gStyle->SetTitleX(0);
   gStyle->SetTitleW(0);
   gStyle->SetTitleBorderSize(0);

   // To make title box transparent https://docs.google.com/Doc?docid=0AaltnJcgAgafZGNoOThnaDZfMzEzZjlxOTZ4Zmg&hl=en
   // gStyle->SetStatStyle(0);
   // gStyle->SetTitleStyle(0);
   // gROOT->ForceStyle();

   // gStyle->SetStatStyle(3001);

   // text sizes
   gStyle->SetTitleFontSize(.07);
   gStyle->SetLabelSize(.05,"xyz");
   gStyle->SetTitleSize(.05,"xyz");
   // and fonts
   gStyle->SetTitleFont(42,"xyz");  // if "xyz" set all 3 axes, any other value of axis will set the pad title font
   gStyle->SetTitleFont(42,"a");
   gStyle->SetLabelFont(42,"xyz");
   gStyle->SetLabelFont(42,"a");
   gStyle->SetStatFont(42);
   gStyle->SetTextFont(42);
   
   // gStyle->SetNdivisions(509,"XYZ");
   gStyle->SetNdivisions(508,"XYZ");
}
