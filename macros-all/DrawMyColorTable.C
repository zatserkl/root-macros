#include <TROOT.h>
#include <TColor.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TText.h>
#include <TBox.h>

#include <TColor.h>
#include <iostream>

using std::cout;     using std::endl;

Int_t Color(Int_t number);
TCanvas* newcan(Int_t width=0, Int_t height=0, const char* name="", const char* title="");

// Int_t Color0(Int_t number)
// {
//    if (number >= 0 and number < 10) return number;
//    if (number == 10) return 30;
//    if (number == 11) return 33;
//    if (number == 12) return 38;
//    if (number == 13) return 40;
//    if (number == 14) return 41;
//    if (number == 15) return 42;
//    if (number == 16) return 43;
//    if (number == 17) return 46;
//    if (number == 18) return TColor::GetColorDark(2);
//    if (number == 19) return TColor::GetColorDark(3);
//    if (number == 20) return TColor::GetColorDark(4);
//    if (number == 21) return TColor::GetColorDark(5);
//    if (number == 22) return TColor::GetColorDark(6);
//    if (number == 23) return TColor::GetColorDark(7);
//    if (number == 24) return TColor::GetColorDark(8);
//    if (number == 25) return TColor::GetColorDark(30);
//    if (number == 26) return TColor::GetColorDark(42);
//    if (number == 27) return TColor::GetColorDark(46);
//    if (number == 28) return TColor::GetColorBright(28);
//    return number;
// }

Int_t Color_test(Int_t number)
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

void DrawMyColorTable()
{
   // Static function to Display Color Table in a pad.

   /////////////////////////////////////////////////////////////////////////
   //                                                                      //
   // TColor                                                               //
   //                                                                      //
   // Color defined by RGB or HLS.                                         //
   // At initialization time, a table of colors is generated. This linked  //
   // list can be accessed from the ROOT object                            //
   // (see TROOT::GetListOfColors()). When a color is defined in the range //
   // of [1,50], two "companion" colors are also defined:                  //
   //    - the dark version (color_index + 100)                            //
   //    - the bright version (color_index + 150)                          //
   // The dark and bright color are used to give 3-D effects when drawing  //
   // various boxes (see TWbox, TPave, TPaveText, TPaveLabel,etc).         //
   //                                                                      //
   // This is the list of currently supported basic colors (here dark and  //
   // bright colors are not shown).                                        //
   //

   Int_t i, j;
   Int_t color;
   Double_t xlow, ylow, xup, yup, hs, ws;
   Double_t x1, y1, x2, y2;
   x1 = y1 = 0;
   x2 = y2 = 20;

   TText *text = new TText(0,0,"");
   text->SetTextFont(61);
   text->SetTextSize(0.07);
   text->SetTextAlign(22);

   TBox *box = new TBox();

   // Draw color table boxes.
   hs = (y2-y1)/Double_t(5);
   ws = (x2-x1)/Double_t(10);

   Int_t icolor = 0;

   cout<< "Default (utils.C) colors" <<endl;

   // new TCanvas;
   newcan(0,0,"", "Default (utils.C) colors");

   gPad->SetFillColor(0);
   gPad->Clear();
   gPad->Range(x1,y1,x2,y2);

   icolor = 0;
   for (i=0;i<10;i++) {
      xlow = x1 + ws*(Double_t(i)+0.1);
      xup  = x1 + ws*(Double_t(i)+0.9);
      for (j=0;j<5;j++) {
         ylow = y1 + hs*(Double_t(j)+0.1);
         yup  = y1 + hs*(Double_t(j)+0.9);
         color = 10*j + i;
         box->SetFillStyle(1001);
         color = icolor <= 30? Color(icolor): 1;                     // icolor <= 30
         // if (icolor <= 30) cout<< "my color \t " << icolor <<endl;   // icolor <= 30
         box->SetFillColor(color);
         box->DrawBox(xlow, ylow, xup, yup);
         box->SetFillStyle(0);
         box->SetLineColor(1);
         box->DrawBox(xlow, ylow, xup, yup);
         if (color == 1) text->SetTextColor(0);
         else            text->SetTextColor(1);
         text->DrawText(0.5*(xlow+xup), 0.5*(ylow+yup), Form("%d",icolor));
         ++icolor;
      }
   }

   cout<< "My colors" <<endl;

   // new TCanvas;
   newcan(0,0,"", "My colors");

   gPad->SetFillColor(0);
   gPad->Clear();
   gPad->Range(x1,y1,x2,y2);

   icolor = 0;
   for (i=0;i<10;i++) {
      xlow = x1 + ws*(Double_t(i)+0.1);
      xup  = x1 + ws*(Double_t(i)+0.9);
      for (j=0;j<5;j++) {
         ylow = y1 + hs*(Double_t(j)+0.1);
         yup  = y1 + hs*(Double_t(j)+0.9);
         color = 10*j + i;
         box->SetFillStyle(1001);
         color = icolor <= 30? Color_test(icolor): 1;                     // icolor <= 30
         // if (icolor <= 30) cout<< "my color \t " << icolor <<endl;   // icolor <= 30
         box->SetFillColor(color);
         box->DrawBox(xlow, ylow, xup, yup);
         box->SetFillStyle(0);
         box->SetLineColor(1);
         box->DrawBox(xlow, ylow, xup, yup);
         if (color == 1) text->SetTextColor(0);
         else            text->SetTextColor(1);
         text->DrawText(0.5*(xlow+xup), 0.5*(ylow+yup), Form("%d",icolor));
         ++icolor;
      }
   }
}
