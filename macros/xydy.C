/*
root -l xydy.C

Data file xydy.dat:
        this is a comment
 1   1   .5     part of line after 3 numbers is ignored
2   20  2--- this line is correct as well
3   30  3
4   25  3
5   35  4
100
# NB: last three lines: line with spaces " " and two empty lines
   


*/
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <sstream>

using std::cout;     using std::endl;

TGraphErrors* xydy(const char* ifname="xydy.dat");
TGraphErrors* xydy(const char* ifname)
{
   cout<< "Using data from " << ifname <<endl;

   std::ifstream ifile(ifname);

   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return 0;
   }

   const Int_t dim = 1000;
   Float_t x[dim], y[dim], dy[dim], dx[dim];
   Int_t np = 0;

   std::string line;
   while (std::getline(ifile,line))
   {
      std::stringstream ss;  ss.str("");    ss << line;
      ss >> x[np] >> y[np] >> dy[np];
      if (ss.fail()) continue;         // conversion error: this can be a comment line

      // successfully read 3 numbers from the line
      cout<< np << "\t " << x[np] <<" "<< y[np] <<" "<< dy[np] <<endl;
      dx[np] = 0.;
      np++;
   }
   ifile.close();
   cout<< "Read " << np << " data lines" <<endl;

   new TCanvas;

   gPad->SetFillColor(0);
   gPad->SetGrid(1,1);

   TGraphErrors *gr = new TGraphErrors(np,x,y,dx,dy);
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(2);
   gr->SetLineColor(2);
   gr->SetName("gr");
   gr->SetTitle("y vs x;abscissa;ordinate");
   gr->Draw("ap");

   // gPad->SaveAs(Form("%s.png", ifname));

   return gr;
}

/*
NB: inverse quotes ``!
g++ `$ROOTSYS/bin/root-config --cflags --glibs` xydy.C
*/
#if !defined(__CINT__)
#include <TRint.h>         // include TRint.h to build stand-along ROOT program
#include <TROOT.h>
#include <cstdarg>

int main(int argc, char *argv[])                      // w/o int ACLiC warning ~ "no type"
{
   TRint* theApp = new TRint("Rint", 0, 0, 0, 0, 1);  // do not show splash screen

   const char* ifname = "xydy.dat";
   if (argc > 1) ifname = argv[1];

   // xydy(ifname); // NB: CINT will not see the return value

   // work in command line mode
   cout<< argv[0] << ": load macro xydy.C" <<endl;
   gROOT->LoadMacro("xydy.C");                        // add "+" to compile with ACLiC
   if (*ifname != '\0') gROOT->ProcessLine(Form("xydy(\"%s\")", ifname));
   else                 gROOT->ProcessLine("xydy()");

   theApp->Run();
   delete theApp;
}
#endif
