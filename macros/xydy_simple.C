#include <iostream>
#include <fstream>

using namespace std;

void xydy_simple(const char* fname="xydy.dat")
{
   float x[1000], y[1000], dy[1000];
   int np = 0;
   
   ifstream ifile(fname);
   if (!ifile.is_open()) {
      cout<< "File not found: " << fname <<endl;
      return;
   }
   
   while (ifile >> x[np] >> y[np] >> dy[np]) {
      cout<< np << " x = " << x[np] << " y = " << y[np] << " dy = " << dy[np] <<endl;
      ++np;
   }
   ifile.close();
}

#ifndef CINT
main()
{
   xydy_simple("xydy.dat");
}
#endif
