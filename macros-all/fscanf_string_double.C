#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>

using std::cout;      using std::endl;

/*
   Content of the file fscanf_string_double.dat:

"Record Length"	10000	"Points"	-2.07000000e-008	-1.87507643e-005
"Sample Interval"	5E-12	s	-2.06950000e-008	-7.62939517e-010
*/

void fscanf_string_double(const char* ifname="fscanf_string_double.dat")
{
  FILE* ifile = fopen(ifname, "r");
  if (!ifile) {
     cout<< "Could not open file " << ifname <<endl;
     return;
  }

  int RecordLength = -1.;
  double SampleInterval = -1.;
  double x[1000], y[1000];
  int np = 0;

  //--
  //-- Without \n I cannot read the next line correctly
  //
  // I think that this is because I'm using explicit fields like "Record Length"
  // NB: you can add/remove spaces in the format string of fscanf even in \"Record Length\"
  //
  //--
  //--ERROR-- fscanf(ifile, "\"Record Length\" %d \"Points\" %lf %lf", &RecordLength, &x[np], &y[np]);
  fscanf(ifile, "\"Record Length\" %d \"Points\" %lf %lf \n", &RecordLength, &x[np], &y[np]);
  np++;
  cout<< "Record Length = " << RecordLength << " Points" <<endl;

  fscanf(ifile, "\"Sample Interval\"	%lf	s	%lf %lf \n", &SampleInterval, &x[np], &y[np]);
  np++;
  fclose(ifile);

  cout<< "SampleInterval = " << SampleInterval <<endl;
  for (int i=0; i<np; ++i) cout<< i <<"\t "<< x[i] <<"\t "<< y[i] <<endl;
}

#ifndef __CINT__
int main()
{
   fscanf_string_double();

   return 0;
}
#endif
