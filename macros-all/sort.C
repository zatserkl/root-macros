/*
root -l -q sort.C
*/
#include <TMath.h>
#include <iostream>

using std::cout;        using std::endl;

void sort()
{
   // qsort algorithm?

   Float_t x[5] = {2,  1,  5,  4,  3};
   Float_t y[5] = {20, 10, 50, 40, 30};
   
   Long64_t np = 5;
   // Long64_t index[np]; // ROOT: doesn't work with interpretor because np isn't const
   Long64_t* index = new Long64_t[np];
   
   TMath::Sort(np, x, index, kFALSE);  // kFALSE to sort in ascending order
   
   for (int i=0; i<np; ++i) {
      cout<< "index[i] = " << index[i] << " x[index[i]] = " << x[index[i]] << " y[index[i]] = " << y[index[i]] <<endl;
   }

   delete[] index;
}
