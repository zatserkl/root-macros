#include <iostream>

// doesn't work: #if defined(__MAKECINT__)
//-- #if !defined(__CINT__) || defined(__MAKECINT__)
   #define ACRed  "\e[0;31m"
   #define ACPlain  "\e[0m"
#else
   #define ACRed ""
   #define ACPlain ""
#endif

using std::cout;        using std::endl;

void rootcolor()
{
   cout<<ACRed<< "Color" <<ACPlain<<endl;
}
