#include <iostream>

Double_t tdcu(Int_t rtd, Int_t Vpc=0);

Double_t tdcu(Int_t rtd, Int_t Vpc) {
   static Int_t Vportcard;
   if (Vpc > 0) Vportcard = Vpc;
   if (Vportcard == 0) {
      std::cout<< "Usage: tdcu(Int_t rtd, Int_t Vpc)" <<std::endl;
      return 0;
   }
   Double_t T = (2/3850.e-6) * (1. - Double_t(rtd)/Double_t(Vportcard));
   return T;
}
