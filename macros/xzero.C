#include <TF1.h>

#include <iostream>
using std::cout;     using std::endl;

Double_t xzero(TF1* f, Double_t xmin=0)
{
   // intersection with x-axis tangent to the leading edge

   Double_t xmaximum = f->GetMaximumX();
   Double_t ymaximum = f->Eval(xmaximum);
   Double_t y_hmax = 0.5*ymaximum;
   Double_t x_hmax = f->GetX(y_hmax, xmin,xmaximum);
   Double_t dfdx = f->Derivative(x_hmax);

   // y = kx + b
   // Double_t b = y_hmax - dfdx*x_hmax;
   // Double_t x = -b / dfdx;

   Double_t x = x_hmax - y_hmax/dfdx;
   return x;
}
