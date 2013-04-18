#include <iostream>

using std::cout;      using std::endl;

void uniona()
{
  union Par {
    int ar[10];
    struct {
      int a1;
      int a2;
      int a3;
    } el;
  };

  Par par;
  par.el.a1 = 1;
  par.el.a2 = 2;
  par.el.a3 = 3;
  cout<< "par.el.a1 = " << par.el.a1 << " par.el.a2 = " << par.el.a2 << " par.el.a3 = " << par.el.a3 <<endl;
  for (int i=0; i<3; ++i) cout<< "par.ar[" << i << "] = " << par.ar[i] <<endl;
}
