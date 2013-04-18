#include <TROOT.h>
#include <TTree.h>

void ali()
{
   TTree* pulse = (TTree*) gDirectory->Get("pulse");
   if (pulse) {
      // cout<< "Set aliases for the tree pulse" <<endl;

      pulse->SetAlias("t1", "b1_t");
      pulse->SetAlias("t2", "b2_t");

      pulse->SetAlias("c1", "b1_c1");
      pulse->SetAlias("c2", "b1_c2");
      pulse->SetAlias("c3", "b1_c3");
      pulse->SetAlias("c4", "b1_c4");
      pulse->SetAlias("c5", "b2_c1");
      pulse->SetAlias("c6", "b2_c2");
      pulse->SetAlias("c7", "b2_c3");
      pulse->SetAlias("c8", "b2_c4");
   }
}
