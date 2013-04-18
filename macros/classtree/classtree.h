// #include <Rtypes.h>
#include <TObject.h>

#ifndef classtree_h
#define classtree_h

class Disk: public TObject {
public:
   Float_t  vpc;        // Vpc for this disk
   Float_t  rtd[8];     // rtd
   Int_t    bp[8];      // blade/panel for rtd
   Disk(): TObject() {
      // set blades and panels labels
      bp[0]   = 11;
      bp[1]   = 12;
      bp[2]   = 61;
      bp[3]   = 62;
      bp[4]   = 71;
      bp[5]   = 72;
      bp[6]   = 121;
      bp[7]   = 122;
      clear();
   }
   void clear() {
      vpc = 0;
      for (int i=0; i<8; ++i) {
         rtd[i] = 0;
      }
   }
   ClassDef(Disk, 1)
};

class HCyl: public TObject {
public:
   Float_t  vccu;       // Vccu for the whole detector
   Disk     disk[2];    //
   Int_t    d[2];       // disk #
   HCyl(): TObject() {
      // set disk labels
      d[0] = 1;
      d[1] = 2;
      clear();
   }
   void clear() {
      vccu = 0;
      for (int i=0; i<2; ++i) {
         disk[i].clear();
      }
   }
   ClassDef(HCyl, 1)
};

#endif    // ifndef classtree_h
