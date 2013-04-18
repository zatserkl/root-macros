/*
work/root> root -l -q ClassDef_b.C+
*-- Default rootlogon
// 
Processing ClassDef_b.C+...
Info in <TUnixSystem::ACLiC>: creating shared library /srv/zatserkl/work/root/./ClassDef_b_C.so
In file included from /srv/zatserkl/work/root/ClassDef_b_C_ACLiC_dict.h:34,
                 from /srv/zatserkl/work/root/ClassDef_b_C_ACLiC_dict.cxx:17:
/srv/zatserkl/work/root/./ClassDef_b.C: In member function ‘void class_with_b::StreamerNVirtual(TBuffer&)’:
/srv/zatserkl/work/root/./ClassDef_b.C:18: warning: declaration of ‘b’ shadows a member of 'this'
/srv/zatserkl/work/root/./ClassDef_b.C: In member function ‘void BadClass::fun(int)’:
/srv/zatserkl/work/root/./ClassDef_b.C:30: warning: declaration of ‘a’ shadows a member of 'this'
*/

#include <TObject.h>

class class_with_b: public TObject {
public:
   Float_t b;
   ClassDef(class_with_b, 1);
};

void ClassDef_b() {
   class_with_b obj;
   obj.b = 0;
}

// example of the same mistake
class BadClass {
public:
   int a;
   int dummy;
   void fun(int a) {dummy = a;}
   BadClass(): a(0), dummy(0) {}
};

#include <iostream>
using std::cout;        using std::endl;

void useBadClass() {
   BadClass badClass;
   cout<< "badClass.a = " << badClass.a <<endl;
   badClass.fun(1);
   cout<< "badClass.a = " << badClass.a <<endl;
}
