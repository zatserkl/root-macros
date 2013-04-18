{
gROOT->Macro(gSystem->ExpandPathName("$(HOME)/macros/rootlogon.C"));
cout<< "*-- Local rootlogon" << endl;

// cout<< "Local rootlogon.C: Defined const char s[] = \"sames\"" <<endl;
// const char s[] = "sames";

// cout<< "Load FWLite" <<endl;
// gSystem->Load("libFWCoreFWLite.so");
// AutoLibraryLoader::enable();
// gSystem->Load("libDataFormatsFWLite.so");
// gROOT->ProcessLine("namespace edm {typedef edm::Wrapper<vector<float> > Wrapper<vector<float,allocator<float> > >; }");
// gROOT->ProcessLine("namespace edm {typedef edm::Wrapper<vector<double> > Wrapper<vector<double,allocator<double> > >; }");

// cout<< "Load eloop.C+" <<endl;
// gROOT->LoadMacro("eloop.C+");
}
