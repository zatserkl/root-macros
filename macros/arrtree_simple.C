namespace Tree
{
   Int_t    n[3];
   Float_t  a[3];
   Float_t  x;

   void clear()
   {
      for (int i=0; i<3; ++i) {
	 n[i] = 0;
	 a[i] = 0;
	 x = 0;
      }
   }

   void book(TTree* tree) {
      tree->Branch("n",           &n,           "n[3]/I");
      tree->Branch("a",           &a,           "a[3]/F");
      tree->Branch("x",           &x,           "x/F");
   }
   void connect(TTree* tree)                                // need for event-by-event analysis
   {  
      // connects tree buffers with variables to use for event-by-event analysis
      tree->SetBranchAddress("n",             &n);
      tree->SetBranchAddress("a",             &a);
      tree->SetBranchAddress("x",             &x);
   }
}  // namespace Tree

void arrtree_simple()
{
  const char ofname[] = "test.root";
  
  TFile* ofile = TFile::Open(ofname, "recreate");
  TTree* tree = new TTree("t", "test tree");
  
  Tree::book(tree);
  
  // event 0
  Tree::clear();
    Tree::n[0] = 1;
    Tree::a[0] = 11;
  
    Tree::n[1] = 2;
    Tree::a[1] = 12;
  
    Tree::n[2] = 3;
    Tree::a[2] = 13;
  
    Tree::x = 100;
  tree->Fill();
  
  // event 1
  Tree::clear();
    Tree::n[0] = 1;
    Tree::a[0] = 21;
  
    Tree::n[1] = 2;
    Tree::a[1] = 22;
  
    Tree::n[2] = 0;
    Tree::a[2] = 23;
  
    Tree::x = 100;
  tree->Fill();
  
  ofile->Write();
  
  new TCanvas;
  t->Draw("a");
  
  new TCanvas;
  t->Draw("a","n==2");
  
  new TCanvas;
  t->Draw("a","n==3");
}

/*
.q
root -l arrtree_simple.C
*/

