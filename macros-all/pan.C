//namespace root_utils {

//-- remarks

TString rem;

void remark(const char* text = 0);
void remark(const char* text) {
  if (text == 0) {
    if (rem.IsNull()) {
      rem = GetStringDialog("Enter remark", rem.Data());
    }  
  }
  else rem = text;
  //gPad->DrawTextNDC(.01,.01, rem); // not seen on ps
  TText txt;
  txt.DrawTextNDC(.01,.01, rem.Data());
}

void remark_all() {
  if (rem.IsNull()) rem = GetStringDialog("Enter remark", rem.Data());
  TText txt;
  TSeqCollection* clist =  gROOT->GetListOfCanvases();
  TCanvas* canvas_curr = gPad;
  for (int i=0; i<clist->GetSize(); i++) {
    TCanvas* canvas = (TCanvas*) clist->At(i);
    canvas->cd();
    txt.DrawTextNDC(.01,.01, rem.Data());
    canvas->Update();
  }
  canvas_curr->cd();
}

void new_remark() {
  rem = GetStringDialog("Enter remark", rem.Data());
  TText txt;
  txt.DrawTextNDC(.01,.01, rem);
}

void remfile() {
  if (gFile == 0) return;
  TText txt;
  txt.DrawTextNDC(.01,.01, gFile->GetName());
}

void remfile_all() {
  if (gFile == 0) return;
  TText txt;
  TSeqCollection* clist =  gROOT->GetListOfCanvases();
  TCanvas* canvas_curr = gPad;
  for (int i=0; i<clist->GetSize(); i++) {
    TCanvas* canvas = (TCanvas*) clist->At(i);
    canvas->cd();
    txt.DrawTextNDC(.01,.01, gFile->GetName());
    canvas->Update();
  }
  canvas_curr->cd();
}

//-- PostScript

TPostScript* postscript;
TControlBar* ps_bar = 0;
bool ps_opened = false;

// Int_t ps_type = 111;   // portrait  ps
Int_t ps_type = 112;      // landscape ps
// Int_t ps_type = 113;   // eps

void ps_save_selected() {
  const char* psfile = ps_open();
  if (!ps_opened) return;
  ps_bar = new TControlBar("vertical", psfile, 600, 0);
  ps_bar->AddButton("save current canvas (select it using middle click)", "ps_save_current();", "");
    TString mess_close("close ");
    mess_close += psfile;
  ps_bar->AddButton(mess_close, "ps_close()", "");
  ps_bar->Show();
}

void ps_save_current() {
  postscript->NewPage();
  gPad->Modified(1);
  gPad->Update();
}

void ps_save_all() {
  const char *gOpenAsTypes[] = {
    "PostScript",   "*.ps",
    "Macro files",  "*.C",
    "ROOT files",   "*.root",
    "Encapsulated PostScript", "*.eps",
    "Gif files",    "*.gif",
    "All files",    "*",
     0,              0
  };

  static TGFileInfo fi;
  fi.fFileTypes = gOpenAsTypes;
  new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen, &fi);
  
  const char* psfile = fi.fFilename;
  if (psfile == 0) return;
  
  postscript = new TPostScript(psfile, ps_type);
  postscript->Clear();
  
  TSeqCollection* clist =  gROOT->GetListOfCanvases();
  for (int i=0; i<clist->GetSize(); i++) {
    TCanvas* canvas = (TCanvas*) clist->At(i);
    postscript->NewPage();
    canvas->Modified(1);  // mark modified to apply Update()
    canvas->Update();     // this command writes to ps, actually
  }
  postscript->Close();
  cout<< "Output file: " << psfile <<endl;
}

const char* ps_open() {
  const char *gOpenAsTypes[] = {
    "PostScript",   "*.ps",
    "Macro files",  "*.C",
    "ROOT files",   "*.root",
    "Encapsulated PostScript", "*.eps",
    "Gif files",    "*.gif",
    "All files",    "*",
     0,              0
  };

  static TGFileInfo fi;
  fi.fFileTypes = gOpenAsTypes;
  new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen, &fi);
  
  const char* psfile = fi.fFilename;
  if (psfile == 0) return;
  
  postscript = new TPostScript(psfile, ps_type);
  postscript->Clear();
  ps_opened = true;
  return psfile;
}

void ps_close() {
  postscript->Close();
  cout<< "Output file: " << psfile <<endl;
  ps_opened = false;
  delete ps_bar;      ps_bar = 0;
}

// //-- canvas
// 
// TCanvas* SquareCanvas(const char* name="square");
// TCanvas* SquareCanvas(const char* name)
// {
//    TCanvas *squire = new TCanvas(name, "square",0,0,400,400);
//    squire->SetGrid();
//    squire->SetFillColor(0);
//    squire->SetBorderSize(2);
//    squire->SetFrameBorderMode(0);
//    squire->SetFrameBorderMode(0);
//    squire->cd();
//    return squire;
// }
// 
// TCanvas* BigCanvas(char const* name="big");
// TCanvas* BigCanvas(char const* name)
// {
//    char name_def[] = "big";
//    if (!name) name = name_def;
//    TCanvas *big_canv = new TCanvas(name, name, 800,1000);
//    big_canv->SetGrid();
//    big_canv->SetFillColor(0);
//    big_canv->SetBorderSize(2);
//    big_canv->SetFrameBorderMode(0);
//    big_canv->SetFrameBorderMode(0);
//    big_canv->cd();
//    return big_canv;
// }
// 
// TCanvas* WideCanvasSmall(char const* name="widesmall");
// TCanvas* WideCanvasSmall(char const* name)
// {
//    char name_def[] = "widesmall";
//    if (!name) name = name_def;
//    TCanvas *wide_canv = new TCanvas(name, name, 800,600);
//    wide_canv->SetGrid();
//    wide_canv->SetFillColor(0);
//    wide_canv->SetBorderSize(2);
//    wide_canv->SetFrameBorderMode(0);
//    wide_canv->SetFrameBorderMode(0);
//    wide_canv->cd();
//    return wide_canv;
// }
// 
// TCanvas* WideCanvas(char const* name="wide");
// TCanvas* WideCanvas(char const* name)
// {
//    char name_def[] = "wide";
//    if (!name) name = name_def;
//    TCanvas *wide_canv = new TCanvas(name, name, 1000,700);
//    wide_canv->SetGrid();
//    wide_canv->SetFillColor(0);
//    wide_canv->SetBorderSize(2);
//    wide_canv->SetFrameBorderMode(0);
//    wide_canv->SetFrameBorderMode(0);
//    wide_canv->cd();
//    return wide_canv;
// }
// 
// TCanvas* c1x1_old(char const* name="c1x1");
// TCanvas* c1x1_old(char const* name)
// {
//    char name_def[] = "c1x1";
//    if (!name) name = name_def;
//    TCanvas *c1x1_canv = new TCanvas(name, name, 500,350);
//    c1x1_canv->SetGrid(0);
//    c1x1_canv->SetFillColor(0);
//    c1x1_canv->SetBorderSize(2);
//    c1x1_canv->SetFrameBorderMode(0);
//    c1x1_canv->SetFrameBorderMode(0);
//    c1x1_canv->cd();
//    return c1x1_canv;
// }
// 
// TCanvas* c1x2_old(char const* name="c1x2");
// TCanvas* c1x2_old(char const* name)
// {
//    char name_def[] = "c1x2";
//    if (!name) name = name_def;
//    TCanvas *c1x2_canv = new TCanvas(name, name, 500,725);
//    c1x2_canv->SetGrid();
//    c1x2_canv->SetFillColor(0);
//    c1x2_canv->SetBorderSize(2);
//    c1x2_canv->SetFrameBorderMode(0);
//    c1x2_canv->SetFrameBorderMode(0);
//    c1x2_canv->cd();
//    return c1x2_canv;
// }
// 
// TCanvas* bigcan(const char* caname=0, const char* cantit=0) {
//    if (!caname) {
//       cout<< "--> Usage: bigcan(name, title=name)" <<endl;
//       return 0;
//    }
//    if (!cantit) cantit = caname;
//    TCanvas* can = new TCanvas(caname,cantit,420,600);
//    return can;
// }
// 
// TCanvas* c1x2(const char* caname=0, const char* cantit=0) {
//    if (!caname) {
//       cout<< "--> Usage: c1x2(name, title=name)" <<endl;
//       return 0;
//    }
//    if (!cantit) cantit = caname;
//    TCanvas* can = new TCanvas(caname,cantit,420,600);
//    can->Divide(1,2);
//    return can;
// }
// 
// TCanvas* c1x3(const char* caname=0, const char* cantit=0) {
//    if (!caname) {
//       cout<< "--> Usage: c1x3(name, title=name)" <<endl;
//       return 0;
//    }
//    if (!cantit) cantit = caname;
//    TCanvas* can = new TCanvas(caname,cantit,420,600);
//    can->Divide(1,3);
//    return can;
// }

void DrawTextFile() {
  const char *gOpenAsTypes[] = {
    "Text",         "*.txt",
    "Macro files",  "*.C",
    "All files",    "*",
     0,              0
  };

  static TGFileInfo fi;
  fi.fFileTypes = gOpenAsTypes;
  new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen, &fi);
  
  const char* txtfile = fi.fFilename;
  if (txtfile == 0) return;
  
  TPaveText* pavetext = new TPaveText(0.1,0.1,0.9,0.9);
  
//   ifstream fin(txtfile);   
//   char line[256];
//   Int_t len = sizeof(line);
//   while (fin.getline(line,len)) pavetext->AddText(line);
//   fin.close();

  pavetext->DrawFile(txtfile);
  pavetext->SetTextAlign(12);
  pavetext->Draw();
}

void pan()
{
  //gROOT->LoadMacro("$ROOTSYS/tutorials/dialogs.C");
  gROOT->LoadMacro("dialogs.C");
  
  bar = new TControlBar("vertical", "utils", 500, 0);
  bar->AddButton("save all canv",    "ps_save_all();",      "");
  bar->AddButton("save selected",    "ps_save_selected();", "");
  bar->AddButton("remark",      "remark();",                "");
  bar->AddButton("remark all",  "remark_all();",            "");
  bar->AddButton("new remark", "new_remark();",             "");
  bar->AddButton("remfile",     "remfile();",               "");
  bar->AddButton("remfile_all", "remfile_all();",           "");
  bar->AddButton("ct", "ct();", "");         // from $(HOME)/macros/rootalias.C
  bar->AddButton("nul", "nul();", "");       // from $(HOME)/macros/rootalias.C
  bar->AddButton("nedit utils.C", "gROOT->ProcessLine(\".!nedit -background DarkKhaki ~/macros/utils.C ~/macros/rootalias.C &\");", "");
  bar->Show();
}

//}  // namespace root_utils
