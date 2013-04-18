Stat_t hint(char* hname, Float_t x1=0, Float_t x2=0);
Stat_t hint(char* hname, Float_t x1, Float_t x2)
{
  TObject* obj = (TObject*) gDirectory->Get(hname);
  if (obj && strcmp(obj->IsA()->GetName(),TH1F::Class()->GetName())==0)
  {
    TH1F* h = (TH1F*) obj;
    Int_t bin1 = h->GetXaxis()->FindBin(x1);
    if (x2 == 0) x2 = h->GetXaxis()->GetXmax();
    Int_t bin2 = h->GetXaxis()->FindBin(x2);
    Stat_t sum = h->Integral(bin1,bin2);
    Stat_t tot = h->Integral(h->GetXaxis()->GetFirst(),h->GetXaxis()->GetLast());
    //cout<< "Intergal is " << sum << " ==> " << 100.*Float_t(sum)/Float_t(tot) << " %" <<endl;
    cout<< "Intergal from " << x1 << " to " << x2 << " is " << sum << " ==> " << 100.*Float_t(sum)/Float_t(tot) << " %" <<endl;
    return sum;
  }
  else {
    cout<< "Could not found histogram " << hname <<endl;
    return 0;
  }
}

Stat_t hint(TH1F* h, Float_t x1=0, Float_t x2=0);
Stat_t hint(TH1F* h, Float_t x1, Float_t x2)
{
  Int_t bin1 = h->GetXaxis()->FindBin(x1);
  if (x2 == 0) x2 = h->GetXaxis()->GetXmax();
  Int_t bin2 = h->GetXaxis()->FindBin(x2);
  Stat_t sum = h->Integral(bin1,bin2);
  Stat_t tot = h->Integral(h->GetXaxis()->GetFirst(),h->GetXaxis()->GetLast());
  cout<< "Intergal is " << sum << " ==> " << 100.*Float_t(sum)/Float_t(tot) << " %" <<endl;
  return sum;
}
