void FitAttenuation()
{
  TFile* inFile = new TFile("FiberSpread.root");
  TH2D* hSpread = (TH2D*)inFile->Get("VEff_Dz");
  TH1D* hFrac = hSpread->ProjectionX("hFrac");
  hFrac->Sumw2();
  hFrac->Scale(1.0/hFrac->GetBinContent(1));

  hFrac->Draw("hist");
  
  Double_t par[12];
  TF1* g1 = new TF1("e1", "expo", 0, 20);
  TF1* g2 = new TF1("e2", "expo", 30, 100);
  TF1* g3 = new TF1("e3", "expo", 125, 200);
  TF1* g4 = new TF1("e4", "expo", 225, 600);
  TF1* g5 = new TF1("e4", "expo", 650, 1000);
  TF1* g6 = new TF1("e4", "expo", 2000, 3500);

  TF1* total = new TF1("total", "expo(0)+expo(2)+expo(4)+expo(6)+expo(8)+expo(10)",0,3500);
  
  hFrac->Fit(g1, "R WW");
  hFrac->Fit(g2, "R WW");
  hFrac->Fit(g3, "R WW");
  hFrac->Fit(g4, "R WW");
  hFrac->Fit(g5, "R WW");
  hFrac->Fit(g6, "R WW");
  g1->GetParameters(&par[0]);
  g2->GetParameters(&par[2]);
  g3->GetParameters(&par[4]);
  g4->GetParameters(&par[6]);
  g5->GetParameters(&par[8]);
  g5->GetParameters(&par[10]);
  total->SetParameters(par);

  total->FixParameter(1, par[1]);
  total->FixParameter(3, par[3]);
  total->FixParameter(5, par[5]);
  total->FixParameter(7, par[7]);
  total->FixParameter(9, par[9]);
  total->FixParameter(11, par[11]);
  hFrac->Fit(total, "R WW B");

  total->ReleaseParameter(1);
  total->ReleaseParameter(3);
  total->ReleaseParameter(5);
  total->ReleaseParameter(7);
  total->ReleaseParameter(9);
  total->ReleaseParameter(11);
  hFrac->Fit(total, "");
  
}
