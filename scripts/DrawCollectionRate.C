{
  gROOT->ForceStyle();
  
  TFile* file = new TFile("dT_dZ_CollectionRate.root");
  TH2D*  histo = (TH2D*)file->Get("dT_dZ_CollectionRate");
  
  TCanvas* canv = new TCanvas("cdT_dZ_CollectionRate","cdT_dZ_CollectionRate");
  CenterTitles(histo);
  histo->SetStats(false);
  histo->SetTitle(";|Z^{collection} - Z^{emission}| (cm);T^{collection} - T^{emission} (ns);");
  histo->SetMinimum(1.0e-9);
  canv->SetLogz();
  histo->Draw("colz");

  
  Simulation();

}
