void DrawBetaDz()
{
  TFile* f = new TFile("FiberSpread.root");
  TH2D* h = (TH2D*)f->Get("Beta_Dz");

  int nX = h->GetNbinsX();
  int nY = h->GetNbinsY();

  for (int iX = 1; iX <= nX; ++iX)
    {
      TH1D* hProj = h->ProjectionY();    
      double norm = hProj->Integral();
      cout << "Column " << iX << " norm = " << norm << endl;
      for (int iY = 1; iY <= nY; ++iY)
	{
	  double val = h->GetBinContent(iX, iY);
	  if (norm > 0)
	    {
	      val /= norm;
	      h->SetBinContent(iX, iY, val);
	    }	  
	}
    }
  h->Draw("colz");
}
