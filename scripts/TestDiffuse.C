void TestDiffuse()
{
  TVector3 nhat(1, 1, 0);
  nhat *= 1.0/nhat.Mag();
  
  TRandom3 gen;
  gen.SetSeed(0);

  TH3D* hPoints = new TH3D("hPoints", ";X;Y", 100, -1.5, 1.5, 100, -1.5, 1.5, 100, -1.5, 1.5);
  TH1D* hDotProd = new TH1D("hDotProd", ";rhat*nhat", 300, -1.5, 1.5);
    
  for (int i = 0; i < 100000; ++i)
    {
      //double Th = 0.5*TMath::Pi()*gen.Rndm();
      //double sinTh = sin(Th);
      //double cosTh = cos(Th);

      double cosTh = gen.Uniform(0.5, 1.0);
      double Th = acos(cosTh);
      double sinTh = sin(Th);
      
      //double sinTh = sqrt( gen.Rndm() );
      //double cosTh = sqrt( 1 - pow(sinTh,2) );
      double phi = 2*TMath::Pi()*gen.Rndm();
      
      TVector3 tHat1 = nhat.Orthogonal();
      tHat1 *= 1.0/tHat1.Mag();
      
      TVector3 tHat2 = nhat.Cross(tHat1);
      tHat2 *= 1.0/tHat2.Mag();
      
      TVector3 rhat = cosTh*nhat + sinTh*cos(phi)*tHat1 + sinTh*sin(phi)*tHat2;

      hPoints->Fill( rhat.X(), rhat.Y(), rhat.Z() );
      hDotProd->Fill( rhat*nhat );
    }

  TCanvas* cPoints = new TCanvas("cPoints", "cPoints");
  hPoints->Draw();

  TCanvas* cDotProd = new TCanvas("cDotProd", "cDotProd");
  hDotProd->Draw();
  
}
