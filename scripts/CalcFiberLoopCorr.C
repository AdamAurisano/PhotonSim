#include <TRandom3.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TMath.h>
#include <TFile.h>

#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace TMath;

bool inTorus(double x, double y, double z, double r, double R)
{
  if (x > 0)
    {
      //double det = pow( x*x + y*y + R*R - r*r, 2) - 4.0*R*R*(x*x + y*y);
      double det = pow( R - sqrt(x*x + y*y),2) + z*z - r*r;
      if (det <= 0.0) return true;
      else return false;
    }
  else
    {
      double det1 = pow(y - R,2) + z*z - r*r;
      double det2 = pow(y + R,2) + z*z - r*r;
      if (det1 <= 0.0 || det2 <= 0.0) return true;
      else return false;
    }
}

void CalcFiberLoopCorr()
{
  //ring diameter = 61.5 mm
  //groove 1.6 mm deep
  //fiber radius = 0.7 mm

  double ringD = 6.15 - 0.16;
  double r = 0.07;  
   
  double R = 0.5*ringD+r;

  TRandom3 gen;
  gen.SetSeed(0);
  
  double xmin = -2*R;
  //double xmin = 0;
  double xmax = (R+r);
  
  double ymin = (-R-r);
  double ymax = (R+r);
  
  double zmin = (-r);
  double zmax = (r);

  double Abox = (ymax-ymin)*(zmax-zmin);
  double Afiber = 2.0*Pi()*r*r;

  TFile* outfile = new TFile("FiberLoopCorr.root", "recreate");
  TH1D* histo = new TH1D("hFiberLoopCorr", ";Distance from start of the bend;Correction Factor", 100, xmin, xmax);
  histo->SetDirectory(outfile);
  histo->SetStats(false);
  histo->GetXaxis()->CenterTitle();
  histo->GetYaxis()->CenterTitle();

  double m = R+r;

  //TH3D* hBox = new TH3D("hBox", "", 50, -m, m, 50, -m, m, 50, -m, m);
  
  int ntot = 10000000;
  /*
  for (int iter = 0; iter < ntot; ++iter)
    {
      double x = gen.Uniform(-m, m);
      double y = gen.Uniform(-m, m);
      double z = gen.Uniform(-m, m);
      if (inTorus(x, y, z, r, R)) hBox->Fill(x, y, z);
    }
  */
  
  for (int iX = 1; iX <= histo->GetXaxis()->GetNbins(); ++iX)
    {
      cout << "iX = " << iX << endl;
      double successes = 0.0;
      double x = histo->GetBinCenter(iX);
      for (int iter = 0; iter < ntot; ++iter)
	{
	  //if (iter % 2000 == 0) cout << "iter = " << iter << endl;
	  double y = gen.Uniform(ymin, ymax);
	  double z = gen.Uniform(zmin, zmax);
	  //assume retaining ring blocks about half the light...
	  if (inTorus(x, y, z, r, R))
	    {
	      //hBox->Fill(x, y, z);
	      successes += 1.0;
	    }
	}
      histo->SetBinContent(iX, (Abox/Afiber)*(successes/ntot));
    }
    
  outfile->cd();
  histo->Write();
  histo->Draw();

  //TCanvas* cBox = new TCanvas();
  //hBox->Draw();
  
}


