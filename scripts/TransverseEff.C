#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TRandom2.h>
#include <TVector3.h>
#include <TString.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TMath.h>

#include "include/RootUtils.hh"
#include "include/PhotonTracer.hh"

using namespace std;
using namespace TMath;

void TransverseEff()
{
  double z = 750.0;
  TRandom2 generator(0);

  //fiberSeparation = 6.15 cm at the retaining ring
  //PhotonTracer tracer(6.15);
  PhotonTracer tracer(5.0);

  TH2D* hTransverseEff = new TH2D("hTransverseEff", ";X (cm);Y (cm)", 100, -0.5, 4.1, 100, -0.5, 6.14);
  int nbins = hTransverseEff->GetXaxis()->GetNbins()*hTransverseEff->GetYaxis()->GetNbins();
  
  Long64_t iterations = 1e4;
  Long64_t totalIter(0);
  cout.setf(ios::left);
  for (int iX = 1; iX <= hTransverseEff->GetXaxis()->GetNbins(); ++iX)
    {
      for (int iY = 1; iY <= hTransverseEff->GetYaxis()->GetNbins(); ++iY)
	{
	  //cout << "Bin (" << iX << ", " << iY << ")" << endl;
	  //cout << endl;
	  double x = hTransverseEff->GetXaxis()->GetBinCenter(iX);
	  double y = hTransverseEff->GetYaxis()->GetBinCenter(iY);
	  
	  Long64_t successes(0);
	  for (Long64_t i = 0; i < iterations; ++i)
	    {
	      totalIter += 1;
	      if (i%2000 == 0) 
		{
		  double percent_done = 100.0*((double)totalIter/(double)(iterations*nbins));
		  cout << "percent done: " << setprecision(2) << fixed << percent_done << "%\r" << flush;
		}
	      if ( tracer.TraceOnePhoton( z, x, y, false ) ) successes += 1.0;
	    }
	  hTransverseEff->SetBinContent(iX, iY, successes/(double)iterations);
	}
      //cout << endl;
    }

  hTransverseEff->Draw("colz");
}

