#include <cmath>
#include <cstdlib>
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
#include <TThread.h>
#include "TStopwatch.h"

#include "RootUtils.hh"
#include "SimpleFiberSim.hh"

using namespace std;
using namespace TMath;

void GetFiberTime(double zCut);

int main(int argc, char* argv[])
{
  if (argc > 1)
    {
      double zCut = atof(argv[1]);
      GetFiberTime(zCut);
      return 0;
    }
  else
    {
      return -1;
    }
}

void GetFiberTime(double zCut)
{
  TRandom2 generator(0);
  
  SimpleFiberSim tracer(zCut);
  double iterations = 1e5;
  double percentDone(0);
  
  //TH1D* hFracVDiff = new TH1D("hFracVDiff", ";FracVDiff", 1000, -1, 3);
  TH1D* hDtDz = new TH1D("hDtDz", ";FracTDiff", 500, -0.02, 1.0);
  TH1D* hWLS  = new TH1D("hWLS", ";N WLS", 20, 0, 20);
  TH1D* hT    = new TH1D("hT",  ";Time (ns)", 1000, 0, 1000);
  double v0 = 29.9792458/1.59;
  
  Long64_t counter = 0;
  Long64_t successes = 0;
  cout.setf(ios::left);
  for (Long64_t i = 0; i < iterations; ++i)
    {
      double z = 0;
      counter++;
      percentDone = 100.0*(double)counter/(double)iterations;
      if (counter%10000 == 0) cout << percentDone << "% Done " << endl;
      if ( tracer.TraceOnePhoton( z ) )
	{
	  z = tracer.GetPhoton().GetZ0();
	  double time = tracer.GetPhoton().GetTime();
	  double expTime = z/v0;
	  //double v = z/time;
	  //hFracVDiff->Fill( (v - v0)/v0 );
	  //hFracTDiff->Fill( (time - expTime)/expTime );
	  hDtDz->Fill( (time - expTime)/z );
	  hWLS->Fill( tracer.GetWLSIterations() );
	  hT->Fill( time );
	  successes++;
	  //cout << "Straightline distance = " << z << " path length = " << tracer.GetPhoton().GetPathLength() << " exp. time = " << z/v0 << " true time = " << time << endl;	  
	}
    }
  cout << "Fraction captured = " << successes/iterations << endl;
  TCanvas* c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);
  hDtDz->Draw();
  c1->cd(2);
  hWLS->Draw();
  c1->cd(3);
  hT->Draw();
  c1->cd(4);
  tracer.GetBeta_Dz()->Draw("colz");
  
  TCanvas* c2 = new TCanvas();
  c2->cd();
  tracer.GetDtDz_Dz()->Draw("colz");
  
}
