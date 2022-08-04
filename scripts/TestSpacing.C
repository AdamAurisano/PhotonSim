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

vector<Double_t> MakeBins(bool isAbs)
{
  vector<Double_t> bins;

  Double_t highEdge, binWidth;

  if (isAbs) {
    Double_t currEdge = 0.0;
    bins.push_back(currEdge);
    highEdge = 50.0; binWidth = 2.5;
    while( currEdge < highEdge - 1e-3 ) {currEdge += binWidth; bins.push_back(currEdge);}
    highEdge = 100.0; binWidth = 5.0;
    while( currEdge < highEdge - 1e-3 ) {currEdge += binWidth; bins.push_back(currEdge);}
    highEdge = 200.0; binWidth = 10.0;
    while( currEdge < highEdge - 1e-3 ) {currEdge += binWidth; bins.push_back(currEdge);}
    highEdge = 500.0; binWidth = 20.0;
    while( currEdge < highEdge - 1e-3 ) {currEdge += binWidth; bins.push_back(currEdge);}
  }
  else {
    Double_t currEdge = -500.0;
    bins.push_back(currEdge);
    highEdge = -200.0; binWidth = 20.0;
    while( currEdge < highEdge - 1e-3 ) {currEdge += binWidth; bins.push_back(currEdge);}    
    highEdge = -100.0; binWidth = 10.0;
    while( currEdge < highEdge - 1e-3 ) {currEdge += binWidth; bins.push_back(currEdge);}    
    highEdge = -50.0; binWidth = 5.0;
    while( currEdge < highEdge - 1e-3 ) {currEdge += binWidth; bins.push_back(currEdge);}    
    highEdge = 0.0; binWidth = 2.5;
    while( currEdge < highEdge - 1e-3 ) {currEdge += binWidth; bins.push_back(currEdge);}    
    highEdge = 50.0; binWidth = 2.5;
    while( currEdge < highEdge - 1e-3 ) {currEdge += binWidth; bins.push_back(currEdge);}
    highEdge = 100.0; binWidth = 5.0;
    while( currEdge < highEdge - 1e-3 ) {currEdge += binWidth; bins.push_back(currEdge);}
    highEdge = 200.0; binWidth = 10.0;
    while( currEdge < highEdge - 1e-3 ) {currEdge += binWidth; bins.push_back(currEdge);}
    highEdge = 500.0; binWidth = 20.0;
    while( currEdge < highEdge - 1e-3 ) {currEdge += binWidth; bins.push_back(currEdge);}
  }
  
  return bins;
}


void TestSpacing()
{
  double z = 750.0;
  TRandom2 generator(0);

  TH1D* hEffVSpacing = new TH1D("hEffVSpacing", ";Fiber Spacing (cm);Collection Efficiency", 100, 0.1, 6.15);

  Long64_t iterations = 1e6; 
  cout.setf(ios::left);
  for (int i = 1; i <= hEffVSpacing->GetXaxis()->GetNbins(); ++i)
    {
      double spacing = hEffVSpacing->GetXaxis()->GetBinCenter(i);
      PhotonTracer tracer(spacing);      
      Long64_t successes(0);
      for (Long64_t i = 0; i < iterations; ++i)
	{
	  if (i%2000 == 0) 
	    {
	      double percent_done = 100.0*((double)i/(double)iterations);
	      cout << "percent done: " << setprecision(2) << fixed << percent_done << "%\r" << flush;
	    }
	  if ( tracer.TraceOnePhoton( z ) ) successes += 1.0;
	}
      hEffVSpacing->SetBinContent(i, successes/(double)iterations);
    }
  
  hEffVSpacing->Draw();
}

