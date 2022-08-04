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

void PhotonSim();

int main(int argc, char* argv[])
  {
    PhotonSim();
    return 1;
  }

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


void PhotonSim()
{
  double z = 750.0;
  TRandom2 generator(0);

  //double fiberSeparation = 6.15;
  double fiberSeparation = 5.0;
  PhotonTracer tracer(fiberSeparation);

  vector<Double_t> Zbins = MakeBins(false);
  vector<Double_t> ZbinsAbs = MakeBins(true);

  TFile* file = new TFile("dT_dZ_CollectionRate.root","recreate");
  TH1D* hCollectionZ    = new TH1D("dZ_CollectionRate", ";#Delta Z;Collection Rate", Zbins.size()-1, &Zbins[0]);
  TH1D* hCollectionT    = new TH1D("dT_CollectionRate", ";#Delta T;Collection Rate", 100, 0, 100);
  TH2D* hCollectionRate = new TH2D("dT_dZ_CollectionRate",";#Delta Z;#Delta T;Collection Rate", ZbinsAbs.size()-1, &ZbinsAbs[0], 20, 0, 100);

  double fracCollected(0.0);
  double frac100(0.0);
  double frac150(0.0);
  double frac200(0.0);
  double frac250(0.0);
  double frac300(0.0);
  double frac350(0.0);
  double frac400(0.0);
  double frac450(0.0);
  double frac500(0.0);

  double trapEff = 0.054; //From Kuraray documents

  Long64_t iterations = 1e5;
  double weight = 1.0/(double)iterations;

  cout.setf(ios::left);
  for (Long64_t i = 0; i < iterations; ++i)
    {
      if (i%2000 == 0) 
	{
	  double percent_done = 100.0*((double)i/(double)iterations);
	  cout << "percent done: " << setprecision(2) << fixed << percent_done << "%\r" << flush;
	}
      //if ( tracer.TraceOnePhoton( z, 1.8, 2.82 ) )
      if ( tracer.TraceOnePhoton( z, 0, 0) )
	{
	  double phoZ = z - tracer.GetPhoton().GetZ0();
	  double phoT = tracer.GetPhoton().GetTime();
	  //double phoW = tracer.GetPhoton().GetWeight();
	  //smear phoT by emission time
	  phoT += generator.Exp(9.0);
	  hCollectionZ->Fill( phoZ, weight );
	  hCollectionT->Fill( phoT, weight );
	  hCollectionRate->Fill( fabs(phoZ), phoT, weight);
	  
	  fracCollected += weight;
	  if (fabs(phoZ) > 100) frac100 += weight;
	  if (fabs(phoZ) > 150) frac150 += weight;
	  if (fabs(phoZ) > 200) frac200 += weight;
	  if (fabs(phoZ) > 250) frac250 += weight;
	  if (fabs(phoZ) > 300) frac300 += weight;
	  if (fabs(phoZ) > 350) frac350 += weight;
	  if (fabs(phoZ) > 400) frac400 += weight;
	  if (fabs(phoZ) > 450) frac450 += weight;
	  if (fabs(phoZ) > 500) frac500 += weight;
	}
    }

  cout << "Percent collected = " << setprecision(4) << 100.0*fracCollected << "%" << endl;
  cout << "Percent past 100 = " << setprecision(4) << 100.0*frac100/fracCollected << "%" << endl;
  cout << "Percent past 150 = " << setprecision(4) << 100.0*frac150/fracCollected << "%" << endl;
  cout << "Percent past 200 = " << setprecision(4) << 100.0*frac200/fracCollected << "%" << endl;
  cout << "Percent past 250 = " << setprecision(4) << 100.0*frac250/fracCollected << "%" << endl;
  cout << "Percent past 300 = " << setprecision(4) << 100.0*frac300/fracCollected << "%" << endl;
  cout << "Percent past 350 = " << setprecision(4) << 100.0*frac350/fracCollected << "%" << endl;
  cout << "Percent past 400 = " << setprecision(4) << 100.0*frac400/fracCollected << "%" << endl;
  cout << "Percent past 450 = " << setprecision(4) << 100.0*frac450/fracCollected << "%" << endl;
  cout << "Percent past 500 = " << setprecision(4) << 100.0*frac500/fracCollected << "%" << endl;

  hCollectionRate->SetDirectory(file);
  file->cd();
  hCollectionRate->Write();

  TCanvas* c1 = new TCanvas("CollectionRate", "CollectionRate");
  c1->cd();
  hCollectionRate->Draw("colz");
  c1->Write();
  
  TCanvas* c2 = new TCanvas("CollectionZ", "CollectionZ");
  c2->cd();
  hCollectionZ->Draw();
  c2->Write();

  TCanvas* c3 = new TCanvas("CollectionT", "CollectionT");
  c3->cd();
  hCollectionT->Draw();
  c3->Write();
  
  TCanvas* c4 = new TCanvas("cRefXY", "cRefXY");
  c4->cd();
  tracer.GetRefXY()->Draw("colz");
  c4->Write();
  
  TCanvas* c5 = new TCanvas("cAttenXY", "cAttenXY");
  c5->cd();
  tracer.GetAttenXY()->Draw("colz");
  c5->Write();
  
  TCanvas* c6 = new TCanvas("cAbsXY", "cAbsXY");
  c6->cd();
  tracer.GetAbsXY()->Draw("colz");
  c6->Write();

  file->Close();

  return;
}

