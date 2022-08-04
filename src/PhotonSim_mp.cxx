#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <unistd.h>
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

const int nThreads = 10;
int      threadIdx[nThreads];
Long64_t iterations[nThreads];
TH1D*    hCollectionZ[nThreads];
TH1D*    hCollectionT[nThreads];
TH2D*    hCollectionRate[nThreads];
TH2D*    hCollectionZ_NBounces[nThreads];
TH1D*    hFiberAngle[nThreads];
TH1D*    hNBounces[nThreads];
TH1D*    hPhi[nThreads];
TH1D*    hTheta[nThreads];
TH2D*    hTheta_Phi[nThreads];

void PhotonSim_mp();

int main(int argc, char* argv[])
  {
    PhotonSim_mp();
    return 1;
  }


void* RunOneTrace(void* arg)
{
  TThread::Lock();

  sleep(1);
  double z = 750.0;
  TRandom2 generator(0);

  //double fiberSeparation(6.15);
  double fiberSeparation(5.0);
  PhotonTracer tracer(fiberSeparation);

  int tID = *(int*)arg;

  //cout << "tID = " << tID << endl;

  vector<Double_t> Zbins = MakeBins(false);
  vector<Double_t> ZbinsAbs = MakeBins(true);

  TString hName = "dZ_CollectionRate";
  hName += tID;
  hCollectionZ[tID]    = new TH1D(hName, ";#Delta Z;Collection Rate", 500, -500, 500);
  hCollectionZ[tID]->SetDirectory(0);

  hName = "dT_CollectionRate";
  hName += tID;
  hCollectionT[tID]    = new TH1D(hName, ";#Delta T;Collection Rate", 100, 0, 100);
  hCollectionT[tID]->SetDirectory(0);

  hName = "dT_dZ_CollectionRate";
  hName += tID;
  //hCollectionRate[tID] = new TH2D(hName,";#Delta Z;#Delta T;Collection Rate", ZbinsAbs.size()-1, &ZbinsAbs[0], 30, 0, 60);
  hCollectionRate[tID] = new TH2D(hName,";#Delta Z;#Delta T;Collection Rate", 101, 0, 500, 30, 0, 30);
  hCollectionRate[tID]->SetDirectory(0);

  hName = "CollectionZ_NBounces";
  hName += tID;
  hCollectionZ_NBounces[tID] = new TH2D(hName, ";N Bounces;#DeltaZ", 100, 0, 100, 500, 0, 500);
  hCollectionZ_NBounces[tID]->SetDirectory(0);

  hName = "NBounces";
  hName += tID;
  hNBounces[tID] = new TH1D(hName, ";N Bounces", 100, 0, 100);
  hNBounces[tID]->SetDirectory(0);
  
  hName = "FiberAngle";
  hName += tID;
  hFiberAngle[tID] = new TH1D(hName, ";Angle to Fiber", 500, -0.1, TMath::Pi()+0.1);
  hFiberAngle[tID]->SetDirectory(0);

  hName = "Phi";
  hName += tID;
  hPhi[tID] = new TH1D(hName, ";Phi", 500, -0.1, 2*TMath::Pi()+0.1);
  hPhi[tID]->SetDirectory(0);

  hName = "Theta";
  hName += tID;
  hTheta[tID] = new TH1D(hName, ";Theta", 500, -0.1, TMath::Pi()+0.1);
  hTheta[tID]->SetDirectory(0);

  hName = "Theta_Phi";
  hName += tID;
  hTheta_Phi[tID] = new TH2D(hName, ";Phi;Theta", 500, -0.1, 2*TMath::Pi()+0.1, 500, -0.1, TMath::Pi()+0.1);
  hTheta_Phi[tID]->SetDirectory(0);

  TThread::UnLock();

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

  //iterations[tID] = 1e7;
  iterations[tID] = 1e4;
  //double weight = 1.0/(double)iterations[tID];
  double weight = 1.0;

  cout.setf(ios::left);
  for (Long64_t i = 0; i < iterations[tID]; ++i)
    {
      if (i%2000 == 0 && tID == 0) 
	{
	  double percent_done = 100.0*((double)i/(double)iterations[tID]);
	  cout << "percent done: " << setprecision(2) << fixed << percent_done << "%\r" << flush;
	}
      //if ( tracer.TraceOnePhoton( z, -999, -999, false ) )
      if ( tracer.TraceOnePhoton( z, 3.6/2, 5.64/2, false ) )
	{
	  double phoZ = z - tracer.GetPhoton().GetZ0();
	  double phoT = tracer.GetPhoton().GetTime();

	  double phoX = tracer.GetPhoton().GetX0();
	  double phoY = tracer.GetPhoton().GetY0();

	  double x01 = 2.63518; double y01 = 5.18318;
	  double x02 = 0.934823; double y02 = 0.406816;
	  
	  double phoX01 = phoX - x01;  double phoX02 = phoX - x02;
	  double phoY01 = phoY - y01;  double phoY02 = phoY - y02;
	  double phoR01 = sqrt( phoX01*phoX01 + phoY01*phoY01 );
	  double phoR02 = sqrt( phoX02*phoX02 + phoY02*phoY02 );

	  double phoX0 = (phoR01 < phoR02) ? phoX01 : phoX02;
	  double phoY0 = (phoR01 < phoR02) ? phoY01 : phoY02;
	  
	  double phoXhat = tracer.GetPhoton().GetXHat();
	  double phoYhat = tracer.GetPhoton().GetYHat();
	  double phoZhat = tracer.GetPhoton().GetZHat();
	  
	  double rhat2D = sqrt( phoXhat*phoXhat + phoYhat*phoYhat );
	  double r2D = sqrt( phoX0*phoX0 + phoY0*phoY0 );
	  double phi = atan2(phoY0, phoX0);

	  double phoXRot = r2D*phoX0*cos(phi) + r2D*phoY0*sin(phi);
	  double phoYRot = -r2D*phoX0*sin(phi) + r2D*phoY0*cos(phi);

	  double phoXhatRot = rhat2D*phoXhat*cos(phi) + rhat2D*phoYhat*sin(phi);
	  double phoYhatRot = -rhat2D*phoXhat*sin(phi) + rhat2D*phoYhat*cos(phi);
	  double phiHat = atan2(phoYhatRot, phoXhatRot);

	  while (phiHat > 2*TMath::Pi() || phiHat < 0.0)
	    {
	      if (phiHat < 0.0) phiHat += TMath::TwoPi();
	      if (phiHat > TMath::TwoPi()) phiHat -= TMath::TwoPi();
	    }
	  double thetaHat = acos(phoZhat);

	  //cout << "x = " << phoXRot*1000 << " y = " << phoYRot*1000 << " xhat = " << sin(thetaHat)*cos(phiHat) << " yhat = " << sin(thetaHat)*sin(phiHat) << endl;

	  hTheta_Phi[tID]->Fill( phiHat, thetaHat );
	  hTheta[tID]->Fill(thetaHat);
	  hPhi[tID]->Fill(phiHat);
	  
	  //double phoW = tracer.GetPhoton().GetWeight();
	  //smear phoT by emission time
	  //phoT += generator.Exp(9.0);
	  hCollectionZ[tID]->Fill( phoZ, weight );
	  hCollectionT[tID]->Fill( phoT, weight );
	  hCollectionRate[tID]->Fill( fabs(phoZ), phoT, weight);
	  hCollectionZ_NBounces[tID]->Fill( tracer.GetPhoton().GetNReflections(), fabs(phoZ), weight );
	  hNBounces[tID]->Fill( tracer.GetPhoton().GetNReflections(), weight );
	  hFiberAngle[tID]->Fill( tracer.GetPhoton().GetFiberAngle() );
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

  /*
  if (tID == 0)
    {
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
    }
  */
  return 0;
}

void PhotonSim_mp()
{
  TH1::AddDirectory(false);

  TThread* threads[nThreads];
  for (int iThread = 0; iThread < nThreads; ++iThread)
    {
      TString threadName = "thread_";
      threadName += iThread;
      threadIdx[iThread] = iThread;
      threads[iThread] = new TThread(threadName, RunOneTrace, (void*)&threadIdx[iThread]);
      threads[iThread]->Run();
    }
  TThread::Ps();
  for (int iThread = 0; iThread < nThreads; ++iThread)
    {
      threads[iThread]->Join();
    }

  TH1D* _hCollectionZ = (TH1D*)hCollectionZ[0]->Clone("dZ_CollectionRate");
  _hCollectionZ->Reset();

  TH1D* _hCollectionT = (TH1D*)hCollectionT[0]->Clone("dT_CollectionRate");
  _hCollectionT->Reset();

  TH2D* _hCollectionRate = (TH2D*)hCollectionRate[0]->Clone("dT_dZ_CollectionRate");
  _hCollectionRate->Reset();

  TH2D* _hCollectionZ_NBounces = (TH2D*)hCollectionZ_NBounces[0]->Clone("CollectionZ_NBounces");
  _hCollectionZ_NBounces->Reset();

  TH1D* _hNBounces = (TH1D*)hNBounces[0]->Clone("CollectionZ");
  _hNBounces->Reset();

  TH1D* _hFiberAngle = (TH1D*)hFiberAngle[0]->Clone("FiberAngle");
  _hFiberAngle->Reset();  

  TH1D* _hPhi = (TH1D*)hPhi[0]->Clone("Phi");
  _hPhi->Reset();
  TH1D* _hTheta = (TH1D*)hTheta[0]->Clone("Theta");
  _hTheta->Reset();
  TH2D* _hTheta_Phi = (TH2D*)hTheta_Phi[0]->Clone("Theta_Phi");
  _hTheta_Phi->Reset();    
  
  Long64_t totIter = 0;
  for (int iThread = 0; iThread < nThreads; ++iThread)
    {
      totIter += iterations[iThread];
      _hCollectionZ->Add(hCollectionZ[iThread], 1);
      _hCollectionT->Add(hCollectionT[iThread], 1);
      _hCollectionRate->Add(hCollectionRate[iThread], 1);
      _hCollectionZ_NBounces->Add(hCollectionZ_NBounces[iThread], 1);
      _hNBounces->Add(hNBounces[iThread], 1);
      _hFiberAngle->Add(hFiberAngle[iThread], 1);
      _hPhi->Add(hPhi[iThread], 1);
      _hTheta->Add(hTheta[iThread], 1);
      _hTheta_Phi->Add(hTheta_Phi[iThread], 1);
      cout << "Integral of hCollectionZ[" << iThread << "] = " << hCollectionZ[iThread]->Integral() << endl; 
    }
  cout << "Integral of _hCollectionZ  = " << _hCollectionZ->Integral() << endl; 
  cout << "Total Collection Rate = " << 100.0*(_hCollectionZ->Integral()/totIter) << "%" <<  endl;

  _hCollectionZ->Scale(1.0/totIter);
  _hCollectionT->Scale(1.0/totIter);
  _hCollectionRate->Scale(1.0/totIter);
  _hCollectionZ_NBounces->Scale(1.0/totIter);

  TFile* file = new TFile("dT_dZ_CollectionRate.root","recreate");

  _hCollectionRate->SetDirectory(file);
  file->cd();
  _hCollectionRate->Write();
  _hTheta_Phi->Write();
  
  TCanvas* c1 = new TCanvas("CollectionRate", "CollectionRate");
  c1->cd();
  _hCollectionRate->Draw("colz");

  TCanvas* c2 = new TCanvas("CollectionZ", "CollectionZ");
  c2->cd();
  _hCollectionZ->Draw();


  TCanvas* c3 = new TCanvas("CollectionT", "CollectionT");
  c3->cd();
  _hCollectionT->Draw();

  TCanvas* c4 = new TCanvas("CollectionZ_NBounces","CollectionZ_NBounces");
  c4->cd();
  _hCollectionZ_NBounces->Draw("colz");

  TCanvas* c5 = new TCanvas("FiberAngle", "FiberAngle");
  c5->Divide(2,2);
  c5->cd(1);
  _hFiberAngle->Draw();
  c5->cd(2);
  _hPhi->Draw();
  c5->cd(3);
  _hTheta->Draw();
  c5->cd(4);
  _hTheta_Phi->Draw();

  TCanvas* c6 = new TCanvas("NBounces","NBounces");
  c6->cd();
  _hNBounces->Draw();
}
