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
#include <TThread.h>
#include "TStopwatch.h"

#include "include/RootUtils.hh"
#include "include/SimpleFiberSim.hh"

using namespace std;
using namespace TMath;

const int nThreads = 10;
const double v0 = 29.9792458/1.59;

double   _zCut;
int      threadIdx[nThreads];
Long64_t iterations[nThreads];
Long64_t successes[nThreads];
TH1D*    _hDtDz[nThreads];
TH1D*    _hWLS[nThreads];
TH1D*    _hT[nThreads];
TH1D*    _hvEff[nThreads];
TH2D*    _hDtDz_Dz[nThreads];
TH2D*    _hBeta_Dz[nThreads];
TH2D*    _hvEff_Dz[nThreads];
double   percentDone[nThreads];
bool     isDone[nThreads];

void* GetFiberTime(void* arg)
{
  TThread::Lock();

  sleep(1);
  
  TRandom2 generator(0);
  
  SimpleFiberSim tracer(_zCut);

  int tID = *(int*)arg;

  isDone[tID] = false;
  percentDone[tID] = 0.0;
  iterations[tID] = 1e9;
  //iterations[tID] = 1e5;
  successes[tID] = 0;

  TString hName;

  hName = "hDtDz";
  hName += tID;
  _hDtDz[tID] = new TH1D(hName, ";Dt/z (ns/cm)", 500, -0.02, 1.0);

  hName = "hWLS";
  hName += tID;
  _hWLS[tID]  = new TH1D(hName, ";N WLS", 20, 0, 20);

  hName = "hT";
  hName += tID;
  _hT[tID]    = new TH1D(hName,  ";Time (ns)", 600, -100, 500);

  hName = "hvEff";
  hName += tID;
  _hvEff[tID]    = new TH1D(hName,  ";v_{eff} (cm/ns)", 120, 0, 30);

  TThread::UnLock();
  
  Long64_t counter = 0;
  cout.setf(ios::left);
  for (Long64_t i = 0; i < iterations[tID]; ++i)
    {
      double z = 0;
      counter++;
      percentDone[tID] = 100.0*(double)counter/(double)iterations[tID];
      if ( tracer.TraceOnePhoton( z ) )
	{
	  z = fabs(tracer.GetPhoton().GetZ0());
	  double time = tracer.GetPhoton().GetTime();
	  double weight = tracer.GetPhoton().GetWeight();
	  double expTime = z/v0;
	  double v = z/time;
	  //hFracVDiff->Fill( (v - v0)/v0 );
	  //hFracTDiff->Fill( (time - expTime)/expTime );
	  //_hDtDz[tID]->Fill( (time - expTime)/z );
	  _hDtDz[tID]->Fill( (time - expTime)/expTime, weight );
	  _hWLS[tID]->Fill( tracer.GetWLSIterations(), weight );
	  _hT[tID]->Fill( time - expTime, weight);
	  _hvEff[tID]->Fill( v, weight );
	  successes[tID]++;
	}
    }

  hName = "hDtDz_Dz";
  hName += tID;
  _hDtDz_Dz[tID]    = (TH2D*)tracer.GetDtDz_Dz()->Clone(hName);

  hName = "hBeta_Dz";
  hName += tID;
  _hBeta_Dz[tID]    = (TH2D*)tracer.GetBeta_Dz()->Clone(hName);

  hName = "hvEff_Dz";
  hName += tID;
  _hvEff_Dz[tID]    = (TH2D*)tracer.GetvEff_Dz()->Clone(hName);

  isDone[tID] = true;
  return 0;
}

void GetFiberTime_mp(double zCut)
{
  _zCut = zCut;

  TH1::AddDirectory(false);
  TStopwatch stopwatch;

  stopwatch.Start();
  TThread* threads[nThreads];
  for (int iThread = 0; iThread < nThreads; ++iThread)
    {
      TString threadName = "thread_";
      threadName += iThread;
      threadIdx[iThread] = iThread;
      threads[iThread] = new TThread(threadName, GetFiberTime, (void*)&threadIdx[iThread]);
      threads[iThread]->Run();
    }
  TThread::Ps();

  while (true)
    {
      sleep(5);
      bool allDone = true;
      for (int iThread = 0; iThread < nThreads; ++iThread)
	{
	  allDone = allDone && isDone[iThread];
	  cout << iThread << ": " << setprecision(4) << fixed << percentDone[iThread] << "%     ";
	}
      cout << "\r" << flush;
      if (allDone)
	{
	  cout << "\nAll threads seem to be finished!" << endl;
	  break;
	}
    }

  for (int iThread = 0; iThread < nThreads; ++iThread)
    {
      threads[iThread]->Join();
    }

  stopwatch.Stop();
  
  cout << "Total time to run fiber simulation " << stopwatch.RealTime() << "s" << endl;

  Long64_t iterationsTot(0);
  Long64_t successesTot(0);
  
  TH1D* hDtDz = (TH1D*)_hDtDz[0]->Clone("hDtDz");
  hDtDz->Reset();

  TH1D* hWLS = (TH1D*)_hWLS[0]->Clone("hWLS");
  hWLS->Reset();
  
  TH1D* hT = (TH1D*)_hT[0]->Clone("hT");
  hT->Reset();

  TH1D* hvEff = (TH1D*)_hvEff[0]->Clone("hvEff");
  hvEff->Reset();
  
  TH2D* hDtDz_Dz = (TH2D*)_hDtDz_Dz[0]->Clone("hDtDz_Dz");
  hDtDz_Dz->Reset();

  TH2D* hBeta_Dz = (TH2D*)_hBeta_Dz[0]->Clone("hBeta_Dz");
  hBeta_Dz->Reset();

  TH2D* hvEff_Dz = (TH2D*)_hvEff_Dz[0]->Clone("hvEff_Dz");
  hvEff_Dz->Reset();

  for (int iThread = 0; iThread < nThreads; ++iThread)
    {
      iterationsTot += iterations[iThread];
      successesTot += successes[iThread];
      hDtDz->Add(   _hDtDz[iThread],    1);
      hWLS->Add(    _hWLS[iThread],     1);
      hT->Add(      _hT[iThread],       1);
      hvEff->Add(   _hvEff[iThread],    1);
      hDtDz_Dz->Add(_hDtDz_Dz[iThread], 1);
      hBeta_Dz->Add(_hBeta_Dz[iThread], 1);
      hvEff_Dz->Add(_hvEff_Dz[iThread], 1);
    }
  
  
  cout << "Fraction captured = " << successesTot/(double)iterationsTot << endl;
  TCanvas* c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);
  hDtDz->Draw();
  c1->cd(2);
  hWLS->Draw();
  c1->cd(3);
  hT->Draw();
  c1->cd(4);
  hDtDz_Dz->Draw("colz");

  TCanvas* c2 = new TCanvas("c2", "c2");
  c2->Divide(2,1);
  c2->cd(1);
  hBeta_Dz->Draw("colz");
  c2->cd(2);
  hvEff_Dz->Draw("colz");

  TCanvas* c3 = new TCanvas("c3", "c3");
  hvEff->Draw();
  
  TFile* file = new TFile("FiberSpread.root", "recreate");
  file->WriteTObject(hDtDz_Dz, "FiberSpread");
  file->WriteTObject(hBeta_Dz, "Beta_Dz");
  file->WriteTObject(hvEff_Dz, "VEff_Dz");
}
