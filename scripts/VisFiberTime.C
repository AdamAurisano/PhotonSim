#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

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
#include <TStyle.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TPolyLine3D.h>
#include <TThread.h>

#include "include/RootUtils.hh"
#include "include/SimpleFiberSim.hh"

using namespace std;
using namespace TMath;

void VisFiberTime(double zCut)
{
  gStyle->SetCanvasPreferGL(true);
  TCanvas* c1 = new TCanvas("c1");
  TGeoManager* manager = new TGeoManager("world", "world");
  TGeoMaterial *mat = new TGeoMaterial("Vacuum", 0, 0, 0);
  TGeoMedium   *med = new TGeoMedium("Vacuum", 1, mat);

  double D = 0.07;
  double claddingPercent = 0.02;
  double rCore = D*(0.5 - 2*claddingPercent);
  double rInner = rCore + claddingPercent*D;
  double rOuter = rInner + claddingPercent*D;

  double scale = 10000;
  
  //TGeoVolume *top = manager->MakeTube("Top", med, 0.0, 40, zCut);
  TGeoVolume *top = manager->MakeBox("Top", med, zCut, zCut, zCut);
  manager->SetTopVolume(top);
  TGeoVolume *core  = manager->MakeTube("Core",  med, 0.0*scale,    rCore*scale, zCut);
  TGeoVolume *inner = manager->MakeTube("Inner", med, rCore*scale,  rInner*scale, zCut);
  TGeoVolume *outer = manager->MakeTube("Outer", med, rInner*scale, rOuter*scale, zCut);

  top->AddNode(core,1);
  top->AddNode(inner,1);
  top->AddNode(outer,1);

  manager->CloseGeometry();

  core->SetLineColor(kYellow);
  inner->SetLineColor(kBlue);
  outer->SetLineColor(kRed);

  core->SetTransparency(75);
  inner->SetTransparency(75);
  outer->SetTransparency(75);

  top->Draw("ogl");

  string input;
  bool done = false;

  TRandom gen;
  TPolyLine3D* poly = new TPolyLine3D();
  poly->SetLineWidth(1);
  poly->SetLineColor(kGreen);
  TTimer *timer = new TTimer("gSystem->ProcessEvents();",50,false);

  TFile* inFile = new TFile("FiberSpread_tree.root");
  TTree* tree(0);
  inFile->GetObject("TrajTree",tree);
  vector<double>* x = 0;
  vector<double>* y = 0;
  vector<double>* z = 0;
  int wls = 0;
  double dT = 0;
  double dTdZ = 0;
  double vEff = 0;
  double beta = 0;

  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("z", &z);
  tree->SetBranchAddress("wls", &wls);
  tree->SetBranchAddress("dT", &dT);
  tree->SetBranchAddress("dTdZ", &dTdZ);
  tree->SetBranchAddress("vEff", &vEff);
  tree->SetBranchAddress("beta", &beta);
  //SimpleFiberSim tracer(zCut);

  Long64_t ientry = 0;
  
  while (ientry < tree->GetEntries()) {
    timer->TurnOn();
    timer->Reset();
    cout << "Type <return> for next event, <q> to exit: " << endl;
    getline(cin, input);
    timer->TurnOff();

    if (input == "q") break;

    /*
  label:
   
    while( !tracer.TraceOnePhoton(0.0) ) {};

    double v0 = 29.9792458/1.59;
    double z = fabs(tracer.GetPhoton().GetZ0());
    double time = tracer.GetPhoton().GetTime();
    double expTime = z/v0;
    double vEff = z/time;
    double wls = tracer.GetWLSIterations();
    cout << "z = " << z << " time = " << time << " expTime  = " << expTime << " n wls =  " << wls << " vEff = " << vEff << endl;
    
    if (vEff < 17) goto label;
    //if (time < 100) continue;
    
    //tracer.TraceOnePhoton(0.0);
    
    poly->SetPolyLine(-1);
    int npoints = tracer.GetX().size();
    for (int i = 0; i < npoints; ++i)
      {
	poly->SetPoint(i,
		       1000*tracer.GetX()[i],
		       1000*tracer.GetY()[i],
		       tracer.GetZ()[i]
		       );
      }
    poly->Draw();
    */

    tree->GetEntry(ientry);
    poly->SetPolyLine(-1);
    int npoints = x->size();
    int sign = (z->at(npoints-1) < 0) ? -1 : 1;
    for (int i = 0; i < npoints; ++i)
      {
	poly->SetPoint(i,
		       scale*x->at(i),
		       scale*y->at(i),
		       sign*z->at(i)
		       );
      }
    poly->Draw();

    cout << "N WLS: " << wls << " v effective: " << vEff << " excess time: " << dT << " excess time/distance: " << dTdZ << endl;
    
    c1->Modified();
    c1->Update();
    ++ientry;
  }

}
