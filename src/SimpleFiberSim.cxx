#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TFile.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TRandom2.h>
#include <TString.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TMath.h>
#include <TThread.h>

#include "TRandom3.h"
#include "include/SurfaceContainer.hh"
#include "include/Trajectory.hh"
#include "include/SimpleFiberSim.hh"

using namespace std;
using namespace TMath;

/*
double PSTIndex(double lambda)
{
  double x = pow(lambda, -2);
  return 1.5718 + 8412*x + 2.35e8*x*x;
}

double PMMAIndex(double lambda)
{
  double x = pow(lambda, -2);
  return 1.48024 + 4131*x + 3.86e7*x*x;
}
*/

double PSTIndex(double lambda)
{
  double x = pow(lambda, -2);
  return 1.56 + 11184*x - 2.3e8*x*x;
}

double PMMAIndex(double lambda)
{
  double x = pow(lambda, -2);
  return 1.47 + 4544*x - 5.8e7*x*x;
}

double FPIndex(double lambda)
{
  double x = pow(lambda, -2);
  return 1.40 + 4549*x - 6.2e7*x*x;
}

double OilIndex(double lambda)
{
  double x = pow(lambda, -2);
  return 1.45689 + 4362.*x;
}

//implement Fresnel equations
double SimpleFiberSim::Reflectivity(double idx1, double idx2, double theta, bool isS)
{
  double det = 1.0 - pow(idx1*sin(theta)/idx2,2);
  if (det < 0.0) return 1.0;

  if (isS) {
    double Rs_num = idx1*cos(theta) - idx2*sqrt(1.0 - pow(idx1*sin(theta)/idx2,2));
    double Rs_den = idx1*cos(theta) + idx2*sqrt(1.0 - pow(idx1*sin(theta)/idx2,2));
    double Rs = pow( Rs_num/Rs_den, 2 );  
    return Rs;
  }
  else {
    double Rp_num = idx1*sqrt(1.0 - pow(idx1*sin(theta)/idx2,2)) - idx2*cos(theta);
    double Rp_den = idx1*sqrt(1.0 - pow(idx1*sin(theta)/idx2,2)) + idx2*cos(theta);
    double Rp = pow( Rp_num/Rp_den, 2 );
    return Rp;
  }
  //double R = 0.5*(Rs + Rp);
  //return R;
}

double SimpleFiberSim::DrawWLS(double startWavelength)
{
  //cout << "In DrawWLS" << endl;
  double tmpWavelength = -999;
  while (tmpWavelength < startWavelength)
    {
      tmpWavelength = _emission->GetRandom();
    }
  return tmpWavelength;
}

SimpleFiberSim::SimpleFiberSim(double zCut)
  : _container("data/wallReflectivity.dat")
{
  _zCut = zCut;
  //_wlsTimeConstant = 7.4; //ns
  _wlsTimeConstant = 11.8; //ns
  _gen.SetSeed(0);

  _hDtDz_Dz = new TH2D("hDtDz_Dz", ";#DeltaZ;#DeltaT",   350, 0, 3500, 200, 0, 1);
  _hBeta_Dz = new TH2D("hBeta_Dz", ";#DeltaZ;v/c",       350, 0, 3500, 200, 0, 0.63);
  _hvEff_Dz = new TH2D("hvEff_Dz", ";#DeltaZ;v_{eff}",   350, 0, 3500, 200, 0.0, 20.5);
  
  double D = 0.07;
  double claddingPercent = 0.02;
  double thickness = claddingPercent*D;
  double rCore = 0.5*D - 2*thickness;
  double rInner = rCore + thickness;
  double rOuter = rInner + thickness;

  _container.AddFiber( 0.0, 0.0, rOuter  );
  _container.AddFiber( 0.0, 0.0, rInner );
  _container.AddFiber( 0.0, 0.0, rCore );

  TGraph specGraph("data/novaSpectrum.dat");
  specGraph.SetBit(TGraph::kIsSortedX);
  _spectrum = new TH1D("spectrum", "", 100, 360, 520);
  for (int i = 1; i <= _spectrum->GetXaxis()->GetNbins(); ++i)
    {
      double bc = _spectrum->GetXaxis()->GetBinCenter(i);
      double value = specGraph.Eval(bc);
      _spectrum->SetBinContent(i, value);
    }
  _spectrum->SetDirectory(0);

  _absWLS  = new TGraph("data/y11_abs_length.dat");
  _absPS   = new TGraph("data/fiberPSTabsorb.dat");
  _absPMMA = new TGraph("data/PMMABulkAbsorb.dat");
  _absFP   = new TGraph("data/fluorinatedPolymerAbs.dat");

  _absWLS->SetBit(TGraph::kIsSortedX); 
  _absPS->SetBit(TGraph::kIsSortedX); 
  _absPMMA->SetBit(TGraph::kIsSortedX); 
  _absFP->SetBit(TGraph::kIsSortedX); 
  
  TFile* angleFile = new TFile("data/dT_dZ_CollectionRate.root");
  _hTheta_Phi = (TH2D*)angleFile->Get("Theta_Phi");
  _hTheta_Phi->SetDirectory(0);
  angleFile->Close();
  
  TGraph emissionGraph("data/y11Emission.dat");
  _emission = new TH1D("emission", "", 100, 420, 636);
  for (int i = 1; i <= _emission->GetXaxis()->GetNbins(); ++i)
    {
      double bc = _emission->GetXaxis()->GetBinCenter(i);
      double value = emissionGraph.Eval(bc);
      _emission->SetBinContent(i, value);
    }
  _emission->SetDirectory(0);

  _wlsIterations = 0;

  /*
  _tree = new TTree("TrajCoord", "Tree with trajectory coordinates");
  _tree->Branch("x", &_x);
  _tree->Branch("y", &_y);
  _tree->Branch("z", &_z);
  _tree->Branch("wls", &_wlsIterations);
  _tree->Branch("dT", &_dT);
  _tree->Branch("dTdZ", &_dTdZ);
  _tree->Branch("vEff", &_vEff);
  _tree->Branch("beta", &_beta);
  */
}

void SimpleFiberSim::GeneratePhoton( double x, double y, double z, double t, double wavelength)
{
  double rx(-1),ry(-1),rz(-1); //,angle(0.0);
  _gen.Sphere(rx, ry, rz, 1.0);

  double u = _gen.Rndm();
  bool isS = (u < 0.5);
  
  _photon.Init(x, y, z, rx, ry, fabs(rz), t, wavelength, isS); //some default wavelength      
}

bool SimpleFiberSim::TraceOnePhoton( double z )
{
  double currIdx(0);

  _x.clear();
  _y.clear();
  _z.clear();
  
  //cout << "-----------------New Photon---------------------" << endl;

  double D = 0.07;
  double claddingPercent = 0.02;
  double thickness = claddingPercent*D;
  double rCore = D*0.5 - 2*thickness;
  double rInner = rCore + thickness;
  double rOuter = rInner + thickness;  
  double eps = 1e-10;
    
  //double rStart = _gen.Uniform(0, rOuter);
  double x(0), y(0), r(0);
  int iterations = 0;
  _wlsIterations = 0;

  //Assume a scintillation ray has made it just under the outer surface
  //without loss of generality, assume the ray contacts at (x, y, z) = (rOuter - epsilon, 0, 0)
  //ray must have a direction pointing inwards in x, but y and z directions are unconstrained
  //this is the ending condition for the cell simulation
  _photon.R0(rOuter - 1e-8, 0, 0);
  double rx(-1),ry(-1),rz(-1);
  double theta(0), phi(0);
  _hTheta_Phi->GetRandom2(phi, theta);
  rx = sin(theta)*cos(phi);
  ry = sin(theta)*sin(phi);
  rz = cos(theta);
  _photon.RHat(rx, ry, rz);
  //randomly pick S or P polarization
  bool isS = (_gen.Rndm() < 0.5);
  _photon.IsS(isS);
  //draw a wavelength from the scintillator spectrum
  double scintWavelength = _spectrum->GetRandom();
  _photon.Wavelength(scintWavelength);
  _photon.Time(0.0);
  //--------------------------------------------------------------------------
  _x.push_back(_photon.GetX0());
  _y.push_back(_photon.GetY0());
  _z.push_back(_photon.GetZ0());
   
  //int iterations(0);
  while ( true )
    {
      //if (iterations > 100000) return false;
      ++iterations;
      double u(0.0);

      x = _photon.GetX0();
      y = _photon.GetY0();
      r = sqrt( x*x + y*y );
      
      //double cellMediumIdx = 1.46;
      double stepToAbs = 9999.0;
      double stepToWLS = 9999.0;

      double absLength = -9999.0;
      if (r < rCore)
	{
	  //currIdx = 1.59;
	  currIdx = PSTIndex(_photon.GetWavelength());
	  absLength = 100*_absPS->Eval(_photon.GetWavelength());
	  stepToAbs = _gen.Exp(absLength);
	  //wavelength shifting only happens in the core
	  double wlsLength = 100*_absWLS->Eval(_photon.GetWavelength());
	  stepToWLS = _gen.Exp(wlsLength);
	}
      else if (r < rInner)
	{
	  //currIdx = 1.49;
	  currIdx = PMMAIndex(_photon.GetWavelength());
	  absLength = 100*_absPMMA->Eval(_photon.GetWavelength());
	  stepToAbs = _gen.Exp(absLength);
	}
      else if (r < rOuter)
	{
	  //currIdx = 1.42;
	  currIdx = FPIndex(_photon.GetWavelength());
	  absLength = 100*_absFP->Eval(_photon.GetWavelength());
	  stepToAbs = _gen.Exp(absLength);
	}
      else
	{
	  //currIdx = cellMediumIdx;
	  currIdx = OilIndex(_photon.GetWavelength());
	}

      //turn off absorption, do weighting instead
      //stepToAbs = 9999.0;
      
      double stepToFiber(-9999.0);
      double angleToFiber(-9999.0);
      bool hitsFiber = _container.BestFiberIntersection( _photon, stepToFiber, angleToFiber, _fiberNx, _fiberNy, _fiberNz );

      //cout << "stepToAbs = " << stepToAbs
      //	   << " stepToWLS = " << stepToWLS
      //   << " stepToFiber = " << stepToFiber
      //   << endl;
      
      if (hitsFiber)
	{
	  //which process happens first?
	  //absorption happens first -> kill the photon
	  if ( (stepToAbs <= stepToWLS) && (stepToAbs <= stepToFiber) )
	    {
	      //cout << "Absorption - kill photon" << endl;
	      return false;
	    }
	  //wls happens first.  Throw a shifted photon
	  if ( (stepToWLS <= stepToAbs) && (stepToWLS <= stepToFiber) )
	    {
	      u = _gen.Uniform();
	      //QE for k27 is 70%
	      if (u > 0.7) return false;

	      //scale the weight to account for 70% QE for k27
	      //you lose 30% of the weight due to the chance of non-radiative relaxation of k27
	      //_photon.ScaleWeight( 0.7 );

	      //new photon would be going backwards
	      /*
	      if (u < 0.5)
		{
		  //cout << "Backward WLS - kill photon" << endl;
		  return false;
		}
	      */
	      //cout << "Forward WLS - continue" << endl;
	      //new random direction, new wavelength, and additional time
	      double nextWavelength = DrawWLS(_photon.GetWavelength());
	      double extraTime = _gen.Exp(_wlsTimeConstant);
	      double rx(-1),ry(-1),rz(-1); 
	      _gen.Sphere(rx, ry, rz, 1.0);
	      _photon.Step(stepToWLS, currIdx, extraTime);
	      _photon.Wavelength(nextWavelength);
	      _photon.RHat( rx, ry, rz );
	      _wlsIterations++;
	      continue;
	    }
	  //if we make it to here, stepToFiber was the smallest
	  //cout << "Step to next volume" << endl;
	  
	  double nextIdx(0.0);
	  double x1 = _photon.GetX0() + _photon.GetXHat()*(stepToFiber*(1.0 + 1e-6));
	  double y1 = _photon.GetY0() + _photon.GetYHat()*(stepToFiber*(1.0 + 1e-6));

	  double r1 = sqrt( x1*x1 + y1*y1 );
	  if (r1 < rCore)                            
	    {
	      nextIdx = PSTIndex(_photon.GetWavelength());
	      //nextIdx = 1.59;
	    }
	  else if (r1 < rInner)              
	    {
	      nextIdx = PMMAIndex(_photon.GetWavelength());
	      //nextIdx = 1.49;
	    }
	  else if (r1 < rOuter)
	    {
	      nextIdx = FPIndex(_photon.GetWavelength());
	      //nextIdx = 1.42;
	    }
	  else
	    {
	      //nextIdx = cellMediumIdx;
	      nextIdx = OilIndex(_photon.GetWavelength());
	    }
	  
	  //check for reflection
	  double rTh = Reflectivity(currIdx, nextIdx, angleToFiber, _photon.GetIsS());
	  
	  u = _gen.Uniform();
	  if (u < rTh) 
	    {
	      //do the reflection
	      //_photon.Step( stepToFiber - eps, currIdx);
	      _photon.Step( stepToFiber, currIdx);
	      _photon.Reflect( _fiberNx, _fiberNy, _fiberNz );
	      _photon.Step( eps, currIdx);
	      //still in the same medium, so don't change the idx of refraction
	    }	  
	  else 
	    {
	      //refract
	      //_photon.Step( stepToFiber + eps, currIdx);
	      _photon.Step( stepToFiber, currIdx);
	      _photon.Refract( _fiberNx, _fiberNy, _fiberNz, currIdx, nextIdx );
	      currIdx = nextIdx;
	      _photon.Step( eps, currIdx);
	      //medium changed so set currIdx to nextIdx
	    }
	  //scale the weight to account for absorption probability
	  //_photon.ScaleWeight( exp(-stepToFiber/absLength) );

	  //Next check for exit conditions:
	  //if r > 0.5*0.07, we are not in the fiber any more
	  //if z >= _zCut, it made it to the APD
	  double phoW = _photon.GetWeight();
	  double v0 = 29.9792458/1.59;
	  double phoZ = _photon.GetZ0();
	  double phoT = _photon.GetTime();
	  double expT = fabs(phoZ/v0);
	  //double DtDz = (phoT - expT)/phoZ;
	  //if (phoZ > 0)

	  //cout << "PhoWeight = " << phoW << endl;
	  _x.push_back(_photon.GetX0());
	  _y.push_back(_photon.GetY0());
	  _z.push_back(_photon.GetZ0());
	  _dT = phoT - expT;
	  _dTdZ = _dT/fabs(phoZ);
	  _vEff = fabs(phoZ)/phoT;
	  _beta = _vEff/29.9792458;

	  _hDtDz_Dz->Fill( fabs(phoZ), _dTdZ, phoW);
	  _hBeta_Dz->Fill( fabs(phoZ), _beta, phoW );
	  _hvEff_Dz->Fill( fabs(phoZ), _vEff, phoW );

	  //photon makes it to either end of the fiber
	  if ( fabs(_photon.GetZ0()) >= _zCut) 
	    {
	      //if (_vEff > 18 || _vEff < 10) _tree->Fill();
	      //double veff = fabs(phoZ)/phoT;
	      //if (veff > 18.5 || veff < 5) _tree->Fill();
	      return true;
	    }
	}
      else
	{
	  //cout << "Exited the fiber" << endl;
	  return false;
	}
    }
  
  //if it made it here, exceeded the path length limit
  return false;
}

