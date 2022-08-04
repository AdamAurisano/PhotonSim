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
#include "include/PhotonTracer.hh"

using namespace std;
using namespace TMath;

//implement Fresnel equations
double PhotonTracer::Reflectivity(double idx1, double idx2, double theta, bool isS)
{
  double det = 1.0 - pow(idx1*sin(theta)/idx2,2);
  if (det < 0.0) return 1.0;

  if (isS) {
    double Rs_num = idx1*cos(theta) - idx2*sqrt(1.0 - pow(idx1*sin(theta)/idx2,2));
    double Rs_den = idx1*cos(theta) + idx2*sqrt(1.0 - pow(idx1*sin(theta)/idx2,2));
    return pow( Rs_num/Rs_den, 2 );
  }
  else {
    double Rp_num = idx1*sqrt(1.0 - pow(idx1*sin(theta)/idx2,2)) - idx2*cos(theta);
    double Rp_den = idx1*sqrt(1.0 - pow(idx1*sin(theta)/idx2,2)) + idx2*cos(theta);
    return pow( Rp_num/Rp_den, 2 );
  }
  //double R = 0.5*(Rs + Rp);
  //return R;
}


//0.877 reflectivity
//beta attenuation was 50 m
PhotonTracer::PhotonTracer(double fiberSeparation) : _container("data/wallReflectivity.dat" ),
						     _beta_atten( 900.0 ), _beta_emission( 9.0 ),
						     _wallReflectivity( "data/wallReflectivity.dat" ),
						     _endcapReflectivity( "data/endReflectivity.dat" ),
						     _width( 3.57 ) , _height( 5.59 ), _length( 1550.0 ),
						     _fiberSeparation( fiberSeparation )
{
  //old: width = 3.57, height = 5.59
  _gen.SetSeed(0);

  _refXY = new TH2D("hRefXY", ";X (cm); Y (cm)", 1000, -0.5*_width, 1.5*_width, 1000, -0.5*_height, 1.5*_height);
  _refXY->SetDirectory(0);

  _attenXY = new TH2D("hAttenXY", ";X (cm); Y (cm)", 1000, -0.5*_width, 1.5*_width, 1000, -0.5*_height, 1.5*_height);
  _attenXY->SetDirectory(0);

  _absXY = new TH2D("hAbsXY", ";X (cm); Y (cm)", 1000, -0.5*_width, 1.5*_width, 1000, -0.5*_height, 1.5*_height);
  _absXY->SetDirectory(0);

  //TVector3 p0, p1, p2;

  double cornerR = 0.79;
  double cornerRTop = 0.813;
  double cornerRBot = 0.770;
  
  RectSurface* end1 = new RectSurface(0.0, 0.0, 0.0,
				      _width, 0.0, 0.0,
				      0.0, _height, 0.0,
				      _endcapReflectivity);
  _container.AddSurface( end1 );
  
  RectSurface* end2 = new RectSurface(0.0, 0.0, _length,
				      _width, 0.0, _length,
				      0.0, _height, _length,
				      _endcapReflectivity);
  _container.AddSurface( end2 );
  
  RectSurface* side1 = new RectSurface(0.0, _height, 0.0,
				       _width, _height, 0.0,
				       0.0, _height, _length,
				       _wallReflectivity);
  _container.AddSurface( side1 );

  RectSurface* side2 = new RectSurface(_width, _height, 0.0,
				       _width, 0.0, 0.0,
				       _width, _height, _length,
				       _wallReflectivity);
  _container.AddSurface( side2 );

  RectSurface* side3 = new RectSurface(0.0, 0.0, 0.0,
				       _width, 0.0, 0.0,
				       _width, 0.0, _length,
				       _wallReflectivity);
  _container.AddSurface( side3 );

  RectSurface* side4 = new RectSurface(0.0, 0.0, 0.0,
				       0.0, _height, 0.0,
				       0.0, _height, _length,
				       _wallReflectivity);
  _container.AddSurface( side4 );

  /*
  _container.AddCorner(0,      0,       cornerR, false, true);
  _container.AddCorner(_width, 0,       cornerR, false, false);
  _container.AddCorner(0,      _height, cornerR, true, true);
  _container.AddCorner(_width, _height, cornerR, true, false);
  */

  _container.AddCorner(0,      0,       cornerRBot, false, true);
  _container.AddCorner(_width, 0,       cornerRBot, false, false);
  _container.AddCorner(0,      _height, cornerRTop, true, true);
  _container.AddCorner(_width, _height, cornerRTop, true, false);
  
  //Adding fibers:
  /*
  double a = 1 + pow(_width/_height,2);
  double b = -(_height + _width*_width/_height);
  double c = pow(_height/2,2) + pow(_width/2,2) - pow(0.5*(_fiberSeparation + 0.07),2);
  double yf1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
  double xf1 = _width*yf1/_height;
  double yf2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
  double xf2 = _width*yf2/_height;
  */

  double m = _height/(_width - 2*cornerR);
  double xc = _width/2;
  //double yc = _height/2;
  double a = 1 + pow(m,2);
  double b = -2*xc*a;
  double c = a*pow(xc,2) - pow(0.5*(_fiberSeparation + 0.07),2);
  double xf1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
  double xf2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
  double yf1 = m*(xf1 - cornerR);
  double yf2 = m*(xf2 - cornerR);
  
  cout << "Fiber1: x = " << xf1 << " y = " << yf1 << endl;
  cout << "Fiber2: x = " << xf2 << " y = " << yf2 << endl;

  //0.07 is diameter of fiber
  //_container.AddFiber(0.07/2, xf1, yf1);
  //_container.AddFiber(0.07/2, xf2, yf2);
  _container.AddFiber(xf1, yf1, 0.07/2);
  _container.AddFiber(xf2, yf2, 0.07/2);
  
  //_container.SetFiber(_length, 0.07, 0.5*(6.15 + 0.07));

  //TGraph specGraph("data/bc517p_spectrum.dat");
  TGraph specGraph("data/novaSpectrum.dat");

  _spectrum = new TH1D("spectrum", "", 100, 360, 520);
  _spectrum->SetDirectory(0);
  for (int i = 1; i <= _spectrum->GetXaxis()->GetNbins(); ++i)
    {
      double bc = _spectrum->GetXaxis()->GetBinCenter(i);
      double value = specGraph.Eval(bc);
      _spectrum->SetBinContent(i, value);
    }

}

void PhotonTracer::GeneratePhoton( double z, double x, double y, bool doRandomXY )
{
  double rx,ry,rz,rt;
  //generate start position
  _gen.Sphere(rx,ry,rz,1.0);
  rt = 0.0;
  //generate start time
  //if (doTime) rt = _gen.Exp(_beta_emission);
  //else        rt = 0.0;

  double epsilon = 1e-3;
  if (doRandomXY)
    {
      x = _gen.Uniform(epsilon, _width - epsilon);
      y = _gen.Uniform(epsilon, _height  - epsilon);
    }

  //cout << "x = " << x << " y = " << y << endl;
  //if (x < 0) x = _width/2;
  //if (y < 0) y = _height/2;
  
  //TThread::Lock();
  double wavelength = _spectrum->GetRandom();
  //TThread::UnLock();


  double u = _gen.Rndm();
  bool isS = (u < 0.5); //s polarization if u < 0
    
  //TVector3 r0, rhat;
  //r0.SetXYZ( 0.5*_height, 0.5*_width, z );
  //r0.SetXYZ( x, y, z );
  //rhat.SetXYZ( rx, ry, rz );
  _photon.Init(x, y, z, rx, ry, rz, rt, wavelength, isS);      
}

bool PhotonTracer::TraceOnePhoton( double z, double x, double y, bool doRandomXY )
{
  double pvcIdx = 1.53;
  double oilIdx = 1.42;
  
  GeneratePhoton( z, x, y, doRandomXY );
  double eps = 1e-3;

  if (!_container.IsInside(_photon))
    {
      //cout << "Photon outside the volume" << endl;
      return false;
    }

  cout << "Start Tracing" << endl;
  //double step_size(2000.0);
  while ( true )
    {
      cout << "New iteration" << endl;
      //_photon.Print();
      //stop tracing if weight drops too low
      //if (_photon.GetWeight() < 1.0e-5) return false;
      //stop tracing if pathlength gets too long
      //if (_photon.GetPathLength() > 600) return false;
      
      double u(0.0);
      double stepToInt(-9999.0), reflectivity(-9999.0);
      double surfNx(0.0), surfNy(0.0), surfNz(0.0);
      
      _container.GetBestIntersection( _photon, _photon.GetWavelength(), stepToInt, reflectivity, surfNx, surfNy, surfNz, false );

      //check to see if collected
      double stepToFiber(-9999.0);
      double angleToFiber(-9999.0);
      double stepToAbs(-9999.0);
      double fiberNx(0.0), fiberNy(0.0), fiberNz(0.0);
      bool hitsFiber = _container.BestFiberIntersection( _photon, stepToFiber, angleToFiber, fiberNx, fiberNy, fiberNz );

      stepToAbs = _gen.Exp(_beta_atten);
      
      //if hits fiber and wall, check which one is closer
      //if hits only fiber, do fiber check
      //if hits only wall, only do wall checks


      if (hitsFiber)
	{
	  //absorption happens first -> kill the photon
	  if (stepToAbs <= stepToFiber)
	    {
	      _photon.Step(stepToAbs);
	      _attenXY->Fill( _photon.GetX0(), _photon.GetY0() );
	      return false;
	    }
	  //hits fiber
	  else 
	    {
	      _photon.FiberAngle(angleToFiber);
	      _photon.Step( stepToFiber );
	      _absXY->Fill( _photon.GetX0(), _photon.GetY0() );
	      return true;
	    }
	}
      else if (stepToInt > 0)
	{
	  //absorption happens first -> kill the photon
	  if (stepToAbs <= stepToFiber)
	    {
	      _photon.Step(stepToAbs);
	      _attenXY->Fill( _photon.GetX0(), _photon.GetY0() );
	      return false;
	    }
	  //hits wall
	  else 
	    {
	      //check for specular Fresnel reflection first
	      double angleToSurf = acos(fabs(surfNx*_photon.GetXHat() + surfNy*_photon.GetYHat() + surfNz*_photon.GetZHat()));
	      double probReflect = Reflectivity(oilIdx, pvcIdx, angleToSurf, _photon.GetIsS());
	      u = _gen.Uniform();
	      if (u < probReflect)
		{
		  _photon.Step( stepToInt );
		  _photon.Reflect( surfNx, surfNy, surfNz );
		  _refXY->Fill( _photon.GetX0(), _photon.GetY0() );
		  _photon.Step( eps );
		  continue;
		}
	      //if it didn't reflect specularly, check for diffuse reflection
	      u = _gen.Uniform();
	      probReflect = reflectivity;
	      if (u < probReflect ) 
		{
		  _photon.Step( stepToInt );
		  _photon.ReflectDiffuse( surfNx, surfNy, surfNz );
		  _refXY->Fill( _photon.GetX0(), _photon.GetY0() );
		  _photon.Step( eps );
		  continue;
		}		  
	      else
		{
		  //Failed at reflecting
		  _photon.Step( stepToInt );
		  _absXY->Fill( _photon.GetX0(), _photon.GetY0() );
		  return false;
		}	      
	    }
	}
      else
	{
	  //oops, no reflecting surfaces
	  //generate new photon 
	  cout << "No reflection found! - reset" << endl;
	  _container.IsInside(_photon, true);
	  //_photon.Print();
	  _container.GetBestIntersection( _photon, _photon.GetWavelength(), stepToInt, reflectivity, surfNx, surfNy, surfNz, true );
	  cout << "----------------------------" << endl;
	  return false;
	  GeneratePhoton( z, x, y, doRandomXY);
	  continue;
	}

      /*
	  
      if (hitsFiber && (stepToInt <= 0 || stepToFiber < stepToInt))
	{
	  //does photon attenuate before reaching the fiber?
	  u = _gen.Uniform();
	  double probAtten = 1.0 - exp( -stepToFiber/_beta_atten );
	  if ( u < probAtten ) {

	    return false;
	  }

	  _photon.FiberAngle(angleToFiber);
	  if (stepToFiber > eps) stepToFiber -= eps; 
	  else stepToFiber = 0.0;
	  _photon.Step( stepToFiber );
	  _absXY->Fill( _photon.GetX0(), _photon.GetY0() );

	  

	  
	  //if angle is greater than critical angle, reflect off the fiber
	  //if (angleToFiber > 1.336) {
	  //  _photon.Reflect( fiberNormal );
	  //  continue;
	  // }

	  //otherwise, assume the fiber accepts the photon
	  return true;
	}
      else if (stepToInt > 0)
	{
	  //NA = 0.72 = 1.46*sinTheta (1.46 = index of refraction of oil)
	  //The half-angle at intersection must be less than the angle implied by the numerical aperature to enter (0.5157)
	  //If instead we take critAngle = asin(1.42/1.46) [outer cladding to oil], max angle is 1.336
	  //if (hitsFiber && stepToFiber < stepToInt && angleToFiber <= 1.336)
	  //if (hitsFiber && stepToFiber < stepToInt && angleToFiber <= 0.5157)
	    //if (hitsFiber && stepToFiber < stepToInt)
	    //{
	  //_photon.FiberAngle(angleToFiber);
	  //  _photon.Step( stepToFiber );
	  //  return true;
	      //if (angleToFiber <= 0.5157) return true;
	      //if (angleToFiber <= 1.336) return true;
	      //else return false;
	  //}
	  
	  //bool isIntersection(false);
	  //double step(0.0);
	  //if ( stepToInt < step_size )
	  //{
	  //  isIntersection = true;
	  //  step = stepToInt;
	  //}
	  //else
	  //{
	  //  isIntersection = false;
	  //  step = step_size;
	  //}

	  if (stepToInt > eps) stepToInt -= eps;
	  else stepToInt = 0.0;
	  _photon.Step( stepToInt );
	  _absXY->Fill( _photon.GetX0(), _photon.GetY0() );


	  //check to see if collected
	  //if ( _container.InFiber( _photon ) ) return true;
	  //u = _gen.Uniform();
	  //double probCollect = (step/_beta_collect)*exp( -step/_beta_collect );
	  //if ( u < probCollect ) return true;

	  //check to see if attenuated
	  u = _gen.Uniform();
	  //integral[x_0^s] 1/beta * exp(-x/beta) dx = 1 - exp(-s/beta)
	  //double probAtten = (stepToInt/_beta_atten)*exp( -stepToInt/_beta_atten );
	  double probAtten = 1.0 - exp( -stepToInt/_beta_atten );
	  if ( u < probAtten ) {
	    _attenXY->Fill( _photon.GetX0(), _photon.GetY0() );
	    return false;
	  }
	  _photon.ScaleWeight( probAtten );
	  
	  //if (isIntersection)
	  //{
	  u = _gen.Uniform();
	  //normal*incidentDir is Lambertian weighting of reflectivity for diffuse surfaces
	  //currently neglecting Fresnel enhancement to reflectivity when light is glancing
	  //double probReflect = reflectivity*(fabs(surfNormal*_photon.GetRHat()));

	  //check for specular Fresnel reflection first
	  double angleToSurf = acos(fabs(surfNx*_photon.GetXHat() + surfNy*_photon.GetYHat() + surfNz*_photon.GetZHat()));
	  double probReflect = Reflectivity(oilIdx, pvcIdx, angleToSurf, _photon.GetIsS());
	  if ( u < probReflect )
	    {
	      _photon.Reflect( surfNx, surfNy, surfNz );
	      _refXY->Fill( _photon.GetX0(), _photon.GetY0() );
	      continue;
	    }
	  //if it didn't reflect specularly, check for diffuse reflection
	  u = _gen.Uniform();
	  probReflect = reflectivity;
	  _photon.ScaleWeight( probReflect );
	  if (u < probReflect ) 
	    {
	      _photon.ReflectDiffuse( surfNx, surfNy, surfNz );
	      _refXY->Fill( _photon.GetX0(), _photon.GetY0() );
	      continue;
	    }
	  //failed at reflecting
	  else
	    {
	      return false;
	    }
	}
      else
	{
	  //oops, no reflecting surfaces
	  //generate new photon 
	  cout << "No reflection found! - reset" << endl;
	  _container.IsInside(_photon, true);
	  //_photon.Print();
	  _container.GetBestIntersection( _photon, _photon.GetWavelength(), stepToInt, reflectivity, surfNx, surfNy, surfNz, true );
	  cout << "----------------------------" << endl;
	  return false;
	  GeneratePhoton( z, x, y, doRandomXY);
	  continue;
	}
      */
    }
  
  //if it made it here, exceeded the path length limit
  return false;
}

