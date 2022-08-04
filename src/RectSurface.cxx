#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <TGraph.h>

#include "include/RectSurface.hh"

using namespace std;
using namespace TMath;

RectSurface::RectSurface(double p0x, double p0y, double p0z,
			 double p1x, double p1y, double p1z,
			 double p2x, double p2y, double p2z,
			 string reflectivityFile) : _p0x(p0x), _p0y(p0y), _p0z(p0z)
{
  _v1x = p1x - p0x;
  _v1y = p1y - p0y;
  _v1z = p1z - p0z;
  _v2x = p2x - p0x;
  _v2y = p2y - p0y;
  _v2z = p2z - p0z;
    
  _v1mag = sqrt( _v1x*_v1x + _v1y*_v1y + _v1z*_v1z );
  _v2mag = sqrt( _v2x*_v2x + _v2y*_v2y + _v2z*_v2z );

  _nhatx = _v1y*_v2z - _v1z*_v2y;
  _nhaty = _v1z*_v2x - _v1x*_v2z;
  _nhatz = _v1x*_v2y - _v1y*_v2x;

  double nhatMag = sqrt(_nhatx*_nhatx + _nhaty*_nhaty + _nhatz*_nhatz);
  //_nhat = _v1.Cross(_v2);
  //_nhat *= 1.0/_nhat.Mag();
  _nhatx /= nhatMag;
  _nhaty /= nhatMag;
  _nhatz /= nhatMag;
  
  _reflectivity = new TGraph(reflectivityFile.c_str());
}

double RectSurface::GetReflectivity(double wavelength)
{
  return _reflectivity->Eval(wavelength); 
}


bool RectSurface::IsInBounds( double px, double py, double pz )
{
  double pvecx = px - _p0x;
  double pvecy = py - _p0y;
  double pvecz = pz - _p0z;
  double u1 = (pvecx*_v1x + pvecy*_v1y + pvecz*_v1z)/(_v1mag*_v1mag);
  //TVector3 pvec = point - _point0;
  //double u1 = (pvec*_v1)/(_v1mag*_v1mag);
  bool inU1 = ( u1 >= 0 && u1 <= 1 );
  //double u2 = (pvec*_v2)/(_v2mag*_v2mag);
  double u2 = (pvecx*_v2x + pvecy*_v2y + pvecz*_v2z)/(_v2mag*_v2mag);
  bool inU2 = ( u2 >= 0 && u2 <= 1 );
  return inU1 && inU2;
}

double RectSurface::StepToIntersect( Trajectory traj )
{
  double step = -9999;
  //double den = traj.GetRHat()*_nhat;
  double den = traj.GetXHat()*_nhatx + traj.GetYHat()*_nhaty + traj.GetZHat();
  if( fabs(den) > 0.0 ) 
    {
      double diffx = _p0x - traj.GetX0();
      double diffy = _p0y - traj.GetY0();
      double diffz = _p0z - traj.GetZ0();
      //double tmp_step = ( (_point0 - traj.GetR0())*_nhat )/den;
      double tmp_step = (diffx*_nhatx + diffy*_nhaty + diffz*_nhatz)/den;
      //TVector3 newPoint = traj.GetR0() + tmp_step*traj.GetRHat();
      //if ( IsInBounds( newPoint ) ) step = tmp_step;
      double npx = traj.GetX0() + tmp_step*traj.GetXHat();
      double npy = traj.GetY0() + tmp_step*traj.GetYHat();
      double npz = traj.GetZ0() + tmp_step*traj.GetZHat();
      if ( IsInBounds( npx, npy, npz ) ) step = tmp_step;
    }
  return step;
}

