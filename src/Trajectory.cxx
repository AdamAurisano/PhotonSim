#include <iostream>
#include <iomanip>

#include <TVector3.h>

#include "include/Trajectory.hh"

using namespace std;
using namespace TMath;

Trajectory::Trajectory() : _x0(0.0), _y0(0.0), _z0(0.0),
			   _xhat(0.0), _yhat(0.0), _zhat(0.0),
			   _pathlength(0.0), _time(0.0), _weight(1.0),
			   _wavelength(420.0), _nReflections(0),
			   _fiberAngle(-999.0), _gen(0), _isS(true)
{}

Trajectory::Trajectory(double x0, double y0, double z0,
		       double xhat, double yhat, double zhat,
		       double time, double wavelength, bool isS)
{
  Init(x0, y0, z0, xhat, yhat, zhat,  time, wavelength, isS);
}

void Trajectory::Print()
{
  cout << "x0 = " << _x0 << " y0 = " << _y0 << " z0 = " << _z0 << endl;
  cout << "xhat = " << _xhat << " yhat = " << _yhat << " zhat = " << _zhat << endl;
  cout << "Pathlength = " << _pathlength << " N Reflections = " << _nReflections << endl;
}

void Trajectory::Step(double step)
{
  this->Step(step, 1.46);
}

void Trajectory::Step(double step, double idxRefraction)
{
  this->Step(step, idxRefraction, 0.0);
}

void Trajectory::Step(double step, double idxRefraction, double extraTime)
{
  //updates _r0 after taking a step
  _x0 += step*_xhat;
  _y0 += step*_yhat;
  _z0 += step*_zhat;
  _pathlength += step;

  const double c = 29.9792458;
  const double cn = c/idxRefraction;

  _time += step/cn;
  _time += extraTime;

  return;
}

void Trajectory::Reflect(double nx, double ny, double nz)
{
  //double dp = _rhat.Dot(nhat);
  double dp = _xhat*nx + _yhat*ny + _zhat*nz;
  //int sign = 1;
  //if (dp > 0) sign = -1;

  //double c1 = -sign*nhat*_rhat;
  //_rhat += 2.0*c1*sign*nhat;

  _xhat += -2.0*dp*nx;
  _yhat += -2.0*dp*ny;
  _zhat += -2.0*dp*nz;
  
  ++_nReflections;
  return;
}

void Trajectory::ReflectDiffuse(double nx, double ny, double nz)
{
  //double Th = _gen.Uniform(0.0, 0.5*TMath::Pi() - 1.0e-3);
  //double cosTh = cos(Th);
  //double sinTh = sin(Th);

  double epsilon = 1.0e-3;
  double cosTh = _gen.Uniform(0.0, 1.0 - epsilon);
  double sinTh = sqrt( 1 - pow(cosTh,2) );

  //double sinTh = sqrt( _gen.Rndm() );
  //double cosTh = sqrt( 1 - pow(sinTh,2) );
  double phi = 2*TMath::Pi()*_gen.Rndm();

  //find first orthogonal to nhat
  double tHat1_x(0.0), tHat1_y(0.0), tHat1_z(0.0);
  double xx = nx < 0.0 ? -nx : nx;
  double yy = ny < 0.0 ? -ny : ny;
  double zz = nz < 0.0 ? -nz : nz;
  if (xx < yy) {
    if (xx < zz) {tHat1_x = 0; tHat1_y = nz; tHat1_z = -ny;}
    else         {tHat1_x = ny; tHat1_y = -nx; tHat1_z = 0;}
  }
  else {
    if (yy < zz) {tHat1_x = -nz; tHat1_y = 0; tHat1_z = nx;}
    else         {tHat1_x = ny;  tHat1_y = -nx; tHat1_z = 0;}
  }
  
  double tHat1Mag = sqrt( tHat1_x*tHat1_x + tHat1_y*tHat1_y + tHat1_z*tHat1_z );
  tHat1_x /= tHat1Mag;
  tHat1_y /= tHat1Mag;
  tHat1_z /= tHat1Mag;

  //TVector3 tHat1 = nhat.Orthogonal();
  //tHat1 *= 1.0/tHat1.Mag();
  
  //TVector3 tHat2 = nhat.Cross(tHat1);
  //tHat2 *= 1.0/tHat2.Mag();
  double tHat2_x = ny*tHat1_z - nz*tHat1_y;
  double tHat2_y = nz*tHat1_x - nx*tHat1_z;
  double tHat2_z = nx*tHat1_y - ny*tHat1_x;

  double tHat2Mag = sqrt( tHat2_x*tHat2_x + tHat2_y*tHat2_y + tHat2_z*tHat2_z );
  tHat2_x /= tHat2Mag;
  tHat2_y /= tHat2Mag;
  tHat2_z /= tHat2Mag;
  
  //if reflecting, rhat*nhat should be negative
  //double dp = _rhat.Dot(nhat);
  double dp = _xhat*nx + _yhat*ny + _zhat*nz;
  int sign = 1;
  if (dp > 0) sign = -1;
  
  //_rhat = sign*cosTh*nhat + sinTh*cos(phi)*tHat1 + sinTh*sin(phi)*tHat2;
  double cosp = cos(phi);
  double sinp = sin(phi);
  _xhat = sign*cosTh*nx + sinTh*cosp*tHat1_x + sinTh*sinp*tHat2_x;
  _yhat = sign*cosTh*ny + sinTh*cosp*tHat1_y + sinTh*sinp*tHat2_y;
  _zhat = sign*cosTh*nz + sinTh*cosp*tHat1_z + sinTh*sinp*tHat2_z;
  
  ++_nReflections;
  return;
}

void Trajectory::Refract(double nx, double ny, double nz, double idxRef1, double idxRef2)
{
  //double dp = _rhat.Dot(nhat);
  double dp = _xhat*nx + _yhat*ny + _zhat*nz;
  int sign = 1;
  if (dp > 0) sign = -1;

  const double n = idxRef1/idxRef2;
  //const double c1 = fabs(nhat*_rhat);
  const double c1 = fabs(dp);
  const double c2 = sqrt(1.0 - n*n*(1.0 - c1*c1));

  //_rhat = (n*_rhat) + (n*c1 - c2)*sign*nhat;
  //_rhat *= 1.0/_rhat.Mag();

  _xhat = n*_xhat + (n*c1 - c2)*sign*nx;
  _yhat = n*_yhat + (n*c1 - c2)*sign*ny;
  _zhat = n*_zhat + (n*c1 - c2)*sign*nz;
  double rHatMag = sqrt( _xhat*_xhat + _yhat*_yhat + _zhat*_zhat );
  _xhat /= rHatMag;
  _yhat /= rHatMag;
  _zhat /= rHatMag;
  
  return;
}



