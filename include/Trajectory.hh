#ifndef TRAJECTORY_HH
#define TRAJECTORY_HH

#include <TRandom3.h>

class Trajectory
{
public:
  Trajectory();
  Trajectory(double x0, double y0, double z0, double xhat, double yhat, double zhat, double time, double wavelength, bool isS);
  void Init(double x0, double y0, double z0, double xhat, double yhat, double zhat, double time, double wavelength, bool isS)
  {
    _x0 = x0; _y0 = y0; _z0 = z0;
    _xhat = xhat; _yhat = yhat; _zhat = zhat;
    _pathlength = 0.0; _time = time;
    _wavelength = wavelength;
    _weight = 1.0; _nReflections = 0;
    _isS = isS;
  }
  void Print();
  void Step(double step);
  void Step(double step, double idxRefraction);
  void Step(double step, double idxRefraction, double extraTime);
  void Reflect(double nx, double ny, double nz);
  void ReflectDiffuse(double nx, double ny, double nz);
  void Refract(double nx, double ny, double nz, double idxRef1, double idxRef2);
  void ScaleWeight(double factor) { _weight *= factor; }
  void FiberAngle(double fiberAngle) {_fiberAngle = fiberAngle;}
  void Wavelength(double wavelength) {_wavelength = wavelength;}
  void R0(double x, double y, double z) {_x0 = x, _y0 = y; _z0 = z;}
  void RHat(double rx, double ry, double rz) {_xhat = rx; _yhat = ry; _zhat = rz;}
  void IsS(bool isS) {_isS = isS;}
  void Time(double time) {_time = time;}
  
  double   GetX0() {return _x0;}
  double   GetY0() {return _y0;}
  double   GetZ0() {return _z0;}
  double   GetXHat() {return _xhat;}
  double   GetYHat() {return _yhat;}
  double   GetZHat() {return _zhat;}
  double   GetPathLength() {return _pathlength;}
  double   GetTime() {return _time;}
  double   GetWeight() {return _weight;}
  double   GetWavelength() {return _wavelength;}
  int      GetNReflections() {return _nReflections;}
  double   GetFiberAngle() {return _fiberAngle;}
  bool     GetIsS() {return _isS;}
  
private:
  double   _x0;
  double   _y0;
  double   _z0;
  double    _xhat;
  double    _yhat;
  double    _zhat;
  double   _pathlength;
  double   _time;
  double   _weight;
  double   _wavelength;
  int      _nReflections;
  double   _fiberAngle;
  TRandom3 _gen;
  bool     _isS;
};

#endif
