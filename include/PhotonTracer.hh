#ifndef PHOTONTRACER_HH
#define PHOTONTRACER_HH

#include <TRandom3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <string>
#include "include/SurfaceContainer.hh"
#include "include/Trajectory.hh"

class PhotonTracer
{
public:
  PhotonTracer(double fiberSeparation);
  void GeneratePhoton( double z, double x = -999, double y = -999, bool doRandomXY=false );
  bool TraceOnePhoton( double z, double x = -999, double y = -999, bool doRandomXY=false );
  Trajectory& GetPhoton() { return _photon; }
  double Reflectivity(double idx1, double idx2, double theta, bool isS);
  
  //double GetWallReflectivity()   { return _wallReflectivity; }
  //double GetEndcapReflectivity() { return _endcapReflectivity; }
  double GetWidth()  { return _width; }
  double GetHeight() { return _height; }
  double GetLength() { return _length; }
  TH2D*   GetRefXY()  { return _refXY; }
  TH2D*   GetAttenXY()  { return _attenXY; }
  TH2D*   GetAbsXY()  { return _absXY; }
  
private:
  TRandom3          _gen;
  SurfaceContainer  _container;
  Trajectory        _photon;
  TH1D*             _spectrum;
  const double      _beta_atten;
  const double      _beta_emission;
  const std::string _wallReflectivity;
  const std::string _endcapReflectivity;
  const double      _width;
  const double      _height;
  const double      _length;
  const double      _rtop;
  const double      _rbottom;
  const double      _fiberSeparation;
  TH2D*             _refXY;
  TH2D*             _attenXY;
  TH2D*             _absXY;
};

#endif

