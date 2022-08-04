#ifndef RECTSURFACE_HH
#define RECTSURFACE_HH

#include <TGraph.h>
#include "include/Trajectory.hh"
#include <string>

//define a rectangular surface
class RectSurface
{
public: 
  RectSurface(double p0x, double p0y, double p0z,
	      double p1x, double p1y, double p1z,
	      double p2x, double p2y, double p2z,
	      std::string reflectivityFile);
  double GetNormalX() { return _nhatx; }
  double GetNormalY() { return _nhaty; }
  double GetNormalZ() { return _nhatz; }

  bool IsInBounds( double px, double py, double pz );
  double GetReflectivity(double wavelength);
  double StepToIntersect( Trajectory traj );

private:
  double  _p0x; //initial point
  double  _p0y; //initial point
  double  _p0z; //initial point
  double  _v1x; //p1 - p0
  double  _v1y; //p1 - p0
  double  _v1z; //p1 - p0
  double  _v2x; //p2 - p0
  double  _v2y; //p2 - p0
  double  _v2z; //p2 - p0
  TGraph* _reflectivity;
  double  _nhatx;
  double  _nhaty;
  double  _nhatz;
  double  _v1mag;
  double  _v2mag;
};

#endif
