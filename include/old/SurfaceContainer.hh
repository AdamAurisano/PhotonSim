#ifndef SURFACECONTAINER_HH
#define SURFACECONTAINER_HH

#include <vector>
#include <string>
#include "TVector3.h"
#include "TGraph.h"
#include "include/RectSurface.hh"
#include "include/Trajectory.hh"

class SurfaceContainer
{
public:
  SurfaceContainer(std::string reflectivityFile);
  ~SurfaceContainer();
  void AddSurface( RectSurface* rect ) { _surfaces.push_back( rect ); }
  void AddFiber( double length, double fiberRadius, double fiberX, double fiberY)
  { _length.push_back(length); _fiberRadius.push_back(fiberRadius); _fiberX.push_back(fiberX); _fiberY.push_back(fiberY);}

  //x0 and y0 are the bounding surfaces for the corner that are tangent to circle
  void AddCorner( double x0, double y0, double r, bool isTop, bool isLeft);
  
  bool BestFiberIntersection( Trajectory traj, double& step, double& angle, TVector3& normal );

  //best wall or corner
  void GetBestIntersection( Trajectory traj, double wavelength, double& step, double& reflectivity, TVector3& normal, bool debug=false );

  bool IsInside(Trajectory traj, bool debug = false);
  
private:
  //normal of a cylinder
  TVector3 FindNormal(TVector3 intersection, double fiberX, double fiberY);
  //finds intersections for cylinders - used for both fibers and corners
  void     FindCylinderIntersections(Trajectory traj, double xc, double yc, double r,
				     double& step1, double& angle1, TVector3& normal1,
				     double& step2, double& angle2, TVector3& normal2);
  
  void     FindFiberIntersection(Trajectory traj, double fiberX, double fiberY, double fiberR, double& step, double& angle, TVector3& normal);
  
  void     FindCornerIntersection(Trajectory traj, double x0, double xc, double y0, double yc, double r, double& step, TVector3& normal);
  
  std::vector< RectSurface* > _surfaces;
  std::vector<double> _length;
  std::vector<double> _fiberRadius;
  std::vector<double> _fiberX;
  std::vector<double> _fiberY;

  //corners
  std::vector<double> _x0;
  std::vector<double> _y0;
  std::vector<double> _xc;
  std::vector<double> _yc;
  std::vector<double> _rc;

  TGraph* _reflectivity;
  
};

#endif
