#ifndef SURFACECONTAINER_HH
#define SURFACECONTAINER_HH

#include <vector>
#include <string>
#include "TGraph.h"
#include "include/RectSurface.hh"
#include "include/Trajectory.hh"

class SurfaceContainer
{
public:
  SurfaceContainer(std::string reflectivityFile);
  ~SurfaceContainer();
  void AddSurface( RectSurface* rect ) { _surfaces.push_back( rect ); }
  void AddFiber( double x0, double y0, double r  );
  
  //x0 and y0 are the bounding surfaces for the corner that are tangent to circle
  void AddCorner( double x0, double y0, double r, bool isTop, bool isLeft);
  
  bool BestFiberIntersection( Trajectory& traj, double& step, double& angle, double& nx, double& ny, double& nz);

  //best wall or corner
  void GetBestIntersection( Trajectory& traj, double wavelength, double& step, double& reflectivity, double& nx, double& ny, double& nz, bool debug=false );

  //bool IsInside(Trajectory& traj, bool debug = false);
  bool IsInside(double x, double y, double W, double H, double RTop, double RBottom);
  
private:
  //normal of a cylinder
  void FindNormal(double interX, double interY, double interZ, double fiberX, double fiberY, double& nx, double& ny, double& nz);
  //finds intersections for cylinders - used for both fibers and corners
  void     FindCylinderIntersections(Trajectory& traj, double xc, double yc, double r,
				     double& step1, double& angle1, double& nx1, double& ny1, double& nz1,
				     double& step2, double& angle2, double& nx2, double& ny2, double& nz2);
  
  void     FindFiberIntersection(Trajectory& traj, double fiberX, double fiberY, double fiberR, double& step, double& angle, double& nx, double& ny, double& nz);
  
  void     FindCornerIntersection(Trajectory& traj, double x0, double xc, double y0, double yc, double r, double& step, double& n1, double& ny, double& nz);
  
  std::vector< RectSurface* > _surfaces;
  std::vector<double> _fiberX0;
  std::vector<double> _fiberY0;
  std::vector<double> _fiberR;
  
  //corners
  std::vector<double> _x0;
  std::vector<double> _y0;
  std::vector<double> _xc;
  std::vector<double> _yc;
  std::vector<double> _rc;

  TGraph* _reflectivity;
  
};

#endif
