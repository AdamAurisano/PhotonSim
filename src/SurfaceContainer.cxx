#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <array>

#include "include/Trajectory.hh"
#include "include/RectSurface.hh"
#include "include/SurfaceContainer.hh"

using namespace std;
using namespace TMath;

SurfaceContainer::SurfaceContainer(string reflectivityFile)
{
  _reflectivity = new TGraph(reflectivityFile.c_str());
}

SurfaceContainer::~SurfaceContainer()
{
  vector< RectSurface* >::iterator it;
  for (it = _surfaces.begin(); it != _surfaces.end(); ++it)
    {
      RectSurface* surf = *it;
      delete surf;
    }
  delete _reflectivity;
}

void SurfaceContainer::AddFiber( double x0, double y0, double r  )
{
  _fiberX0.push_back(x0);
  _fiberY0.push_back(y0);
  _fiberR.push_back(r);
}

bool SurfaceContainer::IsInside(Trajectory& traj, bool debug )
{
  vector<array<double,3> > unitVec;
  unitVec.push_back( {1, 0, 0} );
  unitVec.push_back( {-1, 0, 0} );
  unitVec.push_back( {0, 1, 0} );
  unitVec.push_back( {0, -1, 0} );
  unitVec.push_back( {0, 0, 1} );
  unitVec.push_back( {0, 0, -1} );

  bool goodSurface(false), goodCorner(false);
  for (int i = 0; i < (int)unitVec.size(); ++i)
    {
      Trajectory testPho = Trajectory(traj.GetX0(), traj.GetY0(), traj.GetZ0(),
				      unitVec[i][0], unitVec[i][1], unitVec[i][2],
				      0.0, 420.,true);
      double tmp_step = -9999;
      //double tmp_refl = -9999;
      //double tmp_nx(0.0), tmp_ny(0.0), tmp_nz(0.0);

      for (auto it = _surfaces.begin(); it != _surfaces.end(); ++it)
	{
	  RectSurface* tmp_surf = *it;
	  double tmp_step = tmp_surf->StepToIntersect( traj );
	  if (tmp_step > 0) goodSurface = true;
	}

      for (int iC = 0; iC < (int)_x0.size(); ++iC)
	{
	  double tmp_step(0);
	  double tmp_nx(0), tmp_ny(0), tmp_nz(0);
	  FindCornerIntersection(traj, _x0[iC], _xc[iC], _y0[iC], _yc[iC], _rc[iC], tmp_step, tmp_nx, tmp_ny, tmp_nz);
	  if (tmp_step > 0) goodCorner = true;
	}
      
      //GetBestIntersection(testPho, 420., tmp_step, tmp_refl, tmp_nx, tmp_ny, tmp_nz);
      if (debug) cout << "RHat: ( " << unitVec[i][0] << ", " << unitVec[i][1] << ", " << unitVec[i][2] << "): " << tmp_step << endl;
      //if (tmp_step < 0) return false;
    }
  //return true;
  return goodSurface && goodCorner;
}

void SurfaceContainer::AddCorner(double x0, double y0, double r, bool isTop, bool isLeft)
{
  _x0.push_back(x0);
  _y0.push_back(y0);
  _rc.push_back(r);
  double xc(0), yc(0);
  if (isTop) yc = y0 - r;
  else yc = y0 + r;
  if (isLeft) xc = x0 + r;
  else xc = x0 - r;
  _xc.push_back(xc);
  _yc.push_back(yc);
}

void SurfaceContainer::FindCornerIntersection(Trajectory& traj, double x0, double xc, double y0, double yc, double r,
					      double& step, double& nx, double& ny, double& nz)
{
  double step1(-999), step2(-999), angle1(-999), angle2(-999);
  double nx1(0.0), ny1(0.0), nz1(0.0);
  double nx2(0.0), ny2(0.0), nz2(0.0);
  
  FindCylinderIntersections(traj, xc, yc, r, step1, angle1, nx1, ny, nz1, step2, angle2, nx2, ny2, nz2);

  if (step1 <= 0 && step2 <= 0)
    {
      step = step1;
      nx = nx1;
      ny = ny1;
      nz = nz1;
      return;
    }

  bool isValid1(false), isValid2(false);
  double xint1 = traj.GetX0() + step1*traj.GetXHat();
  double xint2 = traj.GetX0() + step2*traj.GetXHat();

  double yint1 = traj.GetY0() + step1*traj.GetYHat();
  double yint2 = traj.GetY0() + step2*traj.GetYHat();

  //if xint1 is between x0 and xc  and yint1 is between y0 and yc
  if ( ((xint1 >= x0 && xint1 <= xc) || (xint1 <= x0 && xint1 >= xc))
       && ((yint1 >= y0 && yint1 <= yc) || (yint1 <= y0 && yint1 >= yc)) ) isValid1 = true;
  if (step1 <= 0) isValid1 = false;
  
  //if xint2 is between x0 and xc  and yint2 is between y0 and yc
  if ( ((xint2 >= x0 && xint2 <= xc) || (xint2 <= x0 && xint2 >= xc))
       && ((yint2 >= y0 && yint2 <= yc) || (yint2 <= y0 && yint2 >= yc)) ) isValid2 = true;
  if (step2 <= 0) isValid2 = false;
  
  if (isValid1 && !isValid2)
    {
      step = step1;
      nx = nx1;
      ny = ny1;
      nz = nz1;
    }
  else if (!isValid1 && isValid2)
    {
      step = step2;
      nx = nx2;
      ny = ny2;
      nz = nz2;
    }
  else if (isValid1 && isValid2)
    {
      if (step1 < step2)
	{
	  step = step1;
	  nx = nx1;
	  ny = ny1;
	  nz = nz1;
	}
      else
	{
	  step = step2;
	  nx = nx2;
	  ny = ny2;
	  nz = nz2;
	}
    }
  else
    {
      step = -999;
      nx = nx1;
      ny = ny1;
      nz = nz1;
    }
  return;
}

void SurfaceContainer::GetBestIntersection( Trajectory& traj,  double wavelength, double& step, double& reflectivity, double& nx, double& ny, double& nz, bool debug )
{

  RectSurface* surf(0);
  double best_step(9999.0);
  double best_reflectivity(9999.0);
  double best_nx(0.0), best_ny(0.0), best_nz(0.0);
  //TVector3 best_normal;
  
  vector< RectSurface* >::iterator it;
  if (debug)
    {
      cout << "-----Debug GetBestIntersection-----------" << endl;
    }

  //cout << "--------------------------" << endl;
  //cout << "Wall Intersections" << endl;
  for (it = _surfaces.begin(); it != _surfaces.end(); ++it)
    {
      RectSurface* tmp_surf = *it;
      double tmp_step = tmp_surf->StepToIntersect( traj );
      //  cout << tmp_step << endl;
      //TVector3 tmp_normal = tmp_surf->GetNormal();
      double tmp_nx = tmp_surf->GetNormalX();
      double tmp_ny = tmp_surf->GetNormalY();
      double tmp_nz = tmp_surf->GetNormalZ();
      if (debug)
	{
	  cout << "Step to intersect = " << tmp_step << endl;
	}
      //if ( tmp_step > 1e-8 )
      //{
      //  if ( (best_step < 1e-8) || (tmp_step < best_step) )
      if (tmp_step >= 0 && tmp_step < best_step)
	{
	  best_step = tmp_step;
	  best_nx = tmp_nx;
	  best_ny = tmp_ny;
	  best_nz = tmp_nz;
	  surf = tmp_surf;	      
	  best_reflectivity = surf->GetReflectivity(wavelength);
	}
      //}
    }

  //cout << "Corner Intersections" << endl;
  for (int iC = 0; iC < (int)_x0.size(); ++iC)
    {
      double tmp_step(0);
      double tmp_nx(0), tmp_ny(0), tmp_nz(0);
      FindCornerIntersection(traj, _x0[iC], _xc[iC], _y0[iC], _yc[iC], _rc[iC], tmp_step, tmp_nx, tmp_ny, tmp_nz);
      //cout << tmp_step << endl;
      //if ( tmp_step > 1e-8 )
      //{
      //  if ( (best_step < 1e-8) || (tmp_step < best_step) )
      if (tmp_step >= 0 && tmp_step < best_step)
	{
	  best_step = tmp_step;
	  //normal for corner is opposite what you would have for a cylinder
	  best_nx = -tmp_nx;
	  best_ny = -tmp_ny;
	  best_nz = -tmp_nz;
	  best_reflectivity = _reflectivity->Eval(wavelength);
	}
      //}
    }
  //cout << "best_step = " << best_step << endl;
  
  step = best_step;
  nx = best_nx;
  ny = best_ny;
  nz = best_nz;
  reflectivity = best_reflectivity;
}

//fiberX and fiberY are the center coordinate - this still works for the corner
void SurfaceContainer::FindNormal(double intersectionX, double intersectionY, double intersectionZ,
				  double fiberX, double fiberY,
				  double& nx, double& ny, double& nz)
{
  double cx = fiberX;
  double cy = fiberY;
  double cz = intersectionZ;
  //TVector3 center(fiberX, fiberY, intersection.Z());
  nx = cx - intersectionX;
  ny = cy - intersectionY;
  nz = cz - intersectionZ;
  double nMag = sqrt( nx*nx + ny*ny + nz*nz);
  //TVector3 normal = center - intersection;
  nx /= nMag;
  ny /= nMag;
  nz /= nMag;
  //normal *= 1.0/normal.Mag();
  //return normal;
}

void SurfaceContainer::FindCylinderIntersections(Trajectory& traj, double xc, double yc, double r,
						 double& step1, double& angle1, double& nx1, double& ny1, double& nz1,
						 double& step2, double& angle2, double& nx2, double& ny2, double& nz2)
{
  double x = traj.GetX0();
  double y = traj.GetY0();
  double z = traj.GetZ0();
  
  double dx = traj.GetXHat();  
  double dy = traj.GetYHat();  
  double dz = traj.GetZHat();

  double a = dx*dx + dy*dy;
  double b = 2.0*(x - xc)*dx + 2.0*(y-yc)*dy;
  double c = pow(x - xc,2) + pow(y - yc,2) - pow(r,2);
  double det = b*b - 4.0*a*c;  

  if (det < 0) 
    {
      step1 = -9999.0;
      step2 = -9999.0;
    }
  else if (det == 0)
    {
      step1 = -b/(2.0*a);
      step2 = -9999.0;
    }
  else 
    {
      step1 = (-b - sqrt(det))/(2.0*a);
      step2 = (-b + sqrt(det))/(2.0*a);
    }

  if (step1 <= 0) angle1 = -9999.0;
  else
    {
      double intersectionX = x + step1*dx;
      double intersectionY = y + step1*dy;
      double intersectionZ = z + step1*dz;
      //TVector3 intersection( x + step1*dx, y + step1*dy, z + step1*dz );
      FindNormal( intersectionX, intersectionY, intersectionZ, xc, yc, nx1, ny1, nz1);
      angle1 = acos(fabs( traj.GetXHat()*nx1 + traj.GetYHat()*ny1 + traj.GetZHat()*nz1 ) );      
    }

  if (step2 <= 0) angle2 = -9999.0;
  else
    {
      double intersectionX = x + step2*dx;
      double intersectionY = y + step2*dy;
      double intersectionZ = z + step2*dz;
      //TVector3 intersection( x + step2*dx, y + step2*dy, z + step2*dz );
      FindNormal( intersectionX, intersectionY, intersectionZ, xc, yc, nx2, ny2, nz2 );
      angle2 = acos(fabs(traj.GetXHat()*nx2 + traj.GetYHat()*ny2 + traj.GetZHat()*nz2));
      //angle2 = acos(fabs(traj.GetRHat()*normal2));
    }
  return;
}

void SurfaceContainer::FindFiberIntersection(Trajectory& traj, double fiberX, double fiberY, double fiberR, double& step, double& angle, double& nx, double& ny, double& nz)
{
  double step1(-9999.0);
  double step2(-9999.0);
  
  double angle1(-9999.0);
  double angle2(-9999.0);

  double nx1(0.0), ny1(0.0), nz1(0.0);
  double nx2(0.0), ny2(0.0), nz2(0.0);

  FindCylinderIntersections(traj, fiberX, fiberY, fiberR,
			    step1, angle1, nx1, ny1, nz1,
			    step2, angle2, nx2, ny2, nz2);
  
  if (step1 <= 0 && step2 <= 0)
    {
      step = -9999.0;
      angle = -9999.0;
      nx = nx1;
      ny = ny1;
      nz = nz1;
    }
  else if (step1 <= 0 && step2 > 0)
    {
      step = step2;
      angle = angle2;
      nx = nx2;
      ny = ny2;
      nz = nz2;
    }
  else if (step1 > 0 && step2 <= 0)
    {
      step = step1;
      angle = angle1;
      nx = nx1;
      ny = ny1;
      nz = nz1;
    }
  else
    {
      if (step1 < step2)
	{
	  step = step1;
	  angle = angle1;
	  nx = nx1;
	  ny = ny1;
	  nz = nz1;
	}
      else
	{
	  step = step2;
	  angle = angle2;
	  nx = nx2;
	  ny = ny2;
	  nz = nz2;
	}
    }
  return;
}


bool SurfaceContainer::BestFiberIntersection(Trajectory& traj, double& step, double& angle, double& nx, double& ny, double& nz)
{
  double best_step(-9999.0);
  double best_angle(-9999.0);

  double best_nx(0.0), best_ny(0.0), best_nz(0.0);
  double tmp_nx(0.0), tmp_ny(0.0), tmp_nz(0.0);
  
  bool foundIntersection(false);
  for (size_t i = 0; i < _fiberX0.size(); ++i)
    {
      double tmp_step(-9999.0), tmp_angle(-9999.0);
      //TVector3 tmp_normal;
      FindFiberIntersection(traj, _fiberX0[i], _fiberY0[i], _fiberR[i], tmp_step, tmp_angle, tmp_nx, tmp_ny, tmp_nz);
      if (tmp_step > 0.0)
	{
	  if ( (best_step < 0) || (tmp_step < best_step) )
	    {
	      foundIntersection = true;
	      best_step = tmp_step;
	      best_angle = tmp_angle;
	      best_nx = tmp_nx;
	      best_ny = tmp_ny;
	      best_nz = tmp_nz;
	    }
	}
    }

  step = best_step;
  angle = best_angle;
  nx = best_nx;
  ny = best_ny;
  nz = best_nz;

  return foundIntersection;
}

