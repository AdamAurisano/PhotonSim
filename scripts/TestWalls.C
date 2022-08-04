#include "include/SurfaceContainer.hh"
#include "include/Trajectory.hh"
#include "include/PhotonTracer.hh"

#include <TRandom3.h>
#include <TH2D.h>

#include <iostream>
#include <iomanip>

using namespace std;
using namespace TMath;


void TestWalls()
{
  double width(3.57), height(5.59), length(1550.0);
  double cornerR = 0.79;
  double cornerRTop = 0.813;
  double cornerRBot = 0.770;
  string wallReflectivity = "data/wallReflectivity.dat";
  string endcapReflectivity = "data/endReflectivity.dat";
  double fiberSeparation = 5.0;
  
  SurfaceContainer container("data/wallReflectivity.dat");

  RectSurface* end1 = new RectSurface(0.0, 0.0, 0.0,
                                      width, 0.0, 0.0,
                                      0.0, height, 0.0,
                                      endcapReflectivity);
  container.AddSurface( end1 );

  RectSurface* end2 = new RectSurface(0.0, 0.0, length,
                                      0.0, height, length,
                                      width, 0.0, length,
                                      endcapReflectivity);
  container.AddSurface( end2 );

  RectSurface* side1 = new RectSurface(0.0, height, 0.0,
                                       width, height, 0.0,
                                       0.0, height, length,
                                       wallReflectivity);
  container.AddSurface( side1 );

  RectSurface* side2 = new RectSurface(width, height, 0.0,
                                       width, 0.0, 0.0,
                                       width, height, length,
                                       wallReflectivity);
  container.AddSurface( side2 );

  RectSurface* side3 = new RectSurface(0.0, 0.0, 0.0,
				       0.0, 0.0, length,
				       width, 0.0, 0.0,
				       wallReflectivity);
  container.AddSurface( side3 );

  RectSurface* side4 = new RectSurface(0.0, 0.0, 0.0,
                                       0.0, height, 0.0,
                                       0.0, 0.0, length,
                                       wallReflectivity);
  container.AddSurface( side4 );

  container.AddCorner(0,      0,       cornerRBot, false, true);
  container.AddCorner(width, 0,       cornerRBot, false, false);
  container.AddCorner(0,      height, cornerRTop, true, true);
  container.AddCorner(width, height, cornerRTop, true, false);

  double m = height/(width - 2*cornerR);
  double xc = width/2;
  double a = 1 + pow(m,2);
  double b = -2*xc*a;
  double c = a*pow(xc,2) - pow(0.5*(fiberSeparation + 0.07),2);
  double xf1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
  double xf2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
  double yf1 = m*(xf1 - cornerR);
  double yf2 = m*(xf2 - cornerR);
  
  cout << "Fiber1: x = " << xf1 << " y = " << yf1 << endl;
  cout << "Fiber2: x = " << xf2 << " y = " << yf2 << endl;
  container.AddFiber(xf1, yf1, 0.07/2);
  container.AddFiber(xf2, yf2, 0.07/2);


  TRandom3 gen;

  double z = 750;
  double x(0), y(0);

  Trajectory photon;

  TH2D* hCell = new TH2D("hCell", ";X;Y", 1000, -0.1*width, 1.1*width, 1000, -0.1*height, 1.1*height);

 
  for (Long64_t i = 0; i < 1e6; ++i)
    {
      x = gen.Uniform(-0.1*width, 1.1*width);
      y = gen.Uniform(-0.1*height, 1.1*height);
      photon.Init(x, y, z, 0, 0, 0, 0, 430, 1);
      if (container.IsInside(photon))
	{
	  hCell->Fill(x, y);
	}
    }

  hCell->Draw("colz");
}
