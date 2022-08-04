#ifndef SIMPLEFIBERSIM_HH
#define SIMPLEFIBERSIM_HH

#include <TRandom2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TGraph.h>
#include <TString.h>
#include <string>
#include <vector>
#include "include/SurfaceContainer.hh"
#include "include/Trajectory.hh"

class SimpleFiberSim
{
public:
  SimpleFiberSim(double zCut);
  void GeneratePhoton( double x, double y, double z, double t, double wavelength);
  bool TraceOnePhoton( double z );
  Trajectory& GetPhoton() { return _photon; }
  double Reflectivity(double idx1, double idx2, double theta, bool isS);
  double DrawWLS(double startWavelength);
  double GetWLSIterations() {return _wlsIterations;}
  TH2D*  GetDtDz_Dz() {return _hDtDz_Dz;}
  TH2D*  GetBeta_Dz() {return _hBeta_Dz;}
  TH2D*  GetvEff_Dz() {return _hvEff_Dz;}
  std::vector<float> GetX() {return _x;}
  std::vector<float> GetY() {return _y;}
  std::vector<float> GetZ() {return _z;}
  TTree* GetTree() {return _tree;}
  
private:
  TRandom2          _gen;
  SurfaceContainer  _container;
  Trajectory        _photon;

  int               _wlsIterations;
  TH1D*             _spectrum;
  TH1D*             _emission;
  TGraph*           _absWLS;
  TGraph*           _absPS;
  TGraph*           _absPMMA;
  TGraph*           _absFP;
  double            _zCut;
  double            _wlsTimeConstant;

  TH2D*             _hTheta_Phi;
  
  TH2D*             _hDtDz_Dz;
  TH2D*             _hBeta_Dz;
  TH2D*             _hvEff_Dz;

  double            _fiberNx;
  double            _fiberNy;
  double            _fiberNz;

  std::vector<float> _x;
  std::vector<float> _y;
  std::vector<float> _z;

  double             _dT;
  double             _dTdZ;
  double             _vEff;
  double             _beta;
  
  TTree*             _tree;
};

#endif

