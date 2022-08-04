#ifndef ROOTUTILS_HH
#define ROOTUTILS_HH

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

#include <vector>

void PadCoords2NDC(TVirtualPad* pad, Double_t& x, Double_t& y);

//finds the upper and lower bounds of the y axis
void PadBounds(TVirtualPad* pad, Double_t& xmin, Double_t& xmax, Double_t& ymin, Double_t& ymax);

//rounds value to the place specified in pow10 (ie if pow10 = 2 --> round to the nearest hundred, -2 --> nearest hundreth)
Double_t Round(Double_t value, Int_t pow10);

//Creates a legend with the right starting properties
TLegend* CreateLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2); 

//Updates a running mean and RMS
void UpdateStats(Double_t n, Double_t x, Double_t& mean, Double_t& sigma);

//Updates histos with mean and RMS
void FillStatsHistos(TH2D* histo2D, TH1D* meanhisto, TH1D* rmshisto, Double_t ymin, Double_t ymax);

//Returns a style good for publications
TStyle* GetStyle();

#endif
