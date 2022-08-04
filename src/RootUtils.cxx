#include <TFrame.h>
#include <TMath.h>
#include <TH1.h>

#include <cmath>
#include <iostream>
#include <iomanip>

#include "include/RootUtils.hh"

using namespace std;
using namespace TMath;

void PadCoords2NDC(TVirtualPad* pad, Double_t& x, Double_t& y)
{
  Double_t dpx  = pad->GetX2() - pad->GetX1();
  Double_t dpy  = pad->GetY2() - pad->GetY1();
  Double_t xp1  = pad->GetX1();
  Double_t yp1  = pad->GetY1();
  x = (x-xp1)/dpx;
  y = (y-yp1)/dpy; 
  return;
}

void PadBounds(TVirtualPad* pad, Double_t& xmin, Double_t& xmax, Double_t& ymin, Double_t& ymax)
{
  pad->Modified();
  pad->Update();
  if (pad->GetLogx())
    {
      xmin = TMath::Power(10, pad->GetFrame()->GetX1());
      xmax = TMath::Power(10, pad->GetFrame()->GetX2());
    }
  else
    {
      xmin = pad->GetFrame()->GetX1();
      xmax = pad->GetFrame()->GetX2();
    }
  if (pad->GetLogy())
    {
      ymin = TMath::Power(10, pad->GetFrame()->GetY1());
      ymax = TMath::Power(10, pad->GetFrame()->GetY2());
    }
  else
    {
      ymin = pad->GetFrame()->GetY1();
      ymax = pad->GetFrame()->GetY2();
    }
  return;
}

Double_t Round(Double_t value, Int_t pow10)
{
  Double_t sign = (value > 0) ? 1.0 : -1.0;
  Double_t tmp_value = abs(value)*pow(10.0, -pow10);
  Double_t tmp_round = floor( tmp_value + 0.45 );
  Double_t rounded = sign*tmp_round*pow(10.0, pow10);
  return rounded;
}

TLegend* CreateLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
  TLegend* legend = new TLegend(x1, y1, x2, y2);
  legend->SetTextSize(0.04);
  legend->SetFillColor(kWhite);
  legend->SetLineColor(kWhite);
  legend->SetShadowColor(kWhite);
  return legend;
}

void UpdateStats(Double_t n, Double_t x, Double_t& mean, Double_t& sigma)
{
  Double_t M2;
  if (n < 2) M2 = 0.0;
  else M2 = (n - 2.0)*pow(sigma, 2);

  Double_t delta = x - mean;
  mean = mean + delta/n;
  M2   = M2 + delta*( x - mean );

  if (n < 2) sigma = 0.0;
  else sigma = sqrt( M2/(n - 1.0) );

  return;
}

void FillStatsHistos(TH2D* histo2D, TH1D* meanhisto, TH1D* rmshisto, Double_t ymin, Double_t ymax)
{
  meanhisto->SetStats(false);
  meanhisto->GetXaxis()->CenterTitle();
  meanhisto->SetMinimum(ymin);
  meanhisto->SetMaximum(ymax);
  rmshisto->SetStats(false);
  rmshisto->GetXaxis()->CenterTitle();
  rmshisto->SetMinimum(ymin);
  rmshisto->SetMaximum(ymax);

  for (int ibin = 0; ibin < histo2D->GetNbinsX(); ++ibin)
    {
      TH1D* projection = histo2D->ProjectionY("proj", ibin+1, ibin+1, "e");
      meanhisto->SetBinContent(ibin+1, projection->GetMean(1));
      meanhisto->SetBinError(ibin+1, projection->GetMeanError(1));
      rmshisto->SetBinContent(ibin+1, projection->GetRMS(1));
      rmshisto->SetBinError(ibin+1, projection->GetRMSError(1));
      projection->Delete();
    }

  return;
}

TStyle* GetStyle()
{
  TStyle *RootStyle = new TStyle("Root-Style","The Perfect Style for Plots ;-)");

  Int_t font_num(62); //helvetica bold

  //gSystem->Load("$ROOTSYS/lib/libPhysics.so");                                   
  //gSystem->Load("$ROOTSYS/lib/libMatrix.so");                                  
  //gSystem->Load("/usr/lib/libgsl.so");
  //gDebug = 1;

  RootStyle->SetPalette(1);
  //gROOT->SetStyle("Modern");

  RootStyle->SetPaperSize(TStyle::kUSLetter);
  
  //Canvas
  RootStyle->SetCanvasColor(0);
  RootStyle->SetCanvasBorderMode(0);
  RootStyle->SetCanvasBorderSize(10);
  RootStyle->SetCanvasBorderMode(0);
  RootStyle->SetCanvasDefH(600);
  RootStyle->SetCanvasDefW(600);
  RootStyle->SetCanvasDefX(10);
  RootStyle->SetCanvasDefY(10);
  //Canvas Background colour:
  //RootStyle->SetCanvasDefH(555);
  //RootStyle->SetCanvasDefW(450);
  //RootStyle->SetCanvasDefH(370);
  //RootStyle->SetCanvasDefW(300);

  //Pads
  RootStyle->SetPadColor       (0);
  RootStyle->SetPadBorderMode  (0);
  RootStyle->SetPadBorderSize  (2);
  RootStyle->SetPadBottomMargin(0.15);
  RootStyle->SetPadLeftMargin  (0.15);
  RootStyle->SetPadTopMargin   (0.10);
  RootStyle->SetPadRightMargin (0.07);
  RootStyle->SetPadGridX       (0);
  RootStyle->SetPadGridY       (0);
  RootStyle->SetPadTickX       (1);
  RootStyle->SetPadTickY       (1);

  // Frames                                       

  RootStyle->SetFrameBorderMode(0);
  RootStyle->SetFrameFillStyle ( 0);
  RootStyle->SetFrameFillColor ( 0);
  RootStyle->SetFrameLineColor ( 1);
  RootStyle->SetFrameLineStyle ( 0);
  RootStyle->SetFrameLineWidth ( 2);
  RootStyle->SetFrameBorderSize(10);

  //Histograms
  RootStyle->SetHistLineWidth(2);
  RootStyle->SetHistLineColor(1);  

  //Functions
  RootStyle->SetFuncColor(1);
  RootStyle->SetFuncStyle(0);
  RootStyle->SetFuncWidth(2);

  RootStyle->SetLineWidth(2);

  //Legends
  RootStyle->SetStatBorderSize(1);
  RootStyle->SetStatFont      (font_num);
  //n(name), e(entries), mM(mean), rR(rms), o(overflow), u(underflow), i(integral), sS(skewness), kK(kurtosis)
  //RootStyle->SetOptStat   ("iRMen");
  RootStyle->SetOptStat   ("rmen");
  RootStyle->SetStatColor (0);
  RootStyle->SetStatBorderSize(1);
  RootStyle->SetStatX     (0.990);
  RootStyle->SetStatY     (0.993);
  RootStyle->SetStatW     (0.25);
  RootStyle->SetStatH     (0.09);
  RootStyle->SetLegendBorderSize(0);
  RootStyle->SetLegendFillColor(0);

  // Labels, Ticks, and Titles 
  RootStyle->SetTickLength ( 0.015    ,"X");
  RootStyle->SetTitleSize  ( 0.050    ,"X");
  RootStyle->SetTitleColor ( 1        ,"X");
  //RootStyle->SetTitleOffset(1.2,"X");
  RootStyle->SetTitleOffset( 1.100    ,"X");
  RootStyle->SetLabelOffset( 0.015    ,"X");
  RootStyle->SetLabelSize  ( 0.045    ,"X");
  RootStyle->SetLabelFont  ( font_num ,"X");
  RootStyle->SetTitleFont  ( font_num ,"X");
  RootStyle->SetNdivisions ( 505      ,"X");

  RootStyle->SetTickLength ( 0.015    ,"Y");
  RootStyle->SetTitleSize  ( 0.050    ,"Y");
  RootStyle->SetTitleOffset( 1.350    ,"Y");
  RootStyle->SetLabelOffset( 0.015    ,"Y");
  RootStyle->SetLabelSize  ( 0.045    ,"Y");
  RootStyle->SetLabelFont  ( font_num ,"Y");
  RootStyle->SetTitleFont  ( font_num ,"Y");
  RootStyle->SetNdivisions ( 505      ,"Y");

  RootStyle->SetTickLength ( 0.015    ,"Z");
  RootStyle->SetTitleSize  ( 0.050    ,"Z");
  RootStyle->SetTitleOffset( 1.350    ,"Z");
  RootStyle->SetLabelOffset( 0.015    ,"Z");
  RootStyle->SetLabelSize  ( 0.045    ,"Z");
  RootStyle->SetLabelFont  ( font_num ,"Z");
  RootStyle->SetTitleFont  ( font_num ,"Z");
  RootStyle->SetNdivisions ( 505      ,"Z");
  //RootStyle->SetTitleOffset(1.5,"Z");
  //RootStyle->SetLabelFont  ( 42   ,"Z");
  //RootStyle->SetTitleFont  ( 42   ,"Z");

  RootStyle->SetTitleBorderSize  (0);
  RootStyle->SetTitleFillColor   (0);
  RootStyle->SetTitleFont        (font_num);
  RootStyle->SetTitleColor       (1);

  RootStyle->SetTextAlign(11);
  RootStyle->SetTextFont (font_num);
  RootStyle->SetTextSize (0.04);
  // Options                                                                         

  RootStyle->SetOptFit(112);
  RootStyle->SetOptTitle(0);
  //RootStyle->SetMarkerSize(0.9);
  //RootStyle->SetMarkerStyle(8);

  // global Variable: gROOT,gFile,gDirectory, gPad,gRandom,gEnv
  // gROOT->GetListOfGlobalFunctions();
  //Int_t colours[6]={1,13,9,28,35,40};   
  // 6 colours for s/w print;
  // colours are not good for LEGO2-plots
  //Int_t colours[6]={13,9,28,35,40,19};  // 6 colours for s/w print;
  //gSystem->SetIncludePath(" -I$HOME/Physics/CDF/root/include ");
  //gSystem->SetLinkedLibs("");
  // converting to HTML objects 
  // (see http://root.cern.ch/root/html303/THtml.html ):
  //gHtml = new THtml();                                                        
  //gHtml->SetOutputDir("html");
  //RootStyle->SetStyle("Plain");
  //---Style defintions: ------
  //RootStyle->SetPalette(6,colours);
  //report->SetNumberContours(6); // = no. of colours

  /*
  // get different angle for LEGO plot:
  gPad->SetTheta(theta); // default is 30
  gPad->SetPhi(phi); // default is 30
  gPad->Update();
  */

  return RootStyle;
}
