#ifndef __PULSESHAPER__H
#define __PULSESHAPER__H

#include <TRandom.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <list>
#include <vector>

//Clocktick:                    62.5
//NumClockticksInSpill:         8800
//NumClockticksPerDigitization: 8
//ADCMaxPE:                     1807
//ASICIntegratorTime:           83300
//ASICRiseTime:                 432
//ASICFallTime:                 7000
//VaryBaseline:                 true
//BaselineMean:                 411.5
//BaselineSigma:                88.40
//NoiseSigma:
//Phase:                        [0-7]
//MaxStep:

class PulseShaper 
{
public:
  PulseShaper(float IInt, float ISlope, float RInt, float RSlope, float FInt, float FSlope);

  void Reset();
  double Step(float xNext, float dt);
  double NextPreAmp()  {return fPre[4];}
  double NextRise()    {return fRise[2];}
  double NextRiseDot() {return fRiseDot[2];}

  void SetNumClockticksInSpill(int NumClockticksInSpill)                 { fNumClockticksInSpill = NumClockticksInSpill; }
  void SetClocktick(float Clocktick)                                     { fClocktick = Clocktick; }
  void SetNumClockticksPerDigitization(int NumClockticksPerDigitization) { fNumClockticksPerDigitization = NumClockticksPerDigitization; }
  void SetPhase( int Phase )                                             { fPhase = Phase; }
  void SetMaxStep( float MaxStep )                                       { fMaxStep = MaxStep; }
  void SetBaseline( float Baseline )                                     { fBaseline = Baseline; }

  std::vector<float> CreateTrace( std::list< std::pair<float,float> >& ADCPulses ); 

private:
  float fIInt;
  float fISlope;
  float fRInt;
  float fRSlope;
  float fFInt;
  float fFSlope;

  int    fNumClockticksInSpill;
  float  fClocktick;
  int    fNumClockticksPerDigitization;
  int    fPhase;
  float  fMaxStep;
  float  fBaseline;

  double fnsPerDigitization;
  double fnsPerSpill;

  float fDt;
  float fTime;
  std::vector<float> fInput;
  std::vector<float> fPre;
  std::vector<float> fRise;
  std::vector<float> fRiseDot;
  std::vector<float> fShape;
};

#endif
