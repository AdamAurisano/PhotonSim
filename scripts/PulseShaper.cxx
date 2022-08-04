#include <algorithm>
#include <TRandom.h>
#include <vector>
#include <TMath.h>
#include <TRandom3.h>
#include <TF1.h>
#include <cmath>
#include <TCanvas.h>
#include <TGraph.h>
#include <iostream>
#include <iomanip>
#include <TLine.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>

#include "PulseShaper.hh"

using namespace std;
using namespace TMath;

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

bool CompPulses( pair<float, float> first, pair<float, float> second)
{
  if (first.first < second.first) return true;
  else if (first.first == second.first && first.second > second.second) return true;
  else return false;
}

bool IsSagged( pair<float, float> first, pair<float, float> second ) 
{
  return first.first == second.first;
}

PulseShaper::PulseShaper(float IInt, float ISlope, float RInt, float RSlope, float FInt, float FSlope)
  : fIInt(IInt), fISlope(ISlope), fRInt(RInt), fRSlope(RSlope), fFInt(FInt), fFSlope(FSlope), 
    fNumClockticksInSpill(8800), fClocktick(62.5), fNumClockticksPerDigitization(8),
    fPhase(0), fMaxStep(6.25), fBaseline(250.0), fTime(0.0)
{
  fInput.resize(2);
  fPre.resize(7);
  fRise.resize(4);
  fRiseDot.resize(3);
  fShape.resize(2);
  Reset();
}

void PulseShaper::Reset()
{
  for (int i = 0; i < 2; ++i) fInput[i]   = 0.0;
  for (int i = 0; i < 7; ++i) fPre[i]     = 0.0;
  for (int i = 0; i < 4; ++i) fRise[i]    = 0.0;
  for (int i = 0; i < 3; ++i) fRiseDot[i] = 0.0;
  for (int i = 0; i < 2; ++i) fShape[i]   = 0.0;
  //trace.resize(fNumClockticksInSpill);
  //for (int i = fPhase; i < fNumClockticksInSpill; i += fNumClockticksPerDigitization)
  //{
  //  if (i%fNumClockticksPerDigitization == fPhase) trace[i] = fBaseline;
  //  else trace[i] = 0.0;
  //}
}

double PulseShaper::Step(float xNext, float dt)
{
  fInput[0]   = fInput[1];
  fPre[0]     = fPre[4];
  fRise[0]    = fRise[2];
  fRiseDot[0] = fRiseDot[2];
  fShape[0]   = fShape[1];
  //set the new time step and the next impulse
  fDt = dt;
  fTime += fDt;
  fInput[1]   = xNext;

  double IInv(0.0);
  double w1(0.0), w2(0.0), w3(0.0), w4(0.0);
  IInv = (fPre[0] > 0) ? 1.0/(fIInt + fISlope*fPre[0]) : 1.0/fFInt;  
  fPre[1] = fPre[0] - 0.25*fDt*fPre[0]*IInv;
  IInv = (fPre[1] > 0) ? 1.0/(fIInt + fISlope*fPre[1]) : 1.0/fFInt;  
  fPre[2] = fPre[1] - 0.25*fDt*fPre[1]*IInv;
  IInv = (fPre[2] > 0) ? 1.0/(fIInt + fISlope*fPre[2]) : 1.0/fFInt;  
  fPre[3] = fPre[2] - 0.25*fDt*fPre[2]*IInv;
  IInv = (fPre[3] > 0) ? 1.0/(fIInt + fISlope*fPre[3]) : 1.0/fFInt;  
  fPre[4] = fPre[3] - 0.25*fDt*fPre[3]*IInv + fInput[1];
  IInv = (fPre[4] > 0) ? 1.0/(fIInt + fISlope*fPre[4]) : 1.0/fFInt;  
  fPre[5] = fPre[4] - 0.25*fDt*fPre[4]*IInv;
  IInv = (fPre[5] > 0) ? 1.0/(fIInt + fISlope*fPre[5]) : 1.0/fFInt;  
  fPre[6] = fPre[5] - 0.25*fDt*fPre[5]*IInv;

  double RInv(0.0);
  RInv = (fPre[0] > 0) ? 1.0/(fRInt + fRSlope*fPre[0]) : 1.0/fRInt;    
  w1 = 0.5*fDt*(fPre[0] - fRise[0])*RInv;
  RInv = (fPre[1] > 0) ? 1.0/(fRInt + fRSlope*fPre[1]) : 1.0/fRInt;    
  w2 = 0.5*fDt*(fPre[1] - fRise[0] - 0.5*w1)*RInv;
  w3 = 0.5*fDt*(fPre[1] - fRise[0] - 0.5*w2)*RInv;
  RInv = (fPre[2] > 0) ? 1.0/(fRInt + fRSlope*fPre[2]) : 1.0/fRInt;    
  w4 = 0.5*fDt*(fPre[2] - fRise[0] - w3)*RInv;
  fRise[1] = fRise[0] + (1.0/6.0)*(w1 + 2*w2 + 2*w3 + w4);

  RInv = (fPre[2] > 0) ? 1.0/(fRInt + fRSlope*fPre[2]) : 1.0/fRInt;    
  w1 = 0.5*fDt*(fPre[2] - fRise[1])*RInv;
  RInv = (fPre[3] > 0) ? 1.0/(fRInt + fRSlope*fPre[3]) : 1.0/fRInt;    
  w2 = 0.5*fDt*(fPre[3] - fRise[1] - 0.5*w1)*RInv;
  w3 = 0.5*fDt*(fPre[3] - fRise[1] - 0.5*w2)*RInv;
  RInv = (fPre[4] > 0) ? 1.0/(fRInt + fRSlope*fPre[4]) : 1.0/fRInt;    
  w4 = 0.5*fDt*(fPre[4] - fRise[1] - w3)*RInv;
  fRise[2] = fRise[1] + (1.0/6.0)*(w1 + 2*w2 + 2*w3 + w4);

  RInv = (fPre[4] > 0) ? 1.0/(fRInt + fRSlope*fPre[4]) : 1.0/fRInt;
  w1 = 0.5*fDt*(fPre[4] - fRise[2])*RInv;
  RInv = (fPre[5] > 0) ? 1.0/(fRInt + fRSlope*fPre[5]) : 1.0/fRInt;
  w2 = 0.5*fDt*(fPre[5] - fRise[2] - 0.5*w1)*RInv;
  w3 = 0.5*fDt*(fPre[5] - fRise[2] - 0.5*w2)*RInv;
  RInv = (fPre[6] > 0) ? 1.0/(fRInt + fRSlope*fPre[6]) : 1.0/fRInt;
  w4 = 0.5*fDt*(fPre[6] - fRise[2] - w3)*RInv;
  fRise[3] = fRise[2] + (1.0/6.0)*(w1 + 2*w2 + 2*w3 + w4);

  fRiseDot[1] = (fRise[2] - fRise[0])/fDt;
  fRiseDot[2] = (fRise[3] - fRise[1])/fDt;

  double FInv0 = (fPre[0] > 0) ? 1.0/(fFInt + fFSlope*fPre[0]) : 1.0/fFInt;
  double FInv1 = (fPre[2] > 0) ? 1.0/(fFInt + fFSlope*fPre[2]) : 1.0/fFInt;
  double FInv2 = (fPre[4] > 0) ? 1.0/(fFInt + fFSlope*fPre[4]) : 1.0/fFInt;
  w1 = fDt*(fRiseDot[0] - fShape[0]*FInv0);
  w2 = fDt*(fRiseDot[1] - (fShape[0] + 0.5*w1)*FInv1);
  w3 = fDt*(fRiseDot[1] - (fShape[0] + 0.5*w2)*FInv1);
  w4 = fDt*(fRiseDot[2] - (fShape[0] + w3)*FInv2);
  fShape[1] = fShape[0] + (1.0/6.0)*(w1 + 2*w2 + 2*w3 + w4);
  
  //cout << "Time: " << fTime << " Value: " << fShape[1] << endl;
  
  return fShape[1];
}

std::vector<float> PulseShaper::CreateTrace( std::list< std::pair<float,float> >& ADCPulses )
{
  ADCPulses.sort(CompPulses);
  ADCPulses.unique(IsSagged);

  std::vector<float> trace(fNumClockticksInSpill, 0.0);
  for (int i = fPhase; i < fNumClockticksInSpill; i += fNumClockticksPerDigitization) trace[i] = fBaseline;
  
  fnsPerDigitization = fNumClockticksPerDigitization*fClocktick;
  fnsPerSpill = fNumClockticksInSpill*fClocktick;
  
  Reset();
  if (ADCPulses.size() == 0) return trace;

  int iDigitization(0);
  double firstPulseTime        = ADCPulses.front().first;
  double firstDigitizationTime = fClocktick*fPhase;

  fTime = (firstPulseTime < firstDigitizationTime) ? firstPulseTime : firstDigitizationTime;
  fTime -= fClocktick;
  
  while (fTime < fnsPerSpill) {
    double timeToNextInput(9999999.0);
    if (ADCPulses.size() > 0) timeToNextInput = ADCPulses.front().first - fTime;
    double timeToNextDigitization = (iDigitization*fNumClockticksPerDigitization + fPhase)*fClocktick - fTime;

    double maxStep = fMaxStep;
    if (std::fabs(fInput[1]) == 0)
      {
	double delta = std::fabs(NextRiseDot());
	double step(0.0);
	if (delta > 0)
	  {
	    step = 1.0/delta;
	    if (step < fMaxStep) maxStep = step;
	  }
      }

    //there is neither a clocktick nor an input before we hit the maximum step
    if ( (maxStep < timeToNextInput) && (maxStep < timeToNextDigitization) ) {
      Step(0.0, maxStep);
    }

    //a clocktick and the input are at about the same time
    else if ( std::fabs(timeToNextInput - timeToNextDigitization) < 1e-8 ) {
      double input = ADCPulses.front().second;
      double output = Step(input, timeToNextDigitization);
      ADCPulses.pop_front();
      if ( iDigitization*fNumClockticksPerDigitization + fPhase < fNumClockticksInSpill) {
	trace[iDigitization*fNumClockticksPerDigitization + fPhase] += output;
      }
      ++iDigitization;
    }

    //there is an input before the next clocktick
    else if (timeToNextInput < timeToNextDigitization) {
      double input = ADCPulses.front().second;
      Step(input, timeToNextInput);
      ADCPulses.pop_front();
    }

    //the clocktick is before the input
    else {
      double output = Step(0.0, timeToNextDigitization);
      if ( iDigitization*fNumClockticksPerDigitization + fPhase < fNumClockticksInSpill) {
	trace[iDigitization*fNumClockticksPerDigitization + fPhase] += output;
      }
      ++iDigitization;
    }
  }

  //cout << "Trace: " << endl;
  //for (int i = 0; i < fNumClockticksInSpill; ++i) cout << trace[i] << endl;
  return trace;
}
