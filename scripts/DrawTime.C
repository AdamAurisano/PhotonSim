#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"

#include <vector>
#include <algorithm>
#include <numeric>

using namespace std;

#include "PulseShaper.hh"

template<class T> inline T sqr(T x){return x*x;}

// Number of points we expect to receive, to find the best peak match with
const unsigned int kNumFineTimingADCPoints = 4; // Including the DCS origin
// Empirically, the best fit is never found outside this range. +3 is an
// upper limit on the positive side because it means the pulse starts beyond
// the last ADC observation.
const double kZeroOffsetSamples = -1;
// Have the capacity to consider 65,536 offsets, which is way too many to
// give good performance.
const int kOffsetStep = 16;
const double kSamplesPerOffset = 1./(256*64);
// With 4096 steps that means the maximum offset is +3 and step size is
// 0.49ns (in the FD).

double GetExpectations(double t0,
		       double riseTime, double fallTime, double preAmp,
		       const int16_t* obs,
		       double* exps, double& base, int mode)
{
  // There are only a few thousand distinct inputs to this function per
  // combination of riseTime, fallTime and preAmp. Cache all possible
  // outputs.  This saves about 40% of the run time of CalHit for ND and FD
  // jobs, for which this cache needs to be built only once.  For TB jobs,
  // two kinds of FEBs are present, so the cache is rebuilt many times as we
  // see hits from each kind.  This is *still* faster than not having a
  // cache, although obviously a little dumb since we could save both sets of
  // cached data with a little more work.
  // The calculation of cache_size is slightly fragile because
  // kSamplesPerOffset and kZeroOffsetSamples are defined as doubles even
  // though they are representing integer concepts (also the meaning of
  // kSamplesPerOffset is the reciprocal of its name).  But I think the rest
  // of the code also breaks if kZeroOffsetSamples isn't an integer or if
  // kSamplesPerOffset isn't the reciprocal of an integer.
  const int steps_per_offset = 1 / kSamplesPerOffset / kOffsetStep;
  const int num_offsets = (kNumFineTimingADCPoints - 1) - kZeroOffsetSamples;
  const int cache_size = steps_per_offset * num_offsets;
  
  static double cache_ideal[cache_size] = {0};
  static double cache_3param[cache_size] = {0};
  static double cache_riseTime = -1, cache_fallTime = -1, cache_preAmp = -1;
  if(cache_ideal[0] == 0 || riseTime != cache_riseTime ||
     fallTime != cache_fallTime || preAmp != cache_preAmp){
    for(int i = 0; i < cache_size; i++){
      const double t_minus_t0 = double(i+1)/steps_per_offset;
      // This is the ideal ASIC formula from SimpleReadout
      cache_ideal[i] = exp(-t_minus_t0/fallTime)*(1-exp(-t_minus_t0/riseTime));
      //new three parameters function
      cache_3param[i] = (fallTime * preAmp)*(
					     (exp(-t_minus_t0/fallTime) / ((preAmp-fallTime)*(fallTime-riseTime))) -
					     (exp(-t_minus_t0/preAmp)   / ((preAmp-fallTime)*(preAmp - riseTime))) +
					     (exp(-t_minus_t0/riseTime) / ((preAmp-riseTime)*(riseTime-fallTime))) );
    }
  }
  cache_riseTime = riseTime;
  cache_fallTime = fallTime;
  cache_preAmp = preAmp;
  
  for(unsigned int t = 0; t < kNumFineTimingADCPoints; ++t){
    if(t <= t0){
      exps[t] = 0;
    }
    else{
      const int cachei = (t-t0) * steps_per_offset - 1;
      exps[t] = mode == 0?cache_ideal[cachei]: cache_3param[cachei];
    }
  }
  
  // We can calculate the best fit baseline and normalization at this offset
  // analytically (just write the chisq expression, and take derivatives with
  // norm or base), like so:
  double U, V, W, X, Z;
  U = V = W = X = Z = 0;
  for(unsigned int i = 0; i < kNumFineTimingADCPoints; ++i){
    U += exps[i];
    V += 1;
    W += obs[i];
    X += exps[i]*exps[i];
    Z += obs[i]*exps[i];
  }
  
  const double norm = (W*U-V*Z)/(U*U-X*V);
  base = (W*X-Z*U)/(X*V-U*U);
  
  // Correct all the expectations for the best guess baseline and norm
  for(unsigned int n = 0; n < kNumFineTimingADCPoints; ++n){
    exps[n] = norm*exps[n]+base;
  }
  return norm;
}
  
//......................................................................
uint16_t ADCShapeFit(int16_t adc1, int16_t adc2, int16_t adc3,
		     double riseTime, double fallTime, double preAmp,
		     bool& goodTime, int fMode)
{
  // The observations, in chronological order
  const int16_t obs[kNumFineTimingADCPoints] = {0, adc1, adc2, adc3};
  double bestchisq = 1e10; // infinity
  
  // Sane default
  
  uint16_t bestoffset = -kZeroOffsetSamples/kSamplesPerOffset;
  goodTime = false;
  
  // Scan through possible offsets of the peak, looking for the best match
  uint16_t offset = 0;
  do{
    // In TDC
    const double t0 = (double(offset)*kSamplesPerOffset+kZeroOffsetSamples);
    
    // Expectations
    double exps[kNumFineTimingADCPoints];
    double junk;
    const double norm = GetExpectations(t0, riseTime, fallTime, preAmp, obs, exps, junk, fMode);
    
    if(norm < 0) continue;
    
    // Assuming equal, uncorrelated, errors on all points. (Including
    // correlations was not found to improve resolution).
    double chisq = 0;
    for(unsigned int i = 0; i < kNumFineTimingADCPoints; ++i)
      chisq += sqr(obs[i]-exps[i]);
    assert(chisq >= 0);
    if(chisq < bestchisq){
      bestchisq = chisq;
      bestoffset = offset;
      goodTime = true;
    }
  } while((offset += kOffsetStep) != 0);
  
  if(bestoffset == 0 || bestoffset == 256*256-kOffsetStep){
    bestoffset = -kZeroOffsetSamples/kSamplesPerOffset;
    goodTime = false;
  }
  return bestoffset;
}

double TNS(double tns0, int16_t adc1, int16_t adc2, int16_t adc3,
	   double& adcpeak, bool& goodTime, double& base)
{
  uint16_t offset;
  double riseTime ;
  double fallTime ;
  double preAmp;
  bool fIsND = true;
  bool fIsMC = true;
  double fMode = 0;
  
  const double sampleTime = fIsND ? 125 : 500; // ns
  
  if (fMode == 0){
    riseTime = (fIsND ?  140 :  380)/sampleTime;
    fallTime = (fIsND ? 4500 : 7000)/sampleTime;
    preAmp = 0;
    if(!fIsND && !fIsMC) riseTime = 460/sampleTime;
  } else{
    riseTime = (fIsND ? 116 : 424)/sampleTime;
    fallTime = (fIsND ? 6729 : 10449)/sampleTime;
    preAmp = ( fIsND ? 147000 : 110000)/sampleTime;
    
    if(!fIsND && !fIsMC) {
      riseTime = 410/sampleTime;
      fallTime = 9000/sampleTime;
    }
  }
  offset = ADCShapeFit(adc1, adc2, adc3, riseTime, fallTime, preAmp, goodTime, fMode);

  if(goodTime){
    // The observations, in chronological order
    const int16_t obs[kNumFineTimingADCPoints] = {0, adc1, adc2, adc3};
    const double t0 = (double(offset)*kSamplesPerOffset+kZeroOffsetSamples);
    double junk[kNumFineTimingADCPoints]; // Expectations
    const double norm = GetExpectations(t0, riseTime, fallTime, preAmp, obs, junk, base, fMode);
    
    // The value the peak of the curve would have before scaling by norm
    // Analytic calculation from the functional form in ReadoutSim
    static const double peak = pow(riseTime/(riseTime+fallTime), riseTime/fallTime)-
      pow(riseTime/(riseTime+fallTime), 1+riseTime/fallTime);
    
    adcpeak = norm*peak;
  }
  else{
    // Unless the fit failed, in which case we're probably safest just doing
    // this.
    adcpeak = adc3;
    return tns0;
  }
  
  return tns0+(offset*kSamplesPerOffset+kZeroOffsetSamples)*sampleTime;
}

/*std::vector<float>*/

float DCSTrigger(std::vector<float> trace,
		 int firstDigitizationClocktick,
		 int fNumClockticksInSpill,
		 int fNumClockticksPerDigitization,
		 double& adcpeak)
{
  double threshold = 45;
  int nSamples = 4;
  int nPretrig = 3;
  int fLookbackSamples = 3;
  double fRetriggerThreshold = 0.5;

  bool triggered = false;
  double peakDCSValue = 0;
  int peakDCSTime = 0;

  vector<float> triggers(fNumClockticksInSpill, 0);

  double DCS_Space[fNumClockticksInSpill];
  
  //const int iMin = firstDigitizationClocktick + std::max(nPretrig-fLookbackSamples, 0)*fNumClockticksPerDigitization;
  //const int iMax = fNumClockticksInSpill - (nSamples-nPretrig+fLookbackSamples-1)*fNumClockticksPerDigitization;

  const int iMin = firstDigitizationClocktick + (nSamples-nPretrig+fLookbackSamples-1)*fNumClockticksPerDigitization;
  const int iMax = fNumClockticksInSpill;

  int ipoint = 0;
  for (int iClocktick = iMin; iClocktick < iMax; iClocktick += fNumClockticksPerDigitization) {
    //DCS_Space[iClocktick] = trace[iClocktick + fLookbackSamples*fNumClockticksPerDigitization] - trace[iClocktick];
    DCS_Space[iClocktick] = trace[iClocktick] - trace[iClocktick - fLookbackSamples*fNumClockticksPerDigitization];
    ++ipoint;
    double value = DCS_Space[iClocktick];
    //if (trace[iClocktick] > 0) cout << trace[iClocktick] << endl;;
    //cout << trace[iClocktick] << endl;;

    if (value > threshold) triggered = true;
    else if (triggered && value < peakDCSValue)
      {
	double adc[nSamples];
	for (int iSample = 0; iSample < nSamples; ++iSample) {
	  //cout << "adc" << iSample << ": " << trace[peakDCSTime - (nPretrig - iSample)*fNumClockticksPerDigitization] << endl;
	  adc[iSample] = trace[peakDCSTime - (nPretrig - iSample)*fNumClockticksPerDigitization];
	}
	
	bool goodTime = false;
	double base(0);
	adcpeak = 0;
	
	base = adc[0];
	double tns0 = peakDCSTime*62.5;
	
	double tns = TNS(tns0, adc[1] - adc[0], adc[2] - adc[0], adc[3] - adc[0], adcpeak, goodTime, base);
	//cout << "tns0 = " << tns0 << " tns = " << tns << " adcpeak = " << adcpeak << " base = " << base << endl;

	if (goodTime) return tns;
	else return -999.0;
	
       	triggered = false;
	peakDCSValue = 0;
	peakDCSTime  = 0;
      }
    // looking for peak?
    if (triggered) {
      //cout << "tick: " << iClocktick << " value: " << value << endl;
      if (value > peakDCSValue) {
	peakDCSTime  = iClocktick;
	peakDCSValue = value;
      }
    }
  }
  adcpeak = 0;
  return -999;
  //return triggers;
}

double Atten(double x)
{
  //return 0.6245*exp(-x/322.3) + 0.3755*exp(-x/902.9);
  return 0.4624*exp(-x/175.1) + 0.5376*exp(-x/687.2);
}

void DrawTime(double startZ, int ifile)
{
  TRandom3 gen;
  gen.SetSeed(0);

  /*
  const float clocktick = 62.5;
  const float step  = 1.0; //ns
  const float RInt   = 420.;
  const float RSlope = 0.0;
  const float FInt   = 10373;
  const float FSlope = 2.4;
  const float IInt   = 112557;
  const float ISlope = 0.0;
  int phase = 0;
  const int nclockticks = 8800;
  const int nperdigitization = 8;
  const double baseline = 400;
  */
  
  const float clocktick = 62.5;
  const float step  = 1.0; //ns
  const float RInt   = 107.;
  const float RSlope = 0.0;
  const float FInt   = 6728;
  const float FSlope = 0.9;
  const float IInt   = 119278;
  const float ISlope = 0.0;
  int phase = 0;
  const int nclockticks = 8800;
  const int nperdigitization = 2;
  const double baseline = 250;

  double pe2adc = 1.0/1705.0*( (1<<12) - 1 );
  
  PulseShaper shaper(IInt, ISlope, RInt, RSlope, FInt, FSlope);
  shaper.SetPhase(phase);
  shaper.SetClocktick(62.5);
  shaper.SetMaxStep(step);
  shaper.SetBaseline(baseline);
  shaper.SetNumClockticksInSpill(nclockticks);
  shaper.SetNumClockticksPerDigitization(nperdigitization);

  //const double length(1500); //length of cell
  const double length(399.28 - 0.5); //length of cell

  TH1D* hOldTime = new TH1D("hOldTime", ";Photon Time (ns)", 500, 0, 2000);
  TH1D* hFiberTime = new TH1D("hFiberTime", ";Photon Time (ns)", 500, 0, 2000);
  TH1D* hFiberScintTime = new TH1D("hFiberScintTime", ";Photon Time (ns)", 500, 0, 2000);
  
  TH1D* hFiberPulseRMS      = new TH1D("hFiberPulseRMS", ";Pulse RMS (ns)", 400, 0, 5000);
  TH1D* hFiberScintPulseRMS = new TH1D("hFiberScintPulseRMS", ";Pulse RMS (ns)", 400, 0, 5000);
  TH1D* hOldPulseRMS        = new TH1D("hOldPulseRMS", ";Pulse RMS (ns)", 400, 0, 5000);
  
  TH1D* hFiberPulseWidth      = new TH1D("hFiberPulseWidth", ";Pulse Width (ns)", 400, 0, 5000);
  TH1D* hFiberScintPulseWidth = new TH1D("hFiberScintPulseWidth", ";Pulse Width (ns)", 400, 0, 5000);
  TH1D* hOldPulseWidth        = new TH1D("hOldPulseWidth", ";Pulse Width (ns)", 400, 0, 5000);

  TH1D* hFiberPulseTime      = new TH1D("hFiberPulseTime", ";Pulse Time (ns)", 200, -5000, 5000);
  TH1D* hFiberScintPulseTime = new TH1D("hFiberScintPulseTime", ";Pulse Time (ns)", 200, -5000, 5000);
  TH1D* hOldPulseTime        = new TH1D("hOldPulseTime", ";Pulse Time (ns)", 200, -5000, 5000);

  TFile* fCollection = new TFile("dT_dZ_CollectionRate.root");
  TH2D* hCollection = (TH2D*)fCollection->Get("dT_dZ_CollectionRate");
  
  TFile* fSpread = new TFile("FiberSpread.root");
  TH2D* hSpread = (TH2D*)fSpread->Get("FiberSpread");

  TFile* fOldSpread = new TFile("Dt_per_z_distribution.root");
  TH1D* hOldSpread = (TH1D*)fOldSpread->Get("Dt_per_z");
  
  double v0 = TMath::Ccgs()*1e-9/1.59;
  int nSlice = 100;
  //make individual traces
  for (int i = 0; i < 100; ++i)
    {
      //if (i%10 == 0) cout << "iteration: " << i << endl;
      cout << "iteration: " << i << endl;
      vector<double> vecOldTNS;
      vector<double> vecFiberTNS;
      vector<double> vecFiberScintTNS;
      //draw 50 hits and find the delta T between the earliest and latest hits
      //nSlice = (int)gen.Uniform(50, 140);
      //double sliceExtent = gen.Uniform(20, 240);
      for (int iSlice = 0; iSlice < nSlice; ++iSlice)
	{
	  //startZ = gen.Uniform(750 - 0.5*sliceExtent, 750 + 0.5*sliceExtent);
	  startZ = 0.5*length;
	  phase = gen.Integer(8);
	  shaper.SetPhase( phase );
	  
	  list< pair<float, float> > inputsFiber;
	  list< pair<float, float> > inputsFiberScint;
	  list< pair<float, float> > inputsOld;

	  int view = gen.Integer(1);
	  double factor = (view == 0) ? 0.57861 : 0.57023;
	  factor *= 0.13*3151.04*1.05*0.85;
	  int nPho = gen.Poisson( gen.Gaus(17.956, 0.002)*factor );
	  inputsFiber.clear();
	  inputsFiberScint.clear();
	  inputsOld.clear();
	  for (int iPho = 0; iPho < nPho; ++iPho)
	    {
	      double dZ(0.0), dT(0.0);
	      hCollection->GetRandom2(dZ, dT);
	      vector<double> lengths;
	      double Dsp = startZ + dZ;
	      double Dlp = 2*length - Dsp;
	      if (Dsp > 0 && Dsp < length)
		{
		  lengths.push_back(Dsp);
		  lengths.push_back(Dlp);
		}
	      double Dsm = startZ - dZ;
	      double Dlm = 2*length - Dsp;
	      if (Dsm > 0 && Dsm < length)
		{
		  lengths.push_back(Dsm);
		  lengths.push_back(Dlm);
		}
		      
	      double q[]   = {0.933, 0.024, 0.022, 0.0213};
	      double tau[] = {3.95, 23.56, 78.86, 546.39};
	      
	      double prob(0);
	      double scintTime(0.0);
	      double u = gen.Rndm();
	      if (u < 0.05) scintTime = gen.Exp(3.6);
	      else
		{
		  u = gen.Rndm();
		  for (int i = 0; i < 4; ++i)
		    {
		      prob += q[i];
		      if (u < prob)
			{
			  scintTime = gen.Exp(tau[i]);
			  break;
			}
		    }
		}
	      
	      for (unsigned int iL = 0; iL < lengths.size(); ++iL)
		{
		  double attenWeight = Atten(lengths[iL])*0.25;
		  //double time = 2000 + dT;
		  double time = 100 + dT;
		  time += lengths[iL]/v0;
		  int iBin = hSpread->GetXaxis()->FindBin(lengths[iL]);
		  TH1D* hDelay = hSpread->ProjectionY("hDelay", iBin, iBin);

		  double fiberTime = time + hDelay->GetRandom() + gen.Exp(3.6);
		  double fiberScintTime  = time + hDelay->GetRandom() + scintTime;
		  double oldTime = time + hOldSpread->GetRandom()*lengths[iL] + gen.Exp(9.0);
		  
		  inputsFiber.push_back( make_pair<float, float>( fiberTime, attenWeight*pe2adc ) );

		  inputsFiberScint.push_back( make_pair<float, float>( fiberScintTime, attenWeight*pe2adc ) );
		  
		  inputsOld.push_back( make_pair<float, float>( oldTime, attenWeight*pe2adc ) );

		  hOldTime->Fill(oldTime, attenWeight*pe2adc);
		  hFiberTime->Fill(fiberTime, attenWeight*pe2adc);
		  hFiberScintTime->Fill(fiberScintTime, attenWeight*pe2adc);
		  
		  //cout << setw(15) << "old time" << setw(15) << "fiber time" << setw(15) << "fiber+scint" << endl;
		  //cout << setw(15) << oldTime << setw(15) << fiberTime << setw(15) << fiberScintTime << endl;
		}
	    }

	  inputsFiber.sort();

	  double totPho = 0;
	  for (list<pair<float,float> >::iterator it = inputsFiber.begin(); it != inputsFiber.end(); it++) totPho += it->second;
	  //cout << "totPho = " << totPho << endl;


	  vector<float> traceFiber = shaper.CreateTrace( inputsFiber );
	  //cout << "traceFiber.size() = " << traceFiber.size() << endl;
	  double adcFiber(0);
	  double tnsFiber = DCSTrigger(traceFiber, phase, nclockticks, nperdigitization, adcFiber);
	  if (tnsFiber > 0) vecFiberTNS.push_back(tnsFiber);
	  
	  inputsFiberScint.sort();
	  vector<float> traceFiberScint = shaper.CreateTrace( inputsFiberScint );
	  //cout << "traceFiberScint.size() = " << traceFiberScint.size() << endl;
	  double adcFiberScint(0);
	  double tnsFiberScint = DCSTrigger(traceFiberScint, phase, nclockticks, nperdigitization, adcFiberScint);
	  if (tnsFiberScint > 0) vecFiberScintTNS.push_back(tnsFiberScint);
	  
	  inputsOld.sort();
	  vector<float> traceOld = shaper.CreateTrace( inputsOld );
	  //cout << "traceOld.size() = " << traceOld.size() << endl;
	  double adcOld(0);
	  double tnsOld = DCSTrigger(traceOld, phase, nclockticks, nperdigitization, adcOld);
	  if (tnsOld > 0) vecOldTNS.push_back(tnsOld);

	  /*
	  cout << setw(15) << "   " << setw(15) << "Old"  << setw(15) << "Fiber"  << setw(15) << "Fiber+Scint" << endl;
	  cout << setw(15) << "ADC" << setw(15) << adcOld << setw(15) << adcFiber << setw(15) << adcFiberScint << endl;
	  cout << setw(15) << "TNS" << setw(15) << tnsOld << setw(15) << tnsFiber << setw(15) << tnsFiberScint << endl;
	  */
	}
      
      sort(vecFiberTNS.begin(), vecFiberTNS.end());
      sort(vecFiberScintTNS.begin(), vecFiberScintTNS.end());
      sort(vecOldTNS.begin(), vecOldTNS.end());

      double meanFiber = accumulate(vecFiberTNS.begin(), vecFiberTNS.end(), 0.0)/vecFiberTNS.size();
      double SqSumFiber = inner_product(vecFiberTNS.begin(), vecFiberTNS.end(), vecFiberTNS.begin(), 0.0);
      double StDevFiber = sqrt( SqSumFiber/vecFiberTNS.size() - meanFiber*meanFiber);
	
      double meanFiberScint = accumulate(vecFiberScintTNS.begin(), vecFiberScintTNS.end(), 0.0)/vecFiberScintTNS.size();
      double SqSumFiberScint = inner_product(vecFiberScintTNS.begin(), vecFiberScintTNS.end(), vecFiberScintTNS.begin(), 0.0);
      double StDevFiberScint = sqrt( SqSumFiberScint/vecFiberScintTNS.size() - meanFiberScint*meanFiberScint);
      
      double meanOld = accumulate(vecOldTNS.begin(), vecOldTNS.end(), 0.0)/vecOldTNS.size();
      double SqSumOld = inner_product(vecOldTNS.begin(), vecOldTNS.end(), vecOldTNS.begin(), 0.0);
      double StDevOld = sqrt( SqSumOld/vecOldTNS.size() - meanOld*meanOld);
      
      for (int iSlice = 0; iSlice < vecFiberTNS.size(); ++iSlice) hFiberPulseTime->Fill( vecFiberTNS[iSlice] - meanFiber);
      for (int iSlice = 0; iSlice < vecFiberScintTNS.size(); ++iSlice) hFiberScintPulseTime->Fill( vecFiberScintTNS[iSlice] - meanFiberScint);
      for (int iSlice = 0; iSlice < vecOldTNS.size(); ++iSlice) hOldPulseTime->Fill( vecOldTNS[iSlice] - meanOld);

      hFiberPulseWidth->Fill( vecFiberTNS[vecFiberTNS.size() - 1] - vecFiberTNS[0]);
      hFiberScintPulseWidth->Fill( vecFiberScintTNS[vecFiberScintTNS.size() - 1] - vecFiberScintTNS[0]);
      hOldPulseWidth->Fill( vecOldTNS[vecOldTNS.size() - 1] - vecOldTNS[0]);

      hFiberPulseRMS->Fill(StDevFiber);
      hFiberScintPulseRMS->Fill(StDevFiberScint);
      hOldPulseRMS->Fill(StDevOld);

      cout << setw(15) << ""
	   << setw(15) << "Old"
	   << setw(15) << "Fiber"
	   << setw(15) << "Fiber+Scint"
	   << endl;

      cout << setw(15) << "Width"
	   << setw(15) << vecOldTNS[vecOldTNS.size() - 1] - vecOldTNS[0]
	   << setw(15) << vecFiberTNS[vecFiberTNS.size() - 1] - vecFiberTNS[0]
	   << setw(15) <<  vecFiberScintTNS[vecFiberScintTNS.size() - 1] - vecFiberScintTNS[0]
	   << endl;

      cout << setw(15) << "RMS"
	   << setw(15) << StDevOld
	   << setw(15) << StDevFiber
	   << setw(15) << StDevFiberScint
	   << endl;
      
    }

  /*
  TCanvas* cPulseTime = new TCanvas("cPulseTime","cPulseTime");
  cPulseTime->Divide(2,1);
  cPulseTime->cd(1);
  hOldPulseTime->SetLineColor(kBlue);
  hOldPulseTime->SetLineWidth(2);
  hOldPulseTime->Draw("hist");

  cPulseTime->cd(2);
  hPulseTime->SetLineColor(kRed);
  hPulseTime->SetLineWidth(2);
  hPulseTime->Draw("hist");

  TCanvas* cPulseWidth = new TCanvas("cPulseWidth","cPulseWidth");
  cPulseWidth->Divide(2,1);
  cPulseWidth->cd(1);
  hOldPulseWidth->SetLineColor(kBlue);
  hOldPulseWidth->SetLineWidth(2);
  hOldPulseWidth->Draw("hist");

  cPulseWidth->cd(2);
  hPulseWidth->SetLineColor(kRed);
  hPulseWidth->SetLineWidth(2);
  hPulseWidth->Draw("hist");
  */
  
  TString fname = "TimeSim_";
  fname += ifile;
  fname += ".root";
  TFile* file = new TFile(fname, "recreate");
  file->WriteTObject(hFiberPulseTime);
  file->WriteTObject(hFiberScintPulseTime);
  file->WriteTObject(hOldPulseTime);
  file->WriteTObject(hFiberPulseWidth);
  file->WriteTObject(hFiberScintPulseWidth);
  file->WriteTObject(hOldPulseWidth);
  file->WriteTObject(hFiberPulseRMS);
  file->WriteTObject(hFiberScintPulseRMS);
  file->WriteTObject(hOldPulseRMS);
  file->WriteTObject(hOldTime);
  file->WriteTObject(hFiberTime);
  file->WriteTObject(hFiberScintTime);
  file->Close();
  exit(0);
}
  
void OneTrace()
{
  TRandom3 gen;
  gen.SetSeed(0);

  double startZ = 750;

  const float clocktick = 62.5;
  const float step  = 1.0; //ns
  const float RInt   = 420.;
  const float RSlope = 0.0;
  const float FInt   = 10373;
  const float FSlope = 2.4;
  const float IInt   = 112557;
  const float ISlope = 0.0;
  const int phase = 0;
  const int nclockticks = 8800;
  const int nperdigitization = 8;
  const double baseline = 400;
  
  PulseShaper shaper(IInt, ISlope, RInt, RSlope, FInt, FSlope);
  shaper.SetPhase(phase);
  shaper.SetClocktick(62.5);
  shaper.SetMaxStep(step);
  shaper.SetBaseline(baseline);
  shaper.SetNumClockticksInSpill(nclockticks);
  shaper.SetNumClockticksPerDigitization(nperdigitization);

  list< pair<float, float> > inputsNewMuon;
  list< pair<float, float> > inputsNewProton;
  list< pair<float, float> > inputsOld;
  
  const double length(1500); //length of cell
  
  TH1D* hPhotonTimeMuon   = new TH1D("hPhotonTimeMuon",   ";Photon Arrival Time (ns)", 5000, 0, 10000);
  TH1D* hPhotonTimeProton = new TH1D("hPhotonTimeProton", ";Photon Arrival Time (ns)", 5000, 0, 10000);
  TH1D* hOldPhotonTime    = new TH1D("hOldPhotonTime",    ";Photon Arrival Time (ns)", 5000, 0, 10000);

  TGraph* gNewTime = new TGraph();
  TGraph* gNewTimeA = new TGraph();
  TGraph* gOldTime = new TGraph();
  
  TFile* fCollection = new TFile("dT_dZ_CollectionRate.root");
  TH2D* hCollection = (TH2D*)fCollection->Get("dT_dZ_CollectionRate");
  
  TFile* fSpread = new TFile("FiberSpread.root");
  TH2D* hSpread = (TH2D*)fSpread->Get("FiberSpread");

  TFile* fOldSpread = new TFile("Dt_per_z_distribution.root");
  TH1D* hOldSpread = (TH1D*)fOldSpread->Get("Dt_per_z");
  
  double v0 = TMath::Ccgs()*1e-9/1.59;

  Long64_t niterations = 1e4;
  double weight = 1000./niterations;
  
  for (Long64_t iterations = 0; iterations < niterations; ++iterations)
    {
      if (iterations%10000 == 0) cout << "iteration: " << iterations << endl;
      double dZ(0.0), dT(0.0);
      hCollection->GetRandom2(dZ, dT);
      double Dsp = startZ + dZ;
      double Dlp = 2*length - Dsp;
      double Dsm = startZ - dZ;
      double Dlm = 2*length - Dsp;
      vector<double> lengths = {Dsp, Dlp, Dsm, Dlm};

      double q[]   = {0.933, 0.024, 0.022, 0.0213};
      double tau[] = {3.95, 23.56, 78.86, 546.39};
      
      double qA[]   = {0.679, 0.144, 0.102, 0.075};
      double tauA[] = {4.15, 19.90, 99.91, 617.96};
      double u = gen.Rndm();

      double prob(0);
      double scintTimeMuon(0.0);
      if (u < 0.05) scintTimeMuon = gen.Exp(3.6);
      else
	{
	  u = gen.Rndm();
	  for (int i = 0; i < 4; ++i)
	    {
	      prob += q[i];
	      if (u < prob)
		{
		  scintTimeMuon = gen.Exp(tau[i]);
		  break;
		}
	    }
	}

      double scintTimeProton = gen.Exp(3.6);
      
      /*
      double probProton(0);
      double scintTimeProton(0.0);
      for (int i = 0; i < 4; ++i)
	{
	  probProton += q[i];
	  if (u < probProton)
	    {
	      scintTimeProton = gen.Exp(tau[i]);
	      break;
	    }
	}
      

	double probA(0);
	double scintTimeA(0.0);
	for (int i = 0; i < 4; ++i)
	{
	probA += qA[i];
	if (u < probA)
	{
	scintTimeA = gen.Exp(tauA[i]);
	break;
	}
	}
      */
      
      for (unsigned int iL = 0; iL < lengths.size(); ++iL)
	{
	  double attenWeight = Atten(lengths[iL]);
	  double time = 2000 + dT;
	  time += lengths[iL]/v0;
	  int iBin = hSpread->GetXaxis()->FindBin(lengths[iL]);
	  TH1D* hDelay = hSpread->ProjectionY("hDelay", iBin, iBin);
	  hPhotonTimeMuon->Fill(time + hDelay->GetRandom() + scintTimeMuon, attenWeight);
	  
	  hPhotonTimeProton->Fill(time + hDelay->GetRandom() + scintTimeProton, attenWeight);
	  
	  hOldPhotonTime->Fill(time + hOldSpread->GetRandom()*lengths[iL] + gen.Exp(9.0), attenWeight);
	}
    }

  for (int iBin = 0; iBin < hPhotonTimeMuon->GetNbinsX(); ++iBin)
    {
      double time = hPhotonTimeMuon->GetBinCenter(iBin+1);
      double nPho = hPhotonTimeMuon->GetBinContent(iBin+1)*weight;
      inputsNewMuon.push_back( make_pair<float, float>( time, nPho ) );

      time = hPhotonTimeProton->GetBinCenter(iBin+1);
      nPho = hPhotonTimeProton->GetBinContent(iBin+1)*weight;
      inputsNewProton.push_back( make_pair<float, float>( time, nPho ) );
      
      time = hOldPhotonTime->GetBinCenter(iBin+1);
      nPho = hOldPhotonTime->GetBinContent(iBin+1)*weight;
      inputsOld.push_back( make_pair<float, float>( time, nPho ) );
    }

  vector<float> traceNewMuon = shaper.CreateTrace( inputsNewMuon );
  //for (unsigned int i = 0; i < nclockticks/nperdigitization; ++i) gNewTime->SetPoint(i, i, traceNew[phase + i*nperdigitization]);
  
  vector<float> traceNewProton = shaper.CreateTrace( inputsNewProton );
  //for (unsigned int i = 0; i < nclockticks/nperdigitization; ++i) gNewTimeA->SetPoint(i, i, traceNewA[phase + i*nperdigitization]);
  
  vector<float> traceOld = shaper.CreateTrace( inputsOld );
  //for (unsigned int i = 0; i < nclockticks/nperdigitization; ++i) gOldTime->SetPoint(i, i, traceOld[phase + i*nperdigitization]);
  
  //TGraph* gDCSNew = new TGraph();
  //TGraph* gDCSNewA = new TGraph();
  //TGraph* gDCSOld = new TGraph();
  
  cout << "Triggers for new, muon: " << endl;
  double adcMuon(0);
  double tnsMuon = DCSTrigger(traceNewMuon, phase, nclockticks, nperdigitization, adcMuon);

  cout << "Triggers for new, proton: " << endl;
  double adcProton(0);
  double tnsProton = DCSTrigger(traceNewProton, phase, nclockticks, nperdigitization, adcProton);

  cout << "Triggers for old: " << endl;
  double adcOld(0);
  double tnsOld = DCSTrigger(traceOld, phase, nclockticks, nperdigitization, adcOld);

  cout << setw(15) << "Type:"    << setw(15) <<  "TNS:"   << setw(15) << "ADC:"    << endl;
  cout << setw(15) << "Old:"     << setw(15) << tnsOld    << setw(15) << adcOld    << endl;
  cout << setw(15) << "Muon:"    << setw(15) << tnsMuon   << setw(15) << adcMuon   << endl;
  cout << setw(15) << "Proton:"  << setw(15) << tnsProton << setw(15) << adcProton << endl;
  
  //TCanvas* cTriggers = new TCanvas("cTriggers", "cTriggers");
  //gDCSNew->Draw("AL");
  //gDCSNewA->Draw("sameL");
  //gDCSOld->Draw("sameL");
  
  //TCanvas* cTraces = new TCanvas("cTraces", "cTraces");
  //gOldTime->SetLineColor(kBlack);
  //gOldTime->Draw("AL");
  //gNewTime->SetLineColor(kRed);
  //gNewTime->Draw("sameL");
  //gNewTimeA->SetLineColor(kBlue); 
  //gNewTimeA->Draw("sameL");
 
  TCanvas* cPhotonSignals = new TCanvas("cPhotonSignals", "cPhotonSignals");
  hOldPhotonTime->SetLineColor(kBlack);
  hOldPhotonTime->Draw("hist");
  hPhotonTimeMuon->SetLineColor(kRed);
  hPhotonTimeMuon->Draw("same hist");
  hPhotonTimeProton->SetLineColor(kBlue);
  hPhotonTimeProton->Draw("same hist");
}
