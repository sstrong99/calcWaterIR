#ifndef CALCIR_TAA_H
#define CALCIR_TAA_H

#include "calculation.h"
#include "traj.h"
#include "calcW.h"
#include "exactDiag.h"
#include "timer.h"
#include "initTraj.h"
#include "input.h"

#include <omp.h>
#include <cstdio>
#include <vector>
#include <string>

class CalcIR_TAA : public Calculation {
public:
  CalcIR_TAA(const Input &inp);
  ~CalcIR_TAA();

  void printResults(string postfix) const;

private:
  const float T1;
  int nTavg;
  int N;  //number of frequency points

  const string xtcfile;
  int model;

  float *spectrum;
  float *w;

  float timestep;   //in ps
  int dt_skip;
  int nH,nT;

  void init(const float &Tavg);
  void loopSamples(const int nSample,const int step);
  void calcAvg(const int start, float *spec1);
  inline float lorentz(const float &w);

  //for debugging
  PrintDebug prnt;
};

#endif
