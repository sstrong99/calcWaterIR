#ifndef CALCIR_H
#define CALCIR_H

#include "calculation.h"
#include "traj.h"
#include "calcW.h"
#include "exactDiag.h"
#include "adamsBashforth.h"
#include "vecManip.h"
#include "timer.h"
#include "mycomplex.h"
#include "initTraj.h"
#include "input.h"

#include <omp.h>
#include <fftw3.h>
#include <cstdio>
#include <vector>
#include <string>

class CalcIR : public Calculation {
public:
  CalcIR(const Input &inp);
  ~CalcIR();

  void printResults(string postfix) const;

private:
  const float T1;

  const char *xtcfile;
  int integrator,model;

  float *spectrum;
  float *w;
  cpx    *avgCorr;

  float timestep;   //in ps
  int dt_skip;
  int nH,nT;
  int nTCF;        //number of timesteps in TCF
  int N;           //size of FFT
  float avgF;

  void init();
  void loopSamples(const int nSample,const int step);
  void calcTCF(const int start, cpx *corr1);
  void calcFFT(cpx *y);
  static int nextPow2(int v);

  cpx sumMFM(const rvec *m0,const rvec *m,const cpx *F);

  //for debugging
  float norm(const cpx *mat);
  PrintDebug prnt;
  //TODO: make PrintDebug static, only one instance in program
};

#endif
