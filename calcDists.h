#ifndef CALCDISTS_H
#define CALCDISTS_H

#include "calculation.h"
#include "traj.h"
#include "calcW.h"
#include "histogram.h"
#include "exactDiag.h"
#include "timer.h"
#include "input.h"
#include "vecManip.h"
#include "initTraj.h"

#include <omp.h>
#include <string>
#include <cmath>

class CalcDists : public Calculation {
public:
  CalcDists(const Input &inp);
  ~CalcDists() {};

  void printResults(string postfix) const;

private:
  void loopSamples();
  Histogram Pu_tot,Pc_tot,spdn_tot,Pintra_tot;
};
#endif
