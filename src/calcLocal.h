#ifndef CALCLOCAL_H
#define CALCLOCAL_H

#include "calculation.h"
#include "traj.h"
#include "calcW.h"
#include "timer.h"
#include "input.h"
#include "initTraj.h"
#include "vecManip.h"

#include <omp.h>
#include <string>
#include <vector>
#include <cmath>

class CalcLocal : public Calculation {
public:
  CalcLocal(const Input &inp);
  ~CalcLocal() {};

  void printResults(string postfix) const;

private:
  vector<vector<float>> localFreqs;
  vector<vector<float>> transStrength;
};
#endif
