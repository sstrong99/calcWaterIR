#ifndef CALCULATION_H
#define CALCULATION_H

#include "input.h"
#include <cmath>
#include <cstdio>

class Calculation {
public:
  Calculation() {};
  virtual ~Calculation() {};
  virtual void printResults(string postfix) const = 0;

protected:
  static int getNsample(const int trajL,const int tsL,int &nSample);
  static int ompNumThreads();
  static void print(const char *filename, const float *x, const float *y, const int sz);
  static void print(const char *filename, const float *x, const float *y1, const float *y2, const int sz);
};
#endif
