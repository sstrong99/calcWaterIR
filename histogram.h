#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <utility>  //for std::swap
#include <cstdlib>
#include <cstdio>
#include <cmath>

class Histogram {
public:
  Histogram(const int N,const float &lo,const float &hi);
  Histogram(const Histogram &hist);
  ~Histogram();

  void addData(const float &pt,const float &weight=1.0);
  void print(const char *filename) const;

  Histogram& operator+=(const Histogram &rhs);
  friend Histogram operator+(const Histogram &a,const Histogram &b);
  bool operator==(const Histogram &rhs) const;
  bool operator!=(const Histogram &rhs) const;

private:
  const int N;  //number of bins
  const float lo,hi;  //range of histogram
  const float dx;     //width of bin
  float *binCenters;
  float *binEdges;
  float *counts;  //not int, for weighted histograms

  void init();
};
#endif
