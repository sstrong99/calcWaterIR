#include "calculation.h"

//must include these here, not in header, because otherwise Traj isn't defined
//so children won't know what parent is
#include "calcIR.h"
#include "calcDists.h"
#include "calcIR_TAA.h"
#include "calcLocal.h"

Calculation *Calculation::createCalc(const Input &inp) {
  int whichCalc = inp.getWhichCalc();
  if (whichCalc==0)
    return new CalcIR(inp);
  else if (whichCalc==1)
    return new CalcDists(inp);
  else if (whichCalc==2)
    return new CalcIR_TAA(inp);
  else if (whichCalc==3)
    return new CalcLocal(inp);
  else {
    printf("ERROR: calc keyword accepts 0=IR,1=dists,2=TAA,3=localFreqs\n");
    exit(EXIT_FAILURE);
  }
}

int Calculation::getNsample(const int trajL,const int tsL,int &nSample) {
  //for a trajectory of length trajL and a timeseries of length tsL
  //find the number of samples given the requested nSample
  //input nSample=0: maximum number, -1: don't reuse trajectory data
  //nSample is changed to the correct value, and the sampleStep is returned
  int maxSample=trajL-tsL+1;
  int sampleStep;
  if ( nSample > maxSample ) {
    printf("WARNING: reducing the number of requested samples.\n");
    nSample=0;
  }

  if (nSample > 0) {
    sampleStep = floor((trajL-tsL)/nSample);
  } else if (nSample == 0) {
    sampleStep = 1;
    nSample = maxSample;
  } else if (nSample == -1) {
    sampleStep = tsL;
    nSample = floor(trajL/tsL);
  } else {
    printf("Undefined behavior for nSample = %d\n",nSample);
    exit(EXIT_FAILURE);
  }

  return sampleStep;
}

//get thread count even in serial part of program
int Calculation::ompNumThreads() {
  int n = 0;
#pragma omp parallel reduction(+:n)
  n += 1;
  return n;
}

void Calculation::print(const char *filename, const float *x, const float *y, const int sz)
{
  FILE *myfile = fopen(filename,"w");
  for (int ii=0; ii<sz; ii++)
    fprintf(myfile,"%.5e %.5e\n",x[ii],y[ii]);
  fclose(myfile);
}

void Calculation::print(const char *filename, const float *x, const float *y1, const float *y2,const int sz)
{
  FILE *myfile = fopen(filename,"w");
  for (int ii=0; ii<sz; ii++)
    fprintf(myfile,"%.5e %.5e %.5e\n",x[ii],y1[ii],y2[ii]);
  fclose(myfile);
}
