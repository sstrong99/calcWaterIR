#ifndef INTEGRATEF_H
#define INTEGRATEF_H

#include "traj.h"
#include "mycomplex.h"

using namespace std;
class IntegrateF
{
public:
  IntegrateF();
  virtual ~IntegrateF() {};

  virtual void next(cpx *F,const float* wMat) = 0;
  static void initF(cpx *F,const int nH);

  const cpx pureIm;  //pure imaginary number i
};

#endif
