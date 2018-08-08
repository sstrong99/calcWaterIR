#ifndef INTEGRATEF_H
#define INTEGRATEF_H

#include "traj.h"
#include "printDebug.h"
#include "mycomplex.h"

using namespace std;
class IntegrateF
{
public:
  IntegrateF() {};
  virtual ~IntegrateF() {};

  virtual void next(cpx *F,const float* wMat) = 0;
  static void initF(cpx *F,const int nH);
  PrintDebug prnt;
};

#endif
