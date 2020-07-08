#ifndef ADAMSBASHFORTH_H
#define ADAMSBASHFORTH_H

#include "integrateF.h"
#include "mycomplex.h"

#include <vector>

//integrate F matrix using Adams-Bashforth method
using namespace std;
class AdamsBashforth : public IntegrateF
{
public:
  AdamsBashforth(const int nprev,const int nH,const float &dt);
  ~AdamsBashforth();

  void next(cpx *F,const float *wMat);

private:
  const int nprev;       //number of previous timesteps of Fw matrix to store
  int nextInd;     //index of next Fw to overwrite
  const int nH;
  bool trunc;      //true if not enough results to do full Adams-Bashforth yet
  vector<cpx*> iFw;
  const float dt;

  void assignMat(cpx *Fnew,const cpx *Fold);
  void addMatAB(cpx *F,const cpx *Fadd);
  void addMatAB(cpx *F,const cpx *F1,const cpx *F2);
  void addMatAB(cpx *F,const cpx *F1,const cpx *F2,const cpx *F3);
  void addMatAB(cpx *F,const cpx *F1,const cpx *F2,const cpx *F3,const cpx *F4);
};
#endif
