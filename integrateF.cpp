#include "integrateF.h"

//initialize Fw to identity
void IntegrateF::initF(cpx *F,const int nH)
{
  int ii,jj;
  for (ii=0; ii<nH; ii++)
  {
    for (jj=0; jj<nH; jj++)
      F[jj+ii*nH]=0.0;

    F[ii+ii*nH]=1.0;
  }
}
