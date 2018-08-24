#include "adamsBashforth.h"

AdamsBashforth::AdamsBashforth(const int nprev,const int nH,const float &dt) : nprev(nprev),nH(nH),iFw(nprev),dt(dt)
{
//declare iFw matrix
  for (int ii=0; ii<nprev; ii++)
    iFw[ii]=new cpx[nH*nH];

  nextInd=0;
  trunc=true;

  if (nprev>4)
  {
    printf("Need to add code to do Adams-Bashforth greater than 4th order\n");
    exit(EXIT_FAILURE);
  }
}

AdamsBashforth::~AdamsBashforth()
{
  for (int ii=0; ii<nprev; ii++)
    delete[] iFw[ii];
}

void AdamsBashforth::next(cpx *F,const float *wMat)
{
  //calculate next iFw, and store it
  //multipy Flast with wMat
  cpx tmpout;
  int ii,jj,kk;
  for (ii=0; ii<nH; ii++)
    for (jj=0; jj<nH; jj++)
    {
      tmpout=0.0;
      for (kk=0; kk<nH; kk++) {
	tmpout += wMat[kk+ii*nH] * F[jj+kk*nH];
	//tmpout += F[kk+ii*nH] * wMat[jj+kk*nH];
      }
      iFw[nextInd][jj+ii*nH]=tmpout*pureIm*dt;
    }
  //check if have enough results to do full order Adams-Bashforth
  if (trunc && nextInd==nprev-1)
    trunc=false;

  //compute truncated Adams-Bashforth
  if (trunc) {
    switch(nextInd)
    {
    case 0:
      addMatAB(F,iFw[nextInd]);
      break;
    case 1:
      addMatAB(F,iFw[nextInd],iFw[nextInd-1]);
      break;
    case 2:
      addMatAB(F,iFw[nextInd],iFw[nextInd-1],iFw[nextInd-2]);
      break;
    case 3:
      addMatAB(F,iFw[nextInd],iFw[nextInd-1],iFw[nextInd-2],iFw[nextInd-3]);
      break;
    }
  } else { //compute untruncated Adams-Bashforth
    //add nprev to indicies to prevent negative indicies
    addMatAB(F,iFw[nextInd],iFw[(nextInd-1+nprev) % nprev ],
	     iFw[(nextInd-2+nprev) % nprev],iFw[(nextInd-3+nprev) % nprev]);
  }

  //iterate nextInd, better than modulo
  nextInd = (nextInd+1 == nprev ? 0 : nextInd+1);
}

inline void AdamsBashforth::assignMat(cpx *Fnew,const cpx *Fold)
{
  int ii,jj;
  for (ii=0; ii<nH; ii++)
    for (jj=0; jj<nH; jj++)
      Fnew[jj+ii*nH]=Fold[jj+ii*nH];
};

inline void AdamsBashforth::addMatAB(cpx *F, const cpx *Fadd)
{
  int ii,jj;
  for (ii=0; ii<nH; ii++)
    for (jj=0; jj<nH; jj++)
      F[jj+ii*nH]+=Fadd[jj+ii*nH];
};

inline void AdamsBashforth::addMatAB(cpx *F, const cpx *F1, const cpx *F2)
{
  int ii,jj;
  for (ii=0; ii<nH; ii++)
    for (jj=0; jj<nH; jj++)
      F[jj+ii*nH]+=  1.5f  * F1[jj+ii*nH]
		   - 0.5f  * F2[jj+ii*nH];
};

inline void AdamsBashforth::addMatAB(cpx *F, const cpx *F1, const cpx *F2, const cpx *F3)
{
  int ii,jj;
  for (ii=0; ii<nH; ii++)
    for (jj=0; jj<nH; jj++)
      F[jj+ii*nH]+=  1.9166666666666667f  * F1[jj+ii*nH]
		   - 1.3333333333333333f  * F2[jj+ii*nH]
		   + 0.41666666666666667f * F3[jj+ii*nH];
};

inline void AdamsBashforth::addMatAB(cpx *F, const cpx *F1, const cpx *F2, const cpx *F3,const cpx *F4)
{
  int ii,jj;
  for (ii=0; ii<nH; ii++)
    for (jj=0; jj<nH; jj++)
      F[jj+ii*nH]+=  2.2916666666666667f * F1[jj+ii*nH]
		   - 2.4583333333333333f * F2[jj+ii*nH]
		   + 1.5416666666666667f * F3[jj+ii*nH]
		   - 0.375f              * F4[jj+ii*nH];
};
