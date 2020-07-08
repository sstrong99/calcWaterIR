#ifndef EXACTDIAG_H
#define EXACTDIAG_H

#include "integrateF.h"
#include "mycomplex.h"

#ifdef USEGPU
#include "magma_v2.h"
#include <cuda_runtime_api.h>  //to get free GPU memory
#include <omp.h> //to compute GPU memory use
typedef magma_int_t myint;
typedef magmaFloatComplex mcpx;
#else
typedef int myint;
#endif

class ExactDiag : public IntegrateF
{
public:
  ExactDiag(int nH,const float &dt,const bool vecflag=true);
  ~ExactDiag();

  void next(cpx *F,const float *wMat);
  void spdn(const float *wMat,const rvec *m,float *out_weights);
  void diag(const float *wMat);
  void mult_taa(const rvec *m,float *d);
  float getEigenvalue(const int ii) const {return w[ii];};

private:
  void propF(cpx *F);
  void testmult(const float* wMat);

  myint nH,nH2; //can't make const b/c of fortran interfaces
  const float dt;
  myint lwork;
  myint liwork;

  float *work;
  myint *iwork;
  float *evs;
  cpx *cevs;
  float *w;
  cpx *prop;
  cpx *tmp;

#ifdef USEGPU
  magma_queue_t queue;

  //variables on GPU memory
  mcpx *cevs_d;
  mcpx *prop_d;
  mcpx *tmp_d;
  mcpx *F_d;
  float *tmpm_d;
  float *tmpmOut_d;
  float *evs_d;

  void printMat(const mcpx *mat_d) const;
  void allocateGPUmem();
  void freeGPUmem();
  void testmult_gpu(const float* wMat);

  magma_vec_t vec;  //calc eigen vecs or not?
  magma_uplo_t uplo; //upper or lower triangular dsyev
#else
  char vec;
  char uplo;
  bool warnflag;
#endif
};
#endif
