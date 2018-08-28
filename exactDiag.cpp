#include "exactDiag.h"
//use lapack to diagonalize symmetric float matricies
#ifndef USEGPU
extern "C" {
  extern int ssyevd_(char*,char*,int*,float*,int*,float*,float*,int*,int*,int*,int*);
  extern int cgemm_(char*,char*,int*,int*,int*,cpx*,cpx*,int*,cpx*,int*,cpx*,cpx*,int*);
  extern int sgemm_(char*,char*,int*,int*,int*,float*,float*,int*,float*,int*,float*,float*,int*);
  extern int ccopy_(int *,cpx*,int*,cpx*,int*);
}
#endif

ExactDiag::ExactDiag(int nH,const float &dt,const bool vecflag) :
  nH(nH),nH2(nH*nH),dt(dt)
{

//get workspace size for diagonalization
  float tmpwork;
  myint tmpiwork;
  myint info;

  //upper triangular part in row-major format is
  //equivalent to lower triangular in column-major
  bool upperflag=false;

#ifdef USEGPU
  if (vecflag)   vec = MagmaVec;
  else           vec = MagmaNoVec;
  if (upperflag) uplo=MagmaUpper;
  else           uplo=MagmaLower;

  magma_init();
  magma_queue_create(0,&queue);

  magma_ssyevd(vec,uplo,nH,NULL,nH,NULL,
	       &tmpwork,-1,&tmpiwork,-1,&info);
//  magma_ssyevd_gpu(vec,uplo,nH,NULL,nH,NULL,NULL,nH,
//  	       &tmpwork,-1,&tmpiwork,-1,&info);

  //4 matricies stored on
  int gpuMemMM = 4*sizeof(cpx)*nH2 * omp_get_num_threads();
  // retrieve and print info about gpu
  //cudaDeviceProp prop;
  //int totalGPUmem = prop.totalGlobalMem;
  size_t free,total;
  cudaMemGetInfo(&free,&total);
  float warnFrac=0.95;
  //printf("%d %d %d\n",gpuMemMM,(int) free,(int) total);
  if (gpuMemMM > warnFrac*free) {
    printf("Warning: %.1f%% of the GPU's free memory will be used for matrix multiplication\n",100.0 * (float) gpuMemMM / (float) free);
  }
#else
  if (vecflag)   vec='V';
  else           vec='N';
  if (upperflag) uplo='U';
  else           uplo='L';

  lwork=-1; liwork=-1;
  ssyevd_(&vec,&uplo,&nH,NULL,&nH,NULL,
	  &tmpwork,&lwork,&tmpiwork,&liwork,&info);
#endif
  //printf("lwork = %d\t liwork = %d\n",(myint) tmpwork,tmpiwork);

  lwork  = (myint) tmpwork;
  liwork = tmpiwork;
  work = new float[lwork];
  iwork= new myint[liwork];
  evs = new float[nH2];
  cevs= new cpx[nH2];
  w   = new float[nH];
  prop= new cpx[nH2];
  tmp = new cpx[nH2];
}

ExactDiag::~ExactDiag()
{
#ifdef USEGPU
  magma_queue_destroy(queue);
#endif
  delete[] work;
  delete[] iwork;
  delete[] evs;
  delete[] cevs;
  delete[] w;
  delete[] prop;
  delete[] tmp;
}

void ExactDiag::next(cpx *F,const float *wMat)
{
  //diagonalize w matrix
  diag(wMat);

  //test O^T W O = D
  //testmult(wMat);
//#ifdef USEGPU
    //testmult_gpu(wMat);
//#endif

  //propagate F
  propF(F);
}

void ExactDiag::spdn(const float *wMat,const rvec *m,float *out_weights)
{
  //diagonalize w matrix
  diag(wMat);

  //zero weights
  int jj;
  for (jj=0; jj<nH; jj++)
    out_weights[jj]=0.0;

  //set up variables on GPUs for multiplication
#ifdef USEGPU
  magma_smalloc(&tmpm_d,nH);
  magma_smalloc(&tmpmOut_d,nH);
  magma_smalloc(&evs_d,nH2);
  magma_ssetmatrix(nH,nH,evs,nH,evs_d,nH,queue);
  float *tmpm = new float[nH];
  float *tmpm_eb = new float[nH];
  int ii;

  //project m to eigenbasis
  for (ii=0; ii<DIM; ii++) {
    for (jj=0; jj<nH; jj++)
      tmpm[jj]=(float) m[jj][ii];
    magma_ssetvector(nH,tmpm,1,tmpm_d,1,queue);
    magma_sgemv( MagmaTrans,nH,nH,1.0,evs_d,nH,tmpm_d,1,0.0,tmpmOut_d,1,queue);
    magma_sgetvector(nH,tmpmOut_d,1,tmpm_eb,1,queue);
    for (jj=0; jj<nH; jj++)
      out_weights[jj]+=tmpm_eb[jj]*tmpm_eb[jj];
  }

  //cleanup
  magma_free(tmpm_d);
  magma_free(tmpmOut_d);
  magma_free(evs_d);

  delete[] tmpm;
  delete[] tmpm_eb;
#else
  //This function doesn't do anything if compiled without GPUs yet
  printf("ERROR: Code has not been written to compute the spectral density without GPUs\n");
  exit(EXIT_FAILURE);
#endif
}

void ExactDiag::mult_taa(const rvec *m,float *d) {
  int ii,pp;
  float sum;
  float tmpE;
  for (int kk=0; kk<nH; kk++) {
    sum=0.0;
    for (ii=0; ii<nH; ii++) {
      tmpE=evs[ii+kk*nH];
      for (pp=0; pp<DIM; pp++)
	sum+=tmpE*m[ii][pp];
    }
    d[kk]=sum/DIM;
  }
}

void ExactDiag::diag(const float *wMat)
{
  //only copy upper triangular part of W
  int jj;
  for (int ii=0; ii<nH; ii++)
    for (jj=ii; jj<nH; jj++)
      evs[jj+ii*nH] = wMat[jj+ii*nH];
  //stored upper triangular part in row-major format
  //equivalent to lower triangular in column-major
  //magma/lapack use column-major (fortran-style)

  //printMat(evs);

  //diagonalize wMat
  //maybe should leave evs on GPU memory, for multiplication next
  int info;
#ifdef USEGPU
  magma_ssyevd(vec,uplo,nH,evs,nH,w,work,lwork,iwork,liwork,&info);
  //float *wA = new float[nH*nH];
  //float *evs_d;
  //magma_smalloc(&evs_d,nH2);
  //magma_ssetmatrix(nH,nH,evs,nH,evs_d,nH,queue);
  //magma_ssyevd_gpu(vec,uplo,nH,evs_d,nH,w,wA,nH,work,lwork,iwork,liwork,&info);
  //magma_sgetmatrix(nH,nH,evs_d,nH,evs,nH,queue);
  //delete[] wA;
  //magma_free(evs_d);
#else
  ssyevd_(&vec,&uplo,&nH,evs,&nH,w,work,&lwork,iwork,&liwork,&info);
#endif
  //for some reason, the top row of the eigenvector matrix has the opposite
  //sign as in nick's code, but this doesn't change the results
}

void ExactDiag::propF(cpx *F)
{
  //multiply diagonal prop matrix by ev^T
  cpx diagEl;
  int jj;
  for (int ii=0; ii<nH; ii++)
  {
    diagEl=exp(pureIm*w[ii]*dt);  //dt in ps, w in 1/ps
    //diagEl=1.0;  //SES tested 6/4/18 that this gives O*O^T = I
    for (jj=0; jj<nH; jj++)
    {
      tmp[jj+ii*nH]  = evs[jj+ii*nH] * diagEl; //O*e^D
      cevs[jj+ii*nH] = evs[jj+ii*nH]; //automatically assigns to real part
    }
  }

#ifdef USEGPU
  //copy to GPU memory
  allocateGPUmem(); //when allocated in constructor, dsyev nullified these ptrs
  magma_csetmatrix(nH,nH,(mcpx*) F   ,nH,F_d   ,nH,queue);
  magma_csetmatrix(nH,nH,(mcpx*) cevs,nH,cevs_d,nH,queue);
  magma_csetmatrix(nH,nH,(mcpx*) tmp ,nH,tmp_d ,nH,queue);

  //perform multiplications
  //magmablas_zlascl_diag only works with triangluar cevs, not full
  //magmablas_zlascl2 only works with real diagonal matrix, not cpx
  magma_cgemm(MagmaNoTrans,MagmaTrans,nH,nH,nH,
	      MAGMA_C_ONE,cevs_d,nH,tmp_d,nH,
	      MAGMA_C_ZERO,prop_d,nH,queue); //(O*e^D)^T = e^D O^T
  magma_cgemm(MagmaNoTrans,MagmaNoTrans,nH,nH,nH,
	      MAGMA_C_ONE,prop_d,nH,F_d,nH,
	      MAGMA_C_ZERO,tmp_d,nH,queue);

  //copy F back to CPU memory
  //magma_cgetmatrix(nH,nH,tmp_d,nH,(mcpx*) F,nH,queue);
  magma_cgetvector(nH2,tmp_d,1,(mcpx*) F,1,queue);

  magma_queue_sync(queue);  //is this necessary?
  freeGPUmem();
#else
  char notrans='N'; char trans='T';
  cpx one  = 1.0;
  cpx zero = 0.0;
  int ione=1;

  cgemm_(&notrans,&trans,&nH,&nH,&nH,&one,cevs,&nH,tmp,&nH,&zero,prop,&nH);
  cgemm_(&notrans,&notrans,&nH,&nH,&nH,&one,F,&nH,prop,&nH,&zero,tmp,&nH);
  ccopy_(&nH2,tmp,&ione,F,&ione);
#endif
}

#ifdef USEGPU
void ExactDiag::allocateGPUmem()
{
  magma_cmalloc(&cevs_d,nH2);
  magma_cmalloc(&prop_d,nH2);
  magma_cmalloc(&tmp_d,nH2);
  magma_cmalloc(&F_d,nH2);
}

void ExactDiag::freeGPUmem()
{
  magma_free(cevs_d);
  magma_free(prop_d);
  magma_free(tmp_d);
  magma_free(F_d);
}
void ExactDiag::testmult_gpu(const float* wMat)
{
  float *tmp1 = new float[nH2];
  float *evs1_d,*mat1_d,*tmp1_d;

  for (int ii=0; ii<nH2; ii++)
    tmp1[ii]  = wMat[ii];

  magma_smalloc(&evs1_d,nH2);
  magma_smalloc(&mat1_d,nH2);
  magma_smalloc(&tmp1_d,nH2);

  magma_ssetmatrix(nH,nH,evs  ,nH, evs1_d,nH,queue);
  magma_ssetmatrix(nH,nH,tmp1 ,nH, mat1_d,nH,queue);

  magma_sgemm(MagmaNoTrans,MagmaNoTrans,nH,nH,nH,
	      1.0,mat1_d,nH,evs1_d,nH,
	      0.0,tmp1_d,nH,queue);
  magma_sgemm(MagmaTrans,MagmaNoTrans,nH,nH,nH,
	      1.0,evs1_d,nH,tmp1_d,nH,
	      0.0,mat1_d,nH,queue);

  magma_sgetmatrix(nH,nH,mat1_d,nH,tmp1,nH,queue);

  magma_queue_sync(queue);

  magma_free(evs1_d);
  magma_free(mat1_d);
  magma_free(tmp1_d);
  delete[] tmp1;
}
#else
void ExactDiag::testmult(const float* wMat)
{
  char notrans='N'; char trans='T';
  float one  = 1.0;
  float zero = 0.0;

  float *tmp1 = new float[nH2];
  float *tmp2 = new float[nH2];
  int jj;

  float tmpel;
  int kk;
  for (int ii=0; ii<nH; ii++)
    for (jj=0; jj<nH; jj++) {
      tmpel=0.0;
      for (kk=0; kk<nH; kk++)
	tmpel += evs[kk+ii*nH]*wMat[jj+kk*nH];
      tmp1[jj+ii*nH] = tmpel;
    }

  for (int ii=0; ii<nH; ii++)
    for (jj=0; jj<nH; jj++) {
      tmpel=0.0;
      for (kk=0; kk<nH; kk++)
	tmpel += tmp1[kk+ii*nH]*evs[kk+jj*nH];
      tmp2[jj+ii*nH] = tmpel;
    }

  //dgemm_(&trans,&notrans,&nH,&nH,&nH,&one,evs,&nH,wMat,&nH,&zero,tmp1,&nH);
  //dgemm_(&notrans,&notrans,&nH,&nH,&nH,&one,tmp1,&nH,evs,&nH,&zero,tmp2,&nH);

  sgemm_(&trans,&notrans,&nH,&nH,&nH,&one,evs,&nH,evs,&nH,&zero,tmp1,&nH);

  delete[] tmp1;
  delete[] tmp2;
}
#endif
