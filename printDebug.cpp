#include "printDebug.h"

PrintDebug::PrintDebug() : nH(-1),max(10) {};

void PrintDebug::setNH(const int nH1) {
    nH=nH1;
    max = max>nH ? nH : max;
  };

void PrintDebug::mat(const float *mat) const {
  float tmpval;
  int jj;
  for (int ii=0; ii<max; ii++)
  {
    for (jj=0; jj<max; jj++)
    {
      tmpval=mat[jj+ii*nH];
      tmpval = fabs(tmpval)<1e-13 ? 0.0 : tmpval;
      printf("%9.2e ",tmpval);
    }
    printf("\n");
  }
  printf("\n");
}

void PrintDebug::mat(const cpx *mat) const {
  float tmpval;
  for (int ii=0; ii<max; ii++)
  {
    for (int jj=0; jj<max; jj++)
    {
      tmpval=real(mat[jj+ii*nH]);
      tmpval = fabs(tmpval)<1e-13 ? 0.0 : tmpval;
      printf("%9.2e ",tmpval);
    }

    printf("\n");
  }
  printf("\n");
}

void PrintDebug::mat(const char *filename, const cpx *mat) const {
  FILE *myfile = fopen(filename,"w");
  for (int ii=0; ii<max; ii++)
  {
    for (int jj=0; jj<max; jj++)
      fprintf(myfile,"%5.1f ",real(mat[jj+ii*nH]));

    fprintf(myfile,"\n");
  }
  fclose(myfile);
}

void PrintDebug::vec(const float *vec) const {
  for (int ii=0; ii<max; ii++)
    printf("%8.4f\n",vec[ii]);

  printf("\n\n");
}

void PrintDebug::vec(const rvec *vec) const {
  for (int ii=0; ii<max; ii++) {
    for (int jj=0; jj<DIM; jj++)
      printf("%12.4e ",vec[ii][jj]);

    printf("\n");
  }
  printf("\n\n");
}
