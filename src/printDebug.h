#ifndef PRINTDEBUG_H
#define PRINTDEBUG_H

#include "mycomplex.h"
#include <xdrfile_xtc.h>
#include <cstdio>

class PrintDebug {
public:
  PrintDebug();
  ~PrintDebug() {};

  void setNH(const int nH1);

  void mat(const float *mat) const;
  void mat(const cpx *mat) const;
  void mat(const char *filename, const cpx *mat) const;
  void vec(const float *vec) const;
  void vec(const rvec *vec) const;
private:
  int nH;
  int max;
};

#endif
