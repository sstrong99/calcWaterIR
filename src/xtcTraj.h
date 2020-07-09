#ifndef XTCTRAJ_H
#define XTCTRAJ_H

#include "vecManip.h"
#include "traj.h"

#include <xdrfile_xtc.h>
#include <xdrfile.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>

using namespace std;
class XTCTraj : public Traj {
public:
  XTCTraj(const char *xtcfile);
  ~XTCTraj();
  int next(const bool convertFlag=true);
  float allT();
  void skip(const int n);

private:
  XDRFILE *trj;   //pointer to file
  float prec;     //precision
  matrix boxT;    //simulation box in raw output

  void convert();
};

#endif
