#ifndef CALCW_H
#define CALCW_H

#include "h2xind.h"
#include "map.h"
#include "traj.h"
#include "vecManip.h"
#include "printDebug.h"

#include <cmath>

#define CM2PS 0.18836516 //cm/ps = 1/2pi*c
#define HART2CM 2.1947463e5 //convert hartree to wavenumber

class CalcW {
public:
  CalcW(const int model,const int natoms,float avef=0.0);
  ~CalcW();
  void compute(Traj &traj,rvec *m);

  int getNH() const {return nH;};

  float calcAveF();
  void getW(float* out) const ;
  const float* getW() const { return wMat; };
  inline float getW(int nn,int mm) const
  { return wMat[mm+nn*nH]; };

  float maxInterFreq() const;
  float avgIntraFreq() const;

private:
  void calcE(const Traj &traj);
  void calcE_wrong(const Traj &traj);
  void mapE2W(rvec *m);

  inline void setNN(float *M, const float val, const int nn, const int mm)
  { M[mm+nn*nH]=val; };
  inline float getDipdip(const int nn,const int mm)
  { return dipdip[mm+nn*nH]; };

  float *charges;
  int aPerM; //atoms per molecule
  int nO;
  int nH;
  PrintDebug prnt; //need to init

  float *E;      //scalar Efield at each H, along OH
  rvec *OH;       //OH vectors
  float *dipdip; //dipole dipole term
  Map *mymap;
  H2Xind *h2x;
  float *wMat;
  float avef;

  float moveMfrac;   //how much to move M by
  int mInc;       //increment O ind by this much to get M ind
};


#endif
