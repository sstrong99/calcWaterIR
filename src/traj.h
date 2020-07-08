#ifndef TRAJ_H
#define TRAJ_H

#include "vecManip.h"

#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#define A0INV 18.89726125    //inverse bohr radius (nm/a0)
#define PI 3.14159265359

using namespace std;
class Traj {
public:
  static Traj *getTraj(const string &filename);
  virtual ~Traj() {};
  
  virtual int next(const bool convertFlag=true) = 0;
  virtual float allT() = 0;
  virtual void skip(const int n) = 0;

  int getNatoms() const {return natoms;};
  int getNT() const {return nT;};
  float getT() const {return t;};

  const rvec* getCoords() const { return x; };
  void getBox(rvec &box) const;
  void moveM(const float &frac,const int aPerM);

  int getModel() const;

protected:
  int natoms;     //number of atoms
  int nT;         //number of timesteps processed

  int step;       //integer timestep
  float t;        //time in ps (units correct?)
  rvec box;       //simulation box
  rvec *x;        //atom positions

  float dt;       //timestep between current step and last
};

#endif
