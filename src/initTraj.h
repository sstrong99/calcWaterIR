#ifndef INITTRAJ_H
#define INITTRAJ_H

#include "traj.h"
#include "calcW.h"

#include <string>

class InitTraj {
public:
  InitTraj(const string &trajfile);
  ~InitTraj() {};

  void printModel();
  string modelString();
  float getAvgF() {return avgF;};
  int getNT() {return nT;};
  int getModel() {return model;};
  int getNH() {return nH;};
  float getDT() {return dt;};
  int adjustTimestep(const float &new_dt);

private:
  int nT; //number of timesteps
  float dt; //timestep
  int nH;
  int model;
  float avgF;
};
#endif
