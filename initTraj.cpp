#include "initTraj.h"

InitTraj::InitTraj(const char *xtcfile) {
  printf("Reading trajectory file to find number of timesteps...\n");
  Traj traj(xtcfile);
  dt=traj.allT();
  nT=traj.getNT();
  model = traj.getModel();

  CalcW calcW(model,traj.getNatoms());
  nH=calcW.getNH();

  rvec *junk = new rvec[nH];
  calcW.compute(traj,junk);
  avgF = calcW.calcAveF();
  delete[] junk;
}

string InitTraj::modelString() {
  string str;
  switch (model) {
  case 0 : str="SPC/E"; break;
  case 1 : str="TIP4P"; break;
  case 2 : str="TIP4P/2005"; break;
  default :
    printf("ERROR: Water model was not identified\n");
    exit(EXIT_FAILURE);
  }

  return str;
}

int InitTraj::adjustTimestep(const float &new_dt) {
  //if pass 0, use the same timestep as trajectory
  if (new_dt == 0.0)
    return 1;

  float tRatio = new_dt/dt;
  int dt_skip;
  if ( fabs(tRatio-1) < 1e-3 ) {
    return 1;
  } else if ( tRatio < 1 ) {
    printf("ERROR: timestep is less than trajectory timestep.\n");
    exit(EXIT_FAILURE);
  } else if (fabs(round(tRatio)-tRatio) < 1e-3 ) {
    dt_skip=round(tRatio);
  } else {
    printf("ERROR: timestep is not a multiple of the trajectory timestep\n");
    exit(EXIT_FAILURE);
  }
  dt*=dt_skip;

  nT=floor((float) nT/(float) dt_skip);
  return dt_skip;
}

void InitTraj::printModel() {
  printf("%d %s water molecules\n",nH/2,modelString().c_str());
}
