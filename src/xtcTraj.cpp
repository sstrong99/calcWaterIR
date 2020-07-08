#include "xtcTraj.h"

using namespace std;
XTCTraj::XTCTraj(const char *xtcfile)
{
  //get number of atoms
  read_xtc_natoms((char*)xtcfile,&natoms);

  //open xtc file
  trj = xdrfile_open(xtcfile,"r");
  if (trj == NULL) {
    printf("Can't open %s\n",xtcfile);
    exit(EXIT_FAILURE);
  }

  //initialize vars
  nT = 0;
  x  = new rvec[natoms];
}

XTCTraj::~XTCTraj()
{
  delete[] x;
  xdrfile_close(trj);
}

//read next step of trajectory
int XTCTraj::next(const bool convertFlag)
{
  //read one timestep, returns 0 if success
  int stat=read_xtc(trj,natoms,&step,&t,boxT,x,&prec);

  //return if finished reading file
  if (stat) return stat;

  if (convertFlag) convert();

  nT++;

  return stat;
}

float XTCTraj::allT()
{
  if (nT != 0 ) { printf("WARNING: this Traj object has already been used\n");}

  //read one step to get 1st time
  int stat=next(false);
  float tlast=t;

  //loop rest of timesteps
  float dtsum=0.0;
  while (!stat)
  {
    stat=next(false);
    dtsum += t-tlast;
    tlast=t;
  }

  //convert last step for future processing
  convert();

  nT--; //for last timestep that wasn't read

  //round dt to nearest tenth of a fs
  float dtavg=dtsum*10000/((float) nT-1);
  if (dtavg < 1.0) {
    printf("ERROR: the code only works with dt>0.1 fs.\n");
    exit(EXIT_FAILURE);
  }

  return round(dtavg)/10000;
}

void XTCTraj::skip(const int n)
{
  if (n<0) {
    printf("ERROR: n must be non-negative\n");
    exit(EXIT_FAILURE);
  }

  //loop will never execute if n is zero
  int stat;
  for (int ii=0; ii<n; ii++) {
    stat=next(false);

    if (stat) {
      printf("Tried to go past last timestep: %d\n",nT);
      exit(EXIT_FAILURE);
    }
  }
}

void XTCTraj::convert()
{
  for (int jj=0; jj<natoms; jj++)
    for (int kk=0; kk<DIM; kk++)
      x[jj][kk]*=A0INV;

  for (int ii=0; ii<DIM; ii++)
    box[ii]=boxT[ii][ii]*A0INV;
}
