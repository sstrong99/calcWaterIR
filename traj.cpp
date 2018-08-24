#include "traj.h"

using namespace std;
Traj::Traj(const char *xtcfile)
{
  //get number of atoms
  read_xtc_natoms((char*)xtcfile,&natoms);
  //prnt.setNH(natoms);

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

Traj::~Traj()
{
  delete[] x;
  xdrfile_close(trj);
}

//read next step of trajectory
int Traj::next(const bool convertFlag)
{
  //read one timestep, returns 0 if success
  int stat=read_xtc(trj,natoms,&step,&t,boxT,x,&prec);

  //return if finished reading file
  if (stat) return stat;

  if (convertFlag) convert();

  nT++;

  return stat;
}

float Traj::allT()
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

void Traj::skip(const int n)
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

void Traj::convert()
{
  //for (int jj=0; jj<natoms; jj++)
  //  for (int kk=0; kk<DIM; kk++)
  //    x[jj][kk]*=A0INV;

  for (int ii=0; ii<DIM; ii++)
    box[ii]=boxT[ii][ii];//*A0INV;
}

void Traj::getBox(rvec &out) const {
  for (int kk=0; kk<DIM; kk++)
    out[kk]=box[kk];
}

void Traj::moveM(const float &frac,const int aPerM) {
  //assume that O is first and M is last atom
  int nO=natoms/aPerM;
  int Oind,Mind;
  rvec O,M,vec;
  for (int ii=0; ii<nO; ii++) {
    Oind=ii*aPerM;
    Mind=Oind+aPerM-1;
    getRvec(x[Oind],O);
    getRvec(x[Mind],M);
    addRvec(M,O,vec,-1);
    multRvec(vec,frac);
    addRvec(O,vec,x[Mind],+1);
  }
}

//This is inefficient and not meant to be used more than once on a file
//figure out which water model the trajectory is using
//0=SPCE, 1=TIP4P, 2=TIP4P2005
int Traj::getModel() const {
  //1st atom should be O
  //next 2 atoms should be Hs
  rvec O,H1,H2;
  getRvec(x[0],O);
  getRvec(x[1],H1);
  getRvec(x[2],H2);

  rvec OH1,OH2;
  addRvec(H1,O,OH1,-1);
  addRvec(H2,O,OH2,-1);
  pbc(OH1,box);
  pbc(OH2,box);
  float dOH1=sqrt(norm2vec(OH1));
  float dOH2=sqrt(norm2vec(OH2));
  multRvec(OH1,1.0/dOH1);
  multRvec(OH2,1.0/dOH2);

  //distinguish between 3-site and 4-site
  float angle = acos(dot(OH1,OH2))*180/PI;
  int model;
  const float SPCEangle = 109.47;
  const float TIP4Pangle = 104.52;
  const float angleTol = 0.1;
  if ( fabs(angle - SPCEangle) < angleTol ) {
    model = 0;
  } else if ( fabs(angle - TIP4Pangle ) < angleTol ) {
    model = 1;
  } else {
    printf("ERROR: bond angle does not match any known model.\n");
    exit(EXIT_FAILURE);
  }

  //distinguish between TIP4P and TIP4P/2005
  if (model == 1) {
    rvec M,OM;
    getRvec(x[3],M);
    addRvec(M,O,OM,-1);
    pbc(OM,box);
    float dOM=sqrt(norm2vec(OM));

    const float dOM_t4p=0.015; //*A0INV;
    const float dOM_2005=0.01546; //*A0INV;
    const float dOMtol = 0.001/A0INV; //in A0 units
    if (fabs(dOM-dOM_t4p) < dOMtol)
      model=1;
    else if (fabs(dOM-dOM_2005) < dOMtol)
      model=2;
    else {
    printf("ERROR: OM bond length does not match any known model.\n");
    exit(EXIT_FAILURE);
    }
  }

  //float check OH bond length
  const float dOH_spce = 0.1; //*A0INV;
  const float dOH_t4p = 0.09572; //*A0INV;
  const float dOHtol = 0.01/A0INV; //in A0 units
  if (model == 0) {
    if (fabs(dOH1 - dOH_spce) > dOHtol ) {
      printf("ERROR: OH bond length does not match SPC/E.\n");
      exit(EXIT_FAILURE);
    }
    if (fabs(dOH2 - dOH_spce) > dOHtol ) {
      printf("ERROR: OH bond length does not match SPC/E.\n");
      exit(EXIT_FAILURE);
    }
  } else {
    if (fabs(dOH1 - dOH_t4p) > dOHtol ) {
      printf("ERROR: OH bond length does not match TIP4P-type model.\n");
      exit(EXIT_FAILURE);
    }
    if (fabs(dOH2 - dOH_t4p) > dOHtol ) {
      printf("ERROR: OH bond length does not match TIP4P-type model.\n");
      exit(EXIT_FAILURE);
    }
  }

  return model;
}
