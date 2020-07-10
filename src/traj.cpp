#include "traj.h"

//must include these here, not in header, because otherwise Traj isn't defined
//so children won't know what parent is
#include "xtcTraj.h" 
#include "groTraj.h"

using namespace std;
//factory method
Traj *Traj::getTraj(const string &filename) {
  uint extPos = filename.find_last_of('.');
  if (extPos==filename.length()) {
    printf("ERROR: trajectory file does not have an extension\n");
    exit(EXIT_FAILURE);
  }

  string ext=filename.substr(extPos+1);
  if (ext.compare("xtc") == 0) {
    return new XTCTraj(filename.c_str());
//} else if (ext.compare("lammpstrj") == 0) {
    //  return new lammpsTraj(filename);
  } else if (ext.compare("gro") == 0) {
    return new groTraj(filename);
  } else {
    printf("ERROR: trajectory extension \"%s\" is not supported\n",ext.c_str());
    exit(EXIT_FAILURE);
  }
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

  //the below comparisons need high tolerances for the low precision .gro files

  //distinguish between 3-site and 4-site
  float angle = acos(dot(OH1,OH2))*180/PI;
  int model;
  const float SPCEangle = 109.47;
  const float TIP4Pangle = 104.52;
  const float angleTol = 1.5;
  //printf("angle = %f\nOH1 = %f\nOH2 = %f\n",angle,dOH1,dOH2);
  if ( fabs(angle - SPCEangle) < angleTol ) {
    model = 0;
  } else if ( fabs(angle - TIP4Pangle ) < angleTol ) {
    model = 1;
  } else {
    printf("ERROR: bond angle %f does not match any known model.\n",angle);
    exit(EXIT_FAILURE);
  }

  //distinguish between TIP4P and TIP4P/2005
  if (model == 1) {
    rvec M,OM;
    getRvec(x[3],M);
    addRvec(M,O,OM,-1);
    pbc(OM,box);
    float dOM=sqrt(norm2vec(OM));
    //printf("OM = %f\n",dOM);
    const float dOM_t4p=0.015*A0INV;
    const float dOM_2005=0.01546*A0INV;
    const float dOMtol = 0.02*A0INV;
    if (fabs(dOM-dOM_t4p) < dOMtol)
      model=1;
    else if (fabs(dOM-dOM_2005) < dOMtol)
      model=2;
    else {
      printf("ERROR: OM bond length %f does not match any known model.\n",dOM/A0INV);
      exit(EXIT_FAILURE);
    }
  }

  //float check OH bond length
  const float dOH_spce = 0.1*A0INV;
  const float dOH_t4p = 0.09572*A0INV;
  const float dOHtol = 0.002*A0INV;
  if (model == 0) {
    if (fabs(dOH1 - dOH_spce) > dOHtol ) {
      printf("ERROR: OH bond length %f does not match SPC/E.\n",dOH1/A0INV);
      exit(EXIT_FAILURE);
    }
    if (fabs(dOH2 - dOH_spce) > dOHtol ) {
      printf("ERROR: OH bond length %f does not match SPC/E.\n",dOH2/A0INV);
      exit(EXIT_FAILURE);
    }
  } else {
    if (fabs(dOH1 - dOH_t4p) > dOHtol ) {
      printf("ERROR: OH bond length %f does not match TIP4P-type model.\n",dOH1/A0INV);
      exit(EXIT_FAILURE);
    }
    if (fabs(dOH2 - dOH_t4p) > dOHtol ) {
      printf("ERROR: OH bond length %f does not match TIP4P-type model.\n",dOH2/A0INV);
      exit(EXIT_FAILURE);
    }
  }

  return model;
}
