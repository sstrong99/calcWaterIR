#include "groTraj.h"

groTraj::groTraj(const std::string &trajfile)
{
  //open gro file
  ifstream infile(trajfile);
  std::string line;
  
  //skip comment line
  if (!std::getline(infile, line)) {
    printf("ERROR: gro file is not formatted correctly\n");
    exit(EXIT_FAILURE);
  }

  //get natoms
  if (!std::getline(infile, line)) {
    printf("ERROR: gro file is not formatted correctly\n");
    exit(EXIT_FAILURE);
  }
  std::istringstream iss(line);
  iss >> natoms;
  x  = new rvec[natoms];

  for (int ii=0; ii<natoms; ii++) {
    if (!std::getline(infile, line) || line.length() < 44) {
      printf("ERROR: not enough atoms in gro file\n");
      exit(EXIT_FAILURE);
    }
    x[ii][0] = stof(line.substr(20,27))*A0INV;
    x[ii][1] = stof(line.substr(28,35))*A0INV;
    x[ii][2] = stof(line.substr(36,43))*A0INV;
  }

  //get box
  if (!std::getline(infile, line)) {
    printf("ERROR: gro file is not formatted correctly\n");
    exit(EXIT_FAILURE);
  }
  std::istringstream iss1(line);
  for (int ii=0; ii<DIM; ii++) {
    iss1 >> box[ii];
    box[ii] *= A0INV;
  }

  //initialize vars
  nT   = 1;
  step = 0;
  t    = 0.0;
  dt   = 0.0;
  nextFlag=false;

  //close gro file
  infile.close();
}

groTraj::~groTraj()
{
  delete[] x;
}

//read next step of trajectory
int groTraj::next(const bool convertFlag)
{
  if (nextFlag) {
    printf("WARNING: attempt to read more than one timestep from gro file\n");
    return 1;
  }
  nextFlag=true;
  return 0;
}

float groTraj::allT() { return dt; }

void groTraj::skip(const int n)
{
  if (n>0) 
    printf("WARNING: attempt to read more than one timestep from gro file\n");
}
