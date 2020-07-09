#ifndef INPUT_H
#define INPUT_H

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;
class Input {
public:
  Input(int argc, const char *argv[]);
  ~Input() {};

  void readInputFile(const string &inputfile);
  string getTrajFile() const {return trajFile;};
  string getOutPostfix() const {return outPostfix;};
  int getIntMethod() const {return intMethod;};
  int getNsample() const {return nSample;};
  float getTimestep() const {return timestep;};
  float getAvgF() const {return avgF;};
  int getNtcf() const {return nTCF;};
  int getWhichCalc() const {return whichCalc;};
  float getT1() const {return T1;};
  float getTavg() const {return Tavg;};

private:
  string trajFile;
  string outPostfix;
  int intMethod;  //0 = exactDiag,+integer=Adams-Bashforth nth order
  int nSample;    //0 = max,-1=don't reuse data,+integer=literal
  float timestep; //0 = use traj timestep
  float avgF;    //0 = estimate from one snapshot
  int nTCF;       //0 = default of 6*T1 relaxation time
  int whichCalc; //0=IR, 1=dists, 2=TAA
  float T1;
  float Tavg;     //for TAA approach
};
#endif
