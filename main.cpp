/**
    calcWaterIR
    Purpose: calculate the IR spectra of water according to Yang&Skinner 2010 PCCP
    
    @author Steven E. Strong
*/
#include "calculation.h"
#include "calcIR.h"
#include "calcDists.h"
#include "calcIR_TAA.h"
#include "timer.h"
#include "input.h"

#include <cstdio>
#include <string>

//using namespace std;
int main(int argc, const char *argv[])
{
  string inputfile="in.spec";
  if (argc == 2)
    inputfile=argv[1];
  Input input(inputfile);

  //TODO: make Traj class flexible to read lammps input
  Timer time_entire;

  //perform calculation, according to whichCalc
  Calculation *mycalc;
  switch(input.getWhichCalc()) {
  case 0 : mycalc = new CalcIR(input); break;
  case 1 : mycalc = new CalcDists(input); break;
  case 2 : mycalc = new CalcIR_TAA(input); break;
  default :
    printf("ERROR: calc keyword accepts 0=IR,1=dists,2=TAA\n");
    return EXIT_FAILURE;
  }
  mycalc->printResults(input.getOutPostfix());
  delete mycalc;

  string time = time_entire.getDHMS();
  printf("Completed in %s\n",time.c_str());

  return EXIT_SUCCESS;
};
