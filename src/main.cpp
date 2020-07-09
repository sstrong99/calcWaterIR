/**
    calcWaterIR
    Purpose: calculate the IR spectra of water according to Yang&Skinner 2010 PCCP
    
    @author Steven E. Strong
*/
#include "calculation.h"

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
  Calculation *mycalc = Calculation::createCalc(input);
  mycalc->printResults(input.getOutPostfix());
  delete mycalc;

  string time = time_entire.getDHMS();
  printf("Completed in %s\n",time.c_str());

  return EXIT_SUCCESS;
};
