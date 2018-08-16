#include "input.h"

Input::Input(const string &inputfile) : outPostfix(""),intMethod(0),nSample(-1),
					timestep(0.0),avgF(0.0),nTCF(0),
					whichCalc(0),T1(0.26),Tavg(0.076)
{
  ifstream file(inputfile);
  if (!file.good()) {
    printf("ERROR: Input file %s cannot be read.\n",inputfile.c_str());
    exit(EXIT_FAILURE);
  }

  string   line;
  string   key;

  while(getline(file, line))
  {
    //check if line is empty
    if (line.empty())
      continue;

    stringstream   linestream(line);

    linestream >> key;
    if (key.compare("calc")==0)
      linestream >> whichCalc;
    else if (key.compare("trajFile")==0)
      linestream >> trajFile;
    else if (key.compare("outPostfix")==0) {
      linestream >> outPostfix;
      outPostfix.insert(0,"_");
    }
    else if (key.compare("intMethod")==0)
      linestream >> intMethod;
    else if (key.compare("nSample")==0)
      linestream >> nSample;
    else if (key.compare("timestep")==0)
      linestream >> timestep;
    else if (key.compare("avgF")==0)
      linestream >> avgF;
    else if (key.compare("nTCF")==0)
      linestream >> nTCF;
    else if (key.compare("T1")==0)
      linestream >> T1;
    else if (key.compare("Tavg")==0)
      linestream >> Tavg;
    else if (key.compare(0,1,"#")==0)
      continue; //skip comment
    else {
      printf("ERROR: unrecognized keyword %s in %s\n",
	     key.c_str(),inputfile.c_str());
      exit(EXIT_FAILURE);
    }
  }

  //TODO: add ability to skip comments

  file.close();

  if (whichCalc < 0 && whichCalc > 2) {
    printf("ERROR: calc = %d is invalid.\n",whichCalc);
    exit(EXIT_FAILURE);
  }

  if (Tavg<=0.0) {
    printf("ERROR: Tavg = %f is invalid.\n",Tavg);
    exit(EXIT_FAILURE);
  }

  if (intMethod < 0 || intMethod > 4) {
    printf("ERROR: intMethod = %d is invalid.\n",intMethod);
    exit(EXIT_FAILURE);
  }

  if (nSample < -1) {
    printf("ERROR: nSample = %d is invalid.\n",nSample);
    exit(EXIT_FAILURE);
  }

  if (timestep < 0.0) {
    printf("ERROR: timestep = %f is invalid.\n",timestep);
    exit(EXIT_FAILURE);
  }

  if (avgF < 0.0) {
    printf("ERROR: avgF = %f is invalid.\n",avgF);
    exit(EXIT_FAILURE);
  }

  if (nTCF < 0) {
    printf("ERROR: nTCF = %d is invalid.\n",nTCF);
    exit(EXIT_FAILURE);
  }

  if (intMethod>0 && timestep>0.005)
    printf("WARNING: Adams-Bashforth might not work well with a large timestep\n");

}
