#include "calcLocal.h"

CalcLocal::CalcLocal(const Input &inp) {
  InitTraj init(inp.getTrajFile());
  int dt_skip=init.adjustTimestep(inp.getTimestep());
  int nT = init.getNT();
  int nSample = inp.getNsample();
  nSample = nSample < nT && nSample > 0 ? nSample : nT;
  int model = init.getModel();
  init.printModel();
  int maxThreads=1;
#ifdef USEOMP
  int nThreads = ompNumThreads();
  maxThreads = nSample < nThreads ? nSample : nThreads;
#endif
  int tPerThread=floor(nSample/maxThreads);
  int nFinished=0;

  Timer timer;

  localFreqs.resize(nT);
  transStrength.resize(nT);
  printf("Computing local frequencies using %d threads...\n",maxThreads);
#pragma omp parallel num_threads(maxThreads)
  {
    int startInd = omp_get_thread_num()*tPerThread;
    Traj *traj = Traj::getTraj(inp.getTrajFile());
    traj->skip(startInd);

    CalcW calcW(model,traj->getNatoms(),0.0);
    int nH=calcW.getNH();

    //ExactDiag diag(nH,0.0);

    rvec *m = new rvec[nH];

    for (int tt=startInd; tt<startInd+tPerThread; tt++) {
      traj->next();
      calcW.compute(traj,m);

      localFreqs[tt].resize(nH);
      transStrength[tt].resize(nH);
      for (int ii=0; ii<nH; ii++) {
	localFreqs[tt][ii] = calcW.getW(ii,ii)/CM2PS;
	transStrength[tt][ii] = sqrt(norm2vec(m[ii]));

	//tmpW=diag.getEigenvalue(ii)/CM2PS;
      }

#pragma omp critical
      {
	nFinished++;
	if ( nFinished % 1000 == 0 ) {
	  printf("\rCompleted %d of %d timesteps (%.1f%%) in %s - %.2fs/timestep",
		 nFinished,nT,100*(float)nFinished/(float)nT,
		 (timer.getDHMS()).c_str(),timer.getTime()/nFinished);
	  fflush(stdout);
	}
      }

      if (dt_skip>1)
	traj->skip(dt_skip-1);
    } //end of parallel for
    delete[] m;
    delete traj;

  }
  printf("\n");
}

void CalcLocal::printResults(string postfix) const {
  string tmp=postfix;
  string fname=tmp.insert(0,"localFreq");

  FILE * pFile;
  pFile = fopen (fname.c_str(),"w");
  fprintf(pFile,"#local frequency (cm-1)\ttransition dipole magnitude (atomic units)\n");
  
  for (uint tt=0; tt<localFreqs.size(); tt++) {
    fprintf(pFile,"#step = %d\n",tt);
    for (uint ii=0; ii<localFreqs[tt].size(); ii++) {
      fprintf(pFile,"%f\t%f\n",localFreqs[tt][ii],transStrength[tt][ii]);
    }
  }
   fclose (pFile);
}
