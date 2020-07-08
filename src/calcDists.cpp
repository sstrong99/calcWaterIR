#include "calcDists.h"

CalcDists::CalcDists(const Input &inp) :
  Pu_tot(1000,2900,3900),Pc_tot(Pu_tot),spdn_tot(Pu_tot),Pintra_tot(200,-60,10)
{
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

  printf("Computing frequency distributions using %d threads...\n",maxThreads);
#pragma omp parallel num_threads(maxThreads)
  {
    Traj *traj = Traj::getTraj(inp.getTrajFile());
    traj->skip(omp_get_thread_num()*tPerThread);

    CalcW calcW(model,traj->getNatoms(),0.0);
    int nH=calcW.getNH();

    ExactDiag diag(nH,0.0);

    rvec *m = new rvec[nH];
    float *weights = new float[nH];
    float tmpW;

    Histogram Pu(Pu_tot);
    Histogram Pc(Pc_tot);
    Histogram spdn(spdn_tot);
    Histogram Pintra(Pintra_tot);

    for (int tt=0; tt<tPerThread; tt++) {
      traj->next();

      calcW.compute(traj,m);
      diag.spdn(calcW.getW(),m,weights);

      for (int ii=0; ii<nH; ii++) {
	Pu.addData(calcW.getW(ii,ii)/CM2PS);

	tmpW=diag.getEigenvalue(ii)/CM2PS;
	Pc.addData(tmpW);
	spdn.addData(tmpW,weights[ii]); //spectral density
      }

      //histogram intramolecular couplings
      for (int ii=0; ii<nH/2; ii++) {
	tmpW=calcW.getW(ii*2,ii*2+1)/CM2PS;
	Pintra.addData(tmpW);
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
    delete[] weights;
    delete traj;

#pragma omp critical
    {
      Pu_tot  += Pu;
      Pc_tot  += Pc;
      spdn_tot += spdn;
      Pintra_tot += Pintra;
    }
  }
  printf("\n");
}

void CalcDists::printResults(string postfix) const {
  string tmp=postfix;
  Pu_tot.print(tmp.insert(0,"Pu").c_str());

  tmp=postfix;
  Pc_tot.print(tmp.insert(0,"Pc").c_str());

  tmp=postfix;
  Pintra_tot.print(tmp.insert(0,"Pintra").c_str());

  tmp=postfix;
  spdn_tot.print(tmp.insert(0,"spdn").c_str());
}
