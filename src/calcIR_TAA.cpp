#include "calcIR_TAA.h"

CalcIR_TAA::CalcIR_TAA(const Input &inp) :
  T1(inp.getT1()),xtcfile(inp.getTrajFile()),timestep(inp.getTimestep())
{
  init(inp.getTavg());

  int nSample = inp.getNsample();
  int sampleStep = getNsample(nT,nTavg,nSample);

  //loop over samples
  loopSamples(nSample,sampleStep);
}

CalcIR_TAA::~CalcIR_TAA()
{
  delete[] spectrum;
  delete[] w;
}

void CalcIR_TAA::loopSamples(const int nSample, const int step)
{
  int nFinished=0;
  int maxThreads=1;

#ifdef USEOMP
  int nThreads = ompNumThreads();
  maxThreads = nSample < nThreads ? nSample : nThreads;
#endif
  //omp_set_num_threads(maxThreads);
  printf("Looping over %d samples using %d threads...\n",
	 nSample,maxThreads);

  Timer totalTime;

#pragma omp parallel for num_threads(maxThreads)
  for (int ii=0; ii<nSample; ii++)
  {
    Timer timer;

    float *spec=new float[N];
    calcAvg(ii*step,spec);

    //this is not good parallelism, but this section is very fast relative to
    //the rest of the loop, so should be ok
#pragma omp critical
    {
      //reduction
      for (int kk=0; kk<N; kk++)
	spectrum[kk]+=spec[kk];
      delete[] spec;

      nFinished++;
      if ( nFinished % (maxThreads*5) == 0 ) {
	printf("\rCompleted %d of %d samples (%.1f%%) in %s - %.2fs/timestep",
	       nFinished,nSample,100*(float)nFinished/(float)nSample,
	       (totalTime.getDHMS()).c_str(),timer.getTime()/nTavg);
	fflush(stdout);
      }
    }
  }
  printf("\n");

  for (int jj=0; jj<N; jj++)
    spectrum[jj] /= ((float) nSample);
}

//preloading the traj would only work for very short and
//small sims, so i removed that code
void CalcIR_TAA::calcAvg(const int start,float* spec1)
{
  //initialize
  Traj *traj = Traj::getTraj(xtcfile);
  traj->skip(start);

  CalcW calcW(model,traj->getNatoms(),0.0);

  int nH2=nH*nH;
  rvec *m = new rvec[nH];
  rvec *m0 = new rvec[nH];
  float *avgW = new float[nH2];
  int ii;
  for (ii=0; ii<nH2; ii++)
    avgW[ii]=0.0;

  //loop through time
  for (int tt=0; tt<nTavg; tt++) {
    //get next timestep of coords
    if (traj->next()) {
      printf("ERROR: no more coordinates in file\n");
      exit(EXIT_FAILURE);
    }

    //compute w matrix and m vectors at this timestep
    calcW.compute(traj,m);
    if (tt==0)
      for (ii=0; ii<nH; ii++)
	getRvec(m[ii],m0[ii]);

    //add W matrix to average
    for (ii=0; ii<nH2; ii++)
      avgW[ii]+=calcW.getW()[ii];

    if (dt_skip > 1)
      traj->skip(dt_skip-1);
  }
  //printf("inter,intra= %.1f\t%.1f\n",maxInter,avgIntra);

  for (ii=0; ii<nH2; ii++)
    avgW[ii]/=nTavg;

  //diagonalize average kappa
  ExactDiag diag(nH,timestep);
  diag.diag(avgW);
  float *d0 = new float[nH];
  float *d = new float[nH];
  diag.mult_taa(m0,d0);
  diag.mult_taa(m,d);

  delete[] m;
  delete[] m0;
  delete[] avgW;
  delete traj;
  
  float *eig = new float[nH];
  for (ii=0; ii<nH; ii++)
    eig[ii]=diag.getEigenvalue(ii)/CM2PS;

  int jj;
  float sum,wi,lor;
  for (ii=0; ii<N; ii++) {
    sum=0.0;
    wi=w[ii];
    for (jj=0; jj<nH; jj++) {
      lor=lorentz(wi-eig[jj]);
      sum+=d[jj]*d0[jj]*lor;
    }
    spec1[ii]=sum;
  }

  delete[] eig;
  delete[] d;
  delete[] d0;
}

inline float CalcIR_TAA::lorentz(const float &w) {
  return 0.5*PI*T1/(w*w + 0.25*T1*T1);
}

void CalcIR_TAA::printResults(string postfix) const {
  print(postfix.insert(0,"spec").c_str(),w,spectrum,N);
}

void CalcIR_TAA::init(const float &Tavg) {
  //initialize trajectory variables
  InitTraj traj(xtcfile);
  dt_skip=traj.adjustTimestep(timestep);
  timestep=traj.getDT();
  nT=traj.getNT();
  model = traj.getModel();
  nH=traj.getNH();

  traj.printModel();

  //initialize spectrum
  N=1000;
  float wStart=2900; //cm-1
  float wStop=3900;
  float dw=(wStop-wStart)/(N-1);
  w=new float[N];
  spectrum=new float[N];
  for (int ii=0; ii<N; ii++) {
    w[ii]=wStart+ii*dw;
    spectrum[ii]=0.0;
  }

  //initialize averaging time
  float tmpNT=Tavg/timestep;
  if (fabs(round(tmpNT)-tmpNT) > 1e-3)
    printf("WARNING: The averaging time has been adjusted to be a multiple of the timestep.\n");

  nTavg=(int) tmpNT;

  printf("There are %d steps with a timestep of %.1f fs.\n",nT,1000*timestep);
  printf("The averaging time is %.1f fs.\n",1000*nTavg*timestep);
}
