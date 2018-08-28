#include "calcIR.h"

extern "C" {
  extern int sger_(int*,int*,float*,float*,int*,float*,int*,float*,int*);
}

CalcIR::CalcIR(const Input &inp) :
  T1(inp.getT1()),xtcfile((inp.getTrajFile())),
  integrator(inp.getIntMethod()),timestep(inp.getTimestep()),
  nTCF(inp.getNtcf()),avgF(inp.getAvgF())
{
  init();

  int nSample = inp.getNsample();
  int sampleStep = getNsample(nT,nTCF,nSample);

  //loop over samples to get TCF
  avgCorr=new cpx[nTCF];
  loopSamples(nSample,sampleStep);

  //FFT TCF to get I(w)
  printf("Computing FFT...\n");
  calcFFT(avgCorr);
}

CalcIR::~CalcIR()
{
  delete[] spectrum;
  delete[] w;
  delete[] avgCorr;
}

void CalcIR::loopSamples(const int nSample, const int step)
{
  vector<cpx*> corr(nSample);  //initialize size of corr vector
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

    corr[ii]=new cpx[nTCF];
    calcTCF(ii*step,corr[ii]);

    //this is not good parallelism, but this section is very fast relative to
    //the rest of the loop, so should be ok
#pragma omp critical
    {
      nFinished++;
      if ( nFinished % maxThreads == 0 ) {
	printf("\rCompleted %d of %d samples (%.1f%%) in %s - %.2fs/timestep",
	       nFinished,nSample,100*(float)nFinished/(float)nSample,
	       (totalTime.getDHMS()).c_str(),timer.getTime()/nTCF);
	fflush(stdout);
      }
    }
  }
  printf("\n");

  //average correlation function over samples
  const float T1x2=2*T1;
  cpx sum;
  int ii;
  for (int jj=0; jj<nTCF; jj++) {
    sum=0.0;
    for (ii=0; ii<nSample; ii++)
      sum+=corr[ii][jj];

    //NOTE: Nick's code doesnt normalize by nSample until after printing TCF
    avgCorr[jj]=sum * exp(-jj*timestep/T1x2) / ((float) nSample);
  }

  for (uint kk=0; kk<corr.size(); kk++)
    delete[] corr[kk];
}

//preloading the traj would only work for very short and
//small sims, so i removed that code
void CalcIR::calcTCF(const int start,cpx* corr1)
{
  //initialize
  Traj traj(xtcfile.c_str());
  traj.skip(start);

  CalcW calcW(model,traj.getNatoms(),avgF);

  IntegrateF *intF;
  if (integrator > 0)
    intF = new AdamsBashforth(integrator,nH,timestep);
  else
    intF = new ExactDiag(nH,timestep);

  cpx *F = new cpx[nH*nH];
  intF->initF(F,nH);
  rvec *m = new rvec[nH];
  rvec *m0= new rvec[nH];

  //loop through time
  float maxInter=0.0;
  float avgIntra=0.0;
  int ii;
  for (int tt=0; tt<nTCF; tt++) {
    //get next timestep of coords
    if (traj.next()) {
      printf("ERROR: no more coordinates in file\n");
      exit(EXIT_FAILURE);
    }

    //printf("norm(F) = %f\n",norm(F));

    //compute w matrix and m vectors at this timestep
    calcW.compute(traj,m);
    maxInter+=calcW.maxInterFreq();
    avgIntra+=calcW.avgIntraFreq();
    if (tt==0)
      for (ii=0; ii<nH; ii++)
	getRvec(m[ii],m0[ii]);
    else
      intF->next(F,calcW.getW());  //integrate F matrix forward one step
    //integrate F with W matrix at current time, not previous time
    //Fdot(t) = i F(t) W(t)

    //compute correlation function
    corr1[tt]=sumMFM(m0,m,F);

    if (dt_skip > 1)
      traj.skip(dt_skip-1);
  }
  maxInter/=nTCF;
  avgIntra/=nTCF;
  //printf("inter,intra= %.1f\t%.1f\n",maxInter,avgIntra);

  delete[] m0;
  delete[] m;
  delete[] F;
  delete intF;
}

void CalcIR::calcFFT(cpx *y)
{
  int ii;
  //pad y with zeros
  //mult input by (-1)^ii to put zero freq at center
  N = nextPow2(nTCF)*4; //can add more zeros to interpolate more finely
  cpx fact = 1;
  cpx *yPad = new cpx[N];
  for (ii=0; ii<nTCF; ii++) {
    yPad[ii]=y[ii] * fact;
    fact*=-1;
  }
  for (ii=nTCF; ii<N; ii++)
    yPad[ii]=0.0;

  //c2r effectively mirrors yPad about t=0, and returns abs of FFT
  spectrum = new float[N];
  fftwf_plan p = fftwf_plan_dft_c2r_1d(N, reinterpret_cast<fftwf_complex*>(yPad),
				     spectrum,FFTW_ESTIMATE);
  fftwf_execute(p);
  fftwf_destroy_plan(p);
  delete[] yPad;

  //output frequencies are in cycle/ps, so convert to radians, then wavenumber
  w=new float[N];
  float freqFactor=2*PI/(N*timestep*CM2PS);
  float intFactor = 1.0/(freqFactor*N);
  for (ii=0; ii<N; ii++) {
    spectrum[ii] *= intFactor;  //convert to energy
    w[ii]         = avgF - (ii-N/2)*freqFactor; //negate frequency, b/c c2r fft is actually inverse fft
  }
}

//from https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
int CalcIR::nextPow2(int v)
{
v--;
v |= v >> 1;
v |= v >> 2;
v |= v >> 4;
v |= v >> 8;
v |= v >> 16;
return ++v;
}

cpx CalcIR::sumMFM(const rvec *m0,const rvec *m,const cpx *F) {
  cpx sum=0.0;
  int ii,jj;
  //could do as 1d sum, maybe faster
  for (ii=0; ii<nH; ii++)
    for (jj=0; jj<nH; jj++) {
      sum += F[jj+ii*nH] * (float) dot(m0[ii],m[jj]);
      //magma handles matricies as col-major, so F is transposed.
      //instead of transposing, switch m0 and m indicies
      //if ODE is Fdot=prop F, then should be m0[jj],m[ii],
      //but I switch the order here
      //TODO: check that adams-bashforth is compatible with this
    }
  //NOTE: Nick's code doesn't divide by 3
  //return sum/((float) DIM);  //average over 3 dims
  return sum;
}

float CalcIR::norm(const cpx *mat)
{
  float sum=0.0;
  int ii,jj;
  for (ii=0; ii<nH; ii++)
    for (jj=0; jj<nH; jj++)
      sum += real( mat[jj+ii*nH]*conj(mat[jj+ii*nH]) );

  return sqrt(sum/nH);
}

void CalcIR::printResults(string postfix) const {
  //print spectrum
  string postfixOrig=postfix;
  print(postfix.insert(0,"spec").c_str(),w,spectrum,N);

  //print TCF
  float *t=new float[nTCF];
  float *realTCF=new float[nTCF];
  float *imagTCF=new float[nTCF];
  for (int ii=0; ii<nTCF; ii++)
  {
    t[ii]=ii*timestep;
    realTCF[ii]=real(avgCorr[ii]);
    imagTCF[ii]=imag(avgCorr[ii]);
  }
  print(postfixOrig.insert(0,"tcf").c_str(), t, realTCF, imagTCF, nTCF);

  delete[] t;
  delete[] realTCF;
}

void CalcIR::init() {
  InitTraj traj(xtcfile.c_str());
  dt_skip=traj.adjustTimestep(timestep);
  timestep=traj.getDT();  //this is the adjusted timestep, not the original
  nT=traj.getNT();
  model = traj.getModel();
  nH=traj.getNH();

  if (avgF==0.0)
    avgF = traj.getAvgF();

  traj.printModel();

  //chose nTCF based on damping time T1
  //damping term e^(-t/2*T1) ~ 0.05 for t=6*T1
  if (nTCF == 0)
    nTCF=(int) (T1 * 6 / timestep);
  if (nTCF < 2) {
    printf("ERROR: T1 = %f is bad\n",T1);
    exit(EXIT_FAILURE);
  }

  printf("There are %d steps with a timestep of %.1f fs.\n",nT,1000*timestep);
  printf("The time correlation function will be %.1f ps long.\n",nTCF*timestep);
  printf("Removing average frequency %.1f cm^-1\n",avgF);
}
