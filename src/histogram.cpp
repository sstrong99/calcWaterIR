#include "histogram.h"

Histogram::Histogram(const int N,const float &lo,const float &hi) :
  N(N),lo(lo),hi(hi),dx((hi-lo)/N) {
  if (hi<=lo) {
    printf("ERROR: Range of histogram is backwards\n");
    exit(EXIT_FAILURE);
  }
  if (N<1) {
    printf("ERROR: the number of bins must be positive\n");
    exit(EXIT_FAILURE);
  }

  init();
}

Histogram::Histogram(const Histogram &hist) :
   N(hist.N),lo(hist.lo),hi(hist.hi),dx(hist.dx) {
  init();

  for (int ii=0; ii<N; ii++)
    counts[ii] = hist.counts[ii];
}

Histogram::~Histogram() {
  delete[] binCenters;
  delete[] binEdges;
  delete[] counts;
}

void Histogram::init() {
  binCenters = new float[N];
  binEdges = new float[N+1];
  counts = new float[N];

  for (int ii=0; ii<N; ii++) {
    binEdges[ii]   = lo + ii*dx;
    binCenters[ii] = binEdges[ii] + dx/2;
    counts[ii]     = 0.0;
  }
  binEdges[N] = hi;
}

//if no weight supplied, defaults to 1.0
void Histogram::addData(const float &pt,const float &weight) {
  if ( pt >= lo && pt < hi ) {
    int bin = floor( (pt-lo)/dx );
    counts[bin]+=weight;
  }
}

Histogram& Histogram::operator+=(const Histogram &rhs) {
  if (*this!=rhs) {
    printf("ERROR: Histograms being added do not have the same bins\n");
    exit(EXIT_FAILURE);
  }

  for (int ii=0; ii<N; ii++)
    counts[ii] += rhs.counts[ii];

  return *this;
}

Histogram operator+(const Histogram &a,const Histogram &b) {
  Histogram newhist(a);
  newhist+=b;
  return newhist;
}

bool Histogram::operator==(const Histogram &rhs) const {
  if (lo != rhs.lo || hi != rhs.hi || N != rhs.N )
    return false;
  return true;
}

bool Histogram::operator!=(const Histogram &rhs) const {
  return !(*this == rhs);
}

//to add an assignment operator,
//see the following explanation of copy and swap:
//https://stackoverflow.com/questions/3279543/what-is-the-copy-and-swap-idiom

void Histogram::print(const char *filename) const {
  float sum=0.0;
  for (int ii=0; ii<N; ii++) sum+=counts[ii];

  if (sum==0.0) {
    printf("WARNING: Printing an empty histogram\n");
    sum=1.0; //to prevent divide-by-zero
    //exit(EXIT_FAILURE);
  }

  FILE *myfile = fopen(filename,"w");
  for (int ii=0; ii<N; ii++)
    fprintf(myfile,"%8e %8e\n",binCenters[ii],counts[ii]/(sum*dx));
  fclose(myfile);
}
