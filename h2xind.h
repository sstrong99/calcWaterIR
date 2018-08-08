#ifndef H2XIND_H
#define H2XIND_H

#include <cmath>

//convert index from 0->nH-1 to 0->natoms-1
struct H2Xind {
  virtual ~H2Xind() {};

  virtual int convert(int ind) = 0;
  //functor doesn't work well, b/c have to call like (*h2x)()
  //virtual int operator()(int ind) = 0;
};

struct SPCEind : H2Xind {
  int convert(int ind)
    { return ind + floor( (ind + 2)/2 ); };
};

struct TIP4Pind : H2Xind {
  int convert(int ind)
    { return 1-(ind%2)+2*ind; };
};

#endif
