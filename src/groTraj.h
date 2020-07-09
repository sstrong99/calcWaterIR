#ifndef GROTRAJ_H
#define GROTRAJ_H

#include "vecManip.h"
#include "traj.h"

#include <fstream>
#include <sstream>
#include <string>

class groTraj : public Traj {
 public:
  groTraj(const std::string &trajfile);
  ~groTraj();
  int next(const bool convertFlag=true);
  float allT();
  void skip(const int n);

 private:
  bool nextFlag; //can call next once, no more
};

#endif
