#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <string>
#include <cmath>

using namespace std;
using namespace chrono;
class Timer {
public:
  Timer() { t1=high_resolution_clock::now(); };
  ~Timer() {};

  float getTime() const;
  string getDHMS(float secs) const;
  string getDHMS() const;
  void restart() { t1=high_resolution_clock::now(); };

private:
  high_resolution_clock::time_point t1;
};

#endif
