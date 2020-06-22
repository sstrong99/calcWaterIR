#include "timer.h"

float Timer::getTime() const {
    high_resolution_clock::time_point t2=high_resolution_clock::now();
    duration<float> time_span;
    time_span=duration_cast<duration<float>>(t2 - t1);
    return time_span.count();
  };

string Timer::getDHMS() const {
  float time = getTime();
  return getDHMS(time);
};

string Timer::getDHMS(float secs) const {
  int mins=0;
  int hrs=0;
  int days=0;
  string str;
  char tmp[4]; //extra char for null terminator
  mins=floor(secs/60);
  if (mins > 0) {
    secs = secs - mins*60;
    hrs = floor(mins/60);
    if (hrs > 0) {
      mins = mins - hrs*60;
      days = floor(hrs/24);
      if (days > 0) {
	hrs = hrs - days*24;
	str.append(to_string(days));
	str.append("-");
      }
      //only print hrs if non-zero, but print mins:secs no matter what
      str.append(to_string(hrs));
      str.append(":");
    }
  }
  sprintf(tmp,"%02d:",mins);
  str.append(tmp);

  sprintf(tmp,"%02d",(int) round(secs));
  str.append(tmp);

  return str;
}
