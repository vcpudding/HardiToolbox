#include "utils.h"

timeval UtilStopWatch::_tStart = {0,0};

void UtilStopWatch::tic ()
{
  gettimeofday(&(UtilStopWatch::_tStart), NULL);
}

long UtilStopWatch::toc ()
{
  long mtime, seconds, useconds;    

  struct timeval tEnd;
  gettimeofday(&tEnd, NULL);

  seconds  = tEnd.tv_sec  - UtilStopWatch::_tStart.tv_sec;
  useconds = tEnd.tv_usec - UtilStopWatch::_tStart.tv_usec;

  mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;

  return mtime;
}
