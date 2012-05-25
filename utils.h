#ifndef STELLA_UTILS_H
#define STELLA_UTILS_H

#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>

class UtilStopWatch
{
  static struct timeval _tStart;

public:
static void tic ();
static long  toc ();
};

#endif
