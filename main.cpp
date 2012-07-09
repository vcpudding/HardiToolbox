#include <iostream>
#include <math.h>
#include <stdio.h>
#include "HardiToolbox.h"

using namespace std;
using namespace HardiToolbox;

#include "testProc.h"

int main (int argc, char **argv)
{
  // testSparsityTerm();
  // return 0;
  struct timeval tStart, tEnd;

  int bVal = 3000;
  double s0 = 1500;
  double snr = 40;
  double angle = M_PI/4;
  int nFibers = 2;
  int nTrials = 10;
  mat fibDirs;
  vec weights;
  mat diffMat;
  vec diffVec;
  double d1 = 2.0e-3;
  double d0 = 5e-4;

  if (nFibers == 1)
  {
    fibDirs <<1 <<endr <<0 <<endr <<0 <<endr;
    weights <<1  <<endr;
    diffMat <<d1 <<endr <<d0 <<endr <<d0 <<endr;
    diffVec <<d1-d0 <<endr <<d0 <<endr;
  } else
  {
    fibDirs <<1 <<cos(angle) <<endr <<0 <<sin(angle) <<endr <<0 <<0 <<endr;
    weights <<0.5 <<endr <<0.5  <<endr;
    diffMat <<d1 <<d1 <<endr <<d0 <<d0 <<endr <<d0 <<d0 <<endr;
    diffVec <<d1-d0 <<endr <<d1-d0 <<endr <<d0 <<endr;
  }

  for (int i=0; i<nFibers; ++i) {
  }

  mat R(3,3);
  initRandom(R, 3);
  //fibDirs = R*fibDirs;

  DeconvOption options;
  options.order= 16; // polynomial order
  options.delta= 20; //29.2; % estimated from the data, b-value multiplied by the dominant diffusivity
  options.lambda=1e-6; //initial damping parameter value
  options.tol= 1e-6; // convergence criteria
  options.maxIt= 3000;
  options.step = 5e-4;
  options.useAccurateIntegral = false;
  options.init = 0;
  options.useLineSearch = true;
  options.optMethod = 1;
  options.weightThres = 0.2;
  options.sparse = 0;


  // MultiTensorOption fgOptions;
  // fgOptions.maxIt = 50000;
  // fgOptions.init = 1;
  // fgOptions.step = 1e-6;
  // fgOptions.tolerance = 1e-6;

  // if (argc<4) {
  //   cout <<"input x, y and nFibers" <<endl;
  //   return 0;
  // } else {
  //   int x = atoi(argv[1]);
  //   int y = atoi(argv[2]);
  //   int nFibers = atoi(argv[3]);
  //   //testTomsAlgorithmOneVoxel(x, y, nFibers, stickOptions);
  // }
  testTomsAlgorithmOnBrainOneVoxel(argc, argv);
  // testWeightEstimation(argc, argv);
  // testBASRicianEM();
  // estimateCCDiffusivities();
  return 0;
}
