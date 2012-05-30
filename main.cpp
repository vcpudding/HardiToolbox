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
  double s0 = 220;
  double d = 1.7e-3;
  double snr = 20;
  double angle = M_PI/3;
  int nFibers = 2;
  int nTrials = 10;
  mat fibDirs;
  vec weights;
  vec diffusivities;

  if (nFibers == 1)
  {
    fibDirs <<cos(angle) <<endr <<sin(angle) <<endr <<0 <<endr;
    weights <<1  <<endr;
    diffusivities <<1.4e-3 <<endr <<3e-4 <<endr;
  } else
  {
    fibDirs <<1 <<cos(angle) <<endr <<0 <<sin(angle) <<endr <<0 <<0 <<endr;
    weights <<0.5 <<endr <<0.5  <<endr;
    diffusivities <<1.8e-3 <<endr <<1.4e-3 <<endr <<3e-4 <<endr;
  }

  DeconvOption options;
  options.order= 16; // polynomial order
  options.delta= 20; //29.2; % estimated from the data, b-value multiplied by the dominant diffusivity
  options.lambda=1e-5; //initial damping parameter value
  options.tol= 1e-12; // convergence criteria
  options.maxIt= 3000;
  options.step = 5e-4;
  options.useAccurateIntegral = false;
  options.init = 0;
  options.useLineSearch = true;
  options.optMethod = 1;
  options.weightThres = 0.2;
  options.sparse = 0;

  StickEstimateOption stickOptions;
  stickOptions.maxIt = 550000;
  stickOptions.maxInnerIt = 5;
  stickOptions.init = 1;
  stickOptions.useManifold = true;
  stickOptions.step = 1e-8;
  stickOptions.kappaStep = 1e-13;
  stickOptions.kappa0Step = 1e-12;
  stickOptions.weightStep = 1e-12;
  stickOptions.tolerance = 1e-12;
  stickOptions.innerTolerance = 1e-6;
  stickOptions.isEstWeights = true;
  stickOptions.isEstDiffusivities = false;

  MultiTensorOption fgOptions;
  fgOptions.maxIt = 5000;
  fgOptions.step = 1e-2;
  fgOptions.tolerance = 1e-6;

  // if (argc<4) {
  //   cout <<"input x, y, and nFibers" <<endl;
  //   return 0;
  // } else {
  //   int x = atoi(argv[1]);
  //   int y = atoi(argv[2]);
  //   int nFibers = atoi(argv[3]);
  //   testTomsAlgorithmOneVoxel(x, y, nFibers, stickOptions);
  // }
  mat gradientOrientations = loadGradientOrientations("gradients.txt");
  vec S = simulateMultiTensor(bVal, 1.0, gradientOrientations, fibDirs, weights);
  vec S1 = addRicianNoise(S, 1.0/snr);
  FiberComposition fibComp;
  estimateMultiTensor(fibComp, S1, gradientOrientations, bVal, 1.0, nFibers, fgOptions);
  fibComp.fibDirs.print("estimated directions");
  fibDirs.print("true directions");
  cout <<"deviation: " <<directionDeviation(fibComp.fibDirs, fibDirs) <<endl;
  return 0;
}
