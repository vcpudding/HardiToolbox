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
  double snr = 40;
  double angle = M_PI/6;
  int nFibers = 1;
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

  mat R(3,3);
  initRandom(R, 3);
  fibDirs = R*fibDirs;

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

  StickEstimateOption stickOptions;
  stickOptions.maxIt = 5500;
  stickOptions.maxInnerIt = 5;
  stickOptions.init = 0;
  stickOptions.useManifold = true;
  stickOptions.step = 1e-7;
  stickOptions.kappaStep = 1e-13;
  stickOptions.kappa0Step = 1e-12;
  stickOptions.weightStep = 0;
  stickOptions.tolerance = 1e-5;
  stickOptions.innerTolerance = 1e-5;
  stickOptions.isEstWeights = false;
  stickOptions.isEstDiffusivities = true;
  stickOptions.useLineSearch = false;

  MultiTensorOption fgOptions;
  fgOptions.maxIt = 5000;
  fgOptions.init = 1;
  fgOptions.step = 1e-6;
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
  // testTomsAlgorithmParams(stickOptions);
  mat gradientOrientations = loadGradientOrientations("gradients.txt");
  vec S = simulateMultiTensor(bVal, s0, gradientOrientations, fibDirs, weights);
  vec S1 = addRicianNoise(S, s0/snr);
  FiberComposition fibComp;
  UtilStopWatch::tic();
  estimateMultiTensor(fibComp, S1, gradientOrientations, bVal, s0, nFibers, fgOptions);
  cout <<"time: " <<UtilStopWatch::toc() <<"ms" <<endl;
  fibComp.fibDiffs.print("estimated diffusivities");
  fibComp.fibDirs.print("estimated directions");
  fibComp.fibWeights.print("estimated weights");
  fibDirs.print("true directions");
  vec devAngle = directionDeviation(fibComp.fibDirs, fibDirs)*180/M_PI;
  devAngle.print("direction deviation");
  return 0;
}
