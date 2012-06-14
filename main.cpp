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
  double s0 = 1;
  double d = 1.7e-3;
  double snr = 20;
  double angle = M_PI/6;
  int nFibers = 2;
  int nTrials = 10;
  mat fibDirs;
  vec weights;
  mat diffusivities;

  if (nFibers == 1)
  {
    fibDirs <<1 <<endr <<0 <<endr <<0 <<endr;
    weights <<1  <<endr;
    diffusivities <<1.4e-3 <<endr <<3e-4 <<endr <<3e-4 <<endr;
  } else
  {
    fibDirs <<1 <<cos(angle) <<endr <<0 <<sin(angle) <<endr <<0 <<0 <<endr;
    weights <<0.5 <<endr <<0.5  <<endr <<0 <<endr;
    diffusivities <<1.8e-3 <<1.6e-3 <<endr <<0.3e-3 <<0.3e-3 <<endr <<3e-4 <<3e-4 <<endr;
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
  stickOptions.maxIt = 15000;
  stickOptions.maxInnerIt = 5;
  stickOptions.init = 0;
  stickOptions.useManifold = true;
  stickOptions.step = 1e-7;
  stickOptions.kappaStep = 1e-12;
  stickOptions.kappa0Step = 1e-11;
  stickOptions.weightStep = 1e-7;
  stickOptions.tolerance = 1e-6;
  stickOptions.innerTolerance = 1e-5;
  stickOptions.isEstWeights = false;
  stickOptions.isEstDiffusivities = true;
  stickOptions.useLineSearch = false;

  MultiTensorOption fgOptions;
  fgOptions.maxIt = 50000;
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
  // testTomsAlgorithmOnPhantom(stickOptions, "results/phantom/weights");
  // testWeightEstimation(stickOptions);
  mat gradientOrientations = loadGradientOrientations("gradients.txt");
  vec S = simulateBallAndStick(bVal, s0, gradientOrientations, fibDirs, weights);
  vec S1 = addRicianNoise(S, s0/snr);
  FiberComposition fibComp;
  UtilStopWatch::tic();
  //estimateMultiTensor(fibComp, S1, gradientOrientations, bVal, s0, nFibers, fgOptions);
  estimateBallAndSticks(fibComp, S1, gradientOrientations, bVal, s0, snr, nFibers, stickOptions);
  cout <<"time: " <<UtilStopWatch::toc() <<"ms" <<endl;
  fibComp.fibDiffs.print("estimated diffusivities");
  fibComp.fibDirs.print("estimated directions");
  fibComp.fibWeights.print("estimated weights");
  fibDirs.print("true directions");
  vec devAngle = directionDeviation(fibComp.fibDirs, fibDirs)*180/M_PI;
  devAngle.print("direction deviation");
  return 0;
}
