#include "../HardiToolbox.h"
#include "../utils.h"

using namespace std;
using namespace HardiToolbox;

int main (int argc, char **argv)
{
  int i=1;
  string method = string(argv[i++]);
  double snr = atof(argv[i++]);
  int degAngle = atoi(argv[i++]);
  int nFibers = atoi(argv[i++]);

  vec diffus = zeros(nFibers+1);
  for (int iDiffus=0; iDiffus<nFibers+1; ++iDiffus) {
    diffus(iDiffus) = atof(argv[i++]);
  }

  vec weights = zeros(method=="bas"?nFibers+1:nFibers);
  for (int iWeight=0; iWeight<weights.n_elem; ++iWeight) {
    if (iWeight<weights.n_elem-1)
      weights(iWeight) = atof(argv[i++]);
    else
      weights(iWeight) = 1-sum(weights);
  }

  bool bEstWeights = atoi(argv[i++]);
  bool bEstDiffus = atoi(argv[i++]);
  int nEst = atoi(argv[i++]);
  int nRep = atoi(argv[i++]);

  double angle = degAngle*M_PI/180;
  int bVal = 3000;
  double s0 = 1500;
  mat fibDirs;
  // vec weights;
  // mat diffMat;
  // vec diffVec;

  if (nFibers == 1) {
    fibDirs <<1 <<endr <<0 <<endr <<0 <<endr;
    // if (method=="bas") {
    //   weights <<w <<endr <<1-w <<endr;
    //   diffVec <<d1 <<endr <<d0 <<endr;
    // } else {
    //   weights <<1  <<endr;
    //   diffVec <<d1-d0 <<endr <<d0 <<endr;
    // }
    // diffMat <<d1 <<endr <<d0 <<endr <<d0 <<endr;
  } else if (nFibers==2) {
    fibDirs <<1 <<cos(angle) <<endr <<0 <<sin(angle) <<endr <<0 <<0 <<endr;
    // if (method=="bas") {
    //   weights <<w <<endr <<w <<endr  <<1-2*w <<endr;
    //   diffVec <<d1 <<endr <<d1 <<endr <<d0 <<endr;
    // } else {
    //   weights <<w <<endr <<1-w <<endr;
    //   diffVec <<d1-d0 <<endr <<d1-d0 <<endr <<d0 <<endr;
    // }
    // diffMat <<d1 <<d1 <<endr <<d0 <<d0 <<endr <<d0 <<d0 <<endr;
  } else {
    initFixed(fibDirs, 3);
  }

  StickEstimateOption stickOptions;
  stickOptions.maxIt = 25000;
  stickOptions.maxInnerIt = 5;
  stickOptions.init = 1;
  stickOptions.useManifold = true;
  stickOptions.step = 1e-9;
  stickOptions.kappaStep = 1e-12;
  stickOptions.kappa0Step = 1e-10;
  stickOptions.weightStep = 1e-10;
  stickOptions.tolerance = 1e-5;
  stickOptions.innerTolerance = 1e-6;
  stickOptions.isEstWeights = bEstWeights;
  stickOptions.isEstDiffusivities = bEstDiffus;
  stickOptions.useLineSearch = false;
  stickOptions.sparseFactor = 1e5;

  FiberComposition trueFibComp;
  trueFibComp.nFibers = nFibers;
  trueFibComp.fibDirs = fibDirs;
  trueFibComp.fibDiffs = diffus;
  trueFibComp.fibWeights = weights;

  mat gradientOrientations = loadGradientOrientations("gradients.txt");
  vec S;
  if (method=="bas")
    S = simulateBallAndStick(bVal, s0, gradientOrientations, fibDirs, weights, diffus);
  else
    S = simulateModBAS(bVal, s0, diffus, weights, gradientOrientations, fibDirs);
  vec S1 = addRicianNoise(S, s0/snr);

  // char path [1024];
  // if (bEstWeights && bEstDiffus)
  //   sprintf(path, "results/synthetic/%s/weights_diffus", method.c_str());
  // else if (bEstWeights)
  //   sprintf(path, "results/synthetic/%s/weights", method.c_str());
  // else
  //   sprintf(path, "results/synthetic/%s/diffus", method.c_str());

  // char fileName [2048];
  // sprintf(fileName, "%s/n=%d__s=%d__a=%d__d1=%0.1e__d0=%0.1e__w=%0.1f.txt", path, nFibers, (int)snr, degAngle, d1, d0, w);
  // ofstream outFile(fileName);
  
  for (int iRep = 0; iRep<nRep; ++iRep) {
    FiberComposition fibComp;
    double l;
    UtilStopWatch::tic();
    if (method=="bas") 
      l = estimateBallAndSticks(fibComp, S1, gradientOrientations, bVal, s0, snr, nEst, stickOptions);
    else
      l = estimateModifiedBAS(fibComp, S1, gradientOrientations, bVal, s0, snr, nEst, stickOptions);

    long timeElapsed = UtilStopWatch::toc();
    // vec dirDev, weightDev, diffusDev;
    // double err = fiberDeviation(fibComp, trueFibComp, dirDev, weightDev, diffusDev);
    // cout <<"mean error: " <<err <<endl;
    // vec resVec = join_cols(dirDev, join_cols(weightDev, diffusDev));
    fibComp.fibDirs.reshape(3*nEst, 1);
    vec resVec = join_cols(fibComp.fibDirs, join_cols(fibComp.fibWeights, fibComp.fibDiffs));
    for (int i=0; i<resVec.n_elem; ++i) {
      // outFile <<resVec(i) <<"\t";
      cout <<resVec(i) <<"\t";
    }
    // outFile <<l <<"\t" <<timeElapsed <<endl;
    cout <<l <<"\t" <<timeElapsed <<endl;
  }
  // outFile.close();
 
  return 0;
}
