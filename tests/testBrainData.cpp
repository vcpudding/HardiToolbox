#include "../HardiToolbox.h"
#include "../testProc.h"

using namespace std;
using namespace HardiToolbox;

int main (int argc, char **argv)
{
  int iArg = 1;
  int x = atoi(argv[iArg++]);
  int y = atoi(argv[iArg++]);
  int z = atoi(argv[iArg++]);
  string method = string(argv[iArg++]); 
  bool bEstWeights = atoi(argv[iArg++]);
  bool bEstDiffus = atoi(argv[iArg++]);
  int nEst = atoi(argv[iArg++]);
  double sparsity = atof(argv[iArg++]);
  int nRep = atoi(argv[iArg++]);
  
  StickEstimateOption options;
  options.maxIt = 50000;
  options.maxInnerIt = 5;
  options.init = 1;
  options.useManifold = true;
  options.step = 1e-9;
  options.kappaStep = 1e-12;
  options.kappa0Step = 1e-10;
  options.weightStep = 3e-10;
  options.tolerance = 1e-5;
  options.innerTolerance = 1e-6;
  options.isEstWeights = bEstWeights;
  options.isEstDiffusivities = bEstDiffus;
  options.useLineSearch = false;
  options.sparseFactor = sparsity;
  options.isPrintDebugInfo = false;

  mat gradMat = loadGradientOrientations("assets/brain/gradients.txt");
  double snr = 9.4377;
  int bVal = 2000;

  double s0;
  vec S = readVoxelFromBrainData(x-1, y-1, z-1, s0);
  vec diffusis = -log(S/s0)/bVal;
    
  FiberComposition fibComp;

  FiberComposition fullTensor;
  estimateTensor(fullTensor, S, gradMat, bVal, s0);
  double d1 = 2.0e-3;
  double d0 = 1.3e-3;

  double d2 = diffusis.max();
  double smax = sqrt(M_PI/(4*bVal*d2))*erf(sqrt(bVal*d2));
  double smin = exp(-bVal*d2);
  double smean = mean(S)/s0;
  double w0 = 1-(smean-smin)/(smax-smin);
  
  for (int i=0; i<nRep; ++i) {
    fibComp.nFibers = nEst;
    //initRandom(fibComp.fibDirs, nEst);
    fibComp.fibDirs = fullTensor.fibDirs.cols(0,nEst-1);
    fibComp.fibDiffs = d1*ones(nEst+1);
    fibComp.fibDiffs(nEst) = d0;
  
    if (method=="bas") {
      // fibComp.fibWeights = (1.0-w0)/nEst*ones(nEst+1);
      // fibComp.fibWeights(nEst) = w0;
      fibComp.fibWeights = 1.0/(nEst+1)*ones(nEst+1);
    } else
      fibComp.fibWeights = 1.0/nEst*ones(nEst);
    
    double l ;
    if (method=="bas")
      l = estimateBallAndSticks (fibComp, S, gradMat, bVal, s0, snr, nEst, options);
    else
      l = estimateModifiedBAS (fibComp, S, gradMat, bVal, s0, snr, nEst, options);

    fibComp.fibDirs.reshape(3*nEst,1);

    mat resMat = join_cols(fibComp.fibDirs,
  			   join_cols(fibComp.fibWeights, fibComp.fibDiffs));
    //cout <<x <<"\t" <<y <<"\t" <<z <<"\t" <<nEst <<"\t";
    for (int j=0; j<resMat.n_elem; ++j)
      cout <<resMat(j,0) <<"\t";
    cout <<l <<endl;
  }
  return 0;
}
