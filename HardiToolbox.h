#ifndef HARDITOOLBOX_H
#define HARDITOOLBOX_H

#include <armadillo>
using namespace arma;



namespace HardiToolbox
{

  mat       loadGradientOrientations (const char * fileName);
  vec       addRicianNoise (const vec &cleanSignal, double sigma);
  vec       directionDeviation (const mat &dirs, const mat &trueDirs);

  /**************************************/
  /******direct fiber deconvolution******/
  struct    DeconvOption
  {
    int     order;
    int     delta;
    int     maxIt;
    int     init;
    int     optMethod;
    bool    useAccurateIntegral;
    bool    useLineSearch;
    double  step;
    double  tol;
    double  lambda;
    double  sparse;
    double  weightThres;
  };

  struct    FiberComposition
  {
    mat     fibDirs;
    vec     fibWeights;
    vec     fibDiffs;
    int     nFibers;
  };

  vec       fiberDeviation (const FiberComposition &fibComp, const FiberComposition &trueFibComp);

  vec       simulateMultiTensor (int bVal, double s0, const mat &gradientOrientations, const vec &fibDirs,
  			                            const vec &weights=vec(), bool isAnisotropic=true);
  vec       simulateMultiTensor (int bVal, double s0, const mat &gradientOrientations, const mat &fibDirs,
			       const vec &weights=vec(), bool isAnisotropic=true);
  mat       simulateMultiTensorByComponent (int bVal, double s0, const mat &gradientOrientations, const mat &eulerAngles, const mat &diffusivities = mat());

  void      loadSphereVecs (mat &triangleVecs, rowvec &triangleAreas, bool isAccurate = false);
  void      initRandom (mat &fibDirs, int nFibers);
  void      initFixed (mat &fibDirs, int nFibers);
  vec       signalError (const mat &fibDirs, const mat &kernel, int order, const mat &triangleVecs, const rowvec &triangleAreas, const vec &dwSignal);
  double    energyFunc (const mat &fibDirs, const mat &kernel, int order, const mat &triangleVecs, const rowvec &triangleAreas, const vec &dwSignal, double sparse);
  mat       dirGradient (const mat &fibDirs, const mat &kernel, int order, const mat &triangleVecs, const rowvec &triangleAreas, const vec &dwSignal, double sparse, int updateFibIdx);
  mat       dirGradientLM (const mat &fibDirs, const mat &kernel, int order, const mat &triangleVecs, const rowvec &triangleAreas, const vec &dwSignal, const vec &lambdas, int updateFibIdx);
  void      deconvolveFibers (FiberComposition &fibComp, const vec &dwSignal,
                              const mat &gradientOrientations,
			      int bVal, int nFibers, const DeconvOption &options);

  void      estimateTensor (FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrienations, int bVal, double s0);

  /**************************************************/
  /******finite mixture of Gaussian [Tuch2002]*******/
  struct    MultiTensorOption
  {
    int     maxIt;
    int     init;
    double  step;
    double  tolerance;
  };

  mat       eulerAngleToMatrix (const vec &eulerAngle);
  void      estimateMultiTensor (FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations, int bVal, double s0, int nFibers, const MultiTensorOption &options);

  /**************************************/
  /******ball and stick estimation*******/

  vec       simulateStick (int bVal, double s0, double d, const mat &gradientOrientations, const vec &fibDirs, const vec &weights=vec());
  mat       simulateStickByComponent (int bVal, double s0, double d, const mat &gradientOrientations, const mat &fibDirs, const vec &weights=vec());

  vec       simulateStick (int bVal, double s0, const vec &diffusivities, const mat &gradientOrientations, const vec &fibDirs);
  mat       simulateStickByComponent (int bVal, double s0, const vec &diffusivities, const mat &gradientOrientations, const mat &fibDirs);

  vec       simulateBallAndStick (int bVal, double s0, const vec &diffusivities, const vec &weights, const mat &gradientOrientations, const vec &fibDirs);
  mat       simulateBallAndStickByComponent (int bVal, double s0, const vec &diffusivities, const vec &weights, const mat &gradientOrientations, const mat &fibDirs);

  vec       simulateModBAS (int bVal, double s0, const vec &diffusivities, const vec &weights, const mat &gradientOrientations, const vec &fibDirs);
  mat       simulateModBASByComponent (int bVal, double s0, const vec &diffusivities, const vec &weights, const mat &gradientOrientations, const mat &fibDirs);


  struct    StickEstimateOption
  {
    int     maxIt;
    int     maxInnerIt;
    int     init;
    bool    useManifold;
    double  step;
    double  kappaStep;
    double  kappa0Step;
    double  weightStep;
    double  tolerance;
    double  innerTolerance;
    bool    isEstDiffusivities;
    bool    isEstWeights;
    bool    useLineSearch;
  };

  void      estimateSticks (FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations,
                            int bVal, double s0, double d, double snr, int nFibers, const StickEstimateOption &options);
  void      estimateSticksWeights (FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations,
                            int bVal, double s0, double d, double snr, int nFibers, const StickEstimateOption &options);
  void      estimateSticksDiffs (FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations,
                            int bVal, double s0, double snr, int nFibers, const StickEstimateOption &options);

  void      estimateBallAndSticks (FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations,
                            int bVal, double s0, double snr, int nFibers, const StickEstimateOption &options);
  void      estimateModifiedBAS (FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations,
                            int bVal, double s0, double snr, int nFibers, const StickEstimateOption &options);
  void      estimateModifiedBASByStick (FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations,
                            int bVal, double s0, double snr, int nFibers, const StickEstimateOption &options);

  double    stickLogLikelihood (const mat &fibDirs, const mat &gradientOrientations, int bVal, double s0, double d, const vec &A, const vec &dwSignal);

  double    stickLogLikelihood (const mat &estimatedSignal, const mat &oldSignal, const vec &dwSignal, double sigma);

  vec       sphereLog (const vec &p, const vec &q);
  vec       sphereExp (const vec &p, const vec &v);
}

#endif
