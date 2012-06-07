#include "HardiToolbox.h"
#include <fstream>
#include <math.h>
#include "besselFunction.h"

using namespace std;

mat HardiToolbox::loadGradientOrientations (const char *fileName)
{
  ifstream gradFile (fileName);
  if (!gradFile) {
    cerr <<"Can't open specified file!" <<endl;
    return mat ();
  }

  double tmpGradBuf [3000];

  char line [1280];
  int i=0;
  while (!gradFile.eof()) {
    gradFile.getline(line, 1280, '\n');
    stringstream ss (line);
    for (int j=0; j<3; ++j) {
      double dat;
      ss >>dat;
      tmpGradBuf[i*3+j] = dat;
    }
    ++i;
  }
  mat gradMat (tmpGradBuf, 3, i-1);
  return gradMat.t();    
}

vec HardiToolbox::simulateMultiTensor (int bVal, double s0, const mat &gradientOrientations, const vec &fibDirs,
				  const vec &weights, bool isAnisotropic)
{
  int nGrads = gradientOrientations.n_rows;
  int nFibers = fibDirs.n_elem;

  vec dwSignal (nGrads);
  dwSignal.fill(0);
  vec dwSignal1 (nGrads);
  dwSignal1.fill(0);

  mat D (3,3);
  D.eye();
  if (isAnisotropic) {
    D(0,0) = 1.7e-3;
    D(1,1) = 0.3e-3;
    D(2,2) = 0.3e-3;
  } else {
    D(0,0) = 7e-4;
    D(1,1) = 7e-4;
    D(2,2) = 7e-4;
  }

  for (int i=0; i<nGrads; ++i) {
    for (int j=0; j<nFibers; ++j) {
      double angle = fibDirs(j);
      mat R (3,3);
      R <<cos(angle) <<-sin(angle) <<0 <<endr
	<<sin(angle) <<cos(angle) <<0 <<endr
	<<0 <<0 <<1 <<endr;

      rowvec g = gradientOrientations.row(i);
      double w = weights.empty()?1.0:weights(j);
      dwSignal(i) += s0*w*exp(-bVal*dot(g*R*D*R.t(),g));
    }

  }
  return dwSignal;
}

vec HardiToolbox::simulateMultiTensor (int bVal, double s0, const mat &gradientOrientations, const mat &fibDirs,
				       const vec &weights, const mat &diffus)
{
  int nGrads = gradientOrientations.n_rows;
  int nFibers = fibDirs.n_cols;

  vec dwSignal = zeros(nGrads);

  for (int i=0; i<nGrads; ++i) {
    for (int j=0; j<nFibers; ++j) {
      vec eulerAngle (3);
      eulerAngle <<0 <<endr <<asin(fibDirs(2,j)) <<endr <<atan2(-fibDirs(1,j), fibDirs(0,j)) <<endr;
      mat R = eulerAngleToMatrix(eulerAngle);

      mat D (3,3);
      if (diffus.is_empty()) {
	D <<1.7e-3  <<0      <<0   <<endr
	  <<0       <<3e-4   <<0   <<endr
	  <<0       <<0      <<3e-4   <<endr;
      } else {
	D = diagmat(diffus.col(j));
      }

      rowvec g = gradientOrientations.row(i);
      double w = weights.empty()?1.0/nFibers:weights(j);
      dwSignal(i) += s0*w*exp(-bVal*dot(g*R*D*R.t(),g));
    }

  }
  return dwSignal;
}

mat HardiToolbox::simulateMultiTensorByComponent(int bVal, double s0, const mat &gradientOrientations, const mat &eulerAngles, const mat &diffusivities)
{
  mat D (3,3);

  int nGrads = gradientOrientations.n_rows;
  int nFibers = eulerAngles.n_cols;

  mat dwSignal = zeros(nGrads, nFibers);
  for (int i=0; i<nFibers; ++i)
  {
    mat R = eulerAngleToMatrix(eulerAngles.col(i));

    if (diffusivities.is_empty()) {
      D <<1.7e-3  <<0      <<0   <<endr
	<<0       <<3e-4   <<0   <<endr
	<<0       <<0      <<3e-4   <<endr;
    } else {
      D = diagmat(diffusivities.col(i));
    }

    mat tensor = R*D*R.t();

    for (int j=0; j<nGrads; ++j)
    {
      dwSignal(j, i) = s0*exp(-bVal*dot(gradientOrientations.row(j)*tensor, gradientOrientations.row(j)));
    }
  }

  return dwSignal;
}

vec HardiToolbox::addRicianNoise(const vec &cleanSignal, double sigma)
{
  int nGrads = cleanSignal.n_elem;
  vec noisySignal = sqrt(pow(cleanSignal+sigma*randn<vec>(nGrads), 2) + pow(sigma*randn<vec>(nGrads), 2));
  return noisySignal;
}

void HardiToolbox::loadSphereVecs (mat &triangleVecs, rowvec &triangleAreas, bool isAccurate)
{
  ifstream sphereFile;
  if (isAccurate) {
    sphereFile.open("assets/Sphere2562.txt");
  } else {
    sphereFile.open("assets/Sphere1280.txt");
  }

  int nVecs = isAccurate?2562:1280;
  triangleVecs = mat(3, nVecs);
  triangleAreas = rowvec(nVecs);

  char numCh [128];
  int i = 0;
  while (!sphereFile.eof()&&i<nVecs*4) {
    memset(numCh, 0, 128);
    if (i%4==3) {
      sphereFile.getline(numCh, 128, '\n');
      triangleAreas(i/4) = atof(numCh); //1 x nVecs
    } else {
      sphereFile.getline(numCh, 128, ',');
      triangleVecs(i%4, i/4) = atof(numCh); //3 x nVecs
    }

    ++i;
  }
}

void HardiToolbox::initRandom(mat &fibDirs, int nFibers)
{
  srand(time(0));
  vec r = randu(3);
  double phi = r(0)*2.0*M_PI;
  double theta = r(1)*M_PI;
  double psi = r(2)*2.0*M_PI;

  mat randMat(3,3);
  randMat <<cos(theta)*cos(psi)   <<cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)  <<sin(phi)*sin(psi)-cos(phi)*sin(theta)*cos(psi)  <<endr
          <<-cos(theta)*sin(psi)  <<cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi)  <<sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)  <<endr
          <<sin(theta)            <<-sin(phi)*cos(theta)                            <<cos(phi)*cos(theta)                             <<endr;
  fibDirs = randMat.cols(0, nFibers-1);
}

void HardiToolbox::initFixed(mat &fibDirs, int nFibers)
{
  mat e = eye(3,3);
  fibDirs = e.cols(0, nFibers-1);
}

vec HardiToolbox::signalError(const mat &fibDirs, const mat &kernel, int order, const mat &triangleVecs, const rowvec &triangleAreas, const vec &dwSignal)
{
  int nFibers = fibDirs.n_cols;
  mat triangleFODs = ones(1, nFibers)*pow(fibDirs.t()*triangleVecs, order); // 1 x nVecs
  vec simSignal = kernel*(triangleFODs.t()%triangleAreas.t()); // nGrads x 1
  vec e = dwSignal - simSignal;
  return e;
}

double HardiToolbox::energyFunc (const mat &fibDirs, const mat &kernel, int order, const mat &triangleVecs, const rowvec &triangleAreas, const vec &dwSignal, double sparse)
{
  vec sigErr = signalError(fibDirs, kernel, order, triangleVecs, triangleAreas, dwSignal);
  double sparsity = 0;
  for (int i=0; i<fibDirs.n_cols; ++i) {
    sparsity += pow(norm(fibDirs.col(i),2), order/2);
  }
  return accu(pow(sigErr, 2))-sparse*sparsity;
}

mat HardiToolbox::dirGradient (const mat &fibDirs, const mat &kernel, int order, const mat &triangleVecs, const rowvec &triangleAreas, const vec &dwSignal, double sparse, int updateFibIdx)
{
  int nFibers = fibDirs.n_cols;
  mat dDir = zeros(3, nFibers);
  vec sigErr = signalError(fibDirs, kernel, order, triangleVecs, triangleAreas, dwSignal); // nGrads x 1
  vec fibDir = fibDirs.col(updateFibIdx);
  mat dFODs = pow(fibDir.t()*triangleVecs, order-1); // 1 x nVecs
  mat dFODxArea = dFODs%triangleAreas; // 1 x nVecs
  dFODxArea = join_cols(dFODxArea, join_cols(dFODxArea, dFODxArea));
  mat c = -order*(dFODxArea%triangleVecs)*kernel.t(); // 3 x nGrads
  dDir.col(updateFibIdx) = c*sigErr;

  double n = norm(fibDir,2), d;
  double m = 0;
  for (int i=0; i<nFibers; ++i) {
    m += pow(norm(fibDirs.col(i), 2), d);
  }
  //dDir.col(updateFibIdx) -= sparse*order*pow(n, order-2)*(pow(n,order)-m/nFibers)*fibDir;
  dDir.col(updateFibIdx) += sparse*order/4.0*pow(n, order/2.0-2.0)*fibDir;
  
  return dDir;
}

mat HardiToolbox::dirGradientLM (const mat &fibDirs, const mat &kernel, int order, const mat &triangleVecs, const rowvec &triangleAreas, const vec &dwSignal, const vec &lambdas, int updateFibIdx)
{
  int nFibers = fibDirs.n_cols;
  mat dDir = zeros(3, nFibers);
  vec sigErr = signalError(fibDirs, kernel, order, triangleVecs, triangleAreas, dwSignal); // nGrads x 1
  vec fibDir = fibDirs.col(updateFibIdx);
  mat dFODs = pow(fibDir.t()*triangleVecs, order-1); // 1 x nVecs
  mat dFODxArea = dFODs%triangleAreas; // 1 x nVecs
  dFODxArea = join_cols(dFODxArea, join_cols(dFODxArea, dFODxArea));
  mat J = -order*kernel*(dFODxArea.t()%triangleVecs.t()); // nGrads x 3
  dDir.col(updateFibIdx) = solve(J.t()*J+lambdas(updateFibIdx)*diagmat(J.t()*J), J.t()*sigErr);
  return dDir;

}

void HardiToolbox::deconvolveFibers (FiberComposition &fibComp, const vec &dwSignal,
				    const mat &gradientDirs,
				    int bVal, int nFibers, const DeconvOption &options)
{
  mat triangleVecs;
  rowvec triangleAreas;

  loadSphereVecs(triangleVecs, triangleAreas, options.useAccurateIntegral);

  //initialize
  mat fibDirs (3, nFibers);

  switch (options.init) {
  default:
  case 0:
    //random initialization
    initRandom(fibDirs, nFibers);
    //fibDirs.print("random init");
    break;
  }
  mat kernel = exp(-options.delta*pow(gradientDirs*triangleVecs, 2)); //nGrads x nVecs
  vec lambdas (nFibers);
  lambdas.fill(options.lambda);
  double c = 5;

  //optimize
  double energy = energyFunc(fibDirs, kernel, options.order, triangleVecs, triangleAreas, dwSignal, options.sparse);
  double lastEnergy = 0;
  for (int it=0; it< options.maxIt; ++it)
  {
    lastEnergy = energy;
    //update
    for (int i=0; i<nFibers; ++i)
    {
      if (options.optMethod==0)
      {
        mat dDir = dirGradient(fibDirs, kernel, options.order, triangleVecs, triangleAreas, dwSignal, options.sparse, i);
        double step = options.step;
        if (options.useLineSearch)
        {
          step = 1;
          double a = 0.05;
          double b = 0.3;
          double n = dot(dDir.col(i), dDir.col(i));

          while (energyFunc(fibDirs-step*dDir, kernel, options.order, triangleVecs, triangleAreas, dwSignal, options.sparse)
		 >energyFunc(fibDirs, kernel, options.order, triangleVecs, triangleAreas, dwSignal, options.sparse) - a*n*step)
          {
            step *= b;
            if (step<1e-6)
            {
              break;
            }
          }
        }
        fibDirs -= step*dDir;
      } else
      {

        double e1 = energy;
        while (true)
        {
          mat dDir = dirGradientLM(fibDirs, kernel, options.order, triangleVecs, triangleAreas, dwSignal, lambdas, i);
          mat tmpFibDirs = fibDirs - dDir;
          energy = energyFunc(tmpFibDirs, kernel, options.order, triangleVecs, triangleAreas, dwSignal, 0);
          if (energy<=e1 || lambdas(i)>1e5)
          {
            fibDirs = tmpFibDirs;
            lambdas(i) /= c;
            break;
          } else
          {
            lambdas(i) *= c;
          }
        }
      }
    }

    //check for termination
    energy = energyFunc(fibDirs, kernel, options.order, triangleVecs, triangleAreas, dwSignal, options.sparse);
    cout <<"#" <<it <<":\te=" <<energy <<endl;
    if (fabs(lastEnergy-energy)<=options.tol)
    {
      break;
    }
  }

  //output
  //fibComp.nFibers = nFibers;
  fibComp.fibDirs = fibDirs;
  fibComp.fibWeights = zeros(nFibers);
  double sumWeights = 0;
  for (int i=0; i<nFibers; ++i)
  {
    double n = norm(fibDirs.col(i), 2);
    fibComp.fibDirs.col(i) /= n;
    n = std::pow(n, options.order);
    fibComp.fibWeights(i) = n;
    sumWeights += n;
  }
  fibComp.fibWeights /= sumWeights;

  for (int i=nFibers-1; i>=0; --i) {
    if (fibComp.fibWeights(i) < options.weightThres) {
      fibComp.fibWeights.shed_row(i);
      fibComp.fibDirs.shed_col(i);
    }
  }
  fibComp.nFibers = fibComp.fibWeights.n_elem;
  fibComp.fibDiffs = 1.7e-3*ones(fibComp.nFibers+1);
  fibComp.fibDiffs(fibComp.nFibers) = 0.3e-3;
}

vec HardiToolbox::directionDeviation(const mat &dirs, const mat &trueDirs)
{
  vec devs = zeros(trueDirs.n_cols);
  for (int i=0; i<(int)trueDirs.n_cols; ++i)
  {
    double minDev = M_PI/2;
    for (int j=0; j<(int)dirs.n_cols; ++j)
    {
      double c = fabs(dot(dirs.col(j), trueDirs.col(i)));
      c = c<0?0:(c>1?1:c);
      double dev = std::acos(c);
      minDev = min(minDev, dev);
    }
    devs(i) = minDev;
  }
  return devs;
}

vec HardiToolbox::fiberDeviation(const FiberComposition &fibComp, const FiberComposition &trueFibComp)
{
  double sumDev = 0;
  double sumWeightErr = 0;
  const mat &dirs = fibComp.fibDirs;
  const mat &trueDirs = trueFibComp.fibDirs;

  for (int i=0; i<(int)trueDirs.n_cols; ++i)
  {
    double minDev = M_PI/2;
    double weightErr = 0;
    for (int j=0; j<(int)dirs.n_cols; ++j)
    {
      double dev = std::acos(dot(dirs.col(j), trueDirs.col(i)));
      dev = dev>M_PI/2?M_PI-dev:dev;
      if (dev<minDev) {
	minDev = dev;
	weightErr = fibComp.fibWeights(j) - trueFibComp.fibWeights(i);
      }
    }
    sumDev += minDev;
    sumWeightErr += fabs(weightErr);
  }

  vec retDev (2);
  retDev <<sumDev/trueDirs.n_cols*180/M_PI <<endr <<sumWeightErr/trueDirs.n_cols;
  return retDev;
}

vec HardiToolbox::simulateStick(int bVal, double s0, double d, const mat &gradientOrientations, const vec &fibDirs, const vec &weights)
{
  int nFibers = fibDirs.n_elem;
  int nGrads = gradientOrientations.n_rows;
  vec dwSignal (nGrads);
  dwSignal.fill(0);

  for (int j=0; j<nFibers; ++j)
  {
    vec dir (3);
    dir <<cos(fibDirs(j)) <<endr <<sin(fibDirs(j)) <<endr <<0 <<endr;
    dir.print("synthesized dir");
    double w = weights.is_empty()?1.0 : weights(j);
    for (int i=0; i<nGrads; ++i)
    {
      dwSignal(i) += w*s0*exp(-bVal*d*pow(dot(gradientOrientations.row(i), dir),2));
    }
  }
  return dwSignal;
}

vec HardiToolbox::simulateStick(int bVal, double s0, const vec &d, const mat &gradientOrientations, const vec &fibDirs)
{
  int nFibers = fibDirs.n_elem;
  int nGrads = gradientOrientations.n_rows;
  vec dwSignal (nGrads);
  dwSignal.fill(0);

  for (int j=0; j<nFibers; ++j)
  {
    vec dir (3);
    dir <<cos(fibDirs(j)) <<endr <<sin(fibDirs(j)) <<endr <<0 <<endr;
    for (int i=0; i<nGrads; ++i)
    {
      dwSignal(i) += s0*exp(-bVal*d(j)*pow(dot(gradientOrientations.row(i), dir),2));
    }
  }
  return dwSignal;
}

vec HardiToolbox::simulateBallAndStick(int bVal, double s0, const vec &d, const vec &w, const mat &gradientOrientations, const vec &fibDirs)
{
  int nFibers = fibDirs.n_elem;
  int nGrads = gradientOrientations.n_rows;
  vec dwSignal (nGrads);
  dwSignal.fill(0);

  for (int j=0; j<nFibers; ++j)
  {
    vec dir (3);
    dir <<cos(fibDirs(j)) <<endr <<sin(fibDirs(j)) <<endr <<0 <<endr;

    double w1 = w.is_empty()?1.0/(nFibers+1):w(j);
    for (int i=0; i<nGrads; ++i)
    {
      dwSignal(i) += w1*s0*exp(-bVal*d(j)*pow(dot(gradientOrientations.row(i), dir),2));
    }
  }

  for (int i=0; i<nGrads;++i)
  {
    dwSignal(i) += (w.is_empty()?1.0:w(nFibers))*s0*exp(-bVal*d(nFibers));
  }

  return dwSignal;
}

vec HardiToolbox::simulateModBAS(int bVal, double s0, const vec &d, const vec &w, const mat &gradientOrientations, const vec &fibDirs)
{
  int nFibers = fibDirs.n_elem;
  int nGrads = gradientOrientations.n_rows;
  vec dwSignal (nGrads);
  dwSignal.fill(0);


  for (int i=0; i<nGrads; ++i)
  {
    for (int j=0; j<nFibers; ++j)
    {
      vec dir (3);
      dir <<cos(fibDirs(j)) <<endr <<sin(fibDirs(j)) <<endr <<0 <<endr;
      double weight = w.is_empty()?1.0/nFibers:w(j);
      dwSignal(i) += weight*s0*exp(-bVal*d(nFibers))*exp(-bVal*d(j)*pow(dot(gradientOrientations.row(i), dir),2));
    }
  }
  return dwSignal;
}

mat HardiToolbox::simulateStickByComponent(int bVal, double s0, double d, const mat &gradientOrientations, const mat &fibDirs, const vec &weights)
{
  int nFibers = fibDirs.n_cols;
  int nGrads = gradientOrientations.n_rows;

  mat dwSignal (nGrads, nFibers);
  dwSignal.fill(0);
  for (int i=0; i<nGrads; ++i)
  {
    for (int j=0; j<nFibers; ++j)
    {
      double w = weights.is_empty()?1.0:weights(j);
      dwSignal(i,j) = w*s0*exp(-bVal*d*pow(dot(gradientOrientations.row(i), fibDirs.col(j)),2));
    }
  }
  return dwSignal;
}

mat HardiToolbox::simulateStickByComponent(int bVal, double s0, const vec &d, const mat &gradientOrientations, const mat &fibDirs)
{
  int nFibers = fibDirs.n_cols;
  int nGrads = gradientOrientations.n_rows;

  mat dwSignal (nGrads, nFibers);
  dwSignal.fill(0);
  for (int i=0; i<nGrads; ++i)
  {
    for (int j=0; j<nFibers; ++j)
    {
      dwSignal(i,j) = s0*exp(-bVal*d(j)*pow(dot(gradientOrientations.row(i), fibDirs.col(j)),2));
    }
  }
  return dwSignal;
}

mat HardiToolbox::simulateBallAndStickByComponent(int bVal, double s0, const vec &d, const vec &w, const mat &gradientOrientations, const mat &fibDirs)
{
  int nFibers = fibDirs.n_cols;
  int nGrads = gradientOrientations.n_rows;

  mat dwSignal (nGrads, nFibers+1);
  dwSignal.fill(0);
  for (int i=0; i<nGrads; ++i)
  {
    for (int j=0; j<nFibers; ++j)
    {
      double w1 = w.is_empty()?1.0/(nFibers+1):w(j);
      dwSignal(i,j) = w1*s0*exp(-bVal*d(j)*pow(dot(gradientOrientations.row(i), fibDirs.col(j)),2));
    }
  }

  for (int i=0; i<nGrads;++i)
  {
    dwSignal(i, nFibers) = (w.is_empty()?1.0/(nFibers+1):w(nFibers))*s0*exp(-bVal*d(nFibers));
  }

  return dwSignal;
}

mat HardiToolbox::simulateModBASByComponent(int bVal, double s0, const vec &d, const vec &w,  const mat &gradientOrientations, const mat &fibDirs)
{
  int nFibers = fibDirs.n_cols;
  int nGrads = gradientOrientations.n_rows;

  mat dwSignal (nGrads, nFibers);
  dwSignal.fill(0);
  for (int i=0; i<nGrads; ++i)
  {
    for (int j=0; j<nFibers; ++j)
    {
      double weight = w.is_empty()?1.0/nFibers:w(j);
      dwSignal(i,j) = weight*s0*exp(-bVal*d(nFibers))*exp(-bVal*d(j)*pow(dot(gradientOrientations.row(i), fibDirs.col(j)),2));
    }
  }

  return dwSignal;
}

void HardiToolbox::estimateSticks(FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations,
                                  int bVal, double s0, double d, double snr, int nFibers, const StickEstimateOption &options)
{
  //initialize
  int nGrads = gradientOrientations.n_rows;
  mat fibDirs(3, nFibers);
  initRandom(fibDirs, nFibers);
  fibDirs.print("random init");
  double sigma = s0/snr;

  //optimize
  mat estimatedSignal = simulateStickByComponent(bVal, s0, d, gradientOrientations, fibDirs);
  vec A (nGrads);
  bool waitForInput = false;
  double e = 1e10;
  for (int it=0; it<options.maxIt; ++it)
  {
    mat oldDirs = fibDirs;
    mat oldSignal = estimatedSignal;
    vec oldSignalSum = sum(oldSignal, 1);
    for (int i=0; i<nGrads; ++i)
    {
      double x = oldSignalSum(i) * dwSignal(i)/sigma/sigma;
      A(i) = bessiRatio(x);
    }

    for (int iit=0; iit<options.maxInnerIt; ++iit)
    {

      estimatedSignal = simulateStickByComponent(bVal, s0, d, gradientOrientations, fibDirs);
      mat dDirs = zeros(3, nFibers);
      for (int j=0; j<nFibers; ++j)
      {
        for (int i=0; i<nGrads; ++i)
        {
          vec g = gradientOrientations.row(i).t();
          dDirs.col(j) += 4*(estimatedSignal(i,j)-oldSignal(i,j)+oldSignalSum(i)/nFibers-dwSignal(i)/nFibers*A(i))*estimatedSignal(i,j)*(bVal*d*dot(g,fibDirs.col(j)))*g;
        }
      }

      double step = options.step;
//      double a = 0.2;
//      double b = 0.8;
//      double n = accu(pow(dDirs, 2));
//      while (stickLogLikelihood(fibDirs-step*dDirs, gradientOrientations, bVal, A, dwSignal) > stickLogLikelihood(fibDirs, gradientOrientations, bVal, A, dwSignal)-a*n*step)
//      {
//        step *= b;
//        if (step<1e-3)
//        {
//          break;
//        }
//      }

      fibDirs += step*dDirs;

      if (sqrt(accu(pow(dDirs,2)))<options.innerTolerance)
      {
        break;
      }
    }

    e = sqrt(accu(pow(fibDirs-oldDirs, 2)));
    cout <<"it #" <<it <<": d=" <<e <<endl;
    if (e<options.tolerance)
    {
      break;
    }

    if(waitForInput)
    {
      std::string s;
      std::getline(std::cin, s);
      if (strcmp(s.c_str(), "c")==0)
      {
        waitForInput = false;
      }
    }
  }

  //output
  fibComp.fibDiffs = vec(nFibers);
  fibComp.fibDirs = mat(3, nFibers);
  for (int i=0; i<nFibers; ++i)
  {
    double n = norm(fibDirs.col(i), 2);
    fibComp.fibDirs.col(i) = fibDirs.col(i)/n;
    fibComp.fibDiffs(i) = n*n*d;
  }
  fibComp.fibWeights = 1.0/nFibers*ones(nFibers);
  fibComp.nFibers = nFibers;
}

void HardiToolbox::estimateSticksDiffs(FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations,
                                  int bVal, double s0, double snr, int nFibers, const StickEstimateOption &options)
{
  //initialize
  int nGrads = gradientOrientations.n_rows;
  mat fibDirs(3, nFibers);
  initRandom(fibDirs, nFibers);

  fibDirs.print("random init");
  double sigma = s0/snr;

  vec kappas = 1.0e-3*ones(nFibers);

  //optimize
  mat estimatedSignal = simulateStickByComponent(bVal, s0, kappas, gradientOrientations, fibDirs);
  vec A (nGrads);
  bool waitForInput = false;
  double e = stickLogLikelihood(estimatedSignal, estimatedSignal, dwSignal, sigma);
  for (int it=0; it<options.maxIt; ++it)
  {
    mat oldDirs = fibDirs;
    mat oldSignal = estimatedSignal;
    vec oldSignalSum = sum(oldSignal, 1);
    for (int i=0; i<nGrads; ++i)
    {
      double x = oldSignalSum(i) * dwSignal(i)/sigma/sigma;
      A(i) = bessiRatio(x);
    }

    for (int iit=0; iit<options.maxInnerIt; ++iit)
    {

      estimatedSignal = simulateStickByComponent(bVal, s0, kappas, gradientOrientations, fibDirs);
      mat dDirs = zeros(3, nFibers);
      vec dKappas = zeros(nFibers);
      for (int j=0; j<nFibers; ++j)
      {
        vec fibDir = fibDirs.col(j);
        for (int i=0; i<nGrads; ++i)
        {
          vec g = gradientOrientations.row(i).t();
          double b = 2*(estimatedSignal(i,j)-oldSignal(i,j)+oldSignalSum(i)/nFibers-dwSignal(i)/nFibers*A(i));

          vec logUG = sphereLog(fibDir, g);
          double theta = norm(logUG, 2);
          if (theta>1e-12)
          {
            dDirs.col(j) += b*estimatedSignal(i,j)*bVal*kappas(j)*sin(2*theta)/theta*logUG;
          }

          dKappas(j) -= b*estimatedSignal(i,j)/s0*pow(dot(g, fibDir), 2);
        }
      }

      for (int i=0; i<nFibers; ++i)
      {
        fibDirs.col(i) = sphereExp(fibDirs.col(i), options.step*dDirs.col(i));
      }

      kappas -= options.kappaStep*dKappas;

      if (sqrt(accu(pow(dDirs,2)))<options.innerTolerance)
      {
        break;
      }
    }

    double lastE = e;
    e = stickLogLikelihood(estimatedSignal, oldSignal, dwSignal, sigma);
    cout <<"it #" <<it <<": e=" <<e <<endl;
    if (fabs(e-lastE)<options.tolerance)
    {
      break;
    }

    if(waitForInput)
    {
      std::string s;
      std::getline(std::cin, s);
      if (strcmp(s.c_str(), "c")==0)
      {
        waitForInput = false;
      }
    }
  }

  //output
  fibComp.fibDirs = fibDirs;
  for (int i=0; i<nFibers; ++i)
  {
    fibComp.fibDirs.col(i) /= norm(fibComp.fibDirs.col(i), 2);
  }
  fibComp.fibDiffs = kappas;
  fibComp.fibWeights = 1.0/nFibers*ones(nFibers);
  fibComp.nFibers = nFibers;
}


void HardiToolbox::estimateBallAndSticks(FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations,
                                  int bVal, double s0, double snr, int nFibers, const StickEstimateOption &options)
{
  //initialize
  int nGrads = gradientOrientations.n_rows;
  mat fibDirs(3, nFibers);
  if (options.init ==0 ) {
    initRandom(fibDirs, nFibers);
  } else {
    fibDirs = fibComp.fibDirs;
  }

  double sigma = s0/snr;

  vec kappas = 2.3e-3*ones(nFibers+1);
  kappas(nFibers) = 1.6e-3;

  vec weights = 1.0/(nFibers+1)*ones(nFibers+1);

  //optimize
  mat estimatedSignal = simulateBallAndStickByComponent(bVal, s0, kappas, weights, gradientOrientations, fibDirs);
  vec A (nGrads);
  bool waitForInput = false;
  double e = stickLogLikelihood(estimatedSignal, estimatedSignal, dwSignal, sigma);
  for (int it=0; it<options.maxIt; ++it)
  {
    mat oldDirs = fibDirs;
    mat oldSignal = estimatedSignal;
    vec oldSignalSum = sum(oldSignal, 1);
    for (int i=0; i<nGrads; ++i)
    {
      double x = oldSignalSum(i) * dwSignal(i)/sigma/sigma;
      A(i) = bessiRatio(x);
    }

    for (int iit=0; iit<options.maxInnerIt; ++iit)
    {

      estimatedSignal = simulateBallAndStickByComponent(bVal, s0, kappas, weights, gradientOrientations, fibDirs);
      mat dDirs = zeros(3, nFibers);
      vec dKappas = zeros(nFibers+1);
      vec dWeights = zeros(nFibers+1);

      for (int i=0; i<nGrads; ++i)
      {
        for (int j=0; j<nFibers+1; ++j)
        {
          vec g = gradientOrientations.row(i).t();
          double b = 2*(estimatedSignal(i,j)-oldSignal(i,j)+oldSignalSum(i)/nFibers-dwSignal(i)/nFibers*A(i));

          if (j<nFibers) /**sticks**/
          {
            vec fibDir = fibDirs.col(j);
            vec logUG = sphereLog(fibDir, g);
            double theta = norm(logUG, 2);
            if (theta>1e-12)
            {
              dDirs.col(j) += b*estimatedSignal(i,j)*bVal*kappas(j)*sin(2*theta)/theta*logUG;
            }

            dKappas(j) += b*estimatedSignal(i,j)/s0*pow(dot(g, fibDir), 2);
	    dWeights(j) += s0*exp(-bVal*kappas(j)*pow(dot(g, fibDir), 2));
          } else /**ball**/
          {
            dKappas(j) += b*estimatedSignal(i,j)/s0;
	    dWeights(j) += s0*exp(-bVal*kappas(nFibers));
          }
        }
      }

      for (int i=0; i<nFibers; ++i)
      {
        fibDirs.col(i) = sphereExp(fibDirs.col(i), options.step*dDirs.col(i));
      }
      
      if(options.isEstDiffusivities) {
	kappas += options.kappaStep*dKappas;
      }

      if (options.isEstWeights) {
	cout <<dWeights <<endl;
	weights += options.weightStep*dWeights;
	weights -= (sum(weights)-1)/(nFibers+1)*ones(nFibers+1);
	for (int i=0; i<nFibers+1; ++i) {
	  weights(i) = weights(i)<0?0:weights(i);
	}
	weights /= sum(weights);
      }

      if (sqrt(accu(pow(dDirs,2)))<options.innerTolerance)
      {
        break;
      }
    }

    double lastE = e;
    e = stickLogLikelihood(estimatedSignal, oldSignal, dwSignal, sigma);
    cout <<"it #" <<it <<": e=" <<e <<endl;
    if (fabs(e-lastE)<options.tolerance)
    {
      break;
    }

    if(waitForInput)
    {
      std::string s;
      std::getline(std::cin, s);
      if (strcmp(s.c_str(), "c")==0)
      {
        waitForInput = false;
      }
    }
  }

  //output
  fibComp.fibDirs = fibDirs;
  // for (int i=0; i<nFibers; ++i)
  // {
  //   fibComp.fibDirs.col(i) /= norm(fibComp.fibDirs.col(i), 2);
  // }
  fibComp.fibDiffs = kappas;
  fibComp.fibWeights = weights;
  fibComp.nFibers = nFibers;
}

double HardiToolbox::estimateModifiedBAS(FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations,
				       int bVal, double s0, double snr, int nFibers, const StickEstimateOption &options)
{  
  /**initialize**/
  int nGrads = gradientOrientations.n_rows;
  mat fibDirs(3, nFibers);
  vec kappas;
  vec weights;

  if (options.init==0 || fibComp.fibDirs.is_empty()) {
    initRandom(fibDirs, nFibers);
    kappas = 2.0e-3*ones(nFibers+1);
    kappas(nFibers) = 1.0e-3;
    weights = 1.0/nFibers*ones(nFibers);
  } else {
    initRandom(fibDirs, nFibers);
    for (int i=0; i<fibComp.fibDirs.n_cols; ++i) {
      fibDirs.col(i) = fibComp.fibDirs.col(i);
    }

    kappas = fibComp.fibDiffs(0)*ones(nFibers+1);
    kappas(nFibers) = fibComp.fibDiffs(fibComp.nFibers);
    weights = 1.0/nFibers*ones(nFibers);
  }

  fibDirs.print("init dirs");

  double sigma = s0/snr;

  /**optimize**/
  mat estimatedSignal = simulateModBASByComponent(bVal, s0, kappas, weights, gradientOrientations, fibDirs);
  vec A (nGrads);
  bool waitForInput = false;
  // double e = stickLogLikelihood(estimatedSignal, estimatedSignal, dwSignal, sigma);
  vec S = sum(estimatedSignal,1);
  double e = ricianLikelihood(S, dwSignal, sigma);
  double lastE = e;
  bool bDecreasing = true;
  int it = 0;
  for (it=0; it<options.maxIt; ++it)
  {

    /**optimize ball**/
    if(options.isEstDiffusivities) {
      for (int iit=0; iit<options.maxInnerIt; ++iit) {
	double dKappa0 = 0;
	estimatedSignal = simulateModBASByComponent(bVal, s0, kappas, weights, gradientOrientations, fibDirs);
	vec estimatedSignalSum = sum(estimatedSignal, 1);

	for (int i=0; i<nGrads; ++i) {
	  double x = dwSignal(i) * estimatedSignalSum(i)/sigma/sigma;
	  dKappa0 += (estimatedSignalSum(i) - dwSignal(i)*bessiRatio(x))*estimatedSignalSum(i)/pow(sigma, 2);
	}

	kappas(nFibers) += options.kappa0Step*dKappa0;

	if (fabs(dKappa0)<options.innerTolerance) {
	  break;
	}
      }
    }

    /**optimize sticks**/
    mat oldDirs = fibDirs;
    mat oldSignal = estimatedSignal;
    vec oldSignalSum = sum(oldSignal, 1);
    for (int i=0; i<nGrads; ++i)
    {
      double x = oldSignalSum(i) * dwSignal(i)/sigma/sigma;
      A(i) = bessiRatio(x);
    }

    for (int iit=0; iit<options.maxInnerIt; ++iit)
    {
      estimatedSignal = simulateModBASByComponent(bVal, s0, kappas, weights, gradientOrientations, fibDirs);
      mat dDirs = zeros(3, nFibers);
      vec dKappas = zeros(nFibers+1);
      vec dWeights = zeros(nFibers);

      for (int i=0; i<nGrads; ++i)
      {
        for (int j=0; j<nFibers; ++j)
        {
          vec g = gradientOrientations.row(i).t();
          double b = 2*(estimatedSignal(i,j)-oldSignal(i,j)+oldSignalSum(i)/nFibers-dwSignal(i)/nFibers*A(i));

          vec fibDir = fibDirs.col(j);
          vec logUG = sphereLog(fibDir, g);
          double theta = norm(logUG, 2);
          if (theta>1e-12)
          {
            dDirs.col(j) += b*estimatedSignal(i,j)*bVal*kappas(j)*sin(2*theta)/theta*logUG;
          }

	  if (options.isEstDiffusivities) {
	    dKappas(j) += b*estimatedSignal(i,j)*pow(dot(gradientOrientations.row(i), fibDirs.col(j)), 2);
	  }

	  if (options.isEstWeights) {
	    dWeights(j) -= bVal*s0*exp(-bVal*kappas(nFibers))*exp(-bVal*kappas(j)*pow(dot(g, fibDir),2));
	  }
        }
      }

      for (int i=0; i<nFibers; ++i)
      {
	double step = options.step;
	mat tmpFibDirs = fibDirs;
	if (options.useLineSearch) {
	  step = 1.0;
	  double a = 0.5;
	  double b = 0.3;
	  mat s1 = simulateModBASByComponent(bVal, s0, kappas, weights, gradientOrientations, fibDirs);
	  tmpFibDirs.col(i) = sphereExp(fibDirs.col(i), step*dDirs.col(i));
	  mat s2 = simulateModBASByComponent(bVal, s0, kappas, weights, gradientOrientations, tmpFibDirs);
	  double e1 = stickLogLikelihood(s1, oldSignal, dwSignal, sigma);
	  double e2 = stickLogLikelihood(s2, oldSignal, dwSignal, sigma);
	  while (step >= 1e-10 && e2 > e1 - b*step*dot(dDirs.col(i), dDirs.col(i))) {
	    step *= a;
	    s1 = simulateModBASByComponent(bVal, s0, kappas, weights, gradientOrientations, fibDirs);
	    tmpFibDirs.col(i) = sphereExp(fibDirs.col(i), step*dDirs.col(i));
	    s2 = simulateModBASByComponent(bVal, s0, kappas, weights, gradientOrientations, tmpFibDirs);
	    e1 = stickLogLikelihood(s1, oldSignal, dwSignal, sigma);
	    e2 = stickLogLikelihood(s2, oldSignal, dwSignal, sigma);
	  }
	}
        fibDirs.col(i) = sphereExp(fibDirs.col(i), step*dDirs.col(i));
      }

      if (options.isEstDiffusivities) {
	kappas += options.kappaStep*dKappas;
      }

      if (options.isEstWeights) {
	weights += options.weightStep*dWeights;
	weights -= (sum(weights)-1)/nFibers*ones(nFibers);
	for (int i=0; i<nFibers; ++i) {
	  weights(i) = weights(i)<0?0:weights(i);
	}
	weights /= sum(weights);
      }

      if (sqrt(accu(pow(dDirs,2)))<options.innerTolerance)
      {
	cout <<"internal stopped. norm of kappas: " <<norm(dKappas,2) <<endl;
        break;
      }
    }

    lastE = e;
    estimatedSignal = simulateModBASByComponent(bVal, s0, kappas, weights, gradientOrientations, fibDirs);
    S = sum(estimatedSignal,1);
    //e = stickLogLikelihood(estimatedSignal, oldSignal, dwSignal, sigma);
    e = ricianLikelihood(S, dwSignal, sigma);
    if (bDecreasing && e>lastE) {
      cout <<"it #" <<it <<": e=" <<e <<endl;
      // cout <<"energy increased " <<e-lastE <<endl;
      // fibDirs.print("--local minima dirs");
      // kappas.print("--local minima diffus");
      // weights.print("--local minima weights");
      bDecreasing = false;
    }

    if (!bDecreasing && e<lastE) {
      cout <<"it #" <<it <<": e=" <<e <<endl;
      // cout <<"energy decreased" <<endl;
      // fibDirs.print("++local maxima dirs");
      // kappas.print("++local maxima diffus");
      // weights.print("++local maxima weights");
      bDecreasing = true;
    }

    // if (fabs(lastE-e)<options.tolerance)
    if (e-lastE<options.tolerance)
    {
      break;
    }

    if(waitForInput)
    {
      std::string s;
      std::getline(std::cin, s);
      if (strcmp(s.c_str(), "c")==0)
      {
        waitForInput = false;
      }
    }
  }
  cout <<"it #" <<it <<": e=" <<e <<". decreased " <<lastE-e <<endl;

  //output
  fibComp.fibDirs = fibDirs;
  // for (int i=0; i<nFibers; ++i) {
  //   fibComp.fibDirs.col(i) /= norm(fibComp.fibDirs.col(i), 2);
  // }
  fibComp.fibDiffs = kappas;
  fibComp.fibWeights = weights;
  fibComp.nFibers = nFibers;

  return e;
}

void HardiToolbox::estimateModifiedBASByStick(FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations,
				       int bVal, double s0, double snr, int nFibers, const StickEstimateOption &options)
{  
  /**initialize**/
  int nGrads = gradientOrientations.n_rows;
  mat fibDirs(3, nFibers);
  vec kappas;
  vec weights;

  if (options.init==0 || fibComp.fibDirs.is_empty()) {
    initRandom(fibDirs, nFibers);
    kappas = 2.0e-3*ones(nFibers+1);
    kappas(nFibers) = 1.6e-3;
    weights = 1.0/nFibers*ones(nFibers);
  } else {
    initRandom(fibDirs, nFibers);
    for (int i=0; i<fibComp.fibDirs.n_cols; ++i) {
      fibDirs.col(i) = fibComp.fibDirs.col(i);
    }

    kappas = fibComp.fibDiffs(0)*ones(nFibers+1);
    kappas(nFibers) = fibComp.fibDiffs(fibComp.nFibers);
    weights = 1.0/nFibers*ones(nFibers);
  }

  fibDirs.print("init dirs");

  double sigma = s0/snr;

  /**optimize**/
  mat estimatedSignal = simulateModBASByComponent(bVal, s0, kappas, weights, gradientOrientations, fibDirs);
  vec A (nGrads);
  bool waitForInput = false;
  double e = stickLogLikelihood(estimatedSignal, estimatedSignal, dwSignal, sigma);
  double lastE = e;
  bool bDecreasing = true;
  int it = 0;
  for (it=0; it<options.maxIt; ++it)
  {

    /**optimize ball**/
    if(options.isEstDiffusivities) {
      for (int iit=0; iit<options.maxInnerIt; ++iit) {
	double dKappa0 = 0;
	estimatedSignal = simulateModBASByComponent(bVal, s0, kappas, weights, gradientOrientations, fibDirs);
	vec estimatedSignalSum = sum(estimatedSignal, 1);

	for (int i=0; i<nGrads; ++i) {
	  double x = dwSignal(i) * estimatedSignalSum(i)/sigma/sigma;
	  dKappa0 += (estimatedSignalSum(i) - dwSignal(i)*bessiRatio(x))*estimatedSignalSum(i)/pow(sigma, 2);
	}

	kappas(nFibers) += options.kappa0Step*dKappa0;
	cout <<dKappa0 <<endl;
	if (fabs(dKappa0)<options.innerTolerance) {
	  break;
	}
      }
    }

    /**optimize sticks**/
    mat oldDirs = fibDirs;
    mat oldSignal = estimatedSignal;
    vec oldSignalSum = sum(oldSignal, 1);
    for (int i=0; i<nGrads; ++i)
    {
      double x = oldSignalSum(i) * dwSignal(i)/sigma/sigma;
      A(i) = bessiRatio(x);
    }

    for (int j=0; j<nFibers; ++j) {
      for (int iit=0; iit<options.maxInnerIt; ++iit) {
	estimatedSignal = simulateModBASByComponent(bVal, s0, kappas, weights, gradientOrientations, fibDirs);
	mat dDirs = zeros(3, nFibers);
	vec dKappas = zeros(nFibers+1);
	vec dWeights = zeros(nFibers);

	for (int i=0; i<nGrads; ++i) {
          vec g = gradientOrientations.row(i).t();
          double b = 2*(estimatedSignal(i,j)-oldSignal(i,j)+oldSignalSum(i)/nFibers-dwSignal(i)/nFibers*A(i));

          vec fibDir = fibDirs.col(j);
          vec logUG = sphereLog(fibDir, g);
          double theta = norm(logUG, 2);
          if (theta>1e-12) {
            dDirs.col(j) += b*estimatedSignal(i,j)*bVal*kappas(j)*sin(2*theta)/theta*logUG;
          }

	  if (options.isEstDiffusivities) {
	    dKappas(j) += b*estimatedSignal(i,j)*pow(dot(gradientOrientations.row(i), fibDirs.col(j)), 2);
	    // dKappas(0) += b*estimatedSignal(i,j)*pow(dot(gradientOrientations.row(i), fibDirs.col(j)), 2)/100;
	  }

	  if (options.isEstWeights) {
	    dWeights(j) -= bVal*s0*exp(-bVal*kappas(nFibers))*exp(-bVal*kappas(j)*pow(dot(g, fibDir),2));
	  }
        }

	for (int i=0; i<nFibers; ++i) {
	  fibDirs.col(i) = sphereExp(fibDirs.col(i), options.step*dDirs.col(i));
	}

	// cout <<"dkappa: " <<dKappas(0) <<endl;
	if (options.isEstDiffusivities) {
	  kappas += options.kappaStep*dKappas;
	  // for (int i=0; i<nFibers; ++i) {
	  //   kappas(i) += options.kappaStep*dKappas(0);
	  // }
	}

	if (options.isEstWeights) {
	  weights += options.weightStep*dWeights;
	  weights -= (sum(weights)-1)/nFibers*ones(nFibers);
	  for (int i=0; i<nFibers; ++i) {
	    weights(i) = weights(i)<0?0:weights(i);
	  }
	  weights /= sum(weights);
      }

	if (sqrt(accu(pow(dDirs,2)))<options.innerTolerance) {
	  cout <<"internal stopped. norm of kappas: " <<norm(dKappas,2) <<endl;
	  break;
	}
      }
    }

    lastE = e;
    e = stickLogLikelihood(estimatedSignal, oldSignal, dwSignal, sigma);
    if (bDecreasing && e>lastE) {
      cout <<"it #" <<it <<": e=" <<e <<endl;
      cout <<"energy increased " <<e-lastE <<endl;
      fibDirs.print("local minima");
      bDecreasing = false;
    }

    if (!bDecreasing && e<lastE) {
      cout <<"it #" <<it <<": e=" <<e <<endl;
      cout <<"energy decreased" <<endl;
      bDecreasing = true;
    }
    if (e-lastE<options.tolerance)
      // if (fabs(lastE-e)<options.tolerance)
    {
      break;
    }

    if(waitForInput)
    {
      std::string s;
      std::getline(std::cin, s);
      if (strcmp(s.c_str(), "c")==0)
      {
        waitForInput = false;
      }
    }
  }
  cout <<"it #" <<it <<": e=" <<e <<endl;

  //output
  fibComp.fibDirs = fibDirs;
  // for (int i=0; i<nFibers; ++i) {
  //   fibComp.fibDirs.col(i) /= norm(fibComp.fibDirs.col(i), 2);
  // }
  fibComp.fibDiffs = kappas;
  fibComp.fibWeights = weights;
  fibComp.nFibers = nFibers;
}
 
void HardiToolbox::estimateSticksWeights(FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations,
                                  int bVal, double s0, double d, double snr, int nFibers, const StickEstimateOption &options)
{
  /**initialize**/
  int nGrads = gradientOrientations.n_rows;
  mat fibDirs(3, nFibers);
  initRandom(fibDirs, nFibers);

//  fibDirs <<1 <<0 <<0 <<endr <<0 <<1 <<0 <<endr;
//  fibDirs = fibDirs.t();

  fibDirs.print("random init");
  vec weights = 1.0/nFibers*ones(nFibers);
//  weights <<0.6 <<endr <<0.4 <<endr;
  double sigma = s0/snr;

  /**optimize**/
  mat estimatedSignal = simulateStickByComponent(bVal, s0, d, gradientOrientations, fibDirs, weights);
  vec A (nGrads);
  bool waitForInput = false;
  double e = stickLogLikelihood(estimatedSignal, estimatedSignal, dwSignal, sigma);
  for (int it=0; it<options.maxIt; ++it)
  {
    mat oldDirs = fibDirs;
    mat oldSignal = estimatedSignal;
    vec oldSignalSum = sum(estimatedSignal,1);
    for (int i=0; i<nGrads; ++i)
    {
      double x = oldSignalSum(i) * dwSignal(i)/sigma/sigma;
      A(i) = bessiRatio(x);
    }

    for (int iit=0; iit<options.maxInnerIt; ++iit)
    {

      estimatedSignal = simulateStickByComponent(bVal, s0, d, gradientOrientations, fibDirs, weights);
      mat dDirs = zeros(3, nFibers);
      vec dWeights = zeros(nFibers);
      for (int j=0; j<nFibers; ++j)
      {
        vec fibDir = fibDirs.col(j);
        for (int i=0; i<nGrads; ++i)
        {
          vec g = gradientOrientations.row(i).t();
          double b = 2*(estimatedSignal(i,j)-oldSignal(i,j)+oldSignalSum(i)/nFibers-dwSignal(i)/nFibers*A(i));
          if (!options.useManifold)
          {
            dDirs.col(j) += 2*b*estimatedSignal(i,j)*(bVal*d*dot(g,fibDirs.col(j)))*g;
          } else
          {
            /**on sphere**/
            vec logUG = sphereLog(fibDir, g);
            double theta = norm(logUG, 2);
//            logUG.print("logUG");
//            cout <<"theta=" <<theta <<endl;
            if (theta>1e-12)
            {
              dDirs.col(j) += b*estimatedSignal(i,j)*bVal*d*sin(2*theta)/theta*logUG;
            }

          }

          dWeights(j) += b*s0*exp(-bVal*d*pow(dot(g, fibDirs.col(j)),2));
        }
      }

//      dDirs.print("dDirs");
      if (!options.useManifold)
      {
        fibDirs += options.step*dDirs;
      } else
      {
        for (int i=0; i<nFibers; ++i)
        {
          fibDirs.col(i) = sphereExp(fibDirs.col(i), options.step*dDirs.col(i));
        }
      }


      weights -= options.step*dWeights;
      weights -= (sum(weights)-1)/nFibers*ones(nFibers);
      for (int i=0; i<nFibers; ++i)
      {
        weights(i) = weights(i)<0?0:weights(i);
      }
      weights /= sum(weights);

      if (sqrt(accu(pow(dDirs,2)))<options.innerTolerance)
      {
        break;
      }
    }

    double lastE = e;
    e = stickLogLikelihood(estimatedSignal, oldSignal, dwSignal, sigma);
    cout <<"it #" <<it <<": e=" <<e <<endl;
    if (fabs(e-lastE)<options.tolerance)
    {
      break;
    }

    if(waitForInput)
    {
      std::string s;
      std::getline(std::cin, s);
      if (strcmp(s.c_str(), "c")==0)
      {
        waitForInput = false;
      }
    }
  }

  //output
  for (int i=0; i<nFibers; ++i)
  {
    fibDirs.col(i) /= norm(fibDirs.col(i), 2);
  }
  fibComp.fibDirs = fibDirs;
  fibComp.fibWeights = weights;
  fibComp.nFibers = nFibers;
}

double HardiToolbox::stickLogLikelihood(const mat &fibDirs, const mat &gradientOrientations,
                                        int bVal, double s0, double d, const vec &A, const vec &dwSignal)
{
  mat estimatedSignal = simulateStickByComponent(bVal, s0, d, gradientOrientations, fibDirs);
  vec estimatedSignalSum = sum(estimatedSignal, 1);
  int nFibers = fibDirs.n_cols;
  int nGrads = gradientOrientations.n_rows;
  double logLike = 0;
  for (int i=0; i<nGrads; ++i)
  {
    for (int j=0; j<nFibers; ++j)
    {
      logLike += estimatedSignal(i,j)*estimatedSignal(i,j) - 2.0/nFibers*dwSignal(i)*A(i)*estimatedSignal(i,j);

    }
  }
  return logLike;
}

double HardiToolbox::stickLogLikelihood(const mat &estimatedSignal, const mat &oldSignal, const vec &dwSignal, double sigma)
{
  int nFibers = estimatedSignal.n_cols;
  int nGrads = estimatedSignal.n_rows;
  vec oldSignalSum = sum(oldSignal,1);
  double logLike = 0;
  for (int i=0; i<nGrads; ++i)
  {
    for (int j=0; j<nFibers; ++j)
    {
      logLike += pow(estimatedSignal(i,j),2)-2*estimatedSignal(i,j)*(oldSignal(i,j)-oldSignalSum(i)/nFibers+dwSignal(i)/nFibers*bessiRatio(oldSignalSum(i)*dwSignal(i)/sigma/sigma));
    }
  }
  return logLike;
}

mat HardiToolbox::eulerAngleToMatrix(const vec &eulerAngle)
{
  mat R (3,3);
  double phi = eulerAngle(0);
  double theta = eulerAngle(1);
  double psi = eulerAngle(2);
  R <<cos(theta)*cos(psi)   <<cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)  <<sin(phi)*sin(psi)-cos(phi)*sin(theta)*cos(psi)  <<endr
    <<-cos(theta)*sin(psi)  <<cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi)  <<sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)  <<endr
    <<sin(theta)            <<-sin(phi)*cos(theta)                            <<cos(phi)*cos(theta)                             <<endr;

  return R;
}

void HardiToolbox::estimateMultiTensor (FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations, int bVal, double s0, int nFibers, const MultiTensorOption &options)
{
  //initialize
  int nGrads = gradientOrientations.n_rows;
  srand(time(0));
  mat eulerAngles = randu(3, nFibers);
  eulerAngles.row(0) *= 2*M_PI;
  eulerAngles.row(1) *= M_PI;
  eulerAngles.row(2) *= 2*M_PI;

  mat diffusivities (3, nFibers);
  for (int i=0; i<nFibers; ++i) {
    diffusivities(0,i) = 1.7e-3;
    diffusivities(1,i) = 0.3e-3;
    diffusivities(2,i) = 0.3e-3;
  }

  if (options.init==1) {
    /**initialize with full tensor**/
    FiberComposition fullTensor;
    estimateTensor (fullTensor, dwSignal, gradientOrientations, bVal, s0);
    mat fibDirs = fullTensor.fibDirs;
    for (int i=0; i<nFibers; ++i) {
      fibDirs.swap_cols(0,i);
      eulerAngles(1,i) = asin(fibDirs(2,0)) + M_PI/2;
      eulerAngles(2,i) = atan2(-fibDirs(1,0), fibDirs(0,0)) + M_PI;
      eulerAngles(0,i) = atan2(-fibDirs(2,1), fibDirs(2,2)) + M_PI;
    }

    // for (int i=0; i<min(2,nFibers); ++i) {
    //   diffusivities(0,i) = 2*fullTensor.fibDiffs(i) - fullTensor.fibDiffs(2);
    //   diffusivities(1,i) = fullTensor.fibDiffs(2);
    //   diffusivities(2,i) = fullTensor.fibDiffs(2);
    // }
  }

  /////////////////test//////////////////////
  //  eulerAngles <<0 <<0 <<0 <<endr <<0 <<0 <<M_PI/2 <<endr;
  //  eulerAngles = eulerAngles.t();

  cube rotMatrices (3,3,nFibers);

  vec etas = ones(nFibers);
  vec weights = 1.0/nFibers*ones(nFibers);

  bool waitForInput = false;
  //optimize
  double e = 1e10;
  for (int it = 0; it<options.maxIt; ++it)
  {
    weights = exp(etas)/accu(exp(etas));
    for (int i=0; i<nFibers; ++i)
    {
      rotMatrices.slice(i) = eulerAngleToMatrix(eulerAngles.col(i));
      vec fibDir = rotMatrices.slice(i).col(0);
    }

    mat estimatedSignal = simulateMultiTensorByComponent(bVal, s0, gradientOrientations, eulerAngles, diffusivities);
    //estimatedSignal.print("estimated signal");
    vec estimatedSignalSum = zeros(nGrads);
    for (int i=0; i<nGrads; ++i)
    {
      for (int j=0; j<nFibers; ++j)
      {
        estimatedSignalSum(i) += weights(j)*estimatedSignal(i,j);
      }
    }

    //estimatedSignal.print("estimated signal by component");
    //estimatedSignalSum.print("estimated signal");
    //dwSignal.print("true signal");

    double lastE = e;
    e = norm(estimatedSignalSum-dwSignal,2);
    // cout <<"it #" <<it <<":\t" <<e <<endl;
    if (lastE - e<options.tolerance)
    {
      break;
    }


    if(waitForInput)
    {
      std::string s;
      std::getline(std::cin, s);
      if (strcmp(s.c_str(), "c")==0)
      {
        waitForInput = false;
      }
    }

    mat dEulerAngles = zeros<mat>(3, nFibers);
    vec dEtas = zeros(nFibers);
    double sqSumExpEta = pow(accu(exp(etas)),2);

    for (int i=0; i<nFibers; ++i)
    {
      //directions
      double phi = eulerAngles(0,i);
      double theta = eulerAngles(1,i);
      double psi = eulerAngles(2,i);
      mat R = rotMatrices.slice(i);
      mat D = diagmat(diffusivities.col(i));

      for (int j=0; j<3; ++j)
      {
        mat dR (3,3);
        switch(j)
        {
        case 0: //d_phi
          dR  <<0   <<-sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)   <<cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)  <<endr
              <<0   <<-sin(phi)*cos(psi)-cos(phi)*sin(theta)*sin(psi)   <<cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi)  <<endr
              <<0   <<-cos(phi)*cos(theta)                              <<-sin(phi)*cos(theta)                            <<endr;
          break;
        case 1: //d_theta
          dR  <<-sin(theta)*cos(psi)  <<sin(phi)*cos(theta)*cos(psi)  <<-cos(phi)*cos(theta)*cos(psi)   <<endr
              <<sin(theta)*sin(psi)   <<-sin(phi)*cos(theta)*sin(psi) <<cos(phi)*cos(theta)*sin(psi)    <<endr
              <<cos(theta)            <<sin(phi)*sin(theta)           <<-cos(phi)*sin(theta)            <<endr;
          break;
        case 2: //d_psi
          dR  <<-cos(theta)*sin(psi)  <<cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi)    <<sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)    <<endr
              <<-cos(theta)*cos(psi)  <<-cos(phi)*sin(psi)-sin(phi)*sin(theta)*cos(psi)   <<-sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)   <<endr
              <<0                     <<0                                                 <<0                                                 <<endr;
          break;
        }

        mat B = dR*D*R.t()+R*D*dR.t();

        for (int k=0; k<nGrads; ++k)
        {
          rowvec g = gradientOrientations.row(k);
          dEulerAngles(j,i) += -(estimatedSignalSum(k)-dwSignal(k))*weights(i)*estimatedSignal(k,i)*dot(g*B, g)*bVal;
        }
      }

      //weights
      for (int k=0; k<nGrads; ++k)
      {
        double b=0;
        for (int j=0; j<nFibers; ++j)
        {
          b += i==j? 0: (estimatedSignal(k,i)-estimatedSignal(k,j))*exp(etas(j));
        }
        dEtas(i) += (estimatedSignalSum(k)-dwSignal(k))*b;
      }
      dEtas(i) *= exp(etas(i))/sqSumExpEta;
    }

//    dEulerAngles.print("dEulerAngles");
//    dEtas.print("dEtas");
    eulerAngles -= options.step*dEulerAngles;
    etas -= options.step*dEtas;
  }

  //output
  fibComp.fibDirs = mat(3, nFibers);
  for (int i=0; i<nFibers; ++i)
  {
    double theta = eulerAngles(1,i);
    double psi = eulerAngles(2,i);
    fibComp.fibDirs(0,i) = cos(theta)*cos(psi);
    fibComp.fibDirs(1,i) = -cos(theta)*sin(psi);
    fibComp.fibDirs(2,i) = sin(theta);
  }
  fibComp.fibWeights = weights;
  fibComp.nFibers = nFibers;
}

vec HardiToolbox::sphereLog(const vec &p, const vec &q)
{
  vec q1 = q;
  if (dot(p,q1)<0)
  {
    q1 *= -1;
  }
  double cosTheta = dot(p,q1);
  if (fabs(1-cosTheta)<1e-12)
  {
    return zeros(3);
  }
  vec l = (q1-cosTheta*p)/sqrt(1.0-cosTheta*cosTheta)*acos(cosTheta);
  return l;
}

vec HardiToolbox::sphereExp(const vec &p, const vec &v)
{
  double theta = norm(v, 2);
  if (theta<1e-12)
  {
    return p;
  }
  vec e = cos(theta)*p+sin(theta)*v/theta;
  e /= norm(e,2);
  return e;
}

void HardiToolbox::estimateTensor (FiberComposition &fibComp, const vec &dwSignal, const mat &gradientOrientations, int bVal, double s0)
{
  int nGrads = gradientOrientations.n_rows;
  mat G (nGrads, 6);
  for (int i=0; i<nGrads; ++i) {
    rowvec g = gradientOrientations.row(i);
    G(i,0) = g(0)*g(0);
    G(i,1) = g(1)*g(1);
    G(i,2) = g(2)*g(2);
    G(i,3) = 2*g(0)*g(1);
    G(i,4) = 2*g(1)*g(2);
    G(i,5) = 2*g(2)*g(0);
  }

  vec S = -log(dwSignal/s0)/bVal;
  vec d = solve(G, S);
  mat D (3,3);
  D <<d(0) <<d(3) <<d(5) <<endr
    <<d(3) <<d(1) <<d(4) <<endr
    <<d(5) <<d(4) <<d(2) <<endr;

  vec eigVal;
  mat eigVec;

  // use standard algorithm by default
  eig_sym(eigVal, eigVec, D);

  fibComp.nFibers = 1;
  fibComp.fibDirs.resize(3, 3);
  fibComp.fibDiffs.resize(3);
  for (int i=0; i<3; ++i) {
    fibComp.fibDirs.col(i) = eigVec.col(2-i);
    fibComp.fibDiffs(i) = eigVal(2-i);
  }
  fibComp.fibWeights = ones(1);
}

double HardiToolbox::ricianLikelihood (const vec &estimatedSignal, const vec &realSignal, double sigma)
{
  double l = 0;
  int nGrads = estimatedSignal.n_rows;
  double s = sigma*sigma;
  for (int i=0; i<nGrads; ++i) {
    l += log(realSignal(i)) - 2*log(s) - (realSignal(i)*realSignal(i)+estimatedSignal(i)*estimatedSignal(i))/(2*s) + log(bessi0(realSignal(i)*estimatedSignal(i)/s));
  }
  return l;
}
