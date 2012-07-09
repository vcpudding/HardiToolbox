#include "testProc.h"
#include "erf.h"

#include <fstream>

using namespace std;
using namespace HardiToolbox;

void testBASRicianEM()
{

  struct timeval tStart, tEnd;

  int bVal = 2000;
  double s0 = 1500;
  double snr = 2000;
  double angle = M_PI/6;
  int nFibers = 2;
  int nEstFibers = 2;
  int nTrials = 10;
  mat fibDirs;
  vec weights;
  mat diffMat;
  vec diffVec;
  double d1 = 1.7e-3;
  double d0 = 3e-4;

  if (nFibers == 1)
  {
    double w = 0.2;
    fibDirs <<1 <<endr <<0 <<endr <<0 <<endr;
    weights <<w  <<endr <<1-w <<endr;
    diffMat <<d1 <<endr <<d0 <<endr <<d0 <<endr;
    diffVec <<d1 <<endr <<d0 <<endr;
  } else
  {
    double w1 = 1.0/3.0, w2 = 1.0/3.0;
    fibDirs <<1 <<cos(angle) <<endr <<0 <<sin(angle) <<endr <<0 <<0 <<endr;
    weights <<w1 <<endr <<w2  <<endr <<1-w1-w2 <<endr;
    diffMat <<d1 <<d1 <<endr <<d0 <<d0 <<endr <<d0 <<d0 <<endr;
    diffVec <<d1 <<endr <<d1 <<endr <<d0 <<endr;
  }

  for (int i=0; i<nFibers; ++i) {
  }

  mat R(3,3);
  initRandom(R, 3);
  //fibDirs = R*fibDirs;

  StickEstimateOption stickOptions;
  stickOptions.maxIt = 25000;
  stickOptions.maxInnerIt = 5;
  stickOptions.init = 1;
  stickOptions.useManifold = true;
  stickOptions.step = 1e-9;
  stickOptions.kappaStep = 1e-12;
  stickOptions.kappa0Step = 1e-10;
  stickOptions.weightStep = 1e-10;
  stickOptions.tolerance = 1e-6;
  stickOptions.innerTolerance = 1e-6;
  stickOptions.isEstWeights = false;
  stickOptions.isEstDiffusivities = true;
  stickOptions.useLineSearch = false;
  stickOptions.sparseFactor = 1e5;
  stickOptions.isPrintDebugInfo = true;

  mat gradientOrientations = loadGradientOrientations("gradients.txt");
  vec S = simulateBallAndStick(bVal, s0, gradientOrientations, fibDirs, weights, diffVec);
  vec S1 = addRicianNoise(S, s0/snr);
  FiberComposition fibComp;
  UtilStopWatch::tic();
  fibComp.nFibers = nEstFibers;
  initRandom(fibComp.fibDirs, nEstFibers);
  fibComp.fibDiffs = 2.0e-3*ones(nEstFibers+1);
  fibComp.fibDiffs(nEstFibers) = 0.3e-4;
  fibComp.fibWeights = 1.0/(nEstFibers+1)*ones(nEstFibers+1);
  estimateBallAndSticks(fibComp, S1, gradientOrientations, bVal, s0, snr, nEstFibers, stickOptions);
  cout <<"time: " <<UtilStopWatch::toc() <<"ms" <<endl;
  fibComp.fibDiffs.print("estimated diffusivities");
  fibComp.fibDirs.print("estimated directions");
  fibComp.fibWeights.print("estimated weights");
  FiberComposition trueFibComp;
  trueFibComp.nFibers = nFibers;
  trueFibComp.fibDirs = fibDirs;
  trueFibComp.fibDiffs = diffVec;
  trueFibComp.fibWeights = weights;
  
  fibDirs.print("true directions");
  vec dirDev, weightDev, diffusDev;
  double meanErr = fiberDeviation(fibComp, trueFibComp, dirDev, weightDev, diffusDev);
  dirDev.print("direction error");
  weightDev.print("weight error");
  diffusDev.print("diffusion error");
}

void testModBASRicianEM()
{
  struct timeval tStart, tEnd;

  int bVal = 3000;
  double s0 = 1500;
  double snr = 20;
  double angle = M_PI/4;
  int nFibers = 2;
  int nEst = 3;
  int nTrials = 10;
  mat fibDirs;
  vec weights;
  mat diffMat;
  vec diffVec;
  double d1 = 2.0e-3;
  double d0 = 3e-4;

  if (nFibers == 1)
  {
    fibDirs <<1 <<endr <<0 <<endr <<0 <<endr;
    weights <<1  <<endr;
    diffMat <<d1 <<endr <<d0 <<endr <<d0 <<endr;
    diffVec <<d1-d0 <<endr <<d0 <<endr;
  } else
  {
    fibDirs <<1 <<cos(angle) <<endr <<0 <<sin(angle) <<endr <<0 <<0 <<endr;
    weights <<0.3 <<endr <<0.7 <<endr;
    diffMat <<d1 <<d1 <<endr <<d0 <<d0 <<endr <<d0 <<d0 <<endr;
    diffVec <<d1-d0 <<endr <<d1-d0 <<endr <<d0 <<endr;
  }

  mat R(3,3);
  initRandom(R, 3);
  //fibDirs = R*fibDirs;

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
  stickOptions.isEstWeights = true;
  stickOptions.isEstDiffusivities = false;
  stickOptions.useLineSearch = false;
  stickOptions.sparseFactor = 1e5;

  mat gradientOrientations = loadGradientOrientations("gradients.txt");
  vec S = simulateMultiTensor(bVal, s0, gradientOrientations, fibDirs, weights, diffMat);
  vec S1 = addRicianNoise(S, s0/snr);
  FiberComposition fibComp;
  fibComp.nFibers = nEst;
  initRandom(fibComp.fibDirs, nEst);
  fibComp.fibDiffs = 1.7e-3*ones(nEst+1);
  fibComp.fibDiffs(nEst) = 3e-4;
  fibComp.fibWeights = 1.0/nEst*ones(nEst);

  UtilStopWatch::tic();
  estimateModifiedBAS(fibComp, S1, gradientOrientations, bVal, s0, snr, nEst, stickOptions);
  cout <<"time: " <<UtilStopWatch::toc() <<"ms" <<endl;
  fibComp.fibDiffs.print("estimated diffusivities");
  fibComp.fibDirs.print("estimated directions");
  fibComp.fibWeights.print("estimated weights");
  FiberComposition trueFibComp;
  trueFibComp.nFibers = nFibers;
  trueFibComp.fibDirs = fibDirs;
  trueFibComp.fibDiffs = diffVec;
  trueFibComp.fibWeights = weights;
  
  fibDirs.print("true directions");
  vec dirDev, weightDev, diffusDev;
  double meanErr = fiberDeviation(fibComp, trueFibComp, dirDev, weightDev, diffusDev);
  dirDev.print("direction error");
  weightDev.print("weight error");
  diffusDev.print("diffusion error");
}

ImageType::Pointer readImage (const char *fileName)
{
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fileName);
  reader->Update();

  ImageType::Pointer image = reader->GetOutput();
  ImageType::RegionType region = image->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();
 
  std::cout <<"dimensions: \t" <<size <<endl;
  return image;
}

vec getVoxelValue (ImageType::Pointer img, int x, int y, int z)
{
  ImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();
  int nDim = size[3];
  vec dat (nDim);
  for (int i=0; i<nDim; ++i) {
    ImageType::IndexType currentIndex;
    currentIndex[0] = x;
    currentIndex[1] = y;
    currentIndex[2] = z;
    currentIndex[3] = i;
			   
    dat(i) = img->GetPixel(currentIndex);
  }

  return dat;    
}

void trimSignalAndGrads (vec &S, mat &gradMat, double &s0)
{
  int nGrads = S.n_elem;
  s0 = 0;
  int nRep = 0;
  for (int i=nGrads-1; i>=0; --i) {
    if (norm(gradMat.row(i), 2)<1e-3) {
      gradMat.shed_row(i);
      s0 += S(i);
      S.shed_row(i);
      ++nRep;
    }
  }
  s0 /= nRep;
}



void estimate3DVolume (const char *dwiFileName, const char *maskFileName, const char *gradientFileName, 
		       const char *saveFileName, int bVal, double snr, int nFibers, void *opt)
{
  ImageType::Pointer pDWImage = readImage(dwiFileName);
  // ImageType::Pointer pMaskImage = readImage(maskFileName);

  ImageType::SizeType imgSize = pDWImage->GetLargestPossibleRegion().GetSize();

  mat oriGradMat = loadGradientOrientations(gradientFileName);
  mat estimateResults = zeros(imgSize[0]*imgSize[1]*imgSize[2], nFibers*3+3);
  DeconvOption *options = static_cast<DeconvOption*>(opt);
  FiberComposition fibComp;
  int nEstimated = 0;

  for (int x=0; x<imgSize[0]; ++x) {
    for (int y=0; y<imgSize[1]; ++y) {
      for (int z=0; z<1/*imgSize[2]*/; ++z) {
	cout <<"(" <<x <<", " <<y <<", " <<z <<"): ";
	ImageType::IndexType currentIndex;
	currentIndex[0] = x;
	currentIndex[1] = y;
	currentIndex[2] = z;
	currentIndex[3] = 0;
	double m = pDWImage->GetPixel(currentIndex);
	if (m<300.0) {
	  cout <<"skipped" <<endl;
	  continue;
	}

	cout <<"estimating" <<endl;
	mat gradMat = oriGradMat;
	double s0 = 0;
	vec S = getVoxelValue(pDWImage, x, y, z);
	trimSignalAndGrads(S, gradMat, s0);
	deconvolveFibers(fibComp, S, gradMat, bVal, nFibers, *options);

	estimateResults(nEstimated,0) = x;
	estimateResults(nEstimated,1) = y;
	estimateResults(nEstimated,2) = z;
	for (int i=0; i<nFibers; ++i) {
	  for (int j=0; j<3; ++j) {
	    estimateResults(nEstimated,3+3*i+j) = fibComp.fibDirs(j, i);
	  }
	}
	++nEstimated;
      }
    }
  }

  estimateResults.resize(nEstimated, nFibers*3+3);
  estimateResults.save(saveFileName, arma_ascii);
}




void testSparsityTerm () 
{
  int nRep = 20;
  int bVal = 3000;
	    
  DeconvOption options;
  options.order= 20; // polynomial order
  options.delta= 3.4; //29.2; % estimated from the data, b-value multiplied by the dominant diffusivity
  options.lambda=1e-5; //initial damping parameter value
  options.tol= 1e-6; // convergence criteria
  options.maxIt= 2000;
  options.step = 5e-4;
  options.useAccurateIntegral = false;
  options.init = 0;
  options.useLineSearch = true;
  options.optMethod = 0;
  options.weightThres = 0.1;
  options.sparse = 0.02;

  for (int snr = 10; snr<=30; snr+=10) {
    for (int nFibers = 2; nFibers<=2; ++nFibers) {
      for (int useSparse=0; useSparse<=1; ++useSparse) {
	if (!useSparse) {
	  options.sparse = 0;
	} else {
	  options.sparse = 0.02;
	  /* switch (snr) { */
	  /* case 10: */
	  /*   options.sparse = 0.1; */
	  /*   break; */
	  /* default: */
	  /*   options.sparse = 0.05; */
	  /* } */

	}

	char fileName [1024];
	sprintf(fileName, "rec_snr=%d_M=%d_s=%0.2f.txt\0", snr, nFibers, options.sparse);

	ofstream recFile (fileName);

	for (int iRep = 0; iRep<nRep; ++iRep) {

	  FiberComposition trueFibComp;
	  trueFibComp.nFibers = nFibers;
	  trueFibComp.fibDirs = mat(3, nFibers);
	  trueFibComp.fibWeights = vec(nFibers);
	  vec simFibDirs (nFibers);
	
	  srand(time(NULL));
	  double sepAngle = (30+rand()%61)*M_PI/180.0;
	  double weight = (rand()*1.0/RAND_MAX)*0.2+0.3;
	
	  if (nFibers==1) {
	    trueFibComp.fibDirs <<cos(sepAngle) <<endr <<sin(sepAngle) <<endr <<0 <<endr;
	    trueFibComp.fibWeights <<1.0 <<endr;
	    simFibDirs <<sepAngle <<endr;
	  } else {
	    trueFibComp.fibDirs <<1 <<cos(sepAngle) <<endr <<0 <<sin(sepAngle) <<endr <<0 <<0 <<endr;
	    trueFibComp.fibWeights <<weight <<endr <<1.0-weight <<endr;
	    simFibDirs <<0 <<endr <<sepAngle <<endr;
	  }
	  
	  mat gradMat = loadGradientOrientations("assets/Gradient64.txt");
	  vec S = simulateMultiTensor(bVal, 1, gradMat, simFibDirs, trueFibComp.fibWeights);
	  vec S1 = addRicianNoise(S, 1.0/snr);
	  FiberComposition fibComp;
	  UtilStopWatch::tic();
	  deconvolveFibers(fibComp, S1, gradMat, bVal, 3, options);
	  long elapsedTime = UtilStopWatch::toc();
	  vec err = fiberDeviation(fibComp, trueFibComp);

	  recFile <<sepAngle*180/M_PI <<"," <<weight <<"," <<fibComp.nFibers-nFibers <<"," <<err(0) <<"," <<err(1) <<"," <<elapsedTime <<endl;
	}
	recFile.close();
      }
    }
  }
}

void testTomsAlgorithm (const StickEstimateOption &options)
{
  ifstream inFile ("Matlab/log.txt");
  ofstream outFile ("tomsalgorithm.txt");
  int nTrials = 10;
  char line [1280];

  mat gradMat0 = loadGradientOrientations("diffusion_directions.txt");
  ImageType::Pointer img = readImage("dwi-b1500.nii");
  double snr = 17.6;
  int bVal = 1500;

  while (!inFile.eof()) {
    inFile.getline(line, 1280, '\n');
    if (inFile.gcount()<1) {
      break;
    }

    stringstream ss (line);
    int x, y, nFibers;
    ss >>x >>y >>nFibers;
    cout <<"(" <<x <<"," <<y <<"): " <<nFibers <<endl;
    mat trueDirs (3, nFibers);
    for (int i=0; i<nFibers; ++i) {
      ss >>trueDirs(0,i) >>trueDirs(1,i) >>trueDirs(2,i);
    }

    mat gradMat = gradMat0;
    vec S = getVoxelValue(img, x-1, y-1, 0);
    double s0;
    trimSignalAndGrads(S, gradMat, s0);
    FiberComposition fibComp;
    vec devs = zeros(nTrials);
    for (int iTrial=0; iTrial<nTrials; ++iTrial) {
      estimateModifiedBAS (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);
      devs(iTrial) = mean(directionDeviation(fibComp.fibDirs, trueDirs));
      cout <<"Trial #" <<iTrial <<": " <<devs(iTrial) <<endl;
    }
    
    char logStr [1280];
    sprintf(logStr, "%d\t%d\t%d\t%0.4f\t%0.4f\n", x, y, nFibers, mean(devs), stddev(devs));
    outFile <<logStr;
  }

  inFile.close();
  outFile.close();
}

void testTomsAlgorithmOneVoxel (int x, int y, int nFibers, const StickEstimateOption &options)
{
  mat gradMat = loadGradientOrientations("diffusion_directions.txt");
  ImageType::Pointer img = readImage("dwi-b1500.nii");
  double snr = 17.6;
  int bVal = 1500;

  vec S = getVoxelValue(img, x-1, y-1, 0);
  double s0;
  trimSignalAndGrads(S, gradMat, s0);
  FiberComposition fibComp;

  ifstream inFile ("Matlab/log.txt");
  mat trueDirs;
  bool bTrueDirsFound = false;
  char line [1280];

  while (!inFile.eof()) {
    inFile.getline(line, 1280, '\n');
    if (inFile.gcount()<1) {
      break;
    }

    stringstream ss (line);
    int x0, y0, nFibers0;
    ss >>x0 >>y0 >>nFibers0;
    if (x0 == x && y0 == y) {
      trueDirs.resize(3, nFibers0);
      for (int i=0; i<nFibers0; ++i) {
	ss >>trueDirs(0,i) >>trueDirs(1,i) >>trueDirs(2,i);
      }
      bTrueDirsFound = true;
      break;      
    }
  }

  FiberComposition fullTensor;
  estimateTensor(fullTensor, S, gradMat, bVal, s0);

  fibComp.nFibers = nFibers;
  initRandom(fibComp.fibDirs, nFibers);
  fibComp.fibDiffs = 2.7e-3*ones(nFibers+1);
  fibComp.fibDiffs(nFibers) = 1.2e-3;
  
  // fibComp.fibDirs = fullTensor.fibDirs.cols(0, nFibers-1);
  // fibComp.fibDiffs = fullTensor.fibDiffs(0)*ones(nFibers+1);
  // fibComp.fibDiffs(nFibers) = (fullTensor.fibDiffs(1)+fullTensor.fibDiffs(2))/2;
  fibComp.fibWeights = 1.0/nFibers*ones(nFibers);

  estimateModifiedBAS (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);
  //estimateBallAndSticks (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);
  fibComp.fibDirs.print("estimated dirs");
  fibComp.fibDiffs.print("estimated diffusivities");
  fibComp.fibWeights.print("estimated weights");

  if (bTrueDirsFound) {
    trueDirs.print("true dirs");
    double devAngle = mean(directionDeviation(fibComp.fibDirs, trueDirs))*180/M_PI;
    cout <<"deviation to true dirs: " <<devAngle <<endl;
  } else {
    cout <<"true dirs not found" <<endl;
  }
}

void testTomsAlgorithmParams (const StickEstimateOption &options)
{
  ifstream inFile ("Matlab/log.txt");
  ofstream outFile ("tomsalgorithmparams.txt");
  int nTrials = 10;
  char line [1280];

  mat gradMat0 = loadGradientOrientations("diffusion_directions.txt");
  ImageType::Pointer img = readImage("dwi-b1500.nii");
  double snr = 17.6;
  int bVal = 1500;

  /******parameter settings*********/
  mat diffus;
  diffus <<2.0e-3 <<1.6e-3 <<endr
	 <<2.0e-3 <<1.5e-3 <<endr
	 <<1.9e-3 <<1.6e-3 <<endr
	 <<2.0e-3 <<1.0e-3 <<endr
	 <<2.0e-3 <<0.5e-3 <<endr;

  while (!inFile.eof()) {
    inFile.getline(line, 1280, '\n');
    if (inFile.gcount()<1) {
      break;
    }

    stringstream ss (line);
    int x, y, nFibers;
    ss >>x >>y >>nFibers;
    if (nFibers<2) {
      continue;
    }

    cout <<"(" <<x <<"," <<y <<"): " <<nFibers <<endl;
    mat trueDirs (3, nFibers);
    for (int i=0; i<nFibers; ++i) {
      ss >>trueDirs(0,i) >>trueDirs(1,i) >>trueDirs(2,i);
    }

    outFile <<"(" <<x <<"," <<y <<"): " <<nFibers <<endl;

    mat gradMat = gradMat0;
    vec S = getVoxelValue(img, x-1, y-1, 0);
    double s0;
    trimSignalAndGrads(S, gradMat, s0);
    FiberComposition fibComp;
    fibComp.nFibers = nFibers;

    for (int iDif = 0; iDif<diffus.n_rows; ++iDif) {
      for (int iInit = 0; iInit<=1; ++iInit) {

	if (iInit==0) {
	  initRandom(fibComp.fibDirs, nFibers);
	} else {
	  fibComp.fibDirs = trueDirs;
	}
	fibComp.fibDiffs = diffus(iDif,0)*ones(nFibers+1);
	fibComp.fibDiffs(nFibers) = diffus(iDif,1);
	fibComp.fibWeights = 1.0/nFibers*ones(nFibers);
  
	vec devs = zeros(nTrials);
	for (int iTrial=0; iTrial<nTrials; ++iTrial) {
	  estimateModifiedBAS (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);
	  devs(iTrial) = mean(directionDeviation(fibComp.fibDirs, trueDirs));
	  cout <<"Trial #" <<iTrial <<": " <<devs(iTrial) <<endl;
	}
    
	char logStr [1280];
	sprintf(logStr, "\t%0.4f\t%0.4f\t%d\t%0.4f\t%0.4f\t%0.4f\n", diffus(iDif,0), diffus(iDif,1), iInit, mean(devs), stddev(devs), devs.min());
	outFile <<logStr;
      }
    }
  }

  inFile.close();
  outFile.close();
}

void testTomsAlgorithmOnPhantom (const StickEstimateOption &options, const char *folderName)
{
  mat gradMat0 = loadGradientOrientations("assets/phantom/diffusion_directions.txt");
  ImageType::Pointer img = readImage("assets/phantom/dwi-b1500.nii");
  double snr = 17.6;
  int bVal = 1500;
  int nRep = 10;
  char str [1024];
  char line [1024];

  ifstream inFile ("assets/phantom/mask.txt");
  sprintf(str, "%s/likelihood.txt\0", folderName);
  ofstream likeFile (str);

  while (!inFile.eof()) {
    inFile.getline(line, 1280, '\n');
    if (inFile.gcount()<1) {
      break;
    }

    stringstream ss (line);
    int x, y, nFibers;
    ss >>x >>y >>nFibers;

    cout <<"(" <<x <<"," <<y <<"): " <<nFibers <<endl;

    mat gradMat = gradMat0;
    vec S = getVoxelValue(img, x-1, y-1, 0);
    double s0;
    trimSignalAndGrads(S, gradMat, s0);
    FiberComposition fibComp;

    vec likeBuf (nRep);
    sprintf(str, "%s/res_%d_%d_%d.txt", folderName, x, y, nFibers);
    ofstream resFile (str);

    for (int i=0; i<nRep; ++i) {
      fibComp.nFibers = nFibers;
      initRandom(fibComp.fibDirs, nFibers);
      fibComp.fibDiffs = 2.5e-3*ones(nFibers+1);
      fibComp.fibDiffs(nFibers) = 1.2e-3;
      fibComp.fibWeights = 1.0/nFibers*ones(nFibers);
    
      likeBuf(i) = estimateModifiedBAS (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);
      mat resMat = join_cols(fibComp.fibDirs, fibComp.fibWeights.t());
      for (int j=0; j<resMat.n_rows; ++j) {
	for (int k=0; k<resMat.n_cols; ++k)
	  resFile <<resMat(j,k) <<"\t";
	resFile <<endl;
      }
      for (int j=0; j<nFibers+1; ++j)
	resFile <<fibComp.fibDiffs(j) <<"\t";
      resFile <<endl;
    }
    resFile.close();
    sprintf(str, "%d\t%d\t%d\t%0.2f\t%0.2f\n", x, y, nFibers, mean(likeBuf), stddev(likeBuf));
    likeFile <<str;
  }

  inFile.close();
  likeFile.close();
}

void testWeightEstimation (int argc, char **argv)
{

  if (argc<3) {
    cout <<"insufficient input" <<endl;
    exit (1);
  }

  StickEstimateOption options;
  options.maxIt = 50000;
  options.maxInnerIt = 5;
  options.init = 1;
  options.useManifold = true;
  options.step = 1e-8;
  options.kappaStep = 1e-10;
  options.kappa0Step = 1e-11;
  options.weightStep = 1e-8;
  options.tolerance = 1e-6;
  options.innerTolerance = 1e-5;
  options.isEstWeights = atoi(argv[1]);
  options.isEstDiffusivities = atoi(argv[2]);
  options.useLineSearch = false;

  int nRep = 20;
  int bVal = 2000;
  double s0 = 200;

  mat gradMat = loadGradientOrientations("gradients.txt");

  rowvec snrs;
  snrs <<20 <<40 <<endr;
  rowvec stickDiffus;
  stickDiffus <<1.3e-3 <<1.5e-3 <<1.7e-3 <<1.9e-3 <<2.1e-3 <<endr;
  rowvec ballDiffus;
  ballDiffus <<0.2e-3 <<0.3e-3 <<0.4e-3 <<0.5e-3 <<endr;

  for (int nFibers=2; nFibers<=2; ++nFibers) {
    rowvec sepAngles;
    rowvec weights;
    if (nFibers == 1) {
      sepAngles <<0 <<endr;
      weights <<1.0 <<endr;
    } else {
      sepAngles <<30 <<45 <<60 <<90 <<endr;
      weights <<0.3 <<0.4 <<0.5 <<endr;
    }

    for (int iSnr=0; iSnr<snrs.n_elem; ++iSnr) {
    for (int iSepAngle = 0; iSepAngle<sepAngles.n_elem; ++iSepAngle) {     
    for (int iStickDiffus=0; iStickDiffus<stickDiffus.n_elem; ++iStickDiffus) {
    for (int iBallDiffus=0; iBallDiffus<ballDiffus.n_elem; ++iBallDiffus) {
    for (int iWeight=0; iWeight<weights.n_elem; ++iWeight) {
      char str [1024];

      if (options.isEstWeights && options.isEstDiffusivities)
	sprintf(str, "results/synthetic/weights_diffus/weight_est__n=%d__s=%d__a=%d__d1=%0.1e__d0=%0.1e__w=%0.1f.txt", 
		nFibers, (int)snrs(iSnr), (int)sepAngles(iSepAngle), stickDiffus(iStickDiffus), ballDiffus(iBallDiffus), weights(iWeight));
      else if (options.isEstWeights)
	sprintf(str, "results/synthetic/weights/weight_est__n=%d__s=%d__a=%d__d1=%0.1e__d0=%0.1e__w=%0.1f.txt", 
	      nFibers, (int)snrs(iSnr), (int)sepAngles(iSepAngle), stickDiffus(iStickDiffus), ballDiffus(iBallDiffus), weights(iWeight));
      else
	sprintf(str, "results/synthetic/diffus/weight_est__n=%d__s=%d__a=%d__d1=%0.1e__d0=%0.1e__w=%0.1f.txt", 
	      nFibers, (int)snrs(iSnr), (int)sepAngles(iSepAngle), stickDiffus(iStickDiffus), ballDiffus(iBallDiffus), weights(iWeight));

      ofstream outFile (str);

      srand(time(0));
      double angle = sepAngles(iSepAngle)*M_PI/180;
      mat trueDirs (3, nFibers);
      if (nFibers==1) {
	trueDirs <<1 <<endr <<0 <<endr <<0 <<endr;
      } else {
	trueDirs <<1 <<cos(angle) <<endr <<0 <<sin(angle) <<endr <<0 <<0 <<endr;
      }

      mat R (3,3);
      initRandom(R, 3);
      trueDirs = R*trueDirs;

      vec trueWeights = ones(nFibers);
      if (nFibers>1)
	trueWeights <<weights(iWeight) <<endr <<1-weights(iWeight) <<endr;

      mat trueDiffus (3, nFibers);
      for (int i=0; i<nFibers; ++i) {
	trueDiffus(0,i) = stickDiffus(iStickDiffus);
	trueDiffus(1,i) = ballDiffus(iBallDiffus);
	trueDiffus(2,i) = ballDiffus(iBallDiffus);
      }
      FiberComposition trueFibComp;
      trueFibComp.nFibers = nFibers;
      trueFibComp.fibDirs = trueDirs;
      trueFibComp.fibWeights = trueWeights;
      trueFibComp.fibDiffs = stickDiffus(iStickDiffus)*ones(nFibers+1);
      trueFibComp.fibDiffs(nFibers) = ballDiffus(iBallDiffus);

      vec S = simulateMultiTensor(bVal, s0, gradMat, trueDirs, trueWeights, trueDiffus);
      vec S1 = addRicianNoise(S, s0/snrs(iSnr));

      for (int iRep=0; iRep<nRep; ++iRep) {	
	FiberComposition fibComp;
	fibComp.nFibers = nFibers;
	initRandom(fibComp.fibDirs, nFibers);
	fibComp.fibDiffs = 1.7e-3*ones(nFibers+1);
	fibComp.fibDiffs(nFibers) = 0.3e-3;
	fibComp.fibWeights = 1.0/nFibers*ones(nFibers);

	double l = estimateModifiedBAS (fibComp, S, gradMat, bVal, s0, snrs(iSnr), nFibers, options);
	vec dirDev, weightDev, diffusDev;
	double err = fiberDeviation(fibComp, trueFibComp, dirDev, weightDev, diffusDev);
	cout <<"mean error: " <<err <<endl;
	vec resVec = join_cols(dirDev, join_cols(weightDev, diffusDev));
	for (int i=0; i<resVec.n_elem; ++i) {
	  outFile <<resVec(i) <<"\t";
	}
	outFile <<l <<endl;
      }
      outFile.close();
    }}}}}
  }
}

vec readVoxelFromBrainData (int x, int y, int z, double &s0)
{
  ifstream brainFile ("assets/brain/im_b2000.raw");
  int sx = 106, sy = 106, sz = 76;
  int sVol = sx*sy*sz;

  vec S = zeros(64);

  for (int i=0; i<65; ++i) {
    brainFile.seekg((sVol*i+sx*sy*z+sx*y+x)*sizeof(float), ios::beg);
    float val;
    brainFile.read((char*)&val, sizeof(float));
    if (i==0)
      s0 = val;
    else
      S(i-1) = val;
  }
  return S;
}


void testTomsAlgorithmOnBrainOneVoxel (int argc, char **argv)
{
  int iArg = 1;
  int x= atoi(argv[iArg++]);
  int y = 40; //atoi(argv[iArg++]);
  int z = atoi(argv[iArg++]);
  int nFibers = atoi(argv[iArg++]);
  double sparse = atof(argv[iArg++]);
  testTomsAlgorithmOnBrainOneVoxel(x,y,z,nFibers,sparse);
}

void testTomsAlgorithmOnBrainOneVoxel (int x, int y, int z, int nFibers, double sparse)
{
  StickEstimateOption stickOptions;
  stickOptions.maxIt = 50000;
  stickOptions.maxInnerIt = 5;
  stickOptions.init = 1;
  stickOptions.useManifold = true;
  stickOptions.step = 1e-9;
  stickOptions.kappaStep = 1e-12;
  stickOptions.kappa0Step = 1e-10;
  stickOptions.weightStep = 5e-10;
  stickOptions.tolerance = 1e-5;
  stickOptions.innerTolerance = 1e-6;
  stickOptions.isEstWeights = true;
  stickOptions.isEstDiffusivities = false;
  stickOptions.useLineSearch = false;
  stickOptions.sparseFactor = sparse;
  stickOptions.isPrintDebugInfo = true;

  double snr = 9.4377;
  int bVal = 2000;

  double s0;
  vec S = readVoxelFromBrainData(x-1, y-1, z-1, s0);
  vec diffusis = -log(S/s0)/bVal;
  mat gradMat = loadGradientOrientations("assets/brain/gradients.txt");

  FiberComposition fibComp;
  
  FiberComposition fullTensor;
  estimateTensor(fullTensor, S, gradMat, bVal, s0);
  fullTensor.fibDirs.print("full tensor dirs");
  fullTensor.fibDiffs.print("full tensor diffus");

  fibComp.nFibers = nFibers;
  //initRandom(fibComp.fibDirs, nFibers);
  fibComp.fibDirs = fullTensor.fibDirs.cols(0,nFibers-1);
  double d1 = 2.0e-3;
  double d0 = 1.3e-3;
  fibComp.fibDiffs = d1*ones(nFibers+1);
  fibComp.fibDiffs(nFibers) = d0;

  d0 = diffusis.max();
  double smax = sqrt(M_PI/(4*bVal*d0))*erf(sqrt(bVal*d0));
  double smin = exp(-bVal*d0);
  double smean = mean(S)/s0;
  double w0 = 1-(smean-smin)/(smax-smin);
  cout <<"ball weight: " <<w0 <<endl;
  fibComp.fibWeights = (1.0-w0)/nFibers*ones(nFibers+1);
  fibComp.fibWeights(nFibers) = w0;

  fibComp.fibWeights = 1.0/(nFibers+1)*ones(nFibers+1);

  //estimateModifiedBAS (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);
  estimateBallAndSticks (fibComp, S, gradMat, bVal, s0, snr, nFibers, stickOptions);
  fibComp.fibDirs.print("estimated dirs");
  fibComp.fibDiffs.print("estimated diffusivities");
  fibComp.fibWeights.print("estimated weights");

  vec angDev = directionDeviation(fibComp.fibDirs, fullTensor.fibDirs.col(0))*180/M_PI;
  angDev.print("direction deviation");
}



void testTomsAlgorithmOnBrain (int argc, char **argv)
{
  if (argc<3) {
    cout <<"insufficient input" <<endl;
    exit (1);
  }

  StickEstimateOption options;
  options.maxIt = 50000;
  options.maxInnerIt = 5;
  options.init = 1;
  options.useManifold = true;
  options.step = 1e-8;
  options.kappaStep = 1e-10;
  options.kappa0Step = 1e-11;
  options.weightStep = 1e-8;
  options.tolerance = 1e-6;
  options.innerTolerance = 1e-5;
  options.isEstWeights = atoi(argv[1]);
  options.isEstDiffusivities = atoi(argv[2]);
  options.useLineSearch = false;

  // cout <<options.isEstWeights <<endl;
  // cout <<options.isEstDiffusivities <<endl;
  // return;

  mat gradMat = loadGradientOrientations("assets/brain/gradients.txt");
  double snr = 20;
  int bVal = 2000;
  int nRep = 20;
  char str [1024];
  char line [1024];
  
  double baselineMask = 100;
  const char *folderName;
  if (options.isEstWeights && options.isEstDiffusivities)
    folderName = "results/brain/weights_diffus";
  else if (options.isEstWeights)
    folderName = "results/brain/weights";
  else if (options.isEstDiffusivities)
    folderName = "results/brain/diffus";

  ifstream inFile ("assets/brain/mask.txt");
  sprintf(str, "%s/likelihood.txt\0", folderName);
  ofstream likeFile (str);

  while (!inFile.eof()) {
    inFile.getline(line, 1024, '\n');
    if (inFile.gcount()<1) {
      break;
    }

    stringstream ss (line);
    int x, y, z, nFibers;
    ss >>x >>y >>z >>nFibers;

    cout <<"(" <<x <<"," <<y <<"," <<z <<"): " <<nFibers <<endl;

    double s0;
    vec S = readVoxelFromBrainData(x,y,z, s0);

    if (s0<baselineMask)
      continue;
    
    FiberComposition fibComp;

    vec likeBuf (nRep);
    sprintf(str, "%s/res_%d_%d_%d_%d.txt", folderName, x, y, z, nFibers);
    ofstream resFile (str);

    for (int i=0; i<nRep; ++i) {
      fibComp.nFibers = nFibers;
      initRandom(fibComp.fibDirs, nFibers);
      fibComp.fibDiffs = 1.7e-3*ones(nFibers+1);
      fibComp.fibDiffs(nFibers) = 0.3e-3;
      fibComp.fibWeights = 1.0/(nFibers+1)*ones(nFibers+1);
    
      likeBuf(i) = estimateModifiedBAS (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);
      mat resMat = join_cols(fibComp.fibDirs, fibComp.fibWeights.t());
      for (int j=0; j<resMat.n_rows; ++j) {
	for (int k=0; k<resMat.n_cols; ++k)
	  resFile <<resMat(j,k) <<"\t";
	resFile <<endl;
      }
      for (int j=0; j<nFibers+1; ++j)
	resFile <<fibComp.fibDiffs(j) <<"\t";
      resFile <<endl;
    }
    resFile.close();
    sprintf(str, "%d\t%d\t%d\t%d\t%0.2f\t%0.2f\n", x, y, z, nFibers, mean(likeBuf), stddev(likeBuf));
    likeFile <<str;
  }

  inFile.close();
  likeFile.close();
}

void estimateCCDiffusivities ()
{
  StickEstimateOption stickOptions;
  stickOptions.maxIt = 50000;
  stickOptions.maxInnerIt = 5;
  stickOptions.init = 1;
  stickOptions.useManifold = true;
  stickOptions.step = 1e-9;
  stickOptions.kappaStep = 1e-13;
  stickOptions.kappa0Step = 1e-10;
  stickOptions.weightStep = 3e-10;
  stickOptions.tolerance = 1e-5;
  stickOptions.innerTolerance = 1e-6;
  stickOptions.isEstWeights = false;
  stickOptions.isEstDiffusivities = true;
  stickOptions.useLineSearch = false;
  stickOptions.sparseFactor = 0;
  stickOptions.isPrintDebugInfo = false;

  double snr = 9.4377;
  int bVal = 2000;
  mat gradMat = loadGradientOrientations("assets/brain/gradients.txt");

  ifstream posFile ("assets/brain/cc_b2000.txt");
  char line [1024];
  int x, y, z;
  int nFibers = 1;

  while(!posFile.eof()) {
    posFile.getline(line, 1024, '\n');
    if (posFile.gcount()<3)
      break;

    stringstream ss (line);
    ss >>x >>y >>z;
    //cout <<x <<"\t" <<y <<"\t" <<z <<endl;

    double s0;
    vec S = readVoxelFromBrainData(x-1, y-1, z-1, s0);
    vec diffusis = -log(S/s0)/bVal;
     
    FiberComposition fibComp;
    
    FiberComposition fullTensor;
    estimateTensor(fullTensor, S, gradMat, bVal, s0);
    
    fibComp.nFibers = nFibers;
    //initRandom(fibComp.fibDirs, nFibers);
    fibComp.fibDirs = fullTensor.fibDirs.cols(0,nFibers-1);
    double d1 = 2.0e-3;
    double d0 = 1.3e-3;
    fibComp.fibDiffs = d1*ones(nFibers+1);
    fibComp.fibDiffs(nFibers) = d0;
    
    d0 = diffusis.max();
    double smax = sqrt(M_PI/(4*bVal*d0))*erf(sqrt(bVal*d0));
    double smin = exp(-bVal*d0);
    double smean = mean(S)/s0;
    double w0 = 1-(smean-smin)/(smax-smin);
    
    // fibComp.fibWeights = (1.0-w0)/nFibers*ones(nFibers+1);
    // fibComp.fibWeights(nFibers) = w0;
    fibComp.fibWeights = 1.0/(nFibers+1)*ones(nFibers+1);

    estimateBallAndSticks (fibComp, S, gradMat, bVal, s0, snr, nFibers, stickOptions);
    for (int i=0; i<nFibers+1; ++i) {
      cout <<fibComp.fibWeights(i) <<"\t";
    }
    for (int i=0; i<nFibers+1; ++i) {
      cout <<fibComp.fibDiffs(i) <<"\t";
    }
    cout <<endl;
  }
}
