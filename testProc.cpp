#include "testProc.h"

#include <fstream>

using namespace std;
using namespace HardiToolbox;

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
      trueDirs.resize(3, nFibers);
      for (int i=0; i<nFibers; ++i) {
	ss >>trueDirs(0,i) >>trueDirs(1,i) >>trueDirs(2,i);
      }
      fibComp.nFibers = nFibers;
      fibComp.fibDirs = trueDirs;
      fibComp.fibDiffs = 2.0e-3*ones(nFibers+1);
      fibComp.fibDiffs(nFibers) = 1.0e-3;
      fibComp.fibWeights = 1.0/nFibers*ones(nFibers);
      bTrueDirsFound = true;
      break;      
    }
  }

  // estimateModifiedBASByStick (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);
  estimateModifiedBAS (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);
  fibComp.fibDirs.print("estimated dirs");
  fibComp.fibWeights.print("estimated weights");
  fibComp.fibDiffs.print("estimated diffusivities");

  if (bTrueDirsFound) {
    trueDirs.print("true dirs");
    vec devAngle = directionDeviation(fibComp.fibDirs, trueDirs)*180/M_PI;

    devAngle.print("deviation to true dirs");
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
  //double snr = 17.6;
  int bVal = 1500;

  /******parameter settings*********/
  mat diffus = 2.0*ones(15, 2);
  for (int i=0; i<diffus.n_rows; ++i) {
    diffus(i,1) = i*0.1e-3+0.4e-3;
  }

  vec snrs;
  snrs <<17.6 <<endr;

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
    for (int iSnr = 0; iSnr<snrs.n_elem; ++iSnr) {
      for (int iInit = 0; iInit<=0; ++iInit) {

	if (iInit==0) {
	  initRandom(fibComp.fibDirs, nFibers);
	} else {
	  fibComp.fibDirs = trueDirs;
	}
	fibComp.fibDiffs = diffus(iDif,0)*ones(nFibers+1);
	fibComp.fibDiffs(nFibers) = diffus(iDif,1);
	fibComp.fibWeights = 1.0/nFibers*ones(nFibers);
  
	double snr = snrs(iSnr);

	vec devs = zeros(nTrials);
    
	char logStr [1280];
	sprintf(logStr, "%0.4f\t%0.4f\t%0.0f\n", diffus(iDif,0), diffus(iDif,1), snr);
	outFile <<logStr;
	for (int iTrial=0; iTrial<nTrials; ++iTrial) {
	  estimateModifiedBAS (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);
	  for (int i=0; i<3; ++i) {
	    outFile <<"\t" <<fibComp.fibDirs(i,0) <<"\t" <<fibComp.fibDirs(i,1) <<endl;
	  }
	  outFile <<"\t---------------------" <<endl;
	  devs(iTrial) = mean(directionDeviation(fibComp.fibDirs, trueDirs));
	  cout <<"Trial #" <<iTrial <<": " <<devs(iTrial) <<endl;
	}
      }
    }
    }
  }

  inFile.close();
  outFile.close();
}
