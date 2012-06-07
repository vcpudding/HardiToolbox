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
      trueDirs.resize(3, nFibers0);
      for (int i=0; i<nFibers0; ++i) {
	ss >>trueDirs(0,i) >>trueDirs(1,i) >>trueDirs(2,i);
      }
      bTrueDirsFound = true;
      break;      
    }
  }

  fibComp.nFibers = nFibers;
  initRandom(fibComp.fibDirs, nFibers);
  fibComp.fibDiffs = 2.5e-3*ones(nFibers+1);
  fibComp.fibDiffs(nFibers) = 1.2e-3;
  fibComp.fibWeights = 1.0/nFibers*ones(nFibers);

  estimateModifiedBAS (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);
  // estimateBallAndSticks (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);
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

void testWeightEstimation (const StickEstimateOption &options)
{
  int nRep = 10;
  int bVal = 2000;
  double s0 = 1000;

  mat gradMat = loadGradientOrientations("gradients.txt");

  rowvec snrs;
  snrs <<20 <<40 <<endr;
  rowvec stickDiffus;
  stickDiffus <<1.3e-3 <<1.5e-3 <<1.7e-3 <<1.9e-3 <<2.1e-3 <<endr;
  rowvec ballDiffus;
  ballDiffus <<0.2e-3 <<0.3e-3 <<0.4e-3 <<0.5e-3 <<endr;

  for (int nFibers=1; nFibers<=2; ++nFibers) {
    rowvec sepAngles;
    if (nFibers == 1)
      sepAngles <<0 <<endr;
    else
      sepAngles <<30 <<45 <<60 <<90 <<endr;

    for (int iSnr=0; iSnr<snrs.n_elem; ++iSnr) {
    for (int iSepAngle = 0; iSepAngle<sepAngles.n_elem; ++iSepAngle) {     
    for (int iStickDiffus=0; iStickDiffus<stickDiffus.n_elem; ++iStickDiffus) {
    for (int iBallDiffus=0; iBallDiffus<ballDiffus.n_elem; ++iBallDiffus) {
      char str [1024];
      sprintf(str, "results/synthetic/weights/weight_est__n=%d__s=%d__a=%d_d1=%0.1e_d0=%0.1e.txt", 
	      nFibers, (int)snrs(iSnr), (int)sepAngles(iSepAngle), stickDiffus(iStickDiffus), ballDiffus(iBallDiffus));
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

      vec trueWeights = 1.0/nFibers*ones(nFibers);
      if (nFibers>1) {
	vec r = randu(nFibers)-0.5;
	r -= mean(r);
	trueWeights += 0.2*r;
      }

      mat trueDiffus (3, nFibers);
      for (int i=0; i<nFibers; ++i) {
	trueDiffus(0,i) = stickDiffus(iStickDiffus);
	trueDiffus(1,i) = ballDiffus(iBallDiffus);
	trueDiffus(2,i) = ballDiffus(iBallDiffus);
      }

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
	vec dev = directionDeviation(fibComp.fibDirs, trueDirs)*180/M_PI;
	vec resVec = join_cols(trueWeights, join_cols(fibComp.fibWeights, dev));
	for (int i=0; i<resVec.n_elem; ++i) {
	  outFile <<resVec(i) <<"\t";
	}
	outFile <<l <<endl;
      }
      outFile.close();
    }}}}
  }
}
