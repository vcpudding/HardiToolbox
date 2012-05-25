#include <iostream>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include "HardiToolbox.h"

using namespace std;
using namespace HardiToolbox;

#include "itkVectorImage.h"
#include "itkImageFileReader.h"

#include "testProc.h"


typedef itk::Image<double, 4>  ImageType; 
typedef itk::ImageFileReader<ImageType> ReaderType;

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

int main (int argc, char **argv)
{
  // testSparsityTerm();
  // return 0;
  struct timeval tStart, tEnd;

  int bVal = 3000;
  double s0 = 220;
  double d = 1.7e-3;
  double snr = 10;
  double angle = M_PI/3;
  int nFibers = 2;
  int nTrials = 10;
  vec fibDirs;
  vec weights;
  vec diffusivities;

  if (nFibers == 1)
  {
    fibDirs <<angle <<endr;
    weights <<1  <<endr;
    diffusivities <<1.4e-3 <<endr <<3e-4 <<endr;
  } else
  {
    fibDirs <<0 <<endr <<angle <<endr;
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
  stickOptions.maxIt = 55000;
  stickOptions.maxInnerIt = 5;
  stickOptions.init = 0;
  stickOptions.useManifold = true;
  stickOptions.step = 1e-8;
  stickOptions.kappaStep = 1e-12;
  stickOptions.kappa0Step = 1e-12;
  stickOptions.tolerance = 1e-3;
  stickOptions.innerTolerance = 1e-6;
  stickOptions.isEstWeights = false;
  stickOptions.isEstDiffusivities = true;

  FiniteGaussianOption fgOptions;
  fgOptions.maxIt = 5000;
  fgOptions.step = 1e-6;
  fgOptions.tolerance = 1e-6;

  FiberComposition fibComp1;
  FiberComposition fibComp2;
  FiberComposition fibComp3;

  // mat gradMat = loadGradientOrientations("assets/Gradient64.txt");
  // vec S = simulateMultiTensor(bVal, s0, gradMat, fibDirs, weights);
  // vec S1 = addRicianNoise(S, s0/snr);

  mat gradMat = loadGradientOrientations("diffusion_directions.txt");
  ImageType::Pointer img = readImage("dwi-b1500.nii");
  if (argc<4) {
    cout <<"input voxel locations x and y, number of fibers" <<endl;
    return 0;
  }
  int x = atoi(argv[1]);
  int y = atoi(argv[2]);
  nFibers = atoi(argv[3]);
  vec S1 = getVoxelValue(img, x-1, y-1, 1);
  trimSignalAndGrads(S1, gradMat, s0);
  snr = 15;
  bVal = 1500;

  UtilStopWatch::tic();

  // estimateSticksDiffs(fibComp1, S1, gradMat, bVal, s0, snr, nFibers, stickOptions);
//  estimateSticks(fibComp2, S1, gradMat, bVal, s0, d, snr, nFibers, stickOptions);
  // estimateBallAndSticksDiffs(fibComp2, S1, gradMat, bVal, s0, snr, nFibers, stickOptions);
  // estimateFiniteGaussian(fibComp2, S1, gradMat, bVal, s0, nFibers, fgOptions);

  deconvolveFibers(fibComp1, S1/s0, gradMat, bVal, nFibers, options);
  // estimateTensor (fibComp1, S1, gradMat, bVal, s0);
  fibComp2 = fibComp1;
  estimateModifiedBAS (fibComp2, S1, gradMat, bVal, s0, snr, nFibers, stickOptions);
  // estimateBallAndSticks(fibComp2, S1, gradMat, bVal, s0, snr, nFibers, stickOptions);

//  cout <<endl;
//  cout <<"elapsed time 1: " <<elapsedTime1 <<endl;
//  cout <<"elapsed time 2: " <<elapsedTime2 <<endl;
//
  // weights.t().print("true weights");
  fibComp1.fibWeights.print("estimated weights 1");
  fibComp2.fibWeights.print("estimated weights 2");
  // fibComp3.fibWeights.print("estimated weights 3");

  fibComp1.fibDirs.print("estimated directions 1");
  fibComp2.fibDirs.print("estimated directions 2");
  // fibComp3.fibDirs.print("estimated directions 3");

  fibComp1.fibDiffs.print("estimated diffusivities 1");
  fibComp2.fibDiffs.print("estimated diffusivities 2");
  // fibComp3.fibDiffs.print("estimated diffusivities 3");

  // mat trueDirs;
  // // trueDirs.t().print("true dirs");
  // double dev1 = directionDeviation(fibComp1.fibDirs, trueDirs.t());
  // // double dev2 = directionDeviation(fibComp2.fibDirs, trueDirs.t());
  // // double dev3 = directionDeviation(fibComp3.fibDirs, trueDirs.t());
  cout <<"direction deviation: " <<directionDeviation(fibComp1.fibDirs, fibComp2.fibDirs) <<endl;
  // cout <<"sparsity: " <<options.sparse <<endl;
  cout <<"timer: " <<UtilStopWatch::toc() <<endl;
  return 0;
}