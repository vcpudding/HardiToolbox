#include "../HardiToolbox.h"

using namespace std;
using namespace HardiToolbox;


#include "itkVectorImage.h"
#include "itkImageFileReader.h"

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

int main (int argc, char **argv)
{
  int lineNum = atoi(argv[1]);
  bool bEstWeights = atoi(argv[2]);
  bool bEstDiffus = atoi(argv[3]);
  string method = string(argv[4]);  

  StickEstimateOption options;
  options.maxIt = 50000;
  options.maxInnerIt = 5;
  options.init = 1;
  options.useManifold = true;
  options.step = 1e-9;
  options.kappaStep = 1e-12;
  options.kappa0Step = 1e-10;
  options.weightStep = 1e-11;
  options.tolerance = 1e-5;
  options.innerTolerance = 1e-6;
  options.isEstWeights = bEstWeights;
  options.isEstDiffusivities = bEstDiffus;
  options.useLineSearch = false;

  mat gradMat = loadGradientOrientations("assets/phantom/diffusion_directions.txt");
  double snr = 17.6;
  int bVal = 1500;
  int nRep = 20;
  char str [1024];
  char line [1024];
  
  double baselineMask = 100;
  char folderName [1024];
  if (options.isEstWeights && options.isEstDiffusivities)
    sprintf(folderName, "results/phantom/%s/weights_diffus", method.c_str());  
  else if (options.isEstWeights)
    sprintf(folderName, "results/phantom/%s/weights", method.c_str());
  else if (options.isEstDiffusivities)
    sprintf(folderName, "results/phantom/%s/diffus", method.c_str());

  ifstream inFile ("assets/phantom/mask.txt");

  //while (!inFile.eof()) {
  for (int i=0; i<=lineNum; ++i) {
    inFile.getline(line, 1024, '\n');
    if (inFile.gcount()<1) {
      return 0;
    }
  }

  stringstream ss (line);
  int x, y, nFibers;
  ss >>x >>y >>nFibers;

  cout <<"(" <<x <<"," <<y <<"): " <<nFibers <<endl;

  ImageType::Pointer img = readImage("assets/phantom/dwi-b1500.nii");
  vec S = getVoxelValue(img, x-1, y-1, 0);
  double s0;
  trimSignalAndGrads(S, gradMat, s0);
    
  FiberComposition fibComp;

sprintf(str, "%s/res_%d_%d_%d_%d.txt", folderName, x, y, 0, nFibers);
  ofstream resFile (str);

  for (int i=0; i<nRep; ++i) {
    fibComp.nFibers = nFibers;
    initRandom(fibComp.fibDirs, nFibers);
    fibComp.fibDiffs = 2.5e-3*ones(nFibers+1);
    fibComp.fibDiffs(nFibers) = 1.2e-3;
    if (method=="bas")
      fibComp.fibWeights = 1.0/(nFibers+1)*ones(nFibers+1);
    else
      fibComp.fibWeights = 1.0/nFibers*ones(nFibers);
    
    double l ;
    if (method=="bas")
      l = estimateBallAndSticks (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);
    else
      l = estimateModifiedBAS (fibComp, S, gradMat, bVal, s0, snr, nFibers, options);

    fibComp.fibDirs.reshape(3*nFibers,1);
    mat resMat = join_cols(fibComp.fibDirs,
			   join_cols(fibComp.fibWeights, fibComp.fibDiffs));
    for (int j=0; j<resMat.n_elem; ++j)
      resFile <<resMat(j,0) <<"\t";
    resFile <<l <<endl;
  }
  resFile.close();

  inFile.close();
  return 0;
}
