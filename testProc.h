#ifndef TESTPROC_H
#define TESTPROC_H

#include <time.h>
#include <armadillo>
using namespace arma;

#include "HardiToolbox.h"
#include "utils.h"

#include "itkVectorImage.h"
#include "itkImageFileReader.h"

using namespace HardiToolbox;

typedef itk::Image<double, 4>  ImageType; 
typedef itk::ImageFileReader<ImageType> ReaderType;

ImageType::Pointer readImage (const char *fileName);
vec getVoxelValue (ImageType::Pointer img, int x, int y, int z);
void estimate3DVolume (const char *dwiFileName, const char *maskFileName, const char *gradientFileName, 
		       const char *saveFileName, int bVal, double snr, int nFibers, void *opt);

void trimSignalAndGrads (vec &S, mat &gradMat, double &s0);
vec readVoxelFromBrainData (int x, int y, int z, double &s0);

void testBASRicianEM ();
void testModBASRicianEM ();

void testSparsityTerm ();
void testTomsAlgorithm (const StickEstimateOption &options);
void testTomsAlgorithmOneVoxel (int x, int y, int nFibers, const StickEstimateOption &options);
void testTomsAlgorithmParams (const StickEstimateOption &options);
void testTomsAlgorithmOnPhantom (const StickEstimateOption &options, const char *folderName);
void testWeightEstimation (int argc, char **argv);
void testTomsAlgorithmOnBrain (int argc, char **argv);
void testTomsAlgorithmOnBrainOneVoxel (int argc, char **argv);
void testTomsAlgorithmOnBrainOneVoxel (int x, int y, int z, int nFibers, double sparse);
void estimateCCDiffusivities();
  
#endif
