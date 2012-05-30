#include <iostream>
#include <armadillo>

#include "HardiToolbox.h"

using namespace std;
using namespace arma;
using namespace HardiToolbox;

int main (int argc, char **argv)
{
  int nRep = 200;
  int bVal = 3000;
  int nFibers = 2;
  double s0 = 1.0;

  rowvec sepAngles;
  rowvec snrs;
  rowvec weights;

  sepAngles <<30 <<45 <<60 <<75 <<90 <<endr;
  snrs <<20 <<40 <<endr;
  weights <<0.5 <<endr;

  const char * gradFileName = "gradients.txt";
  if (argc>1) {
    gradFileName = argv[1];
  }
  cout <<gradFileName <<endl;

  mat gradientOrientations = loadGradientOrientations (gradFileName);
  int nGrads = gradientOrientations.n_rows;

  gradientOrientations.print("gradients");
  return 0;

  MultiTensorOption options;
  options.maxIt = 5000;
  options.step = 1e-2;
  options.tolerance = 1e-6;

  for (int iSepAngle = 0; iSepAngle<sepAngles.n_elem; ++iSepAngle) {
    for (int iSnr = 0; iSnr<snrs.n_elem; ++iSnr) {
      for (int iWeight = 0; iWeight<weights.n_elem; ++iWeight) {

	double angle = sepAngles(iSepAngle)*M_PI/180.0;
	double snr = snrs(iSnr);
	double weight = weights(iWeight);

	char dirFileName [1024];
	char weightFileName [1024];
	sprintf(dirFileName, "results/g=%d__a=%0.0f__s=%0.0f__w=%0.1f_d.txt", nGrads, sepAngles(iSepAngle), snr, weight);
	sprintf(weightFileName, "results/g=%d__a=%0.0f__s=%0.0f__w=%0.1f_w.txt", nGrads, sepAngles(iSepAngle), snr, weight);
	ofstream dirFile (dirFileName);
	ofstream weightFile (weightFileName);
	
	for (int i=0; i<nRep; ++i) {

	  mat trueDirs (3,nFibers);
	  vec trueWeights;

	  if (nFibers==1) {
	    trueDirs <<1 <<endr <<0 <<endr <<0 <<endr;
	    trueWeights <<1 <<endr;
	  } else {
	    trueDirs <<1 <<cos(angle) <<endr <<0 <<sin(angle) <<endr <<0 <<0 <<endr;
	    trueWeights <<weight <<endr <<1.0-weight <<endr;
	  }
	  mat R (3,3);
	  initRandom(R, 3);
	  trueDirs = R*trueDirs;

	  trueDirs.print("true dirs");

	  vec S = simulateMultiTensor(bVal, s0, gradientOrientations, trueDirs, trueWeights, true);
	  vec S1 = addRicianNoise (S, 1.0/snr);

	  FiberComposition fibComp;
	  estimateMultiTensor(fibComp, S1, gradientOrientations, bVal, s0, nFibers, options);
	  
	  for (int k=0; k<3; ++k) {
	    for (int j=0; j<nFibers; ++j) {
	      dirFile <<fibComp.fibDirs(k,j) <<" ";
	    }
	    for (int j=0; j<nFibers; ++j) {
	      dirFile <<trueDirs(k,j) <<" ";
	    }
	    dirFile <<endl;
	  }
	  
	  for (int j=0; j<nFibers; ++j) {
	    weightFile <<fibComp.fibWeights(j) <<" ";
	  }
	  weightFile <<endl;

	}

	dirFile.close();
	weightFile.close();
      }
    }
  }
  return 0;
}
