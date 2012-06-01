#include <iostream>
#include <armadillo>
#include <string>

#include "HardiToolbox.h"

using namespace std;
using namespace arma;
using namespace HardiToolbox;


void parseArgs (int argc, char ** argv, int &nFibers, int &bVal, int &nRep, const char **gradFileName)
{
  for (int i = 1; i < argc; i++) {
    if (i + 1 != argc) { // Check that we haven't finished parsing already
      if (string(argv[i]) == "-n") {
	nFibers = atoi(argv[++i]);
      } else if (string(argv[i]) == "-b") {
	bVal = atoi(argv[++i]);
      } else if (string(argv[i]) == "-r") {
	nRep = atoi(argv[++i]);
      } else if (string(argv[i]) == "-g") {
	*gradFileName = argv[++i];
      } else {
	std::cout << "Not enough or invalid arguments, please try again.\n";
	exit(0);
      }
    }
  }
}

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
  
  parseArgs(argc, argv, nFibers, bVal, nRep, &gradFileName);

  cout <<"Number of fibers: " <<nFibers <<endl;
  cout <<"b-value: " <<bVal <<endl;
  cout <<"Number of repititions: " <<nRep <<endl;
  cout <<"Gradient file: " <<gradFileName <<endl;

  cout <<"continue? (y/n) ";
  char c;
  cin >>c;
  if (c!='y') {
    return 0;
  }

  mat gradientOrientations = loadGradientOrientations (gradFileName);
  int nGrads = gradientOrientations.n_rows;

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
	sprintf(dirFileName, "results/n=%d__b=%d__g=%d__a=%0.0f__s=%0.0f__w=%0.1f_d.txt", nFibers, bVal, nGrads, sepAngles(iSepAngle), snr, weight);
	sprintf(weightFileName, "results/n=%d__b=%d__g=%d__a=%0.0f__s=%0.0f__w=%0.1f_w.txt", nFibers, bVal, nGrads, sepAngles(iSepAngle), snr, weight);
	ofstream dirFile (dirFileName);
	ofstream weightFile (weightFileName);


	mat trueDirs (3,nFibers);
	vec trueWeights;

	if (nFibers==1) {
	  trueDirs <<1 <<endr <<0 <<endr <<0 <<endr;
	  trueWeights <<1 <<endr;
	} else if (nFibers==2) {
	  trueDirs <<1 <<cos(angle) <<endr <<0 <<sin(angle) <<endr <<0 <<0 <<endr;
	  trueWeights <<weight <<endr <<1.0-weight <<endr;
	} else {
	  initFixed(trueDirs, 3);
	  trueWeights = 1.0/3.0*ones(3);
	}

	mat R (3,3);
	initRandom(R, 3);
	trueDirs = R*trueDirs;

	vec S = simulateMultiTensor(bVal, s0, gradientOrientations, trueDirs, trueWeights, true);
	
	for (int i=0; i<nRep; ++i) {
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

	cout <<dirFileName <<" saved." <<endl;
	cout <<weightFileName <<" saved." <<endl;
      }
    }
  }
  return 0;
}
