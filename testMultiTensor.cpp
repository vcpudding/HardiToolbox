#include <iostream>
#include <armadillo>
#include <string>
#include <math.h>

#include "HardiToolbox.h"

using namespace std;
using namespace arma;
using namespace HardiToolbox;

#define _USE_MATH_DEFINES

void parseArgs (int argc, char ** argv, int &nSimFibers, int &maxEstFibers, int &bVal, int &nRep, int &init, const char **gradFileName)
{
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-help") {
      cout <<"Experiment options:" <<endl;
      cout <<"-nsim:\tnumber of fibers to simulate" <<endl;
      cout <<"-nmaxest:\tmaximum number of fibers to estimate (minimum is set by -nsim)" <<endl;
      cout <<"-bval:\tb-value" <<endl;
      cout <<"-numrep:\tnumber of repetitions per set of parameters" <<endl;
      cout <<"-init:\tinitialization method. 0: random init; 1: init with full tensor" <<endl;
      cout <<"-gradfile:\tpath to gradient direction file" <<endl;
      exit(0);
    }
    if (i+1 != argc) { // Check that we haven't finished parsing already
      if (string(argv[i]) == "-nsim") {
	nSimFibers = atoi(argv[++i]);
      } else if (string(argv[i]) == "-nmaxest") {
	maxEstFibers = atoi(argv[++i]);
	if (maxEstFibers<nSimFibers) {
	  cout <<"Max number of estimated fibers must not be less than number of simulated fibers." <<endl;
	  maxEstFibers = nSimFibers;
	}
      } else if (string(argv[i]) == "-bval") {
	bVal = atoi(argv[++i]);
      } else if (string(argv[i]) == "-numrep") {
	nRep = atoi(argv[++i]);
      } else if (string(argv[i]) == "-gradfile") {
	*gradFileName = argv[++i];
      } else if (string(argv[i]) == "-init") {
	init = atoi(argv[++i]);
      }
    }
  }
}

int main (int argc, char **argv)
{
  int nRep = 200;
  int bVal = 3000;
  int nSimFibers = 2;
  int nMaxEst = 2;
  int init = 1;
  double s0 = 1.0;

  rowvec sepAngles;
  rowvec snrs;
  rowvec weights;

  sepAngles <<30 <<45 <<60 <<75 <<90 <<endr;
  snrs <<20 <<40 <<endr;
  weights <<0.5 <<endr;

  const char * gradFileName = "gradients.txt";
  
  parseArgs(argc, argv, nSimFibers, nMaxEst, bVal, nRep, init, &gradFileName);

  cout <<"Number of simulated fibers: " <<nSimFibers <<endl;
  cout <<"Maximum number of estimated fibers: " <<nMaxEst <<endl;
  cout <<"b-value: " <<bVal <<endl;
  cout <<"Number of repetitions: " <<nRep <<endl;
  cout <<"Initialization: " <<(init==0?"random": "full tensor") <<endl;
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
  options.maxIt = 50000;
  options.init = init;
  options.step = 1e-3;
  options.tolerance = 1e-6;

  for (int iSepAngle = 0; iSepAngle<sepAngles.n_elem; ++iSepAngle) {
    for (int iSnr = 0; iSnr<snrs.n_elem; ++iSnr) {
      for (int iWeight = 0; iWeight<weights.n_elem; ++iWeight) {

	double angle = sepAngles(iSepAngle)*3.1415926/180.0;
	double snr = snrs(iSnr);
	double weight = weights(iWeight);

	char dirFileName [1024];
	char weightFileName [1024];
	sprintf(dirFileName, "results/n=%d__est=%d__b=%d__g=%d__a=%0.0f__s=%0.0f__w=%0.1f__init=%d_d.txt", 
		nSimFibers, nMaxEst, bVal, nGrads, sepAngles(iSepAngle), snr, weight, init);
	sprintf(weightFileName, "results/n=%d__est=%d__b=%d__g=%d__a=%0.0f__s=%0.0f__w=%0.1f__init=%d_w.txt", 
		nSimFibers, nMaxEst, bVal, nGrads, sepAngles(iSepAngle), snr, weight, init);
	ofstream dirFile (dirFileName);
	ofstream weightFile (weightFileName);

	mat trueDirs (3,nSimFibers);
	vec trueWeights;

	if (nSimFibers==1) {
		trueDirs <<1 <<endr <<0 <<endr <<0 <<endr;
		trueWeights <<1 <<endr;
	} else if (nSimFibers==2) {
		trueDirs <<1 <<cos(angle) <<endr <<0 <<sin(angle) <<endr <<0 <<0 <<endr;
		trueWeights <<weight <<endr <<1.0-weight <<endr;
	} else {
		initFixed(trueDirs, 3);
		trueWeights = 1.0/3.0*ones(3);
	}

	mat R (3,3);
	initRandom(R, 3);
	trueDirs = R*trueDirs;

	vec S = simulateMultiTensor(bVal, s0, gradientOrientations, trueDirs, trueWeights);
	
	for (int i=0; i<nRep; ++i) {
	  vec S1 = addRicianNoise (S, 1.0/snr);

	  mat dirMat = trueDirs;
	  mat weightMat = trueWeights.t();

	  for (int n=nSimFibers; n<=nMaxEst; ++n) {
		FiberComposition fibComp;
		estimateMultiTensor(fibComp, S1, gradientOrientations, bVal, s0, n, options);
		dirMat = join_rows(dirMat, fibComp.fibDirs);
		weightMat = join_rows(weightMat, fibComp.fibWeights.t());
	  }
	  
	  for (int k=0; k<3; ++k) {
	    for (int j=0; j<dirMat.n_cols; ++j) {
	      dirFile <<dirMat(k,j) <<" ";
	    }
	    dirFile <<endl;
	  }
	  
	  for (int j=0; j<weightMat.n_cols; ++j) {
	    weightFile <<weightMat(j) <<" ";
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
