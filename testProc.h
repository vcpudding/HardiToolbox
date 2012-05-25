#ifndef TESTPROC_H
#define TESTPROC_H

#include <time.h>
#include <armadillo>
using namespace arma;

#include "HardiToolbox.h"
#include "utils.h"


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
  
#endif
