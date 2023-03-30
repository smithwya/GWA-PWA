#include <iostream>
#include <complex>
#include <vector>
#include <Eigen/Dense>
#include <sstream>
#include <fstream>
#include <map>
#include <chrono>
#include "amplitude.h"
#include "channel.h"
#include "filereader.h"
#include "observable.h"
#include "JPsi.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"


using namespace std;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;


double minfunc(const double *xx){

	return 0;
}



int main()
{
	//TODO: data reading, write chisquared function 'minfunc'

	//reads the file and creates an observable object with the information from the file
	filereader testReader("Data/GWA_dataformat.txt");
	testReader.SetAllCommandLists();
	testReader.ConstructBareAmps();
	testReader.setChebys();
	testReader.setPoles();
	testReader.setKmats();

	//saves the observable object outside of filereader object
	observable testObs = testReader.getObs();

	//make the minimzer
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2","");
	//Set some criteria for the minimzer to stop
	min->SetMaxFunctionCalls(1000000);
	min->SetMaxIterations(10000);
	min->SetTolerance(0.001);
	min->SetPrintLevel(1);
	//get the initial parameters and steps from the constructed observable object
	vector<double> fitparams = testObs.getFitParams();
	vector<double> steps = testObs.getStepSizes();
	int nParams = fitparams.size();

	//make a function wrapper to minimize the function minfunc (=chisquared)
	ROOT::Math::Functor f(&minfunc,nParams);
	min->SetFunction(f);
	//set the initial conditions and step sizes
	for(int i = 0; i < nParams; i++){
		min->SetVariable(i,to_string(i),fitparams[i],steps[i]);
	}

	min->Minimize();

	//store the parameters for the minimum that the minimizer found in xs
	const double *xs = min->X();

	//note to self: need to get rid of 'dumbJ' in amplitude.cpp later when doing non-radJPsi fits
}
