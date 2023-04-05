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
observable testObs = observable();
int nParams = 0;

double minfunc(const double *xx){
	return 1.0+xx[0]*xx[0]+pow(xx[1]-1.5,2)+pow(xx[2]+4,2);

	vector<double> params = {};
	for(int i = 0; i < nParams; i++){
		params.push_back(xx[i]);
	}

	testObs.setFitParams(params);

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
	testReader.loadExpData();


	//saves the observable object outside of filereader object
	testObs = testReader.getObs();

	/*
	JPsi testJ = JPsi();
	
	function<double(double*)> jminf = [&](double *x){
		return x[0]*x[0];
		//return testJ.JPsiminfunc({x[0]});
	};
	*/
	vector<double> sols = {0,1.5,-4.0};
	ROOT::Math::Functor f(&minfunc,sols.size());
	

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
	nParams = fitparams.size();

	min->SetFunction(f);
	/*
	//make a function wrapper to minimize the function minfunc (=chisquared)
	ROOT::Math::Functor f(&minfunc,nParams);
	min->SetFunction(f);
	//set the initial conditions and step sizes
	for(int i = 0; i < nParams; i++){
		min->SetVariable(i,to_string(i),fitparams[i],steps[i]);
	}
	*/
	min->SetVariable(0,"x",5,0.2);
	min->SetVariable(1,"y",32,0.2);
	min->SetVariable(2,"z",20,0.2);
	min->Minimize();

	//store the parameters for the minimum that the minimizer found in xs
	const double *xs = min->X();

	//note to self: need to get rid of 'dumbJ' in amplitude.cpp later when doing non-radJPsi fits
}
