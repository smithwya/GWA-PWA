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
#include "bottomonium.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom3.h"
#include "TError.h"


using namespace std;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;
observable testObs = observable();
int nParams = 0;
vector<comp> temppoles = {};

double minfunc(const double *xx){

	vector<double> params = {};
	for(int i = 0; i < nParams; i++){
		params.push_back(xx[i]);
	}

	testObs.setFitParams(params);
	return testObs.chisq(params);
}

int main(int argc, char *argv[])
{

	//reads the file and creates an observable object with the information from the file
	
	//filereader testReader("Data/simpledat.txt");
	if (argc == 1) { cout << "Missing filename" << endl; exit(0);}
	filereader testReader(argv[1]);
	testReader.SetAllCommandLists();
	testReader.ConstructBareAmps();
	testReader.setChebys();
	testReader.setPoles();
	testReader.setKmats();
	testReader.loadExpData();
	if(testReader.getRandomize()) testReader.randomize(); 

	//saves the observable object outside of filereader object
	testObs = testReader.getObs();

	//print out the observable starting params
	cout << testObs.amplitudes[0] << endl; 
	vector<double> params = testObs.getFitParams();

	//Giorgio's graphing shit
	auto intensityP_BB = [&](double x){
		comp value = testObs.amplitudes[0].getValue(pow(x,2))(0);
		return (value*conj(value)).real();
	};

	
	auto intensityP_BBstar = [&](double x){
		comp value = testObs.amplitudes[0].getValue(pow(x,2))(1);
		return (value*conj(value)).real();
	};
	

	
	auto intensityP_BstarBstar = [&](double x){
		comp value = testObs.amplitudes[0].getValue(pow(x,2))(2);
		return (value*conj(value)).real();
	};
	

	testObs.makePlotGraphWithExp("P", "BB", "BottP_BB_Graph_WithExp", intensityP_BB, 10.6322,11.0208);
	testObs.makePlotGraphWithExp("P", "BBstar", "BottP_BBstar_Graph_WithExp", intensityP_BBstar, 10.6322,11.0208);
	testObs.makePlotGraphWithExp("P", "BstarBstar", "BottP_BstarBstar_Graph_WithExp", intensityP_BstarBstar, 10.6322,11.0208);


if(testReader.getFitFlag()){
	//make the minimzer
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2","");
	//Set some criteria for the minimzer to stop
	min->SetMaxFunctionCalls(1000);
	min->SetMaxIterations(100);
	min->SetTolerance(0.01);
	min->SetPrintLevel(1);
	//get the initial parameters and steps from the constructed observable object
	vector<double> fitparams = testObs.getFitParams();
	vector<double> steps = testObs.getStepSizes();
	nParams = fitparams.size();

	//make a function wrapper to minimize the function minfunc (=chisquared)
	ROOT::Math::Functor f(&minfunc,nParams);
	min->SetFunction(f);
	//set the initial conditions and step sizes
	for(int i = 0; i < nParams; i++){
		min->SetVariable(i,to_string(i),fitparams[i],steps[i]);
	}
	//run the minimization
	min->Minimize();
	min->X();

	//extract the resulting fit parameters
	vector<double> finalParams = {};
	for(int i = 0; i < nParams; i ++){
		finalParams.push_back(min->X()[i]);
	}

	//is this necessary?
	testObs.setFitParams(finalParams);
	
	//print out all the fit parameters
	for(double x: testObs.getFitParams()) cout << x <<endl;

	//store the parameters for the minimum that the minimizer found in xs
	//unnecessary?
	const double *xs = min->X();
	//print out the final params
	cout << testObs.amplitudes[0] << endl;
}	

	testReader.writeOutputFile();
	
	//more graphing shit
	testObs.makePlotGraphWithExp("P", "BB", "testBott_BB", intensityP_BB, 10.6322, 11.0208);
	//testObs.makePlotGraphWithExp("P", "BBstar", "testBott_BBstar", intensityP_BBstar, 10.6322, 11.0208);
	//testObs.makePlotGraphWithExp("P", "BstarBstar", "testBott_BstarBstar", intensityP_BstarBstar, 10.6322, 11.0208);


	return 0;
	
}
