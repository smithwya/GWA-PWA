#include <iostream>
#include <complex>
#include <vector>
#include <Eigen/Dense>
#include <sstream>
#include <fstream>
#include <map>
#include <chrono>
#include <ctime>
#include "amplitude.h"
#include "channel.h"
#include "filereader.h"
#include "observable.h"
#include "bottomonium.h"
#include "polesearcher.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom3.h"
#include "TError.h"
#include <TGraph2D.h>
#include <TF2.h>

using namespace std;
typedef std::chrono::system_clock Clock;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;
observable testObs = observable();
polesearcher ps;
int nParams = 0;
vector<comp> poles = {}; //this vector holds the poles found by polesearcher
double f_val_poles = 0; //this variable holds the value of the minization function for poles

double minfunc(const double *xx){

	vector<double> params = {};
	for(int i = 0; i < nParams; i++){
		params.push_back(xx[i]);
	}

	testObs.setFitParams(params);
	return testObs.chisq();
}

double minfunc_with_InclCrossSec(const double *xx){

	vector<double> params = {};
	for(int i = 0; i < nParams; i++){
		params.push_back(xx[i]);
	}

	testObs.setFitParams(params);
	return testObs.chisq()+testObs.chisq_with_InclCrossSec();
	
}

double minfunc_for_poles(const double *xx){

	vector<double> params = {};
	for(int i = 0; i < nParams; i++){
		params.push_back(xx[i]);
	}
	
	return ps.minfuncforpoles(params);
}

template<typename T>
vector<double> linspace(T start_in, T end_in, int num_in)
{

  	std::vector<double> linspaced;

  	double start = static_cast<double>(start_in);
  	double end = static_cast<double>(end_in);
  	double num = static_cast<double>(num_in);

  	if (num == 0) { return linspaced; }
  	if (num == 1) 
  	{
  	  linspaced.push_back(start);
  	  return linspaced;
  	}

  	double delta = (end - start) / (num - 1);

  	for(int i=0; i < num-1; ++i)
  	{
  	  linspaced.push_back(start + delta * i);
  	}
  	linspaced.push_back(end); // I want to ensure that start and end
  	                          // are exactly the same as the input
  	return linspaced;

}

int main(int argc, char ** argv)
{
	int jobnum = atoi(argv[1]);
	int numfits = atoi(argv[2]);
	string inputfile = (string) argv[3];
	string fitsfolder = (string) argv[4];

	//reads the file and creates an observable object with the information from the file
	
	filereader testReader(inputfile);
	testReader.SetAllCommandLists();
	testReader.ConstructBareAmps();
	testReader.setChebys();
	testReader.setPoles();
	testReader.setKmats();
	testReader.loadExpData();
	if(testReader.getInclCrossSecFlag()) testReader.loadExpInclCrossSec();
	
	//saves current time
	std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
	time_t tt = std::chrono::system_clock::to_time_t(now);
	tm local_tm = *localtime(&tt);
	
	std::stringstream timebuffer;
	auto t = std::chrono::system_clock::now();
	timebuffer << local_tm.tm_year + 1900 << '-'<< local_tm.tm_mon + 1 << '-'<< local_tm.tm_mday;
	
	//selects a seed based off clock + job number
	int seed = now.time_since_epoch().count()+jobnum+numfits;
	testReader.setSeed(seed); 

	//gets chisq cutoff
	double cutoff = testReader.getChi2CutOff();

	//saves the observable object outside of filereader object
	testObs = testReader.getObs();

	//saves original starting parameters
	//vector<double> startparams = testObs.getFitParams();
	vector<double> steps = testObs.getStepSizes();
	
	//gets degrees of freedom
	string fittype = "excl";
	int dof = testObs.getNumData()-steps.size();
	
	if(testReader.getInclCrossSecFlag()){
		dof+=testObs.getNumInclData();
		fittype = "incl";
	}

	for (int j = 0; j < numfits; j++){
		
		if(testReader.getRandomizeFlag()) testReader.randomize(seed);
		testObs = testReader.getObs();
		
		
	if(testReader.getFitFlag()){
		//make the minimzer
		ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2","");
		//Set some criteria for the minimzer to stop
		min->SetMaxFunctionCalls(100000);
		min->SetMaxIterations(10000);
		min->SetTolerance(0.001);
		min->SetPrintLevel(1);
		//get the initial parameters and steps from the constructed observable object
		vector<double> fitparams = testObs.getFitParams();
		nParams = fitparams.size();
		//make a function wrapper to minimize the function minfunc (=chisquared)
		ROOT::Math::Functor f(&minfunc,nParams);
		ROOT::Math::Functor g(&minfunc_with_InclCrossSec,nParams);
		if(testReader.getInclCrossSecFlag()) min->SetFunction(g);
		else min->SetFunction(f);
		//set the initial conditions and step sizes
		for(int i = 0; i < nParams; i++){
			min->SetVariable(i,to_string(i),fitparams[i],steps[i]);
		}
		//run the minimization
		min->Minimize();
		//extract the resulting fit parameters
		vector<double> finalParams = {};
		for(int i = 0; i < nParams; i ++){
			finalParams.push_back(min->X()[i]);
		}

		double chisq = min->MinValue()/dof;

		if(chisq<cutoff){
			testObs.setFitParams(finalParams);
			double excl_chisq = testObs.chisq()/(testObs.getNumData()-steps.size());
			string fname = fitsfolder+timebuffer.str()+"-"+to_string(jobnum)+"-"+to_string(j)+"-"+fittype+"-"+to_string(chisq)+"-"+to_string(excl_chisq);
			testReader.setObs(testObs);
			testReader.writeOutputFile(fname);
			ofstream outputfile(fname,ios::app);
			outputfile<<"chisq = "<<chisq<<endl;
			outputfile<<"excl_chisq = "<<excl_chisq<<endl;
			outputfile.close();
		}
		
	}	
	
}
	return 0;
	
}
