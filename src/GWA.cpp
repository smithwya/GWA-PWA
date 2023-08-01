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

double minfunc_with_InclCrossSec(const double *xx){

	vector<double> params = {};
	for(int i = 0; i < nParams; i++){
		params.push_back(xx[i]);
	}

	testObs.setFitParams(params);
	return testObs.chisq_with_InclCrossSec(params);
}

int main(int argc, char ** argv)
{

	int jobnum = atoi(argv[1]);
	int fitnum = atoi(argv[2]);
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
	//selects a seed based off clock + job number
	int seed = std::chrono::system_clock::now().time_since_epoch().count()+jobnum+fitnum;
	testReader.setSeed(seed);
	if(testReader.getRandomizeFlag()) testReader.randomize(seed); 

	//gets chisq cutoff
	double cutoff = testReader.getChi2CutOff();

	//saves the observable object outside of filereader object
	testObs = testReader.getObs();

	////tests and plots

	testReader.writeMathematicaOutputFile("Data/Math_test2.dat");

	double lower_bound = testObs.amplitudes[0].getFitInterval()[0];
	double upper_bound = testObs.amplitudes[0].getFitInterval()[1];
	//testObs.plotInclCrossSec("InclCrossSec", lower_bound, upper_bound);
	//testObs.plotInclCrossSecVsSumOfExcl("Diff", lower_bound, upper_bound);

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

	/*auto intensityP_B_sstarB_sstar = [&](double x){
		comp value = testObs.amplitudes[0].getValue(pow(x,2))(3);
		return (value*conj(value)).real();
	};*/

	auto ImagPartInt = [&](double x){
		double value = (testObs.amplitudes[0].getIntegral(pow(x,2), 0)).imag();
		return value;
	};

	auto RealPartInt = [&](double x){
		double value = (testObs.amplitudes[0].getIntegral(pow(x,2), 0)).real();
		return value;
	};

	auto AlternImagPartInt = [&](double x){
		double value = (testObs.amplitudes[0].getIntegral(pow(x,2), 0)).imag();
		return value;
	};

	auto AlternRealPartInt = [&](double x){
		double value = (testObs.amplitudes[0].getIntegral(pow(x,2), 0)).real();
		return value;
	};

	comp test = testObs.amplitudes[0].getValue(pow(10.87,2))(0);

	cout << "test " << (test*conj(test)).real() << endl;

	testObs.makePlotGraph("P", "BB", "test2_ImagPartInt", ImagPartInt, 10.6322, 11.0208);
	testObs.makePlotGraph("P", "BB", "test2_AlternImagPartInt", AlternImagPartInt, 10.6322, 11.0208);
	testObs.makePlotGraph("P", "BB", "test2_RealPartInt", RealPartInt, 10.6322, 11.0208);
	testObs.makePlotGraph("P", "BB", "test2_AlternRealPartInt", AlternRealPartInt, 10.6322, 11.0208);
	testObs.makePlotGraphWithExp("P", "BB", "test2_BB", intensityP_BB, 10.6322,11.0208);
	testObs.makePlotGraphWithExp("P", "BBstar", "test2_BBstar", intensityP_BBstar, 10.6322,11.0208);
	testObs.makePlotGraphWithExp("P", "BstarBstar", "test2_BstarBstar", intensityP_BstarBstar, 10.6322,11.0208);
	//testObs.makePlotGraphWithExp("P", "B_sstarB_sstar", "BottB_sstarB_sstar_Graph_WithExp", intensityP_B_sstarB_sstar, 10.6322,11.0208);

	return 0;

	////

	//saves original starting parameters
	vector<double> startparams = testObs.getFitParams();
	vector<double> steps = testObs.getStepSizes();
	
	//gets degrees of freedom
	int dof = testObs.getNumData()-steps.size();

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
			string fname = fitsfolder+"fit"+to_string(jobnum)+"-"+to_string(fitnum);
			testObs.setFitParams(finalParams);
			testReader.setObs(testObs);
			testReader.writeOutputFile(fname);
			ofstream outputfile(fname,ios::app);
			outputfile<<"chisq = "<<chisq<<endl;
			outputfile.close();
		}
	}	

	return 0;
	
}
