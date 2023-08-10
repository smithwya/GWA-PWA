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
#include "polesearcher.h"

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
polesearcher ps;
int nParams = 0;
vector<comp> temppoles = {}; //this vector hold the poles found by polesearcher

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
	return testObs.chisq_with_InclCrossSec();
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

	//make the minimzer
	ps.settestObs(testObs);
	ROOT::Math::Minimizer* minpoles = ROOT::Math::Factory::CreateMinimizer("Minuit2","");
	//Set some criteria for the minimizer to stop
	minpoles->SetMaxFunctionCalls(100000);
	minpoles->SetMaxIterations(10000);
	minpoles->SetTolerance(0.001);
	minpoles->SetPrintLevel(1);
	
	nParams = 2;//real and imaginary part of the pole
	//make a function wrapper to minimize the function minfuncforpoles
	ROOT::Math::Functor f(&minfunc_for_poles,nParams);
	minpoles->SetFunction(f);
	//set the initial conditions and step sizes
	string ampname = "P";
	ps.setAmpIndex(ampname);
	int ampindex = ps.getAmpIndex();
	
	TRandom3 gen;
	ofstream letwrite("Data/poles.txt");
	vector<double> grid_Re = linspace(113, 121, 5);
	vector<double> grid_Im = linspace(-3, 3, 5);
	vector<double> fitparamspoles = {};
	double steppoles[2] = {0.01,0.01};
	for(int i = 0; i < grid_Re.size(); i++){
		for(int j = 0; j < grid_Im.size(); j++){
			
			//vector<double> fitparamspoles = {testObs.amplitudes[ampindex].getResMasses()[0]+1, 1};
			//vector<double> fitparamspoles = {grid_Re[i], grid_Im[j]};
			fitparamspoles = {grid_Re[i], grid_Im[j]};//cout << grid_Re[i] << endl;
			
			//for {118,0} it finds the same pole two time
			//for {121,0} it finds the \Upsilon(11020) maybe
			cout << fitparamspoles[0] << " " << fitparamspoles[1] << endl;
			for(int counter = 0; counter < testReader.getAddPoleList().size(); counter++){
		
				for(int l = 0; l < nParams; l++){
					minpoles->SetVariable(l,to_string(l),fitparamspoles[l],steppoles[l]);
				}
				ps.setTemppoles(temppoles);
				//run the minimization
				minpoles->Minimize();
				//extract the resulting fit parameters
				comp finalParams = comp(minpoles->X()[0], minpoles->X()[1]);
				temppoles.push_back(finalParams);
				minpoles->Clear();

			}

			for(int k = 0; k < temppoles.size(); k++){
				letwrite << temppoles[k].real() << "	" << temppoles[k].imag() << endl; 
				//cout << temppoles[k] << endl;
			}
		}
		//minpoles->Clear();
	}

	TCanvas c1;
	TGraph *gr = new TGraph("Data/poles.txt");
	gr->SetMarkerStyle(21);
	gr->Draw("AP");
	c1.SaveAs("Plots/poles_graph.pdf");

	//testReader.writeMathematicaOutputFile("Data/Math_test2.dat");
	
	/*
	double lower_bound = testObs.amplitudes[0].getFitInterval()[0];
	double upper_bound = testObs.amplitudes[0].getFitInterval()[1];
	testObs.plotInclCrossSec("InclCrossSec", lower_bound, upper_bound);
	testObs.plotInclCrossSecVsSumOfExcl("Diff", lower_bound, upper_bound);

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

	auto intensityP_B_sstarB_sstar = [&](double x){
		comp value = testObs.amplitudes[0].getValue(pow(x,2))(3);
		return (value*conj(value)).real();
	};

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

	cout << "test " << testObs.amplitudes[0].getValue(pow(10.7,2)) << endl;

	testObs.makePlotGraph("P", "BB", "test2_ImagPartInt", ImagPartInt, 10.6322, 11.0208);
	testObs.makePlotGraph("P", "BB", "test2_AlternImagPartInt", AlternImagPartInt, 10.6322, 11.0208);
	testObs.makePlotGraph("P", "BB", "test2_RealPartInt", RealPartInt, 10.6322, 11.0208);
	testObs.makePlotGraph("P", "BB", "test2_AlternRealPartInt", AlternRealPartInt, 10.6322, 11.0208);
	testObs.makePlotGraphWithExp("P", "BB", "test2_BB", intensityP_BB, 10.6322,11.0208);
	testObs.makePlotGraphWithExp("P", "BBstar", "test2_BBstar", intensityP_BBstar, 10.6322,11.0208);
	testObs.makePlotGraphWithExp("P", "BstarBstar", "test2_BstarBstar", intensityP_BstarBstar, 10.6322,11.0208);
	testObs.makePlotGraphWithExp("P", "B_sstarB_sstar", "BottB_sstarB_sstar_Graph_WithExp", intensityP_B_sstarB_sstar, 10.6322,11.0208);
	*/

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
