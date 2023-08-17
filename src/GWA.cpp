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
#include <TGraph2D.h>
#include <TF2.h>

using namespace std;
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
	minpoles->SetPrintLevel(0);
	
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
	vector<double> grid_Im = linspace(-1, 1, 5);
	vector<double> fitparamspoles = {};
	double steppoles[2] = {0.01,0.01};
	for(int i = 0; i < grid_Re.size(); i++){
		for(int j = 0; j < grid_Im.size(); j++){
			
			fitparamspoles = {grid_Re[i], grid_Im[j]};//cout << grid_Re[i] << endl;
			
			//cout << fitparamspoles[0] << " " << fitparamspoles[1] << endl;
		
			for(int l = 0; l < nParams; l++){
				minpoles->SetVariable(l,to_string(l),fitparamspoles[l],steppoles[l]);
			}
			ps.setPoles(poles);
			//run the minimization
			minpoles->Minimize();
			//extract the resulting fit parameters
			comp finalParams = comp(minpoles->X()[0], minpoles->X()[1]);
			f_val_poles = minpoles->MinValue();
			if(f_val_poles < -10.) poles.push_back(finalParams);
			minpoles->Clear();

		}

	}

	for(int k = 0; k < poles.size(); k++){
		letwrite << poles[k].real() << "	" << poles[k].imag() << "	" << f_val_poles << endl; 
		//cout << poles[k] << endl;
	}

	auto abs_det = [](double* x, double* p){
		//return abs(testObs.amplitudes[0].getDenominator(comp(x[0], x[1])).determinant());
		MatrixXcd temp = (comp(x[0], x[1]) - comp(115,0.5)) * (comp(x[0], x[1]) - comp(118,0.7)) * (comp(x[0], x[1]) - comp(115,-0.5)) * (comp(x[0], x[1]) - comp(118,-0.7)) * MatrixXcd({{1}});
		return abs(temp.determinant());
	};

	auto log_abs_det = [](double* x, double* p){
		//return log(abs(testObs.amplitudes[0].getDenominator(comp(x[0], x[1])).determinant()));
		MatrixXcd temp = (comp(x[0], x[1]) - comp(115,0.5)) * (comp(x[0], x[1]) - comp(118,0.7)) * (comp(x[0], x[1]) - comp(115,-0.5)) * (comp(x[0], x[1]) - comp(118,-0.7)) * MatrixXcd({{1}});
		return log(abs(temp.determinant())); 
	};

	TCanvas c1;
	TGraph2D *gr = new TGraph2D("Data/poles.txt");
	gr->SetMarkerStyle(21);
	gr->Draw("pcol");
	c1.SaveAs("Plots/fit47_32_poles_graph2D.pdf");
	TCanvas c2;
	TGraph *gr2 = new TGraph("Data/poles.txt");
	gr2->SetMarkerStyle(21);
	gr2->Draw("AP");
	c2.SaveAs("Plots/fit47_32_poles_graph.pdf");
	TCanvas c3;
	//TF2 *tf = new TF2("tf", detD, 113, 121, -1, 1, 2);
	//TF2 tf("tf", [](double* x, double* p) { return abs(testObs.amplitudes[0].getDenominator(comp(x[0], x[1])).determinant()); }, 113., 121., -1., 1.);
	TF2 tf("tf", abs_det, 113., 121., -1., 1.);
	tf.Draw("surf1");
	c3.SaveAs("Plots/fit47_32_abs_det.pdf");
	TCanvas c4;
	//TF2 *tf = new TF2("tf", detD, 113, 121, -1, 1, 2);
	TF2 tf2("tf2", log_abs_det, 113., 121., -1., 1.);
	tf2.Draw("surf1");
	c4.SaveAs("Plots/fit47_32_log_abs_det.pdf");

	auto IntCh0 = [&](double x){
		comp value = (testObs.amplitudes[0].getIntegral(pow(x,2), 0));
		return value;
	};

	auto IntCh1 = [&](double x){
		comp value = (testObs.amplitudes[0].getIntegral(pow(x,2), 1));
		return value;
	};

	auto IntCh2 = [&](double x){
		comp value = (testObs.amplitudes[0].getIntegral(pow(x,2), 2));
		return value;
	};

	testObs.plotCompGraph("P", "BB", "fit47_32_DispIntCh0", IntCh0, 10.3, 11.0208);
	testObs.plotCompGraph("P", "BBstar", "fit47_32_DispIntCh1", IntCh1, 10.3, 11.0208);
	testObs.plotCompGraph("P", "BstarBstar", "fit47_32_DispIntCh2", IntCh2, 10.3, 11.0208);

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
