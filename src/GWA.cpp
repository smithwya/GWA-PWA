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
vector<double> f_val_poles = {}; //this variable holds the value of the minization function for poles

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
	string inputfile = (string) argv[1]; //just the core e.g. "fit73-36", without the path nor the extension
	string polefile = (string) argv[2];
	string pWave = (string) argv[3];
	double grid_Re_sx = stod(argv[4]);
	double grid_Re_dx = stod(argv[5]);
	int grid_Re_numpts = atoi(argv[6]);
	double grid_Im_sx = stod(argv[7]);
	double grid_Im_dx = stod(argv[8]);
	int grid_Im_numpts = atoi(argv[9]);
	double func_cutoff = stod(argv[10]);

	//reads the file and creates an observable object with the information from the file
	
	filereader testReader("Data/" + inputfile + ".txt");
	testReader.SetAllCommandLists();
	testReader.ConstructBareAmps();
	testReader.setChebys();
	testReader.setPoles();
	testReader.setKmats();
	testReader.loadExpData();

	//saves the observable object outside of filereader object
	testObs = testReader.getObs(); 

//vector<double> steps = testObs.getStepSizes();

//cout << testObs.chisq() << "	" << testObs.chisq()/(testObs.getNumData() - steps.size() + testObs.getNumInclData()) << endl;

	//initialize the polesearcher object
	ps.settestObs(testObs); 
	ps.setAmpIndex(pWave); 
	int ampindex = ps.getAmpIndex(); 
	
	ofstream letwrite(polefile);
	
	//make the minimizer
	ROOT::Math::Minimizer* minpoles = ROOT::Math::Factory::CreateMinimizer("Minuit2","Simplex");
	//Set some criteria for the minimizer to stop
	minpoles->SetMaxFunctionCalls(10000);
	minpoles->SetMaxIterations(1000);
	//minpoles->SetTolerance(0.001);
	//minpoles->SetPrecision(1.e-7);
	minpoles->SetStrategy(1);
	minpoles->SetPrintLevel(0);
	
	nParams = 2;//real and imaginary part of the pole
	//make a function wrapper to minimize the function minfuncforpoles
	ROOT::Math::Functor f(&minfunc_for_poles,nParams);
	minpoles->SetFunction(f);
	//set the initial conditions and step sizes

	//change these numbers to parameters to pass into from command line
	vector<double> grid_Re = linspace(grid_Re_sx, grid_Re_dx, grid_Re_numpts);
	vector<double> grid_Im = linspace(grid_Im_sx, grid_Im_dx, grid_Im_numpts);
	vector<double> fitparamspoles = {};
	double steppoles[2] = {1.e-2, 1.e-2};
	
	
	for(int i = 0; i < grid_Re.size(); i++){
		for(int j = 0; j < grid_Im.size(); j++){
			
			fitparamspoles = {grid_Re[i], grid_Im[j]};//cout << grid_Re[i] << endl;
			//fitparamspoles = {119.,0.8};
			
			//cout << fitparamspoles[0] << " " << fitparamspoles[1] << endl;
		
			for(int l = 0; l < nParams; l++){
				minpoles->SetVariable(l,to_string(l),fitparamspoles[l],steppoles[l]);
			}

			ps.setPoles(poles);

			//run the minimization
			minpoles->Minimize();
			//minpoles->Minimize();

			//cout << "calls = " << ps.calls_counter_poles << endl;

			//extract the resulting fit parameters
			comp finalParams = comp(minpoles->X()[0], -abs(minpoles->X()[1]));
			long double aux = minpoles->MinValue();
			for (int j=0; j < poles.size(); j++) aux += log(abs((finalParams - poles[j])*(finalParams - conj(poles[j]))));
			if(aux < func_cutoff){
				poles.push_back(finalParams);
				f_val_poles.push_back(aux);
			}
			minpoles->Clear();

		}

	}

	for(int k = 0; k < poles.size(); k++){
		letwrite << poles[k].real() << "	" << poles[k].imag() << "	" << f_val_poles[k] << endl; 
		//cout << poles[k] << endl;
	}

	auto abs_det = [&](double* x, double* p){
		return abs(testObs.amplitudes[ampindex].getDenominator(comp(x[0], x[1])).determinant());
		/*MatrixXcd temp = (comp(x[0], x[1]) - comp(115,0.5)) * (comp(x[0], x[1]) - comp(118,0.7)) * (comp(x[0], x[1]) - comp(115,-0.5)) * (comp(x[0], x[1]) - comp(118,-0.7)) * MatrixXcd({{1}});
		return abs(temp.determinant());*/
	};

	auto log_abs_det = [&](double* x, double* p){
		return log10(abs(testObs.amplitudes[ampindex].getDenominator(comp(x[0], x[1])).determinant()));
		/*MatrixXcd temp = (comp(x[0], x[1]) - comp(115,0.5)) * (comp(x[0], x[1]) - comp(118,0.7)) * (comp(x[0], x[1]) - comp(115,-0.5)) * (comp(x[0], x[1]) - comp(118,-0.7)) * MatrixXcd({{1}});
		return log(abs(temp.determinant()));*/
	};

	testObs.PolePlotGraph2D(inputfile, polefile);

	testObs.PolePlotGraph1D(inputfile, polefile);

	testObs.PoleColormapPlotFunc2D(inputfile, abs_det, "abs_det", grid_Re_sx, grid_Re_dx, grid_Im_sx, grid_Im_dx);

	testObs.PoleColormapPlotFunc2D(inputfile, log_abs_det, "log_abs_det", grid_Re_sx, grid_Re_dx, grid_Im_sx, grid_Im_dx);

	return 0;
	
}
