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

double minfunc(vector<double> xx){
	return testObs.chisq(xx);
}

double minfunc2(vector<double> xx){

	vector<double> params = {};
	for(int i = 0; i < nParams; i++){
		params.push_back(xx[i]);
	}

	testObs.setFitParams(params);

	//for(string ampname: testObs.getAmpNames()){
	string ampname = "P";
		
		int amp_index = testObs.getampindex(ampname);
		
		double lower_bound = sqrt(testObs.amplitudes[amp_index].getFitInterval()[0]);
		double upper_bound = sqrt(testObs.amplitudes[amp_index].getFitInterval()[1]);

		int totnumofchans = testObs.getData().size();

		double tot = 0;

		//for(string channame: testObs.amplitudes[amp_index].getChanNames()){
		for(string channame: {"BB"}){
			
			double sum = 0;
			double std = 0;
			double x = 0;
			double y = 0;
			double stat_err = 0;
			double sist_err = 0;
			comp temp = 0;
			
			int chan_index = testObs.getchanindex(ampname ,channame);
			int npts = testObs.getData()[totnumofchans * amp_index + chan_index].amp_expval.size();

			for(int i = 0; i < npts; i++){
				x = testObs.getData()[totnumofchans * amp_index + chan_index].sqrts[i];
				if(x >= lower_bound && x <= upper_bound){	

					temp = testObs.amplitudes[amp_index].getValue(pow(x,2))(chan_index);
					y = (temp*conj(temp)).real();
					stat_err = testObs.getData()[totnumofchans * amp_index + chan_index].amp_expval_stat_err[i];
					sist_err = testObs.getData()[totnumofchans * amp_index + chan_index].amp_expval_sist_err[i];
					std = sqrt(pow(stat_err, 2) + pow(sist_err, 2));
					sum += pow(((y - testObs.getData()[totnumofchans * amp_index + chan_index].amp_expval[i])/std), 2);
					//if(i == 56) cout << x << "	" << y << "	" << testObs.getData()[totnumofchans * numamp + numchan].amp_expval[i] << "	" << std << endl;
		
				}
			}

			tot += sum;
			
		}

	return tot;

}



int main()
{

	//reads the file and creates an observable object with the information from the file
	filereader testReader("Data/testdat.txt");
	testReader.SetAllCommandLists();
	testReader.ConstructBareAmps();
	testReader.setChebys();
	testReader.setPoles();
	testReader.setKmats();
	testReader.loadExpData();


	//saves the observable object outside of filereader object
	testObs = testReader.getObs();

	//print out the observable starting params
	cout << testObs.amplitudes[0] << endl; 
	vector<double> params = testObs.getFitParams();
	for(double x : params) cout<<x<<endl;

	//cout<<"original minfunc: "<< minfunc2(params)<<endl;
	cout<<"my minfunc: "<< minfunc(params)<<endl;
	return 0;
	/*
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

	return 0;
	*/
}
