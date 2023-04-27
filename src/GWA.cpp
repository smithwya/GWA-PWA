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
	
	double sum = 0;
	double std = 0;
	string ampname = "P";
	string channame = "BB";
	double lower_bound = 10632.2;
	double upper_bound = 11020.8;
	double x = 0;
	double y = 0;
	double stat_err = 0;
	double sist_err = 0;
	comp temp = 0;

	int numamp = testObs.getampindex(ampname);
	int numchan = testObs.getchanindex(ampname, channame);
	int totnumofchans = testObs.getData().size();
	int npts = testObs.getData()[totnumofchans * numamp + numchan].amp_expval.size();

	for(int i = 0; i < npts; i++){
		x = testObs.getData()[totnumofchans * numamp + numchan].sqrts[i];
		if(x >= lower_bound && x <= upper_bound){	

			temp = testObs.amplitudes[numamp].getValue(pow(x,2))(numchan);
			y = (temp*conj(temp)).real();
			stat_err = testObs.getData()[totnumofchans * numamp + numchan].amp_expval_stat_err[i];
			sist_err = testObs.getData()[totnumofchans * numamp + numchan].amp_expval_sist_err[i];
			std = sqrt(pow(stat_err, 2) + pow(sist_err, 2));
			sum += pow(((y - testObs.getData()[totnumofchans * numamp + numchan].amp_expval[i])/std), 2);
			//if(i == 56) cout << x << "	" << y << "	" << testObs.getData()[totnumofchans * numamp + numchan].amp_expval[i] << "	" << std << endl;
		
		}
	}

	double tot = sum;

	/*
	sum = 0;
	std = 0;
	
	channame = "KK";
	numchan = testObs.getchanindex(ampname, channame);
	
	x = 0; 
	y = 0;

	temp = 0;
	npts = testObs.getData()[totnumofchans * numamp + numchan].amp_expval.size();

	for(int i = 0; i < npts; i++){
		x = testObs.getData()[totnumofchans * numamp + numchan].sqrts[i];
		if(x >= lower_bound && x <= upper_bound){	
			temp = testObs.amplitudes[numamp].getValue(pow(x,2))(numchan);//this numchan should be wrong because it is no more 0
			y = (temp*conj(temp)).real();
			std = testObs.getData()[totnumofchans * numamp + numchan].amp_expval_stat_err[i];
			sum += pow(((y - testObs.getData()[totnumofchans * numamp + numchan].amp_expval[i])/std), 2);
			//if(i == 56) cout << x << "	" << y << "	" << testObs.getData()[totnumofchans * numamp + numchan].amp_expval[i] << "	" << std << endl;
		}
	}

	tot += sum;
	*/

	//cout<<"param: "<<params[3]<<" chisq "<<tot<<endl;
	return tot;

}


double minfuncforpoles(const double *xx){
	
	int numamp = 0;
	double val = 0;

	MatrixXcd cmat = testObs.amplitudes.at(numamp).getDenominator(comp(xx[0], xx[1]));

	for(int i = 0; i < temppoles.size(); i++){
		cmat = cmat*(comp(xx[0],xx[1])-temppoles[i]);
		cmat = cmat*(comp(xx[0],xx[1])-temppoles[i]);
	}

	if(pow(xx[0], 2) + pow(xx[1], 2) <= 400){
		val = log(abs(cmat.determinant()));
	}
	else{
		val = log(abs(cmat.determinant())) + pow(xx[0], 2) + pow(xx[1], 2);
	}

	return val;

}


int main()
{

	auto intensityP_BB = [&](double x){
		comp value = testObs.amplitudes[0].getValue(pow(x,2))(0);
		return (value*conj(value)).real();
	};


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

	testObs.makePlotGraphWithExp("P", "BB", "BottP_BB_Graph_WithExp", intensityP_BB, 10632.2,11020.8);

	/*
	TRandom3 gen(testReader.getSeed());
	vector<double> randParams(nParams,0);
	for(int i = 0; i < nParams; i++){
		randParams[i] = gen.Uniform(-3000, 3000);
	}
	*/

	/*
	JPsi testJ = JPsi();
	
	function<double(double*)> jminf = [&](double *x){
		return x[0]*x[0];
		//return testJ.JPsiminfunc({x[0]});
	};
	*/



/*for(double x: testObs.getData()[0].sqrts) cout << x << endl; //the first 4 num are broken
cout << testObs.getData().size() << endl;
for(expchan k: testObs.getData()) cout << k.wavename << " " << k.channame << endl;
expchan k = testObs.getData()[0];
for(int i = 0; i < k.amp_expval.size(); i++){
	cout << k.sqrts[i] << " " << k.amp_expval[i] << " " << k.amp_expval_stat_err[i] << endl;
} 
return 0;*/


	//make the minimzer
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2","");
	//Set some criteria for the minimzer to stop
	min->SetMaxFunctionCalls(1000000);
	min->SetMaxIterations(10000);
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
		//min->VariableLimits(i, fitparams[i] - 100, fitparams[i] + 100); // to do with the uncertanties of input file
		//if(i != 3) min->FixVariable(i); 
	}
	
	min->Minimize();
	min->X();

	vector<double> finalParams = {};
	for(int i = 0; i < nParams; i ++){
		finalParams.push_back(min->X()[i]);
	}


	testObs.setFitParams(finalParams);
	
	//cout<<testObs.amplitudes[0]<<endl;

	//store the parameters for the minimum that the minimizer found in xs
	const double *xs = min->X();
	
	testObs.makePlotGraphWithExp("P", "BB", "testBott_BB", intensityP_BB, 10632.2, 11020.8);

	//note to self: need to get rid of 'dumbJ' in amplitude.cpp later when doing non-radJPsi fits

	/*
	//search poles:
	ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
   	minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
   	minimum->SetMaxIterations(10000);  // for GSL
   	minimum->SetTolerance(0.001);
   	minimum->SetPrintLevel(1);
	ROOT::Math::Functor g(&minfuncforpoles, 2);
   	double step[2] = {0.01,0.01};
    double variable[2] = { -1.,1.2};
	//TRandom3 r(3);
    TRandom3 r(testReader.getSeed()); //idk why it doesn't recognize "stoi()" in filereader.cpp/GetSeed
    variable[0] = r.Uniform(-2000,2000);
    variable[1] = r.Uniform(-2000,2000);
	minimum->SetFunction(g);
   	minimum->SetVariable(0,"x",variable[0], step[0]);
   	minimum->SetVariable(1,"y",variable[1], step[1]);
	minimum->SetVariableLimits(0, -20, 20);
	minimum->SetVariableLimits(1, -20, 20);
	
	//temppoles = {comp(0.00631603,0),comp(3.18353,0)};
	for(int i = 0; i <1; i++){
		minimum->Minimize();
		const double *pole = minimum->X();
		temppoles.push_back({pole[0], pole[1]});
	}
	
	for(comp x : temppoles) cout<<x<<endl;
	*/
/*
	expchan tchan = testObs.getData()[totnumofchans * numamp + numchan];
	cout<<tchan.wavename<<" "<<tchan.channame<<endl;
	for(int i = 0; i < tchan.sqrts.size(); i++){
		cout<<tchan.sqrts[i]<<" "<<tchan.amp_expval[i]<<endl;
	}
	*/

	return 0;
}
