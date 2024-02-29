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
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"
#include "TSystem.h"

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
	return testObs.chisq() + testObs.chisq_with_InclCrossSec();
	
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
	
	filereader formatReader(inputfile);
	formatReader.SetAllCommandLists();
	formatReader.ConstructBareAmps();
	formatReader.setChebys();
	formatReader.setPoles();
	formatReader.setKmats();
	formatReader.loadExpData();
	if(formatReader.getInclCrossSecFlag()) formatReader.loadExpInclCrossSec();
	
	//gets chisq cutoff and weights
	double cutoff = formatReader.getChi2CutOff();

	//saves the observable object outside of filereader object
	testObs = formatReader.getObs();

	//saves original starting parameters
	vector<double> startparams = testObs.getFitParams();
	vector<double> steps = testObs.getStepSizes();
	
	//gets degrees of freedom
	string fittype = "excl";
	int dof = testObs.getNumData()-steps.size();
	
	if(formatReader.getInclCrossSecFlag()){
		dof+=testObs.getNumInclData();
		fittype = "incl";
	}
	
	//Opens up the fit files and sets up TTrees
	string fname = fitsfolder+"run"+to_string(jobnum)+".root";

	TFile *f;
	TTree *t1;
	vector<string> cmds = formatReader.getOutputCmds();
	double chi_squared_incl = 1.0;
	double chi_squared_excl = 1.0;
	vector<string> *tree_cmds = &cmds;

	///

	int namp = testObs.amplitudes.size();

	vector<double> resmasses[namp];

	for(int j = 0; j < namp; j++){

		for(int k = 0; k < testObs.amplitudes[j].getResMasses().size(); k++){

			resmasses[j].push_back(testObs.amplitudes[j].getResMasses()[k]);

		}

	}

	int nchs = testObs.amplitudes[0].getChanNames().size();

	int npole = testObs.amplitudes[0].getResMasses().size();

	double couplings[namp][nchs][npole]; 

	for(int j = 0; j < namp; j++){

		for(int l = 0; l < nchs; l++){
			
			for(int k = 0; k < npole; k++){

				couplings[j][l][k] = testObs.amplitudes[j].getChannels()[l].getCoupling(k);

			}

		}

	}

	int polgrade = testObs.amplitudes[0].getChannels()[0].getChebyCoeffs().size();

	double cheby[namp][nchs][polgrade];

	for(int j = 0; j < namp; j++){

		for(int l = 0; l < nchs; l++){
			
			for(int k = 0; k < polgrade; k++){

				cheby[j][l][k] = testObs.amplitudes[j].getChannels()[l].getChebyCoeffs()[k];

			}

		}

	}

	///
	
	//if the file exists, update the ttree on file. otherwise, make it.
	if(!(gSystem->AccessPathName(fname.c_str(),kFileExists))){
		f=TFile::Open(fname.c_str(),"update");
		t1 = (TTree*)f->Get("fits");
		t1->SetBranchAddress("Commands", &tree_cmds);
		t1->SetBranchAddress("Inclusive_chi", &chi_squared_incl);
		t1->SetBranchAddress("Exclusive_chi", &chi_squared_excl);

		for(int j = 0; j < namp; j++){
			for(int k = 0; k < resmasses[j].size(); k++){
				t1->SetBranchAddress(("mJ" + to_string(j)).c_str(), &resmasses[j].at(k));
			}
		}

		for(int j = 0; j < namp; j++){
			for(int l = 0; l < nchs; l++){
				for(int k = 0; k < npole; k++){
					t1->SetBranchAddress(("gJ" + to_string(j) + "CH" + to_string(l) + "P" + to_string(k)).c_str(), &couplings[j][l][k]);
				}		
			}				
		}

		for(int j = 0; j < namp; j++){
			for(int l = 0; l < nchs; l++){
				for(int k = 0; k < polgrade; k++){
					t1->SetBranchAddress(("aJ" + to_string(j) + "CH" + to_string(l) + "G" + to_string(k)).c_str(), &cheby[j][l][k]);
				}
			}
		}
		
	} else {
		f = TFile::Open(fname.c_str(),"recreate");
		t1 = new TTree("fits", ("Fits from code instance "+to_string(jobnum)).c_str());
		t1->Branch("Commands", &cmds);
		t1->Branch("Inclusive_chi", &chi_squared_incl);
		t1->Branch("Exclusive_chi", &chi_squared_excl);
		
		for(int j = 0; j < namp; j++){
			for(int k = 0; k < resmasses[j].size(); k++){
				t1->Branch(("mJ" + to_string(j)).c_str(), &resmasses[j].at(k));
			}
		}

		for(int j = 0; j < namp; j++){
			for(int l = 0; l < nchs; l++){
				for(int k = 0; k < npole; k++){
					t1->Branch(("gJ" + to_string(j) + "CH" + to_string(l) + "P" + to_string(k)).c_str(), &couplings[j][l][k]);
				}		
			}				
		}

		for(int j = 0; j < namp; j++){
			for(int l = 0; l < nchs; l++){
				for(int k = 0; k < polgrade; k++){
					t1->Branch(("aJ" + to_string(j) + "CH" + to_string(l) + "G" + to_string(k)).c_str(), &cheby[j][l][k]);
				}
			}
		}
	}
	
	//does the fitting
	for (int j = 0; j < numfits; j++){
		//resets the observable
		testObs.setFitParams(startparams);
		formatReader.setObs(testObs);
		auto now = std::chrono::system_clock::now();
		int seed = now.time_since_epoch().count()+jobnum;
		//randomizes the parameters
		if(formatReader.getRandomizeFlag()) formatReader.randomize(seed);
		testObs = formatReader.getObs();
		
	if(formatReader.getFitFlag()){
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
		if(formatReader.getInclCrossSecFlag()) min->SetFunction(g);
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
		
		//if the chisquared is less than the cutoff, add it as a leaf to the ttree
		if(chisq<cutoff){
      testObs.setFitParams(finalParams);
			chi_squared_excl= testObs.chisq()/(testObs.getNumData()-steps.size());
			chi_squared_incl=chisq;
			
			testObs.setFitParams(finalParams);
			formatReader.setObs(testObs);
			cmds = formatReader.getOutputCmds();
			t1->Fill();
		}
	}
	
}
	
	t1->Write(0,TObject::kWriteDelete,0);
	f->Close("R");
  
	return 0;
	
}
