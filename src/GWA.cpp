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
vector<double> temp_Re = {};
vector<double> temp_Im = {};
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
	string inputfile = (string) argv[1]; //TO CHANGE: for polesearch is just the core e.g. "fit73-36", without the path nor the extension
	string polefile = (string) argv[2];
	string pWave = (string) argv[3];
	double grid_Re_sx = stod(argv[4]);
	double grid_Re_dx = stod(argv[5]);
	int grid_Re_numpts = atoi(argv[6]);
	double grid_Im_sx = stod(argv[7]);
	double grid_Im_dx = stod(argv[8]);
	int grid_Im_numpts = atoi(argv[9]);
	double func_cutoff = stod(argv[10]);
	string sheet = (string) argv[11]; // NB: "false" = physical sheet; "true" = first unphysical sheet;
	int jobnum = atoi(argv[12]);
	int numfits = atoi(argv[13]);
	//string inputfile = (string) argv[14];
	//string fitsfolder = (string) argv[15];
	string fitsfolder = "Fits/";

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

	/*
	
	//if you want just plotting: add a flag:
	
	double lower_bound = testObs.amplitudes[0].getFitInterval()[0];
	double upper_bound = testObs.amplitudes[0].getFitInterval()[1];

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
	
	auto intensityP_Dummy = [&](double x){
		comp value = testObs.amplitudes[0].getValue(pow(x,2))(2);
		return (value*conj(value)).real();
	};

	testObs.makePlotGraphWithExp("P", "BB", plotname+"_BB", intensityP_BB, 10.6322,11.0208);
	testObs.makePlotGraphWithExp("P", "BBstar", plotname+"_BBstar", intensityP_BBstar, 10.6322,11.0208);
	testObs.makePlotGraphWithExp("P", "BstarBstar", plotname+"_BstarBstar", intensityP_BstarBstar, 10.6322,11.0208);
	testObs.makePlotGraphWithExp("P", "B_sstarB_sstar", plotname+"_B_sstarB_sstar", intensityP_B_sstarB_sstar, 10.6322,11.0208);
	testObs.makePlotGraphDummy(plotname+"_Dummy", intensityP_Dummy, 10.6322,11.0208);
	testObs.plotInclCrossSecVsSumOfExcl(plotname+"_InclCrossSecVsSumOfExcl", 10.6322,11.2062);
	testObs.plotInclCrossSecWithExp(plotname+"_InclCrossSecWithExp", 10.6322,11.2062);

	//testObs.makePlotGraph_ExpOnly("P", "BB", plotname+"_BB_justpts", 10.6322,  11.0208, "#sqrt{s} (GeV)", "#sigma (pb)");
	//testObs.makePlotGraph_ExpOnly("P", "BBstar", plotname+"_BBstar_justpts", 10.6322,  11.0208, "#sqrt{s} (GeV)", "#sigma (pb)");
	//testObs.makePlotGraph_ExpOnly("P", "BstarBstar", plotname+"_BstarBstar_justpts", 10.6322,  11.0208, "#sqrt{s} (GeV)", "#sigma (pb)");
	//testObs.makePlotGraph_ExpOnly("P", "B_sstarB_sstar", plotname+"_B_sstarB_sstar_justpts", 10.84,  11.0208, "#sqrt{s} (GeV)", "#sigma (pb)");

	vector<double> steps = testObs.getStepSizes();

	//cout << testObs.chisq_with_InclCrossSec() << "	" << testObs.chisq_with_InclCrossSec()/(testObs.getNumData() - steps.size() + testObs.getNumInclData()) << endl;

	//cout << testObs.chisq_with_InclCrossSec() - testObs.chisq() << endl;

	//cout << (testObs.chisq_with_InclCrossSec() + testObs.chisq()) / (testObs.getNumData() - steps.size() + testObs.getNumInclData()) << endl;

	
	
	//double totabs = 1.7 * (testObs.getNumData() - steps.size() + testObs.getNumInclData());

	//double exclabs = 24.6 * (testObs.getNumData() - steps.size());

	//double ratio = exclabs/totabs;

	//cout << totabs << "	" << exclabs << "	" << ratio << endl;

	

	cout << testObs.chisq() << endl;

	cout << testObs.chisq() / (testObs.getNumData() - steps.size()) << endl;

	cout << testObs.chisq_with_InclCrossSec() << endl;

	cout << testObs.chisq_with_InclCrossSec() / (testObs.getNumInclData() - steps.size()) << endl;

	cout << (15*testObs.chisq() + testObs.chisq_with_InclCrossSec())/(testObs.getNumData() + testObs.getNumInclData() - steps.size()) << endl; 
	
	return 0;
	
	*/

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
	//for (int j = 0; j < numfits; j++){
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
		//ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2","Simplex");
		ROOT::Math::Minimizer* min[2];
		min[0] = ROOT::Math::Factory::CreateMinimizer("Minuit2","Simplex");
		min[1] = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");

		//get the initial parameters and steps from the constructed observable object
		vector<double> fitparams = testObs.getFitParams();
		nParams = fitparams.size();
		//make a function wrapper to minimize the function minfunc (=chisquared)
		ROOT::Math::Functor f(&minfunc,nParams);
		ROOT::Math::Functor g(&minfunc_with_InclCrossSec,nParams);

		double chisq = 0;

		int numKmatbgcoeffs = formatReader.getKmatList().size() * testObs.numChans * (testObs.numChans + 1)/2.;
		int numremainingparamsperamp = nParams/testObs.getNumAmps() - numKmatbgcoeffs;
		int numresmasses = formatReader.getAddPoleList().size();

		for(int l = 0; l < sizeof(min)/sizeof(ROOT::Math::Minimizer*); l++){

			//Set some criteria for the minimzer to stop
			min[l]->SetMaxFunctionCalls(10);
			min[l]->SetMaxIterations(10);
			min[l]->SetTolerance(1);
			min[l]->SetPrintLevel(1);

			if(formatReader.getInclCrossSecFlag()) min[l]->SetFunction(g);
			else min[l]->SetFunction(f);

			//set the initial conditions and step sizes
			for(int i = 0; i < nParams; i++){
				min[l]->SetVariable(i,to_string(i),fitparams[i],steps[i]);
			} 

			for(int i = 0; i < testObs.getNumAmps(); i++){
				if(testObs.amplitudes[i].getKMatType() == "kmat-CDD"){
					for(int k = numremainingparamsperamp - numresmasses; k < nParams - numresmasses; k++){
						min[l]->SetVariableLowerLimit(k + i * nParams, 0.);
					}
				}
			}

			//run the minimization
			min[l]->Minimize();

			if(min[l]->Status() == 0 /*|| min[l]->IsValidError() == false*/){
				cout << "The fit is not valid:" << endl;
				cout << min[l]->Status() << endl;
				//cout << min[l]->IsValidError() << endl;
				return 0;
			}

			//extract the resulting fit parameters
			for(int i = 0; i < nParams; i++){
				fitparams[i] = min[l]->X()[i];
				//cout << min[l]->X()[i] << endl;
			}

			chisq = min[l]->MinValue()/dof; cout << chisq << endl <<endl;

		}

		vector<double> finalParams = {};

		for(int i = 0; i < nParams; i++){
			finalParams.push_back(fitparams[i]);
		}

		for(int i = 0; i < testObs.getNumAmps(); i++){
			if(testObs.amplitudes[i].getKMatType() == "kmat-CDD"){
				for(int k = numremainingparamsperamp - numresmasses; k < nParams - numresmasses; k++){
					//min[l]->SetVariableLowerLimit(k + i * nParams, 0.);
					cout << finalParams[k + i * nParams] << endl;
				}
			}
		}		

		//exit(0);

		//if the chisquared is less than the cutoff, add it as a leaf to the ttree
		if(chisq<cutoff){
      		testObs.setFitParams(finalParams);
			chi_squared_excl= testObs.chisq()/(testObs.getNumData()-steps.size());
			chi_squared_incl=chisq;
			formatReader.setObs(testObs);
			cmds = formatReader.getOutputCmds();
			t1->Fill();
		}
	//}
	
}


/*

//if you want to add a polesearch flag:

//initialize the polesearcher object
	ps.settestObs(testObs); 
	ps.setAmpIndex(pWave); 
	int ampindex = ps.getAmpIndex(); 
	ps.SetSheet(false);
	if(sheet == "true") ps.SetSheet(true);
	
	ofstream letwrite(polefile);


	TFile *file;
	TTree *t1; 
	string fname = "test.root";

	vector<double> *tree_poles_Re = &temp_Re;
	vector<double> *tree_poles_Im = &temp_Im;
	vector<double> *tree_f_val_poles = &f_val_poles;

	//if the file exists, update the ttree on file. otherwise, make it.
	if(!(gSystem->AccessPathName(fname.c_str(),kFileExists))){
		file=TFile::Open(fname.c_str(),"update");
		t1 = (TTree*)file->Get("fits");
		t1->SetBranchAddress("Poles_Re", &tree_poles_Re);
		t1->SetBranchAddress("Poles_Im", &tree_poles_Im);
		t1->SetBranchAddress("FVAL_Poles", &tree_f_val_poles);
		
	} else {
		file = TFile::Open(fname.c_str(),"recreate");
		t1 = new TTree("fits", ("Fits from code instance "+to_string(jobnum)).c_str());
		t1->Branch("Poles_Re", &tree_poles_Re);
		t1->Branch("Poles_Im", &tree_poles_Im);
		t1->Branch("FVAL_Poles", &tree_f_val_poles);

	}

	
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

	auto abs_det = [&](double* x, double* p){
		return abs(testObs.amplitudes[ampindex].getDenominator(comp(x[0], x[1]),false).determinant());
		//MatrixXcd temp = (comp(x[0], x[1]) - comp(115,0.5)) * (comp(x[0], x[1]) - comp(118,0.7)) * (comp(x[0], x[1]) - comp(115,-0.5)) * (comp(x[0], x[1]) - comp(118,-0.7)) * MatrixXcd({{1}});
		//return abs(temp.determinant());
	};

	auto log_abs_det = [&](double* x, double* p){
		return log10(abs(testObs.amplitudes[ampindex].getDenominator(comp(x[0], x[1]),false).determinant()));
		//MatrixXcd temp = (comp(x[0], x[1]) - comp(115,0.5)) * (comp(x[0], x[1]) - comp(118,0.7)) * (comp(x[0], x[1]) - comp(115,-0.5)) * (comp(x[0], x[1]) - comp(118,-0.7)) * MatrixXcd({{1}});
		//return log(abs(temp.determinant()));
	};

	auto abs_det_II = [&](double* x, double* p){
		return abs(testObs.amplitudes[ampindex].getDenominator(comp(x[0], x[1]),true).determinant());
		//MatrixXcd temp = (comp(x[0], x[1]) - comp(115,0.5)) * (comp(x[0], x[1]) - comp(118,0.7)) * (comp(x[0], x[1]) - comp(115,-0.5)) * (comp(x[0], x[1]) - comp(118,-0.7)) * MatrixXcd({{1}});
		//return abs(temp.determinant());
	};

	auto log_abs_det_II = [&](double* x, double* p){
		return log10(abs(testObs.amplitudes[ampindex].getDenominator(comp(x[0], x[1]),true).determinant()));
		//MatrixXcd temp = (comp(x[0], x[1]) - comp(115,0.5)) * (comp(x[0], x[1]) - comp(118,0.7)) * (comp(x[0], x[1]) - comp(115,-0.5)) * (comp(x[0], x[1]) - comp(118,-0.7)) * MatrixXcd({{1}});
		//return log(abs(temp.determinant()));
	};

	if(ps.GetSheet()){

		testObs.PolePlotGraph2D(inputfile, polefile, sheet, grid_Re_sx, grid_Re_dx, grid_Im_sx, grid_Im_dx);

		testObs.PolePlotGraph1D(inputfile, polefile, sheet, grid_Re_sx, grid_Re_dx, grid_Im_sx, grid_Im_dx);

		testObs.PoleColormapPlotFunc2D(inputfile, abs_det_II, "abs_det_II", grid_Re_sx, grid_Re_dx, grid_Im_sx, grid_Im_dx);

		testObs.PoleColormapPlotFunc2D(inputfile, log_abs_det_II, "log_abs_det_II", grid_Re_sx, grid_Re_dx, grid_Im_sx, grid_Im_dx);

	}else{

		testObs.PolePlotGraph2D(inputfile, polefile, sheet, grid_Re_sx, grid_Re_dx, grid_Im_sx, grid_Im_dx);

		testObs.PolePlotGraph1D(inputfile, polefile, sheet, grid_Re_sx, grid_Re_dx, grid_Im_sx, grid_Im_dx);

		testObs.PoleColormapPlotFunc2D(inputfile, abs_det, "abs_det", grid_Re_sx, grid_Re_dx, grid_Im_sx, grid_Im_dx);

		testObs.PoleColormapPlotFunc2D(inputfile, log_abs_det, "log_abs_det", grid_Re_sx, grid_Re_dx, grid_Im_sx, grid_Im_dx);

	} */

	
	t1->Write(0,TObject::kWriteDelete,0);
	f->Close("R");
	return 0;
	
}
