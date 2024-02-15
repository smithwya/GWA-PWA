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
//polesearch
  
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
	string sheet = (string) argv[11]; // NB: "false" = physical sheet; "true" = first unphysical sheet;
	int jobnum = atoi(argv[12]);
	int numfits = atoi(argv[13]);
	//string inputfile = (string) argv[14];
	//string fitsfolder = (string) argv[15];


	//reads the file and creates an observable object with the information from the file
	
	filereader testReader("Data/" + inputfile + ".txt");
	testReader.SetAllCommandLists();
	testReader.ConstructBareAmps();
	testReader.setChebys();
	testReader.setPoles();
	testReader.setKmats();
	testReader.loadExpData();

	if(testReader.getInclCrossSecFlag()) testReader.loadExpInclCrossSec();
/*	
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
*/
	//saves the observable object outside of filereader object
	testObs = testReader.getObs(); 

//vector<double> steps = testObs.getStepSizes();

//cout << testObs.chisq() << "	" << testObs.chisq()/(testObs.getNumData() - steps.size() + testObs.getNumInclData()) << endl;

	//initialize the polesearcher object
	ps.settestObs(testObs); 
	ps.setAmpIndex(pWave); 
	int ampindex = ps.getAmpIndex(); 
	ps.SetSheet(false);
	if(sheet == "true") ps.SetSheet(true);
	
	ofstream letwrite(polefile);

	///////

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

	///////
	
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

	/*cout << testObs.amplitudes[0].getChannels()[0].getMomentum(comp(118,0.5)) << endl;
	cout << testObs.amplitudes[0].getChannels()[0].getComplexMomentum(comp(118,0.5)) << endl << endl; 

	cout << testObs.amplitudes[0].getRhoN(comp(118,0.5),0,false) << endl;
	cout << testObs.amplitudes[0].getRhoN(comp(118,0.5),0,true) << endl << endl;*/

	//cout << testObs.amplitudes[0].getIntegrand(118.,comp(118,0.5),0) << endl << endl; 

	/*cout << testObs.amplitudes[0].getIntegral(comp(118.,.5),0,false) << endl; 
	cout << testObs.amplitudes[0].getIntegral(comp(118.,.5),0,true) << endl << endl; */

	//testReader.writeMathematicaOutputFile("Data/Math_fit47-32.txt");exit(0);

	auto imagpartintcomplex = [&](double x){
		return testObs.amplitudes[0].getIntegral(comp(x*x,pow(10,-2)),0,false).imag();
	};

	auto imagpartintcomplexII = [&](double x){
		return testObs.amplitudes[0].getIntegral(comp(x*x,-pow(10,-2)),0,true).imag();
	};

	auto repartintcomplex = [&](double x){
		return testObs.amplitudes[0].getIntegral(comp(x*x,pow(10,-2)),0,false).real();
	};

	auto repartintcomplexII = [&](double x){
		return testObs.amplitudes[0].getIntegral(comp(x*x,-pow(10,-2)),0,true).real();
	};

	/*testObs.makePlotGraph("P", "BB", "imagpartintcomplex", imagpartintcomplex, 10.6322, 11.0208);
	testObs.makePlotGraph("P", "BB", "imagpartintcomplexII", imagpartintcomplexII, 10.6322, 11.0208);
	testObs.makePlotGraph("P", "BB", "repartintcomplex", repartintcomplex, 10.6322, 11.0208);
	testObs.makePlotGraph("P", "BB", "repartintcomplexII", repartintcomplexII, 10.6322, 11.0208);*/

	//cout << testObs.amplitudes[0].getKMatrix(comp(118,0.5)) << endl;

	/*cout << testObs.amplitudes[0].getDenominator(comp(118,0.5),false) << endl << endl;
	cout << testObs.amplitudes[0].getDenominator(comp(118,0.5),true) << endl; exit(0);*/
	
	/*cout << testObs.amplitudes[0].getDenominator(comp(118,5.5),false) << endl << endl;
	cout << testObs.amplitudes[0].getDenominator(comp(118,5.5),true) << endl; exit(0);*/

	//cout << log10(abs(testObs.amplitudes[0].getDenominator(comp(113.676, -0.00807683),true).determinant())) << endl; exit(0);
	
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

	//////

	for(comp x : poles) temp_Re.push_back(x.real());
	for(comp x : poles) temp_Im.push_back(x.imag());

	//////

	for(int k = 0; k < poles.size(); k++){
		letwrite << poles[k].real() << "	" << poles[k].imag() << "	" << f_val_poles[k] << endl; 
		cout << temp_Re[k] << temp_Im[k] << endl;
		//cout << poles[k] << endl;
	}

	/*divide re and imag part, but also try with complex in tree, add f_val_poles*/

	/////
	t1->Fill();
	t1->Write(0,TObject::kWriteDelete,0);
	file->Close("R");
	/////

	auto abs_det = [&](double* x, double* p){
		return abs(testObs.amplitudes[ampindex].getDenominator(comp(x[0], x[1]),false).determinant());
		/*MatrixXcd temp = (comp(x[0], x[1]) - comp(115,0.5)) * (comp(x[0], x[1]) - comp(118,0.7)) * (comp(x[0], x[1]) - comp(115,-0.5)) * (comp(x[0], x[1]) - comp(118,-0.7)) * MatrixXcd({{1}});
		return abs(temp.determinant());*/
	};

	auto log_abs_det = [&](double* x, double* p){
		return log10(abs(testObs.amplitudes[ampindex].getDenominator(comp(x[0], x[1]),false).determinant()));
		/*MatrixXcd temp = (comp(x[0], x[1]) - comp(115,0.5)) * (comp(x[0], x[1]) - comp(118,0.7)) * (comp(x[0], x[1]) - comp(115,-0.5)) * (comp(x[0], x[1]) - comp(118,-0.7)) * MatrixXcd({{1}});
		return log(abs(temp.determinant()));*/
	};

	auto abs_det_II = [&](double* x, double* p){
		return abs(testObs.amplitudes[ampindex].getDenominator(comp(x[0], x[1]),true).determinant());
		/*MatrixXcd temp = (comp(x[0], x[1]) - comp(115,0.5)) * (comp(x[0], x[1]) - comp(118,0.7)) * (comp(x[0], x[1]) - comp(115,-0.5)) * (comp(x[0], x[1]) - comp(118,-0.7)) * MatrixXcd({{1}});
		return abs(temp.determinant());*/
	};

	auto log_abs_det_II = [&](double* x, double* p){
		return log10(abs(testObs.amplitudes[ampindex].getDenominator(comp(x[0], x[1]),true).determinant()));
		/*MatrixXcd temp = (comp(x[0], x[1]) - comp(115,0.5)) * (comp(x[0], x[1]) - comp(118,0.7)) * (comp(x[0], x[1]) - comp(115,-0.5)) * (comp(x[0], x[1]) - comp(118,-0.7)) * MatrixXcd({{1}});
		return log(abs(temp.determinant()));*/
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

	}

/*
	cout << "test " << testObs.amplitudes[0].getValue(pow(10.7,2)) << endl;

	testObs.makePlotGraph("P", "BB", "test2_ImagPartInt", ImagPartInt, 10.6322, 11.0208);
	testObs.makePlotGraph("P", "BB", "test2_AlternImagPartInt", AlternImagPartInt, 10.6322, 11.0208);
	testObs.makePlotGraph("P", "BB", "test2_RealPartInt", RealPartInt, 10.6322, 11.0208);
	testObs.makePlotGraph("P", "BB", "test2_AlternRealPartInt", AlternRealPartInt, 10.6322, 11.0208);
	testObs.makePlotGraphWithExp("P", "BB", "test2_BB", intensityP_BB, 10.6322,11.0208);
	testObs.makePlotGraphWithExp("P", "BBstar", "test2_BBstar", intensityP_BBstar, 10.6322,11.0208);
	testObs.makePlotGraphWithExp("P", "BstarBstar", "test2_BstarBstar", intensityP_BstarBstar, 10.6322,11.0208);
	testObs.makePlotGraphWithExp("P", "B_sstarB_sstar", "BottB_sstarB_sstar_Graph_WithExp", intensityP_B_sstarB_sstar, 10.6322,11.0208);
	


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

		if(chisq<cutoff){
			double excl_chisq = testObs.chisq()/(testObs.getNumData()-steps.size());
			string fname = fitsfolder+timebuffer.str()+"-"+to_string(jobnum)+"-"+to_string(j)+"-"+fittype+"-"+to_string(chisq)+"-"+to_string(excl_chisq);
			testObs.setFitParams(finalParams);
			testReader.setObs(testObs);
			testReader.writeOutputFile(fname);
			ofstream outputfile(fname,ios::app);
			outputfile<<"chisq = "<<chisq<<endl;
			outputfile<<"excl_chisq = "<<excl_chisq<<endl;
			outputfile.close();
		}
		
		testReader.setObs(testObs);
	}	
	
}
*/
	return 0;
	
}
