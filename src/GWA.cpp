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
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::seconds;

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
	return testObs.excl_chisq();
}

double minfunc_with_InclCrossSec(const double *xx){

	vector<double> params = {};
	for(int i = 0; i < nParams; i++){
		params.push_back(xx[i]);
	}

	testObs.setFitParams(params);
	return testObs.excl_chisq() + testObs.incl_chisq();
	
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
	string sheet = (string) argv[11]; // NB: "false" = physical sheet; "true" = first unphysical sheet;
	int jobnum = atoi(argv[12]);
	int numfits = atoi(argv[13]);
	//string inputfile = (string) argv[14];
	//string fitsfolder = (string) argv[15];
	string fitsfolder = "Fits/";
	string polesfolder = "Poles/";
	cout<<"[GWA] Initializing GWA..."<<endl;

	auto t_start = high_resolution_clock::now();
	cout<<"[IO] Reading inputfile...";
	filereader formatReader(string("Data/") + inputfile + string(".txt"));

	auto t_end = high_resolution_clock::now();
	auto delta_t = duration_cast<seconds>(t_end-t_start);
	cout<<"done ("<<delta_t.count()<<"s)\n";

	cout<<"[GWA] Initializing objects...";
	t_start = high_resolution_clock::now();
	
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
	
	t_end = high_resolution_clock::now();
        delta_t = duration_cast<seconds>(t_end-t_start);
        cout<<"done ("<<delta_t.count()<<"s)\n";

	
	//if you want just plotting: add a flag:
	if(formatReader.getPlotFlag()){
		cout<<"[PLOT] Plotting..."<<endl;
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
			comp value = testObs.amplitudes[0].getValue(pow(x,2))(4);
			return (value*conj(value)).real();
		};

		testObs.makePlotGraphWithExp("P", "BB", inputfile+"_BB", intensityP_BB, 10.6322,11.0208);
		testObs.makePlotGraphWithExp("P", "BBstar", inputfile+"_BBstar", intensityP_BBstar, 10.6322,11.0208);
		testObs.makePlotGraphWithExp("P", "BstarBstar", inputfile+"_BstarBstar", intensityP_BstarBstar, 10.6322,11.0208);
		testObs.makePlotGraphWithExp("P", "B_sstarB_sstar", inputfile+"_B_sstarB_sstar", intensityP_B_sstarB_sstar, 10.6322,11.0208);
		testObs.makePlotGraphDummy(inputfile+"_Dummy", intensityP_Dummy, 10.6322,11.0208);
		if(formatReader.getInclCrossSecFlag()){
			testObs.plotInclCrossSecVsSumOfExcl(inputfile+"_InclCrossSecVsSumOfExcl", 10.6322,11.2062);
			testObs.plotInclCrossSecWithExp(inputfile+"_InclCrossSecWithExp", 10.6322,11.2062);
		}

		//testObs.makePlotGraph_ExpOnly("P", "BB", inputfile+"_BB_justpts", 10.6322,  11.0208, "#sqrt{s} (GeV)", "#sigma (pb)");
		//testObs.makePlotGraph_ExpOnly("P", "BBstar", inputfile+"_BBstar_justpts", 10.6322,  11.0208, "#sqrt{s} (GeV)", "#sigma (pb)");
		//testObs.makePlotGraph_ExpOnly("P", "BstarBstar", inputfile+"_BstarBstar_justpts", 10.6322,  11.0208, "#sqrt{s} (GeV)", "#sigma (pb)");
		//testObs.makePlotGraph_ExpOnly("P", "B_sstarB_sstar", inputfile+"_B_sstarB_sstar_justpts", 10.84,  11.0208, "#sqrt{s} (GeV)", "#sigma (pb)");

	}
	

	if(formatReader.getFitFlag()){
		cout<<"[GWA] Fit flag detected"<<endl;
		//saves original starting parameters
		vector<double> startparams = testObs.getFitParams();
		vector<double> steps = testObs.getStepSizes();
		//gets degrees of freedom
		string fittype = "excl";
		int dof = testObs.getNumData()-steps.size();		
		cout<<"[FIT] Using inclusive data...";
		if(formatReader.getInclCrossSecFlag()){
			cout<<"YES\n";
			dof+=testObs.getNumInclData();
			fittype = "incl";
		}else{
		cout<<"NO\n";
		}
		
		cout<<"[FIT] Degrees of freedom: "<<dof<<endl;
		
		t_start = high_resolution_clock::now();
		cout<<"[IO] Setting up TTrees...\n";
		//Opens up the fit files and sets up TTrees
		string fname = fitsfolder+"run"+to_string(jobnum)+".root";

		TFile *file;
		TTree *t1;
		vector<string> cmds = formatReader.getOutputCmds();
		double red_chi_squared = 1.0;
		double red_chi_squared_incl = 1.0;
		double red_chi_squared_excl = 1.0;
		vector<string> *tree_cmds = &cmds;

		///

		int namp = testObs.amplitudes.size();

		vector<double> resmasses[namp];

		int nchs = testObs.amplitudes[0].getChanNames().size();//0th element because the number of channels is the same for each wave

		vector<double> couplings[namp][nchs]; 

		vector<double> cheby[namp][nchs];

		vector<double> kmats[namp][nchs][nchs];

		int npolej = 0;
		int polgradej = 0;
		int gradej = 0;
		cout<<"[IO] Checking if .root file exists...";
		//if the file exists, update the ttree on file. otherwise, make it.
		if(!(gSystem->AccessPathName(fname.c_str(),kFileExists))){
			
			cout<<"No. Generating .root file...\n"<<endl;
			file=TFile::Open(fname.c_str(),"update");
			t1 = (TTree*)file->Get("fits");
			t1->SetBranchAddress("Commands", &tree_cmds);
			t1->SetBranchAddress("Red_chi", &red_chi_squared);
			
			cout<<"Setting up TTree branches...\n"<<endl;
			if(formatReader.getInclCrossSecFlag()){
				t1->SetBranchAddress("Inclusive_red_chi", &red_chi_squared_incl);
				t1->SetBranchAddress("Exclusive_red_chi", &red_chi_squared_excl);
			}

			for(int j = 0; j < namp; j++){
				npolej = testObs.amplitudes[j].getResMasses().size();
				resmasses[j].resize(npolej);
				for(int k = 0; k < npolej; k++){
					t1->SetBranchAddress(("mJ" + to_string(j) + "_num" + to_string(k)).c_str(), &resmasses[j].at(k));
				}
			}

			for(int j = 0; j < namp; j++){
				npolej = testObs.amplitudes[j].getResMasses().size();
				for(int l = 0; l < nchs; l++){
					couplings[j][l].resize(npolej);
					for(int k = 0; k < npolej; k++){
						t1->SetBranchAddress(("gJ" + to_string(j) + "CH" + to_string(l) + "P" + to_string(k)).c_str(), &couplings[j][l].at(k));
					}		
				}				
			}

			for(int j = 0; j < namp; j++){
				for(int l = 0; l < nchs; l++){
					polgradej = testObs.amplitudes[j].getChannels()[l].getChebyCoeffs().size();
					cheby[j][l].resize(polgradej);
					for(int k = 0; k < polgradej; k++){
						t1->SetBranchAddress(("aJ" + to_string(j) + "CH" + to_string(l) + "G" + to_string(k)).c_str(), &cheby[j][l].at(k));
					}
				}
			}

			for(int j = 0; j < namp; j++){
				gradej = testObs.amplitudes[j].getkParameters().size();
				for(int l = 0; l < nchs; l++){
					for(int m = l; m < nchs; m++){
						kmats[j][l][m].resize(gradej);
						for(int k = 0; k < gradej; k++){

							t1->SetBranchAddress(("kJ" + to_string(j) + "pow" + to_string(k) + "ch" + to_string(l) + to_string(m)).c_str(), &kmats[j][l][m].at(k));
						
						}
					}
				}
			}	
			
		} else {
			cout<<"Yes. Loading in TTree from .root file..."<<endl;
			file = TFile::Open(fname.c_str(),"recreate");
			t1 = new TTree("fits", ("Fits from code instance "+to_string(jobnum)).c_str());
			t1->Branch("Commands", &cmds);
			t1->Branch("Red_chi", &red_chi_squared);

			if(formatReader.getInclCrossSecFlag()){
				t1->Branch("Inclusive_red_chi", &red_chi_squared_incl);
				t1->Branch("Exclusive_red_chi", &red_chi_squared_excl);
			}
			
			for(int j = 0; j < namp; j++){
				npolej = testObs.amplitudes[j].getResMasses().size();
				resmasses[j].resize(npolej);
				for(int k = 0; k < npolej; k++){
					t1->Branch(("mJ" + to_string(j) + "_num" + to_string(k)).c_str(), &resmasses[j].at(k));//cout << "hello" << endl;
				}
			}

			for(int j = 0; j < namp; j++){
				npolej = testObs.amplitudes[j].getResMasses().size();
				for(int l = 0; l < nchs; l++){
					couplings[j][l].resize(npolej); // Ensure the vector has npolej elements			
					for(int k = 0; k < npolej; k++){
						//cout << couplings[j][l].at(k) << endl;
						t1->Branch(("gJ" + to_string(j) + "CH" + to_string(l) + "P" + to_string(k)).c_str(), &couplings[j][l].at(k));
					}		
				}				
			}

			for(int j = 0; j < namp; j++){
				for(int l = 0; l < nchs; l++){
					polgradej = testObs.amplitudes[j].getChannels()[l].getChebyCoeffs().size(); //cout << polgradej << endl; 
					cheby[j][l].resize(polgradej);
					for(int k = 0; k < polgradej; k++){
						t1->Branch(("aJ" + to_string(j) + "CH" + to_string(l) + "G" + to_string(k)).c_str(), &cheby[j][l].at(k));
					}
				}
			}

			for(int j = 0; j < namp; j++){
				int gradej = testObs.amplitudes[j].getkParameters().size();
				for(int l = 0; l < nchs; l++){
					for(int m = l; m < nchs; m++){
						kmats[j][l][m].resize(gradej);
						for(int k = 0; k < gradej; k++){

							t1->Branch(("kJ" + to_string(j) + "pow" + to_string(k) + "ch" + to_string(l) + to_string(m)).c_str(), &kmats[j][l][m].at(k));
						
						}
					}
				}
			}

		}

		cout<<"[IO] TTrees connected ";
		t_end = high_resolution_clock::now();
		delta_t = duration_cast<seconds>(t_end-t_start);
		cout<<"("<<delta_t.count()<<"s)\n";

		cout<<"[FIT] Getting fit sequence"<<endl;
		vector<string> fitseq = formatReader.readFitSequence(formatReader.getFitSequence());
		
		//does the fitting
		cout<<"[FIT] Fitting "<<numfits<<" times"<<endl;
		for (int j = 0; j < numfits; j++){
			cout<<endl<<endl<<"[FIT] Beginning fit #"<<j+1<<endl;
			t_start = high_resolution_clock::now();
			cout<<"[FIT] Resetting parameters to initial values...";
			//resets the observable
			testObs.setFitParams(startparams);
			formatReader.setObs(testObs);
			auto now = std::chrono::system_clock::now();
			int seed = now.time_since_epoch().count()+jobnum;
			//randomizes the parameters
			if(formatReader.getRandomizeFlag()){
			cout<<"randomizing starting parameters..."<<endl;
			formatReader.randomize(seed);
			}
			testObs = formatReader.getObs();
			
			t_end = high_resolution_clock::now();
			delta_t = duration_cast<seconds>(t_end-t_start);
			cout<<"done ("<<delta_t.count()<<"s)\n";
			//make the minimzer
			ROOT::Math::Minimizer* min[fitseq.size()];

			for(int l = 0; l < fitseq.size(); l++){
				min[l] = ROOT::Math::Factory::CreateMinimizer("Minuit2",fitseq[l]);
			}
			cout<<"[FIT] Wrapping chisquared function for MINUIT...";
			t_start = high_resolution_clock::now();

			//get the initial parameters and steps from the constructed observable object
			vector<double> fitparams = testObs.getFitParams();
			nParams = fitparams.size();
			//make a function wrapper to minimize the function minfunc (=chisquared)
			ROOT::Math::Functor f(&minfunc,nParams);
			ROOT::Math::Functor g(&minfunc_with_InclCrossSec,nParams);

			double chisq[fitseq.size()];
			int numKmatbgcoeffs = formatReader.getKmatList().size() * testObs.numChans * (testObs.numChans + 1)/2.;
			int numremainingparamsperamp = nParams/testObs.getNumAmps() - numKmatbgcoeffs;
			int numresmasses = formatReader.getAddPoleList().size();
			t_end = high_resolution_clock::now();
			delta_t = duration_cast<seconds>(t_end-t_start);
			cout<<"done ("<<delta_t.count()<<"s)\n";
			for(int l = 0; l < sizeof(min)/sizeof(ROOT::Math::Minimizer*); l++){

				//Set some criteria for the minimzer to stop
				cout<<"[FIT] Initializing MINUIT...";
				t_start = high_resolution_clock::now();

				min[l]->SetMaxFunctionCalls(10000);
				min[l]->SetMaxIterations(1000);
				min[l]->SetTolerance(0.1);
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
				t_end = high_resolution_clock::now();
				delta_t = duration_cast<seconds>(t_end-t_start);
				cout<<"done ("<<delta_t.count()<<"s)\n";
				//run the minimization
				cout<<"[FIT] Minimizing..."<<endl;
				t_start = high_resolution_clock::now();
				min[l]->Minimize();
				t_end = high_resolution_clock::now();
				delta_t = duration_cast<seconds>(t_end-t_start);
				cout<<"done ("<<fitseq[l]<<": "<<delta_t.count()<<"s)\n";
				//extract the resulting fit parameters
				for(int i = 0; i < nParams; i++){
					fitparams[i] = min[l]->X()[i];
					//cout << min[l]->X()[i] << endl;
				}

				chisq[l] = min[l]->MinValue()/dof; //cout << chisq[l] << endl << endl;

			}

			vector<double> finalParams = {};

			for(int i = 0; i < nParams; i++){
				finalParams.push_back(fitparams[i]);
			}		

			//for(double x: finalParams) cout << x << endl;	

			//if the chisquared is less than the cutoff, add it as a leaf to the ttree
			if(chisq[fitseq.size() - 1]<cutoff){
				cout<<"[FIT] Reduced chisq is less than cutoff, saving fit...";
				t_start = high_resolution_clock::now();
				testObs.setFitParams(finalParams);

				for(int j = 0; j < namp; j++){
					npolej = testObs.amplitudes[j].getResMasses().size();
					resmasses[j].resize(npolej);
					for(int k = 0; k < npolej; k++){
					
						resmasses[j].push_back(testObs.amplitudes[j].getResMasses()[k]);

					}

				}

				for(int j = 0; j < namp; j++){
					npolej = testObs.amplitudes[j].getResMasses().size();
					for(int l = 0; l < nchs; l++){
						couplings[j][l].resize(npolej); // Ensure the vector has npolej elements
						for(int k = 0; k < npolej; k++){
						
							couplings[j][l].push_back(testObs.amplitudes[j].getChannels()[l].getCoupling(k));
			
						}
			
					}
			
				}

				for(int j = 0; j < namp; j++){

					for(int l = 0; l < nchs; l++){
						polgradej = testObs.amplitudes[j].getChannels()[l].getChebyCoeffs().size();
						cheby[j][l].resize(polgradej);
						for(int k = 0; k < polgradej; k++){
						
							cheby[j][l].push_back(testObs.amplitudes[j].getChannels()[l].getChebyCoeffs()[k]);

						}

					}

				}

				for(int j = 0; j < namp; j++){
					gradej = testObs.amplitudes[j].getkParameters().size();
					for(int l = 0; l < nchs; l++){
						for(int m = l; m < nchs; m++){
							kmats[j][l][m].resize(gradej);
							for(int k = 0; k < gradej; k++){
							
								kmats[j][l][m].push_back(testObs.amplitudes[j].getkParameters()[k](l,m).real());//real because the kmatrix bg coeffs are real
								if(l != m) kmats[j][m][l].push_back(kmats[j][l][m].at(k));

							}
						}
					}
				}

				red_chi_squared=chisq[fitseq.size() - 1]; //cout << "Prova: " << chisq[1] << "VS " << (testObs.excl_chisq() + testObs.incl_chisq())/(testObs.getNumInclData() + testObs.getNumData()-steps.size()) << endl;

				if(formatReader.getInclCrossSecFlag()){
					red_chi_squared_excl= testObs.excl_chisq()/(testObs.getNumData()-steps.size());
					red_chi_squared_incl= testObs.incl_chisq()/(testObs.getNumInclData()-steps.size());
				}
				
				formatReader.setObs(testObs);
				cmds = formatReader.getOutputCmds();
				t1->Fill();
				
				t_end = high_resolution_clock::now();
				delta_t = duration_cast<seconds>(t_end-t_start);
				cout<<"done ("<<delta_t.count()<<"s)\n";
			}
		
		}
		cout<<"[IO] Writing fits to file...";
		t_start = high_resolution_clock::now();

		t1->Write(0,TObject::kWriteDelete,0);		
		file->Close("R");
		t_end = high_resolution_clock::now();
		delta_t = duration_cast<seconds>(t_end-t_start);
		cout<<"done ("<<delta_t.count()<<"s)\n";
	}		


	//if you want to add a polesearch flag:
	if(formatReader.getPolesearchFlag()){

		//initialize the polesearcher object
		ps.settestObs(testObs); 
		ps.setAmpIndex(pWave); 
		int ampindex = ps.getAmpIndex(); 
		ps.SetSheet(false);
		if(sheet == "true") ps.SetSheet(true);
		
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
				
				fitparamspoles = {grid_Re[i], grid_Im[j]};
			
				for(int l = 0; l < nParams; l++){
					minpoles->SetVariable(l,to_string(l),fitparamspoles[l],steppoles[l]);
				}

				ps.setPoles(poles);

				//run the minimization
				minpoles->Minimize();

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

		for(comp x : poles) temp_Re.push_back(x.real());
		for(comp x : poles) temp_Im.push_back(x.imag());

		for(int k = 0; k < poles.size(); k++){
			letwrite << poles[k].real() << "	" << poles[k].imag() << "	" << f_val_poles[k] << endl; 
			//cout << temp_Re[k] << temp_Im[k] << endl;
			//cout << poles[k] << endl;
		}

		//////

		vector<double> ampres_mod[testObs.numChans];
		vector<double> ampres_phase[testObs.numChans];

		vector<double> denomres_mod[testObs.numChans][testObs.numChans];
		vector<double> denomres_phase[testObs.numChans][testObs.numChans];

		MatrixXcd phsp = MatrixXcd::Identity(testObs.numChans,testObs.numChans);

		for(int k = 0; k < testObs.numChans; k++){	
			
			for(int i = 0; i < poles.size(); i++){

				comp x = comp(poles[i].real(),poles[i].imag());

				amplitude tempamp = testObs.amplitudes[testObs.getampindex(pWave)]; 

				int J = tempamp.getJ(); 

				double eps = tempamp.getEpsilon(); 

				for(int i = 0; i < testObs.numChans; i++){
					phsp(i,i)= pow(tempamp.getMomentum(i,x),J+0.5)/pow(x,.25);
				}

				bool sh = true;
				//VectorXcd tempvec = eps * tempamp.getValueForPoles(x, false); // order O(eps)
				VectorXcd tempvec = 0.25 * eps * (tempamp.getValueForPoles(x + eps, sh) - tempamp.getValueForPoles(x - eps, sh) + comp(0.,1.) * tempamp.getValueForPoles(x + comp(0.,1.) * eps, sh) - comp(0.,1.) * tempamp.getValueForPoles(x - comp(0.,1.) * eps, sh)); // order O(eps^2)

				if(abs(poles[i].imag()) < 2 * tempamp.getEpsilon()){
					tempvec = 0.5 * eps * (tempamp.getValueForPoles(x + eps, sh) - tempamp.getValueForPoles(x - eps, sh));
				}
				
				ampres_mod[k].push_back(abs(tempvec(k))); 
				ampres_phase[k].push_back(arg(tempvec(k))); 

				//MatrixXcd tempmat = eps * tempamp.getDenominator(x, false).inverse();
				MatrixXcd tempmat = 0.25 * eps * (tempamp.getDenominator(x + eps, sh).inverse() - tempamp.getDenominator(x - eps, sh).inverse() + comp(0.,1.) * tempamp.getDenominator(x + comp(0.,1.) * eps, sh).inverse() - comp(0.,1.) * tempamp.getDenominator(x - comp(0.,1.) * eps, sh).inverse()); 

				if(abs(poles[i].imag()) < 2 * tempamp.getEpsilon()){
					tempmat = 0.5 * eps * (tempamp.getDenominator(x + eps, sh).inverse() - tempamp.getDenominator(x - eps, sh).inverse());
				}

				for(int j = 0; j < testObs.numChans; j++){
					denomres_mod[k][j].push_back(abs(tempmat(k,j)));
					denomres_phase[k][j].push_back(arg(tempmat(k,j)));
				}

			}

		}	

		//for(VectorXcd h: ampresidues) cout << h << endl << endl;
		//for(MatrixXcd h: denomresidues) cout << h << endl << endl;

		//////

		string fname_poles = polesfolder+inputfile+"_poles.root";

		TFile *file_poles;
		TTree *t1_poles; 

		vector<double> *tree_poles_Re = &temp_Re;
		vector<double> *tree_poles_Im = &temp_Im;
		vector<double> *tree_f_val_poles = &f_val_poles;
		vector<double> *tree_ampres_mod[testObs.numChans];
		vector<double> *tree_ampres_phase[testObs.numChans];
		vector<double> *tree_denomres_mod[testObs.numChans][testObs.numChans];
		vector<double> *tree_denomres_phase[testObs.numChans][testObs.numChans];

		for(int i = 0; i < testObs.numChans; i++){
			tree_ampres_mod[i] = &ampres_mod[i];
			tree_ampres_phase[i] = &ampres_phase[i];
			for(int j = 0; j < testObs.numChans; j++){
				tree_denomres_mod[i][j] = &denomres_mod[i][j];
				tree_denomres_phase[i][j] = &denomres_phase[i][j];
			}
		}

		//if the file exists, update the ttree on file. otherwise, make it.
		if(!(gSystem->AccessPathName(fname_poles.c_str(),kFileExists))){
			//cout << "sono qui! 3" << endl;
			file_poles=TFile::Open(fname_poles.c_str(),"update");
			t1_poles = (TTree*)file_poles->Get("poles");
			t1_poles->SetBranchAddress("Poles_Re", &tree_poles_Re);
			t1_poles->SetBranchAddress("Poles_Im", &tree_poles_Im);
			t1_poles->SetBranchAddress("FVAL_Poles", &tree_f_val_poles);
			for(int i = 0; i < testObs.numChans; i++){
				t1_poles->SetBranchAddress(("amp_res_ch" + to_string(i) + "_mod").c_str(), &tree_ampres_mod[i]);
				t1_poles->SetBranchAddress(("amp_res_ch" + to_string(i) + "_phase").c_str(), &tree_ampres_phase[i]);
				for(int j = 0; j < testObs.numChans; j++){
					t1_poles->SetBranchAddress(("denom_res_ch" + to_string(i) + to_string(j) + "_mod").c_str(), &tree_denomres_mod[i][j]);
					t1_poles->SetBranchAddress(("denom_res_ch" + to_string(i) + to_string(j) + "_phase").c_str(), &tree_denomres_phase[i][j]);				
				}
			}
			
		} else {
			//cout << "sono qui! 4" << endl;
			file_poles = TFile::Open(fname_poles.c_str(),"recreate");
			t1_poles = new TTree("poles", ("Poles from code instance " + inputfile).c_str());
			t1_poles->Branch("Poles_Re", &temp_Re);
			t1_poles->Branch("Poles_Im", &temp_Im);
			t1_poles->Branch("FVAL_Poles", &f_val_poles);
			for(int i = 0; i < testObs.numChans; i++){
				t1_poles->Branch(("amp_res_ch" + to_string(i) + "_mod").c_str(), &ampres_mod[i]);
				t1_poles->Branch(("amp_res_ch" + to_string(i) + "_phase").c_str(), &ampres_phase[i]);
				for(int j = 0; j < testObs.numChans; j++){
					t1_poles->Branch(("denom_res_ch" + to_string(i) + to_string(j) + "_mod").c_str(), &denomres_mod[i][j]);
					t1_poles->Branch(("denom_res_ch" + to_string(i) + to_string(j) + "_phase").c_str(), &denomres_phase[i][j]);				
				}
			}

		}

		//////

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

		} 

		t1_poles->Fill();
		t1_poles->Write(0,TObject::kWriteDelete,0);
		file_poles->Close("R");

	}

	return 0;
	
}
