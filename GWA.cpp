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
#include "observable.h"
#include "GWAFunction.h"
#include "filereader.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "GWAFcn.h"
#include "GWAFunction.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"


using namespace std;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;




//plots a function which takes in a single double and returns a double
void makePlot(string pdfname, function<double(double)> func){
	double lower_bound = 1.0;
	double upper_bound = 2.5;
	int num_bins = 300;
	double delta = (upper_bound - lower_bound)/num_bins;

	TH1D *plotter = new TH1D(pdfname.c_str(),pdfname.c_str(), num_bins, lower_bound, upper_bound);
	//plotter->SetMinimum(-60.0);
	//plotter->SetMaximum(90.0);
	double sqrtS = 1.0;

	for(int i = 0; i < num_bins; i++){
		sqrtS = lower_bound + delta/2 + i * delta; //centroid
		double val = func(sqrtS);
		if(isnan(val)) plotter->SetBinContent(i + 1,  0);
		else plotter->SetBinContent(i + 1,  val);
	}


	TFile file("pdf_folder.root", "recreate");
	TCanvas canv;
	plotter->Write();
	plotter->Draw();
	canv.SaveAs(("Plots/"+pdfname+".pdf").c_str());
	file.Close();
	return;
}

void makeTable(string filename, function<double(double)> func){
	ofstream outfile("Tables/"+filename+".txt");
	double lower_bound = 1.0;
	double upper_bound = 2.5;
	int num_bins = 300;
	double delta = (upper_bound - lower_bound)/num_bins;
	double sqrtS = 1.0;

	for(int i = 0; i < num_bins; i++){
		sqrtS = lower_bound + delta/2 + i * delta; //centroid
		outfile<< sqrtS<<" "<<func(sqrtS) << endl;
	}
	outfile.close();
	return;
}



//plots a function which takes in a double and returns a complex number
void plotComp(string pdfname,function<comp(double)> func){

	auto realFunc = [&](double x){
        return func(x).real();
	};
	auto imagFunc = [&](double x){
        return func(x).imag();
	};

	makePlot("Re "+pdfname, realFunc);
	makePlot("Im "+pdfname, imagFunc);

}

vector<amplitude> readInput(string filename){

//do something

return {};
}



int main()
{

	double smin = pow(0.997,2);
	double smax = pow(2.5,2);

	channel c1 = channel("PiPi",{0.13498, 0.13498});
	channel c2 = channel("KK", {0.49761, 0.49761});
	channel c3 = channel("RhoRho", {0.762, 0.762});

	amplitude S_wave = amplitude(0,0.6,smin,smax,{c1,c2,c3});

	S_wave.setChebyCoeffs("PiPi",3, 1.0, {-2173.73, 3272.05, -1553.73, 361.79});
	S_wave.setChebyCoeffs("KK",3, 1.0, {1927.5, -2996.84, 1354.33 , -287.158});
	S_wave.setChebyCoeffs("RhoRho", 3,1.0, {0});

	S_wave.addPole(0.00631603, {"PiPi","KK","RhoRho"},{-6.05845,-7.24665,1.95572});
	S_wave.addPole(3.18353,{"PiPi","KK","RhoRho"},{-0.235874,-1.09189,-0.232508});
	S_wave.addPole(7.54762,{"PiPi","KK","RhoRho"},{8.29518,12.9082,-1.04908});

	S_wave.setKParams(0,{{16.1903, 18.3081 ,0}, {19.2388, 0}, {-15.811}});
	S_wave.setKParams(1, {{-6.27393,-9.19686,0}, {-13.3706,0}, {6.03082}});
	
	amplitude D_wave = amplitude(2,0.6,smin,smax,{c1,c2,c3});


	D_wave.setChebyCoeffs("PiPi", 3, 1.0, {109.759,-2.8444, 71.095,0 });
	D_wave.setChebyCoeffs("KK", 3, 1.0 , {-26.964, -135.324, -14.887,0});
	D_wave.setChebyCoeffs("RhoRho",3,1.0, {0});

	D_wave.addPole( 2.332, {"PiPi","KK","RhoRho"}, {0.155, 1.027, -0.344});
	D_wave.addPole( 1.667, {"PiPi","KK","RhoRho"}, {0.908, -0.250, 0});
	D_wave.addPole( 2.894, {"PiPi","KK","RhoRho"}, {0.852, 0.698, 1.134});

	D_wave.setKParams(0, {{-1.170,-1.012, 0 }, {0.395, 0}, {-11.365}});
	D_wave.setKParams(1, {{0.211, 0.168, 0}, {0.005,0}, {3.923}});
		
	
	int ch = 0;
	auto intensityS = [&](double x){
		comp value = S_wave.getValue(pow(x,2))(ch);
		return (value*conj(value)).real();
	};

	auto intensityD = [&](double x){
		comp value = D_wave.getValue(pow(x,2))(ch);
		return (value*conj(value)).real();
	};


	//tested: readCh
	filereader testReader("GWA_dataformat.txt");
	
	string testSeed = testReader.getCommand(0);
	string testFitreg = testReader.getCommand(1);
	string testCh = testReader.getCommand(2);
	string testWav = testReader.getCommand(5);
	string testCheb = testReader.getCommand(6);
	string testPole = testReader.getCommand(9);
	string testKmat = testReader.getCommand(10);

	cout<<testSeed<<" -> "<<testReader.readSeed(testSeed)<<endl;

	vector<double> fitreg = testReader.readFitReg(testFitreg);
	cout<<testFitreg<<" -> "<<fitreg[0]<<", "<<fitreg[1]<<endl;

	chanDat testChanDat = testReader.readCh(testCh);
	cout<<testCh<<" -> "<<testChanDat.name<<" ";
	for(double x : testChanDat.ch_masses){
		cout<<x<<" ";
	}
	cout<<endl;

	ampDat testS = testReader.readWave(testWav);

	cout<<testWav<<" -> "<< testS.wname<<" "<< testS.kmatflag<< " "<<testS.rhoN<<" "<<testS.J<<" "<<testS.sL<<endl;


	/*
	for(ch = 0; ch <3; ch ++){
	makePlot("S_intensity_"+std::to_string(ch), intensityS);
	makePlot("D_intensity_"+std::to_string(ch), intensityD);
	}
	*/
	return 0;
}
