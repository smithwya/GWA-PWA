#include <iostream>
#include<complex>
#include<vector>
#include "amplitude.h"
#include "channel.h"
#include "TCanvas.h"
#include <Eigen/Dense>
#include "TH1D.h"
#include "TFile.h"
#include <sstream>
#include <fstream>
#include "TString.h"
#include <chrono>
#include <map>

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

//if you read AddChannel("pipi", 0.13498, 0.13498)
channel c1 = channel("PiPi",{0.13498, 0.13498});

//if you read AddWave("S","kmat","nominal", 0, 1.0)

amplitude a1 = amplitude(0,0,1.0,{c1});

//if you read ChebyCoeffs("S", "PiPi","p+s=1.0", {-2173.73, 3272.05, -1553.73, 361.79})

a1.setChebyCoeffs(0,3,1.0,{-2173.73, 3272.05, -1553.73, 361.79});

//if you read AddPole("S",0.00631603,{-6.05845,-7.24665,1.95572} )

a1.addPole(0.00631603,{"PiPi","KK","RhoRho"},{-6.05845,-7.24665,1.95572});

//if you read AddKmatBackground("S", 0, {{16.1903, 18.3081 ,0}, {19.2388, 0}, {-15.811}})

a1.setKParams(0,{{16.1903, 18.3081 ,0}, {19.2388, 0}, {-15.811}});

return {a1};
}



int main()
{

	double smin = pow(0.9975,2);
	double smax = pow(2.5,2);

	channel c1 = channel("PiPi",{0.13498, 0.13498});
	channel c2 = channel("KK", {0.49761, 0.49761});
	channel c3 = channel("RhoRho", {0.762, 0.762});

	amplitude S_wave = amplitude(0,smin,smax,{c1,c2,c3});

	S_wave.setChebyCoeffs("PiPi",3, 1.0, {-2173.73, 3272.05, -1553.73, 361.79});
	S_wave.setChebyCoeffs("KK",3, 1.0, {1927.5, -2996.84, 1354.33 , -287.158});
	S_wave.setChebyCoeffs("RhoRho", 3,1.0, {0});

	S_wave.addPole(0.00631603, {"PiPi","KK","RhoRho"},{-6.05845,-7.24665,1.95572});
	S_wave.addPole(3.18353,{"PiPi","KK","RhoRho"},{-0.235874,-1.09189,-0.232508});
	S_wave.addPole(7.54762,{"PiPi","KK","RhoRho"},{8.29518,12.9082,-1.04908});

	S_wave.setKParams(0,{{16.1903, 18.3081 ,0}, {19.2388, 0}, {-15.811}});
	S_wave.setKParams(1, {{-6.27393,-9.19686,0}, {-13.3706,0}, {6.03082}});
	
	amplitude D_wave = amplitude(2,smin,smax,{c1,c2,c3});


	D_wave.setChebyCoeffs("PiPi", 3, 1.0, {109.759,-2.84445, 71.0952});
	D_wave.setChebyCoeffs("KK", 3, 1.0 , {-26.9639, -135.324, -14.8872});
	D_wave.setChebyCoeffs("RhoRho",3,1.0, {0});

	D_wave.addPole( 2.32993, {"PiPi","KK","RhoRho"}, {0.154701, 0.908119, 0.852267});
	D_wave.addPole( 1.66739, {"PiPi","KK","RhoRho"}, {1.02689, -0.249618, 0.697538});
	D_wave.addPole( 2.89362, {"PiPi","KK","RhoRho"}, {-0.343538, 0, 1.13407});

	D_wave.setKParams(0, {{-1.17046,-1.01166, 0 }, {0.394645, 0}, {-11.3652}});
	D_wave.setKParams(1, {{0.211384, 0.168251, 0}, {4.86693780658242758 * pow(10, -3),0}, {3.92328}});

	cout<<"S-wave: "<<S_wave<<endl;


	cout<<"D-wave:"<<D_wave<<endl;
		
	
	int ch = 0;
	auto intensity = [&](double x){
		comp value = S_wave.getValue(pow(x,2))(ch);
		return (value*conj(value)).real();

	};

	auto start = chrono::high_resolution_clock::now();
	makePlot("S_intensity_"+std::to_string(ch), intensity);
	auto stop = chrono::high_resolution_clock::now();

	auto duration1 = chrono::duration_cast<chrono::microseconds>(stop-start);
	cout<<"First plot generation time: "<< duration1.count()<<" microseconds"<<endl;

	start = chrono::high_resolution_clock::now();
	makePlot("S_intensity_"+std::to_string(ch), intensity);
	stop = chrono::high_resolution_clock::now();

	auto duration2 = chrono::duration_cast<chrono::microseconds>(stop-start);
	cout<<"First plot generation time: "<< duration2.count()<<" microseconds"<<endl;

	
	return 0;
}
