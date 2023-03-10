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
channel c1 = channel({0.13498, 0.13498});

//if you read AddWave("S","kmat","nominal", 0, 1.0)

amplitude a1 = amplitude(0,0,1.0,{c1});

//if you read ChebyCoeffs("S", "PiPi","p+s=1.0", {-2173.73, 3272.05, -1553.73, 361.79})

a1.setChebyCoeffs(0,3,1.0,{-2173.73, 3272.05, -1553.73, 361.79});

//if you read AddPole("S",0.00631603,{-6.05845,-7.24665,1.95572} )

a1.addPole(0.00631603,{-6.05845,-7.24665,1.95572});

//if you read AddKmatBackground("S", 0, {{16.1903, 18.3081 ,0}, {19.2388, 0}, {-15.811}})

a1.setKParams(0,{{16.1903, 18.3081 ,0}, {19.2388, 0}, {-15.811}});

}



int main()
{
	int J1 = 0;
	int J2 = 2;
	double alpha = 1;
	double sL = 1.0;
	double s0 = 1.0;
	double smin = comp(pow(0.997, 2), 0);
	comp smax = comp(pow(2.5, 2), 0);
	double mass1 = 0.13498;
	double mass2 = 0.49761;
	double mass3 = 0.762;
	// 1st wave
	channel chan1_1= channel({-6.0584496703077093,-0.23587420969943196,8.2951842306611070}, {-2173.7290154401403,3272.0517220458223,-1553.7255278262101, 361.78951904504214}, {mass1,mass1});
	channel chan2_1 = channel({ comp(-7.2466452617409232, 0), comp(-1.0918936513389781, 0), comp(12.908191623140738, 0)}, { comp(1927.4952345586573, 0.), comp(-2996.8378548561377, 0.), comp(1354.3272596278623, 0.), comp(-287.15761634573755, 0.)}, mass2);
	channel chan3_1 = channel({ comp(1.9557195508523364, 0), comp(-0.23250771304628870, 0), comp(-1.0490788402175895, 0)}, { comp(0, 0) },mass3);
	vector<channel> chans_1 = { chan1_1, chan2_1, chan3_1 };
	vector<double> rmasses_1 = {6.31602747547477250 * pow(10, -3), 3.1835288340590329, 7.5476155809082135};
	MatrixXcd c1 = MatrixXcd::Zero(3, 3);
	MatrixXcd d1 = MatrixXcd::Zero(3, 3);
	c1 << comp(16.190285088168821, 0), comp(18.308149119757218, 0), comp(0, 0), comp(18.308149119757218, 0), comp(19.238824892849152, 0), comp(0, 0), comp(0, 0), comp(0, 0), comp(-15.811018219003017, 0);
	d1 << comp(-6.2739277814634988, 0), comp(-9.1968627819824178, 0), comp(0, 0), comp(-9.1968627819824178, 0), comp(-13.370648723526756, 0), comp(0, 0), comp(0, 0), comp(0, 0), comp(6.0308208392725646, 0);
	vector<MatrixXcd> kparams_1 = { c1 , d1 };

	// 2nd wave
	channel chan1_2 = channel({ comp(0.15470117985842080, 0), comp(0.90811894619037048 , 0), comp(0.85226737492303073, 0)}, { comp(109.75915967336491,0.), comp(-2.8444511791990625,0.), comp(71.095211523264609,0.)}, mass1);
	channel chan2_2 = channel({ comp(1.0268866010574129, 0), comp(-0.24961757745677460, 0), comp(0.69753803634648648, 0)}, { comp(-26.963937424908881 , 0.), comp(-135.32369663565123, 0.), comp(-14.887225035532870, 0.)}, mass2);
	channel chan3_2 = channel({ comp(-0.34353840245057654, 0), comp(0, 0), comp(1.1340731281916305, 0)}, { comp(0, 0) },mass3);
	vector<channel> chans_2 = { chan1_2, chan2_2, chan3_2 };
	vector<double> rmasses_2 = {2.3299346099112341, 1.6673851773367765, 2.8936189128830971};
	MatrixXcd c2 = MatrixXcd::Zero(3, 3);
	MatrixXcd d2 = MatrixXcd::Zero(3, 3);
	c2 << comp(-1.1704615912349254, 0), comp(-1.0116586232616100, 0), comp(0, 0), comp(-1.0116586232616100, 0), comp(0.39464452011270623, 0), comp(0, 0), comp(0, 0), comp(0, 0), comp(-11.365225200233908, 0);
	d2 << comp(0.21138364931721298, 0), comp(0.16825103893825144, 0), comp(0, 0), comp(0.16825103893825144, 0), comp(4.86693780658242758 * pow(10, -3), 0), comp(0, 0), comp(0, 0), comp(0, 0), comp(3.9232752975622134, 0);
	vector<MatrixXcd> kparams_2 = { c2 , d2 };

	//Create the amplitude by passing it the channels, masses and k-matrix parameters
	amplitude wave_1 = amplitude(J1, alpha, sL, chans_1, kparams_1,rmasses_1,s0,smin,smax);
	amplitude wave_2 = amplitude(J2, alpha, sL, chans_2, kparams_2,rmasses_2,s0,smin,smax);

	int ch = 0;
	auto intensity = [&](double x){
		comp value = wave_1.getValue(pow(x,2))(ch);
		return (value*conj(value)).real();

	};
	for(ch = 0; ch<3;ch++){
	makePlot("intensity_"+std::to_string(ch), intensity);
	}

	ofstream outfile("wave_2.txt");
	outfile<<wave_2<<endl;
	outfile.close();

	return 0;
}
