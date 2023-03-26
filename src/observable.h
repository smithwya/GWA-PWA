#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "amplitude.h"
#include "channel.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
using namespace std;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;


class observable {


private:
	vector<channel> channels;
	vector<amplitude> amplitudes;
	vector<pair<double,double>> data;

public:
	observable(){
		channels = {};
		amplitudes = {};
		data = {};
	};

	observable(vector<channel> c, vector<amplitude> a){
		channels = c;
		amplitudes = a;
		data = {};
	};

	vector<channel> getChannels(){
		return channels;
	};

	vector<amplitude> getAmps(){
		return amplitudes;
	};

	amplitude getAmp(int i){
		return amplitudes[i];
	};

	void setData(vector<pair<double,double>> dat){
		data = dat;
	}

	void readData(string filename){

		return;
	}

	void makePlot(string pdfname, function<double(double)> func, double lower_bound, double upper_bound, int num_bins){

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

	void plotComp(string pdfname,function<comp(double)> func, double lower_bound, double upper_bound, int num_bins){

		auto realFunc = [&](double x){
			return func(x).real();
		};
		auto imagFunc = [&](double x){
			return func(x).imag();
		};

		makePlot("Re "+pdfname, realFunc, lower_bound, upper_bound, num_bins);
		makePlot("Im "+pdfname, imagFunc, lower_bound, upper_bound, num_bins);

	}

};