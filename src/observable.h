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

struct expchan{

		string wavename;
		string channame;
		vector<double> sqrts;
		vector<double> amp_expval;

		expchan(string wn, string cn, vector<double> x, vector<double> y){
			wavename = wn;
			channame = cn;
			sqrts = x;
			amp_expval = y;
		}


};

class observable {


private:
	
	vector<vector<expchan>> data;
	int numAmps;


public:
	vector<amplitude> amplitudes;
	int numChans;
	observable(){
		amplitudes = {};
		data = {};
		numAmps = 0;
		numChans =0;
	};

	observable(vector<amplitude> a){
		amplitudes = a;
		data = {};
		numAmps = a.size();
		numChans = a[0].getChanNames().size();
	};

	vector<amplitude> getAmps(){
		return amplitudes;
	};

	int getNumAmps(){
		return numAmps;
	};

	amplitude getAmp(int i){
		return amplitudes[i];
	};

	void addData(vector<expchan> dat){
		data.push_back(dat);
	};

	vector<vector<expchan>> getData(){
		return data;
	}

	void readData(string filename){

		return;
	};

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
	};

	void plotComp(string pdfname,function<comp(double)> func, double lower_bound, double upper_bound, int num_bins){

		auto realFunc = [&](double x){
			return func(x).real();
		};
		auto imagFunc = [&](double x){
			return func(x).imag();
		};

		makePlot("Re "+pdfname, realFunc, lower_bound, upper_bound, num_bins);
		makePlot("Im "+pdfname, imagFunc, lower_bound, upper_bound, num_bins);

	};


	vector<double> getFitParams(){
		vector<double> params = {};
		for(amplitude a : amplitudes){
			vector<double> temp = a.getFittedParamList();
			params.insert(params.end(),temp.begin(),temp.end());
		}
		return params;
	};

	void setFitParams(vector<double> newparams){
		int index = 0;
		for(int i = 0; i < amplitudes.size(); i++){
			int nparams = amplitudes[i].getFittedParamList().size();

			vector<double> ampparams(newparams.begin()+index,newparams.begin()+index+nparams);
			amplitudes[i].setFittedParamList(ampparams);
			index+=nparams;
		}

	}

	vector<double> getStepSizes(){
		vector<double> steps = {};
		for(amplitude a : amplitudes){
			vector<double> temp = a.getFittedSteps();
			steps.insert(steps.end(),temp.begin(),temp.end());
		}
		return steps;
	};


	friend ostream& operator<<(std::ostream& os, observable const& m){
		for(amplitude a: m.amplitudes) cout<<a<<endl<<endl;
		return os;
	};


};