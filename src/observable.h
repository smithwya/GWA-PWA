#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "amplitude.h"
#include "channel.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TString.h"
#include <Math/ParamFunctor.h>
using namespace std;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;


struct expchan{//this has to be modify according to the format of the exp data file

		string wavename;
		string channame;
		vector<double> sqrts;
		vector<double> amp_expval, amp_expval_stat_err;

		expchan(string wn, string cn, vector<double> x, vector<double> y, vector<double> y_stat_err){
			wavename = wn;
			channame = cn;
			sqrts = x;
			
			amp_expval = y;
			amp_expval_stat_err = y_stat_err;
		}


};


class observable {


private:
	
	vector<expchan> data;
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

	vector<string> getAmpNames(){
		vector<string> temp;
		for(amplitude amp: getAmps()){
			temp.push_back(amp.getName());
		}
		return temp;
	}

	int getNumAmps(){
		return numAmps;
	};

	amplitude getAmp(int i){
		return amplitudes[i];
	};

	void setData(vector<expchan> dat){
		data=dat;
	};

	vector<expchan> getData(){
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

	/*
	void makePlotGraph(string pdfname, function<double(double)> func, double lower_bound, double upper_bound){	

		auto myFunc = [&](double x){return func(x);};

		TF1 *f = new TF1("f", "myFunc(x)", lower_bound, upper_bound, 0);

		auto gr = new TGraph(f);

		TFile file("pdf_folder.root", "recreate");
		TCanvas canv;
		//gr->Write();
		//gr->Draw();
		f->Draw();
		canv.SaveAs(("Plots/"+pdfname+".pdf").c_str());
		file.Close();
		return;

	}; //not working!
	*/

	int getampindex(string ampname){

		vector<string> amp_names = getAmpNames();
		auto it = find(amp_names.begin(),amp_names.end(),ampname);

		if(it==amp_names.end()){
			cout<<"amplitude doesn't exist"<<endl;
			return -1;
		}
		return it-amp_names.begin();

	}

	int getchanindex (string ampname, string channame){

		int amp_index = getampindex(ampname);

		vector<string> chan_names = amplitudes[amp_index].getChanNames();
		auto it = find(chan_names.begin(),chan_names.end(),channame);

		if(it==chan_names.end()){
			cout<<"channel doesn't exist"<<endl;
			return -1;
		}
		return it-chan_names.begin();


	}
	void makePlotGraph(string ampname, string channame, string pdfname, function<double(double)> func, double lower_bound, double upper_bound){	
		
		int numamp = getampindex(ampname);
		int totnumofchans = data.size(); //in this way is always consistent with the case in which we consider or not the dummy channels
		int numchan = getchanindex(ampname,channame);

		int numpts = data[totnumofchans * numamp + numchan].amp_expval.size();

		double delta = (upper_bound - lower_bound)/ numpts;

		double x[numpts], y[numpts], ex[numpts], ey[numpts];

		for(int i = 0; i < numpts; i++){

			x[i] = 0;
			y[i] = 0;
			ex[i] = 0;
			ey[i] = 0;

		}

		double val = 0;

		for(int i = 0; i < numpts; i++){

			val = data[totnumofchans * numamp + numchan].sqrts[i];
			if(isnan(val)) x[i] = 0;
			else x[i] = val;

			val = func(x[i]);
			if(isnan(val)) y[i] = 0;
			else y[i] = val;

		}

		auto gr = new TGraphErrors(numpts,x,y,ex,ey);

   		//gr->SetTitle("TGraphErrors Example");
   		gr->SetMarkerColor(4);
   		gr->SetMarkerStyle(21);
		gr->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		gr->GetYaxis()->SetRangeUser(0, 15500);
		gr->SetLineWidth(1);

		TFile file("pdf_folder.root", "recreate");
		TCanvas canv;
		gr->Write();
		gr->Draw("ALP");
		canv.SaveAs(("Plots/"+pdfname+".pdf").c_str());
		file.Close();
		return;

	}; 

	void makePlotWithExp(string ampname, string channame, string pdfname, function<double(double)> func, double lower_bound, double upper_bound, int num_bins){

		int numamp = getampindex(ampname);
		int totnumofchans = data.size();
		int numchan = getchanindex(ampname,channame);
		
		double delta = (upper_bound - lower_bound)/num_bins;
		int num_data_points = data[totnumofchans * numamp + numchan].amp_expval.size();

		TH1D *plotter = new TH1D(pdfname.c_str(),pdfname.c_str(), num_bins, lower_bound, upper_bound);
		//plotter->SetMinimum(-60.0);
		//plotter->SetMaximum(90.0);
		double sqrtS = 1.0;
		TH1D *plotter_exp = new TH1D("","", num_data_points, lower_bound, upper_bound);

		for(int i = 0; i < num_bins; i++){
			sqrtS = lower_bound + delta/2 + i * delta; //centroid
			double val = func(sqrtS);
			if(isnan(val)) plotter->SetBinContent(i + 1,  0);
			else plotter->SetBinContent(i + 1,  val);
		}

		for(int i = 0; i < num_data_points; i++){
			double val_exp = data.at(totnumofchans * numamp + numchan).amp_expval.at(i);
			if(isnan(val_exp))plotter_exp->SetBinContent(i + 1,  0);
			else plotter_exp->SetBinContent(i + 1,  val_exp);
		}

		


		TFile file("pdf_folder.root", "recreate");
		TCanvas canv;
		plotter->Write();
		plotter_exp->Write();
		plotter->Draw();
		plotter_exp->Draw("same");
		canv.SaveAs(("Plots/"+pdfname+".pdf").c_str());
		file.Close();
		return;
	};

	void makePlotGraphWithExp(string ampname, string channame, string pdfname, function<double(double)> func, double lower_bound, double upper_bound){

		int numamp = getampindex(ampname);
		int totnumofchans = data.size();
		int numchan = getchanindex(ampname,channame);
		
		int num_exp_pts = data[totnumofchans * numamp + numchan].amp_expval.size();

		double x1[num_exp_pts], y1[num_exp_pts], ex1[num_exp_pts], ey1[num_exp_pts];

		double x2[num_exp_pts], y2[num_exp_pts], ex2[num_exp_pts], ey2[num_exp_pts];

		for(int i = 0; i < num_exp_pts; i++){

			x1[i] = 0;
			y1[i] = 0;
			ex1[i] = 0;
			ey1[i] = 0;

			x2[i] = 0;
			y2[i] = 0;
			ex2[i] = 0;
			ey2[i] = 0;

		}

		double val = 0;

		for(int i = 0; i < num_exp_pts; i++){

			val = data[totnumofchans * numamp + numchan].sqrts[i];
			if(isnan(val)) x1[i] = 0;
			else x1[i] = val;

			val = data[totnumofchans * numamp + numchan].sqrts[i];
			if(isnan(val)) x2[i] = 0;
			else x2[i] = val;

			val = data[totnumofchans * numamp + numchan].amp_expval[i];
			if(isnan(val)) y1[i] = 0;
			else y1[i] = val;

			val = func(x2[i]);
			if(isnan(val)) y2[i] = 0;
			else y2[i] = val;

			val = data[totnumofchans * numamp + numchan].amp_expval_stat_err[i];
			if(isnan(val)) ey1[i] = 0;
			else ey1[i] = val;

		}

		auto gr1 = new TGraphErrors(num_exp_pts,x1,y1,ex1,ey1);
		auto gr2 = new TGraphErrors(num_exp_pts,x2,y2,ex2,ey2);

   		//gr1->SetTitle("TGraphErrors Example");
   		gr1->SetMarkerColor(4);
   		gr1->SetMarkerStyle(21);
		gr1->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		gr1->GetYaxis()->SetRangeUser(0, 15500);
		gr1->SetLineWidth(1);

		//gr2->SetTitle("TGraphErrors Example");
   		gr2->SetMarkerSize(0);
		gr2->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		gr2->GetYaxis()->SetRangeUser(0, 15500);
		gr2->SetLineWidth(1);
		gr2->SetLineColor(kRed);

		TFile file("pdf_folder.root", "recreate");
		TCanvas canv;
		gr1->Write();
		gr2->Write();
		gr1->Draw("AP");
		gr2->Draw("same");
		canv.SaveAs(("Plots/"+pdfname+".pdf").c_str());
		file.Close();
		return;
		
	};

	void makePlotExpOnly(string ampname, string channame, string pdfname, double lower_bound, double upper_bound){

		int numamp = getampindex(ampname);
		int totnumofchans = data.size();
		int numchan = getchanindex(ampname,channame);

		int num_data_points = data[totnumofchans * numamp + numchan].amp_expval.size();

		double sqrtS = 1.0;
		TH1D *plotter_exp = new TH1D("","", num_data_points, lower_bound, upper_bound);

		for(int i = 0; i < num_data_points; i++){
			double val_exp = data.at(totnumofchans * numamp + numchan).amp_expval.at(i);
			if(isnan(val_exp))plotter_exp->SetBinContent(i + 1,  0);
			else plotter_exp->SetBinContent(i + 1,  val_exp);
		}

		TFile file("pdf_folder.root", "recreate");
		TCanvas canv;
	
		plotter_exp->Write();
		plotter_exp->Draw();
		canv.SaveAs(("Plots/"+pdfname+".pdf").c_str());
		file.Close();
		return;
	};

	void makePlotGraph_ExpOnly(string ampname, string channame, string pdfname, double lower_bound, double upper_bound){	

		int numamp = getampindex(ampname);
		int totnumofchans = data.size();
		int numchan = getchanindex(ampname,channame);

		int num_exp_pts = data[totnumofchans * numamp + numchan].amp_expval.size();

		double delta = (upper_bound - lower_bound)/ num_exp_pts;

		double x[num_exp_pts], y[num_exp_pts], ex[num_exp_pts], ey[num_exp_pts];

		for(int i = 0; i < num_exp_pts; i++){

			x[i] = 0;
			y[i] = 0;
			ex[i] = 0;
			ey[i] = 0;

		}

		double val = 0;

		for(int i = 0; i < num_exp_pts; i++){

			val = data[totnumofchans * numamp + numchan].sqrts[i];
			if(isnan(val)) x[i] = 0;
			else x[i] = val;

			val = data[totnumofchans * numamp + numchan].amp_expval[i];
			if(isnan(val)) y[i] = 0;
			else y[i] = val;

			ex[i] = 0;

			val = data[totnumofchans * numamp + numchan].amp_expval_stat_err[i];
			if(isnan(val)) ey[i] = 0;
			else ey[i] = val;

		}

		auto gr = new TGraphErrors(num_exp_pts,x,y,ex,ey);

   		//gr->SetTitle("TGraphErrors Example");
   		gr->SetMarkerColor(4);
   		gr->SetMarkerStyle(21);
		gr->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		gr->GetYaxis()->SetRangeUser(0, 15500);
		gr->SetLineWidth(1);

		TFile file("pdf_folder.root", "recreate");
		TCanvas canv;
		gr->Write();
		gr->Draw("ALP");
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

	void plotCompGraph(string ampname, string channame, string pdfname,function<comp(double)> func, double lower_bound, double upper_bound){

		auto realFunc = [&](double x){
			return func(x).real();
		};
		auto imagFunc = [&](double x){
			return func(x).imag();
		};

		makePlotGraph(ampname, channame, "Re "+pdfname, realFunc, lower_bound, upper_bound);
		makePlotGraph(ampname, channame, "Im "+pdfname, imagFunc, lower_bound, upper_bound);

	};

	void plotCompWithExp(string ampname, string channame, string pdfname,function<comp(double)> func, double lower_bound, double upper_bound, int num_bins){

		auto realFunc = [&](double x){
			return func(x).real();
		};
		auto imagFunc = [&](double x){
			return func(x).imag();
		};

		makePlotWithExp(ampname, channame, "Exp_plus_Re "+pdfname, realFunc, lower_bound, upper_bound, num_bins);
		makePlotWithExp(ampname, channame, "Exp_plus_Im "+pdfname, imagFunc, lower_bound, upper_bound, num_bins);

	};

	void plotCompGraphWithExp(string ampname, string channame, string pdfname,function<comp(double)> func, double lower_bound, double upper_bound, int num_bins){

		auto realFunc = [&](double x){
			return func(x).real();
		};
		auto imagFunc = [&](double x){
			return func(x).imag();
		};

		makePlotGraphWithExp(ampname, channame, "Exp_plus_Re "+pdfname, realFunc, lower_bound, upper_bound);
		makePlotGraphWithExp(ampname, channame, "Exp_plus_Im "+pdfname, imagFunc, lower_bound, upper_bound);

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