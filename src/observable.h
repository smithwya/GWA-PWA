#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "amplitude.h"
#include "channel.h"
#include "filereader.h"
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
		vector<double> amp_expval, amp_expval_stat_err, amp_expval_sist_err;

		expchan(string wn, string cn, vector<double> x, vector<double> y, vector<double> y_stat_err, vector<double> y_sist_err){
			wavename = wn;
			channame = cn;
			sqrts = x;
			
			amp_expval = y;
			amp_expval_stat_err = y_stat_err;
			amp_expval_sist_err = y_sist_err;
		}


};

struct expInclCrossSec{

	vector<double> sqrts = {};
	vector<double> amp_expval = {}, amp_expval_stat_err = {}, amp_expval_sist_err = {};

	expInclCrossSec(vector<double> x, vector<double> y, vector<double> y_stat_err, vector<double> y_sist_err){
			sqrts = x;
			
			amp_expval = y;
			amp_expval_stat_err = y_stat_err;
			amp_expval_sist_err = y_sist_err;
	}

};


class observable {


private:
	
	vector<expchan> data;
	int numAmps;
	expInclCrossSec data_InclCrossSec = expInclCrossSec({},{},{},{});


public:
	vector<amplitude> amplitudes;
	int numChans;
	int nParams;
	observable(){
		amplitudes = {};
		data = {};
		data_InclCrossSec;
		numAmps = 0;
		numChans =0;
	};

	observable(vector<amplitude> a){
		amplitudes = a;
		data = {};
		data_InclCrossSec;
		numAmps = a.size();
		numChans = a[0].getNumOfChans();
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

	void setData_InclCrossSec(expInclCrossSec sigma){
		data_InclCrossSec = sigma;
	};

	expInclCrossSec getData_InclCrossSec(){
		return data_InclCrossSec;
	}

	void readData(string filename){

		return;
	};
	
	int getNumData(){
		int sum = 0;
		for(int i = 0; i < data.size(); i++){
			sum+=data[i].amp_expval.size();
		}
		return sum;
	}
	
	int getNumInclData(){
		int sum = 0;
		for(int i = 0; i < data.size(); i++){
			sum+=data_InclCrossSec.amp_expval.size();
		}
		return sum;
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
	};

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
		int totnumofchans = numChans; 
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
		//gr->GetYaxis()->SetRangeUser(0, 15500);
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
		int totnumofchans = numChans;
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
		int totnumofchans = numChans;
		int numchan = getchanindex(ampname,channame);
		
		int num_exp_pts = data[totnumofchans * numamp + numchan].amp_expval.size();

		double x1[num_exp_pts], y1[num_exp_pts], ex1[num_exp_pts], ey1[num_exp_pts];

		const int num_th_pts = 100;
		double x2[num_th_pts], y2[num_th_pts], ex2[num_th_pts], ey2[num_th_pts];

		for(int i = 0; i < num_exp_pts; i++){

			x1[i] = 0;
			y1[i] = 0;
			ex1[i] = 0;
			ey1[i] = 0;

		}

		double lb = data[totnumofchans * numamp + numchan].sqrts[0];

		double ub = data[totnumofchans * numamp + numchan].sqrts[num_exp_pts - 1];

		for (int i = 0; i < num_th_pts; i++){

			x2[i] = lb + (ub - lb) * i / ((double)num_th_pts - 1.);

		}

		double val = 0;

		for(int i = 0; i < num_exp_pts; i++){

			val = data[totnumofchans * numamp + numchan].sqrts[i];
			if(isnan(val)) x1[i] = 0;
			else x1[i] = val;


			val = data[totnumofchans * numamp + numchan].amp_expval[i];
			if(isnan(val)) y1[i] = 0;
			else y1[i] = val;


			val = pow(data[totnumofchans * numamp + numchan].amp_expval_stat_err[i], 2);
			val += pow(data[totnumofchans * numamp + numchan].amp_expval_sist_err[i], 2);
			val = pow(val, 0.5);
			if(isnan(val)) ey1[i] = 0;
			else ey1[i] = val;

		}

		for(int i = 0; i < num_th_pts; i++){

			val = func(x2[i]);
			if(isnan(val)) y2[i] = 0;
			else y2[i] = val;

		}

		auto gr1 = new TGraphErrors(num_exp_pts,x1,y1,ex1,ey1);
		auto gr2 = new TGraph(num_th_pts,x2,y2);

   		//gr1->SetTitle("TGraphErrors Example");
   		gr1->SetMarkerColor(4);
   		gr1->SetMarkerStyle(21);
		gr1->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		//gr1->GetYaxis()->SetRangeUser(0, 0.2);
		gr1->SetLineWidth(1);

		//gr2->SetTitle("TGraphErrors Example");
   		gr2->SetMarkerSize(0);
		gr2->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		//gr2->GetYaxis()->SetRangeUser(0, 0.2);
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
		int totnumofchans = numChans;
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
		int totnumofchans = numChans;
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
		//gr->GetYaxis()->SetRangeUser(0, 15500);
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

	void plotInclCrossSec(string pdfname, double lower_bound, double upper_bound){
		
		int num_exp_pts = data_InclCrossSec.sqrts.size();

		double x[num_exp_pts], y[num_exp_pts], ex[num_exp_pts], ey[num_exp_pts];

		for(int i = 0; i < num_exp_pts; i++){

			x[i] = 0;
			y[i] = 0;
			ex[i] = 0;
			ey[i] = 0;

		}

		double val = 0;

		for(int i = 0; i < num_exp_pts; i++){

			val = data_InclCrossSec.sqrts[i];
			if(isnan(val)) x[i] = 0;
			else x[i] = val;

			val = data_InclCrossSec.amp_expval[i];
			if(isnan(val)) y[i] = 0;
			else y[i] = val;

			ex[i] = 0;

			val = data_InclCrossSec.amp_expval_stat_err[i];
			if(isnan(val)) ey[i] = 0;
			else ey[i] = val;

		}

		auto gr = new TGraphErrors(num_exp_pts,x,y,ex,ey);

   		//gr->SetTitle("TGraphErrors Example");
   		gr->SetMarkerColor(4);
   		gr->SetMarkerStyle(21);
		gr->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		//gr->GetYaxis()->SetRangeUser(0, 15500);
		gr->SetLineWidth(1);

		TFile file("pdf_folder.root", "recreate");
		TCanvas canv;
		gr->Write();
		gr->Draw("AP");
		canv.SaveAs(("Plots/"+pdfname+".pdf").c_str());
		file.Close();
		return;

	}

	void plotInclCrossSecVsSumOfExcl(string pdfname, double lower_bound, double upper_bound){

		int totnumofmeasuredchans = data.size();
		
		int num_excl_pts = data[0].amp_expval.size();
		int num_Bsstar_pts = data[getchanindex("P", "B_sstarB_sstar")].amp_expval.size();
		//cout << getchanindex("P", "B_sstarB_sstar") << endl;
		int num_incl_pts = data_InclCrossSec.amp_expval.size();

		//cout << num_excl_pts << " " << num_Bsstar_pts << " " << num_incl_pts << endl; //I have
		//to exclude B_sstarB_sstar channel from gr1 because it has a different number of points
		//with respect of the other exclusive sigma!!! 

		double x1[num_excl_pts], y1[num_excl_pts], ex1[num_excl_pts], ey1[num_excl_pts];
		double x2[num_Bsstar_pts], y2[num_Bsstar_pts], ex2[num_Bsstar_pts], ey2[num_Bsstar_pts];
		double x3[num_incl_pts], y3[num_incl_pts], ex3[num_incl_pts], ey3[num_incl_pts];

		double val = 0;

		for(int i = 0; i < num_excl_pts; i++){

			val = data[0].sqrts[i];
			if(isnan(val)) x1[i] = 0;
			else x1[i] = val;
			ex1[i] = 0;

			val = 0;
			for(int j = 0; j < 3; j++) val += data[j].amp_expval[i];
			if(isnan(val)) y1[i] = 0;
			else y1[i] = val;


			val = 0;
			for(int j = 0; j < 3; j++) val += pow(data[j].amp_expval_stat_err[i], 2);
			val = sqrt(val);
			if(isnan(val)) ey1[i] = 0;
			else ey1[i] = val;

		}

		for(int i = 0; i < num_Bsstar_pts; i++){

			int index = getchanindex("P", "B_sstarB_sstar");

			val = data[index].sqrts[i];
			if(isnan(val)) x2[i] = 0;
			else x2[i] = val;
			//cout << val << endl;
			ex2[i] = 0;
			cout << val << " ";

			val = data[index].amp_expval[i];
			if(isnan(val)) y2[i] = 0;
			else y2[i] = val;
			cout << val << " ";


			val = pow(data[index].amp_expval_stat_err[i], 2);
			val += pow(data[index].amp_expval_sist_err[i], 2);
			val = pow(val, 0.5);
			if(isnan(val)) ey2[i] = 0;
			else ey2[i] = val;
			cout << val << endl;

		}

		for(int i = 0; i < num_incl_pts; i++){

			val = data_InclCrossSec.sqrts[i];
			if(isnan(val)) x3[i] = 0;
			else x3[i] = val;
			ex3[i] = 0;

			val += data_InclCrossSec.amp_expval[i];
			if(isnan(val)) y3[i] = 0;
			else y3[i] = val;

			val = pow(data_InclCrossSec.amp_expval_stat_err[i], 2);
			val += pow(data_InclCrossSec.amp_expval_sist_err[i], 2);
			val = pow(val, 0.5);
			if(isnan(val)) ey3[i] = 0;
			else ey3[i] = val;

		}

		auto gr1 = new TGraphErrors(num_excl_pts,x1,y1,ex1,ey1);
		auto gr2 = new TGraphErrors(num_Bsstar_pts,x2,y2,ex2,ey2);
		auto gr3 = new TGraphErrors(num_incl_pts,x3,y3,ex3,ey3);

   		//gr1->SetTitle("TGraphErrors Example");
   		gr1->SetMarkerColor(kBlue);
   		gr1->SetMarkerStyle(21);
		gr1->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		gr1->GetYaxis()->SetRangeUser(0,450.);
		gr1->SetLineWidth(1);
		gr1->SetLineColor(kBlue);

		//gr2->SetTitle("TGraphErrors Example");
   		gr2->SetMarkerColor(kRed);
		gr2->SetMarkerStyle(21);
		gr2->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		//gr2->GetYaxis()->SetRangeUser(0, 0.2);
		gr2->SetLineWidth(1);
		gr2->SetLineColor(kRed);

		//gr3->SetTitle("TGraphErrors Example");
   		gr3->SetMarkerColor(kMagenta);
		gr3->SetMarkerStyle(21);
		gr3->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		//gr3->GetYaxis()->SetRangeUser(0, 0.2);
		gr3->SetLineWidth(1);
		gr3->SetLineColor(kMagenta);

		TFile file("pdf_folder.root", "recreate");
		TCanvas canv;
		gr1->Write();
		gr2->Write();
		gr3->Write();
		gr1->Draw("AP");
		gr2->Draw("Psame");
		gr3->Draw("Psame");
		canv.SaveAs(("Plots/"+pdfname+".pdf").c_str());
		file.Close();
		return;

	}


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


	double chisq(){

		double result = 0;

		//for every amplitude in the list
		for(string ampname : getAmpNames()){

			//locate the amplitude
			int amp_index = getampindex(ampname);
			amplitude amp = amplitudes[amp_index];
			//read out the start and end of the fit interval
			double lower_bound = sqrt(amp.getFitInterval()[0]);
			double upper_bound = sqrt(amp.getFitInterval()[1]);	

			//for every channel
			for(string channame : amp.getChanNames()){

				int chan_index = getchanindex(ampname, channame);
				double sum = 0;
				vector<double> sqrts_vals = data[amp_index * numChans + chan_index].sqrts;

				if(sqrts_vals.size() > 0){//i.e. if the j-th channel is NOT a dummy channel

					for(int i = 0; i < sqrts_vals.size(); i++){

						double x = sqrts_vals[i];

						//for sqrts in the fitting range, calculate relevant quantities
						if(x >= lower_bound && x <= upper_bound){

							comp val = amp.getValue(pow(x, 2))(chan_index);
							double y = (val*conj(val)).real();
							double stat_err = data[numChans * amp_index + chan_index].amp_expval_stat_err[i];
							double sist_err = data[numChans * amp_index + chan_index].amp_expval_sist_err[i];
							double std = sqrt(pow(stat_err, 2) + pow(sist_err, 2));
							sum += pow(((y - data[numChans * amp_index + chan_index].amp_expval[i])/std), 2);

						}

					}

					result += sum;

				}

			}

		}

		return result;

	}


	double chisq_with_InclCrossSec(){

		double sum = chisq();
		double lower_bound = sqrt(amplitudes[0].getFitInterval()[0]);
		double upper_bound = sqrt(amplitudes[0].getFitInterval()[1]);
		comp temp = 0;

		for(int i = 0; i < data_InclCrossSec.sqrts.size(); i++){

			double x = data_InclCrossSec.sqrts[i];
			
			double aux = 0;
			
			if(x >= lower_bound && x <= upper_bound){

				for(string ampname : getAmpNames()){

					int amp_index = getampindex(ampname);
					amplitude amp = amplitudes[amp_index];

					for(string channame : amp.getChanNames()){

						int chan_index = getchanindex(ampname, channame);
						temp = amp.getValue(pow(x,2))(chan_index);
						aux += (temp*conj(temp)).real();

					}

				}

				double y = data_InclCrossSec.amp_expval[i];
				double stat_err = data_InclCrossSec.amp_expval_stat_err[i];
				double sist_err = data_InclCrossSec.amp_expval_sist_err[i];
				//uncorrelated inclusive data
				double std = stat_err;
				//correlated inclusive data
				//double std = sqrt(pow(stat_err, 2) + pow(sist_err, 2));
				sum += pow(((aux - y)/std), 2);

			}

		}

		return sum;

	}

	friend ostream& operator<<(std::ostream& os, observable const& m){
		for(amplitude a: m.amplitudes) cout<<a<<endl<<endl;
		return os;
	};


};
