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
#include <TGraph2D.h>
#include <TF2.h>
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
	double incl_weight = 1.;
	double excl_weight = 1.;

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

	void makePlotGraphDummy(string pdfname, function<double(double)> func, double lower_bound, double upper_bound){

		const int num_th_pts = 1000;
		double x[num_th_pts], y[num_th_pts], ex[num_th_pts], ey[num_th_pts];

		double lb = lower_bound;
		double ub = upper_bound;

		for (int i = 0; i < num_th_pts; i++){

			x[i] = lb + (ub - lb) * i / ((double)num_th_pts - 1.);

		}

		double val = 0;	

		for(int i = 0; i < num_th_pts; i++){

			val = func(x[i]);
			if(isnan(val)) y[i] = 0;
			else y[i] = val;

		}

		auto gr1 = new TGraph(num_th_pts,x,y);

		int mymarkerstyle=20;
  		float mymarkersize=0;
  		int mymarkercolor = 4;
  		float mytextsize=0.04;
  		int mytextfont=132;

		TCanvas *Scatola = new TCanvas("Scatola","Scatola",600,500); //costruttore 600pt x 550 pt
  		//gStyle->SetOptStat(0); //non voglio che mi metti il riquadro con la statistica
  		Scatola->SetFillColor(0);//il fondo del grafico con 0 è bianco...in teoria lo potete cambiare
  		Scatola->SetBorderMode(0);//mette dei bordi attorno alla figura...0 nessun bordo
  		Scatola->SetBorderSize(2); //spessore del bordo
  		Scatola->SetLeftMargin(0.18); //spazio a sinistra della figura ...20% della larghezza
  		Scatola->SetRightMargin(0.11);// 5% della larghezza a destra
  		Scatola->SetTopMargin(0.07); //3% della altezza lato superiore
  		Scatola->SetBottomMargin(0.14); //12% dell'altezza lato inferiore
  		Scatola->SetTickx(1);
  		Scatola->SetTicky(1);

		//gr1->SetTitle("Dummy channel");
		gr1->GetYaxis()->SetTitleSize(mytextsize); //controllo sulla dimension del titolo dell'asse
  		gr1->GetXaxis()->SetTitleSize(mytextsize);
  		gr1->GetXaxis()->SetLabelSize(mytextsize);//cotrollo sulla dimensione dei numeretti dell'asse
  		gr1->GetYaxis()->SetLabelSize(mytextsize);
  		gr1->GetXaxis()->SetTitleFont(mytextfont);//controllo sul carattere usato per il titolo dell'asse
  		gr1->GetYaxis()->SetTitleFont(mytextfont);
  		gr1->GetXaxis()->SetLabelFont(mytextfont);//controllo sul carattere usato per i numeretti dell'asse
  		gr1->GetYaxis()->SetLabelFont(mytextfont);
  		gr1->GetXaxis()->SetNdivisions(908); //suddivisione dei numeri sull'asse x---es 0 a 10 a passo di 1, e ogni passo diviso in 5
  		gr1->GetXaxis()->CenterTitle(1);//che il titolo dell'asse lo voglio quindi 1, se non lo volessi metterei 0
  		gr1->GetYaxis()->CenterTitle(1);
  		gr1->GetXaxis()->SetTitleOffset(1.15);//definisce la distanza del titolo dell'asse dall'asse stesso
  		gr1->GetYaxis()->SetTitleOffset(1.15);
  		gr1->SetTitle("");
  		gr1->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  		gr1->GetYaxis()->SetTitle("#sigma (pb)");
  		gr1->SetMarkerColor(mymarkercolor);
  		gr1->SetMarkerSize(mymarkersize);
  		gr1->SetMarkerStyle(mymarkerstyle);
  		gr1->SetLineColor(kBlue);
  		gr1->SetLineWidth(1);
		gr1->GetXaxis()->SetRangeUser(lb,ub);
		//gr1->GetYaxis()->SetRangeUser(0,450.);
		gr1->SetLineWidth(1);
		//gr1->SetLineColor(kMagenta);

		TFile file("pdf_folder.root", "recreate");
		
		gr1->Write();
		gr1->Draw("AL");
		Scatola->SaveAs(("Plots/"+pdfname+".pdf").c_str());
		file.Close();
		return;

	}

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

		const int num_th_pts = 1000;
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

   		int mymarkerstyle=20;
  		float mymarkersize=1.;
  		int mymarkercolor = 4;
  		float mytextsize=0.04;
  		int mytextfont=132;

		TCanvas *Scatola = new TCanvas("Scatola","Scatola",600,500); //costruttore 600pt x 550 pt
  		//gStyle->SetOptStat(0); //non voglio che mi metti il riquadro con la statistica
  		Scatola->SetFillColor(0);//il fondo del grafico con 0 è bianco...in teoria lo potete cambiare
  		Scatola->SetBorderMode(0);//mette dei bordi attorno alla figura...0 nessun bordo
  		Scatola->SetBorderSize(2); //spessore del bordo
  		Scatola->SetLeftMargin(0.18); //spazio a sinistra della figura ...20% della larghezza
  		Scatola->SetRightMargin(0.11);// 5% della larghezza a destra
  		Scatola->SetTopMargin(0.07); //3% della altezza lato superiore
  		Scatola->SetBottomMargin(0.14); //12% dell'altezza lato inferiore
  		Scatola->SetTickx(1);
  		Scatola->SetTicky(1);

   		gr1->SetTitle(("Wave " + ampname + " channel " + channame).c_str());
		gr1->GetYaxis()->SetTitleSize(mytextsize); //controllo sulla dimension del titolo dell'asse
  		gr1->GetXaxis()->SetTitleSize(mytextsize);
  		gr1->GetXaxis()->SetLabelSize(mytextsize);//cotrollo sulla dimensione dei numeretti dell'asse
  		gr1->GetYaxis()->SetLabelSize(mytextsize);
  		gr1->GetXaxis()->SetTitleFont(mytextfont);//controllo sul carattere usato per il titolo dell'asse
  		gr1->GetYaxis()->SetTitleFont(mytextfont);
  		gr1->GetXaxis()->SetLabelFont(mytextfont);//controllo sul carattere usato per i numeretti dell'asse
  		gr1->GetYaxis()->SetLabelFont(mytextfont);
  		gr1->GetXaxis()->SetNdivisions(908); //suddivisione dei numeri sull'asse x---es 0 a 10 a passo di 1, e ogni passo diviso in 5
  		gr1->GetXaxis()->CenterTitle(1);//che il titolo dell'asse lo voglio quindi 1, se non lo volessi metterei 0
  		gr1->GetYaxis()->CenterTitle(1);
  		gr1->GetXaxis()->SetTitleOffset(1.15);//definisce la distanza del titolo dell'asse dall'asse stesso
  		gr1->GetYaxis()->SetTitleOffset(1.15);
  		gr1->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  		gr1->GetYaxis()->SetTitle("#sigma (pb)");
  		gr1->SetMarkerColor(mymarkercolor);
  		gr1->SetMarkerSize(mymarkersize);
  		gr1->SetMarkerStyle(mymarkerstyle);
  		gr1->SetLineColor(1);
  		gr1->SetLineWidth(1);
		gr1->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		//gr1->GetYaxis()->SetRangeUser(0,450.);
		gr1->SetLineWidth(1);
		//gr1->SetLineColor(kMagenta);

		//gr2->SetTitle("TGraphErrors Example");
   		gr2->SetMarkerSize(0);
		gr2->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		//gr2->GetYaxis()->SetRangeUser(0, 0.2);
		gr2->SetLineWidth(1);
		gr2->SetLineColor(kRed);

		TFile file("pdf_folder.root", "recreate");
		
		gr1->Write();
		gr2->Write();
		gr1->Draw("AP");
		gr2->Draw("same");
		Scatola->SaveAs(("Plots/"+pdfname+".pdf").c_str());
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

		void plotInclCrossSecWithExp(string pdfname, double lower_bound, double upper_bound){
		
		int num_exp_pts = data_InclCrossSec.sqrts.size();

		double x[num_exp_pts], y[num_exp_pts], ex[num_exp_pts], ey[num_exp_pts];

		const int num_th_pts = 1000;
		double x2[num_th_pts], y2[num_th_pts], ex2[num_th_pts], ey2[num_th_pts];

		for(int i = 0; i < num_exp_pts; i++){

			x[i] = 0;
			y[i] = 0;
			ex[i] = 0;
			ey[i] = 0;

		}

		double lb = data_InclCrossSec.sqrts[0];

		double ub = data_InclCrossSec.sqrts[num_exp_pts - 1];

		for (int i = 0; i < num_th_pts; i++){

			x2[i] = lb + (ub - lb) * i / ((double)num_th_pts - 1.);

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

		comp temp = 0;
		double aux = 0;

		for(int i = 0; i < num_th_pts; i++){

			aux = 0;

			for(string ampname : getAmpNames()){

				int amp_index = getampindex(ampname);
				amplitude amp = amplitudes[amp_index];

				for(string channame : amp.getChanNames()){

					int chan_index = getchanindex(ampname, channame);
					temp = amp.getValue(pow(x2[i],2))(chan_index);
					aux += (temp*conj(temp)).real();

				}

			}

			val = aux;
			if(isnan(val)) y2[i] = 0;
			else y2[i] = val;

		}

		auto gr1 = new TGraphErrors(num_exp_pts,x,y,ex,ey);
		auto gr2 = new TGraph(num_th_pts,x2,y2);

		int mymarkerstyle=20;
  		float mymarkersize=1.;
  		int mymarkercolor = 4;
  		float mytextsize=0.04;
  		int mytextfont=132;

		TCanvas *Scatola = new TCanvas("Scatola","Scatola",600,500); //costruttore 600pt x 550 pt
  		//gStyle->SetOptStat(0); //non voglio che mi metti il riquadro con la statistica
  		Scatola->SetFillColor(0);//il fondo del grafico con 0 è bianco...in teoria lo potete cambiare
  		Scatola->SetBorderMode(0);//mette dei bordi attorno alla figura...0 nessun bordo
  		Scatola->SetBorderSize(2); //spessore del bordo
  		Scatola->SetLeftMargin(0.18); //spazio a sinistra della figura ...20% della larghezza
  		Scatola->SetRightMargin(0.11);// 5% della larghezza a destra
  		Scatola->SetTopMargin(0.07); //3% della altezza lato superiore
  		Scatola->SetBottomMargin(0.14); //12% dell'altezza lato inferiore
  		Scatola->SetTickx(1);
  		Scatola->SetTicky(1);

		//gr1->SetTitle("TGraphErrors Example");
		gr1->GetYaxis()->SetTitleSize(mytextsize); //controllo sulla dimension del titolo dell'asse
  		gr1->GetXaxis()->SetTitleSize(mytextsize);
  		gr1->GetXaxis()->SetLabelSize(mytextsize);//cotrollo sulla dimensione dei numeretti dell'asse
  		gr1->GetYaxis()->SetLabelSize(mytextsize);
  		gr1->GetXaxis()->SetTitleFont(mytextfont);//controllo sul carattere usato per il titolo dell'asse
  		gr1->GetYaxis()->SetTitleFont(mytextfont);
  		gr1->GetXaxis()->SetLabelFont(mytextfont);//controllo sul carattere usato per i numeretti dell'asse
  		gr1->GetYaxis()->SetLabelFont(mytextfont);
  		gr1->GetXaxis()->SetNdivisions(908); //suddivisione dei numeri sull'asse x---es 0 a 10 a passo di 1, e ogni passo diviso in 5
  		gr1->GetXaxis()->CenterTitle(1);//che il titolo dell'asse lo voglio quindi 1, se non lo volessi metterei 0
  		gr1->GetYaxis()->CenterTitle(1);
  		gr1->GetXaxis()->SetTitleOffset(1.15);//definisce la distanza del titolo dell'asse dall'asse stesso
  		gr1->GetYaxis()->SetTitleOffset(1.15);
  		gr1->SetTitle("");
  		gr1->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  		gr1->GetYaxis()->SetTitle("#sigma (pb)");
  		gr1->SetMarkerColor(mymarkercolor);
  		gr1->SetMarkerSize(mymarkersize);
  		gr1->SetMarkerStyle(mymarkerstyle);
  		gr1->SetLineColor(1);
  		gr1->SetLineWidth(1);
		gr1->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		gr1->GetYaxis()->SetRangeUser(0.,450.);
		gr1->SetLineWidth(1);
		//gr1->SetLineColor(kMagenta);

		//gr2->SetTitle("TGraphErrors Example");
   		gr2->SetMarkerSize(0);
		gr2->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		//gr2->GetYaxis()->SetRangeUser(0, 0.2);
		gr2->SetLineWidth(1);
		gr2->SetLineColor(kRed);

		TFile file("pdf_folder.root", "recreate");
		
		gr1->Write();
		gr2->Write();
		gr1->Draw("AP");
		gr2->Draw("same");
		Scatola->SaveAs(("Plots/"+pdfname+".pdf").c_str());
		file.Close();
		return;

	}

	void plotInclCrossSecVsSumOfExcl(string pdfname, double lower_bound, double upper_bound){

		int totnumofmeasuredchans = data.size();
		
		int num_excl_pts = data[0].amp_expval.size();
		int num_Bsstar_pts = data[getchanindex("P", "B_sstarB_sstar")].amp_expval.size();
		//cout << getchanindex("P", "B_sstarB_sstar") << endl;
		int num_incl_pts = data_InclCrossSec.amp_expval.size();

		double test1 = 0;
		double test2 = 0;

		for(int i = 0; i < num_excl_pts; i++){
			for(int j = 0; j < 3; j++) test1 += 1./(pow(data[j].amp_expval_sist_err[i], 2) + pow(data[j].amp_expval_stat_err[i], 2));
		}

		test2 = test1;

		for(int i = 0; i < num_Bsstar_pts; i++){
			int index = getchanindex("P", "B_sstarB_sstar");
			//cout << 1./pow(data[index].amp_expval_sist_err[i], 2) + pow(data[index].amp_expval_stat_err[i], 2) << endl;
			test2 += 1./pow(data[index].amp_expval_sist_err[i], 2) + pow(data[index].amp_expval_stat_err[i], 2);
		}

		//cout << num_excl_pts << "	" << num_incl_pts << "	" << test1 << "	" << test2 << endl; 

		//cout << num_excl_pts << " " << num_Bsstar_pts << " " << num_incl_pts << endl; //I have
		//to exclude B_sstarB_sstar channel from gr1 because it has a different number of points
		//with respect of the other exclusive sigma!!! 

		double x1[num_excl_pts], y1[num_excl_pts], ex1[num_excl_pts], ey1[num_excl_pts];
		double x2[num_Bsstar_pts], y2[num_Bsstar_pts], ex2[num_Bsstar_pts], ey2[num_Bsstar_pts];
		double x3[num_incl_pts], y3[num_incl_pts], ex3[num_incl_pts], ey3[num_incl_pts];
		double y4[num_incl_pts], ey4[num_incl_pts];

		double val = 0;

		for(int i = 0; i < num_excl_pts; i++){

			int index = getchanindex("P", "B_sstarB_sstar");

			val = data[0].sqrts[i];
			if(isnan(val)) x1[i] = 0;
			else x1[i] = val;
			ex1[i] = 0;

			val = 0;
			for(int j = 0; j < 3; j++) val += data[j].amp_expval[i];
			if(isnan(val)) y1[i] = 0;
			else y1[i] = val;

			/*int tempdiff = num_excl_pts - num_Bsstar_pts;
			cout << tempdiff << endl; 
			val = 0;
			if(i > tempdiff){cout << i << endl;exit(0);
				val += data[index].amp_expval[i - tempdiff];
				cout << data[index].amp_expval[i - tempdiff] << endl;
				if(isnan(val)) y4[i] = y1[i];
				else y4[i] = val;
			}*/


			val = 0;
			for(int j = 0; j < 3; j++) val += pow(data[j].amp_expval_stat_err[i], 2) + pow(data[j].amp_expval_sist_err[i], 2);
			val = sqrt(val);
			if(isnan(val)) ey1[i] = 0;
			else ey1[i] = val;
			/*if(i > tempdiff){
				val += sqrt(pow(data[index].amp_expval_stat_err[i - tempdiff], 2) + pow(data[index].amp_expval_sist_err[i - tempdiff], 2));
				if(isnan(val)) y4[i] = 0;
				else y4[i] = val;
			}*/

		}

		for(int i = 0; i < num_Bsstar_pts; i++){

			int index = getchanindex("P", "B_sstarB_sstar");

			val = data[index].sqrts[i];
			if(isnan(val)) x2[i] = 0;
			else x2[i] = val;
			//cout << val << endl;
			ex2[i] = 0;
			//cout << val << " ";

			val = data[index].amp_expval[i];
			if(isnan(val)) y2[i] = 0;
			else y2[i] = val;
			//cout << val << " ";


			val = pow(data[index].amp_expval_stat_err[i], 2);
			val += pow(data[index].amp_expval_sist_err[i], 2);
			val = sqrt(val);
			if(isnan(val)) ey2[i] = 0;
			else ey2[i] = val;
			//cout << val << endl;

		}

		int tempdiff = num_excl_pts - num_Bsstar_pts;

		for(int i = 0; i < num_excl_pts; i++){
			y4[i] = y1[i];
			ey4[i] = ey1[i];

			for(int j = 0; j < num_Bsstar_pts; j++){
				if(abs(x1[i] - x2[j]) <= 0.001){
					y4[i] += y2[j];
					//ey4[i] += ey2[j]; 
					ey4[i] = sqrt(pow(ey1[i], 2) + pow(ey2[j], 2));
				}
			}

		}

		for(int i = 0; i < num_incl_pts; i++){

			val = data_InclCrossSec.sqrts[i];
			if(isnan(val)) x3[i] = 0;
			else x3[i] = val;
			ex3[i] = 0;

			val += data_InclCrossSec.amp_expval[i];
			if(isnan(val)) y3[i] = 0;
			else y3[i] = val;

			/*val = pow(data_InclCrossSec.amp_expval_stat_err[i], 2);
			val += pow(data_InclCrossSec.amp_expval_sist_err[i], 2);
			val = pow(val, 0.5);*/
			val = data_InclCrossSec.amp_expval_stat_err[i];
			if(isnan(val)) ey3[i] = 0;
			else ey3[i] = val;

		}

		const int num_th_pts = 1000;
		int num_exp_pts = data_InclCrossSec.sqrts.size();

		double x6[num_th_pts], y6[num_th_pts];

		double lb = data_InclCrossSec.sqrts[0];

		double ub = data_InclCrossSec.sqrts[num_exp_pts - 1];

		for (int i = 0; i < num_th_pts; i++){

			x6[i] = lb + (ub - lb) * i / ((double)num_th_pts - 1.);

		}

		comp temp = 0;
		double aux = 0;

		for(int i = 0; i < num_th_pts; i++){

			aux = 0;

			for(string ampname : getAmpNames()){

				int amp_index = getampindex(ampname);
				amplitude amp = amplitudes[amp_index];

				for(string channame : amp.getChanNames()){

					int chan_index = getchanindex(ampname, channame);
					temp = amp.getValue(pow(x6[i],2))(chan_index);
					aux += (temp*conj(temp)).real();

				}

			}

			val = aux;
			if(isnan(val)) y6[i] = 0;
			else y6[i] = val;

		}

		auto gr1 = new TGraphErrors(num_excl_pts,x1,y1,ex1,ey1);
		auto gr2 = new TGraphErrors(num_Bsstar_pts,x2,y2,ex2,ey2);
		auto gr3 = new TGraphErrors(num_incl_pts,x3,y3,ex3,ey3);
		auto gr4 = new TGraphErrors(num_excl_pts,x1,y4,ex1,ey4);
		auto gr5 = new TGraphErrors("Data/RED_incl_with_weight.txt");
		auto gr6 = new TGraph(num_th_pts,x6,y6);

		int mymarkerstyle=20;
  		float mymarkersize=1.;
  		int mymarkercolor = 4;
  		float mytextsize=0.04;
  		int mytextfont=132;

		TCanvas *Scatola = new TCanvas("Scatola","Scatola",600,500); //costruttore 600pt x 550 pt
  		//gStyle->SetOptStat(0); //non voglio che mi metti il riquadro con la statistica
  		Scatola->SetFillColor(0);//il fondo del grafico con 0 è bianco...in teoria lo potete cambiare
  		Scatola->SetBorderMode(0);//mette dei bordi attorno alla figura...0 nessun bordo
  		Scatola->SetBorderSize(2); //spessore del bordo
  		Scatola->SetLeftMargin(0.18); //spazio a sinistra della figura ...20% della larghezza
  		Scatola->SetRightMargin(0.11);// 5% della larghezza a destra
  		Scatola->SetTopMargin(0.07); //3% della altezza lato superiore
  		Scatola->SetBottomMargin(0.14); //12% dell'altezza lato inferiore
  		Scatola->SetTickx(1);
  		Scatola->SetTicky(1);

   		//gr1->SetTitle("TGraphErrors Example");
		gr1->GetYaxis()->SetTitleSize(mytextsize); //controllo sulla dimension del titolo dell'asse
  		gr1->GetXaxis()->SetTitleSize(mytextsize);
  		gr1->GetXaxis()->SetLabelSize(mytextsize);//cotrollo sulla dimensione dei numeretti dell'asse
  		gr1->GetYaxis()->SetLabelSize(mytextsize);
  		gr1->GetXaxis()->SetTitleFont(mytextfont);//controllo sul carattere usato per il titolo dell'asse
  		gr1->GetYaxis()->SetTitleFont(mytextfont);
  		gr1->GetXaxis()->SetLabelFont(mytextfont);//controllo sul carattere usato per i numeretti dell'asse
  		gr1->GetYaxis()->SetLabelFont(mytextfont);
  		gr1->GetXaxis()->SetNdivisions(908); //suddivisione dei numeri sull'asse x---es 0 a 10 a passo di 1, e ogni passo diviso in 5
  		gr1->GetXaxis()->CenterTitle(1);//che il titolo dell'asse lo voglio quindi 1, se non lo volessi metterei 0
  		gr1->GetYaxis()->CenterTitle(1);
  		gr1->GetXaxis()->SetTitleOffset(1.15);//definisce la distanza del titolo dell'asse dall'asse stesso
  		gr1->GetYaxis()->SetTitleOffset(1.15);
  		gr1->SetTitle("");
  		gr1->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  		gr1->GetYaxis()->SetTitle("#sigma (pb)");
  		gr1->SetMarkerColor(kMagenta);
  		gr1->SetMarkerSize(mymarkersize);
  		gr1->SetMarkerStyle(mymarkerstyle);
  		gr1->SetLineColor(1);
  		gr1->SetLineWidth(0);
		gr1->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		gr1->GetYaxis()->SetRangeUser(0,450.);
		gr1->SetLineWidth(1);
		gr1->SetLineColor(kMagenta);

		//gr2->SetTitle("TGraphErrors Example");
		gr3->GetYaxis()->SetTitleSize(mytextsize); //controllo sulla dimension del titolo dell'asse
  		gr3->GetXaxis()->SetTitleSize(mytextsize);
  		gr3->GetXaxis()->SetLabelSize(mytextsize);//cotrollo sulla dimensione dei numeretti dell'asse
  		gr3->GetYaxis()->SetLabelSize(mytextsize);
  		gr3->GetXaxis()->SetTitleFont(mytextfont);//controllo sul carattere usato per il titolo dell'asse
  		gr3->GetYaxis()->SetTitleFont(mytextfont);
  		gr3->GetXaxis()->SetLabelFont(mytextfont);//controllo sul carattere usato per i numeretti dell'asse
  		gr3->GetYaxis()->SetLabelFont(mytextfont);
  		gr3->GetXaxis()->SetNdivisions(908); //suddivisione dei numeri sull'asse x---es 0 a 10 a passo di 1, e ogni passo diviso in 5
  		gr3->GetXaxis()->CenterTitle(1);//che il titolo dell'asse lo voglio quindi 1, se non lo volessi metterei 0
  		gr3->GetYaxis()->CenterTitle(1);
  		gr3->GetXaxis()->SetTitleOffset(1.15);//definisce la distanza del titolo dell'asse dall'asse stesso
  		gr3->GetYaxis()->SetTitleOffset(1.15);
  		gr3->SetTitle("");
  		gr3->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  		gr3->GetYaxis()->SetTitle("#sigma (pb)");
  		gr3->SetMarkerColor(mymarkercolor);
  		gr3->SetMarkerSize(mymarkersize);
  		gr3->SetMarkerStyle(mymarkerstyle);
		gr3->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		gr3->GetYaxis()->SetRangeUser(0,450.);
		gr3->SetLineWidth(1);
		gr3->SetLineColor(mymarkercolor);
		gr3->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		//gr2->GetYaxis()->SetRangeUser(0, 0.2);
		gr3->SetLineWidth(1);


		//gr2->SetTitle("TGraphErrors Example");
   		gr2->SetMarkerColor(kGreen);
		gr2->SetMarkerStyle(21);
		gr2->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		//gr3->GetYaxis()->SetRangeUser(0, 0.2);
		gr2->SetLineWidth(1);
		gr2->SetLineColor(kGreen);

		//gr4->SetTitle("TGraphErrors Example");
		gr4->GetYaxis()->SetTitleSize(mytextsize); //controllo sulla dimension del titolo dell'asse
  		gr4->GetXaxis()->SetTitleSize(mytextsize);
  		gr4->GetXaxis()->SetLabelSize(mytextsize);//cotrollo sulla dimensione dei numeretti dell'asse
  		gr4->GetYaxis()->SetLabelSize(mytextsize);
  		gr4->GetXaxis()->SetTitleFont(mytextfont);//controllo sul carattere usato per il titolo dell'asse
  		gr4->GetYaxis()->SetTitleFont(mytextfont);
  		gr4->GetXaxis()->SetLabelFont(mytextfont);//controllo sul carattere usato per i numeretti dell'asse
  		gr4->GetYaxis()->SetLabelFont(mytextfont);
  		gr4->GetXaxis()->SetNdivisions(908); //suddivisione dei numeri sull'asse x---es 0 a 10 a passo di 1, e ogni passo diviso in 5
  		gr4->GetXaxis()->CenterTitle(1);//che il titolo dell'asse lo voglio quindi 1, se non lo volessi metterei 0
  		gr4->GetYaxis()->CenterTitle(1);
  		gr4->GetXaxis()->SetTitleOffset(1.15);//definisce la distanza del titolo dell'asse dall'asse stesso
  		gr4->GetYaxis()->SetTitleOffset(1.15);
  		gr4->SetTitle("");
  		gr4->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  		gr4->GetYaxis()->SetTitle("#sigma (pb)");
  		gr4->SetMarkerSize(mymarkersize);
  		gr4->SetMarkerStyle(mymarkerstyle);
   		gr4->SetMarkerColor(kRed);
		gr4->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		gr4->GetYaxis()->SetRangeUser(0,450.);
		gr4->SetLineWidth(1);
		gr4->SetLineColor(kRed);
		gr4->SetMarkerStyle(21);
		gr4->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		//gr4->GetYaxis()->SetRangeUser(0, 0.2);
		gr4->SetLineWidth(1);

		//gr4->SetTitle("TGraphErrors Example");
		gr5->GetYaxis()->SetTitleSize(mytextsize); //controllo sulla dimension del titolo dell'asse
  		gr5->GetXaxis()->SetTitleSize(mytextsize);
  		gr5->GetXaxis()->SetLabelSize(mytextsize);//cotrollo sulla dimensione dei numeretti dell'asse
  		gr5->GetYaxis()->SetLabelSize(mytextsize);
  		gr5->GetXaxis()->SetTitleFont(mytextfont);//controllo sul carattere usato per il titolo dell'asse
  		gr5->GetYaxis()->SetTitleFont(mytextfont);
  		gr5->GetXaxis()->SetLabelFont(mytextfont);//controllo sul carattere usato per i numeretti dell'asse
  		gr5->GetYaxis()->SetLabelFont(mytextfont);
  		gr5->GetXaxis()->SetNdivisions(908); //suddivisione dei numeri sull'asse x---es 0 a 10 a passo di 1, e ogni passo diviso in 5
  		gr5->GetXaxis()->CenterTitle(1);//che il titolo dell'asse lo voglio quindi 1, se non lo volessi metterei 0
  		gr5->GetYaxis()->CenterTitle(1);
  		gr5->GetXaxis()->SetTitleOffset(1.15);//definisce la distanza del titolo dell'asse dall'asse stesso
  		gr5->GetYaxis()->SetTitleOffset(1.15);
  		gr5->SetTitle("");
  		gr5->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  		gr5->GetYaxis()->SetTitle("#sigma (pb)");
  		gr5->SetMarkerSize(mymarkersize);
  		gr5->SetMarkerStyle(mymarkerstyle);
   		gr5->SetMarkerColor(kBlack);
		gr5->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		gr5->GetYaxis()->SetRangeUser(0,450.);
		gr5->SetLineWidth(1);
		gr5->SetLineColor(kBlack);
		gr5->SetMarkerStyle(21);
		gr5->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		//gr4->GetYaxis()->SetRangeUser(0, 0.2);
		gr5->SetLineWidth(1);

		//gr2->SetTitle("TGraphErrors Example");
   		gr6->SetMarkerSize(0);
		gr6->GetXaxis()->SetRangeUser(lower_bound, upper_bound);
		//6r2->GetYaxis()->SetRangeUser(0, 0.2);
		gr6->SetLineWidth(1);
		gr6->SetLineColor(kRed);

		TFile file("pdf_folder.root", "recreate");
		
		gr1->Write();
		gr2->Write();
		gr3->Write();
		gr4->Write();
		gr5->Write();
		gr6->Write();
		//gr3->Draw("AP");
		gr5->Draw("AP");
		//gr2->Draw("Psame");
		//gr1->Draw("Psame");
		gr4->Draw("Psame");
		gr6->Draw("same");
		//gr3->Draw("Psame");
		//gr5->Draw("Psame");
		Scatola->SaveAs(("Plots/"+pdfname+".pdf").c_str());
		file.Close();
		return;

	}


	void PolePlotGraph2D(string inputfile, string polefile, string sheet, double Re_sx, double Re_dx, double Im_sx, double Im_dx){

		string temp = "Plots/" + inputfile + "_poles_graph2D.pdf";
		if (sheet == "true") temp = "Plots/" + inputfile + "_poles_graph2D_II.pdf";

		TCanvas *box = new TCanvas("box","box",600,500); //costruttore 600pt x 550 pt
  		//gStyle->SetOptStat(0); //non voglio che mi metti il riquadro con la statistica
  		box->SetFillColor(0);//il fondo del grafico con 0 è bianco...in teoria lo potete cambiare
  		box->SetBorderMode(0);//mette dei bordi attorno alla figura...0 nessun bordo
  		box->SetBorderSize(2); //spessore del bordo
  		box->SetLeftMargin(0.18); //spazio a sinistra della figura ...20% della larghezza
  		box->SetRightMargin(0.11);// 5% della larghezza a destra
  		box->SetTopMargin(0.07); //3% della altezza lato superiore
  		box->SetBottomMargin(0.14); //12% dell'altezza lato inferiore
  		box->SetTickx(1);
  		box->SetTicky(1);
		
		TGraph2D *gr = new TGraph2D(polefile.c_str());

		int mymarkerstyle=20;
  		float mymarkersize=1.;
  		int mymarkercolor = 4;
  		float mytextsize=0.05;
  		int mytextfont=132;

		gr->SetTitle("");
		gr->SetMarkerColor(mymarkercolor);
  		gr->SetMarkerSize(mymarkersize);
  		gr->SetMarkerStyle(mymarkerstyle);
  		gr->SetLineColor(1);
  		gr->SetLineWidth(0);

		gr->GetYaxis()->SetTitleSize(mytextsize); //controllo sulla dimension del titolo dell'asse
  		gr->GetXaxis()->SetTitleSize(mytextsize);
		gr->GetZaxis()->SetTitleSize(mytextsize);
  		gr->GetXaxis()->SetLabelSize(mytextsize);//cotrollo sulla dimensione dei numeretti dell'asse
  		gr->GetYaxis()->SetLabelSize(mytextsize);
		gr->GetZaxis()->SetLabelSize(mytextsize);	
  		gr->GetXaxis()->SetTitleFont(mytextfont);//controllo sul carattere usato per il titolo dell'asse
  		gr->GetYaxis()->SetTitleFont(mytextfont);
		gr->GetZaxis()->SetTitleFont(mytextfont);
  		gr->GetXaxis()->SetLabelFont(mytextfont);//controllo sul carattere usato per i numeretti dell'asse
  		gr->GetYaxis()->SetLabelFont(mytextfont);
		gr->GetZaxis()->SetLabelFont(mytextfont);
  		gr->GetXaxis()->SetNdivisions(908); //suddivisione dei numeri sull'asse x---es 0 a 10 a passo di 1, e ogni passo diviso in 5
  		
		gr->GetXaxis()->CenterTitle(1);//che il titolo dell'asse lo voglio quindi 1, se non lo volessi metterei 0
  		gr->GetYaxis()->CenterTitle(1);
  		gr->GetZaxis()->CenterTitle(1);
		gr->GetXaxis()->SetTitleOffset(1.20);//definisce la distanza del titolo dell'asse dall'asse stesso
  		gr->GetYaxis()->SetTitleOffset(1.70);
  		gr->GetZaxis()->SetTitleOffset(1.15);
		
  		gr->GetXaxis()->SetTitle("Re(s) (GeV^{2})");
  		gr->GetYaxis()->SetTitle("Im(s) (GeV^{2})");
  		gr->GetZaxis()->SetTitle("log_{10}|det(D_{l}(s))|");
		
		gr->GetXaxis()->SetLimits(Re_sx, Re_dx);
		gr->GetYaxis()->SetLimits(Im_sx, 0);

		gr->Draw("p");

		box->SaveAs(temp.c_str());
		
	};

	void PolePlotGraph1D(string inputfile, string polefile, string sheet, double Re_sx, double Re_dx, double Im_sx, double Im_dx){

		string temp = "Plots/" + inputfile + "_poles_graph.pdf";
		if (sheet == "true") temp = "Plots/" + inputfile + "_poles_graph_II.pdf";

		TCanvas *box = new TCanvas("box","box",600,500); //costruttore 600pt x 550 pt
  		//gStyle->SetOptStat(0); //non voglio che mi metti il riquadro con la statistica
  		box->SetFillColor(0);//il fondo del grafico con 0 è bianco...in teoria lo potete cambiare
  		box->SetBorderMode(0);//mette dei bordi attorno alla figura...0 nessun bordo
  		box->SetBorderSize(2); //spessore del bordo
  		box->SetLeftMargin(0.18); //spazio a sinistra della figura ...20% della larghezza
  		box->SetRightMargin(0.11);// 5% della larghezza a destra
  		box->SetTopMargin(0.07); //3% della altezza lato superiore
  		box->SetBottomMargin(0.14); //12% dell'altezza lato inferiore
  		box->SetTickx(1);
  		box->SetTicky(1);
		
		TGraph *gr = new TGraph(polefile.c_str());

		int mymarkerstyle=20;
  		float mymarkersize=1.;
  		int mymarkercolor = 4;
  		float mytextsize=0.05;
  		int mytextfont=132;

		gr->GetYaxis()->SetTitleSize(mytextsize); //controllo sulla dimension del titolo dell'asse
  		gr->GetXaxis()->SetTitleSize(mytextsize);
  		gr->GetXaxis()->SetLabelSize(mytextsize);//cotrollo sulla dimensione dei numeretti dell'asse
  		gr->GetYaxis()->SetLabelSize(mytextsize);
  		gr->GetXaxis()->SetTitleFont(mytextfont);//controllo sul carattere usato per il titolo dell'asse
  		gr->GetYaxis()->SetTitleFont(mytextfont);
  		gr->GetXaxis()->SetLabelFont(mytextfont);//controllo sul carattere usato per i numeretti dell'asse
  		gr->GetYaxis()->SetLabelFont(mytextfont);

  		gr->GetXaxis()->SetNdivisions(908); //suddivisione dei numeri sull'asse x---es 0 a 10 a passo di 1, e ogni passo diviso in 5
  		
		gr->GetXaxis()->CenterTitle(1);//che il titolo dell'asse lo voglio quindi 1, se non lo volessi metterei 0
  		gr->GetYaxis()->CenterTitle(1);
		gr->GetXaxis()->SetTitleOffset(1.15);//definisce la distanza del titolo dell'asse dall'asse stesso
  		gr->GetYaxis()->SetTitleOffset(1.15);
		gr->SetTitle("");
  		gr->GetXaxis()->SetTitle("Re(s) (GeV^{2})");
  		gr->GetYaxis()->SetTitle("Im(s) (GeV^{2})");
		gr->GetXaxis()->SetRangeUser(Re_sx, Re_dx);
		gr->GetYaxis()->SetRangeUser(Im_sx, 0);
		gr->SetMarkerColor(mymarkercolor);
  		gr->SetMarkerSize(mymarkersize);
  		gr->SetMarkerStyle(mymarkerstyle);
  		gr->SetLineColor(1);
  		gr->SetLineWidth(0);
		gr->Draw("AP");

		box->SaveAs(temp.c_str());
		
	};

	void PoleColormapPlotFunc2D(string inputfile, function <double(double*,double*)> func, string funcname, double Re_sx, double Re_dx, double Im_sx, double Im_dx){

		string temp = "Plots/" + inputfile + "_" + funcname + ".pdf";

		TCanvas *box = new TCanvas("box","box",600,500); //costruttore 600pt x 550 pt
  		//gStyle->SetOptStat(0); //non voglio che mi metti il riquadro con la statistica
  		box->SetFillColor(0);//il fondo del grafico con 0 è bianco...in teoria lo potete cambiare
  		box->SetBorderMode(0);//mette dei bordi attorno alla figura...0 nessun bordo
  		box->SetBorderSize(2); //spessore del bordo
  		box->SetLeftMargin(0.18); //spazio a sinistra della figura ...20% della larghezza
  		box->SetRightMargin(0.11);// 5% della larghezza a destra
  		box->SetTopMargin(0.07); //3% della altezza lato superiore
  		box->SetBottomMargin(0.14); //12% dell'altezza lato inferiore
  		box->SetTickx(1);
  		box->SetTicky(1);

		//TF2 *tf = new TF2("tf", detD, 113, 121, -1, 1, 2);
		//TF2 tf("tf", [](double* x, double* p) { return abs(testObs.amplitudes[0].getDenominator(comp(x[0], x[1])).determinant()); }, 113., 121., -1., 1.);
		TF2 tf("tf", func, Re_sx, Re_dx, Im_sx, Im_dx, 1);
		
		int mymarkerstyle=20;
  		float mymarkersize=1.;
  		int mymarkercolor = 4;
  		float mytextsize=0.05;
  		int mytextfont=132;

		tf.GetYaxis()->SetTitleSize(mytextsize); //controllo sulla dimension del titolo dell'asse
  		tf.GetXaxis()->SetTitleSize(mytextsize);
  		tf.GetXaxis()->SetLabelSize(mytextsize);//cotrollo sulla dimensione dei numeretti dell'asse
  		tf.GetYaxis()->SetLabelSize(mytextsize);
  		tf.GetXaxis()->SetTitleFont(mytextfont);//controllo sul carattere usato per il titolo dell'asse
  		tf.GetYaxis()->SetTitleFont(mytextfont);
  		tf.GetXaxis()->SetLabelFont(mytextfont);//controllo sul carattere usato per i numeretti dell'asse
  		tf.GetYaxis()->SetLabelFont(mytextfont);
  		tf.GetXaxis()->SetNdivisions(908); //suddivisione dei numeri sull'asse x---es 0 a 10 a passo di 1, e ogni passo diviso in 5
		tf.GetXaxis()->CenterTitle(1);//che il titolo dell'asse lo voglio quindi 1, se non lo volessi metterei 0
  		tf.GetYaxis()->CenterTitle(1);
		tf.GetXaxis()->SetTitleOffset(1.15);//definisce la distanza del titolo dell'asse dall'asse stesso
  		tf.GetYaxis()->SetTitleOffset(1.15);
		tf.SetTitle("");
  		tf.GetXaxis()->SetTitle("Re(s) (GeV^{2})");
  		tf.GetYaxis()->SetTitle("Im(s) (GeV^{2})");
		tf.SetMarkerColor(mymarkercolor);
  		tf.SetMarkerSize(mymarkersize);
  		tf.SetMarkerStyle(mymarkerstyle);
  		tf.SetLineColor(1);
  		tf.SetLineWidth(0);

		tf.Draw("COLZ");
		box->SaveAs(temp.c_str());

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

	void setinclchi2weight(double temp){
		incl_weight = temp;
	};

	void setexclchi2weight(double temp){
		excl_weight = temp;
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

		return excl_weight * result;

	}


	double chisq_with_InclCrossSec(){

		double sum = 0;
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

		return incl_weight * sum;

	}

	friend ostream& operator<<(std::ostream& os, observable const& m){
		for(amplitude a: m.amplitudes) cout<<a<<endl<<endl;
		return os;
	};


};
