#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm> // std::min_element
#include <iterator>  // std::begin, std::end
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TCanvas.h"
#include <TF2.h>
#include "TH1D.h"

using namespace std;

struct exppts{//this has to be modify according to the format of the exp data file

		vector<double> sqrts;
		vector<double> amp_expval, amp_expval_stat_err, amp_expval_sist_err;

		exppts(vector<double> x, vector<double> y, vector<double> y_stat_err, vector<double> y_sist_err){
			sqrts = x;
			
			amp_expval = y;
			amp_expval_stat_err = y_stat_err;
			amp_expval_sist_err = y_sist_err;
		}

};

int main()
{

    ifstream letsread1("Data/inclusive-cross-section.dat");
    ifstream letsread2("Data/BB-nominal.dat");

    double a[5];
    double b[5];

    vector<double> av = {};
    vector<double> bv = {};

    vector<double> av_y = {};
    vector<double> bv_y = {};

    vector<double> av_y_stat_err = {};
    vector<double> bv_y_stat_err = {};

    vector<double> av_y_sist_err = {};
    vector<double> bv_y_sist_err = {};

    while (letsread1 >> a[0] >> a[1] >> a[2] >> a[3] >> a[4])
    {
        av.push_back(a[0]);
        av_y.push_back(a[1]);
        av_y_stat_err.push_back(a[2]);
        av_y_sist_err.push_back(a[3]);
    }

    while (letsread2 >> b[0] >> b[1] >> b[2] >> b[3] >> b[4])
    {
        bv.push_back(b[0]);
        bv_y.push_back(a[1]);
        bv_y_stat_err.push_back(a[2]);
        bv_y_sist_err.push_back(a[3]);
    }

    vector<double> av_temp[bv.size()];
    vector<double> av_temp_y[bv.size()];
    vector<double> av_temp_y_err[bv.size()];
    for(int i = 0; i < bv.size(); i++){
        av_temp[i] = {};
        av_temp_y[i] = {};
        av_temp_y_err[i] = {};
    }
    vector<double> av_new = bv;

    vector<double> av_def;
    vector<double> av_def_y;
    vector<double> av_def_err;
    vector<double> av_def_y_err;

    double startfit = 10.65;
    double endfit = 11.0208;

    for (int i = 0; i < av.size(); i++)
    {
        if (av[i] >= startfit && av[i] <= endfit)
        {
            vector<double> diff = {};

            for (int j = 0; j < bv.size(); j++)
            {
                diff.push_back(fabs(av[i] - bv[j]));
            }
            auto it = min_element(begin(diff), end(diff));

            int min_index = distance(begin(diff), it);

            av_temp[min_index].push_back(av[i]);
            av_temp_y[min_index].push_back(av_y[i]);
            av_temp_y_err[min_index].push_back(sqrt(pow(av_y_sist_err[i], 2) + pow(av_y_stat_err[i], 2)));
            //av_temp_y_err[min_index].push_back(1);

        }
    }

    double sum = 0;

    /*for (int i = 0; i < bv.size(); i++){
        for (int j = 0; j < av_temp[i].size(); j++){
            cout << av_temp[i].size() << "    " << av_temp[i].at(j) << "   " << av_temp_y[i].at(j) << "   " << av_temp_y_err[i].at(j) << endl;
        }
    }*/

    for (int i = 0; i < bv.size(); i++)
    {
        if(av_temp[i].size() > 0) 
        {
            av_def.push_back(bv[i]);

            sum = 0;
            for (int j = 0; j < av_temp[i].size(); j++) sum += 1./pow(av_temp_y_err[i].at(j), 2);
            double norm = sum;
            av_def_y_err.push_back(1./sqrt(norm));

            sum = 0;
            for(int j = 0; j < av_temp[i].size(); j++){
                sum += av_temp_y[i].at(j)/pow(av_temp_y_err[i].at(j), 2);
            }
            av_def_y.push_back(sum/norm);
        }
    }

    for (int i = 0; i < av_def.size(); i++)
    {
        cout << "   " << av_def[i] << "   " << av_def_y[i] << "   " << av_def_y_err[i] << "   " << endl; //<< 1./sqrt(av_temp[i].size()) << endl;
        av_def_err.push_back(0);
    }

    ofstream letwrite("Data/RED_incl_with_weight.txt");

    for (int i = 0; i < av_def.size(); i++)
    {
        letwrite << "   " << av_def[i] << "   " << av_def_y[i] << "   " << av_def_err[i] << "   " << sqrt(15) * av_def_y_err[i] << "   " << endl; 
    }

    /*auto result = std::min_element(std::begin(diff), std::end(diff));
    //if (std::end(diff)!=result)
    std::cout << *result << '\n';
    cout << distance(begin(diff),result) << endl;*/

    int npts = av_def.size();
    double *aa_def = &av_def[0];
    double *aa_def_y = &av_def_y[0];
    double *aa_def_err = &av_def_err[0];
    double *aa_def_y_err = &av_def_y_err[0];

    auto gr1 = new TGraphErrors(npts, aa_def, aa_def_y, aa_def_err, aa_def_y_err);

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
	gr1->GetXaxis()->SetRangeUser(startfit, endfit);
	gr1->GetYaxis()->SetRangeUser(0,450.);
	gr1->SetLineWidth(1);
	gr1->SetLineColor(mymarkercolor);
	//gr2->GetYaxis()->SetRangeUser(0, 0.2);
	gr1->SetLineWidth(1);

    gr1->Draw("AP");
    Scatola->SaveAs("Plots/norm_incl.pdf");
}

// g++ norm_incl_maker.cpp `root-config --cflags --glibs` && ./a.out