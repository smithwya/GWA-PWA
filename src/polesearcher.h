#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "amplitude.h"
#include "channel.h"
#include "observable.h"
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

class polesearcher {

    private:
    observable testObs;
    vector<comp> poles = {};
    int wvindex;
    bool sh = false;

    public:

    int calls_counter_poles = 0;

    void settestObs(observable obs){
        testObs = obs;
    }

    observable gettestObs(){
        return testObs;
    }

    void setAmpIndex(string ampname){
        wvindex = testObs.getampindex(ampname);
    }

    int getAmpIndex(){
        return wvindex;
    }

    void setPoles(vector<comp> xx){
        poles = xx;
    }

    void SetSheet(bool sheet){
        sh = sheet;
    }

    bool GetSheet(){
        return sh;
    }

    double minfuncforpoles(vector<double> params){

        calls_counter_poles++;
    
        double val = 0; 

        MatrixXcd cmat = testObs.amplitudes.at(wvindex).getDenominator(comp(params[0], params[1]), sh); 
        //MatrixXcd cmat = (comp(params[0],params[1]) - comp(115,0.5)) * (comp(params[0],params[1]) - comp(118,0.7)) * (comp(params[0],params[1]) - comp(115,-0.5)) * (comp(params[0],params[1]) - comp(118,-0.7)) * MatrixXcd({{1}}); //mock example

        comp det = cmat.determinant();

        for(int i = 0; i < poles.size(); i++){
            det /= (comp(params[0],params[1])-poles[i]);
            det /= (comp(params[0],params[1])-conj(poles[i]));
        }

        if(pow(params[0] - 120., 2) + pow(params[1], 2) <= 40.){
            val = log10(abs(det));
        }
        else{
            val = log10(abs(det)) + pow(params[0] - 120., 2) + pow(params[1], 2) - 40.;
        }

            //cout << comp(params[0],params[1]) << ": val = " << val << endl;
        return val;

    }


};