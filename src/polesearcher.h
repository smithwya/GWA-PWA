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
    public:

    double minfuncforpoles(const double *xx){
        
        int numamp = 0;
        double val = 0;

        MatrixXcd cmat = testObs.amplitudes.at(numamp).getDenominator(comp(xx[0], xx[1]));

        for(int i = 0; i < temppoles.size(); i++){
            cmat = cmat*(comp(xx[0],xx[1])-temppoles[i]);
            cmat = cmat*(comp(xx[0],xx[1])-temppoles[i]);
        }

        if(pow(xx[0], 2) + pow(xx[1], 2) <= 400){
            val = log(abs(cmat.determinant()));
        }
        else{
            val = log(abs(cmat.determinant())) + pow(xx[0], 2) + pow(xx[1], 2);
        }

        return val;

    }


}