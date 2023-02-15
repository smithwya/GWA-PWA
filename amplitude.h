#pragma once
#include <iostream>
#include<vector>
#include "channel.h"
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;

class amplitude {

public:
	amplitude();
	amplitude(comp j, comp alpha, comp sl, vector<channel> c, vector<MatrixXcd> kParams, comp s0, comp smin, comp smax);
	comp chebyshev(comp x, int n);
	comp omega_s(comp s);
	comp omega_p(comp s);
	comp omega_ps(comp s);
	comp getValue(int chan, comp s);
    comp omega(comp s, int type);
    VectorXcd getNumerator(comp s, int type);
    MatrixXcd getRhoN(comp s);
	MatrixXcd getDenominator(comp s);
	MatrixXcd getKMatrix(comp s);
	comp getMomentum(int particle, comp s);
	friend ostream& operator<<(std::ostream& os, amplitude const& m);

private:
	vector<channel> channels;
	vector<MatrixXcd> kParameters;
	int numChannels;
	comp J;
	comp alpha;
	comp sL;
	comp s0;
	comp smin;
	comp smax;
};