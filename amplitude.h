#pragma once
#include <iostream>
#include<vector>
#include "channel.h"
#include <Eigen/Dense>
#include "TF1.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;

class amplitude {

public:
	amplitude();
	amplitude(int j, double alpha, double sl, vector<channel> c, vector<MatrixXcd> kParams, vector<double> rmasses,double s0, double smin, double smax);
	amplitude(int J, double smin, double smax, vector<channel> chans);

	comp chebyshev(comp x, int n);
	comp omega_s(comp s);
	comp omega_p(comp s);
	comp omega_ps(comp s);
	VectorXcd getValue(comp s);
    comp omega(comp s, int type);
    VectorXcd getNumerator(comp s, int type);
    comp getRhoN(comp s, int k);
	MatrixXcd getDenominator(comp s);
	MatrixXcd getKMatrix(comp s);
	comp getMomentum(int particle, comp s);
	friend ostream& operator<<(std::ostream& os, amplitude const& m);
	comp getIntegral(comp s, int k);
	comp getIntegrand(double sp, comp s, int k);
	void setChebyCoeffs(string channel_name, int type, double s0, vector<double> coeffs);
	void setKParams(int power, vector<vector<double>> kparamlist);
	void addPole(double mass,vector<string> chan_names, vector<double> couplings);
	vector<string> getChanNames();
private:
	vector<channel> channels;
	vector<MatrixXcd> kParameters;
	vector<double> resmasses;
	vector<string> channel_names;
	int numChannels;
	int J;
	double alpha;
	double sL;
	double s0;
	double smin;
	double smax;
	const double epsilon = 1e-3;
	
};