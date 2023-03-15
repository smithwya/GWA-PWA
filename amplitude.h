#pragma once
#include <iostream>
#include <vector>
#include <map>
#include "channel.h"
#include <Eigen/Dense>
#include "TF1.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;

struct intKey{

	comp s;
	int k;

	intKey(comp x, int y){
		s = x;
		k = y;
	}

	bool operator<(intKey const& other) const{

		if(k!=other.k){
			return k<other.k;
		}
		
		if(s.real()!=other.s.real()) return s.real()<other.s.real();

		return s.imag() < other.s.imag();
		
	}
};



class amplitude {

public:
	amplitude();
	amplitude(int j, double alpha, double sl, vector<channel> c, vector<MatrixXcd> kParams, vector<double> rmasses,double s0, double smin, double smax);
	amplitude(int J, double sl, double smin, double smax, vector<channel> chans);

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
	void calcIntegrals(vector<comp> slist,int k);
	vector<double> getParamList();
	void setParamList(vector<double> params);
	void setResMasses(vector<double> rm);
private:
	vector<channel> channels;
	vector<MatrixXcd> kParameters;
	vector<double> resmasses;
	vector<string> channel_names;
	map<intKey, comp> integralList;
	int numChannels;
	int J;
	double alpha;
	double sL;
	double s0;
	double smin;
	double smax;
	const double epsilon = 1e-3;
	
};