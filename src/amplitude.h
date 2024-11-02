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
	amplitude(string ampName,int J, double sl, double smin, double smax, vector<channel> chans, string kmat, string rhoN);

	comp chebyshev(comp x, int n);
	comp legendreII(comp x, int n);
	comp omega_s(comp s);
	comp omega_p(comp s);
	comp omega_ps(comp s);
	VectorXcd getValue(comp s);
	VectorXcd getValueForPoles(comp s, int sheet);
    comp omega(comp s, int type);
    VectorXcd getNumerator(comp s, int type);
    comp getRhoN(comp s, int k, int sheet);
	MatrixXcd getDenominator(comp s, int sheet);
	MatrixXcd getKMatrix(comp s);
	comp getMomentum(int particle, comp s);
	comp getComplexMomentum(int particle, comp s);
	comp getTrueMomentum(int particle, comp s);
	friend ostream& operator<<(std::ostream& os, amplitude const& m);
	comp getIntegral(comp s, int k, int sh);
	comp getIntegral(double s, int k);
	comp getIntegrand(double sp, comp s, int k);
	void setChebyCoeffs(string channel_name, int type, double s0, vector<double> coeffs);
	void setChebyCoeffs(string channel_name, int type, double s0, vector<double> coeffs,vector<double> csteps);
	void setKParams(int power, vector<vector<double>> kparamlist);
	void setKParams(int power, vector<vector<double>> kparamlist,vector<double> ksteps);
	void addPole(double mass,vector<string> chan_names, vector<double> couplings);
	void addPole(double mass,double mass_step,vector<string> chan_names, vector<double> couplings, vector<double> steps);
	vector<string> getChanNames();
	int getNumOfChans();
	void calcIntegrals(vector<comp> slist,int k,int sheet);
	vector<double> getParamList();
	void setParamList(vector<double> params);
	void setResMasses(vector<double> rm);
	vector<double> getResMasses();
	vector<double> getResMassesSteps();
	vector<double> getFittedParamList();
	void setFittedParamList(vector<double> fittedParams);
	string getName();
	vector<channel> getChannels();
	vector<MatrixXcd> getkParameters();
	vector<double> getKSteps(int i);
	vector<double> getPoleSteps();
	vector<double> getStepSizes();
	vector<double> getFittedSteps();
	vector<double> getFitInterval();
	string getKMatType();
	int getJ();
	double getEpsilon();
	double getSl();

private:
	vector<channel> channels;
	vector<MatrixXcd> kParameters;
	vector<vector<double>> kParameters_steps;
	vector<double> resmasses;
	vector<double> resmasses_steps; 
	vector<string> channel_names;
	map<intKey, comp> integralList;
	vector<bool> fixedParamList;
	int numChannels;
	int J;
	double alpha;
	double sL;
	double s0;
	double smin;
	double smax;
	double epsilon;
	string name;

	string kmattype;
	string rhoNtype;
	
};