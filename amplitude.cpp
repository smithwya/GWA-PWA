#pragma once
#include <iostream>
#include<complex>
#include<algorithm>
#include<math.h>
#include<chrono>
#include<random>
#include<vector>
#include<fstream>
#include<string>
#include "channel.h"
#include "amplitude.h"
#include <Eigen/Dense>
#include "TF1.h"
#include "TMath.h"
#include "TH1.h"
#include "Math/Integrator.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/GSLIntegrator.h"
#include "Math/IFunction.h"


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;
//typdef comp param;

//amplitude constructor
amplitude::amplitude() {
	channels = {};
	kParameters = {};
	numChannels = 0;
	s0 = 0;
	smin = 0;
	smax = 0;
}

amplitude::amplitude(comp j, comp alp, comp ssl, vector<channel> chans, vector<MatrixXcd> kParams,vector<double> rmasses, comp ss0, comp ssmin, comp ssmax) {
	channels = chans;
	kParameters = kParams;
	numChannels = channels.size();
	J = j;
	alpha = alp;
	sL = ssl;
	s0 = ss0;
	smin = ssmin;
	smax = ssmax;
	resmasses =rmasses;
}

//calculate the nth chebyshev polynomial T_n(x) (probably ok to replace with some library)
comp amplitude::chebyshev(comp x, int n) {
	
	if (n == 0) return comp(1.0);
    if (n == 1) return x;
    if (n > 1) return 2. * x * chebyshev(x, n - 1) - chebyshev(x, n - 2);
	return comp(0,0);
};




//calculates omega_s
comp amplitude::omega_p(comp s) {

	return s/(s + s0);
};

//calculates omega_p
comp amplitude::omega_s(comp s) {

	return 2. * (s - smin)/(smax - smin) - 1.0;
};

//calculates omega_ps
comp amplitude::omega_ps(comp s) {

	return 2. * (omega_p(s) - omega_p(smin))/(omega_p(smax) - omega_p(smin)) - 1.;
};

	// return E_gamma * p_i * numerator * denominator.inverse()
VectorXcd amplitude::getValue(comp s) {
	VectorXcd a=VectorXcd::Zero(numChannels);

	double m_JPsi = 3.0969;

	comp Egamma = (pow(m_JPsi, 2) - s) / (2.0 * sqrt(s));
	a=getNumerator(s,3).transpose()*getDenominator(s).inverse();

	for(int i = 0; i < numChannels; i++){
		a(i)*=Egamma * channels[i].getMomentum(s);
	}

	return a;
}

comp amplitude::omega(comp s, int type){

	switch (type)
	{
	case 1: return omega_s(s);
	case 2: return omega_p(s);
	case 3: return omega_ps(s);
	}

	return 0;
}

//returns the numerator function for an amplitude, where each entry is the entry for the kth channel expressed in terms of chebyshev polynomials
//n_k = sum(a_n * T_n(omega(s)))
VectorXcd amplitude::getNumerator(comp s, int type){

	VectorXcd numerator = VectorXcd::Zero(numChannels);

	for(int k = 0; k< numChannels; k++){
		comp n_k = comp(0,0);
		int ncoeffs = channels[k].getChebyCoeffs().size();

		for(int n = 0; n<ncoeffs; n++){
			n_k+=channels[k].getChebyCoeff(n)*chebyshev(omega(s,type),n);
		}

		numerator[k]=n_k;
	}
	return numerator;
}

//calculates rhoN_ki(s') = delta_ki * (2p_i)^{2J+1}/(s'+sL)^{J+alpha}
comp amplitude::getRhoN(comp sprime,int k)
{
		comp x = pow(2.0*channels[k].getMomentum(sprime),2.0*J+1.0)/pow(sprime + sL,J+alpha);

	return x;
}

//overload later for other sheets
comp amplitude::getIntegrand(double sp,comp s, int k){

	return getRhoN(sp,k)*s/(sp*(sp-s-comp(0,1)*epsilon)*TMath::Pi());

}


comp amplitude::getIntegral(comp s,int k){
	//add parameter to select sheet, "+-++" etc
 	auto realIntegrand = [&](double sp)
    {
        return amplitude::getIntegrand(sp,s,k).real();
    };
	auto imagIntegrand = [&](double sp)
    {
        return amplitude::getIntegrand(sp,s,k).imag();
    };

	ROOT::Math::Functor1D re(realIntegrand);
	ROOT::Math::Functor1D im(imagIntegrand);

	ROOT::Math::Integrator intRe(re,ROOT::Math::IntegrationOneDim::kGAUSS,1.E-9,1.E-6);
	ROOT::Math::Integrator intIm(im,ROOT::Math::IntegrationOneDim::kGAUSS,1.E-9,1.E-6);

	double threshold = 4*pow(channels[k].getMass().real(),2);

	double realpart = intRe.IntegralUp(threshold);
	double imagpart = intIm.IntegralUp(threshold);
	

	return comp(realpart,imagpart);
}



MatrixXcd amplitude::getDenominator(comp s)
{

	MatrixXcd Kinv = getKMatrix(s).inverse();

	MatrixXcd M = MatrixXcd::Zero(numChannels,numChannels);

	for(int k = 0; k < numChannels; k++){

		for(int i = 0; i <= k; i++){
			
			M(k,i)=getIntegral(s,k);
			M(i,k) = M(k,i);
		}


	}

	return Kinv-M;
}

//integratoe f(function, lower bound, upper bound)

comp amplitude::getMomentum(int chan, comp s)
{
	comp mass = channels[chan].getMass();
	return 0.5 * sqrt(s - 4.0 * pow(mass, 2));
}

MatrixXcd amplitude::getKMatrix(comp s) {
	MatrixXcd kmat = MatrixXcd::Zero(numChannels, numChannels);

	for (int k = 0; k < numChannels; k++) {
		for (int i = 0; i <= k; i++) {
			comp tempMatrixTerm = 0;

			for (int R = 0; R < numChannels; R++) {
				tempMatrixTerm += channels[k].getCoupling(R) * channels[i].getCoupling(R) / (resmasses[R] - s);
			}

			kmat(k, i) += tempMatrixTerm;
			if(i!=k) kmat(i, k) += tempMatrixTerm;
		}
	}

	for (int j = 0; j < kParameters.size(); j++) {

		kmat += pow(s, j) * kParameters[j];

	}

	return kmat;
}

ostream& operator<<(ostream& os, amplitude const& m) {
	os << "J = " << m.J << ", alpha = " << m.alpha << ", sL = " << m.sL << endl<<endl;
	for (int i = 0; i < m.numChannels; i++) {
		os << "channel " << i << ": " << endl << m.channels[i] << endl<<endl;
	}
	os<<"resmasses: ";
	for (int i = 0; i < m.resmasses.size(); i++) {
		os << m.resmasses[i]<<" ";
	}
	os<<endl;
	os << endl << endl;
	os << "k-matrix parameters: " << endl;
	for (int i = 0; i < m.kParameters.size(); i++) {
		os << "s^" << i << " coefficient matrix = " << endl << m.kParameters[i] << endl;
	}
	os << endl;
	return os;

}