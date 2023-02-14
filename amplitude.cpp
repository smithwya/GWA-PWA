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

amplitude::amplitude(comp j, comp alp, comp ssl, vector<channel> chans, vector<MatrixXcd> kParams, comp ss0, comp ssmin, comp ssmax) {
	channels = chans;
	kParameters = kParams;
	numChannels = channels.size();
	J = j;
	alpha = alp;
	sL = ssl;
	s0 = ss0;
	smin = ssmin;
	smax = ssmax;
}

//calculate the nth chebyshev polynomial T_n(x) (probably ok to replace with some library)
comp amplitude::chebyshev(comp x, int n) {
	
	if (n == 0) return comp(1.0);
    if (n == 1) return x;
    if (n > 1) return 2. * x * chebyshev(x, n - 1) - chebyshev(x, n - 2);
	return comp(0,0);
};


//calculates omega_s
comp amplitude::omega_s(comp s) {

	return s/(s + s0);
};

//calculates omega_p
comp amplitude::omega_p(comp s) {

	return 2. * (s - smin)/(smax - smin) - 1.;
};

//calculates omega_ps
comp amplitude::omega_ps(comp s) {

	return 2. * (omega_p(s) - omega_p(smin))/(omega_p(smin) - omega_p(smax)) - 1.;
};

comp amplitude::getValue(int chan, comp s) {
	//TODO: return E_gamma * p_i * numerator * denominator.inverse()
	return comp(1, 0);
}


	VectorXcd amplitude::getNumerator(comp s, int type)
{
	VectorXcd numerator = VectorXcd::Zero(numChannels);
	return numerator;
}





MatrixXcd amplitude::getRhoN(comp sprime)
{
	//TODO: calculate delta_ki (2p_i)^2J+1/(sprime + sL)^2J + alpha
	return MatrixXcd();
}




MatrixXcd amplitude::getDenominator(comp s)
{
	MatrixXcd D = MatrixXcd::Zero(numChannels, numChannels);

	MatrixXcd Kinv = getKMatrix(s).inverse();
	//TODO: calculate -s/pi * int(rho N/s'(s'-s-iepsilon),4mk^2, infty)

	return D;
}


comp amplitude::getMomentum(int chan, comp s)
{
	comp mass = channels[chan].getMass();
	return sqrt(s - 4.0 * pow(mass, 2));
}

MatrixXcd amplitude::getKMatrix(comp s) {
	MatrixXcd kmat = MatrixXcd::Zero(numChannels, numChannels);

	for (int k = 0; k < numChannels; k++) {
		for (int i = 0; i <= k; i++) {
			comp tempMatrixTerm = 0;

			for (int R = 0; R < numChannels; R++) {
				tempMatrixTerm += channels[k].getCoupling(R) * channels[i].getCoupling(R) / (pow(channels[R].getMass(), 2) - s);
			}

			kmat(k, i) += tempMatrixTerm;
			kmat(i, k) += tempMatrixTerm;
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
	os << endl << endl;
	os << "k-matrix parameters: " << endl;
	for (int i = 0; i < m.kParameters.size(); i++) {
		os << "s^" << i << " coefficient matrix = " << endl << m.kParameters[i] << endl;
	}
	os << endl;
	return os;

}