#pragma once
#include <iostream>
#include <complex>
#include <algorithm>
#include <math.h>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "channel.h"
#include "amplitude.h"
#include "TMath.h"
#include "TH1.h"
#include "Math/Integrator.h"
#include "Math/Functor.h"
#include <map>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;

typedef std::complex<double> comp;



//amplitude constructor
amplitude::amplitude() {
	channels = {};
	kParameters = {};
	numChannels = 0;
	s0 = 0;
	smin = 0;
	smax = 0;
	integralList = {};
}

amplitude::amplitude(int j, double alp, double ssl, vector<channel> chans, vector<MatrixXcd> kParams, vector<double> rmasses, double ss0, double ssmin, double ssmax) {
	channels = chans;
	kParameters = kParams;
	numChannels = channels.size();
	J = j;
	alpha = alp;
	sL = ssl;
	s0 = ss0;
	smin = ssmin;
	smax = ssmax;
	resmasses = rmasses;
	integralList = {};
}

amplitude::amplitude(int Jj, double ssmin, double ssmax, vector<channel> chans){
	numChannels = chans.size();
	J = Jj;
	channels = chans;
	smin = ssmin;
	smax = ssmax;
	alpha = 1.0;
	sL = 1.0;
	s0 = 1.0;
	integralList = {};
	for(channel c: chans){
		channel_names.push_back(c.getName());
	}

}


	// return E_gamma * p_i * numerator * denominator.inverse()
VectorXcd amplitude::getValue(comp s) {
	
	double m_JPsi = 3.0969;

	comp Egamma = pow((s-pow(m_JPsi,2)),2)/(4.0*s);
	
	MatrixXcd phsp = MatrixXcd::Identity(numChannels,numChannels);


	for(int i = 0; i < numChannels; i++){
		phsp(i,i)=sqrt(Egamma*pow(channels[i].getMomentum(s),2.0*J+1.0));
	}
	
	return (phsp*getDenominator(s).inverse())*(getNumerator(s,3));
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
			s0 = channels[k].getS0();
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
	intKey skpair = intKey(s,k);
	auto mapit = integralList.find(skpair);

	if(mapit!=integralList.end()) return mapit->second;

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
	
	double threshold = channels[k].getThreshold();

	double realpart = intRe.IntegralUp(threshold+.0001);
	double imagpart = intIm.IntegralUp(threshold+.0001);
	
	comp result = comp(realpart,imagpart);


	integralList[skpair]=result;

	return result;
}


void amplitude::calcIntegrals(vector<comp> slist,int k){

	/*for(comp s : slist){
	integralList.insert({make_pair(s,k),0});
	}*/
	

	return;
}

MatrixXcd amplitude::getDenominator(comp s)
{
	MatrixXcd Kinv = getKMatrix(s).inverse();

	MatrixXcd M = MatrixXcd::Zero(numChannels,numChannels);

	for(int k = 0; k < numChannels; k++){
		M(k,k) = getIntegral(s,k);
	}

	return Kinv-M;
}

//integratoe f(function, lower bound, upper bound)

comp amplitude::getMomentum(int chan, comp s)
{
	return channels[chan].getMomentum(s);
}

MatrixXcd amplitude::getKMatrix(comp s) {
	MatrixXcd kmat = MatrixXcd::Zero(numChannels, numChannels);

	for (int k = 0; k < numChannels; k++) {
		for (int i = 0; i <= k; i++) {
			comp tempMatrixTerm = 0;

			for (int R = 0; R < numChannels; R++) {
				tempMatrixTerm += channels[k].getCoupling(R) * channels[i].getCoupling(R) / (resmasses[R] - s);
			}

			kmat(k, i) = tempMatrixTerm;
			kmat(i, k) = tempMatrixTerm;
		}
	}

	for (int j = 0; j < kParameters.size(); j++) {

		kmat += pow(s, j) * kParameters[j];

	}

	return kmat.real();
}

void amplitude::setChebyCoeffs(string cname, int poletype, double s0, vector<double> coeffs){

	auto it = find(channel_names.begin(),channel_names.end(),cname);

	if(it==channel_names.end()){
		cout<<"channel not found"<<endl;
		return;
	} 

	int index = it-channel_names.begin();

	channels[index].setChebyCoeffs(poletype,s0,coeffs);

}


void amplitude::setKParams(int power, vector<vector<double>> kparamlist){
	int nParams = kParameters.size();

	while(kParameters.size()<power+1){
		kParameters.push_back(MatrixXcd::Zero(numChannels,numChannels));
	}
	MatrixXcd newKParam = MatrixXcd::Zero(numChannels,numChannels);

	for(int i = 0; i < numChannels; i++){
		for(int j = numChannels-1; j>=i; j--){
			newKParam(i,j) = kparamlist[i][j-i];
			newKParam(j,i) = newKParam(i,j);

		}
	}
	kParameters[power]=newKParam;
}

void amplitude::addPole(double mass, vector<string> chan_names, vector<double> couplings){

	int nCouplings = chan_names.size();

	if(nCouplings != couplings.size()) return;

	resmasses.push_back(mass);

	for(int i = 0; i < numChannels; i++){

		auto it = find(chan_names.begin(),chan_names.end(),channels[i].getName());

		if(it==chan_names.end()){
			channels[i].addCoupling(0.0);
			continue;
		}

		int index = it-chan_names.begin();
		channels[i].addCoupling(couplings[index]);
	} 
	return;
}


ostream& operator<<(ostream& os, amplitude const& m) {
	os << "J = " << m.J << ", alpha = " << m.alpha << ", sL = " << m.sL <<" num_channels = "<<m.numChannels<<" kmat_mat_params = " <<m.kParameters.size() <<endl;
	os<<"numPoles = "<<m.resmasses.size() <<" s0= "<<m.s0<<" smin = "<<m.smin<<" smax ="<<m.smax<<endl<<endl;

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