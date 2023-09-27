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
	resmasses_steps = {};
	kParameters_steps = {};

	epsilon = 1e-3;
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
	resmasses_steps = {};
	kParameters_steps = {};

	epsilon = 1e-3;
	name = "";
}

amplitude::amplitude(string ampName, int Jj,double ssL, double ssmin, double ssmax, vector<channel> chans){
	name = ampName;
	numChannels = chans.size();
	J = Jj;
	channels = chans;
	smin = ssmin;
	smax = ssmax;
	alpha = 1.0;
	sL = ssL;
	s0 = 1.0;
	integralList = {};
	resmasses_steps = {};
	kParameters_steps = {};

	for(channel c: chans){
		channel_names.push_back(c.getName());
	}
	epsilon = 1e-3;
}

vector<channel> amplitude::getChannels(){

	return channels;
}

vector<double> amplitude::getResMasses(){
	return resmasses;
}

	// return E_gamma * p_i * numerator * denominator.inverse()
VectorXcd amplitude::getValue(comp s) {

	bool sheet = false; //because here we are on the real axis to calculate the amplitudes
	
	/*
	//comp Egamma = (pow(m_bottomonium,2)-s)/(2.0*sqrt(s));
	MatrixXcd phsp = MatrixXcd::Identity(numChannels,numChannels);


	for(int i = 0; i < numChannels; i++){
		//phsp(i,i)=Egamma*pow(getMomentum(i,s),J+0.5);
		phsp(i,i)= pow(getMomentum(i,s),J+0.5)/pow(s,.25);
	}
	
	//return (phsp*getDenominator(s).inverse())*(getNumerator(s,3));
	return ((getNumerator(s,3).transpose())*(getDenominator(s).inverse()))*phsp;
	*/

	
	MatrixXcd phsp = MatrixXcd::Identity(numChannels,numChannels);

	//s = comp(118.5,0.8);
	for(int i = 0; i < numChannels; i++){
		phsp(i,i)= pow(getMomentum(i,s),J+0.5)/pow(s,.25);
	}

	comp test = 0;
	//cout << channels[2].getMasses()[0] << "	" << channels[2].getMasses()[1]<<endl;
	//cout << getKMatrix(comp(118.5,0.8)) << endl; exit(0);
	//test = omega_s(comp(118.5,0.8)).imag(); cout << test << endl; exit(0);
	//cout << getNumerator(comp(118.5,0.8), channels[0].getPoleType()) << endl; exit(0);
	//test = pow(2. * channels[2].getMomentum(comp(118.5,0.8)), 2 + 1); cout << test.real() << "	" << test.imag() << endl; exit(0);
	//test = getRhoN(comp(118.5,0.8), 2); cout << test.real() << "	" << test.imag() << endl; exit(0);
	cout << getIntegral(comp(118.5,0.8), 2, sheet) << endl << endl;
	//cout << phsp << endl; exit(0);
	//cout << (getKMatrix(comp(118.5,0.8)).inverse()) << endl; exit(0);
	//cout << getDenominator(comp(118.5,0.8)) << endl; exit(0);
	cout << getDenominator(comp(118.5,0.8), sheet).inverse() << endl; 
	
	MatrixXcd K = getKMatrix(s);
	MatrixXcd I = MatrixXcd::Identity(numChannels,numChannels);

	MatrixXcd DispRhoN = MatrixXcd::Zero(numChannels,numChannels);

	for(int k = 0; k < numChannels; k++){
		DispRhoN(k,k) = getIntegral(s,k,sheet);
	}

	return (getNumerator(s, channels[0].getPoleType()).transpose() * (I - K * DispRhoN).inverse() * K) * phsp;
	//poletype is the same for every wave and every channel, so I take the 0th
	

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

	return 2. * (omega_p(s) - omega_p(smin))/(omega_p(smax) - omega_p(smin)) - 1.0;
};



comp amplitude::omega(comp s, int type){

	switch (type)
	{
	case 1: return omega_p(s);
	case 2: return omega_s(s);
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
comp amplitude::getRhoN(comp sprime,int k,bool sheet)
{
	int dumbJ = J;
	//if(k==2) dumbJ = 0;

	//cout << "dumbJ = " << J << endl;

	comp x = pow(2.0*channels[k].getMomentum(sprime),2.0*dumbJ+1.0)/pow(sprime + sL,dumbJ+alpha);

	if(sheet) x = pow(-1., dumbJ) * pow(2.0*channels[k].getComplexMomentum(sprime),2.0*dumbJ+1.0)/pow(sprime + sL,dumbJ+alpha);

	return x;
}

//overload later for other sheets: no, see below
comp amplitude::getIntegrand(double sp,comp s, int k){

	bool sheet = false; //because we need this function just to calculate the integral on the first sheet (see "GetIntegral()")

	comp z = getRhoN(sp,k,sheet)*s/(sp*(sp-s-comp(0,1)*epsilon)*TMath::Pi());
	//return z.real() + comp(0,1)*getRhoN(s,k);
	return z;
}


comp amplitude::getIntegral(comp s,int k,bool sh){

	comp result = 0;

	if(sh) result += 2. * getRhoN(s,k,sh);

	//cout << "result before = " << result << endl; // .
	
	//add parameter to select sheet, "+-++" etc
	intKey skpair = intKey(s,k,sh);
	auto mapit = integralList.find(skpair);

	//if(mapit!=integralList.end()) return mapit->second; //anyway for the polesearch there's no need to store the result 

	if(abs(s.imag()) < 2 * epsilon){
		double sreal = s.real();
		//result += getIntegral(sreal, k);
		comp result = getIntegral(sreal, k);
		integralList[skpair] = result;
		return result;
	}

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

	double realpart = intRe.IntegralUp(threshold+.001);
	double imagpart = intIm.IntegralUp(threshold+.001);
	
	//comp result = comp(realpart,imagpart);
	result += comp(realpart,imagpart);

	integralList[skpair]=result; // if(sh) cout << "ehilà" << endl; // .

	// cout << "result after = " << result << endl; // .

	//if(sh) result += 2. * getRhoN(s,k,sh);

	return result;
}

comp amplitude::getIntegral(double s,int k){

	bool sheet = false; //because the value of the amplitude (and therefore also the dispersive integral) on the real axis 
	// is the same for every sheet (in fact the unphyhsical sheets fuses with the physical one in the real axis)

	double tempRhoN = getRhoN(s, k, sheet).real();

	auto integrand = [&](double sp)
	{
		return s * (getRhoN(sp, k, sheet).real() - tempRhoN)/(sp * (sp - s)) / TMath::Pi(); //if sp is real getRhoN(sp,k) is real too
	};

	ROOT::Math::Functor1D integrand_fun(integrand);

	ROOT::Math::Integrator integral(integrand_fun,ROOT::Math::IntegrationOneDim::kGAUSS,1.E-9,1.E-6);

	double threshold = channels[k].getThreshold();

	double integral_result = integral.IntegralUp(threshold + epsilon);

	comp temp = integral_result + tempRhoN * log(threshold / abs(s - threshold)) / TMath::Pi();

	if(s > threshold){

		temp += comp(0,tempRhoN); 

	} 

	return temp;

}


void amplitude::calcIntegrals(vector<comp> slist,int k,bool sheet){

	for(comp s : slist){
		integralList[intKey(s,k,sheet)]=getIntegral(s,k,sheet);
	}
	

	return;
}

MatrixXcd amplitude::getDenominator(comp s,bool sheet)
{
	MatrixXcd Kinv = getKMatrix(s).inverse();

	MatrixXcd M = MatrixXcd::Zero(numChannels,numChannels);

	for(int k = 0; k < numChannels; k++){
		M(k,k) = getIntegral(s,k,sheet);
	}

	return Kinv-M;
}

//integratoe f(function, lower bound, upper bound)

comp amplitude::getMomentum(int chan, comp s)
{
	return channels[chan].getMomentum(s);
}

comp amplitude::getComplexMomentum(int chan, comp s)
{
	return channels[chan].getComplexMomentum(s);
}

comp amplitude::getTrueMomentum(int chan, comp s)
{
	return channels[chan].getTrueMomentum(s);
}

void amplitude::setResMasses(vector<double> rm){
	resmasses = rm;
	return;
}

MatrixXcd amplitude::getKMatrix(comp s) {
//return MatrixXcd::Ones(numChannels, numChannels);

	MatrixXcd kmat = MatrixXcd::Zero(numChannels, numChannels);

	for (int k = 0; k < numChannels; k++) {
		for (int i = 0; i <= k; i++) {
			comp tempMatrixTerm = 0;

			for (int R = 0; R < resmasses.size(); R++) {
				tempMatrixTerm += channels[k].getCoupling(R) * channels[i].getCoupling(R) / (resmasses[R] - s);
			}

			kmat(k, i) = tempMatrixTerm;
			kmat(i, k) = tempMatrixTerm;
		}
	}

	for (int j = 0; j < kParameters.size(); j++) {

		kmat += pow(s, j) * kParameters[j];

	}

	//return kmat.real();
	return kmat;
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

void amplitude::setChebyCoeffs(string cname, int poletype, double s0, vector<double> coeffs,vector<double> csteps){

	auto it = find(channel_names.begin(),channel_names.end(),cname);

	if(it==channel_names.end()){
		cout<<"channel not found"<<endl;
		return;
	} 

	int index = it-channel_names.begin();

	channels[index].setChebyCoeffs(poletype,s0,coeffs,csteps);

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

void amplitude::setKParams(int power, vector<vector<double>> kparamlist, vector<double> ksteps){
	int nParams = kParameters.size();

	while(kParameters.size()<power+1){
		kParameters.push_back(MatrixXcd::Zero(numChannels,numChannels));
	}
	while(kParameters_steps.size()<power+1){
		kParameters_steps.push_back({});
	}
	MatrixXcd newKParam = MatrixXcd::Zero(numChannels,numChannels);

	for(int i = 0; i < numChannels; i++){
		for(int j = numChannels-1; j>=i; j--){
			newKParam(i,j) = kparamlist[i][j-i];
			newKParam(j,i) = newKParam(i,j);

		}
	}
	kParameters[power] = newKParam;
	kParameters_steps[power] = ksteps;
}

vector<double> amplitude::getKSteps(int i){
	if(i<0 || i > kParameters_steps.size()) return {};

	return kParameters_steps[i];
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

void amplitude::addPole(double mass,double mass_step, vector<string> chan_names, vector<double> couplings, vector<double> steps){

	int nCouplings = chan_names.size();

	if(nCouplings != numChannels){
		cout << "(in amplitude.cpp) nCouplings != numChannels" << endl;
		return;
	}

	//cout << nCouplings << " " << numChannels << endl;

	resmasses.push_back(mass);
	resmasses_steps.push_back(mass_step);

	for(int i = 0; i < numChannels; i++){

		auto it = find(chan_names.begin(),chan_names.end(),channels[i].getName());

		if(it==chan_names.end()){
			channels[i].addCoupling(0.0,0.0);
			continue;
		}

		int index = it-chan_names.begin();
		channels[i].addCoupling(couplings[index],steps[index]);
	} 
	return;
}

vector<double> amplitude::getStepSizes(){
	vector<double> stepsizes = {};
	//J
	stepsizes.push_back(0);
	//alpha
	stepsizes.push_back(0);
	//sL
	stepsizes.push_back(0);

	//channels
	for(channel c : channels){

		vector<double> templist = c.getCouplingSteps();
		stepsizes.insert(stepsizes.end(),templist.begin(),templist.end());

		templist = c.getChebySteps();
		stepsizes.insert(stepsizes.end(),templist.begin(),templist.end());

		//assuming masses are fixed
		for(double m : c.getMasses()){
			stepsizes.push_back(0);
		}
		//s0
		stepsizes.push_back(0);
	}

	//KMatrix steps
	for(vector<double> kp : kParameters_steps){
	stepsizes.insert(stepsizes.end(),kp.begin(),kp.end());
	}

	stepsizes.insert(stepsizes.end(),resmasses_steps.begin(),resmasses_steps.end());

	return stepsizes;
}

vector<double> amplitude::getParamList(){
	vector<double> params = {};

	params.push_back(J);
	params.push_back(alpha);
	params.push_back(sL);


	for(channel c: channels){
		vector<double> coupls = c.getCouplings();
		vector<double> chebys = c.getChebyCoeffs();
		vector<double> masses = c.getMasses();

		params.insert(params.end(),coupls.begin(),coupls.end());
		params.insert(params.end(),chebys.begin(),chebys.end());
		params.insert(params.end(),masses.begin(),masses.end());
		params.push_back(c.getS0());
	}

	for(MatrixXcd m : kParameters){

		for(int i = 0; i < numChannels; i++){
			for(int j = i; j < numChannels; j++){
				params.push_back(m(i,j).real());
				//params.push_back(m(i,j).imag());
			}
			
		}
	}

	params.insert(params.end(),resmasses.begin(),resmasses.end());

	return params;
}


void amplitude::setParamList(vector<double> params){
	J = params[0];
	alpha = params[1];
	sL = params[2];

	int index = 3;

	for(int i = 0; i < numChannels; i++){
		
		int Ncoupls = channels[i].getCouplings().size();
		int Nchebys = channels[i].getChebyCoeffs().size();
		int Nmasses = channels[i].getMasses().size();
		vector<double> coupls = {};
		vector<double> chebys = {};
		vector<double> masses = {};
		double ss0 = 0;

		for(int start = index; index<start+Ncoupls; index++){
			coupls.push_back(params[index]);
		}
		for(int start = index; index<start+Nchebys; index++){
			chebys.push_back(params[index]);
		}
		for(int start = index; index<start+Nmasses; index++){
			masses.push_back(params[index]);
		}

		ss0 = params[index];
		index++;

		setChebyCoeffs(channels[i].getName(),channels[i].getPoleType(),ss0,chebys);
		channels[i].setCouplings(coupls);
		channels[i].setMasses(masses);
	}

	for(int k = 0; k < kParameters.size(); k++){

		for(int i = 0; i < numChannels; i++){
			for(int j = i; j < numChannels; j++){
				double repart = params[index];
				index++;
				double impart = 0;
				//index++;
				comp matEl = comp(repart,impart);
				kParameters[k](i,j) = matEl;
				kParameters[k](j,i) = matEl;
			}

		}
	}

	int nResMasses = resmasses.size();


	vector<double> rm = {};
	for(int start = index; index<start+nResMasses; index++){
		rm.push_back(params[index]);
	}
	setResMasses(rm);
}

vector<double> amplitude::getFittedParamList(){
	vector<double> allParams = getParamList();
	vector<double> stepsizes = getStepSizes();
	vector<double> fittedParams = {};
	
	for(int i = 0; i < allParams.size(); i++){
		if(stepsizes[i]!=0) fittedParams.push_back(allParams[i]);
	}

	return fittedParams;
}

vector<double> amplitude::getFittedSteps(){
	vector<double> fsteps = {};

	for(double x : getStepSizes()){
		if(x!=0) fsteps.push_back(x);
	}
	return fsteps;
}


void amplitude::setFittedParamList(vector<double> fittedParams){
	vector<double> allParams = getParamList();
	vector<double> stepsizes = getStepSizes();

	int index = 0;
	for(int i = 0; i < allParams.size(); i++){

		if(stepsizes[i] != 0){ 
			allParams[i] = fittedParams[index];
			index++;
		}
	}

	setParamList(allParams);
}

string amplitude::getName(){
	return name;
}

vector<double> amplitude::getPoleSteps(){
	return resmasses_steps;
}

vector<string> amplitude::getChanNames(){
	vector<string> chlist = {};
	for(channel c : channels) chlist.push_back(c.getName());
	return chlist;
}

int amplitude::getNumOfChans(){
	return getChanNames().size();
}

vector<double> amplitude::getFitInterval(){
	return {smin, smax};
}

vector<double> amplitude::getResMassesSteps(){
	return resmasses_steps;
}

vector<MatrixXcd> amplitude::getkParameters(){
	return kParameters;
}

ostream& operator<<(ostream& os, amplitude const& m) {
	os << m.name<<"-wave: J = " << m.J << ", alpha = " << m.alpha << ", sL = " << m.sL <<" num_channels = "<<m.numChannels<<" kmat_mat_params = " <<m.kParameters.size() <<endl;
	os<<"numPoles = "<<m.resmasses.size() <<" s0= "<<m.s0<<" smin = "<<m.smin<<" smax ="<<m.smax<<endl<<endl;

	for (int i = 0; i < m.numChannels; i++) {
		os << m.channels[i] << endl<<endl;
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