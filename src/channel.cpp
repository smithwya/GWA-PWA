#include "channel.h"

channel::channel() {
	masses = {};
	couplings = {};
	chebyCoefficients = {};
	poletype = 1;
	channel_name = "no_name";

}

channel::channel(vector<double> cp, vector<double> ch, vector<double> m) {
	couplings = cp;
	chebyCoefficients = ch;
	masses = m;
	poletype = 3;
	s0 = 1.0;
	channel_name = "no_name";
	chebyCoefficients = {0};
}

channel::channel(string cname, vector<double> m){
	couplings = {};
	chebyCoefficients = {};
	poletype = 1;
	s0 = 1.0;
	channel_name = cname;
	masses = m;
}

vector<double> channel::getCouplings()
{
	return couplings;
}

void channel::setCouplings(vector<double> c)
{
	couplings = c;
}


double channel::getCoupling(int i)
{
	return couplings[i];
}

void channel::setChebyCoeffs(int ptype, double ss0, vector<double> c)
{
	chebyCoefficients = c;

	poletype = ptype;
	s0 = ss0;
}

void channel::setChebyCoeffs(int ptype, double ss0, vector<double> c, vector<double> csteps)
{
	chebyCoefficients = c;
	poletype = ptype;
	s0 = ss0;
	chebyCoeff_steps = csteps;
}

vector<double> channel::getChebyCoeffs()
{
	return chebyCoefficients;
}

double channel::getChebyCoeff(int i)
{
	return chebyCoefficients[i];
}

vector<double> channel::getMasses()
{
	return masses;
}

comp channel::getMomentum(comp s){
	double masssq =0;
	for(double m : masses){
		masssq +=m;
	}
return 0.5*sqrt(s-pow(masssq,2));
}

comp channel::getTrueMomentum(comp s){
	return 0.5 * (sqrt((s-pow(masses[0] + masses[1],2))*(s-pow(masses[0] - masses[1],2))))/s;
}

double channel::getThreshold(){
	double masssq =0;
	for(double m : masses){
		masssq +=m;
	}

	return pow(masssq,2);

}

double channel::getS0(){
	return s0;
}

void channel::setName(string cname){
	channel_name = cname;
}

string channel::getName(){
	return channel_name;
}

void channel::addCoupling(double x){
	
	couplings.push_back(x);

	return;
}

void channel::addCoupling(double x,double step){
	
	couplings.push_back(x);
	couplings_steps.push_back(step);

	return;
}

int channel::getPoleType(){

	return poletype;
}

void channel::setMasses(vector<double> ms){

	masses = ms;
	return;
}

vector<double> channel::getChebySteps(){
	return chebyCoeff_steps;
}

vector<double> channel::getCouplingSteps(){
	return couplings_steps;
}

ostream& operator<<(ostream& os, channel const& m) {
	os<<"channel "<<m.channel_name<<":"<<endl;
	os << "masses =";
	for(double mass : m.masses){
		os<< mass<<" ";
	}
	os<<endl;
	os << "couplings = { ";
	for (int i = 0; i < m.couplings.size(); i++) {
		os << m.couplings[i] << ", ";
	}
	os << "}" << endl; // << " and the couplings.size is " << m.couplings.size() << endl;
	os << "chebyCoeffs = { ";
	for (int i = 0; i < m.chebyCoefficients.size(); i++) {
		os << m.chebyCoefficients[i] << ", ";
	}
	os << "}";
	return os;

}