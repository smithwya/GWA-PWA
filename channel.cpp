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

ostream& operator<<(ostream& os, channel const& m) {
	os << "masses =";
	for(double mass : m.masses){
		os<< mass<<" ";
	}
	os<<endl;
	os << "couplings = { ";
	for (int i = 0; i < m.couplings.size(); i++) {
		os << m.couplings[i] << ", ";
	}
	os << "}" << endl;
	os << "chebyCoeffs = { ";
	for (int i = 0; i < m.chebyCoefficients.size(); i++) {
		os << m.chebyCoefficients[i] << ", ";
	}
	os << "}";
	return os;

}