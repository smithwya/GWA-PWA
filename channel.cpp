#include "channel.h"

channel::channel() {
	masses = {};
	couplings = {};
	chebyCoefficients = {};
	poletype = 1;

}

channel::channel(vector<double> cp, vector<double> ch, vector<double> m) {
	couplings = cp;
	chebyCoefficients = ch;
	masses = m;
	poletype = 3;
	s0 = 1.0;
}

channel::channel(vector<double> masses){
	couplings = {};
	chebyCoefficients = {};
	poletype = 1;
	s0 = 1.0;
	
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