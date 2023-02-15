#include "channel.h"

channel::channel() {
	couplings = {};
}

channel::channel(vector<comp> cp, vector<comp> ch, comp m) {
	couplings = cp;
	chebyCoefficients = ch;
	mass = m;
}

vector<comp> channel::getCouplings()
{
	return couplings;
}

void channel::setCouplings(vector<comp> c)
{
	couplings = c;
}


comp channel::getCoupling(int i)
{
	return couplings[i];
}

void channel::setChebyCoeffs(vector<comp> c)
{
	chebyCoefficients = c;
}

void channel::setChebyCoeff(int i, comp c)
{
	chebyCoefficients[i] = c;
}

vector<comp> channel::getChebyCoeffs()
{
	return chebyCoefficients;
}

comp channel::getChebyCoeff(int i)
{
	return chebyCoefficients[i];
}

comp channel::getMass()
{
	return mass;
}

void channel::setMass(comp m)
{
	mass = m;
}

comp channel::getMomentum(comp s){
return 0.5*sqrt(s-4.0*pow(mass,2));
}

ostream& operator<<(ostream& os, channel const& m) {
	os << "mass =" << m.mass << endl;
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