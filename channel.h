#pragma once
#include <iostream>
#include<vector>
#include <complex>
typedef std::complex<double> comp;
using namespace std;

class channel {

private:
    vector<comp> couplings;
    vector<comp> chebyCoefficients;
    comp mass;
public:
    channel();
    channel(vector<comp> couplings, vector<comp> chebyCoefficients, comp m);
    comp getMomentum(int particle, comp s);
    vector<comp> getCouplings();
    void setCouplings(vector<comp> c);
    comp getCoupling(int i);
    void setChebyCoeffs(vector<comp> c);
    void setChebyCoeff(int i, comp c);
    vector<comp> getChebyCoeffs();
    comp getChebyCoeff(int i);
    friend ostream& operator<<(std::ostream& os, channel const& m);
    comp getMass();
    void setMass(comp m);
};