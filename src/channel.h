#pragma once
#include <iostream>
#include <vector>
#include <complex>
typedef std::complex<double> comp;
using namespace std;

class channel {

private:
    string channel_name;
    vector<double> couplings;
    vector<double> couplings_steps;
    vector<double> chebyCoefficients;
    vector<double> chebyCoeff_steps;

    vector<double> masses;
    int poletype;
    double s0;
public:
    channel();
    channel(vector<double> couplings, vector<double> chebyCoefficients, vector<double> m);
    channel(string cname, vector<double> particlemasses);

    comp getMomentum(comp s);
    comp getComplexMomentum(comp s);   
    comp getTrueMomentum(comp s);
    vector<double> getCouplings();
    void setCouplings(vector<double> c);
    double getCoupling(int i);
    
    void setChebyCoeffs(int ptype, double s0, vector<double> c);
    void setChebyCoeffs(int ptype, double s0, vector<double> c, vector<double> csteps);
    vector<double> getChebyCoeffs();
    double getChebyCoeff(int i);
    friend ostream& operator<<(std::ostream& os, channel const& m);
    vector<double> getMasses();
    double getThreshold();
    double getS0();
    void setName(string cname);
    string getName();
    void addCoupling(double coup);
    void addCoupling(double coup,double step);
    int getPoleType();
    void setMasses(vector<double> masses);
    vector<double> getChebySteps();
    vector<double> getCouplingSteps();
};