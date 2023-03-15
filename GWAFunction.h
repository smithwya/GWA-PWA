#ifndef MN_GWAFunction_H_
#define MN_GWAFunction_H_

#include <math.h>

class GWAFunction {

public:

  GWAFunction(double mean, double sig, double constant) :
    theMean(mean), theSigma(sig), theConstant(constant) {}

  ~GWAFunction() {}

  double m() const {return theMean;}
  double s() const {return theSigma;}
  double c() const {return theConstant;}

  double operator()(double x) const {

    return
       c()*exp(-0.5*(x-m())*(x-m())/(s()*s()))/(sqrt(2.*M_PI)*s());
       
  }

private:

  double theMean;
  double theSigma;
  double theConstant;
};
#endif // MN_GWAFunction_H_