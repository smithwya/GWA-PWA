#include <iostream>
#include <vector>
#include <string>
#include "observable.h"
#include <Eigen/Dense>
#include "amplitude.h"
using namespace std;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;


class JPsi : public observable{
    public:

    using observable::observable;


    void plotIntensity(int wav, int ch){
        auto s_intens = [&](double x){
		    comp val = getAmp(wav).getValue(pow(x,2))(ch);
		    return pow(abs(val),2);
        };
	    makePlot("wave"+to_string(wav)+"-"+to_string(ch),s_intens,1.0,2.5,300);
    };
    
};