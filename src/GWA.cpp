#include <iostream>
#include <complex>
#include <vector>
#include <Eigen/Dense>
#include <sstream>
#include <fstream>
#include <map>
#include <chrono>
#include "amplitude.h"
#include "channel.h"
#include "filereader.h"
#include "observable.h"
#include "JPsi.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"


using namespace std;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;



vector<amplitude> readInput(string filename){

//do something

return {};
}



int main()
{
	double smin = pow(0.997,2);
	double smax = pow(2.5,2);

	channel c1 = channel("PiPi",{0.13498, 0.13498});
	channel c2 = channel("KK", {0.49761, 0.49761});
	channel c3 = channel("RhoRho", {0.762, 0.762});

	amplitude S_wave = amplitude(0,0.6,smin,smax,{c1,c2,c3});

	S_wave.setChebyCoeffs("PiPi",3, 1.0, {-2173.73, 3272.05, -1553.73, 361.79});
	S_wave.setChebyCoeffs("KK",3, 1.0, {1927.5, -2996.84, 1354.33 , -287.158});
	S_wave.setChebyCoeffs("RhoRho", 3,1.0, {0});

	S_wave.addPole(0.00631603, {"PiPi","KK","RhoRho"},{-6.05845,-7.24665,1.95572});
	S_wave.addPole(3.18353,{"PiPi","KK","RhoRho"},{-0.235874,-1.09189,-0.232508});
	S_wave.addPole(7.54762,{"PiPi","KK","RhoRho"},{8.29518,12.9082,-1.04908});

	S_wave.setKParams(0,{{16.1903, 18.3081 ,0}, {19.2388, 0}, {-15.811}});
	S_wave.setKParams(1, {{-6.27393,-9.19686,0}, {-13.3706,0}, {6.03082}});
	
	amplitude D_wave = amplitude(2,0.6,smin,smax,{c1,c2,c3});


	D_wave.setChebyCoeffs("PiPi", 3, 1.0, {109.759,-2.8444, 71.095,0 });
	D_wave.setChebyCoeffs("KK", 3, 1.0 , {-26.964, -135.324, -14.887,0});
	D_wave.setChebyCoeffs("RhoRho",3,1.0, {0});

	D_wave.addPole( 2.332, {"PiPi","KK","RhoRho"}, {0.155, 1.027, -0.344});
	D_wave.addPole( 1.667, {"PiPi","KK","RhoRho"}, {0.908, -0.250, 0});
	D_wave.addPole( 2.894, {"PiPi","KK","RhoRho"}, {0.852, 0.698, 1.134});

	D_wave.setKParams(0, {{-1.170,-1.012, 0 }, {0.395, 0}, {-11.365}});
	D_wave.setKParams(1, {{0.211, 0.168, 0}, {0.005,0}, {3.923}});
		
	JPsi jObs = JPsi({c1,c2,c3},{S_wave,D_wave});
	//jObs.plotIntensity(0,0);
}