#include <iostream>
#include<complex>
#include<vector>
#include "amplitude.h"
#include "channel.h"
#include "TCanvas.h"
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;

int main()
{
	//All stuff needed to make amplitude:
	int x = 3.0;
	//Create list of matrices to parameterize the polynomial in the k-matrix.
	//the matrices must be NxN matrices, where N is the number of channels
	MatrixXcd c = MatrixXcd::Zero(2, 2);
	MatrixXcd d = MatrixXcd::Zero(2, 2);
	c << comp(1, 0), comp(0, 0), comp(0, 0), comp(1, 0);
	d << comp(0, 0), comp(1, 0), comp(1, 0), comp(0, 0);
	vector<MatrixXcd> kparams = { c , d };

	//Create a channel by passing to the channel object a mass, list of couplings and a list of coefficients
	//the list of couplings is the list of {g^R}'s for each channel (make sure you have a g^R for each channel!)
	//the list of coefficients is the list of {a_n}'s in calculating the numerator as a lin combo of chebyshev polynomials 
	//(length of {a_n} = number of chebyshev polynomials to use)
	comp mass1 = comp(1, 0);
	comp mass2 = comp(1, 1);
	channel chan1 = channel({ comp(1.0), comp(1.0) }, { comp(1.0), comp(1.0) },mass1);
	channel chan2 = channel({ comp(1.0), comp(1.0) }, { comp(1.0), comp(1.0) },mass2);
	vector<channel> chans = { chan1, chan2 };

	//specify the angular momentum J, parameter s_L (for left-hand singularities),
	// and parameter alpha for controlling asymptotic behavior of integrand when calculating the K-matrix
	comp J = comp(1.0, 0);
	comp alpha = comp(1.0, 0);
	comp sL = comp(1.0, 0);
	comp s0 = comp(1, 0);
	comp smin = comp(1, 0);
	comp smax = comp(2.5, 0);


	//Create the amplitude by passing it the channels, masses and k-matrix parameters
	amplitude A = amplitude(J, alpha, sL, chans, kparams,s0,smin,smax);

	//print the amplitude like this:
	cout << "Amplitude: "<<endl << A << endl;

	//print the k matrix at s
	comp s = comp(3.4, 0);
	cout << "k-matrix at s = " << s << endl<< A.getKMatrix(s)<<endl;
	//TCanvas *c1 = new TCanvas("canvas","canvas");
	//c1->SaveAs("canvas.pdf");
	return 0;
}