#include <iostream>
#include<complex>
#include<vector>
#include "amplitude.h"
#include "channel.h"
#include "TCanvas.h"
#include <Eigen/Dense>
#include "TH1D.h"
#include "TFile.h"
#include <sstream>
#include "TString.h"

using namespace std;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;

int main()
{
	//All stuff needed to make amplitude:

	//specify the angular momentum J of the considered waves, parameter s_L (for left-hand singularities),
	// and parameter alpha for controlling asymptotic behavior of integrand when calculating the K-matrix
	comp J1 = comp(0, 0);
	comp J2 = comp(2, 0);
	comp alpha = comp(1., 0);
	comp sL = comp(0.6, 0);
	comp s0 = comp(1, 0);
	comp smin = comp(pow(1., 2), 0);
	comp smax = comp(pow(2.5, 2), 0);

	//specify the masses (in GeV) of the particles in output from each channel (we're considering the case of identical particles in output):

	comp mass1 = comp(0.13498, 0.);
	comp mass2 = comp(0.49761, 0.);
	comp mass3 = comp(0.762, 0.);

	///// 1st wave:

	//Create a channel by passing to the channel object a list of masses (that one defined above), list of couplings and a list of coefficients
	//the list of couplings is the list of {g^R}'s for each channel (make sure you have a g^R for each channel!)
	//the list of coefficients is the list of {a_n}'s in calculating the numerator as a lin combo of chebyshev polynomials 
	//(length of {a_n} = number of chebyshev polynomials to use)

	channel chan1_1= channel({ comp(-6.0584496703077093, 0), comp(-0.23587420969943196, 0), comp(8.2951842306611070, 0)}, { comp(-2173.7290154401403,0.), comp(3272.0517220458223,0.), comp(-1553.7255278262101,0.),  comp(361.78951904504214,0.)}, mass1);
	channel chan2_1 = channel({ comp(-7.2466452617409232, 0), comp(-1.0918936513389781, 0), comp(12.908191623140738, 0)}, { comp(1927.4952345586573, 0.), comp(-2996.8378548561377, 0.), comp(1354.3272596278623, 0.), comp(-287.15761634573755, 0.)}, mass2);
	channel chan3_1 = channel({ comp(1.9557195508523364, 0), comp(-0.23250771304628870, 0), comp(-1.0490788402175895, 0)}, { comp(0, 0) },mass3);
	vector<channel> chans_1 = { chan1_1, chan2_1, chan3_1 };
	vector<double> rmasses_1 = {6.31602747547477250 * pow(10, -3), 3.1835288340590329, 7.5476155809082135};
	
	//Create list of matrices to parameterize the polynomial in the k-matrix.
	//the matrices must be NxN matrices, where N is the number of channels
	MatrixXcd c1 = MatrixXcd::Zero(3, 3);
	MatrixXcd d1 = MatrixXcd::Zero(3, 3);
	c1 << comp(16.190285088168821, 0), comp(18.308149119757218, 0), comp(0, 0), comp(18.308149119757218, 0), comp(19.238824892849152, 0), comp(0, 0), comp(0, 0), comp(0, 0), comp(-15.811018219003017, 0);
	d1 << comp(-6.2739277814634988, 0), comp(-9.1968627819824178, 0), comp(0, 0), comp(-9.1968627819824178, 0), comp(-13.370648723526756, 0), comp(0, 0), comp(0, 0), comp(0, 0), comp(6.0308208392725646, 0);
	vector<MatrixXcd> kparams_1 = { c1 , d1 };

	///// 2nd wave:

	//Create a channel by passing to the channel object a list of masses (that one defined above), list of couplings and a list of coefficients
	//the list of couplings is the list of {g^R}'s for each channel (make sure you have a g^R for each channel!)
	//the list of coefficients is the list of {a_n}'s in calculating the numerator as a lin combo of chebyshev polynomials 
	//(length of {a_n} = number of chebyshev polynomials to use)

	channel chan1_2 = channel({ comp(0.15470117985842080, 0), comp(0.90811894619037048 , 0), comp(0.85226737492303073, 0)}, { comp(109.75915967336491,0.), comp(-2.8444511791990625,0.), comp(71.095211523264609,0.)}, mass1);
	channel chan2_2 = channel({ comp(1.0268866010574129, 0), comp(-0.24961757745677460, 0), comp(0.69753803634648648, 0)}, { comp(-26.963937424908881 , 0.), comp(-135.32369663565123, 0.), comp(-14.887225035532870, 0.)}, mass2);
	channel chan3_2 = channel({ comp(-0.34353840245057654, 0), comp(0, 0), comp(1.1340731281916305, 0)}, { comp(0, 0) },mass3);
	vector<channel> chans_2 = { chan1_2, chan2_2, chan3_2 };
	vector<double> rmasses_2 = {2.3299346099112341, 1.6673851773367765, 2.8936189128830971};

	//Create list of matrices to parameterize the polynomial in the k-matrix.
	//the matrices must be NxN matrices, where N is the number of channels
	MatrixXcd c2 = MatrixXcd::Zero(3, 3);
	MatrixXcd d2 = MatrixXcd::Zero(3, 3);
	c2 << comp(-1.1704615912349254, 0), comp(-1.0116586232616100, 0), comp(0, 0), comp(-1.0116586232616100, 0), comp(0.39464452011270623, 0), comp(0, 0), comp(0, 0), comp(0, 0), comp(-11.365225200233908, 0);
	d2 << comp(0.21138364931721298, 0), comp(0.16825103893825144, 0), comp(0, 0), comp(0.16825103893825144, 0), comp(4.86693780658242758 * pow(10, -3), 0), comp(0, 0), comp(0, 0), comp(0, 0), comp(3.9232752975622134, 0);
	vector<MatrixXcd> kparams_2 = { c2 , d2 };

	//Create the amplitude by passing it the channels, masses and k-matrix parameters
	amplitude wave_1 = amplitude(J1, alpha, sL, chans_1, kparams_1,rmasses_1,s0,smin,smax);

	//Create the amplitude by passing it the channels, masses and k-matrix parameters
	amplitude wave_2 = amplitude(J2, alpha, sL, chans_2, kparams_2,rmasses_2,s0,smin,smax);

	//print the amplitude like this:
	std::cout << "1st-wave amplitude: "<<endl << wave_1 << endl;
	std::cout << "2nd-wave amplitude: "<<endl << wave_2 << endl;

	//print the k matrix at s
	comp s = comp(1,0);

	std::cout << chan1_1.getMomentum(s) << endl;
	std::cout << wave_1.getRhoN(s, 0) << endl;
	std::cout << wave_1.getKMatrix(s) << endl;

	// Plotting:

	double lower_sqrt_s = 1;
	double upper_sqrt_s = 2.5;
	int num_bins = 300;
	double delta = (upper_sqrt_s - lower_sqrt_s)/num_bins;

	//1st channel:

	TH1D *ch1_wv1 = new TH1D("ch1_wv1","ch = 1, J = 0", num_bins, lower_sqrt_s, upper_sqrt_s); //without the normalization constant
	TH1D *ch1_wv2 = new TH1D("ch1_wv2","ch = 1, J = 2", num_bins, lower_sqrt_s, upper_sqrt_s);
	TH1D *integral = new TH1D("integral","ch = 1, J = 0", num_bins, lower_sqrt_s, upper_sqrt_s);
	TH1D *ch1_relative_phases = new TH1D("ch1_relative_phases", "ch1 relative phases SD", num_bins, lower_sqrt_s, upper_sqrt_s);
	comp sqrtS = comp(1.0,0);
	for(int i = 0; i < num_bins; i++){
		sqrtS = lower_sqrt_s + delta/2 + i * delta; //centroid
		ch1_wv1->SetBinContent(i + 1,  real(wave_1.getMomentum(0, pow(sqrtS, 2))) * pow(abs(wave_1.getValue(pow(sqrtS, 2))(0)), 2));
		ch1_wv2->SetBinContent(i + 1,  real(wave_2.getMomentum(0, pow(sqrtS, 2))) * pow(abs(wave_2.getValue(pow(sqrtS, 2))(0)), 2));
		integral->SetBinContent(i + 1, imag(wave_1.getIntegral(pow(sqrtS, 2), 1)));
		ch1_relative_phases->SetBinContent(i + 1, arg(wave_1.getValue(pow(sqrtS, 2))(0)) - arg(wave_2.getValue(pow(sqrtS, 2))(0)));
	}

	TFile file("pdf_folder.root", "recreate");
	TCanvas canv;
	ch1_wv1->Write();
	ch1_wv1->Draw();
	canv.SaveAs("ch1_wv1.pdf");
	ch1_wv2->Write();
	ch1_wv2->Draw();
	canv.SaveAs("ch1_wv2.pdf");
	ch1_relative_phases->Write();
	ch1_relative_phases->Draw();
	canv.SaveAs("ch1_relative_phases.pdf");
	integral->Write();
	integral->Draw();
	canv.SaveAs("integral.pdf");
	file.Close();

	/*vector<int> chan_list = {0, 1, 2}; //user has to modify only these three lines writing the number of each channel considered and...
	vector<amplitude> wv_list = {wave_1, wave_2}; 
	vector<int> J_list = {0, 2}; // ..the waves considered for every channel, and the (#) line for specifing the waves which you want the relative phase of.

	vector<TH1D*> amp_graf; 
	vector<TH1D*> phase_graf;

	char amp_name[10];
	char phase_name[20];

	stringstream dynamic_name1;
	stringstream dynamic_name2;

	TFile file("pdf_folder.root", "recreate");
	TCanvas canv;

	for(int i = 1; i <= chan_list.size(); i++){
		for(int j = 1; j <= J_list.size(); j++){
			dynamic_name1 << "ch" << i << "_wv" << j;
			dynamic_name1 >> amp_name;
			std::cout << amp_name << endl;
			amp_graf.push_back(new TH1D(amp_name, Form("ch = %i, J = %i", chan_list.at(i - 1), J_list.at(j - 1)), num_bins, lower_sqrt_s, upper_sqrt_s));

			for(int k = 0; k < num_bins + 1; k++){
				s = lower_sqrt_s + delta/2 + k * delta; //centroid
				amp_graf.at((i - 1) * J_list.size() + j - 1)->SetBinContent(k + 1, real(wv_list.at(j - 1).getMomentum(j - 1, s)) * pow(abs(wv_list.at(j - 1).getValue(s)(j - 1)), 2));
			}

			amp_graf.at((i - 1) * J_list.size() + j - 1)->Draw();
			amp_graf.at((i - 1) * J_list.size() + j - 1)->Write();
			canv.SaveAs(Form("%s.pdf", amp_name));
			dynamic_name1.clear();
		}
		dynamic_name2 << "ch" << i << "_relative_phases";
		dynamic_name2 >> phase_name;
		std::cout << phase_name << endl;
		phase_graf.push_back(new TH1D(phase_name, Form("ch = %i relative phases", chan_list.at(i - 1)), num_bins, lower_sqrt_s, upper_sqrt_s));

		for(int k = 0; k < num_bins + 1; k++){
			s = lower_sqrt_s + delta/2 + k * delta; //centroid
			phase_graf.at(i - 1)->SetBinContent(k + 1, arg(wv_list.at(0).getValue(s)(i - 1)) - arg(wv_list.at(1).getValue(s)(i - 1))); //(#)
		}

		phase_graf.at(i - 1)->Draw();
		phase_graf.at(i - 1)->Write();
		canv.SaveAs(Form("%s.pdf", phase_name));
		dynamic_name2.clear();
	}*/

	return 0;
}