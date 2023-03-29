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




int main()
{
	filereader testReader("Data/GWA_dataformat.txt");
	testReader.SetAllCommandLists();
	testReader.ConstructBareAmps();
	testReader.setChebys();
	testReader.setPoles();
	testReader.setKmats();
	observable testObs = testReader.getObs();


	amplitude S_wave = testObs.amplitudes[0];
	vector<channel> chans = testObs.amplitudes[0].getChannels();

	vector<double> polesteps = S_wave.getPoleSteps();
	for(double x : polesteps) cout<<x<<endl;




	/*testJPsi.plotIntensity(0,0);
	testJPsi.plotIntensity(0,1);
	testJPsi.plotIntensity(0,2);*/
}
