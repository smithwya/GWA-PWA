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
	//reads the file and creates an observable object with the information from the file
	filereader testReader("Data/GWA_dataformat.txt");
	testReader.SetAllCommandLists();
	testReader.ConstructBareAmps();
	testReader.setChebys();
	testReader.setPoles();
	testReader.setKmats();
	observable testObs = testReader.getObs();

	//grabs the amplitudes from the observable object, and recasts as a JPsi object
	//The JPsi object has information needed for doing the fit
	JPsi testJPsi = JPsi(testObs.amplitudes);
	//plots the S-wave and all channels using the JPsi object
	testJPsi.plotIntensity(0,0);
	testJPsi.plotIntensity(0,1);
	testJPsi.plotIntensity(0,2);
}
