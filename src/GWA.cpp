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

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"


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
	/*testJPsi.plotIntensity(0,0);
	testJPsi.plotIntensity(0,1);
	testJPsi.plotIntensity(0,2);*/

	//make the minimzer
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2","");
	min->SetMaxFunctionCalls(1000000);
	min->SetMaxIterations(10000);
	min->SetTolerance(0.001);
	min->SetPrintLevel(1);

	amplitude S_wave = testObs.amplitudes[0];
	
	vector<double> steps = S_wave.getStepSizes();
	vector<double> params =S_wave.getParamList();
	vector<double> fittedParams = S_wave.getFittedParamList();
	
	cout<<S_wave<<endl;
	S_wave.setFittedParamList(fittedParams);
	cout<<S_wave<<endl;
	
}
