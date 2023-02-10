#include "observable.h"
#include<complex>
#include<algorithm>
#include<math.h>
#include<chrono>
#include<random>
#include<fstream>
#include<vector>
#include<string>
#include <Eigen/Dense>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;



observable::observable()
{
}

observable::observable(vector<channel> c, vector<amplitude> a)
{
}

double observable::getIntensity() {
	return 0;
}

double observable::getPhase() {
	return 0.0;
}

void observable::addChannel(channel c) {

}
void observable::addAmp(amplitude a) {

}
void observable::getAmps() {

}
